package main

import (
	"fmt"
	"image"
	"image/color"
	"image/png"
	"math"
	"math/cmplx"
	"os"
	"runtime"
	"sync"
	"time"

	"github.com/hajimehoshi/ebiten/v2"
	"github.com/hajimehoshi/ebiten/v2/inpututil"
	"github.com/hajimehoshi/ebiten/v2/text"
	"github.com/hajimehoshi/ebiten/v2/vector"
	"golang.org/x/image/font/basicfont"
)

// ──────────────────────────────────────────────
// Constants
// ──────────────────────────────────────────────

const (
	screenWidth  = 1200
	screenHeight = 720

	panelSize = 512
	panelX1   = 40
	panelY1   = 40
	panelX2   = 608

	defaultN      = 512
	defaultLambda = 550e-9 // m (green)
	defaultL      = 0.5    // m
	defaultDx     = 1e-6   // m per aperture pixel

	lambdaR = 650e-9
	lambdaG = 550e-9
	lambdaB = 450e-9
)

// ──────────────────────────────────────────────
// FFT (Cooley–Tukey radix-2 DIT)
// ──────────────────────────────────────────────

// twiddleTable caches twiddle factors per N (power of two).
var (
	twiddleCache = map[int][]complex128{}
	twiddleMu    sync.Mutex
)

func twiddles(N int) []complex128 {
	twiddleMu.Lock()
	defer twiddleMu.Unlock()
	if t, ok := twiddleCache[N]; ok {
		return t
	}
	t := make([]complex128, N/2)
	for k := 0; k < N/2; k++ {
		theta := -2 * math.Pi * float64(k) / float64(N)
		t[k] = complex(math.Cos(theta), math.Sin(theta))
	}
	twiddleCache[N] = t
	return t
}

// bitReverse permutes data in place to bit-reversed order.
func bitReverse(data []complex128) {
	N := len(data)
	j := 0
	for i := 1; i < N; i++ {
		bit := N >> 1
		for ; j&bit != 0; bit >>= 1 {
			j ^= bit
		}
		j ^= bit
		if i < j {
			data[i], data[j] = data[j], data[i]
		}
	}
}

// fft1D performs in-place forward FFT. len(data) must be a power of two.
// tw must be twiddles(N).
func fft1D(data []complex128, tw []complex128) {
	N := len(data)
	bitReverse(data)
	for size := 2; size <= N; size <<= 1 {
		half := size >> 1
		step := N / size
		for i := 0; i < N; i += size {
			for k := 0; k < half; k++ {
				t := tw[k*step] * data[i+k+half]
				u := data[i+k]
				data[i+k] = u + t
				data[i+k+half] = u - t
			}
		}
	}
}

// fft2D performs in-place 2D FFT on N×N row-major data using parallel workers.
func fft2D(data []complex128, N int) {
	tw := twiddles(N)
	nWorkers := runtime.NumCPU()
	if nWorkers > N {
		nWorkers = N
	}

	// Phase 1: FFT each row.
	var wg sync.WaitGroup
	rowsPerWorker := (N + nWorkers - 1) / nWorkers
	for w := 0; w < nWorkers; w++ {
		start := w * rowsPerWorker
		end := start + rowsPerWorker
		if end > N {
			end = N
		}
		if start >= end {
			break
		}
		wg.Add(1)
		go func(s, e int) {
			defer wg.Done()
			for r := s; r < e; r++ {
				fft1D(data[r*N:(r+1)*N], tw)
			}
		}(start, end)
	}
	wg.Wait()

	// Phase 2: FFT each column (gather into buffer, FFT, scatter back).
	for w := 0; w < nWorkers; w++ {
		start := w * rowsPerWorker
		end := start + rowsPerWorker
		if end > N {
			end = N
		}
		if start >= end {
			break
		}
		wg.Add(1)
		go func(s, e int) {
			defer wg.Done()
			buf := make([]complex128, N)
			for c := s; c < e; c++ {
				for i := 0; i < N; i++ {
					buf[i] = data[i*N+c]
				}
				fft1D(buf, tw)
				for i := 0; i < N; i++ {
					data[i*N+c] = buf[i]
				}
			}
		}(start, end)
	}
	wg.Wait()
}

// ──────────────────────────────────────────────
// Aperture generators
// ──────────────────────────────────────────────

const (
	apRect      = 0
	apCircle    = 1
	apDoubleSlt = 2
	apGrating   = 3
	apCross     = 4
	apRing      = 5
	apPaint     = 6
)

// buildAperture writes t(x,y)·(-1)^(m+n) into data (length N²).
// dx is the physical step. apW, apH, apD are aperture parameters in meters.
// For grating: apW = slit width, apD = period, apNSlits = number of slits.
func buildAperture(data []complex128, N int, dx float64, apType int,
	apW, apH, apD float64, apNSlits int, paint []float64) {

	nWorkers := runtime.NumCPU()
	if nWorkers > N {
		nWorkers = N
	}
	var wg sync.WaitGroup
	rowsPerWorker := (N + nWorkers - 1) / nWorkers
	for w := 0; w < nWorkers; w++ {
		start := w * rowsPerWorker
		end := start + rowsPerWorker
		if end > N {
			end = N
		}
		if start >= end {
			break
		}
		wg.Add(1)
		go func(s, e int) {
			defer wg.Done()
			cx := float64(N) / 2
			cy := float64(N) / 2
			halfW := apW / 2
			halfH := apH / 2
			for j := s; j < e; j++ {
				y := (float64(j) - cy) * dx
				for i := 0; i < N; i++ {
					x := (float64(i) - cx) * dx
					v := 0.0
					switch apType {
					case apRect:
						if math.Abs(x) <= halfW && math.Abs(y) <= halfH {
							v = 1
						}
					case apCircle:
						if x*x+y*y <= halfW*halfW {
							v = 1
						}
					case apDoubleSlt:
						// two slits of width apW, separation apD, height apH
						if math.Abs(y) <= halfH {
							if math.Abs(math.Abs(x)-apD/2) <= apW/2 {
								v = 1
							}
						}
					case apGrating:
						if math.Abs(y) <= halfH {
							// nSlits slits, period apD, slit width apW
							totalSpan := float64(apNSlits-1) * apD
							leftEdge := -totalSpan / 2
							rel := x - leftEdge
							if rel >= -apW/2 && rel <= totalSpan+apW/2 {
								// nearest slit index
								idx := math.Floor(rel/apD + 0.5)
								if idx >= 0 && idx < float64(apNSlits) {
									slitCenter := leftEdge + idx*apD
									if math.Abs(x-slitCenter) <= apW/2 {
										v = 1
									}
								}
							}
						}
					case apCross:
						// horizontal bar of half-height apH and half-width apW;
						// vertical bar of half-width apH and half-height apW
						if (math.Abs(x) <= halfW && math.Abs(y) <= halfH) ||
							(math.Abs(x) <= halfH && math.Abs(y) <= halfW) {
							v = 1
						}
					case apRing:
						// annulus: outer radius apW/2, inner radius apD/2
						r2 := x*x + y*y
						if r2 <= halfW*halfW && r2 >= (apD/2)*(apD/2) {
							v = 1
						}
					case apPaint:
						if paint != nil {
							v = paint[j*N+i]
						}
					}
					// checkerboard sign to fold in fftshift
					if (i+j)&1 == 1 {
						v = -v
					}
					data[j*N+i] = complex(v, 0)
				}
			}
		}(start, end)
	}
	wg.Wait()
}

// ──────────────────────────────────────────────
// Display transforms
// ──────────────────────────────────────────────

const (
	dispLog    = 0
	dispGamma  = 1
	dispLinear = 2
)

// transformIntensity maps raw |U|² to [0,1] visible range.
func transformIntensity(I, Imax float64, mode int, gamma, logAlpha float64) float64 {
	if Imax <= 0 {
		return 0
	}
	v := I / Imax
	if v < 0 {
		v = 0
	}
	if v > 1 {
		v = 1
	}
	switch mode {
	case dispLog:
		return math.Log1p(logAlpha*v) / math.Log1p(logAlpha)
	case dispGamma:
		return math.Pow(v, 1.0/gamma)
	default:
		return v
	}
}

// ──────────────────────────────────────────────
// Plasma colormap (copied from modeling-1)
// ──────────────────────────────────────────────

type colorStop struct {
	t, r, g, b float64
}

var plasmaStops = []colorStop{
	{0.00, 13, 2, 33},
	{0.05, 26, 10, 78},
	{0.15, 59, 31, 142},
	{0.30, 107, 63, 160},
	{0.45, 160, 81, 149},
	{0.60, 212, 80, 135},
	{0.75, 249, 93, 106},
	{0.88, 255, 158, 59},
	{1.00, 255, 236, 110},
}

func plasmaColor(t float64) (r, g, b uint8) {
	if t <= 0 {
		return 13, 2, 33
	}
	if t >= 1 {
		return 255, 236, 110
	}
	for i := 0; i < len(plasmaStops)-1; i++ {
		if t <= plasmaStops[i+1].t {
			f := (t - plasmaStops[i].t) / (plasmaStops[i+1].t - plasmaStops[i].t)
			return uint8(plasmaStops[i].r + f*(plasmaStops[i+1].r-plasmaStops[i].r)),
				uint8(plasmaStops[i].g + f*(plasmaStops[i+1].g-plasmaStops[i].g)),
				uint8(plasmaStops[i].b + f*(plasmaStops[i+1].b-plasmaStops[i].b))
		}
	}
	return 255, 236, 110
}

// ──────────────────────────────────────────────
// Bilinear sampling on N×N grid (for RGB rescaling)
// ──────────────────────────────────────────────

func bilinear(I []float64, N int, x, y float64) float64 {
	if x < 0 || y < 0 || x > float64(N-1) || y > float64(N-1) {
		return 0
	}
	i := int(x)
	j := int(y)
	fx := x - float64(i)
	fy := y - float64(j)
	i1 := i + 1
	j1 := j + 1
	if i1 >= N {
		i1 = N - 1
	}
	if j1 >= N {
		j1 = N - 1
	}
	a := I[j*N+i]
	b := I[j*N+i1]
	c := I[j1*N+i]
	d := I[j1*N+i1]
	return a*(1-fx)*(1-fy) + b*fx*(1-fy) + c*(1-fx)*fy + d*fx*fy
}

// ──────────────────────────────────────────────
// Game state
// ──────────────────────────────────────────────

type Game struct {
	// Physical parameters
	lambda float64 // wavelength (m)
	L      float64 // distance to screen (m)
	dx     float64 // grid step (m)
	N      int     // grid size (power of 2)

	// Aperture
	apType   int
	apW      float64 // m
	apH      float64 // m
	apD      float64 // m
	apNSlits int

	// Paint mode
	paintBuf  []float64
	painting  bool
	lastPaint image.Point
	hasLast   bool

	// Display
	colorMode   int     // 0 mono, 1 RGB
	displayMode int     // 0 log, 1 gamma, 2 linear
	gamma       float64 // for gamma mode
	logAlpha    float64 // for log mode

	// Internal buffers (sized for current N)
	apertureBuf []complex128
	fftBuf      []complex128
	intensity   []float64

	// Computed images
	apImage   *ebiten.Image
	diffImage *ebiten.Image

	// Stats
	lastRecomputeMs float64
	dirty           bool

	// Input
	holdFrames map[ebiten.Key]int
}

func newGame() *Game {
	g := &Game{
		lambda:      defaultLambda,
		L:           defaultL,
		dx:          defaultDx,
		N:           defaultN,
		apType:      apRect,
		apW:         60e-6,
		apH:         40e-6,
		apD:         20e-6,
		apNSlits:    8,
		colorMode:   0,
		displayMode: dispLog,
		gamma:       2.4,
		logAlpha:    1e3,
		dirty:       true,
		holdFrames:  make(map[ebiten.Key]int),
	}
	g.allocBuffers()
	return g
}

func (g *Game) allocBuffers() {
	n2 := g.N * g.N
	g.apertureBuf = make([]complex128, n2)
	g.fftBuf = make([]complex128, n2)
	g.intensity = make([]float64, n2)
	if g.apType == apPaint {
		if len(g.paintBuf) != n2 {
			g.paintBuf = make([]float64, n2)
		}
	}
}

func (g *Game) keyRepeat(key ebiten.Key) bool {
	if ebiten.IsKeyPressed(key) {
		g.holdFrames[key]++
		f := g.holdFrames[key]
		return f == 1 || (f > 20 && f%4 == 0)
	}
	g.holdFrames[key] = 0
	return false
}

// ──────────────────────────────────────────────
// Recompute pipeline
// ──────────────────────────────────────────────

func (g *Game) recompute() {
	t0 := time.Now()

	// Step 1: build aperture with checkerboard sign
	buildAperture(g.apertureBuf, g.N, g.dx, g.apType, g.apW, g.apH, g.apD, g.apNSlits, g.paintBuf)

	// Step 2: copy to FFT working buffer
	copy(g.fftBuf, g.apertureBuf)

	// Step 3: 2D FFT
	fft2D(g.fftBuf, g.N)

	// Step 4: intensity |U|² (FFT layout: index k maps to physical X = (k - N/2) * DX_λ;
	// thanks to the checkerboard pre-multiplication, the zero-frequency component
	// already lands at index N/2 — no fftshift required.)
	n2 := g.N * g.N
	parallelFor(n2, func(idx int) {
		c := g.fftBuf[idx]
		re, im := real(c), imag(c)
		g.intensity[idx] = re*re + im*im
	})

	// Step 5: max for normalization
	Imax := 0.0
	for _, v := range g.intensity {
		if v > Imax {
			Imax = v
		}
	}
	if Imax <= 0 {
		Imax = 1
	}

	// Step 6: render aperture and diffraction images
	g.renderApertureImage()
	g.renderDiffractionImage(Imax)

	g.lastRecomputeMs = float64(time.Since(t0).Microseconds()) / 1000.0
}

// parallelFor splits [0,n) over CPUs.
func parallelFor(n int, fn func(i int)) {
	nWorkers := runtime.NumCPU()
	if nWorkers > n {
		nWorkers = n
	}
	if nWorkers <= 1 {
		for i := 0; i < n; i++ {
			fn(i)
		}
		return
	}
	chunk := (n + nWorkers - 1) / nWorkers
	var wg sync.WaitGroup
	for w := 0; w < nWorkers; w++ {
		start := w * chunk
		end := start + chunk
		if end > n {
			end = n
		}
		if start >= end {
			break
		}
		wg.Add(1)
		go func(s, e int) {
			defer wg.Done()
			for i := s; i < e; i++ {
				fn(i)
			}
		}(start, end)
	}
	wg.Wait()
}

func (g *Game) renderApertureImage() {
	N := g.N
	pix := make([]byte, N*N*4)
	parallelFor(N, func(j int) {
		for i := 0; i < N; i++ {
			// extract original aperture magnitude (undo checkerboard sign)
			c := g.apertureBuf[j*N+i]
			v := real(c)
			if (i+j)&1 == 1 {
				v = -v
			}
			if v < 0 {
				v = 0
			}
			b := uint8(v * 255)
			idx := (j*N + i) * 4
			pix[idx] = b
			pix[idx+1] = b
			pix[idx+2] = b
			pix[idx+3] = 255
		}
	})
	if g.apImage == nil || g.apImage.Bounds().Dx() != N {
		g.apImage = ebiten.NewImage(N, N)
	}
	g.apImage.WritePixels(pix)
}

func (g *Game) renderDiffractionImage(Imax float64) {
	N := g.N
	pix := make([]byte, N*N*4)

	if g.colorMode == 0 {
		// Mono: 1:1 mapping, plasma colormap
		parallelFor(N, func(j int) {
			for i := 0; i < N; i++ {
				v := transformIntensity(g.intensity[j*N+i], Imax, g.displayMode, g.gamma, g.logAlpha)
				r, gg, b := plasmaColor(v)
				idx := (j*N + i) * 4
				pix[idx] = r
				pix[idx+1] = gg
				pix[idx+2] = b
				pix[idx+3] = 255
			}
		})
	} else {
		// RGB: sample intensity at scale factor lambdaR/lambda for each channel.
		// Display window matches lambdaR (red just fills the panel).
		half := float64(N) / 2
		lambdas := [3]float64{lambdaR, lambdaG, lambdaB}
		parallelFor(N, func(j int) {
			for i := 0; i < N; i++ {
				var ch [3]uint8
				for c := 0; c < 3; c++ {
					scale := lambdaR / lambdas[c]
					kx := half + (float64(i)-half)*scale
					ky := half + (float64(j)-half)*scale
					val := bilinear(g.intensity, N, kx, ky)
					t := transformIntensity(val, Imax, g.displayMode, g.gamma, g.logAlpha)
					ch[c] = uint8(t * 255)
				}
				idx := (j*N + i) * 4
				pix[idx] = ch[0]
				pix[idx+1] = ch[1]
				pix[idx+2] = ch[2]
				pix[idx+3] = 255
			}
		})
	}

	if g.diffImage == nil || g.diffImage.Bounds().Dx() != N {
		g.diffImage = ebiten.NewImage(N, N)
	}
	g.diffImage.WritePixels(pix)
}

// ──────────────────────────────────────────────
// Input handling
// ──────────────────────────────────────────────

func (g *Game) Update() error {
	shift := ebiten.IsKeyPressed(ebiten.KeyShift)

	// Aperture presets
	for key, idx := range map[ebiten.Key]int{
		ebiten.KeyDigit1: apRect,
		ebiten.KeyDigit2: apCircle,
		ebiten.KeyDigit3: apDoubleSlt,
		ebiten.KeyDigit4: apGrating,
		ebiten.KeyDigit5: apCross,
		ebiten.KeyDigit6: apRing,
	} {
		if inpututil.IsKeyJustPressed(key) {
			g.apType = idx
			g.dirty = true
		}
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyM) {
		g.apType = apPaint
		if len(g.paintBuf) != g.N*g.N {
			g.paintBuf = make([]float64, g.N*g.N)
		}
		g.dirty = true
	}

	// Wavelength
	if g.keyRepeat(ebiten.KeyEqual) || g.keyRepeat(ebiten.KeyKPAdd) {
		g.lambda += 10e-9
		if g.lambda > 780e-9 {
			g.lambda = 780e-9
		}
		g.dirty = true
	}
	if g.keyRepeat(ebiten.KeyMinus) || g.keyRepeat(ebiten.KeyKPSubtract) {
		g.lambda -= 10e-9
		if g.lambda < 380e-9 {
			g.lambda = 380e-9
		}
		g.dirty = true
	}

	// Distance L
	if g.keyRepeat(ebiten.KeyL) {
		if shift {
			g.L = math.Max(0.05, g.L-0.05)
		} else {
			g.L = math.Min(5.0, g.L+0.05)
		}
		g.dirty = true
	}

	// Aperture size (W, H, D)
	if g.keyRepeat(ebiten.KeyW) {
		step := 5e-6
		if shift {
			g.apW = math.Max(2e-6, g.apW-step)
		} else {
			g.apW = math.Min(400e-6, g.apW+step)
		}
		g.dirty = true
	}
	if g.keyRepeat(ebiten.KeyH) {
		step := 5e-6
		if shift {
			g.apH = math.Max(2e-6, g.apH-step)
		} else {
			g.apH = math.Min(400e-6, g.apH+step)
		}
		g.dirty = true
	}
	if g.keyRepeat(ebiten.KeyT) {
		step := 5e-6
		if shift {
			g.apD = math.Max(2e-6, g.apD-step)
		} else {
			g.apD = math.Min(400e-6, g.apD+step)
		}
		g.dirty = true
	}

	// Number of slits in grating
	if inpututil.IsKeyJustPressed(ebiten.KeyBracketLeft) {
		if g.apNSlits > 2 {
			g.apNSlits--
		}
		g.dirty = true
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyBracketRight) {
		if g.apNSlits < 64 {
			g.apNSlits++
		}
		g.dirty = true
	}

	// Grid size
	if inpututil.IsKeyJustPressed(ebiten.KeyN) {
		switch g.N {
		case 256:
			g.N = 512
		case 512:
			g.N = 1024
		case 1024:
			g.N = 256
		}
		g.paintBuf = nil
		g.allocBuffers()
		g.dirty = true
	}

	// Color mode toggle
	if inpututil.IsKeyJustPressed(ebiten.KeyC) {
		g.colorMode = 1 - g.colorMode
		g.dirty = true
	}

	// Display mode (log / gamma / linear)
	if inpututil.IsKeyJustPressed(ebiten.KeyD) {
		g.displayMode = (g.displayMode + 1) % 3
		g.dirty = true
	}
	if g.keyRepeat(ebiten.KeyG) {
		if shift {
			g.gamma = math.Max(1.0, g.gamma-0.1)
		} else {
			g.gamma = math.Min(6.0, g.gamma+0.1)
		}
		g.dirty = true
	}

	// Save PNG
	if inpututil.IsKeyJustPressed(ebiten.KeyS) {
		g.savePNG()
	}

	// Reset
	if inpututil.IsKeyJustPressed(ebiten.KeyR) {
		if g.apType == apPaint {
			for i := range g.paintBuf {
				g.paintBuf[i] = 0
			}
			g.dirty = true
		} else {
			g.lambda = defaultLambda
			g.L = defaultL
			g.apW = 60e-6
			g.apH = 40e-6
			g.apD = 20e-6
			g.apNSlits = 8
			g.dirty = true
		}
	}

	// Mouse painting
	if g.apType == apPaint {
		mx, my := ebiten.CursorPosition()
		leftPress := ebiten.IsMouseButtonPressed(ebiten.MouseButtonLeft)
		rightPress := ebiten.IsMouseButtonPressed(ebiten.MouseButtonRight)
		if leftPress || rightPress {
			val := 1.0
			if rightPress {
				val = 0
			}
			if g.paintAt(mx, my, val) {
				g.dirty = true
			}
		} else {
			g.hasLast = false
		}
	}

	if g.dirty {
		g.recompute()
		g.dirty = false
	}
	return nil
}

// paintAt maps a screen coordinate to an aperture index and stamps a small disk.
// Returns true if anything was modified. Also rasterizes a line from the previous
// position so fast drags do not skip pixels.
func (g *Game) paintAt(mx, my int, val float64) bool {
	// Aperture panel rect on screen
	if mx < panelX1 || mx >= panelX1+panelSize || my < panelY1 || my >= panelY1+panelSize {
		g.hasLast = false
		return false
	}
	// Map to aperture pixel
	fx := float64(mx-panelX1) / float64(panelSize) * float64(g.N)
	fy := float64(my-panelY1) / float64(panelSize) * float64(g.N)
	ix := int(fx)
	iy := int(fy)

	radius := g.N / 64
	if radius < 3 {
		radius = 3
	}

	stamp := func(cx, cy int) {
		for dy := -radius; dy <= radius; dy++ {
			y := cy + dy
			if y < 0 || y >= g.N {
				continue
			}
			for dx := -radius; dx <= radius; dx++ {
				x := cx + dx
				if x < 0 || x >= g.N {
					continue
				}
				if dx*dx+dy*dy <= radius*radius {
					g.paintBuf[y*g.N+x] = val
				}
			}
		}
	}

	if g.hasLast {
		// Rasterize line from last (in aperture coords) to current.
		x0 := g.lastPaint.X
		y0 := g.lastPaint.Y
		dx := ix - x0
		dy := iy - y0
		steps := int(math.Max(math.Abs(float64(dx)), math.Abs(float64(dy))))
		if steps < 1 {
			steps = 1
		}
		for s := 0; s <= steps; s++ {
			t := float64(s) / float64(steps)
			stamp(x0+int(float64(dx)*t), y0+int(float64(dy)*t))
		}
	} else {
		stamp(ix, iy)
	}
	g.lastPaint = image.Point{X: ix, Y: iy}
	g.hasLast = true
	return true
}

// ──────────────────────────────────────────────
// Rendering
// ──────────────────────────────────────────────

func (g *Game) Draw(screen *ebiten.Image) {
	screen.Fill(color.RGBA{14, 14, 18, 255})

	face := basicfont.Face7x13
	white := color.RGBA{220, 220, 220, 255}
	dim := color.RGBA{150, 150, 150, 255}
	gray := color.RGBA{130, 130, 130, 255}

	// Aperture panel
	g.drawPanel(screen, g.apImage, panelX1, panelY1, panelSize, panelSize)
	text.Draw(screen, "t(x, y)  -  amplitude mask", face, panelX1, panelY1-10, white)
	text.Draw(screen, fmt.Sprintf("window: %.0f um", float64(g.N)*g.dx*1e6),
		face, panelX1, panelY1+panelSize+18, dim)

	// Diffraction panel
	g.drawPanel(screen, g.diffImage, panelX2, panelY1, panelSize, panelSize)
	colorLbl := "mono - plasma"
	if g.colorMode == 1 {
		colorLbl = "RGB (650/550/450 nm)"
	}
	text.Draw(screen, "I(X, Y)  -  diffraction pattern  ["+colorLbl+"]",
		face, panelX2, panelY1-10, white)
	// physical window of diffraction panel (based on dominant lambda)
	lambdaShown := g.lambda
	if g.colorMode == 1 {
		lambdaShown = lambdaR
	}
	Wobs := lambdaShown * g.L / g.dx
	text.Draw(screen, fmt.Sprintf("window: %s", formatLength(Wobs)),
		face, panelX2, panelY1+panelSize+18, dim)

	// Crosshair at observation center
	cx := float32(panelX2 + panelSize/2)
	cy := float32(panelY1 + panelSize/2)
	cross := color.RGBA{255, 255, 255, 60}
	vector.StrokeLine(screen, cx-6, cy, cx+6, cy, 1, cross, false)
	vector.StrokeLine(screen, cx, cy-6, cx, cy+6, 1, cross, false)

	// Status / parameters block
	g.drawStatus(screen, face, white, dim, gray)

	// Colormap legend (only in mono mode)
	if g.colorMode == 0 {
		g.drawLegend(screen, face, white, dim, gray)
	}
}

func (g *Game) drawPanel(screen, img *ebiten.Image, x, y, w, h int) {
	// Border
	vector.StrokeRect(screen, float32(x-1), float32(y-1),
		float32(w+2), float32(h+2), 1, color.RGBA{90, 90, 100, 255}, false)
	if img == nil {
		return
	}
	bw := img.Bounds().Dx()
	bh := img.Bounds().Dy()
	op := &ebiten.DrawImageOptions{}
	op.GeoM.Scale(float64(w)/float64(bw), float64(h)/float64(bh))
	op.GeoM.Translate(float64(x), float64(y))
	op.Filter = ebiten.FilterLinear
	screen.DrawImage(img, op)
}

func (g *Game) drawStatus(screen *ebiten.Image, face *basicfont.Face,
	white, dim, gray color.RGBA) {

	y := panelY1 + panelSize + 38
	x := panelX1

	apNames := []string{"rectangle", "circle", "double slit", "grating", "cross", "ring", "paint mode"}
	apName := apNames[g.apType]

	// Fresnel number - half-width of the largest feature
	a := math.Max(g.apW/2, g.apH/2)
	NF := a * a / (g.lambda * g.L)

	dispNames := []string{"log", "gamma", "linear"}
	colorName := "mono+plasma"
	if g.colorMode == 1 {
		colorName = "RGB"
	}

	text.Draw(screen, fmt.Sprintf("wl = %.0f nm   L = %.2f m   N = %d   N_F = %.3f",
		g.lambda*1e9, g.L, g.N, NF), face, x, y, white)
	text.Draw(screen, fmt.Sprintf("aperture: %s    W=%s H=%s T=%s   slits=%d",
		apName,
		formatLength(g.apW), formatLength(g.apH), formatLength(g.apD), g.apNSlits),
		face, x, y+18, white)
	text.Draw(screen, fmt.Sprintf("mode: %s - %s - gamma=%.1f    frame: %.1f ms",
		colorName, dispNames[g.displayMode], g.gamma, g.lastRecomputeMs),
		face, x, y+36, dim)

	if NF > 0.1 {
		text.Draw(screen, "! N_F > 0.1: far field broken, pattern is approximate",
			face, x, y+54, color.RGBA{240, 180, 80, 255})
	} else {
		text.Draw(screen, "OK  N_F << 1: Fraunhofer regime",
			face, x, y+54, color.RGBA{120, 200, 130, 255})
	}

	// Controls
	hint := "1-6 apertures | M paint | W/H/T sizes | [/] slits | +/- wavelength | L/Shift+L distance | " +
		"C mono<>RGB | D log/gamma/lin | G gamma | N grid | S save | R reset"
	text.Draw(screen, hint, face, x, screenHeight-10, gray)
}

func (g *Game) drawLegend(screen *ebiten.Image, face *basicfont.Face,
	white, dim, gray color.RGBA) {
	legendX := panelX2 + panelSize - 180
	legendY := panelY1 - 28
	legendW := 180
	legendH := 10
	for px := 0; px < legendW; px++ {
		t := float64(px) / float64(legendW-1)
		r, gg, b := plasmaColor(t)
		vector.StrokeLine(screen,
			float32(legendX+px), float32(legendY),
			float32(legendX+px), float32(legendY+legendH),
			1, color.RGBA{r, gg, b, 255}, false)
	}
	vector.StrokeRect(screen, float32(legendX), float32(legendY),
		float32(legendW), float32(legendH), 1, color.RGBA{220, 220, 220, 180}, false)
	text.Draw(screen, "0", face, legendX-10, legendY+legendH, gray)
	text.Draw(screen, "I_max", face, legendX+legendW+4, legendY+legendH, gray)
}

func (g *Game) Layout(_, _ int) (int, int) {
	return screenWidth, screenHeight
}

// ──────────────────────────────────────────────
// PNG export
// ──────────────────────────────────────────────

func (g *Game) savePNG() {
	if g.diffImage == nil {
		return
	}
	N := g.N
	img := image.NewRGBA(image.Rect(0, 0, N, N))
	pix := make([]byte, N*N*4)
	g.diffImage.ReadPixels(pix)
	copy(img.Pix, pix)

	ts := time.Now().Format("20060102-150405")
	name := fmt.Sprintf("diffraction-%s.png", ts)
	f, err := os.Create(name)
	if err != nil {
		return
	}
	defer f.Close()
	_ = png.Encode(f, img)
}

// ──────────────────────────────────────────────
// Helpers
// ──────────────────────────────────────────────

func formatLength(m float64) string {
	switch {
	case m >= 1e-3:
		return fmt.Sprintf("%.2f mm", m*1e3)
	case m >= 1e-6:
		return fmt.Sprintf("%.1f um", m*1e6)
	default:
		return fmt.Sprintf("%.0f nm", m*1e9)
	}
}

// kept for completeness, lets future code rely on complex math without warning
var _ = cmplx.Abs

func main() {
	ebiten.SetWindowSize(screenWidth, screenHeight)
	ebiten.SetWindowTitle("Дифракция Фраунгофера")
	ebiten.SetWindowResizingMode(ebiten.WindowResizingModeEnabled)
	if err := ebiten.RunGame(newGame()); err != nil {
		panic(err)
	}
}
