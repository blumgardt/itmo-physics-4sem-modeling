package main

import (
	"fmt"
	"image/color"
	"math"
	"sync"

	"github.com/hajimehoshi/ebiten/v2"
	"github.com/hajimehoshi/ebiten/v2/inpututil"
	"github.com/hajimehoshi/ebiten/v2/text"
	"github.com/hajimehoshi/ebiten/v2/vector"
	"golang.org/x/image/font/basicfont"
)

const (
	screenWidth  = 1200
	screenHeight = 800
	mu0          = 4 * math.Pi * 1e-7

	gridNZ = 161 // physics grid resolution (z axis)
	gridNY = 101 // physics grid resolution (y axis)
)

// ──────────────────────────────────────────────
// Elliptic integrals (AGM method)
// ──────────────────────────────────────────────

func ellipticK(m float64) float64 {
	if m < 0 {
		m = 0
	}
	if m >= 1 {
		m = 1 - 1e-15
	}
	a, g := 1.0, math.Sqrt(1-m)
	for math.Abs(a-g) > 1e-15 {
		a, g = (a+g)/2, math.Sqrt(a*g)
	}
	return math.Pi / (2 * a)
}

func ellipticE(m float64) float64 {
	if m < 0 {
		m = 0
	}
	if m >= 1 {
		return 1.0
	}
	if m == 0 {
		return math.Pi / 2
	}
	a, g := 1.0, math.Sqrt(1-m)
	sum, pow2 := m, 1.0
	for {
		aNew := (a + g) / 2
		gNew := math.Sqrt(a * g)
		c := (a - g) / 2
		pow2 *= 2
		sum += pow2 * c * c
		a, g = aNew, gNew
		if math.Abs(c) < 1e-15 {
			break
		}
	}
	return math.Pi / (2 * a) * (1 - sum/2)
}

// ──────────────────────────────────────────────
// Biot-Savart: single current loop
// ──────────────────────────────────────────────

func singleLoopField(R, current float64, r, z []float64) (Br, Bz []float64) {
    n := len(r)
    Br = make([]float64, n)
    Bz = make([]float64, n)
    for i := 0; i < n; i++ {
        ri, zi := r[i], z[i]
        if ri < 1e-12 {
            Br[i] = 0
            denom := math.Pow(R*R+zi*zi, 1.5)
            if denom > 0 {
                Bz[i] = mu0 * current * R * R / (2 * denom)
            }
            continue
        }
        alpha2 := (R+ri)*(R+ri) + zi*zi
        beta2 := (R-ri)*(R-ri) + zi*zi
        k2 := 4 * R * ri / alpha2
        if k2 < 0 {
            k2 = 0
        }
        if k2 > 1-1e-12 {
            k2 = 1 - 1e-12
        }
        K := ellipticK(k2)
        E := ellipticE(k2)
        alpha := math.Sqrt(alpha2)
        c := mu0 * current / (2 * math.Pi)
        Bz[i] = c / alpha * (K + (R*R-ri*ri-zi*zi)/beta2*E)
        Br[i] = c * zi / (ri * alpha) * (-K + (R*R+ri*ri+zi*zi)/beta2*E)
    }
    return
}

// ──────────────────────────────────────────────
// Solenoid field (superposition of N loops)
// ──────────────────────────────────────────────

func turnCenters(L float64, N int) []float64 {
	if N < 1 {
		return nil
	}
	centers := make([]float64, N)
	pitch := L / float64(N)
	for i := range centers {
		centers[i] = -L/2 + pitch*(float64(i)+0.5)
	}
	return centers
}

func solenoidField(D, L float64, N int, current float64, rGrid, zGrid []float64) (Br, Bz []float64) {
	R := D / 2
	n := len(rGrid)
	Br = make([]float64, n)
	Bz = make([]float64, n)
	centers := turnCenters(L, N)
	if len(centers) == 0 {
		return
	}

	nWorkers := 8
	if len(centers) < nWorkers {
		nWorkers = len(centers)
	}

	var mu sync.Mutex
	var wg sync.WaitGroup
	coilsPerWorker := (len(centers) + nWorkers - 1) / nWorkers

	for w := 0; w < nWorkers; w++ {
		start := w * coilsPerWorker
		end := start + coilsPerWorker
		if end > len(centers) {
			end = len(centers)
		}
		wg.Add(1)
		go func(s, e int) {
			defer wg.Done()
			localBr := make([]float64, n)
			localBz := make([]float64, n)
			zShifted := make([]float64, n)
			for coil := s; coil < e; coil++ {
				zc := centers[coil]
				for j := range zShifted {
					zShifted[j] = zGrid[j] - zc
				}
				dBr, dBz := singleLoopField(R, current, rGrid, zShifted)
				for j := 0; j < n; j++ {
					localBr[j] += dBr[j]
					localBz[j] += dBz[j]
				}
			}
			mu.Lock()
			for j := 0; j < n; j++ {
				Br[j] += localBr[j]
				Bz[j] += localBz[j]
			}
			mu.Unlock()
		}(start, end)
	}
	wg.Wait()
	return
}

// ──────────────────────────────────────────────
// Streamline tracing (RK2)
// ──────────────────────────────────────────────

type vec2 struct{ x, y float64 }

func computeStreamlines(zLin, yLin []float64, Bz, By [][]float64, L, R, scale float64, nz, ny, N int) [][]vec2 {
	dz := zLin[1] - zLin[0]
	dy := yLin[1] - yLin[0]
	dt := math.Min(dz, dy) * 0.25

	centers := turnCenters(L, N)
	pitch := L
	if N > 0 {
		pitch = L / float64(N)
	}
	// Do not trace too close to the discrete wire positions: the analytical loop
	// field is singular there, and the interpolated streamline integrator produces
	// visible numerical curls instead of a clean macroscopic solenoid pattern.
	wireGuard := math.Max(dt*1.5, pitch*0.75)
	blocked := func(z, y float64) bool {
		if len(centers) == 0 {
			return false
		}
		if math.Abs(math.Abs(y)-R) > wireGuard {
			return false
		}
		for _, zc := range centers {
			if math.Hypot(z-zc, math.Abs(y)-R) < wireGuard {
				return true
			}
		}
		return false
	}

	interp := func(z, y float64) (bz, by float64, ok bool) {
		fi := (z - zLin[0]) / dz
		fj := (y - yLin[0]) / dy
		i, j := int(fi), int(fj)
		if i < 0 || i >= nz-1 || j < 0 || j >= ny-1 {
			return 0, 0, false
		}
		fx, fy := fi-float64(i), fj-float64(j)
		bz = Bz[j][i]*(1-fx)*(1-fy) + Bz[j][i+1]*fx*(1-fy) +
			Bz[j+1][i]*(1-fx)*fy + Bz[j+1][i+1]*fx*fy
		by = By[j][i]*(1-fx)*(1-fy) + By[j][i+1]*fx*(1-fy) +
			By[j+1][i]*(1-fx)*fy + By[j+1][i+1]*fx*fy
		return bz, by, true
	}

	// trace one direction: dir=+1 forward, dir=-1 backward
	traceHalf := func(seed vec2, dir float64, maxSteps int) []vec2 {
		if blocked(seed.x, seed.y) {
			return nil
		}
		pts := []vec2{seed}
		z, y := seed.x, seed.y
		for step := 0; step < maxSteps; step++ {
			if blocked(z, y) {
				break
			}
			bz1, by1, ok := interp(z, y)
			if !ok {
				break
			}
			mag := math.Hypot(bz1, by1)
			if mag < 1e-20 {
				break
			}
			dzStep := dir * bz1 / mag * dt
			dyStep := dir * by1 / mag * dt
			bz2, by2, ok := interp(z+dzStep*0.5, y+dyStep*0.5)
			if !ok {
				break
			}
			mag2 := math.Hypot(bz2, by2)
			if mag2 < 1e-20 {
				break
			}
			nextZ := z + dir*bz2/mag2*dt
			nextY := y + dir*by2/mag2*dt
			if blocked(nextZ, nextY) {
				break
			}
			z, y = nextZ, nextY
			pts = append(pts, vec2{z, y})
		}
		return pts
	}

	var seeds []vec2

	// Внутри соленоида: несколько радиусов от левого торца
	for _, f := range []float64{-0.9, -0.7, -0.5, -0.3, -0.1, 0.0, 0.1, 0.3, 0.5, 0.7, 0.9} {
		seeds = append(seeds, vec2{-L / 2, f * R * 0.95})
	}

	// Вблизи оси: больше точек по z
	nAxis := 15
	if scale > 3 {
		nAxis = 21
	}
	for k := 0; k < nAxis; k++ {
		f := -1.2 + 2.4*float64(k)/float64(nAxis-1)
		zs := f * L * math.Min(scale*0.9, 3.5)
		seeds = append(seeds, vec2{zs, 0.01 * R}, vec2{zs, -0.01 * R})
	}

	// Снаружи торцов: дальше и с большим радиусом
	for _, f := range []float64{-0.8, -0.5, -0.2, 0.2, 0.5, 0.8} {
		seeds = append(seeds, vec2{-L/2 - 0.2*L, f * R * 1.8})
		seeds = append(seeds, vec2{L/2 + 0.2*L, f * R * 1.8})
	}

	// На поверхности соленоида (вдоль боковой стенки)
	for _, f := range []float64{-0.45, -0.25, 0.0, 0.25, 0.45} {
		seeds = append(seeds, vec2{f * L, R * 1.12}, vec2{f * L, -R * 1.12})
	}

	maxSteps := int(900 * math.Min(scale, 4))
	var lines [][]vec2

	for _, s := range seeds {
		fwd := traceHalf(s, +1, maxSteps)
		bwd := traceHalf(s, -1, maxSteps)

		// combine: reversed backward + forward (skip duplicate seed)
		line := make([]vec2, 0, len(bwd)+len(fwd))
		for i := len(bwd) - 1; i >= 0; i-- {
			line = append(line, bwd[i])
		}
		if len(fwd) > 1 {
			line = append(line, fwd[1:]...)
		}

		if len(line) > 5 {
			lines = append(lines, line)
		}
	}
	return lines
}

// ──────────────────────────────────────────────
// Colormap (plasma-like)
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
// Game
// ──────────────────────────────────────────────

type Game struct {
	D, L, I float64
	N       int
	scale   float64
	dirty   bool

	bgImage     *ebiten.Image
	streamlines [][]vec2
	bTheory     float64
	bCenter     float64
	heatmapMax  float64
	heatmapAuto bool

	extZ, extY float64 // world extents for coordinate mapping

	holdFrames map[ebiten.Key]int // frames each key has been held
}

func newGame() *Game {
	return &Game{
		D:          0.10, // 10 cm
		L:          0.30, // 30 cm
			N:          100,
			I:          1.0,
			scale:      1.5,
			dirty:      true,
			heatmapAuto: true,
			holdFrames: make(map[ebiten.Key]int),
		}
}

// keyRepeat returns true on first press, then after a delay repeats every few frames.
// Initial delay: 20 frames (~333ms), repeat every 4 frames (~67ms).
func (g *Game) keyRepeat(key ebiten.Key) bool {
	if ebiten.IsKeyPressed(key) {
		g.holdFrames[key]++
		f := g.holdFrames[key]
		return f == 1 || (f > 20 && f%4 == 0)
	}
	g.holdFrames[key] = 0
	return false
}

func (g *Game) worldToScreen(wz, wy float64) (float32, float32) {
	sx := float64(screenWidth)/2 + wz/g.extZ*float64(screenWidth)/2
	sy := float64(screenHeight)/2 - wy/g.extY*float64(screenHeight)/2
	return float32(sx), float32(sy)
}

func (g *Game) recompute() {
	// extents with correct aspect ratio
	base := g.scale * math.Max(g.L, g.D)
	g.extY = base
	g.extZ = base * float64(screenWidth) / float64(screenHeight)

	nz, ny := gridNZ, gridNY
	zLin := linspace(-g.extZ, g.extZ, nz)
	yLin := linspace(-g.extY, g.extY, ny)

	// flatten grid for solenoid computation
	size := nz * ny
	rFlat := make([]float64, size)
	zFlat := make([]float64, size)
	for j := 0; j < ny; j++ {
		for i := 0; i < nz; i++ {
			idx := j*nz + i
			rFlat[idx] = math.Abs(yLin[j])
			zFlat[idx] = zLin[i]
		}
	}

	BrFlat, BzFlat := solenoidField(g.D, g.L, g.N, g.I, rFlat, zFlat)

	// build 2D field arrays
	R := g.D / 2
	Bz2D := make([][]float64, ny)
	By2D := make([][]float64, ny)
	bMag2D := make([][]float64, ny)

	for j := 0; j < ny; j++ {
		Bz2D[j] = make([]float64, nz)
		By2D[j] = make([]float64, nz)
		bMag2D[j] = make([]float64, nz)
		for i := 0; i < nz; i++ {
			idx := j*nz + i
			sign := 1.0
			if yLin[j] < 0 {
				sign = -1.0
			}
			if math.Abs(yLin[j]) < 1e-15 {
				sign = 0
			}
			by := BrFlat[idx] * sign
			bz := BzFlat[idx]
			Bz2D[j][i] = bz
			By2D[j][i] = by
			bMag2D[j][i] = math.Hypot(bz, by)
		}
	}

	g.bTheory = mu0 * float64(g.N) / g.L * g.I
	cBr, cBz := solenoidField(g.D, g.L, g.N, g.I, []float64{0}, []float64{0})
	g.bCenter = math.Hypot(cBr[0], cBz[0])

	// render heatmap to a small image (nz × ny), Ebiten scales it up
	currentMax := 0.0
	for j := 0; j < ny; j++ {
		for i := 0; i < nz; i++ {
			if bMag2D[j][i] > currentMax {
				currentMax = bMag2D[j][i]
			}
		}
	}
	if g.heatmapAuto {
		g.heatmapMax = 1.05 * currentMax
	} else if g.heatmapMax < 1e-15 {
		g.heatmapMax = 1.5 * g.bTheory
	}
	if g.heatmapMax < 1e-15 {
		g.heatmapMax = 1e-15
	}

	pix := make([]byte, nz*ny*4)
	for j := 0; j < ny; j++ {
		screenJ := ny - 1 - j // flip y: grid j=0 is bottom, pixel row 0 is top
		for i := 0; i < nz; i++ {
			t := bMag2D[j][i] / g.heatmapMax
			r, g, b := plasmaColor(t)
			idx := (screenJ*nz + i) * 4
			pix[idx] = r
			pix[idx+1] = g
			pix[idx+2] = b
			pix[idx+3] = 255
		}
	}
	img := ebiten.NewImage(nz, ny)
	img.WritePixels(pix)
	g.bgImage = img

	// streamlines
	g.streamlines = computeStreamlines(zLin, yLin, Bz2D, By2D, g.L, R, g.scale, nz, ny, g.N)
}

// ──────────────────────────────────────────────
// Input handling
// ──────────────────────────────────────────────

func (g *Game) Update() error {
	shift := ebiten.IsKeyPressed(ebiten.KeyShift)

	if g.keyRepeat(ebiten.KeyD) {
		if shift {
			g.D = math.Max(0.02, g.D-0.01)
		} else {
			g.D += 0.01
		}
		g.dirty = true
	}
	if g.keyRepeat(ebiten.KeyL) {
		if shift {
			g.L = math.Max(0.05, g.L-0.05)
		} else {
			g.L += 0.05
		}
		g.dirty = true
	}
	if g.keyRepeat(ebiten.KeyN) {
    if shift {
        if g.N > 10 {
            g.N -= 10
        } else {
            g.N = 1
        }
    } else {
        g.N += 10
        if g.N > 2000 {
            g.N = 2000
        }
    }
    g.dirty = true
}
	if g.keyRepeat(ebiten.KeyI) {
		if shift {
			g.I = math.Max(0.1, g.I-0.1)
		} else {
			g.I += 0.1
		}
		g.dirty = true
	}
	if g.keyRepeat(ebiten.KeyEqual) || g.keyRepeat(ebiten.KeyKPAdd) {
		g.scale = math.Min(6, g.scale+0.5)
		g.dirty = true
	}
	if g.keyRepeat(ebiten.KeyMinus) || g.keyRepeat(ebiten.KeyKPSubtract) {
		g.scale = math.Max(1, g.scale-0.5)
		g.dirty = true
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyR) {
		g.D, g.L, g.N, g.I, g.scale = 0.10, 0.30, 100, 1.0, 1.5
		g.heatmapAuto = true
		g.heatmapMax = 0
		g.dirty = true
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyH) {
		if g.heatmapAuto {
			g.heatmapAuto = false
			g.heatmapMax = 1.5 * g.bTheory
		} else {
			g.heatmapAuto = true
		}
		g.dirty = true
	}

	if g.dirty {
		g.recompute()
		g.dirty = false
	}
	return nil
}

// ──────────────────────────────────────────────
// Rendering
// ──────────────────────────────────────────────

func (g *Game) Draw(screen *ebiten.Image) {
	// background heatmap (scaled up from low-res image)
	if g.bgImage != nil {
		op := &ebiten.DrawImageOptions{}
		bw, bh := g.bgImage.Bounds().Dx(), g.bgImage.Bounds().Dy()
		op.GeoM.Scale(float64(screenWidth)/float64(bw), float64(screenHeight)/float64(bh))
		op.Filter = ebiten.FilterLinear
		screen.DrawImage(g.bgImage, op)
	}

	// streamlines
	lineClr := color.RGBA{255, 255, 255, 150}
	for _, line := range g.streamlines {
		for k := 0; k < len(line)-1; k++ {
			x1, y1 := g.worldToScreen(line[k].x, line[k].y)
			x2, y2 := g.worldToScreen(line[k+1].x, line[k+1].y)
			vector.StrokeLine(screen, x1, y1, x2, y2, 1.2, lineClr, true)
		}
	}

	// solenoid outline
	halfL := g.L / 2
	R := g.D / 2
	tlx, tly := g.worldToScreen(-halfL, R)
	trx, try_ := g.worldToScreen(halfL, R)
	blx, bly := g.worldToScreen(-halfL, -R)
	brx, bry := g.worldToScreen(halfL, -R)
	cyan := color.RGBA{0, 229, 255, 200}
	vector.StrokeLine(screen, tlx, tly, trx, try_, 2.5, cyan, true)
	vector.StrokeLine(screen, blx, bly, brx, bry, 2.5, cyan, true)
	vector.StrokeLine(screen, tlx, tly, blx, bly, 2, cyan, true)
	vector.StrokeLine(screen, trx, try_, brx, bry, 2, cyan, true)

	// Рисуем несколько линий витков
	nLines := 9
	if g.N > 0 {
		nLines = min(9, g.N)
	}
	pitch := g.L / float64(nLines)
	for i := 0; i < nLines; i++ {
		zPos := -g.L/2 + pitch*(float64(i)+0.5)
		x1, y1 := g.worldToScreen(zPos, R)
		x2, y2 := g.worldToScreen(zPos, -R)
		lineColor := color.RGBA{0, 200, 200, 80}
		vector.StrokeLine(screen, x1, y1, x2, y2, 1.0, lineColor, true)
	}

	// axis labels (z and y)
	face := basicfont.Face7x13
	gray := color.RGBA{130, 130, 130, 255}
	text.Draw(screen, "z", face, screenWidth-20, screenHeight/2+16, gray)
	text.Draw(screen, "y", face, screenWidth/2+6, 16, gray)

	// parameter info (top-left)
	white := color.RGBA{220, 220, 220, 255}
	modeLabel := "Auto"
	if !g.heatmapAuto {
		modeLabel = "Fixed"
	}
	text.Draw(screen, fmt.Sprintf("D = %.1f cm   [D / Shift+D]", g.D*100), face, 10, 20, white)
	text.Draw(screen, fmt.Sprintf("L = %.0f cm    [L / Shift+L]", g.L*100), face, 10, 36, white)
	text.Draw(screen, fmt.Sprintf("N = %d         [N / Shift+N]", g.N), face, 10, 52, white)
	text.Draw(screen, fmt.Sprintf("I = %.1f A     [I / Shift+I]", g.I), face, 10, 68, white)
	text.Draw(screen, fmt.Sprintf("Scale %.1fx     [+/-]", g.scale), face, 10, 84, white)
	text.Draw(screen, "[R] reset", face, 10, 100, color.RGBA{150, 150, 150, 255})
	text.Draw(screen, fmt.Sprintf("[H] heatmap scale: %s", modeLabel), face, 10, 116, color.RGBA{180, 180, 180, 255})

	// field info (bottom-left)
	yellow := color.RGBA{240, 200, 80, 255}
	orange := color.RGBA{200, 160, 100, 255}
	text.Draw(screen, fmt.Sprintf("B(center) = %s", formatField(g.bCenter)), face, 10, screenHeight-30, yellow)
	text.Draw(screen, fmt.Sprintf("B(theory, inf) = %s", formatField(g.bTheory)), face, 10, screenHeight-14, orange)

	legendX, legendY := screenWidth-220, 18
	legendW, legendH := 180, 14
	for x := 0; x < legendW; x++ {
		t := float64(x) / float64(legendW-1)
		r, g, b := plasmaColor(t)
		vector.StrokeLine(
			screen,
			float32(legendX+x), float32(legendY),
			float32(legendX+x), float32(legendY+legendH),
			1,
			color.RGBA{r, g, b, 255},
			false,
		)
	}
	vector.StrokeRect(screen, float32(legendX), float32(legendY), float32(legendW), float32(legendH), 1, color.RGBA{220, 220, 220, 180}, false)
	text.Draw(screen, "0", face, legendX, legendY+legendH+12, gray)
	text.Draw(screen, formatField(g.heatmapMax), face, legendX+legendW-70, legendY+legendH+12, gray)
	text.Draw(screen, "|B|", face, legendX+legendW/2-8, legendY-4, white)
	text.Draw(screen, modeLabel, face, legendX, legendY-4, gray)
}

func (g *Game) Layout(_, _ int) (int, int) {
	return screenWidth, screenHeight
}

// ──────────────────────────────────────────────
// Helpers
// ──────────────────────────────────────────────

func formatField(val float64) string {
	uT := val * 1e6
	switch {
	case uT >= 1000:
		return fmt.Sprintf("%.2f mT", uT/1000)
	case uT >= 1:
		return fmt.Sprintf("%.1f uT", uT)
	default:
		return fmt.Sprintf("%.1f nT", uT*1000)
	}
}

func linspace(start, end float64, n int) []float64 {
	s := make([]float64, n)
	if n == 1 {
		s[0] = start
		return s
	}
	step := (end - start) / float64(n-1)
	for i := range s {
		s[i] = start + float64(i)*step
	}
	return s
}

// ──────────────────────────────────────────────
// Entry point
// ──────────────────────────────────────────────

func main() {
	ebiten.SetWindowSize(screenWidth, screenHeight)
	ebiten.SetWindowTitle("Магнитное поле соленоида")
	ebiten.SetWindowResizingMode(ebiten.WindowResizingModeEnabled)
	if err := ebiten.RunGame(newGame()); err != nil {
		panic(err)
	}
}
