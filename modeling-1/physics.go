package main

import "math"

const mu0 = 4 * math.Pi * 1e-7 // magnetic permeability of free space, T·m/A

// Params holds solenoid parameters.
type Params struct {
	D     float64 `json:"d"`     // diameter, m
	L     float64 `json:"l"`     // length, m
	N     int     `json:"n"`     // number of turns
	I     float64 `json:"i"`     // current, A
	Scale float64 `json:"scale"` // view scale multiplier (default 1.5)
}

// Point is an (X, Y) coordinate for streamlines.
type Point struct {
	X float64 `json:"x"`
	Y float64 `json:"y"`
}

// FieldResult contains all computed data for visualization.
type FieldResult struct {
	ZAxis       []float64   `json:"z_axis"`      // z coordinates (cm)
	YAxis       []float64   `json:"y_axis"`      // y coordinates (cm)
	BMag        [][]float64 `json:"b_mag"`       // |B| magnitude (T), [ny][nz]
	BTheory     float64     `json:"b_theory"`    // theoretical B for infinite solenoid (T)
	BCenter     float64     `json:"b_center"`    // B at center (T)
	Streamlines [][]Point   `json:"streamlines"` // field lines
	SolR        float64     `json:"sol_r"`       // solenoid radius (cm)
	SolHalfL    float64     `json:"sol_half_l"`  // solenoid half-length (cm)
}

func turnCenters(L float64, N int) []float64 {
	if N < 1 {
		return nil
	}

	centers := make([]float64, N)
	pitch := L / float64(N)
	for coil := range centers {
		centers[coil] = -L/2 + pitch*(float64(coil)+0.5)
	}

	return centers
}

// singleLoopField computes (Br, Bz) for a single current loop at the origin.
func singleLoopField(R, current float64, r, z []float64) (Br, Bz []float64) {
	n := len(r)
	Br = make([]float64, n)
	Bz = make([]float64, n)

	for i := 0; i < n; i++ {
		ri := r[i]
		zi := z[i]

		if math.Abs(ri) < 1e-10 {
			// on-axis
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

		K := EllipticK(k2)
		E := EllipticE(k2)
		alpha := math.Sqrt(alpha2)
		c := mu0 * current / (2 * math.Pi)

		Bz[i] = c / alpha * (K + (R*R-ri*ri-zi*zi)/beta2*E)
		Br[i] = c * zi / (ri * alpha) * (-K + (R*R+ri*ri+zi*zi)/beta2*E)
	}

	return
}

// solenoidField computes total (Br, Bz) from N loops distributed along the axis.
// Uses goroutines to parallelize across coils.
func solenoidField(D, L float64, N int, current float64, rGrid, zGrid []float64) (Br, Bz []float64) {
	R := D / 2
	n := len(rGrid)
	Br = make([]float64, n)
	Bz = make([]float64, n)
	centers := turnCenters(L, N)
	if len(centers) == 0 {
		return Br, Bz
	}

	// Determine number of worker goroutines
	nWorkers := 8
	if len(centers) < nWorkers {
		nWorkers = len(centers)
	}

	type partialResult struct {
		br, bz []float64
	}

	ch := make(chan partialResult, nWorkers)

	coilsPerWorker := (len(centers) + nWorkers - 1) / nWorkers

	for w := 0; w < nWorkers; w++ {
		startCoil := w * coilsPerWorker
		endCoil := startCoil + coilsPerWorker
		if endCoil > len(centers) {
			endCoil = len(centers)
		}
		go func(start, end int) {
			localBr := make([]float64, n)
			localBz := make([]float64, n)
			zShifted := make([]float64, n)
			for coil := start; coil < end; coil++ {
				zCoil := centers[coil]
				for j := range zShifted {
					zShifted[j] = zGrid[j] - zCoil
				}
				dBr, dBz := singleLoopField(R, current, rGrid, zShifted)
				for j := 0; j < n; j++ {
					localBr[j] += dBr[j]
					localBz[j] += dBz[j]
				}
			}
			ch <- partialResult{localBr, localBz}
		}(startCoil, endCoil)
	}

	for w := 0; w < nWorkers; w++ {
		res := <-ch
		for j := 0; j < n; j++ {
			Br[j] += res.br[j]
			Bz[j] += res.bz[j]
		}
	}

	return
}

// ComputeField runs the full computation and returns visualization data.
func ComputeField(p Params) *FieldResult {
	if p.N < 1 {
		p.N = 1
	}
	if p.Scale < 1.0 {
		p.Scale = 1.5
	}

	R := p.D / 2
	nz, ny := 121, 81
	extent := p.Scale * math.Max(p.L, p.D)

	zLin := linspace(-extent, extent, nz)
	yLin := linspace(-extent, extent, ny)

	// flatten grid
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

	BrFlat, BzFlat := solenoidField(p.D, p.L, p.N, p.I, rFlat, zFlat)

	// convert to Cartesian and compute magnitude
	ByFlat := make([]float64, size)
	bMag2D := make([][]float64, ny)
	Bz2D := make([][]float64, ny)
	By2D := make([][]float64, ny)

	for j := 0; j < ny; j++ {
		bMag2D[j] = make([]float64, nz)
		Bz2D[j] = make([]float64, nz)
		By2D[j] = make([]float64, nz)
		for i := 0; i < nz; i++ {
			idx := j*nz + i
			sign := 1.0
			if yLin[j] < 0 {
				sign = -1.0
			}
			if math.Abs(yLin[j]) < 1e-15 {
				sign = 0
			}
			ByFlat[idx] = BrFlat[idx] * sign
			bMag2D[j][i] = math.Sqrt(BzFlat[idx]*BzFlat[idx] + ByFlat[idx]*ByFlat[idx])
			Bz2D[j][i] = BzFlat[idx]
			By2D[j][i] = ByFlat[idx]
		}
	}

	bTheory := mu0 * float64(p.N) / p.L * p.I
	centerBr, centerBz := solenoidField(p.D, p.L, p.N, p.I, []float64{0}, []float64{0})
	bCenter := math.Hypot(centerBr[0], centerBz[0])

	// convert axes to cm
	zCm := make([]float64, nz)
	yCm := make([]float64, ny)
	for i := range zLin {
		zCm[i] = zLin[i] * 100
	}
	for j := range yLin {
		yCm[j] = yLin[j] * 100
	}

	// compute streamlines
	streamlines := computeStreamlines(zLin, yLin, Bz2D, By2D, p.L, R, p.Scale, nz, ny)

	return &FieldResult{
		ZAxis:       zCm,
		YAxis:       yCm,
		BMag:        bMag2D,
		BTheory:     bTheory,
		BCenter:     bCenter,
		Streamlines: streamlines,
		SolR:        R * 100,
		SolHalfL:    p.L / 2 * 100,
	}
}

// computeStreamlines traces field lines through the vector field using RK2 integration.
func computeStreamlines(zLin, yLin []float64, Bz, By [][]float64, L, R, scale float64, nz, ny int) [][]Point {
	dz := zLin[1] - zLin[0]
	dy := yLin[1] - yLin[0]
	dt := math.Min(dz, dy) * 0.5

	// bilinear interpolation of field at (z, y)
	interp := func(z, y float64) (bz, by float64, ok bool) {
		fi := (z - zLin[0]) / dz
		fj := (y - yLin[0]) / dy
		i := int(fi)
		j := int(fj)
		if i < 0 || i >= nz-1 || j < 0 || j >= ny-1 {
			return 0, 0, false
		}
		fx := fi - float64(i)
		fy := fj - float64(j)
		bz = Bz[j][i]*(1-fx)*(1-fy) + Bz[j][i+1]*fx*(1-fy) +
			Bz[j+1][i]*(1-fx)*fy + Bz[j+1][i+1]*fx*fy
		by = By[j][i]*(1-fx)*(1-fy) + By[j][i+1]*fx*(1-fy) +
			By[j+1][i]*(1-fx)*fy + By[j+1][i+1]*fx*fy
		return bz, by, true
	}

	var lines [][]Point

	seeds := make([]Point, 0, 60)

	// seeds along left boundary — scale-adaptive
	seedZ := -(scale - 0.2) * L
	for _, frac := range []float64{-0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9} {
		seeds = append(seeds, Point{seedZ, frac * (scale - 0.1) * (R + L*0.1)})
	}

	// seeds inside solenoid
	for _, frac := range []float64{-0.8, -0.5, -0.2, 0.0, 0.2, 0.5, 0.8} {
		seeds = append(seeds, Point{-L / 2, frac * R * 0.9})
	}

	// seeds along axis — spread with scale
	nAxisSeeds := 7
	if scale > 3 {
		nAxisSeeds = 11
	}
	for k := 0; k < nAxisSeeds; k++ {
		frac := -1.2 + 2.4*float64(k)/float64(nAxisSeeds-1)
		zSeed := frac * L * math.Min(scale*0.8, 3.0)
		seeds = append(seeds, Point{zSeed, 0.01 * R})
		seeds = append(seeds, Point{zSeed, -0.01 * R})
	}

	// extra seeds outside solenoid ends for return field lines at large scale
	if scale > 2 {
		for _, frac := range []float64{-0.6, -0.3, 0.3, 0.6} {
			seeds = append(seeds, Point{-L/2 - 0.1*L, frac * R * 1.5})
			seeds = append(seeds, Point{L/2 + 0.1*L, frac * R * 1.5})
		}
	}

	maxSteps := int(500 * math.Min(scale, 4))

	for _, seed := range seeds {
		line := make([]Point, 0, maxSteps)
		z, y := seed.X, seed.Y
		line = append(line, Point{z * 100, y * 100})

		for step := 0; step < maxSteps; step++ {
			bz1, by1, ok := interp(z, y)
			if !ok {
				break
			}
			mag := math.Sqrt(bz1*bz1 + by1*by1)
			if mag < 1e-20 {
				break
			}
			// normalize and step
			dz := bz1 / mag * dt
			dy := by1 / mag * dt

			// RK2: midpoint
			bz2, by2, ok := interp(z+dz*0.5, y+dy*0.5)
			if !ok {
				break
			}
			mag2 := math.Sqrt(bz2*bz2 + by2*by2)
			if mag2 < 1e-20 {
				break
			}
			z += bz2 / mag2 * dt
			y += by2 / mag2 * dt
			line = append(line, Point{z * 100, y * 100})
		}

		if len(line) > 5 {
			lines = append(lines, line)
		}
	}

	return lines
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
