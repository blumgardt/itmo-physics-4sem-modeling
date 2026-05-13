package main

import (
	"math"
	"math/cmplx"
	"testing"
)

// TestFFT1DImpulse verifies that FFT of a unit impulse is a flat spectrum (all ones).
func TestFFT1DImpulse(t *testing.T) {
	N := 16
	data := make([]complex128, N)
	data[0] = 1
	fft1D(data, twiddles(N))
	for i, v := range data {
		if cmplx.Abs(v-1) > 1e-12 {
			t.Fatalf("impulse FFT bin %d = %v, want 1", i, v)
		}
	}
}

// TestFFT1DConstant verifies FFT of constant is DC-only.
func TestFFT1DConstant(t *testing.T) {
	N := 32
	data := make([]complex128, N)
	for i := range data {
		data[i] = 1
	}
	fft1D(data, twiddles(N))
	for i, v := range data {
		want := complex(0, 0)
		if i == 0 {
			want = complex(float64(N), 0)
		}
		if cmplx.Abs(v-want) > 1e-9 {
			t.Fatalf("constant FFT bin %d = %v, want %v", i, v, want)
		}
	}
}

// TestFFT1DSine verifies FFT of a complex exponential e^{+i 2π k0 n/N} lands at bin k0.
// (Our forward FFT uses kernel e^{-i 2π kn/N}, so positive-frequency input peaks at +k0.)
func TestFFT1DSine(t *testing.T) {
	N := 64
	k0 := 5
	data := make([]complex128, N)
	for n := 0; n < N; n++ {
		theta := 2 * math.Pi * float64(k0) * float64(n) / float64(N)
		data[n] = complex(math.Cos(theta), math.Sin(theta))
	}
	fft1D(data, twiddles(N))
	for i, v := range data {
		want := 0.0
		if i == k0 {
			want = float64(N)
		}
		if math.Abs(cmplx.Abs(v)-want) > 1e-9 {
			t.Fatalf("bin %d = %v, |.|=%g, want %g", i, v, cmplx.Abs(v), want)
		}
	}
}

// TestFFT2DRectSincZero verifies that a rectangular aperture's |FFT|² profile
// has a clear sinc² zero near the predicted bin (with discretization tolerance).
//
// A discrete rect of (2a+1) ones along x ↔ Dirichlet kernel D(k) = sin(π(2a+1)k/N)/sin(πk/N).
// First zero at k = N/(2a+1). For a=8, N=128, that's k ≈ 7.5 from DC, so 64 + ~7 or ~8.
func TestFFT2DRectSincZero(t *testing.T) {
	N := 128
	a := 8
	data := make([]complex128, N*N)
	cx := N / 2
	cy := N / 2
	for j := 0; j < N; j++ {
		for i := 0; i < N; i++ {
			v := 0.0
			if abs(i-cx) <= a && abs(j-cy) <= a {
				v = 1
			}
			if (i+j)&1 == 1 {
				v = -v
			}
			data[j*N+i] = complex(v, 0)
		}
	}
	fft2D(data, N)

	Imax := 0.0
	for _, c := range data {
		I := real(c)*real(c) + imag(c)*imag(c)
		if I > Imax {
			Imax = I
		}
	}

	// Predicted first Dirichlet zero (bin from center)
	predicted := float64(N) / float64(2*a+1) // ≈ 7.53
	rowOffset := cy * N

	// Scan right half of central row, find first deep minimum.
	bestIdx := -1
	bestVal := Imax
	for k := cx + 2; k < N; k++ {
		c := data[rowOffset+k]
		I := real(c)*real(c) + imag(c)*imag(c)
		if I < bestVal {
			bestVal = I
			bestIdx = k
		}
		// stop search after second-order maximum
		if k-cx > int(predicted*1.5) {
			break
		}
	}
	if bestIdx < 0 {
		t.Fatalf("no minimum found in central row")
	}
	off := float64(bestIdx - cx)
	if math.Abs(off-predicted) > 1.5 {
		t.Fatalf("sinc² first zero offset = %.2f bins, predicted %.2f", off, predicted)
	}
	// For a discrete rect of width W=2a+1 with W∤N, the FFT samples the Dirichlet
	// kernel at integer bins straddling the true zero — local minimum reaches
	// at most ~1/W² of the peak rather than exact zero. Loose check is enough
	// to confirm sinc-shape behaviour.
	if bestVal/Imax > 0.01 {
		t.Fatalf("first minimum not deep enough: I/Imax = %g at bin offset %.1f",
			bestVal/Imax, off)
	}
}

// TestFFT2DCheckerboardCentersDC: a constant aperture, after checkerboard sign and FFT,
// must have its peak at index (N/2, N/2).
func TestFFT2DCheckerboardCentersDC(t *testing.T) {
	N := 32
	data := make([]complex128, N*N)
	for j := 0; j < N; j++ {
		for i := 0; i < N; i++ {
			v := 1.0
			if (i+j)&1 == 1 {
				v = -v
			}
			data[j*N+i] = complex(v, 0)
		}
	}
	fft2D(data, N)
	maxIdx := 0
	maxVal := 0.0
	for k, c := range data {
		I := real(c)*real(c) + imag(c)*imag(c)
		if I > maxVal {
			maxVal = I
			maxIdx = k
		}
	}
	wantIdx := (N/2)*N + N/2
	if maxIdx != wantIdx {
		t.Fatalf("DC peak at index %d (row %d col %d), want %d (row %d col %d)",
			maxIdx, maxIdx/N, maxIdx%N, wantIdx, N/2, N/2)
	}
}

func abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}
