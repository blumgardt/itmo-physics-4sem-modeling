package main

import (
	"math"
	"math/cmplx"
	"testing"
)

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

	predicted := float64(N) / float64(2*a+1)
	rowOffset := cy * N

	bestIdx := -1
	bestVal := Imax
	for k := cx + 2; k < N; k++ {
		c := data[rowOffset+k]
		I := real(c)*real(c) + imag(c)*imag(c)
		if I < bestVal {
			bestVal = I
			bestIdx = k
		}
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
	if bestVal/Imax > 0.01 {
		t.Fatalf("first minimum not deep enough: I/Imax = %g at bin offset %.1f",
			bestVal/Imax, off)
	}
}

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
