package main

import (
	"math/rand"
	"testing"
)

func benchFFT2D(b *testing.B, N int) {
	rng := rand.New(rand.NewSource(1))
	src := make([]complex128, N*N)
	for i := range src {
		src[i] = complex(rng.Float64(), rng.Float64())
	}
	buf := make([]complex128, N*N)
	twiddles(N)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		copy(buf, src)
		fft2D(buf, N)
	}
}

func BenchmarkFFT2D_256(b *testing.B)  { benchFFT2D(b, 256) }
func BenchmarkFFT2D_512(b *testing.B)  { benchFFT2D(b, 512) }
func BenchmarkFFT2D_1024(b *testing.B) { benchFFT2D(b, 1024) }
