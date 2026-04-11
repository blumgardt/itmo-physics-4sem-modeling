package main

import (
	"math"
	"testing"
)

func TestEllipticK(t *testing.T) {
	tests := []struct {
		m    float64
		want float64
	}{
		{0.0, math.Pi / 2},
		{0.5, 1.8540746773013719},
		{0.25, 1.6857503548325898},
		{0.75, 2.1565156474996432},
	}
	for _, tt := range tests {
		got := EllipticK(tt.m)
		if math.Abs(got-tt.want) > 1e-10 {
			t.Errorf("EllipticK(%v) = %v, want %v", tt.m, got, tt.want)
		}
	}
}

func TestEllipticE(t *testing.T) {
	tests := []struct {
		m    float64
		want float64
	}{
		{0.0, math.Pi / 2},
		{1.0, 1.0},
		{0.5, 1.3506438810476755},
		{0.25, 1.4674622093394272},
		{0.75, 1.2110560275684594},
	}
	for _, tt := range tests {
		got := EllipticE(tt.m)
		if math.Abs(got-tt.want) > 1e-10 {
			t.Errorf("EllipticE(%v) = %v, want %v", tt.m, got, tt.want)
		}
	}
}
