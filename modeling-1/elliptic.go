package main

import "math"

// EllipticK computes the complete elliptic integral of the first kind K(m),
// where m = k^2 is the parameter. Uses the AGM (arithmetic-geometric mean) method.
func EllipticK(m float64) float64 {
	if m < 0 {
		m = 0
	}
	if m >= 1 {
		m = 1 - 1e-15
	}
	a := 1.0
	g := math.Sqrt(1 - m)
	for math.Abs(a-g) > 1e-15 {
		a, g = (a+g)/2, math.Sqrt(a*g)
	}
	return math.Pi / (2 * a)
}

// EllipticE computes the complete elliptic integral of the second kind E(m),
// where m = k^2 is the parameter. Uses the AGM method with c_n tracking.
func EllipticE(m float64) float64 {
	if m < 0 {
		m = 0
	}
	if m >= 1 {
		return 1.0
	}
	if m == 0 {
		return math.Pi / 2
	}

	a := 1.0
	g := math.Sqrt(1 - m)
	sum := m
	pow2 := 1.0

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

	k := math.Pi / (2 * a)
	return k * (1 - sum/2)
}
