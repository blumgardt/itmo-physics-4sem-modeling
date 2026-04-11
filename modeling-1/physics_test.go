package main

import (
	"math"
	"testing"
)

func TestLongSolenoidApproachesBTheory(t *testing.T) {
	// For a very long solenoid, B at center should approach mu0 * n * I
	p := Params{
		D: 0.1,
		L: 2.0, // long solenoid (L >> D)
		N: 1000,
		I: 1.0,
	}
	result := ComputeField(p)
	ratio := result.BCenter / result.BTheory
	if math.Abs(ratio-1) > 0.02 {
		t.Errorf("B_center/B_theory = %v, want ~1.0 (within 2%%)", ratio)
	}
}

func TestDefaultParams(t *testing.T) {
	p := Params{D: 0.1, L: 0.3, N: 100, I: 1.0}
	result := ComputeField(p)

	if len(result.ZAxis) != 121 {
		t.Errorf("expected 121 z points, got %d", len(result.ZAxis))
	}
	if len(result.YAxis) != 81 {
		t.Errorf("expected 81 y points, got %d", len(result.YAxis))
	}
	if result.BTheory <= 0 {
		t.Error("B_theory should be positive")
	}
	if result.BCenter <= 0 {
		t.Error("B_center should be positive")
	}
	if len(result.Streamlines) == 0 {
		t.Error("expected some streamlines")
	}
}

func TestTurnCentersStayInsidePhysicalLength(t *testing.T) {
	got := turnCenters(0.4, 4)
	want := []float64{-0.15, -0.05, 0.05, 0.15}

	if len(got) != len(want) {
		t.Fatalf("len(turnCenters) = %d, want %d", len(got), len(want))
	}
	for i := range want {
		if math.Abs(got[i]-want[i]) > 1e-12 {
			t.Fatalf("turnCenters[%d] = %v, want %v", i, got[i], want[i])
		}
	}
}

func TestGridContainsExactCenter(t *testing.T) {
	p := Params{D: 0.1, L: 0.3, N: 100, I: 1.0}
	result := ComputeField(p)
	midZ := len(result.ZAxis) / 2
	midY := len(result.YAxis) / 2

	if math.Abs(result.ZAxis[midZ]) > 1e-12 {
		t.Fatalf("z axis midpoint = %v, want 0", result.ZAxis[midZ])
	}
	if math.Abs(result.YAxis[midY]) > 1e-12 {
		t.Fatalf("y axis midpoint = %v, want 0", result.YAxis[midY])
	}

	centerBr, centerBz := solenoidField(p.D, p.L, p.N, p.I, []float64{0}, []float64{0})
	want := math.Hypot(centerBr[0], centerBz[0])
	if math.Abs(result.BCenter-want) > 1e-15 {
		t.Fatalf("B_center = %v, want %v", result.BCenter, want)
	}
	if math.Abs(result.BMag[midY][midZ]-want) > 1e-15 {
		t.Fatalf("B_mag at center = %v, want %v", result.BMag[midY][midZ], want)
	}
}

func TestGridUsesSquareViewExtent(t *testing.T) {
	p := Params{D: 0.1, L: 0.05, N: 100, I: 1.0, Scale: 4.5}
	result := ComputeField(p)
	wantExtentCm := p.Scale * math.Max(p.D, p.L) * 100

	if math.Abs(result.ZAxis[0]+wantExtentCm) > 1e-9 {
		t.Fatalf("z min = %v, want %v", result.ZAxis[0], -wantExtentCm)
	}
	if math.Abs(result.ZAxis[len(result.ZAxis)-1]-wantExtentCm) > 1e-9 {
		t.Fatalf("z max = %v, want %v", result.ZAxis[len(result.ZAxis)-1], wantExtentCm)
	}
	if math.Abs(result.YAxis[0]+wantExtentCm) > 1e-9 {
		t.Fatalf("y min = %v, want %v", result.YAxis[0], -wantExtentCm)
	}
	if math.Abs(result.YAxis[len(result.YAxis)-1]-wantExtentCm) > 1e-9 {
		t.Fatalf("y max = %v, want %v", result.YAxis[len(result.YAxis)-1], wantExtentCm)
	}
}
