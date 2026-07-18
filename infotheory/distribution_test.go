package infotheory

import (
	"errors"
	"testing"
)

func TestUniformDistribution(t *testing.T) {
	if got := UniformDistribution(0); got != nil {
		t.Errorf("UniformDistribution(0) = %v, want nil", got)
	}
	p := UniformDistribution(4)
	for _, v := range p {
		if !approx(v, 0.25, tol) {
			t.Errorf("UniformDistribution(4) entry = %v, want 0.25", v)
		}
	}
}

func TestNormalize(t *testing.T) {
	p, err := NormalizeWeights([]float64{1, 3})
	if err != nil {
		t.Fatal(err)
	}
	if !approx(p[0], 0.25, tol) || !approx(p[1], 0.75, tol) {
		t.Errorf("NormalizeWeights = %v, want [0.25 0.75]", p)
	}
	if _, err := NormalizeWeights([]float64{0, 0}); !errors.Is(err, ErrNotNormalizable) {
		t.Errorf("expected ErrNotNormalizable, got %v", err)
	}
	c, err := NormalizeCounts([]int{1, 1, 2})
	if err != nil {
		t.Fatal(err)
	}
	if !approx(c[2], 0.5, tol) {
		t.Errorf("NormalizeCounts = %v, want last 0.5", c)
	}
	if _, err := NormalizeCounts([]int{0}); !errors.Is(err, ErrNotNormalizable) {
		t.Errorf("expected ErrNotNormalizable, got %v", err)
	}
}

func TestIsProbabilityDistribution(t *testing.T) {
	if !IsProbabilityDistribution([]float64{0.2, 0.3, 0.5}, 1e-12) {
		t.Error("valid distribution rejected")
	}
	if IsProbabilityDistribution([]float64{0.2, 0.3}, 1e-12) {
		t.Error("under-sum distribution accepted")
	}
	if IsProbabilityDistribution([]float64{-0.1, 1.1}, 1e-12) {
		t.Error("negative entry accepted")
	}
}

func TestEmpiricalDistribution(t *testing.T) {
	values, dist, err := EmpiricalDistribution([]float64{1, 2, 2, 3, 3, 3})
	if err != nil {
		t.Fatal(err)
	}
	wantVals := []float64{1, 2, 3}
	for i, v := range wantVals {
		if values[i] != v {
			t.Errorf("values[%d] = %v, want %v", i, values[i], v)
		}
	}
	wantDist := []float64{1.0 / 6, 2.0 / 6, 3.0 / 6}
	for i, d := range wantDist {
		if !approx(dist[i], d, tol) {
			t.Errorf("dist[%d] = %v, want %v", i, dist[i], d)
		}
	}
	if _, _, err := EmpiricalDistribution(nil); !errors.Is(err, ErrNotNormalizable) {
		t.Errorf("expected ErrNotNormalizable, got %v", err)
	}
}
