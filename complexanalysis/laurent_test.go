package complexanalysis

import (
	"math"
	"testing"
)

func TestLaurentCoefficients(t *testing.T) {
	// exp(z)/z^2 = z^-2 + z^-1 + 1/2 + z/6 + ...
	f := func(z complex128) complex128 { return Exp(z) / (z * z) }
	cases := []struct {
		k    int
		want complex128
	}{
		{-2, 1},
		{-1, 1},
		{0, 0.5},
		{1, 1.0 / 6},
		{2, 1.0 / 24},
	}
	for _, c := range cases {
		got := LaurentCoefficient(f, 0, c.k, 1, 600)
		if !closeC(got, c.want, 1e-8) {
			t.Errorf("c_%d = %v, want %v", c.k, got, c.want)
		}
	}
	all := LaurentCoefficients(f, 0, -2, 2, 1, 600)
	if len(all) != 5 {
		t.Fatalf("expected 5 coefficients, got %d", len(all))
	}
	if !closeC(all[0], 1, 1e-8) {
		t.Errorf("LaurentCoefficients[0] (c_-2) = %v", all[0])
	}
}

func TestTaylorCoefficients(t *testing.T) {
	// exp: a_k = 1/k!.
	coeffs := TaylorCoefficients(Exp, 0, 6, 1, 400)
	for k, c := range coeffs {
		want := complex(1/factorialF(k), 0)
		if !closeC(c, want, 1e-8) {
			t.Errorf("a_%d = %v, want %v", k, c, want)
		}
	}
	if TaylorCoefficients(Exp, 0, 0, 1, 100) != nil {
		t.Error("count<=0 should give nil")
	}
}

func factorialF(n int) float64 {
	f := 1.0
	for i := 2; i <= n; i++ {
		f *= float64(i)
	}
	return f
}

func TestPowerSeriesEval(t *testing.T) {
	// 1 + 2(z-1) + 3(z-1)^2 at z=2 = 1+2+3 = 6.
	got := PowerSeriesEval([]complex128{1, 2, 3}, 1, 2)
	if !closeC(got, 6, testTol) {
		t.Errorf("PowerSeriesEval = %v, want 6", got)
	}
}

func TestAnalyticContinuation(t *testing.T) {
	// Continue exp from 0 to 0.5.
	got := AnalyticContinuation(Exp, 0, 0.5, 14, 1, 600)
	if !closeC(got, complex(math.Exp(0.5), 0), 1e-9) {
		t.Errorf("continuation exp(0.5) = %v, want %v", got, math.Exp(0.5))
	}
	// Continue along a path 0 -> 0.5 -> 1 and recover exp(1).
	path := []complex128{0, 0.5, 1}
	gotp := AnalyticContinuationPath(Exp, path, 14, 0.9, 600)
	if !closeC(gotp, complex(math.E, 0), 1e-8) {
		t.Errorf("path continuation exp(1) = %v, want %v", gotp, math.E)
	}
}

func BenchmarkLaurentCoefficient(b *testing.B) {
	f := func(z complex128) complex128 { return Exp(z) / (z * z) }
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		_ = LaurentCoefficient(f, 0, -1, 1, 1024)
	}
}
