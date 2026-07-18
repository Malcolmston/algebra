package autodiff

import (
	"math"
	"testing"
)

func TestSecondDerivative(t *testing.T) {
	// f(x) = x⁴; f''(x) = 12x². At x=2: value 16, f'=32, f''=48.
	f := func(x HyperDual) HyperDual { return x.PowReal(4) }
	v, d1, d2 := Derivatives2(f, 2)
	if !closeTo(v, 16, tol) || !closeTo(d1, 32, tol) || !closeTo(d2, 48, tol) {
		t.Fatalf("got v=%v d1=%v d2=%v", v, d1, d2)
	}
	// f(x) = sin(x); f'' = -sin(x). At x=1.
	g := func(x HyperDual) HyperDual { return x.Sin() }
	if s := SecondDerivative(g, 1); !closeTo(s, -math.Sin(1), tol) {
		t.Fatalf("sin'' = %v want %v", s, -math.Sin(1))
	}
	// f(x) = exp(x)·ln(x): validate via product/second-derivative closed form.
	h := func(x HyperDual) HyperDual { return x.Exp().Mul(x.Log()) }
	x0 := 1.7
	// f = e^x ln x. f' = e^x(ln x + 1/x). f'' = e^x(ln x + 2/x - 1/x²).
	want := math.Exp(x0) * (math.Log(x0) + 2/x0 - 1/(x0*x0))
	if s := SecondDerivative(h, x0); !closeTo(s, want, 1e-8) {
		t.Fatalf("h'' = %v want %v", s, want)
	}
}

func TestHessian(t *testing.T) {
	// f(x,y) = x²y + y³.
	// grad = (2xy, x² + 3y²).
	// H = [[2y, 2x], [2x, 6y]].
	f := func(v []HyperDual) HyperDual {
		return v[0].Mul(v[0]).Mul(v[1]).Add(v[1].PowReal(3))
	}
	x := []float64{1, 2}
	g, h := GradientHessian(f, x)
	wantG := []float64{2 * x[0] * x[1], x[0]*x[0] + 3*x[1]*x[1]}
	for i := range wantG {
		if !closeTo(g[i], wantG[i], tol) {
			t.Fatalf("grad[%d]=%v want %v", i, g[i], wantG[i])
		}
	}
	wantH := [][]float64{{2 * x[1], 2 * x[0]}, {2 * x[0], 6 * x[1]}}
	for i := range wantH {
		for j := range wantH[i] {
			if !closeTo(h[i][j], wantH[i][j], tol) {
				t.Fatalf("H[%d][%d]=%v want %v", i, j, h[i][j], wantH[i][j])
			}
		}
	}
	// Hessian alone must agree.
	h2 := Hessian(f, x)
	for i := range wantH {
		for j := range wantH[i] {
			if !closeTo(h2[i][j], wantH[i][j], tol) {
				t.Fatalf("Hessian[%d][%d]=%v want %v", i, j, h2[i][j], wantH[i][j])
			}
		}
	}
}

func TestHessianVectorProduct(t *testing.T) {
	// Same f; H·v with v=(1,1) at (1,2): H=[[4,2],[2,12]] -> (6,14).
	f := func(v []HyperDual) HyperDual {
		return v[0].Mul(v[0]).Mul(v[1]).Add(v[1].PowReal(3))
	}
	hv := HessianVectorProduct(f, []float64{1, 2}, []float64{1, 1})
	if !closeTo(hv[0], 6, tol) || !closeTo(hv[1], 14, tol) {
		t.Fatalf("Hv=%v want [6 14]", hv)
	}
}
