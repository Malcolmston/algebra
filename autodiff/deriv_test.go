package autodiff

import (
	"math"
	"testing"
)

func TestDerivative(t *testing.T) {
	// f(x) = x³ - 2x; f'(x) = 3x² - 2. At x=2: value 4, deriv 10.
	f := func(x Dual) Dual { return x.Mul(x).Mul(x).Sub(x.Scale(2)) }
	v, d := ValueAndDerivative(f, 2)
	if !closeTo(v, 4, tol) || !closeTo(d, 10, tol) {
		t.Fatalf("got v=%v d=%v", v, d)
	}
	if got := Derivative(f, 2); !closeTo(got, 10, tol) {
		t.Fatalf("Derivative got %v", got)
	}
}

func TestGradient(t *testing.T) {
	// f(x,y,z) = x²y + sin(z). grad = (2xy, x², cos z).
	f := func(v []Dual) Dual {
		return v[0].Mul(v[0]).Mul(v[1]).Add(Sin(v[2]))
	}
	x := []float64{1.5, 2, 0.7}
	g := Gradient(f, x)
	want := []float64{2 * x[0] * x[1], x[0] * x[0], math.Cos(x[2])}
	for i := range want {
		if !closeTo(g[i], want[i], tol) {
			t.Fatalf("grad[%d]=%v want %v", i, g[i], want[i])
		}
	}
	val, g2 := ValueAndGradient(f, x)
	wantVal := x[0]*x[0]*x[1] + math.Sin(x[2])
	if !closeTo(val, wantVal, tol) {
		t.Fatalf("value %v want %v", val, wantVal)
	}
	for i := range want {
		if !closeTo(g2[i], want[i], tol) {
			t.Fatalf("g2[%d]=%v want %v", i, g2[i], want[i])
		}
	}
}

func TestPartialAndDirectional(t *testing.T) {
	f := func(v []Dual) Dual { return v[0].Mul(v[0]).Add(v[1].Mul(v[1])) }
	x := []float64{1, 2}
	if p := PartialDerivative(f, x, 1); !closeTo(p, 4, tol) {
		t.Fatalf("partial=%v want 4", p)
	}
	// grad = (2,4). dir (3,4): directional = 2*3+4*4 = 22.
	if d := DirectionalDerivative(f, x, []float64{3, 4}); !closeTo(d, 22, tol) {
		t.Fatalf("directional=%v want 22", d)
	}
}

func TestJacobian(t *testing.T) {
	// f(x,y) = [x²y, x+sin y]. J = [[2xy, x²], [1, cos y]].
	f := func(v []Dual) []Dual {
		return []Dual{
			v[0].Mul(v[0]).Mul(v[1]),
			v[0].Add(Sin(v[1])),
		}
	}
	x := []float64{2, 0.5}
	j := Jacobian(f, x)
	want := [][]float64{
		{2 * x[0] * x[1], x[0] * x[0]},
		{1, math.Cos(x[1])},
	}
	for i := range want {
		for k := range want[i] {
			if !closeTo(j[i][k], want[i][k], tol) {
				t.Fatalf("J[%d][%d]=%v want %v", i, k, j[i][k], want[i][k])
			}
		}
	}
	// J·v with v=(1,1): row0 = 2xy + x², row1 = 1 + cos y.
	jv := JacobianVectorProduct(f, x, []float64{1, 1})
	if !closeTo(jv[0], want[0][0]+want[0][1], tol) || !closeTo(jv[1], want[1][0]+want[1][1], tol) {
		t.Fatalf("Jv=%v", jv)
	}
}
