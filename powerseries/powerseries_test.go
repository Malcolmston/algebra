package powerseries

import (
	"math"
	"testing"
)

const tol = 1e-9

// approxSlice fails the test unless got and want agree elementwise within tol.
func approxSlice(t *testing.T, name string, got, want []float64) {
	t.Helper()
	if len(got) < len(want) {
		t.Fatalf("%s: got length %d, want at least %d", name, len(got), len(want))
	}
	for i := range want {
		if math.Abs(got[i]-want[i]) > tol {
			t.Errorf("%s[%d] = %.12g, want %.12g", name, i, got[i], want[i])
		}
	}
}

func TestConstructorsAndAccessors(t *testing.T) {
	s := New(3, 0, -2, 5)
	if s.Len() != 4 {
		t.Fatalf("Len = %d, want 4", s.Len())
	}
	if s.Coeff(2) != -2 || s.Coeff(9) != 0 || s.Coeff(-1) != 0 {
		t.Errorf("Coeff lookups wrong: %v", s.Coeffs())
	}
	if got := s.Order(tol); got != 0 {
		t.Errorf("Order = %d, want 0", got)
	}
	if got := Monomial(4, 3, 6).Order(tol); got != 3 {
		t.Errorf("Order of monomial = %d, want 3", got)
	}
	if !Zero(5).IsZero(tol) {
		t.Errorf("Zero(5) should be zero")
	}
	if Ident(5).IsZero(tol) {
		t.Errorf("Ident should not be zero")
	}
	if got := One(4).Evaluate(2.5); math.Abs(got-1) > tol {
		t.Errorf("One.Evaluate = %g, want 1", got)
	}
	// Horner evaluation of 3 - 2x^2 + 5x^3 at x = 2 -> 3 - 8 + 40 = 35.
	if got := s.Evaluate(2); math.Abs(got-35) > tol {
		t.Errorf("Evaluate = %g, want 35", got)
	}
	fr := FromFunc(func(i int) float64 { return float64(i * i) }, 4)
	approxSlice(t, "FromFunc", fr.Coeffs(), []float64{0, 1, 4, 9})
	if !s.Truncate(2).Equal(New(3, 0), tol) {
		t.Errorf("Truncate wrong: %v", s.Truncate(2).Coeffs())
	}
	if s.Extend(6).Len() != 6 {
		t.Errorf("Extend length wrong")
	}
}

func TestArithmetic(t *testing.T) {
	a := New(1, 2, 3)
	b := New(4, 5)
	approxSlice(t, "Add", a.Add(b).Coeffs(), []float64{5, 7, 3})
	approxSlice(t, "Sub", a.Sub(b).Coeffs(), []float64{-3, -3, 3})
	approxSlice(t, "Neg", a.Neg().Coeffs(), []float64{-1, -2, -3})
	approxSlice(t, "Scale", a.Scale(2).Coeffs(), []float64{2, 4, 6})
	// (1+2x+3x^2)(4+5x) truncated to 3 terms = 4 + 13x + 22x^2.
	approxSlice(t, "Mul", a.Mul(b).Coeffs(), []float64{4, 13, 22})
	approxSlice(t, "Hadamard", a.Hadamard(New(2, 2, 2)).Coeffs(), []float64{2, 4, 6})
	// x^2 * (1+2x+3x^2) truncated to 3 terms = x^2.
	approxSlice(t, "Shift", a.Shift(2).Coeffs(), []float64{0, 0, 1})
	// (1+x)^3 = 1 + 3x + 3x^2 + x^3.
	approxSlice(t, "Pow", New(1, 1, 0, 0).Pow(3).Coeffs(), []float64{1, 3, 3, 1})
	if !a.Pow(0).Equal(One(3), tol) {
		t.Errorf("Pow(0) should be one")
	}
}

func TestCalculus(t *testing.T) {
	// d/dx (1 + 2x + 3x^2 + 4x^3) = 2 + 6x + 12x^2.
	s := New(1, 2, 3, 4)
	approxSlice(t, "Derivative", s.Derivative().Coeffs(), []float64{2, 6, 12})
	// Integral of 2 + 6x + 12x^2 with constant 1 -> 1 + 2x + 3x^2 + 4x^3.
	approxSlice(t, "IntegralConst", s.Derivative().IntegralConst(1).Coeffs(), []float64{1, 2, 3, 4})
	if s.Integral().Coeff(0) != 0 {
		t.Errorf("Integral constant term should be zero")
	}
}

func TestInverseAndDiv(t *testing.T) {
	// 1/(1-x) = 1 + x + x^2 + ...
	inv := New(1, -1, 0, 0, 0).Inverse()
	approxSlice(t, "Inverse", inv.Coeffs(), []float64{1, 1, 1, 1, 1})
	// (1)/(1-x-x^2) = Fibonacci shifted: 1,1,2,3,5.
	den := New(1, -1, -1, 0, 0)
	approxSlice(t, "Div", One(5).Div(den).Coeffs(), []float64{1, 1, 2, 3, 5})
}

func TestExpLogSqrt(t *testing.T) {
	n := 12
	// exp(x) coefficients are 1/k!.
	e := Ident(n).Exp()
	want := make([]float64, n)
	for k := 0; k < n; k++ {
		want[k] = 1 / powerseriesFactorial(k)
	}
	approxSlice(t, "Exp", e.Coeffs(), want)
	// log(exp(x)) = x.
	approxSlice(t, "Log∘Exp", e.Log().Coeffs(), Ident(n).Coeffs())
	// exp(log(1+x)) = 1+x.
	oneX := New(1, 1).Extend(n)
	approxSlice(t, "Exp∘Log", oneX.Log().Exp().Coeffs(), oneX.Coeffs())
	// sqrt((1+x)^2) = 1+x.
	sq := oneX.Mul(oneX)
	approxSlice(t, "Sqrt", sq.Sqrt().Coeffs(), oneX.Coeffs())
	// PowReal: (1+x)^(1/3) cubed returns 1+x.
	cbrt := oneX.PowReal(1.0 / 3.0)
	approxSlice(t, "PowReal", cbrt.Pow(3).Coeffs(), oneX.Coeffs())
	// exp with non-zero constant term: exp(2 + x) has constant term e^2.
	shifted := New(2, 1).Extend(4).Exp()
	if math.Abs(shifted.Coeff(0)-math.Exp(2)) > tol {
		t.Errorf("Exp constant term = %g, want %g", shifted.Coeff(0), math.Exp(2))
	}
}

func TestCompose(t *testing.T) {
	n := 8
	// exp(x) composed with -x gives exp(-x): coefficients (-1)^k/k!.
	comp := ExpGF(n).Compose(Ident(n).Neg())
	want := make([]float64, n)
	for k := 0; k < n; k++ {
		want[k] = powerseriesAltFactInv(k)
	}
	approxSlice(t, "Compose", comp.Coeffs(), want)
	// 1/(1-x) composed with 2x = 1/(1-2x) = geometric ratio 2.
	g := GeometricGF(n).Compose(Ident(n).Scale(2))
	approxSlice(t, "Compose geo", g.Coeffs(), GeometricRatioGF(2, n).Coeffs())
}

func TestReversion(t *testing.T) {
	n := 9
	// Reversion of x + x^2 has coefficients (-1)^{k-1} C_{k-1}.
	f := New(0, 1, 1).Extend(n)
	rev := f.Reversion()
	cat := CatalanNumbers(n)
	want := make([]float64, n)
	for k := 1; k < n; k++ {
		s := 1.0
		if (k-1)%2 == 1 {
			s = -1
		}
		want[k] = s * cat[k-1]
	}
	approxSlice(t, "Reversion", rev.Coeffs(), want)
	// f(rev(x)) = x.
	approxSlice(t, "f∘rev", f.Compose(rev).Coeffs(), Ident(n).Coeffs())
	// Reversion of sin is arcsin: 0,1,0,1/6,0,3/40,...
	s := Ident(n).Sin().Reversion()
	approxSlice(t, "arcsin", s.Coeffs(), []float64{0, 1, 0, 1.0 / 6.0, 0, 3.0 / 40.0, 0, 15.0 / 336.0})
}

func TestTrig(t *testing.T) {
	n := 10
	sin := Ident(n).Sin()
	cos := Ident(n).Cos()
	// sin^2 + cos^2 = 1.
	approxSlice(t, "sin^2+cos^2", sin.Mul(sin).Add(cos.Mul(cos)).Coeffs(), One(n).Coeffs())
	// Known Maclaurin coefficients.
	approxSlice(t, "sin", sin.Coeffs(), []float64{0, 1, 0, -1.0 / 6, 0, 1.0 / 120, 0, -1.0 / 5040, 0, 1.0 / 362880})
	approxSlice(t, "cos", cos.Coeffs(), []float64{1, 0, -0.5, 0, 1.0 / 24, 0, -1.0 / 720, 0, 1.0 / 40320, 0})
	// cosh^2 - sinh^2 = 1.
	sh := Ident(n).Sinh()
	ch := Ident(n).Cosh()
	approxSlice(t, "cosh^2-sinh^2", ch.Mul(ch).Sub(sh.Mul(sh)).Coeffs(), One(n).Coeffs())
	// tan = sin/cos known coefficients.
	approxSlice(t, "tan", Ident(n).Tan().Coeffs(), []float64{0, 1, 0, 1.0 / 3, 0, 2.0 / 15, 0, 17.0 / 315})
	// atan known coefficients.
	approxSlice(t, "atan", Ident(n).Atan().Coeffs(), []float64{0, 1, 0, -1.0 / 3, 0, 1.0 / 5, 0, -1.0 / 7})
}
