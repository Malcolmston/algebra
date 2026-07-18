package autodiff

import (
	"math"
	"testing"
)

const tol = 1e-9

func closeTo(a, b, eps float64) bool {
	d := a - b
	if d < 0 {
		d = -d
	}
	return d <= eps*(1+math.Abs(b))
}

func TestDualArithmetic(t *testing.T) {
	// Differentiate h(x) = (x²+1)·(3x) / (x-4) at x = 2 and compare against a
	// closed-form derivative.
	x := Variable(2)
	num := x.Mul(x).AddReal(1).Mul(x.Scale(3))
	den := x.AddReal(-4)
	h := num.Div(den)
	// h(x) = 3x(x²+1)/(x-4); at x=2: 3*2*5 / (-2) = 30/-2 = -15.
	if !closeTo(h.Val, -15, tol) {
		t.Fatalf("value: got %v want -15", h.Val)
	}
	// h(x) = (3x³+3x)/(x-4). h'(x) = [(9x²+3)(x-4) - (3x³+3x)] / (x-4)².
	num2 := 2.0
	numer := (9*num2*num2 + 3) * (num2 - 4)
	numer -= 3*num2*num2*num2 + 3*num2
	want := numer / ((num2 - 4) * (num2 - 4))
	if !closeTo(h.Der, want, tol) {
		t.Fatalf("derivative: got %v want %v", h.Der, want)
	}
}

func TestDualInvNeg(t *testing.T) {
	x := Variable(3.0)
	got := x.Inv()
	if !closeTo(got.Val, 1.0/3.0, tol) || !closeTo(got.Der, -1.0/9.0, tol) {
		t.Fatalf("inv: got %+v", got)
	}
	n := x.Neg()
	if !closeTo(n.Val, -3, tol) || !closeTo(n.Der, -1, tol) {
		t.Fatalf("neg: got %+v", n)
	}
}

// TestElementaryDerivatives validates every dual elementary overload against an
// independently written closed-form derivative at several points.
func TestElementaryDerivatives(t *testing.T) {
	type tc struct {
		name string
		f    func(Dual) Dual
		val  func(float64) float64
		der  func(float64) float64
		xs   []float64
	}
	cases := []tc{
		{"Sin", Sin, math.Sin, math.Cos, []float64{-1, 0, 0.5, 1.3}},
		{"Cos", Cos, math.Cos, func(x float64) float64 { return -math.Sin(x) }, []float64{-1, 0, 0.5, 1.3}},
		{"Tan", Tan, math.Tan, func(x float64) float64 { c := math.Cos(x); return 1 / (c * c) }, []float64{-1, 0.5, 1.3}},
		{"Exp", Exp, math.Exp, math.Exp, []float64{-1, 0, 1, 2}},
		{"Exp2", Exp2, math.Exp2, func(x float64) float64 { return math.Ln2 * math.Exp2(x) }, []float64{-1, 0, 3}},
		{"Expm1", Expm1, math.Expm1, math.Exp, []float64{-0.5, 0, 0.2}},
		{"Log", Log, math.Log, func(x float64) float64 { return 1 / x }, []float64{0.3, 1, 5}},
		{"Log2", Log2, math.Log2, func(x float64) float64 { return 1 / (x * math.Ln2) }, []float64{0.3, 1, 5}},
		{"Log10", Log10, math.Log10, func(x float64) float64 { return 1 / (x * math.Ln10) }, []float64{0.3, 1, 5}},
		{"Log1p", Log1p, math.Log1p, func(x float64) float64 { return 1 / (1 + x) }, []float64{-0.5, 0, 2}},
		{"Sqrt", Sqrt, math.Sqrt, func(x float64) float64 { return 0.5 / math.Sqrt(x) }, []float64{0.25, 1, 9}},
		{"Cbrt", Cbrt, math.Cbrt, func(x float64) float64 { return 1.0 / (3 * math.Cbrt(x*x)) }, []float64{1, 8, 27}},
		{"Sinh", Sinh, math.Sinh, math.Cosh, []float64{-1, 0, 1}},
		{"Cosh", Cosh, math.Cosh, math.Sinh, []float64{-1, 0, 1}},
		{"Tanh", Tanh, math.Tanh, func(x float64) float64 { t := math.Tanh(x); return 1 - t*t }, []float64{-1, 0, 1}},
		{"Asin", Asin, math.Asin, func(x float64) float64 { return 1 / math.Sqrt(1-x*x) }, []float64{-0.5, 0, 0.5}},
		{"Acos", Acos, math.Acos, func(x float64) float64 { return -1 / math.Sqrt(1-x*x) }, []float64{-0.5, 0, 0.5}},
		{"Atan", Atan, math.Atan, func(x float64) float64 { return 1 / (1 + x*x) }, []float64{-2, 0, 2}},
		{"Asinh", Asinh, math.Asinh, func(x float64) float64 { return 1 / math.Sqrt(x*x+1) }, []float64{-1, 0, 2}},
		{"Acosh", Acosh, math.Acosh, func(x float64) float64 { return 1 / math.Sqrt(x*x-1) }, []float64{1.5, 2, 4}},
		{"Atanh", Atanh, math.Atanh, func(x float64) float64 { return 1 / (1 - x*x) }, []float64{-0.5, 0, 0.5}},
		{"Erf", Erf, math.Erf, func(x float64) float64 { return 2 / math.SqrtPi * math.Exp(-x*x) }, []float64{-1, 0, 1}},
		{"Erfc", Erfc, math.Erfc, func(x float64) float64 { return -2 / math.SqrtPi * math.Exp(-x*x) }, []float64{-1, 0, 1}},
		{"Sigmoid", Sigmoid, func(x float64) float64 { return 1 / (1 + math.Exp(-x)) },
			func(x float64) float64 { s := 1 / (1 + math.Exp(-x)); return s * (1 - s) }, []float64{-2, 0, 2}},
		{"Softplus", Softplus, func(x float64) float64 { return math.Log1p(math.Exp(x)) },
			func(x float64) float64 { return 1 / (1 + math.Exp(-x)) }, []float64{-2, 0, 2}},
		{"Cot", Cot, func(x float64) float64 { return math.Cos(x) / math.Sin(x) },
			func(x float64) float64 { s := math.Sin(x); return -1 / (s * s) }, []float64{0.5, 1.3, 2}},
		{"Sec", Sec, func(x float64) float64 { return 1 / math.Cos(x) },
			func(x float64) float64 { return math.Tan(x) / math.Cos(x) }, []float64{0.3, 1, 1.3}},
		{"Csc", Csc, func(x float64) float64 { return 1 / math.Sin(x) },
			func(x float64) float64 { return -math.Cos(x) / (math.Sin(x) * math.Sin(x)) }, []float64{0.5, 1.3, 2}},
	}
	for _, c := range cases {
		for _, x := range c.xs {
			got := c.f(Variable(x))
			if !closeTo(got.Val, c.val(x), tol) {
				t.Errorf("%s value at %v: got %v want %v", c.name, x, got.Val, c.val(x))
			}
			if !closeTo(got.Der, c.der(x), tol) {
				t.Errorf("%s deriv at %v: got %v want %v", c.name, x, got.Der, c.der(x))
			}
		}
	}
}

func TestPowFamily(t *testing.T) {
	// PowReal: d/dx x^2.5 = 2.5 x^1.5 at x=4 -> 2.5*8 = 20; value 32.
	g := PowReal(Variable(4), 2.5)
	if !closeTo(g.Val, 32, tol) || !closeTo(g.Der, 20, tol) {
		t.Fatalf("PowReal: got %+v", g)
	}
	// Pow with x=y=Variable: f(x)=x^x, f'=x^x(ln x + 1). At x=2: 4*(ln2+1).
	x := Variable(2)
	p := Pow(x, x)
	want := 4 * (math.Log(2) + 1)
	if !closeTo(p.Val, 4, tol) || !closeTo(p.Der, want, tol) {
		t.Fatalf("Pow x^x: got %+v want der %v", p, want)
	}
	// Powf: d/dx 3^x = ln3 * 3^x at x=2 -> ln3*9; value 9.
	f := Powf(3, Variable(2))
	if !closeTo(f.Val, 9, tol) || !closeTo(f.Der, math.Log(3)*9, tol) {
		t.Fatalf("Powf: got %+v", f)
	}
}

func TestAtan2Hypot(t *testing.T) {
	// atan2(y,x) with y=Variable, x const: d/dy = x/(x²+y²). At x=3,y=4: 3/25.
	a := Atan2(Variable(4), Constant(3))
	if !closeTo(a.Val, math.Atan2(4, 3), tol) || !closeTo(a.Der, 3.0/25.0, tol) {
		t.Fatalf("Atan2: got %+v", a)
	}
	// hypot(x,y) with x=Variable: d/dx = x/hypot. At x=3,y=4: 3/5; value 5.
	h := Hypot(Variable(3), Constant(4))
	if !closeTo(h.Val, 5, tol) || !closeTo(h.Der, 0.6, tol) {
		t.Fatalf("Hypot: got %+v", h)
	}
}

func TestAbsReluMaxMin(t *testing.T) {
	if g := Abs(Variable(-2)); !closeTo(g.Val, 2, tol) || !closeTo(g.Der, -1, tol) {
		t.Fatalf("Abs(-2): %+v", g)
	}
	if g := Relu(Variable(-1)); g.Val != 0 || g.Der != 0 {
		t.Fatalf("Relu(-1): %+v", g)
	}
	if g := Relu(Variable(3)); !closeTo(g.Val, 3, tol) || !closeTo(g.Der, 1, tol) {
		t.Fatalf("Relu(3): %+v", g)
	}
	if g := Max(Variable(1), NewDual(5, 7)); g.Val != 5 || g.Der != 7 {
		t.Fatalf("Max picked wrong: %+v", g)
	}
	if g := Min(NewDual(1, 9), Constant(5)); g.Val != 1 || g.Der != 9 {
		t.Fatalf("Min picked wrong: %+v", g)
	}
}

func TestDualString(t *testing.T) {
	if s := NewDual(2, 3).String(); s != "2+3ε" {
		t.Fatalf("string: %q", s)
	}
	if s := NewDual(2, -3).String(); s != "2-3ε" {
		t.Fatalf("string: %q", s)
	}
}
