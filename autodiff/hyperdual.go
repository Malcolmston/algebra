package autodiff

import "math"

// HyperDual is a second-order forward-mode number of the form
//
//	Val + E1·ε₁ + E2·ε₂ + E12·ε₁ε₂,
//
// where ε₁ and ε₂ are independent nilpotent infinitesimals with ε₁² = ε₂² = 0
// but ε₁ε₂ ≠ 0. Evaluating a twice-differentiable function f on the number
// a + ε₁ + ε₂ yields
//
//	f(a) + f'(a)·ε₁ + f'(a)·ε₂ + f''(a)·ε₁ε₂,
//
// so the E12 slot delivers the second derivative exactly, with no subtractive
// cancellation. Seeding the two first-order slots with distinct directions
// yields mixed second partial derivatives, which is what [Hessian] exploits.
type HyperDual struct {
	// Val is the primal value.
	Val float64
	// E1 is the coefficient of the first infinitesimal ε₁.
	E1 float64
	// E2 is the coefficient of the second infinitesimal ε₂.
	E2 float64
	// E12 is the coefficient of the mixed term ε₁ε₂, carrying second-order
	// information.
	E12 float64
}

// NewHyperDual returns the hyper-dual number val + e1·ε₁ + e2·ε₂ + e12·ε₁ε₂.
func NewHyperDual(val, e1, e2, e12 float64) HyperDual {
	return HyperDual{Val: val, E1: e1, E2: e2, E12: e12}
}

// HyperConstant returns the hyper-dual number val with all infinitesimal parts
// zero, representing a quantity with vanishing first and second derivatives.
func HyperConstant(val float64) HyperDual {
	return HyperDual{Val: val}
}

// HyperVariable returns the hyper-dual seed val + ε₁ + ε₂ for a single
// independent variable. After evaluating a scalar function, its E12 slot holds
// the second derivative and either first-order slot holds the first derivative.
func HyperVariable(val float64) HyperDual {
	return HyperDual{Val: val, E1: 1, E2: 1}
}

// Real returns the primal part of the hyper-dual number.
func (x HyperDual) Real() float64 { return x.Val }

// Add returns the sum x + y.
func (x HyperDual) Add(y HyperDual) HyperDual {
	return HyperDual{x.Val + y.Val, x.E1 + y.E1, x.E2 + y.E2, x.E12 + y.E12}
}

// Sub returns the difference x - y.
func (x HyperDual) Sub(y HyperDual) HyperDual {
	return HyperDual{x.Val - y.Val, x.E1 - y.E1, x.E2 - y.E2, x.E12 - y.E12}
}

// Mul returns the product x · y using the full second-order product rule,
// including the cross term E1·E2 + E2·E1 that feeds the mixed slot.
func (x HyperDual) Mul(y HyperDual) HyperDual {
	return HyperDual{
		Val: x.Val * y.Val,
		E1:  x.Val*y.E1 + x.E1*y.Val,
		E2:  x.Val*y.E2 + x.E2*y.Val,
		E12: x.Val*y.E12 + x.E1*y.E2 + x.E2*y.E1 + x.E12*y.Val,
	}
}

// Neg returns the negation -x.
func (x HyperDual) Neg() HyperDual {
	return HyperDual{-x.Val, -x.E1, -x.E2, -x.E12}
}

// Scale returns x multiplied by the real scalar k.
func (x HyperDual) Scale(k float64) HyperDual {
	return HyperDual{k * x.Val, k * x.E1, k * x.E2, k * x.E12}
}

// AddReal returns x + k where k is a real constant.
func (x HyperDual) AddReal(k float64) HyperDual {
	return HyperDual{x.Val + k, x.E1, x.E2, x.E12}
}

// Inv returns the reciprocal 1/x, obtained by applying the reciprocal function
// and its derivatives r(a)=1/a, r'(a)=-1/a², r”(a)=2/a³.
func (x HyperDual) Inv() HyperDual {
	a := x.Val
	return autodiffHyperApply(x, 1/a, -1/(a*a), 2/(a*a*a))
}

// Div returns the quotient x / y as x · y⁻¹.
func (x HyperDual) Div(y HyperDual) HyperDual {
	return x.Mul(y.Inv())
}

// autodiffHyperApply lifts a twice-differentiable unary real function f to the
// hyper-duals. Given f(a)=f0, f'(a)=f1 and f”(a)=f2 at a = x.Val, it applies
// the second-order chain rule to all infinitesimal slots.
func autodiffHyperApply(x HyperDual, f0, f1, f2 float64) HyperDual {
	return HyperDual{
		Val: f0,
		E1:  f1 * x.E1,
		E2:  f1 * x.E2,
		E12: f1*x.E12 + f2*x.E1*x.E2,
	}
}

// Sin returns sin(x); it uses sin' = cos and sin” = -sin.
func (x HyperDual) Sin() HyperDual {
	s, c := math.Sincos(x.Val)
	return autodiffHyperApply(x, s, c, -s)
}

// Cos returns cos(x); it uses cos' = -sin and cos” = -cos.
func (x HyperDual) Cos() HyperDual {
	s, c := math.Sincos(x.Val)
	return autodiffHyperApply(x, c, -s, -c)
}

// Tan returns tan(x); it uses tan' = sec² and tan” = 2·sec²·tan.
func (x HyperDual) Tan() HyperDual {
	t := math.Tan(x.Val)
	sec2 := 1 + t*t
	return autodiffHyperApply(x, t, sec2, 2*sec2*t)
}

// Exp returns e^x; every derivative of the exponential is itself.
func (x HyperDual) Exp() HyperDual {
	e := math.Exp(x.Val)
	return autodiffHyperApply(x, e, e, e)
}

// Log returns ln(x); it uses (ln)' = 1/x and (ln)” = -1/x².
func (x HyperDual) Log() HyperDual {
	a := x.Val
	return autodiffHyperApply(x, math.Log(a), 1/a, -1/(a*a))
}

// Sqrt returns √x; it uses derivatives 1/(2√x) and -1/(4·x^(3/2)).
func (x HyperDual) Sqrt() HyperDual {
	r := math.Sqrt(x.Val)
	return autodiffHyperApply(x, r, 0.5/r, -0.25/(r*x.Val))
}

// PowReal returns x^p for a constant real exponent p, with first derivative
// p·x^(p-1) and second derivative p·(p-1)·x^(p-2).
func (x HyperDual) PowReal(p float64) HyperDual {
	a := x.Val
	return autodiffHyperApply(x,
		math.Pow(a, p),
		p*math.Pow(a, p-1),
		p*(p-1)*math.Pow(a, p-2))
}

// Sinh returns sinh(x); it uses sinh' = cosh and sinh” = sinh.
func (x HyperDual) Sinh() HyperDual {
	s, c := math.Sinh(x.Val), math.Cosh(x.Val)
	return autodiffHyperApply(x, s, c, s)
}

// Cosh returns cosh(x); it uses cosh' = sinh and cosh” = cosh.
func (x HyperDual) Cosh() HyperDual {
	s, c := math.Sinh(x.Val), math.Cosh(x.Val)
	return autodiffHyperApply(x, c, s, c)
}

// Tanh returns tanh(x); it uses tanh' = 1-tanh² and tanh” = -2·tanh·(1-tanh²).
func (x HyperDual) Tanh() HyperDual {
	t := math.Tanh(x.Val)
	d := 1 - t*t
	return autodiffHyperApply(x, t, d, -2*t*d)
}

// Sigmoid returns the logistic function σ(x) = 1/(1+e^(-x)); it uses
// σ' = σ(1-σ) and σ” = σ(1-σ)(1-2σ).
func (x HyperDual) Sigmoid() HyperDual {
	s := 1 / (1 + math.Exp(-x.Val))
	d := s * (1 - s)
	return autodiffHyperApply(x, s, d, d*(1-2*s))
}
