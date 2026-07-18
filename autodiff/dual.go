package autodiff

import (
	"math"
	"strconv"
)

// Dual is a forward-mode dual number of the form Val + Der*ε, where ε is a
// nilpotent infinitesimal satisfying ε² = 0. Evaluating a function f on the
// dual number x = a + ε yields f(a) + f'(a)*ε, so the Der field of the result
// holds the derivative of f at a. More generally, if the input carries a
// derivative seed s (its Der field), the output's Der field holds f'(a)*s by
// the chain rule.
type Dual struct {
	// Val is the primal value of the number.
	Val float64
	// Der is the derivative (the coefficient of the infinitesimal ε).
	Der float64
}

// NewDual returns the dual number val + der*ε.
func NewDual(val, der float64) Dual {
	return Dual{Val: val, Der: der}
}

// Constant returns the dual number val + 0*ε, representing a quantity whose
// derivative with respect to the independent variable is zero.
func Constant(val float64) Dual {
	return Dual{Val: val, Der: 0}
}

// Variable returns the dual number val + 1*ε, the seed for the independent
// variable of a forward-mode differentiation. The unit derivative propagates
// the chain rule so that the output's Der field holds a true derivative.
func Variable(val float64) Dual {
	return Dual{Val: val, Der: 1}
}

// Real returns the primal (real) part of the dual number, discarding the
// derivative. It is shorthand for the Val field.
func (x Dual) Real() float64 { return x.Val }

// Add returns the sum x + y.
func (x Dual) Add(y Dual) Dual {
	return Dual{Val: x.Val + y.Val, Der: x.Der + y.Der}
}

// Sub returns the difference x - y.
func (x Dual) Sub(y Dual) Dual {
	return Dual{Val: x.Val - y.Val, Der: x.Der - y.Der}
}

// Mul returns the product x * y using the product rule
// (a*b)' = a'*b + a*b'.
func (x Dual) Mul(y Dual) Dual {
	return Dual{
		Val: x.Val * y.Val,
		Der: x.Der*y.Val + x.Val*y.Der,
	}
}

// Div returns the quotient x / y using the quotient rule
// (a/b)' = (a'*b - a*b')/b². It returns a number with NaN or infinite parts
// if y.Val is zero, matching IEEE-754 division.
func (x Dual) Div(y Dual) Dual {
	inv := 1 / y.Val
	val := x.Val * inv
	return Dual{
		Val: val,
		Der: (x.Der - val*y.Der) * inv,
	}
}

// Neg returns the negation -x.
func (x Dual) Neg() Dual {
	return Dual{Val: -x.Val, Der: -x.Der}
}

// Scale returns x multiplied by the real scalar k, differentiating k as a
// constant so that (k*x)' = k*x'.
func (x Dual) Scale(k float64) Dual {
	return Dual{Val: k * x.Val, Der: k * x.Der}
}

// AddReal returns x + k where k is a real constant, leaving the derivative
// unchanged.
func (x Dual) AddReal(k float64) Dual {
	return Dual{Val: x.Val + k, Der: x.Der}
}

// Inv returns the reciprocal 1/x, whose derivative is -x'/x².
func (x Dual) Inv() Dual {
	inv := 1 / x.Val
	return Dual{Val: inv, Der: -x.Der * inv * inv}
}

// String formats the dual number in the canonical "a+bε" notation.
func (x Dual) String() string {
	b := strconv.FormatFloat(x.Der, 'g', -1, 64)
	sign := "+"
	if x.Der < 0 || math.Signbit(x.Der) {
		sign = "-"
		b = strconv.FormatFloat(-x.Der, 'g', -1, 64)
	}
	return strconv.FormatFloat(x.Val, 'g', -1, 64) + sign + b + "ε"
}

// autodiffDualApply lifts a differentiable unary real function f to the duals.
// Given f(a)=f0 and f'(a)=f1 at a = x.Val, it returns f(x) with the chain rule
// applied to the derivative slot.
func autodiffDualApply(x Dual, f0, f1 float64) Dual {
	return Dual{Val: f0, Der: f1 * x.Der}
}
