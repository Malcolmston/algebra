package autodiff

import "math"

// This file implements reverse-mode automatic differentiation (backpropagation)
// for scalar-valued functions via a small operation [Tape].

// Tape records the elementary operations of a computation in evaluation order so
// that adjoints can be propagated backward from an output to every registered
// input in a single sweep. Reverse mode computes the whole gradient of a scalar
// function at a cost proportional to one function evaluation, independent of the
// number of inputs, which is what makes it the method of choice for
// high-dimensional optimization. A Tape holds no global state; construct one
// per computation with [NewTape].
type Tape struct {
	vals   []float64 // primal value of each node
	p1, p2 []int     // parent node indices, -1 when absent
	w1, w2 []float64 // local partial derivatives dnode/dparent
	inputs []int     // node indices registered as independent variables
}

// Var is a handle to a single node on a [Tape]. Operations on Var record their
// structure onto the tape as a side effect and return a new Var for the result.
// A Var is only meaningful together with the Tape that produced it.
type Var struct {
	t  *Tape
	id int
}

// NewTape returns a new, empty operation tape ready to record a computation.
func NewTape() *Tape {
	return &Tape{}
}

// autodiffPushNode appends a node with the given value and up to two parents to
// the tape and returns its index. Absent parents are passed as -1.
func autodiffPushNode(t *Tape, val float64, p1 int, w1 float64, p2 int, w2 float64) int {
	id := len(t.vals)
	t.vals = append(t.vals, val)
	t.p1 = append(t.p1, p1)
	t.w1 = append(t.w1, w1)
	t.p2 = append(t.p2, p2)
	t.w2 = append(t.w2, w2)
	return id
}

// Constant records a constant with value v on the tape. Its adjoint is computed
// but it is not returned by [Tape.Backward] because it is not an independent
// variable.
func (t *Tape) Constant(v float64) Var {
	return Var{t: t, id: autodiffPushNode(t, v, -1, 0, -1, 0)}
}

// Variable records an independent variable with value v on the tape and
// registers it so that its adjoint appears in the gradient returned by
// [Tape.Backward], in registration order.
func (t *Tape) Variable(v float64) Var {
	id := autodiffPushNode(t, v, -1, 0, -1, 0)
	t.inputs = append(t.inputs, id)
	return Var{t: t, id: id}
}

// NumVars returns the number of independent variables registered on the tape.
func (t *Tape) NumVars() int { return len(t.inputs) }

// Backward propagates adjoints from the output node y backward through the tape
// and returns the gradient of y with respect to every registered independent
// variable, in the order the variables were created. Because nodes are appended
// after their parents, a single reverse sweep over the tape suffices.
func (t *Tape) Backward(y Var) []float64 {
	adj := make([]float64, len(t.vals))
	adj[y.id] = 1
	for i := len(t.vals) - 1; i >= 0; i-- {
		a := adj[i]
		if a == 0 {
			continue
		}
		if p := t.p1[i]; p >= 0 {
			adj[p] += a * t.w1[i]
		}
		if p := t.p2[i]; p >= 0 {
			adj[p] += a * t.w2[i]
		}
	}
	grad := make([]float64, len(t.inputs))
	for k, idx := range t.inputs {
		grad[k] = adj[idx]
	}
	return grad
}

// Value returns the primal value stored at this node.
func (v Var) Value() float64 { return v.t.vals[v.id] }

// autodiffUnary records a unary operation with result value and local
// derivative deriv = dresult/dv onto v's tape.
func (v Var) autodiffUnary(value, deriv float64) Var {
	return Var{t: v.t, id: autodiffPushNode(v.t, value, v.id, deriv, -1, 0)}
}

// Add records and returns the sum v + w.
func (v Var) Add(w Var) Var {
	return Var{t: v.t, id: autodiffPushNode(v.t, v.Value()+w.Value(), v.id, 1, w.id, 1)}
}

// Sub records and returns the difference v - w.
func (v Var) Sub(w Var) Var {
	return Var{t: v.t, id: autodiffPushNode(v.t, v.Value()-w.Value(), v.id, 1, w.id, -1)}
}

// Mul records and returns the product v · w, whose local partials are the
// opposite operand's value.
func (v Var) Mul(w Var) Var {
	return Var{t: v.t, id: autodiffPushNode(v.t, v.Value()*w.Value(), v.id, w.Value(), w.id, v.Value())}
}

// Div records and returns the quotient v / w with partials 1/w and -v/w².
func (v Var) Div(w Var) Var {
	a, b := v.Value(), w.Value()
	return Var{t: v.t, id: autodiffPushNode(v.t, a/b, v.id, 1/b, w.id, -a/(b*b))}
}

// Neg records and returns the negation -v.
func (v Var) Neg() Var { return v.autodiffUnary(-v.Value(), -1) }

// Scale records and returns k·v for a real constant k.
func (v Var) Scale(k float64) Var { return v.autodiffUnary(k*v.Value(), k) }

// AddReal records and returns v + k for a real constant k.
func (v Var) AddReal(k float64) Var { return v.autodiffUnary(v.Value()+k, 1) }

// Recip records and returns the reciprocal 1/v with derivative -1/v².
func (v Var) Recip() Var {
	a := v.Value()
	return v.autodiffUnary(1/a, -1/(a*a))
}

// PowReal records and returns v^p for a constant real exponent p, with
// derivative p·v^(p-1).
func (v Var) PowReal(p float64) Var {
	a := v.Value()
	return v.autodiffUnary(math.Pow(a, p), p*math.Pow(a, p-1))
}

// Exp records and returns e^v.
func (v Var) Exp() Var {
	e := math.Exp(v.Value())
	return v.autodiffUnary(e, e)
}

// Log records and returns ln(v) with derivative 1/v.
func (v Var) Log() Var {
	a := v.Value()
	return v.autodiffUnary(math.Log(a), 1/a)
}

// Sqrt records and returns √v with derivative 1/(2√v).
func (v Var) Sqrt() Var {
	r := math.Sqrt(v.Value())
	return v.autodiffUnary(r, 0.5/r)
}

// Sin records and returns sin(v) with derivative cos(v).
func (v Var) Sin() Var {
	s, c := math.Sincos(v.Value())
	return v.autodiffUnary(s, c)
}

// Cos records and returns cos(v) with derivative -sin(v).
func (v Var) Cos() Var {
	s, c := math.Sincos(v.Value())
	return v.autodiffUnary(c, -s)
}

// Tan records and returns tan(v) with derivative sec²(v).
func (v Var) Tan() Var {
	t := math.Tan(v.Value())
	return v.autodiffUnary(t, 1+t*t)
}

// Tanh records and returns tanh(v) with derivative 1-tanh²(v).
func (v Var) Tanh() Var {
	t := math.Tanh(v.Value())
	return v.autodiffUnary(t, 1-t*t)
}

// Sigmoid records and returns the logistic function σ(v) with derivative
// σ(v)·(1-σ(v)).
func (v Var) Sigmoid() Var {
	s := 1 / (1 + math.Exp(-v.Value()))
	return v.autodiffUnary(s, s*(1-s))
}

// Abs records and returns |v| with derivative sign(v); the derivative at zero
// is taken as zero.
func (v Var) Abs() Var {
	a := v.Value()
	switch {
	case a > 0:
		return v.autodiffUnary(a, 1)
	case a < 0:
		return v.autodiffUnary(-a, -1)
	default:
		return v.autodiffUnary(0, 0)
	}
}

// GradientReverse returns the gradient of the scalar function f: Rⁿ → R at x
// using a single reverse-mode sweep. The function f receives a fresh [Tape] and
// the input variables already registered on it, and must return the scalar
// output built from those variables. This is the reverse-mode counterpart of
// [Gradient] and is preferable when n is large.
func GradientReverse(f func(t *Tape, x []Var) Var, x []float64) []float64 {
	t := NewTape()
	vars := make([]Var, len(x))
	for i := range x {
		vars[i] = t.Variable(x[i])
	}
	y := f(t, vars)
	return t.Backward(y)
}
