package controltheory

import (
	"math"
	"math/cmplx"
)

// TransferFunction is a SISO continuous-time transfer function
// G(s) = Num(s) / Den(s) with real polynomial numerator and denominator stored
// in the ascending-power convention of [Poly].
type TransferFunction struct {
	// Num is the numerator polynomial N(s).
	Num Poly
	// Den is the denominator polynomial D(s).
	Den Poly
}

// NewTransferFunction builds a TransferFunction from ascending-power numerator
// and denominator coefficient slices. Both slices are copied.
func NewTransferFunction(num, den []float64) TransferFunction {
	return TransferFunction{Num: NewPoly(num...), Den: NewPoly(den...)}
}

// Order returns the order of the system, i.e. the degree of the denominator
// polynomial.
func (g TransferFunction) Order() int {
	return g.Den.Degree()
}

// IsProper reports whether the numerator degree does not exceed the
// denominator degree, the condition for a physically realizable system.
func (g TransferFunction) IsProper() bool {
	return g.Num.Degree() <= g.Den.Degree()
}

// IsStrictlyProper reports whether the numerator degree is strictly less than
// the denominator degree.
func (g TransferFunction) IsStrictlyProper() bool {
	return g.Num.Degree() < g.Den.Degree()
}

// Evaluate returns the complex value of G(s) at the complex point s. It returns
// a value with infinite magnitude when s is a pole (the denominator vanishes).
func (g TransferFunction) Evaluate(s complex128) complex128 {
	d := g.Den.EvalComplex(s)
	if d == 0 {
		return cmplx.Inf()
	}
	return g.Num.EvalComplex(s) / d
}

// FrequencyResponse returns G(jω), the complex response at angular frequency
// omega (in radians per second).
func (g TransferFunction) FrequencyResponse(omega float64) complex128 {
	return g.Evaluate(complex(0, omega))
}

// Poles returns the poles of the system, i.e. the roots of the denominator.
func (g TransferFunction) Poles() []complex128 {
	return g.Den.Roots()
}

// Zeros returns the zeros of the system, i.e. the roots of the numerator.
func (g TransferFunction) Zeros() []complex128 {
	return g.Num.Roots()
}

// DCGain returns the steady-state gain G(0) = Num(0)/Den(0). It returns
// +Inf when the denominator constant term is zero (a pole at the origin).
func (g TransferFunction) DCGain() float64 {
	d := g.Den.Eval(0)
	if d == 0 {
		return math.Inf(1)
	}
	return g.Num.Eval(0) / d
}

// IsStable reports whether every pole has a strictly negative real part, the
// condition for asymptotic stability of a continuous-time system.
func (g TransferFunction) IsStable() bool {
	for _, p := range g.Poles() {
		if real(p) >= 0 {
			return false
		}
	}
	return true
}

// Series returns the transfer function of g followed by h connected in series
// (cascade): the product G(s)·H(s).
func (g TransferFunction) Series(h TransferFunction) TransferFunction {
	return Series(g, h)
}

// Series returns the cascade connection G(s)·H(s) of two transfer functions.
func Series(g, h TransferFunction) TransferFunction {
	return TransferFunction{
		Num: g.Num.Mul(h.Num),
		Den: g.Den.Mul(h.Den),
	}
}

// Parallel returns the parallel connection G(s)+H(s) of two transfer
// functions, whose outputs are summed for a common input.
func Parallel(g, h TransferFunction) TransferFunction {
	num := g.Num.Mul(h.Den).Add(h.Num.Mul(g.Den))
	den := g.Den.Mul(h.Den)
	return TransferFunction{Num: num, Den: den}
}

// Feedback returns the closed-loop transfer function of forward path g with
// feedback path h. When sign is -1 (negative feedback) the result is
// G/(1+G·H); when sign is +1 (positive feedback) it is G/(1-G·H).
func Feedback(g, h TransferFunction, sign int) TransferFunction {
	// Closed loop = G.Num*H.Den / (G.Den*H.Den - sign*G.Num*H.Num).
	num := g.Num.Mul(h.Den)
	loop := g.Num.Mul(h.Num)
	den := g.Den.Mul(h.Den).Sub(loop.Scale(float64(sign)))
	return TransferFunction{Num: num, Den: den}
}

// UnityFeedback returns the closed-loop transfer function of forward path g
// with unity negative feedback, G/(1+G).
func UnityFeedback(g TransferFunction) TransferFunction {
	one := TransferFunction{Num: Poly{1}, Den: Poly{1}}
	return Feedback(g, one, -1)
}

// SecondOrderSystem returns the canonical second-order transfer function
// wn^2 / (s^2 + 2·zeta·wn·s + wn^2) for natural frequency wn (rad/s) and
// damping ratio zeta.
func SecondOrderSystem(wn, zeta float64) TransferFunction {
	return TransferFunction{
		Num: Poly{wn * wn},
		Den: Poly{wn * wn, 2 * zeta * wn, 1},
	}
}

// NaturalFrequency returns the undamped natural frequency of a complex pole,
// which is its magnitude |p|.
func NaturalFrequency(pole complex128) float64 {
	return cmplx.Abs(pole)
}

// DampingRatio returns the damping ratio of a complex pole, defined as
// -Re(p)/|p|. It returns 0 for a pole at the origin.
func DampingRatio(pole complex128) float64 {
	m := cmplx.Abs(pole)
	if m == 0 {
		return 0
	}
	return -real(pole) / m
}
