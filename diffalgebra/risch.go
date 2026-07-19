package diffalgebra

import (
	"math"
	"math/big"
)

// expFloat returns e^x.
func expFloat(x float64) float64 { return math.Exp(x) }

// SolveRischDE solves the Risch differential equation R' + a R = b for a
// polynomial R, given polynomials a and b. It returns the polynomial solution
// and true when one exists, or false when no polynomial solution exists. This
// is the core sub-algorithm used to decide elementary integrability of
// exponential integrands.
func SolveRischDE(a, b Poly) (Poly, bool) {
	if b.IsZero() {
		return ZeroPoly(), true
	}
	if a.IsZero() {
		// R' = b  =>  R = integral of b
		return b.Integral(), true
	}
	da := a.Degree()
	db := b.Degree()
	// Determine a degree bound for R.
	var dR int
	if da >= 1 {
		dR = db - da
	} else {
		dR = db // a is a nonzero constant
	}
	if dR < 0 {
		return ZeroPoly(), false
	}
	// Unknowns: coefficients p_0..p_dR of R. Column j is the polynomial
	// (x^j)' + a*x^j = j x^{j-1} + a x^j.
	cols := make([]Poly, dR+1)
	maxDeg := 0
	for j := 0; j <= dR; j++ {
		mon := Monomial(ratInt(1), j)
		cols[j] = mon.Derivative().Add(a.Mul(mon))
		if cols[j].Degree() > maxDeg {
			maxDeg = cols[j].Degree()
		}
	}
	if b.Degree() > maxDeg {
		maxDeg = b.Degree()
	}
	// Build the linear system: sum_j p_j * cols[j] = b, matched coefficient-wise.
	A := make([][]*big.Rat, maxDeg+1)
	rhs := make([]*big.Rat, maxDeg+1)
	for k := 0; k <= maxDeg; k++ {
		A[k] = make([]*big.Rat, dR+1)
		for j := 0; j <= dR; j++ {
			A[k][j] = cols[j].Coeff(k)
		}
		rhs[k] = b.Coeff(k)
	}
	sol, ok := solveRatSystem(A, rhs)
	if !ok {
		return ZeroPoly(), false
	}
	return NewPoly(sol...), true
}

// RischExpIntegrate decides whether the integrand f(x) e^{g(x)} has an
// elementary antiderivative for polynomials f and g, and if so returns the
// polynomial R with integral = R(x) e^{g(x)}. It solves the Risch differential
// equation R' + g' R = f. The boolean reports elementary integrability.
func RischExpIntegrate(f, g Poly) (Poly, bool) {
	return SolveRischDE(g.Derivative(), f)
}

// ExpIntegrand pairs a polynomial coefficient with an exponential argument,
// representing f(x) e^{g(x)}.
type ExpIntegrand struct {
	Coeff Poly // f(x)
	Arg   Poly // g(x)
}

// Integrate returns the polynomial R such that the antiderivative equals
// R(x) e^{Arg(x)}, together with whether the integrand is elementary.
func (e ExpIntegrand) Integrate() (Poly, bool) {
	return RischExpIntegrate(e.Coeff, e.Arg)
}

// EvalFloat evaluates the elementary antiderivative R(x) e^{g(x)} at x, when it
// exists; the second return reports elementary integrability.
func (e ExpIntegrand) EvalFloat(x float64) (float64, bool) {
	R, ok := e.Integrate()
	if !ok {
		return 0, false
	}
	return R.EvalFloat(x) * expFloat(e.Arg.EvalFloat(x)), true
}

// StructureTheoremExp applies the Risch structure-theorem heuristic to decide
// whether f e^{g} is elementary and returns a human-readable description of the
// antiderivative or the reason it is not elementary.
func StructureTheoremExp(f, g Poly) (string, bool) {
	R, ok := RischExpIntegrate(f, g)
	if !ok {
		return "no elementary antiderivative", false
	}
	return "(" + R.String() + ")*exp(" + g.String() + ")", true
}

// LogarithmicDerivativeIsRational reports whether the logarithmic derivative of
// a rational function f is again a proper rational function whose partial
// fraction has only simple poles, the structural signature of a logarithm. It
// returns the logarithmic derivative and the verdict.
func LogarithmicDerivativeIsRational(f RatFunc) (RatFunc, bool) {
	if f.IsZero() {
		return ZeroRatFunc(), false
	}
	ld, err := f.LogDerivative()
	if err != nil {
		return ZeroRatFunc(), false
	}
	// f'/f has square-free denominator exactly when the pole structure is simple.
	den := ld.Den()
	return ld, den.Equal(den.SquareFreePart())
}

// ratConst is a small helper returning the constant rational polynomial with
// value p/q.
func ratConst(p, q int64) Poly { return ConstPoly(big.NewRat(p, q)) }
