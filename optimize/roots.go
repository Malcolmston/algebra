// Package optimize implements one-dimensional root finding and function
// optimization for real-valued functions of a single real variable.
//
// The package gathers the classical numerical methods used to locate the
// zeros and extrema of a scalar function. It is organised into several
// families:
//
//   - Bracketing solvers, which start from an interval [a, b] on which the
//     function changes sign and are guaranteed to converge: Bisection,
//     RegulaFalsi (with the Illinois, Pegasus and Anderson-Björck
//     accelerations), Brent's method (zeroin), Ridders' method and Dekker's
//     method.
//   - Open solvers, which iterate from one or more starting guesses and
//     converge quickly when close to a simple root: Secant, Newton, Halley,
//     Schröder (for multiple roots), Steffensen, inverse quadratic
//     interpolation, Muller's method and fixed-point iteration (optionally
//     accelerated with Aitken's Δ² process).
//   - Polynomial solvers, which return all real (or complex) roots of a
//     polynomial: closed forms for linear, quadratic and cubic equations,
//     and the iterative Durand-Kerner and Aberth-Ehrlich methods together
//     with a companion-matrix constructor.
//   - Minimizers, which locate a local minimum of a unimodal function:
//     golden-section search, Brent's parabolic minimizer, ternary search and
//     Newton/gradient descent on the derivative, plus routines that bracket a
//     minimum.
//   - Supporting numerical utilities: finite-difference derivatives, sign and
//     bracketing predicates, and convergence helpers.
//
// Polynomials are represented throughout as a slice of coefficients in
// ascending order of degree: coeffs[i] is the coefficient of x^i, so that
// coeffs = {c0, c1, c2} denotes c0 + c1*x + c2*x^2.
//
// Every routine is deterministic and depends only on the Go standard library.
package optimize

import (
	"errors"
	"math"
	"math/cmplx"
	"sort"
)

// Func is a real-valued function of a single real variable, the fundamental
// object operated on by the solvers and minimizers in this package.
type Func func(float64) float64

// DefaultTolerance is a reasonable default convergence tolerance for the
// iterative routines when the caller has no specific requirement.
const DefaultTolerance = 1e-10

// DefaultMaxIterations is a reasonable default cap on the number of iterations
// performed by the iterative routines before they give up.
const DefaultMaxIterations = 100

// ErrNoBracket is returned by bracketing solvers when the supplied endpoints
// do not straddle a root, i.e. the function has the same sign at both ends.
var ErrNoBracket = errors.New("optimize: interval does not bracket a root (endpoints share the same sign)")

// ErrMaxIterations is returned when a routine exhausts its iteration budget
// before satisfying the requested tolerance.
var ErrMaxIterations = errors.New("optimize: maximum number of iterations exceeded before convergence")

// ErrZeroDerivative is returned when a method divides by a derivative or a
// finite difference that has collapsed to zero.
var ErrZeroDerivative = errors.New("optimize: derivative (or denominator) evaluated to zero")

// ErrInvalidInterval is returned when an interval or set of starting points is
// degenerate or otherwise unusable.
var ErrInvalidInterval = errors.New("optimize: invalid interval")

// Bracket describes a closed interval [Lo, Hi] believed to enclose a root or a
// minimum.
type Bracket struct {
	Lo float64
	Hi float64
}

// Width returns the (non-negative) length of the bracket.
func (br Bracket) Width() float64 { return math.Abs(br.Hi - br.Lo) }

// Midpoint returns the arithmetic midpoint of the bracket.
func (br Bracket) Midpoint() float64 { return 0.5 * (br.Lo + br.Hi) }

// Contains reports whether x lies within the closed bracket, irrespective of
// whether Lo or Hi is the larger endpoint.
func (br Bracket) Contains(x float64) bool {
	lo, hi := br.Lo, br.Hi
	if lo > hi {
		lo, hi = hi, lo
	}
	return x >= lo && x <= hi
}

// optimizeSameSign reports whether a and b are both strictly positive or both
// strictly negative. A zero operand is treated as a sign change.
func optimizeSameSign(a, b float64) bool {
	return (a > 0 && b > 0) || (a < 0 && b < 0)
}

// optimizeSignf returns the magnitude of a carrying the sign of b, matching the
// FORTRAN/Numerical-Recipes SIGN intrinsic (b >= 0 yields +|a|).
func optimizeSignf(a, b float64) float64 {
	if b >= 0 {
		return math.Abs(a)
	}
	return -math.Abs(a)
}

// optimizeTrimPoly drops trailing (highest-degree) zero coefficients so that
// the returned slice has a non-zero leading coefficient.
func optimizeTrimPoly(c []float64) []float64 {
	n := len(c)
	for n > 0 && c[n-1] == 0 {
		n--
	}
	return c[:n]
}

// optimizePolyEvalC evaluates an ascending-order complex polynomial at x using
// Horner's scheme.
func optimizePolyEvalC(c []complex128, x complex128) complex128 {
	n := len(c)
	if n == 0 {
		return 0
	}
	y := c[n-1]
	for i := n - 2; i >= 0; i-- {
		y = y*x + c[i]
	}
	return y
}

// -----------------------------------------------------------------------------
// Bracketing solvers
// -----------------------------------------------------------------------------

// Bisection finds a root of f in the bracketing interval [a, b] by repeated
// interval halving. It converges linearly and is unconditionally reliable
// whenever f is continuous and f(a) and f(b) have opposite signs. It returns
// ErrNoBracket if the interval does not straddle a root and ErrMaxIterations if
// the width tolerance is not reached within maxIter steps.
func Bisection(f Func, a, b, tol float64, maxIter int) (float64, error) {
	fa := f(a)
	fb := f(b)
	if fa == 0 {
		return a, nil
	}
	if fb == 0 {
		return b, nil
	}
	if optimizeSameSign(fa, fb) {
		return 0, ErrNoBracket
	}
	for i := 0; i < maxIter; i++ {
		m := 0.5 * (a + b)
		fm := f(m)
		if fm == 0 || 0.5*math.Abs(b-a) < tol {
			return m, nil
		}
		if optimizeSameSign(fa, fm) {
			a, fa = m, fm
		} else {
			b, fb = m, fm
		}
	}
	return 0.5 * (a + b), ErrMaxIterations
}

// RegulaFalsi (the method of false position) finds a root of f in [a, b] by
// linear interpolation between the bracketing endpoints. It keeps the root
// bracketed at all times but can converge slowly when one endpoint stagnates;
// see Illinois, Pegasus and AndersonBjorck for accelerated variants.
func RegulaFalsi(f Func, a, b, tol float64, maxIter int) (float64, error) {
	fa := f(a)
	fb := f(b)
	if fa == 0 {
		return a, nil
	}
	if fb == 0 {
		return b, nil
	}
	if optimizeSameSign(fa, fb) {
		return 0, ErrNoBracket
	}
	c := a
	for i := 0; i < maxIter; i++ {
		c = (a*fb - b*fa) / (fb - fa)
		fc := f(c)
		if math.Abs(fc) < tol {
			return c, nil
		}
		if optimizeSameSign(fc, fa) {
			a, fa = c, fc
		} else {
			b, fb = c, fc
		}
	}
	return c, ErrMaxIterations
}

// FalsePosition is an alias for RegulaFalsi, the traditional English name for
// the method of false position.
func FalsePosition(f Func, a, b, tol float64, maxIter int) (float64, error) {
	return RegulaFalsi(f, a, b, tol, maxIter)
}

// Illinois finds a root of f in [a, b] using the Illinois variant of the
// method of false position, which halves the retained (stagnant) endpoint's
// function value to break the one-sided convergence of RegulaFalsi. It keeps
// the root bracketed and converges super-linearly.
func Illinois(f Func, a, b, tol float64, maxIter int) (float64, error) {
	fa := f(a)
	fb := f(b)
	if fa == 0 {
		return a, nil
	}
	if fb == 0 {
		return b, nil
	}
	if optimizeSameSign(fa, fb) {
		return 0, ErrNoBracket
	}
	for i := 0; i < maxIter; i++ {
		c := b - fb*(b-a)/(fb-fa)
		fc := f(c)
		if math.Abs(fc) < tol {
			return c, nil
		}
		if optimizeSameSign(fc, fb) {
			fa *= 0.5
		} else {
			a, fa = b, fb
		}
		b, fb = c, fc
	}
	return b, ErrMaxIterations
}

// Pegasus finds a root of f in [a, b] using the Pegasus variant of the method
// of false position. Like Illinois it damps the stagnant endpoint, but scales
// its function value by fb/(fb+fc), giving faster convergence in practice.
func Pegasus(f Func, a, b, tol float64, maxIter int) (float64, error) {
	fa := f(a)
	fb := f(b)
	if fa == 0 {
		return a, nil
	}
	if fb == 0 {
		return b, nil
	}
	if optimizeSameSign(fa, fb) {
		return 0, ErrNoBracket
	}
	for i := 0; i < maxIter; i++ {
		c := b - fb*(b-a)/(fb-fa)
		fc := f(c)
		if math.Abs(fc) < tol {
			return c, nil
		}
		if optimizeSameSign(fc, fb) {
			if fb+fc != 0 {
				fa *= fb / (fb + fc)
			} else {
				fa *= 0.5
			}
		} else {
			a, fa = b, fb
		}
		b, fb = c, fc
	}
	return b, ErrMaxIterations
}

// AndersonBjorck finds a root of f in [a, b] using the Anderson-Björck variant
// of the method of false position, whose damping factor 1-fc/fb yields
// asymptotic convergence close to that of the secant method while retaining the
// bracketing guarantee.
func AndersonBjorck(f Func, a, b, tol float64, maxIter int) (float64, error) {
	fa := f(a)
	fb := f(b)
	if fa == 0 {
		return a, nil
	}
	if fb == 0 {
		return b, nil
	}
	if optimizeSameSign(fa, fb) {
		return 0, ErrNoBracket
	}
	for i := 0; i < maxIter; i++ {
		c := b - fb*(b-a)/(fb-fa)
		fc := f(c)
		if math.Abs(fc) < tol {
			return c, nil
		}
		if optimizeSameSign(fc, fb) {
			m := 1 - fc/fb
			if m <= 0 {
				m = 0.5
			}
			fa *= m
		} else {
			a, fa = b, fb
		}
		b, fb = c, fc
	}
	return b, ErrMaxIterations
}

// Brent finds a root of f in [a, b] using Brent's method (the classic zeroin
// algorithm), which combines the reliability of bisection with the speed of the
// secant method and inverse quadratic interpolation. It is the recommended
// general-purpose bracketing solver.
func Brent(f Func, a, b, tol float64, maxIter int) (float64, error) {
	fa := f(a)
	fb := f(b)
	if fa == 0 {
		return a, nil
	}
	if fb == 0 {
		return b, nil
	}
	if optimizeSameSign(fa, fb) {
		return 0, ErrNoBracket
	}
	if math.Abs(fa) < math.Abs(fb) {
		a, b = b, a
		fa, fb = fb, fa
	}
	c, fc := a, fa
	d := c
	mflag := true
	for i := 0; i < maxIter; i++ {
		if fb == 0 || math.Abs(b-a) < tol {
			return b, nil
		}
		var s float64
		if fa != fc && fb != fc {
			s = a*fb*fc/((fa-fb)*(fa-fc)) +
				b*fa*fc/((fb-fa)*(fb-fc)) +
				c*fa*fb/((fc-fa)*(fc-fb))
		} else {
			s = b - fb*(b-a)/(fb-fa)
		}
		lo := (3*a + b) / 4
		between := (s > lo && s < b) || (s < lo && s > b)
		if !between ||
			(mflag && math.Abs(s-b) >= math.Abs(b-c)/2) ||
			(!mflag && math.Abs(s-b) >= math.Abs(c-d)/2) ||
			(mflag && math.Abs(b-c) < tol) ||
			(!mflag && math.Abs(c-d) < tol) {
			s = 0.5 * (a + b)
			mflag = true
		} else {
			mflag = false
		}
		fs := f(s)
		d = c
		c, fc = b, fb
		if !optimizeSameSign(fa, fs) {
			b, fb = s, fs
		} else {
			a, fa = s, fs
		}
		if math.Abs(fa) < math.Abs(fb) {
			a, b = b, a
			fa, fb = fb, fa
		}
	}
	return b, ErrMaxIterations
}

// Ridders finds a root of f in [a, b] using Ridders' method, which fits an
// exponential to the two endpoints and their midpoint to obtain a superlinearly
// convergent, always-bracketed iterate. It is robust and requires only two
// function evaluations per step.
func Ridders(f Func, a, b, tol float64, maxIter int) (float64, error) {
	fl := f(a)
	fh := f(b)
	if fl == 0 {
		return a, nil
	}
	if fh == 0 {
		return b, nil
	}
	if optimizeSameSign(fl, fh) {
		return 0, ErrNoBracket
	}
	xl, xh := a, b
	ans := 0.5 * (a + b)
	for i := 0; i < maxIter; i++ {
		xm := 0.5 * (xl + xh)
		fm := f(xm)
		s := math.Sqrt(fm*fm - fl*fh)
		if s == 0 {
			return ans, nil
		}
		xnew := xm + (xm-xl)*(optimizeSignf(1, fl-fh)*fm/s)
		if math.Abs(xnew-ans) <= tol {
			return xnew, nil
		}
		ans = xnew
		fnew := f(ans)
		if fnew == 0 {
			return ans, nil
		}
		if optimizeSignf(fm, fnew) != fm {
			xl, fl = xm, fm
			xh, fh = ans, fnew
		} else if optimizeSignf(fl, fnew) != fl {
			xh, fh = ans, fnew
		} else if optimizeSignf(fh, fnew) != fh {
			xl, fl = ans, fnew
		}
		if math.Abs(xh-xl) <= tol {
			return ans, nil
		}
	}
	return ans, ErrMaxIterations
}

// Dekker finds a root of f in [a, b] using Dekker's method, the secant/bisection
// hybrid that is the historical predecessor of Brent's method. It maintains a
// bracketing contrapoint and falls back to bisection whenever the secant step
// leaves the interval.
func Dekker(f Func, a, b, tol float64, maxIter int) (float64, error) {
	fa := f(a)
	fb := f(b)
	if fa == 0 {
		return a, nil
	}
	if fb == 0 {
		return b, nil
	}
	if optimizeSameSign(fa, fb) {
		return 0, ErrNoBracket
	}
	if math.Abs(fa) < math.Abs(fb) {
		a, b = b, a
		fa, fb = fb, fa
	}
	bprev, fbprev := a, fa
	for i := 0; i < maxIter; i++ {
		if fb == 0 || math.Abs(b-a) < tol {
			return b, nil
		}
		m := 0.5 * (a + b)
		var s float64
		if fb != fbprev {
			s = b - fb*(b-bprev)/(fb-fbprev)
		} else {
			s = m
		}
		if !((s > b && s < m) || (s > m && s < b)) {
			s = m
		}
		bprev, fbprev = b, fb
		b = s
		fb = f(b)
		if optimizeSameSign(fa, fb) {
			a, fa = bprev, fbprev
		}
		if math.Abs(fa) < math.Abs(fb) {
			a, b = b, a
			fa, fb = fb, fa
		}
	}
	return b, ErrMaxIterations
}

// -----------------------------------------------------------------------------
// Open solvers
// -----------------------------------------------------------------------------

// Secant finds a root of f using the secant method, starting from two initial
// guesses x0 and x1. It approximates the derivative by a finite difference of
// successive iterates and converges super-linearly (order ≈ 1.618) near a
// simple root. It returns ErrZeroDerivative if two successive function values
// coincide.
func Secant(f Func, x0, x1, tol float64, maxIter int) (float64, error) {
	f0 := f(x0)
	f1 := f(x1)
	for i := 0; i < maxIter; i++ {
		if f1 == f0 {
			return x1, ErrZeroDerivative
		}
		x2 := x1 - f1*(x1-x0)/(f1-f0)
		if math.Abs(x2-x1) < tol {
			return x2, nil
		}
		x0, f0 = x1, f1
		x1, f1 = x2, f(x2)
	}
	return x1, ErrMaxIterations
}

// Newton finds a root of f using the Newton-Raphson iteration, which requires
// the derivative df and a single starting guess x0. It converges quadratically
// near a simple root. It returns ErrZeroDerivative if df vanishes at an iterate.
func Newton(f, df Func, x0, tol float64, maxIter int) (float64, error) {
	x := x0
	for i := 0; i < maxIter; i++ {
		dfx := df(x)
		if dfx == 0 {
			return x, ErrZeroDerivative
		}
		dx := f(x) / dfx
		x -= dx
		if math.Abs(dx) < tol {
			return x, nil
		}
	}
	return x, ErrMaxIterations
}

// NewtonSafe finds a root of f in the bracketing interval [a, b] using a
// safeguarded Newton iteration (the rtsafe algorithm): it takes a Newton step
// when that step stays inside the current bracket and is decreasing, and falls
// back to bisection otherwise. It combines quadratic convergence with the
// global reliability of bisection.
func NewtonSafe(f, df Func, a, b, tol float64, maxIter int) (float64, error) {
	fl := f(a)
	fh := f(b)
	if fl == 0 {
		return a, nil
	}
	if fh == 0 {
		return b, nil
	}
	if optimizeSameSign(fl, fh) {
		return 0, ErrNoBracket
	}
	var xl, xh float64
	if fl < 0 {
		xl, xh = a, b
	} else {
		xl, xh = b, a
	}
	rts := 0.5 * (a + b)
	dxold := math.Abs(b - a)
	dx := dxold
	fx := f(rts)
	dfx := df(rts)
	for i := 0; i < maxIter; i++ {
		if ((rts-xh)*dfx-fx)*((rts-xl)*dfx-fx) > 0 ||
			math.Abs(2*fx) > math.Abs(dxold*dfx) {
			dxold = dx
			dx = 0.5 * (xh - xl)
			rts = xl + dx
			if xl == rts {
				return rts, nil
			}
		} else {
			dxold = dx
			if dfx == 0 {
				return rts, ErrZeroDerivative
			}
			dx = fx / dfx
			temp := rts
			rts -= dx
			if temp == rts {
				return rts, nil
			}
		}
		if math.Abs(dx) < tol {
			return rts, nil
		}
		fx = f(rts)
		dfx = df(rts)
		if fx < 0 {
			xl = rts
		} else {
			xh = rts
		}
	}
	return rts, ErrMaxIterations
}

// Halley finds a root of f using Halley's method, a third-order iteration that
// uses the first derivative df and second derivative d2f. It converges cubically
// near a simple root, roughly one iteration faster than Newton's method.
func Halley(f, df, d2f Func, x0, tol float64, maxIter int) (float64, error) {
	x := x0
	for i := 0; i < maxIter; i++ {
		fx := f(x)
		dfx := df(x)
		d2fx := d2f(x)
		denom := 2*dfx*dfx - fx*d2fx
		if denom == 0 {
			return x, ErrZeroDerivative
		}
		dx := 2 * fx * dfx / denom
		x -= dx
		if math.Abs(dx) < tol {
			return x, nil
		}
	}
	return x, ErrMaxIterations
}

// Schroeder finds a root of known multiplicity m using the modified Newton (or
// Schröder) iteration x -= m*f/df, which restores quadratic convergence at a
// root where the ordinary Newton method would degrade to linear convergence.
// Passing m = 1 recovers the standard Newton method.
func Schroeder(f, df Func, x0 float64, m int, tol float64, maxIter int) (float64, error) {
	x := x0
	mf := float64(m)
	for i := 0; i < maxIter; i++ {
		fx := f(x)
		dfx := df(x)
		if dfx == 0 {
			// At a root of multiplicity m both f and f' vanish. A zero
			// derivative together with a zero value therefore signals that
			// the iterate already sits on the root; only a zero derivative
			// with a non-zero value is a genuine breakdown.
			if fx == 0 {
				return x, nil
			}
			return x, ErrZeroDerivative
		}
		dx := mf * fx / dfx
		x -= dx
		if math.Abs(dx) < tol {
			return x, nil
		}
	}
	return x, ErrMaxIterations
}

// Steffensen finds a root of f using Steffensen's method, a derivative-free
// iteration that achieves quadratic convergence by estimating the local slope
// from f(x) and f(x+f(x)). It needs only a single starting guess x0.
func Steffensen(f Func, x0, tol float64, maxIter int) (float64, error) {
	x := x0
	for i := 0; i < maxIter; i++ {
		fx := f(x)
		if fx == 0 {
			return x, nil
		}
		denom := f(x+fx) - fx
		if denom == 0 {
			return x, ErrZeroDerivative
		}
		dx := fx * fx / denom
		x -= dx
		if math.Abs(dx) < tol {
			return x, nil
		}
	}
	return x, ErrMaxIterations
}

// InverseQuadratic finds a root of f by inverse quadratic interpolation through
// three starting points x0, x1 and x2: it fits x as a quadratic in f and
// evaluates it at f = 0. Where the three function values are not distinct it
// falls back to a secant step.
func InverseQuadratic(f Func, x0, x1, x2, tol float64, maxIter int) (float64, error) {
	f0, f1, f2 := f(x0), f(x1), f(x2)
	for i := 0; i < maxIter; i++ {
		if f0 == f1 || f0 == f2 || f1 == f2 {
			if f2 == f1 {
				return x2, ErrZeroDerivative
			}
			x3 := x2 - f2*(x2-x1)/(f2-f1)
			if math.Abs(x3-x2) < tol {
				return x3, nil
			}
			x0, f0 = x1, f1
			x1, f1 = x2, f2
			x2, f2 = x3, f(x3)
			continue
		}
		t0 := x0 * f1 * f2 / ((f0 - f1) * (f0 - f2))
		t1 := x1 * f0 * f2 / ((f1 - f0) * (f1 - f2))
		t2 := x2 * f0 * f1 / ((f2 - f0) * (f2 - f1))
		x3 := t0 + t1 + t2
		if math.Abs(x3-x2) < tol {
			return x3, nil
		}
		x0, f0 = x1, f1
		x1, f1 = x2, f2
		x2, f2 = x3, f(x3)
	}
	return x2, ErrMaxIterations
}

// Muller finds a (possibly complex) root of the complex function f using
// Muller's method, which fits a parabola through three starting points and
// takes the root of that parabola nearest the latest iterate. It is well suited
// to polynomials because it can converge to complex roots from real starting
// data.
func Muller(f func(complex128) complex128, x0, x1, x2 complex128, tol float64, maxIter int) (complex128, error) {
	f0 := f(x0)
	f1 := f(x1)
	f2 := f(x2)
	for i := 0; i < maxIter; i++ {
		if x1 == x0 {
			return x2, ErrZeroDerivative
		}
		q := (x2 - x1) / (x1 - x0)
		a := q*f2 - q*(1+q)*f1 + q*q*f0
		b := (2*q+1)*f2 - (1+q)*(1+q)*f1 + q*q*f0
		c := (1 + q) * f2
		disc := cmplx.Sqrt(b*b - 4*a*c)
		denom := b + disc
		if cmplx.Abs(b-disc) > cmplx.Abs(b+disc) {
			denom = b - disc
		}
		if denom == 0 {
			return x2, ErrZeroDerivative
		}
		x3 := x2 - (x2-x1)*2*c/denom
		if cmplx.Abs(x3-x2) < tol {
			return x3, nil
		}
		x0, f0 = x1, f1
		x1, f1 = x2, f2
		x2, f2 = x3, f(x3)
	}
	return x2, ErrMaxIterations
}

// NewtonComplex finds a complex root of the complex function f with derivative
// df using the Newton-Raphson iteration in the complex plane. It converges
// quadratically near a simple root and is the natural tool for tracing the roots
// of complex polynomials.
func NewtonComplex(f, df func(complex128) complex128, x0 complex128, tol float64, maxIter int) (complex128, error) {
	x := x0
	for i := 0; i < maxIter; i++ {
		dfx := df(x)
		if dfx == 0 {
			return x, ErrZeroDerivative
		}
		dx := f(x) / dfx
		x -= dx
		if cmplx.Abs(dx) < tol {
			return x, nil
		}
	}
	return x, ErrMaxIterations
}

// FixedPoint finds a fixed point of g (a value x with g(x) = x) by the direct
// iteration x_{n+1} = g(x_n) starting from x0. It converges linearly when g is a
// contraction near the fixed point.
func FixedPoint(g Func, x0, tol float64, maxIter int) (float64, error) {
	x := x0
	for i := 0; i < maxIter; i++ {
		xn := g(x)
		if math.Abs(xn-x) < tol {
			return xn, nil
		}
		x = xn
	}
	return x, ErrMaxIterations
}

// FixedPointAitken finds a fixed point of g using Aitken's Δ² acceleration
// (Steffensen's fixed-point method): it applies g twice and extrapolates,
// turning linear convergence into quadratic convergence near the fixed point.
func FixedPointAitken(g Func, x0, tol float64, maxIter int) (float64, error) {
	x := x0
	for i := 0; i < maxIter; i++ {
		x1 := g(x)
		x2 := g(x1)
		denom := x2 - 2*x1 + x
		if denom == 0 {
			if math.Abs(x2-x) < tol {
				return x2, nil
			}
			x = x2
			continue
		}
		xa := x - (x1-x)*(x1-x)/denom
		if math.Abs(xa-x) < tol {
			return xa, nil
		}
		x = xa
	}
	return x, ErrMaxIterations
}

// AitkenDelta2 applies a single step of Aitken's Δ² process to three
// consecutive terms x0, x1, x2 of a linearly convergent sequence and returns
// the accelerated estimate of the limit. When the second difference is zero it
// returns x2 unchanged.
func AitkenDelta2(x0, x1, x2 float64) float64 {
	denom := x2 - 2*x1 + x0
	if denom == 0 {
		return x2
	}
	return x0 - (x1-x0)*(x1-x0)/denom
}

// -----------------------------------------------------------------------------
// Bracketing helpers
// -----------------------------------------------------------------------------

// SameSign reports whether a and b are both strictly positive or both strictly
// negative.
func SameSign(a, b float64) bool { return optimizeSameSign(a, b) }

// OppositeSign reports whether a and b have strictly opposite signs.
func OppositeSign(a, b float64) bool {
	return (a > 0 && b < 0) || (a < 0 && b > 0)
}

// SignChange reports whether f takes opposite-signed values at a and b, i.e.
// whether [a, b] is guaranteed (for continuous f) to contain a root.
func SignChange(f Func, a, b float64) bool {
	return !optimizeSameSign(f(a), f(b))
}

// IsBracket reports whether [a, b] strictly brackets a sign change of f.
func IsBracket(f Func, a, b float64) bool {
	return OppositeSign(f(a), f(b))
}

// BracketExpand geometrically expands the interval [a, b] outward (moving the
// endpoint with the larger |f| by a factor of the current width) until f
// changes sign across it or maxIter expansions have been tried. It reports the
// possibly-widened endpoints and whether a bracket was found.
func BracketExpand(f Func, a, b, factor float64, maxIter int) (float64, float64, bool) {
	if a == b {
		return a, b, false
	}
	if factor <= 0 {
		factor = 1.6
	}
	fa := f(a)
	fb := f(b)
	for i := 0; i < maxIter; i++ {
		if !optimizeSameSign(fa, fb) {
			return a, b, true
		}
		if math.Abs(fa) < math.Abs(fb) {
			a += factor * (a - b)
			fa = f(a)
		} else {
			b += factor * (b - a)
			fb = f(b)
		}
	}
	return a, b, !optimizeSameSign(fa, fb)
}

// BracketSubdivide divides [a, b] into n equal sub-intervals and returns every
// sub-interval across which f changes sign, each as a Bracket. It is the tool
// for isolating several roots of a function on a wide interval before refining
// each with a bracketing solver.
func BracketSubdivide(f Func, a, b float64, n int) []Bracket {
	var out []Bracket
	if n < 1 {
		return out
	}
	dx := (b - a) / float64(n)
	x := a
	fp := f(x)
	for i := 0; i < n; i++ {
		xn := x + dx
		fc := f(xn)
		if !optimizeSameSign(fp, fc) {
			out = append(out, Bracket{Lo: x, Hi: xn})
		}
		x, fp = xn, fc
	}
	return out
}

// FindBracket expands the interval [a, b] outward until it encloses a root and
// returns the resulting Bracket together with a flag reporting success.
func FindBracket(f Func, a, b float64) (Bracket, bool) {
	lo, hi, ok := BracketExpand(f, a, b, 1.6, 60)
	return Bracket{Lo: lo, Hi: hi}, ok
}

// -----------------------------------------------------------------------------
// Polynomials
// -----------------------------------------------------------------------------

// PolyEval evaluates the polynomial with the given ascending-order coefficients
// at x using Horner's scheme. An empty slice evaluates to zero.
func PolyEval(coeffs []float64, x float64) float64 {
	n := len(coeffs)
	if n == 0 {
		return 0
	}
	y := coeffs[n-1]
	for i := n - 2; i >= 0; i-- {
		y = y*x + coeffs[i]
	}
	return y
}

// PolyEvalComplex evaluates the real-coefficient polynomial (ascending order) at
// a complex argument x using Horner's scheme.
func PolyEvalComplex(coeffs []float64, x complex128) complex128 {
	n := len(coeffs)
	if n == 0 {
		return 0
	}
	y := complex(coeffs[n-1], 0)
	for i := n - 2; i >= 0; i-- {
		y = y*x + complex(coeffs[i], 0)
	}
	return y
}

// PolyEvalDeriv simultaneously evaluates the polynomial (ascending-order
// coefficients) and its first derivative at x in a single Horner pass, returning
// (value, derivative).
func PolyEvalDeriv(coeffs []float64, x float64) (float64, float64) {
	n := len(coeffs)
	if n == 0 {
		return 0, 0
	}
	p := coeffs[n-1]
	dp := 0.0
	for i := n - 2; i >= 0; i-- {
		dp = dp*x + p
		p = p*x + coeffs[i]
	}
	return p, dp
}

// PolyDerivative returns the coefficients (ascending order) of the derivative of
// the polynomial given by coeffs. A constant polynomial yields an empty slice.
func PolyDerivative(coeffs []float64) []float64 {
	if len(coeffs) <= 1 {
		return []float64{}
	}
	d := make([]float64, len(coeffs)-1)
	for i := 1; i < len(coeffs); i++ {
		d[i-1] = coeffs[i] * float64(i)
	}
	return d
}

// PolyIntegral returns the coefficients (ascending order) of an antiderivative
// of the polynomial given by coeffs, using constant as the value of the
// integration constant (the new degree-zero term).
func PolyIntegral(coeffs []float64, constant float64) []float64 {
	out := make([]float64, len(coeffs)+1)
	out[0] = constant
	for i := 0; i < len(coeffs); i++ {
		out[i+1] = coeffs[i] / float64(i+1)
	}
	return out
}

// PolyDeflate divides the polynomial given by coeffs (ascending order) by the
// linear factor (x - root) using synthetic division. It returns the
// ascending-order quotient coefficients and the scalar remainder, which is zero
// exactly when root is a root of the polynomial.
func PolyDeflate(coeffs []float64, root float64) ([]float64, float64) {
	n := len(coeffs)
	if n == 0 {
		return nil, 0
	}
	if n == 1 {
		return nil, coeffs[0]
	}
	d := n - 1
	b := make([]float64, n)
	b[0] = coeffs[d]
	for k := 1; k <= d; k++ {
		b[k] = coeffs[d-k] + root*b[k-1]
	}
	rem := b[d]
	q := make([]float64, d)
	for i := 0; i < d; i++ {
		q[i] = b[d-1-i]
	}
	return q, rem
}

// CompanionMatrix returns the companion matrix of the polynomial given by
// coeffs (ascending order), a real n×n matrix whose characteristic polynomial
// equals the given polynomial made monic. Its eigenvalues are the roots of the
// polynomial. The matrix is returned in row-major order.
func CompanionMatrix(coeffs []float64) [][]float64 {
	c := optimizeTrimPoly(coeffs)
	n := len(c) - 1
	if n < 1 {
		return nil
	}
	lead := c[n]
	m := make([][]float64, n)
	for i := range m {
		m[i] = make([]float64, n)
	}
	for i := 0; i < n; i++ {
		m[i][n-1] = -c[i] / lead
	}
	for i := 1; i < n; i++ {
		m[i][i-1] = 1
	}
	return m
}

// LinearRoot returns the root of the linear equation a*x + b = 0. The boolean
// result is false when a is zero and no unique root exists.
func LinearRoot(a, b float64) (float64, bool) {
	if a == 0 {
		return 0, false
	}
	return -b / a, true
}

// QuadraticRoots returns both roots of the quadratic a*x^2 + b*x + c = 0 as
// complex numbers, using a numerically stable formulation that avoids
// cancellation. When a is zero it degenerates to the linear (or empty) case.
func QuadraticRoots(a, b, c float64) (complex128, complex128) {
	if a == 0 {
		if b == 0 {
			return cmplx.NaN(), cmplx.NaN()
		}
		r := complex(-c/b, 0)
		return r, r
	}
	disc := b*b - 4*a*c
	if disc >= 0 {
		s := math.Sqrt(disc)
		var q float64
		if b >= 0 {
			q = -0.5 * (b + s)
		} else {
			q = -0.5 * (b - s)
		}
		if q == 0 {
			return complex(0, 0), complex(0, 0)
		}
		return complex(q/a, 0), complex(c/q, 0)
	}
	s := math.Sqrt(-disc)
	re := -b / (2 * a)
	im := s / (2 * a)
	return complex(re, im), complex(re, -im)
}

// QuadraticRealRoots returns the distinct real roots of a*x^2 + b*x + c = 0 in
// ascending order. Complex roots are omitted, so the result has length 0, 1 or
// 2.
func QuadraticRealRoots(a, b, c float64) []float64 {
	if a == 0 {
		if b == 0 {
			return nil
		}
		return []float64{-c / b}
	}
	disc := b*b - 4*a*c
	if disc < 0 {
		return nil
	}
	if disc == 0 {
		return []float64{-b / (2 * a)}
	}
	s := math.Sqrt(disc)
	var q float64
	if b >= 0 {
		q = -0.5 * (b + s)
	} else {
		q = -0.5 * (b - s)
	}
	if q == 0 {
		r := -b / a
		if r < 0 {
			return []float64{r, 0}
		}
		return []float64{0, r}
	}
	r1 := q / a
	r2 := c / q
	if r1 > r2 {
		r1, r2 = r2, r1
	}
	return []float64{r1, r2}
}

// CubicRoots returns all three roots of the cubic a*x^3 + b*x^2 + c*x + d = 0 as
// complex numbers, using Cardano's formula together with the trigonometric
// solution for the three-real-roots case. When a is zero it degenerates to the
// quadratic case.
func CubicRoots(a, b, c, d float64) []complex128 {
	if a == 0 {
		r1, r2 := QuadraticRoots(b, c, d)
		return []complex128{r1, r2}
	}
	bb := b / a
	cc := c / a
	dd := d / a
	p := cc - bb*bb/3
	q := 2*bb*bb*bb/27 - bb*cc/3 + dd
	shift := bb / 3
	disc := q*q/4 + p*p*p/27
	if disc > 0 {
		sq := math.Sqrt(disc)
		u := math.Cbrt(-q/2 + sq)
		v := math.Cbrt(-q/2 - sq)
		x1 := u + v - shift
		re := -(u+v)/2 - shift
		im := (u - v) * math.Sqrt(3) / 2
		return []complex128{complex(x1, 0), complex(re, im), complex(re, -im)}
	}
	if disc == 0 {
		u := math.Cbrt(-q / 2)
		x1 := 2*u - shift
		x2 := -u - shift
		return []complex128{complex(x1, 0), complex(x2, 0), complex(x2, 0)}
	}
	m := 2 * math.Sqrt(-p/3)
	arg := 3 * q / (p * m)
	if arg < -1 {
		arg = -1
	} else if arg > 1 {
		arg = 1
	}
	theta := math.Acos(arg) / 3
	roots := make([]complex128, 3)
	for k := 0; k < 3; k++ {
		t := m * math.Cos(theta-2*math.Pi*float64(k)/3)
		roots[k] = complex(t-shift, 0)
	}
	return roots
}

// CubicRealRoots returns the real roots of a*x^3 + b*x^2 + c*x + d = 0 in
// ascending order (with multiplicity), discarding any genuinely complex roots.
func CubicRealRoots(a, b, c, d float64) []float64 {
	roots := CubicRoots(a, b, c, d)
	var out []float64
	for _, r := range roots {
		if math.Abs(imag(r)) <= 1e-9*(1+math.Abs(real(r))) {
			out = append(out, real(r))
		}
	}
	sort.Float64s(out)
	return out
}

// QuarticRoots returns all four roots of the quartic
// a*x^4 + b*x^3 + c*x^2 + d*x + e = 0 as complex numbers, computed with the
// Durand-Kerner iteration.
func QuarticRoots(a, b, c, d, e float64) []complex128 {
	return DurandKerner([]float64{e, d, c, b, a}, 1e-14, 500)
}

// DurandKerner returns all n roots (real and complex) of the degree-n
// polynomial given by coeffs (ascending order) using the Durand-Kerner
// (Weierstrass) simultaneous iteration. The polynomial is made monic
// internally; roots are refined until the largest update falls below tol or
// maxIter iterations are reached.
func DurandKerner(coeffs []float64, tol float64, maxIter int) []complex128 {
	c := optimizeTrimPoly(coeffs)
	n := len(c) - 1
	if n < 1 {
		return nil
	}
	lead := c[n]
	mon := make([]complex128, n+1)
	for i := 0; i <= n; i++ {
		mon[i] = complex(c[i]/lead, 0)
	}
	roots := make([]complex128, n)
	seed := complex(0.4, 0.9)
	cur := complex(1, 0)
	for i := 0; i < n; i++ {
		cur *= seed
		roots[i] = cur
	}
	for it := 0; it < maxIter; it++ {
		maxDelta := 0.0
		for i := 0; i < n; i++ {
			num := optimizePolyEvalC(mon, roots[i])
			den := complex(1, 0)
			for j := 0; j < n; j++ {
				if j != i {
					den *= roots[i] - roots[j]
				}
			}
			if den == 0 {
				continue
			}
			delta := num / den
			roots[i] -= delta
			if ad := cmplx.Abs(delta); ad > maxDelta {
				maxDelta = ad
			}
		}
		if maxDelta < tol {
			break
		}
	}
	return roots
}

// Aberth returns all n roots (real and complex) of the degree-n polynomial given
// by coeffs (ascending order) using the Aberth-Ehrlich simultaneous iteration,
// which incorporates the derivative for cubic local convergence and typically
// needs fewer iterations than Durand-Kerner. Initial guesses are distributed on
// a circle whose radius bounds the roots.
func Aberth(coeffs []float64, tol float64, maxIter int) []complex128 {
	c := optimizeTrimPoly(coeffs)
	n := len(c) - 1
	if n < 1 {
		return nil
	}
	lead := c[n]
	mon := make([]complex128, n+1)
	for i := 0; i <= n; i++ {
		mon[i] = complex(c[i]/lead, 0)
	}
	dmon := make([]complex128, n)
	for i := 1; i <= n; i++ {
		dmon[i-1] = mon[i] * complex(float64(i), 0)
	}
	radius := 1.0
	for i := 0; i < n; i++ {
		if v := cmplx.Abs(mon[i]) + 1; v > radius {
			radius = v
		}
	}
	roots := make([]complex128, n)
	for i := 0; i < n; i++ {
		theta := 2*math.Pi*float64(i)/float64(n) + 0.5
		roots[i] = complex(radius*math.Cos(theta), radius*math.Sin(theta))
	}
	for it := 0; it < maxIter; it++ {
		maxDelta := 0.0
		for i := 0; i < n; i++ {
			p := optimizePolyEvalC(mon, roots[i])
			dp := optimizePolyEvalC(dmon, roots[i])
			if dp == 0 {
				continue
			}
			ratio := p / dp
			sum := complex(0, 0)
			for j := 0; j < n; j++ {
				if j != i {
					sum += 1 / (roots[i] - roots[j])
				}
			}
			denom := 1 - ratio*sum
			if denom == 0 {
				continue
			}
			offset := ratio / denom
			roots[i] -= offset
			if ad := cmplx.Abs(offset); ad > maxDelta {
				maxDelta = ad
			}
		}
		if maxDelta < tol {
			break
		}
	}
	return roots
}

// PolyComplexRoots returns all roots of the polynomial given by coeffs
// (ascending order) as complex numbers, using Durand-Kerner with default
// settings.
func PolyComplexRoots(coeffs []float64) []complex128 {
	return DurandKerner(coeffs, 1e-14, 1000)
}

// PolyRoots returns the real roots of the polynomial given by coeffs (ascending
// order) in ascending order. A root is considered real when the magnitude of its
// imaginary part is within tol (scaled by the root magnitude) of zero.
func PolyRoots(coeffs []float64, tol float64) []float64 {
	roots := DurandKerner(coeffs, 1e-14, 1000)
	var out []float64
	for _, r := range roots {
		if math.Abs(imag(r)) <= tol*(1+math.Abs(real(r))) {
			out = append(out, real(r))
		}
	}
	sort.Float64s(out)
	return out
}

// -----------------------------------------------------------------------------
// Minimization
// -----------------------------------------------------------------------------

// GoldenSectionMin returns the location of a minimum of the unimodal function f
// on [a, b] using golden-section search, which narrows the bracket by the golden
// ratio at each step and needs a single function evaluation per iteration.
func GoldenSectionMin(f Func, a, b, tol float64, maxIter int) (float64, error) {
	if a > b {
		a, b = b, a
	}
	gr := (math.Sqrt(5) - 1) / 2
	c := b - gr*(b-a)
	d := a + gr*(b-a)
	fc := f(c)
	fd := f(d)
	for i := 0; i < maxIter; i++ {
		if math.Abs(b-a) < tol {
			return 0.5 * (a + b), nil
		}
		if fc < fd {
			b, d, fd = d, c, fc
			c = b - gr*(b-a)
			fc = f(c)
		} else {
			a, c, fc = c, d, fd
			d = a + gr*(b-a)
			fd = f(d)
		}
	}
	return 0.5 * (a + b), ErrMaxIterations
}

// BrentMinimize returns the location and value of a minimum of f on [a, b] using
// Brent's method, which blends golden-section search with parabolic
// interpolation for fast, reliable convergence without requiring derivatives.
// It returns (xmin, f(xmin), error).
func BrentMinimize(f Func, a, b, tol float64, maxIter int) (float64, float64, error) {
	const gold = 0.3819660112501051
	const zeps = 1e-18
	if a > b {
		a, b = b, a
	}
	x := a + gold*(b-a)
	w, v := x, x
	fx := f(x)
	fw, fv := fx, fx
	var d, e float64
	for iter := 0; iter < maxIter; iter++ {
		xm := 0.5 * (a + b)
		tol1 := tol*math.Abs(x) + zeps
		tol2 := 2 * tol1
		if math.Abs(x-xm) <= tol2-0.5*(b-a) {
			return x, fx, nil
		}
		if math.Abs(e) > tol1 {
			r := (x - w) * (fx - fv)
			q := (x - v) * (fx - fw)
			p := (x-v)*q - (x-w)*r
			q = 2 * (q - r)
			if q > 0 {
				p = -p
			}
			q = math.Abs(q)
			etemp := e
			e = d
			if math.Abs(p) >= math.Abs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x) {
				if x >= xm {
					e = a - x
				} else {
					e = b - x
				}
				d = gold * e
			} else {
				d = p / q
				u := x + d
				if u-a < tol2 || b-u < tol2 {
					d = optimizeSignf(tol1, xm-x)
				}
			}
		} else {
			if x >= xm {
				e = a - x
			} else {
				e = b - x
			}
			d = gold * e
		}
		var u float64
		if math.Abs(d) >= tol1 {
			u = x + d
		} else {
			u = x + optimizeSignf(tol1, d)
		}
		fu := f(u)
		if fu <= fx {
			if u >= x {
				a = x
			} else {
				b = x
			}
			v, w, x = w, x, u
			fv, fw, fx = fw, fx, fu
		} else {
			if u < x {
				a = u
			} else {
				b = u
			}
			if fu <= fw || w == x {
				v, w = w, u
				fv, fw = fw, fu
			} else if fu <= fv || v == x || v == w {
				v = u
				fv = fu
			}
		}
	}
	return x, fx, ErrMaxIterations
}

// TernarySearch returns the location of a minimum of the unimodal function f on
// [a, b] using ternary search, which discards one outer third of the interval at
// each step.
func TernarySearch(f Func, a, b, tol float64, maxIter int) (float64, error) {
	if a > b {
		a, b = b, a
	}
	for i := 0; i < maxIter; i++ {
		if b-a < tol {
			return 0.5 * (a + b), nil
		}
		m1 := a + (b-a)/3
		m2 := b - (b-a)/3
		if f(m1) < f(m2) {
			b = m2
		} else {
			a = m1
		}
	}
	return 0.5 * (a + b), ErrMaxIterations
}

// ParabolicMinimum fits a parabola through the three points (x0,f0), (x1,f1) and
// (x2,f2) and returns the abscissa of its vertex. The boolean result is false
// when the points are collinear and no unique vertex exists.
func ParabolicMinimum(x0, x1, x2, f0, f1, f2 float64) (float64, bool) {
	denom := (x1-x0)*(f1-f2) - (x1-x2)*(f1-f0)
	if denom == 0 {
		return 0, false
	}
	num := (x1-x0)*(x1-x0)*(f1-f2) - (x1-x2)*(x1-x2)*(f1-f0)
	return x1 - 0.5*num/denom, true
}

// BracketMinimum searches downhill from the initial points a and b to return a
// triple (a, b, c) with a < b < c (or c < b < a) such that f(b) is less than
// both f(a) and f(c), thereby bracketing a minimum. It implements the classic
// mnbrak algorithm.
func BracketMinimum(f Func, a, b float64) (float64, float64, float64) {
	const gold = 1.618034
	const glimit = 100.0
	const tiny = 1e-20
	fa := f(a)
	fb := f(b)
	if fb > fa {
		a, b = b, a
		fa, fb = fb, fa
	}
	c := b + gold*(b-a)
	fc := f(c)
	for fb > fc {
		r := (b - a) * (fb - fc)
		q := (b - c) * (fb - fa)
		denom := 2 * optimizeSignf(math.Max(math.Abs(q-r), tiny), q-r)
		u := b - ((b-c)*q-(b-a)*r)/denom
		ulim := b + glimit*(c-b)
		var fu float64
		if (b-u)*(u-c) > 0 {
			fu = f(u)
			if fu < fc {
				a, fa = b, fb
				b, fb = u, fu
				return a, b, c
			} else if fu > fb {
				c, fc = u, fu
				return a, b, c
			}
			u = c + gold*(c-b)
			fu = f(u)
		} else if (c-u)*(u-ulim) > 0 {
			fu = f(u)
			if fu < fc {
				b, c = c, u
				u = c + gold*(c-b)
				fb, fc = fc, fu
				fu = f(u)
			}
		} else if (u-ulim)*(ulim-c) >= 0 {
			u = ulim
			fu = f(u)
		} else {
			u = c + gold*(c-b)
			fu = f(u)
		}
		a, fa = b, fb
		b, fb = c, fc
		c, fc = u, fu
	}
	return a, b, c
}

// NewtonMinimize locates a stationary point of a function by applying Newton's
// method to its derivative, using the first derivative df and second derivative
// d2f. Near a non-degenerate minimum it converges quadratically. It returns
// ErrZeroDerivative if the second derivative vanishes.
func NewtonMinimize(df, d2f Func, x0, tol float64, maxIter int) (float64, error) {
	x := x0
	for i := 0; i < maxIter; i++ {
		d2 := d2f(x)
		if d2 == 0 {
			return x, ErrZeroDerivative
		}
		step := df(x) / d2
		x -= step
		if math.Abs(step) < tol {
			return x, nil
		}
	}
	return x, ErrMaxIterations
}

// GradientDescent locates a stationary point of a function by taking fixed-rate
// steps against its derivative df from the starting point x0. It converges for a
// sufficiently small rate and is included as a simple, dependency-free descent
// method.
func GradientDescent(df Func, x0, rate, tol float64, maxIter int) (float64, error) {
	x := x0
	for i := 0; i < maxIter; i++ {
		step := rate * df(x)
		x -= step
		if math.Abs(step) < tol {
			return x, nil
		}
	}
	return x, ErrMaxIterations
}

// -----------------------------------------------------------------------------
// Finite-difference derivatives and utilities
// -----------------------------------------------------------------------------

// Derivative approximates f'(x) by the central difference (f(x+h)-f(x-h))/(2h),
// which has second-order accuracy in h.
func Derivative(f Func, x, h float64) float64 {
	return (f(x+h) - f(x-h)) / (2 * h)
}

// ForwardDifference approximates f'(x) by the forward difference
// (f(x+h)-f(x))/h, a first-order accurate one-sided estimate.
func ForwardDifference(f Func, x, h float64) float64 {
	return (f(x+h) - f(x)) / h
}

// BackwardDifference approximates f'(x) by the backward difference
// (f(x)-f(x-h))/h, a first-order accurate one-sided estimate.
func BackwardDifference(f Func, x, h float64) float64 {
	return (f(x) - f(x-h)) / h
}

// SecondDerivative approximates f”(x) by the central difference
// (f(x+h)-2f(x)+f(x-h))/h^2, which has second-order accuracy in h.
func SecondDerivative(f Func, x, h float64) float64 {
	return (f(x+h) - 2*f(x) + f(x-h)) / (h * h)
}

// RichardsonDerivative approximates f'(x) by Richardson extrapolation of the
// central difference at step sizes h and h/2, cancelling the leading error term
// to yield a fourth-order accurate estimate.
func RichardsonDerivative(f Func, x, h float64) float64 {
	d1 := (f(x+h) - f(x-h)) / (2 * h)
	h2 := h / 2
	d2 := (f(x+h2) - f(x-h2)) / (2 * h2)
	return (4*d2 - d1) / 3
}

// Sign returns +1, -1 or 0 according to the sign of x.
func Sign(x float64) float64 {
	if x > 0 {
		return 1
	}
	if x < 0 {
		return -1
	}
	return 0
}

// Clamp constrains x to the closed interval [lo, hi], swapping the bounds if
// they are given in the wrong order.
func Clamp(x, lo, hi float64) float64 {
	if lo > hi {
		lo, hi = hi, lo
	}
	if x < lo {
		return lo
	}
	if x > hi {
		return hi
	}
	return x
}

// RelativeError returns the relative error of an approximation with respect to an
// exact value. When exact is zero it falls back to the absolute magnitude of the
// approximation.
func RelativeError(approx, exact float64) float64 {
	if exact == 0 {
		return math.Abs(approx)
	}
	return math.Abs((approx - exact) / exact)
}

// AbsoluteError returns the absolute difference |a - b|.
func AbsoluteError(a, b float64) float64 { return math.Abs(a - b) }

// WithinTolerance reports whether a and b differ by at most tol in absolute
// value.
func WithinTolerance(a, b, tol float64) bool { return math.Abs(a-b) <= tol }

// Converged reports whether successive iterates prev and cur satisfy a combined
// relative/absolute convergence test with tolerance tol.
func Converged(prev, cur, tol float64) bool {
	return math.Abs(cur-prev) <= tol*(1+math.Abs(cur))
}
