package rootfind

import "math"

// DefaultTol is the default absolute tolerance used by scalar solvers when the
// caller passes a non-positive tolerance.
const DefaultTol = 1e-12

// DefaultMaxIter is the default iteration budget used by scalar solvers when the
// caller passes a non-positive maximum.
const DefaultMaxIter = 200

// Func is a real scalar function of one real variable.
type Func func(x float64) float64

// Result records the outcome of a scalar root-finding iteration.
type Result struct {
	// Root is the best estimate of the root that was found.
	Root float64
	// Value is the function value at Root (the residual).
	Value float64
	// Iterations is the number of iterations actually performed.
	Iterations int
	// Converged reports whether the requested tolerance was reached.
	Converged bool
}

// resolveTol returns tol if positive, else DefaultTol.
func resolveTol(tol float64) float64 {
	if tol > 0 {
		return tol
	}
	return DefaultTol
}

// resolveMax returns maxIter if positive, else DefaultMaxIter.
func resolveMax(maxIter int) int {
	if maxIter > 0 {
		return maxIter
	}
	return DefaultMaxIter
}

// SignChange reports whether f(a) and f(b) have strictly opposite signs, which
// guarantees a root of a continuous f in (a, b) by the intermediate value
// theorem. A zero at an endpoint counts as a sign change.
func SignChange(f Func, a, b float64) bool {
	fa, fb := f(a), f(b)
	if fa == 0 || fb == 0 {
		return true
	}
	return (fa < 0) != (fb < 0)
}

// Bisection finds a root of f in the bracketing interval [a, b] by repeated
// halving. It requires f(a) and f(b) to have opposite signs and returns
// ErrNoBracket otherwise. Bisection is unconditionally convergent and halves
// the bracket each step.
func Bisection(f Func, a, b, tol float64, maxIter int) (Result, error) {
	tol = resolveTol(tol)
	maxIter = resolveMax(maxIter)
	fa, fb := f(a), f(b)
	if fa == 0 {
		return Result{Root: a, Value: 0, Converged: true}, nil
	}
	if fb == 0 {
		return Result{Root: b, Value: 0, Converged: true}, nil
	}
	if (fa < 0) == (fb < 0) {
		return Result{}, ErrNoBracket
	}
	var m, fm float64
	for i := 1; i <= maxIter; i++ {
		m = 0.5 * (a + b)
		fm = f(m)
		if fm == 0 || 0.5*(b-a) <= tol {
			return Result{Root: m, Value: fm, Iterations: i, Converged: true}, nil
		}
		if (fm < 0) == (fa < 0) {
			a, fa = m, fm
		} else {
			b = m
		}
	}
	return Result{Root: m, Value: fm, Iterations: maxIter, Converged: false}, ErrNoConvergence
}

// FalsePosition (regula falsi) finds a root in [a, b] by interpolating a secant
// line through the bracket endpoints and keeping the subinterval that still
// brackets the root. It requires an initial sign change.
func FalsePosition(f Func, a, b, tol float64, maxIter int) (Result, error) {
	tol = resolveTol(tol)
	maxIter = resolveMax(maxIter)
	fa, fb := f(a), f(b)
	if (fa < 0) == (fb < 0) && fa != 0 && fb != 0 {
		return Result{}, ErrNoBracket
	}
	var c, fc float64
	for i := 1; i <= maxIter; i++ {
		c = (a*fb - b*fa) / (fb - fa)
		fc = f(c)
		if math.Abs(fc) <= tol {
			return Result{Root: c, Value: fc, Iterations: i, Converged: true}, nil
		}
		if (fc < 0) == (fa < 0) {
			a, fa = c, fc
		} else {
			b, fb = c, fc
		}
	}
	return Result{Root: c, Value: fc, Iterations: maxIter, Converged: false}, ErrNoConvergence
}

// Illinois is the Illinois variant of the false-position method: when an
// endpoint is retained across consecutive iterations its function value is
// halved, which cures the slow one-sided convergence of plain regula falsi and
// restores superlinear behaviour while keeping the guaranteed bracket.
func Illinois(f Func, a, b, tol float64, maxIter int) (Result, error) {
	tol = resolveTol(tol)
	maxIter = resolveMax(maxIter)
	fa, fb := f(a), f(b)
	if (fa < 0) == (fb < 0) && fa != 0 && fb != 0 {
		return Result{}, ErrNoBracket
	}
	var c, fc float64
	side := 0
	for i := 1; i <= maxIter; i++ {
		c = (a*fb - b*fa) / (fb - fa)
		fc = f(c)
		if math.Abs(fc) <= tol || 0.5*math.Abs(b-a) <= tol {
			return Result{Root: c, Value: fc, Iterations: i, Converged: true}, nil
		}
		if (fc < 0) == (fb < 0) {
			b, fb = c, fc
			if side == -1 {
				fa *= 0.5
			}
			side = -1
		} else {
			a, fa = c, fc
			if side == 1 {
				fb *= 0.5
			}
			side = 1
		}
	}
	return Result{Root: c, Value: fc, Iterations: maxIter, Converged: false}, ErrNoConvergence
}

// Secant finds a root of f from two initial guesses x0 and x1 using the secant
// iteration, which approximates the derivative by a finite difference. It
// converges superlinearly (order ~1.618) but is not guaranteed to bracket.
func Secant(f Func, x0, x1, tol float64, maxIter int) (Result, error) {
	tol = resolveTol(tol)
	maxIter = resolveMax(maxIter)
	f0, f1 := f(x0), f(x1)
	for i := 1; i <= maxIter; i++ {
		denom := f1 - f0
		if denom == 0 {
			return Result{Root: x1, Value: f1, Iterations: i, Converged: false}, ErrZeroDerivative
		}
		x2 := x1 - f1*(x1-x0)/denom
		f2 := f(x2)
		if math.Abs(f2) <= tol || math.Abs(x2-x1) <= tol*(1+math.Abs(x2)) {
			return Result{Root: x2, Value: f2, Iterations: i, Converged: true}, nil
		}
		x0, f0, x1, f1 = x1, f1, x2, f2
	}
	return Result{Root: x1, Value: f1, Iterations: maxIter, Converged: false}, ErrNoConvergence
}

// Newton finds a root of f using Newton's method, given the derivative df. It
// converges quadratically near a simple root but may diverge from a poor start
// or where df vanishes, in which case ErrZeroDerivative is returned.
func Newton(f, df Func, x0, tol float64, maxIter int) (Result, error) {
	tol = resolveTol(tol)
	maxIter = resolveMax(maxIter)
	x := x0
	for i := 1; i <= maxIter; i++ {
		fx := f(x)
		if math.Abs(fx) <= tol {
			return Result{Root: x, Value: fx, Iterations: i, Converged: true}, nil
		}
		d := df(x)
		if d == 0 {
			return Result{Root: x, Value: fx, Iterations: i, Converged: false}, ErrZeroDerivative
		}
		xn := x - fx/d
		if math.Abs(xn-x) <= tol*(1+math.Abs(xn)) {
			return Result{Root: xn, Value: f(xn), Iterations: i, Converged: true}, nil
		}
		x = xn
	}
	return Result{Root: x, Value: f(x), Iterations: maxIter, Converged: false}, ErrNoConvergence
}

// Halley finds a root of f using Halley's third-order method, given the first
// derivative df and second derivative d2f. It converges cubically near a simple
// root, faster than Newton's method at the cost of a second derivative.
func Halley(f, df, d2f Func, x0, tol float64, maxIter int) (Result, error) {
	tol = resolveTol(tol)
	maxIter = resolveMax(maxIter)
	x := x0
	for i := 1; i <= maxIter; i++ {
		fx := f(x)
		if math.Abs(fx) <= tol {
			return Result{Root: x, Value: fx, Iterations: i, Converged: true}, nil
		}
		d1 := df(x)
		d2 := d2f(x)
		denom := 2*d1*d1 - fx*d2
		if denom == 0 {
			return Result{Root: x, Value: fx, Iterations: i, Converged: false}, ErrZeroDerivative
		}
		xn := x - 2*fx*d1/denom
		if math.Abs(xn-x) <= tol*(1+math.Abs(xn)) {
			return Result{Root: xn, Value: f(xn), Iterations: i, Converged: true}, nil
		}
		x = xn
	}
	return Result{Root: x, Value: f(x), Iterations: maxIter, Converged: false}, ErrNoConvergence
}

// Steffensen finds a root of f using Steffensen's method, a derivative-free
// iteration that achieves quadratic convergence by estimating the derivative
// from f(x) and f(x+f(x)). It needs only a single starting point.
func Steffensen(f Func, x0, tol float64, maxIter int) (Result, error) {
	tol = resolveTol(tol)
	maxIter = resolveMax(maxIter)
	x := x0
	for i := 1; i <= maxIter; i++ {
		fx := f(x)
		if math.Abs(fx) <= tol {
			return Result{Root: x, Value: fx, Iterations: i, Converged: true}, nil
		}
		g := f(x+fx) - fx
		if g == 0 {
			return Result{Root: x, Value: fx, Iterations: i, Converged: false}, ErrZeroDerivative
		}
		xn := x - fx*fx/g
		if math.Abs(xn-x) <= tol*(1+math.Abs(xn)) {
			return Result{Root: xn, Value: f(xn), Iterations: i, Converged: true}, nil
		}
		x = xn
	}
	return Result{Root: x, Value: f(x), Iterations: maxIter, Converged: false}, ErrNoConvergence
}

// FixedPoint iterates x <- g(x) from x0 to find a fixed point of g, which is a
// root of f(x) = g(x) - x. Convergence requires |g'| < 1 near the fixed point.
func FixedPoint(g Func, x0, tol float64, maxIter int) (Result, error) {
	tol = resolveTol(tol)
	maxIter = resolveMax(maxIter)
	x := x0
	for i := 1; i <= maxIter; i++ {
		xn := g(x)
		if math.Abs(xn-x) <= tol*(1+math.Abs(xn)) {
			return Result{Root: xn, Value: xn - x, Iterations: i, Converged: true}, nil
		}
		if math.IsNaN(xn) || math.IsInf(xn, 0) {
			return Result{Root: x, Iterations: i, Converged: false}, ErrNoConvergence
		}
		x = xn
	}
	return Result{Root: x, Iterations: maxIter, Converged: false}, ErrNoConvergence
}

// Ridders finds a root of f in the bracket [a, b] using Ridders' method, which
// applies an exponential correction to the false-position estimate and converges
// quadratically while always maintaining a bracket. It requires a sign change.
func Ridders(f Func, a, b, tol float64, maxIter int) (Result, error) {
	tol = resolveTol(tol)
	maxIter = resolveMax(maxIter)
	fa, fb := f(a), f(b)
	if fa == 0 {
		return Result{Root: a, Converged: true}, nil
	}
	if fb == 0 {
		return Result{Root: b, Converged: true}, nil
	}
	if (fa < 0) == (fb < 0) {
		return Result{}, ErrNoBracket
	}
	ans := 0.5 * (a + b)
	for i := 1; i <= maxIter; i++ {
		m := 0.5 * (a + b)
		fm := f(m)
		s := math.Sqrt(fm*fm - fa*fb)
		if s == 0 {
			return Result{Root: ans, Value: f(ans), Iterations: i, Converged: true}, nil
		}
		sign := 1.0
		if fa < fb {
			sign = -1.0
		}
		xnew := m + (m-a)*sign*fm/s
		fnew := f(xnew)
		if math.Abs(fnew) <= tol || math.Abs(xnew-ans) <= tol*(1+math.Abs(xnew)) {
			return Result{Root: xnew, Value: fnew, Iterations: i, Converged: true}, nil
		}
		ans = xnew
		if (fm < 0) != (fnew < 0) {
			a, fa = m, fm
			b, fb = xnew, fnew
		} else if (fa < 0) != (fnew < 0) {
			b, fb = xnew, fnew
		} else {
			a, fa = xnew, fnew
		}
	}
	return Result{Root: ans, Value: f(ans), Iterations: maxIter, Converged: false}, ErrNoConvergence
}

// Brent finds a root of f in [a, b] using Brent's method, which combines
// bisection, the secant method, and inverse quadratic interpolation. It is the
// recommended general-purpose bracketing solver: robust like bisection yet
// usually superlinear. It requires f(a) and f(b) to have opposite signs.
func Brent(f Func, a, b, tol float64, maxIter int) (Result, error) {
	tol = resolveTol(tol)
	maxIter = resolveMax(maxIter)
	fa, fb := f(a), f(b)
	if fa == 0 {
		return Result{Root: a, Converged: true}, nil
	}
	if fb == 0 {
		return Result{Root: b, Converged: true}, nil
	}
	if (fa < 0) == (fb < 0) {
		return Result{}, ErrNoBracket
	}
	if math.Abs(fa) < math.Abs(fb) {
		a, b = b, a
		fa, fb = fb, fa
	}
	c, fc := a, fa
	mflag := true
	var d float64
	var s, fs float64
	for i := 1; i <= maxIter; i++ {
		if fb == 0 || math.Abs(b-a) <= tol {
			return Result{Root: b, Value: fb, Iterations: i, Converged: true}, nil
		}
		if fa != fc && fb != fc {
			// Inverse quadratic interpolation.
			s = a*fb*fc/((fa-fb)*(fa-fc)) +
				b*fa*fc/((fb-fa)*(fb-fc)) +
				c*fa*fb/((fc-fa)*(fc-fb))
		} else {
			// Secant.
			s = b - fb*(b-a)/(fb-fa)
		}
		cond := (s-(3*a+b)/4)*(s-b) >= 0
		if cond ||
			(mflag && math.Abs(s-b) >= math.Abs(b-c)/2) ||
			(!mflag && math.Abs(s-b) >= math.Abs(c-d)/2) ||
			(mflag && math.Abs(b-c) < tol) ||
			(!mflag && math.Abs(c-d) < tol) {
			s = 0.5 * (a + b)
			mflag = true
		} else {
			mflag = false
		}
		fs = f(s)
		d = c
		c, fc = b, fb
		if (fa < 0) != (fs < 0) {
			b, fb = s, fs
		} else {
			a, fa = s, fs
		}
		if math.Abs(fa) < math.Abs(fb) {
			a, b = b, a
			fa, fb = fb, fa
		}
		if math.Abs(fs) <= tol {
			return Result{Root: s, Value: fs, Iterations: i, Converged: true}, nil
		}
	}
	return Result{Root: b, Value: fb, Iterations: maxIter, Converged: false}, ErrNoConvergence
}

// BracketOutward searches outward from the interval [a, b] by geometrically
// expanding it until f changes sign across the endpoints, returning the widened
// bracket. It is a convenience for seeding bracketing solvers. It returns
// ErrNoBracket if no sign change is found within maxIter expansions.
func BracketOutward(f Func, a, b, factor float64, maxIter int) (lo, hi float64, err error) {
	if a == b {
		return 0, 0, ErrBadInput
	}
	if factor <= 1 {
		factor = 1.6
	}
	maxIter = resolveMax(maxIter)
	fa, fb := f(a), f(b)
	for i := 0; i < maxIter; i++ {
		if (fa < 0) != (fb < 0) || fa == 0 || fb == 0 {
			return a, b, nil
		}
		if math.Abs(fa) < math.Abs(fb) {
			a += factor * (a - b)
			fa = f(a)
		} else {
			b += factor * (b - a)
			fb = f(b)
		}
	}
	return a, b, ErrNoBracket
}

// FindBrackets scans the interval [a, b] on a uniform grid of n subintervals and
// returns every subinterval across which f changes sign. It is a simple way to
// locate multiple roots of a function before refining each with a bracketing
// solver.
func FindBrackets(f Func, a, b float64, n int) [][2]float64 {
	if n < 1 {
		n = 1
	}
	var out [][2]float64
	h := (b - a) / float64(n)
	x0 := a
	f0 := f(x0)
	for i := 1; i <= n; i++ {
		x1 := a + float64(i)*h
		f1 := f(x1)
		if f0 == 0 || (f0 < 0) != (f1 < 0) {
			out = append(out, [2]float64{x0, x1})
		}
		x0, f0 = x1, f1
	}
	return out
}
