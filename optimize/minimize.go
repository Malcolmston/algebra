package optimize

import (
	"errors"
	"math"
	"math/rand"
	"sort"
)

// MultiFunc is a scalar-valued objective function of a real vector argument.
// It is the fundamental object minimized by the multivariate routines in this
// file. Implementations must not mutate the slice they are handed.
type MultiFunc func(x []float64) float64

// GradFunc returns the gradient vector (vector of first partial derivatives)
// of an objective evaluated at x. The returned slice has the same length as x.
type GradFunc func(x []float64) []float64

// HessFunc returns the Hessian matrix (matrix of second partial derivatives)
// of an objective evaluated at x, as an n-by-n row-major slice of slices.
type HessFunc func(x []float64) [][]float64

// VectorFunc is a vector-valued function of a vector argument, used by the
// Jacobian utility. It maps an n-vector to an m-vector.
type VectorFunc func(x []float64) []float64

// ScalarFunc is a real-valued function of a single real variable, operated on
// by the one-dimensional minimizers GoldenSection and BrentParabolic.
type ScalarFunc func(x float64) float64

// GoldenRatio is the golden ratio phi = (1 + sqrt(5)) / 2.
const GoldenRatio = 1.6180339887498948482045868343656381

// InvGoldenRatio is the reciprocal of the golden ratio, 1/phi = phi - 1,
// the shrink factor used by golden-section search.
const InvGoldenRatio = 0.6180339887498948482045868343656381

// DefaultStep is the default finite-difference step used by the numeric
// gradient and directional-derivative routines when the caller passes h <= 0.
const DefaultStep = 1e-6

// DefaultLearningRate is a reasonable default step size for the first-order
// descent methods when the caller has no specific requirement.
const DefaultLearningRate = 0.01

// DefaultMaxIter is a reasonable default cap on the number of iterations
// performed by the iterative multivariate routines.
const DefaultMaxIter = 1000

// DefaultTol is a reasonable default convergence tolerance for the
// multivariate routines.
const DefaultTol = 1e-8

// ErrDimensionMismatch is returned when two vectors or a matrix and a vector
// have incompatible lengths.
var ErrDimensionMismatch = errors.New("optimize: vector or matrix dimensions do not match")

// ErrEmptyPoint is returned when a routine is called with an empty starting
// point.
var ErrEmptyPoint = errors.New("optimize: starting point has zero length")

// ErrNonFinite is returned when a computation produces a NaN or infinite value.
var ErrNonFinite = errors.New("optimize: encountered a non-finite value")

// ErrSingularMatrix is returned by SolveLinearSystem and the Newton solver when
// the coefficient matrix is (numerically) singular.
var ErrSingularMatrix = errors.New("optimize: matrix is singular")

// Result reports the outcome of a multivariate minimization: the located point
// X, the objective value F at that point, the number of Iterations performed,
// and whether the routine Converged to the requested tolerance.
type Result struct {
	X          []float64
	F          float64
	Iterations int
	Converged  bool
}

// ScalarResult reports the outcome of a one-dimensional minimization: the
// located abscissa X, the objective value F there, the number of Iterations
// performed, and whether the routine Converged.
type ScalarResult struct {
	X          float64
	F          float64
	Iterations int
	Converged  bool
}

// Options bundles the common tuning parameters shared by the iterative
// multivariate minimizers. A zero Options is not usable directly; obtain a
// populated value from DefaultOptions and adjust the fields you care about.
type Options struct {
	// Tol is the convergence tolerance on the gradient norm (or on the
	// simplex/step size for derivative-free methods).
	Tol float64
	// MaxIter is the maximum number of iterations.
	MaxIter int
	// Step is the finite-difference step for numeric derivatives.
	Step float64
	// LearningRate is the base step size for first-order methods.
	LearningRate float64
	// Momentum is the momentum coefficient for the momentum methods.
	Momentum float64
}

// DefaultOptions returns an Options value populated with sensible defaults.
func DefaultOptions() Options {
	return Options{
		Tol:          DefaultTol,
		MaxIter:      DefaultMaxIter,
		Step:         DefaultStep,
		LearningRate: DefaultLearningRate,
		Momentum:     0.9,
	}
}

// SAOptions configures the SimulatedAnnealing minimizer.
type SAOptions struct {
	// InitialTemp is the starting temperature; higher values accept more
	// uphill moves early on.
	InitialTemp float64
	// CoolingRate is the geometric cooling factor in (0, 1); the temperature
	// is multiplied by it each iteration.
	CoolingRate float64
	// MinTemp is a floor below which the temperature is not allowed to fall.
	MinTemp float64
	// StepSize scales the Gaussian perturbation used to propose neighbours.
	StepSize float64
	// MaxIter is the number of annealing steps to perform.
	MaxIter int
}

// DefaultSAOptions returns an SAOptions value populated with sensible defaults.
func DefaultSAOptions() SAOptions {
	return SAOptions{
		InitialTemp: 10.0,
		CoolingRate: 0.995,
		MinTemp:     1e-8,
		StepSize:    1.0,
		MaxIter:     10000,
	}
}

// ---------------------------------------------------------------------------
// Vector utilities
// ---------------------------------------------------------------------------

// VecCopy returns a fresh copy of x.
func VecCopy(x []float64) []float64 {
	out := make([]float64, len(x))
	copy(out, x)
	return out
}

// VecZeros returns a new zero vector of length n.
func VecZeros(n int) []float64 { return make([]float64, n) }

// VecFill returns a new vector of length n with every element set to v.
func VecFill(n int, v float64) []float64 {
	out := make([]float64, n)
	for i := range out {
		out[i] = v
	}
	return out
}

// VecAdd returns the element-wise sum a + b. It panics if the lengths differ.
func VecAdd(a, b []float64) []float64 {
	optimizemustMatch(len(a), len(b))
	out := make([]float64, len(a))
	for i := range a {
		out[i] = a[i] + b[i]
	}
	return out
}

// VecSub returns the element-wise difference a - b. It panics if the lengths
// differ.
func VecSub(a, b []float64) []float64 {
	optimizemustMatch(len(a), len(b))
	out := make([]float64, len(a))
	for i := range a {
		out[i] = a[i] - b[i]
	}
	return out
}

// VecScale returns the vector s*a.
func VecScale(a []float64, s float64) []float64 {
	out := make([]float64, len(a))
	for i := range a {
		out[i] = s * a[i]
	}
	return out
}

// VecNegate returns -a.
func VecNegate(a []float64) []float64 { return VecScale(a, -1) }

// VecAxpy returns the combination s*x + y (the classic "a x plus y"). It panics
// if the lengths differ.
func VecAxpy(s float64, x, y []float64) []float64 {
	optimizemustMatch(len(x), len(y))
	out := make([]float64, len(x))
	for i := range x {
		out[i] = s*x[i] + y[i]
	}
	return out
}

// VecLinComb returns the linear combination a*x + b*y. It panics if the lengths
// differ.
func VecLinComb(a float64, x []float64, b float64, y []float64) []float64 {
	optimizemustMatch(len(x), len(y))
	out := make([]float64, len(x))
	for i := range x {
		out[i] = a*x[i] + b*y[i]
	}
	return out
}

// VecDot returns the Euclidean inner product a . b. It panics if the lengths
// differ.
func VecDot(a, b []float64) float64 {
	optimizemustMatch(len(a), len(b))
	var s float64
	for i := range a {
		s += a[i] * b[i]
	}
	return s
}

// VecNormSquared returns the squared Euclidean norm |a|^2.
func VecNormSquared(a []float64) float64 { return VecDot(a, a) }

// VecNorm returns the Euclidean (L2) norm |a|.
func VecNorm(a []float64) float64 { return math.Sqrt(VecDot(a, a)) }

// VecInfNorm returns the maximum-absolute-value (L-infinity) norm of a.
func VecInfNorm(a []float64) float64 {
	var m float64
	for _, v := range a {
		if av := math.Abs(v); av > m {
			m = av
		}
	}
	return m
}

// VecDistance returns the Euclidean distance between a and b.
func VecDistance(a, b []float64) float64 { return VecNorm(VecSub(a, b)) }

// VecClamp returns a copy of x with every element confined to the scalar range
// [lo, hi].
func VecClamp(x []float64, lo, hi float64) []float64 {
	out := make([]float64, len(x))
	for i, v := range x {
		out[i] = optimizeclamp(v, lo, hi)
	}
	return out
}

// ProjectBox returns a copy of x projected onto the axis-aligned box defined by
// the per-coordinate bounds lo and hi. It panics if the lengths differ.
func ProjectBox(x, lo, hi []float64) []float64 {
	optimizemustMatch(len(x), len(lo))
	optimizemustMatch(len(x), len(hi))
	out := make([]float64, len(x))
	for i := range x {
		out[i] = optimizeclamp(x[i], lo[i], hi[i])
	}
	return out
}

// Centroid returns the arithmetic mean of a set of points, all assumed to have
// the same dimension. It returns nil for an empty set.
func Centroid(points [][]float64) []float64 {
	if len(points) == 0 {
		return nil
	}
	n := len(points[0])
	c := make([]float64, n)
	for _, p := range points {
		for i := 0; i < n; i++ {
			c[i] += p[i]
		}
	}
	inv := 1 / float64(len(points))
	for i := range c {
		c[i] *= inv
	}
	return c
}

// ---------------------------------------------------------------------------
// Matrix utilities
// ---------------------------------------------------------------------------

// MatIdentity returns the n-by-n identity matrix.
func MatIdentity(n int) [][]float64 {
	m := make([][]float64, n)
	for i := 0; i < n; i++ {
		m[i] = make([]float64, n)
		m[i][i] = 1
	}
	return m
}

// MatCopy returns a deep copy of the matrix a.
func MatCopy(a [][]float64) [][]float64 {
	m := make([][]float64, len(a))
	for i := range a {
		m[i] = make([]float64, len(a[i]))
		copy(m[i], a[i])
	}
	return m
}

// MatTranspose returns the transpose of the (rectangular) matrix a.
func MatTranspose(a [][]float64) [][]float64 {
	if len(a) == 0 {
		return [][]float64{}
	}
	rows, cols := len(a), len(a[0])
	t := make([][]float64, cols)
	for j := 0; j < cols; j++ {
		t[j] = make([]float64, rows)
		for i := 0; i < rows; i++ {
			t[j][i] = a[i][j]
		}
	}
	return t
}

// MatVec returns the matrix-vector product a*x. It panics if the column count
// of a does not equal len(x).
func MatVec(a [][]float64, x []float64) []float64 {
	out := make([]float64, len(a))
	for i := range a {
		optimizemustMatch(len(a[i]), len(x))
		var s float64
		for j := range x {
			s += a[i][j] * x[j]
		}
		out[i] = s
	}
	return out
}

// OuterProduct returns the outer product a*b^T, an len(a)-by-len(b) matrix
// whose (i, j) entry is a[i]*b[j].
func OuterProduct(a, b []float64) [][]float64 {
	m := make([][]float64, len(a))
	for i := range a {
		m[i] = make([]float64, len(b))
		for j := range b {
			m[i][j] = a[i] * b[j]
		}
	}
	return m
}

// SolveLinearSystem solves the linear system a*x = b for x using Gaussian
// elimination with partial pivoting. The inputs are left unmodified. It returns
// ErrSingularMatrix if the matrix is numerically singular and
// ErrDimensionMismatch if the shapes are inconsistent.
func SolveLinearSystem(a [][]float64, b []float64) ([]float64, error) {
	n := len(a)
	if n == 0 {
		return []float64{}, nil
	}
	if len(b) != n {
		return nil, ErrDimensionMismatch
	}
	m := MatCopy(a)
	rhs := VecCopy(b)
	for k := 0; k < n; k++ {
		if len(m[k]) != n {
			return nil, ErrDimensionMismatch
		}
		// Partial pivot: find the largest magnitude entry in column k.
		p := k
		best := math.Abs(m[k][k])
		for i := k + 1; i < n; i++ {
			if v := math.Abs(m[i][k]); v > best {
				best, p = v, i
			}
		}
		if best == 0 {
			return nil, ErrSingularMatrix
		}
		if p != k {
			m[k], m[p] = m[p], m[k]
			rhs[k], rhs[p] = rhs[p], rhs[k]
		}
		for i := k + 1; i < n; i++ {
			factor := m[i][k] / m[k][k]
			for j := k; j < n; j++ {
				m[i][j] -= factor * m[k][j]
			}
			rhs[i] -= factor * rhs[k]
		}
	}
	sol := make([]float64, n)
	for i := n - 1; i >= 0; i-- {
		s := rhs[i]
		for j := i + 1; j < n; j++ {
			s -= m[i][j] * sol[j]
		}
		sol[i] = s / m[i][i]
	}
	return sol, nil
}

// ---------------------------------------------------------------------------
// Numeric derivatives
// ---------------------------------------------------------------------------

// PartialDerivative approximates the partial derivative of f with respect to
// coordinate i at x by a central finite difference with step h. If h <= 0 the
// DefaultStep is used.
func PartialDerivative(f MultiFunc, x []float64, i int, h float64) float64 {
	if h <= 0 {
		h = DefaultStep
	}
	xp := VecCopy(x)
	orig := xp[i]
	xp[i] = orig + h
	fp := f(xp)
	xp[i] = orig - h
	fm := f(xp)
	return (fp - fm) / (2 * h)
}

// NumericGradientCentral approximates the gradient of f at x using central
// finite differences with step h. If h <= 0 the DefaultStep is used. This is
// the most accurate of the finite-difference gradient estimators.
func NumericGradientCentral(f MultiFunc, x []float64, h float64) []float64 {
	if h <= 0 {
		h = DefaultStep
	}
	n := len(x)
	g := make([]float64, n)
	xp := VecCopy(x)
	for i := 0; i < n; i++ {
		orig := xp[i]
		xp[i] = orig + h
		fp := f(xp)
		xp[i] = orig - h
		fm := f(xp)
		xp[i] = orig
		g[i] = (fp - fm) / (2 * h)
	}
	return g
}

// NumericGradientForward approximates the gradient of f at x using forward
// finite differences with step h. If h <= 0 the DefaultStep is used. It uses
// one fewer evaluation per coordinate than the central estimator at the cost of
// accuracy.
func NumericGradientForward(f MultiFunc, x []float64, h float64) []float64 {
	if h <= 0 {
		h = DefaultStep
	}
	n := len(x)
	g := make([]float64, n)
	f0 := f(x)
	xp := VecCopy(x)
	for i := 0; i < n; i++ {
		orig := xp[i]
		xp[i] = orig + h
		g[i] = (f(xp) - f0) / h
		xp[i] = orig
	}
	return g
}

// NumericGradient approximates the gradient of f at x using central finite
// differences. It is a convenience alias for NumericGradientCentral with the
// DefaultStep.
func NumericGradient(f MultiFunc, x []float64) []float64 {
	return NumericGradientCentral(f, x, DefaultStep)
}

// NumericHessian approximates the Hessian of f at x using central second
// differences with step h. If h <= 0 a step of 1e-4 is used. The returned
// matrix is symmetric by construction.
func NumericHessian(f MultiFunc, x []float64, h float64) [][]float64 {
	if h <= 0 {
		h = 1e-4
	}
	n := len(x)
	hess := make([][]float64, n)
	for i := range hess {
		hess[i] = make([]float64, n)
	}
	xp := VecCopy(x)
	f0 := f(xp)
	for i := 0; i < n; i++ {
		for j := i; j < n; j++ {
			if i == j {
				// Diagonal second derivative: the mixed four-point stencil
				// degenerates when both perturbations act on the same
				// coordinate, so use the standard central formula
				// (f(x+h) - 2 f(x) + f(x-h)) / h².
				oi := xp[i]
				xp[i] = oi + h
				fp := f(xp)
				xp[i] = oi - h
				fm := f(xp)
				xp[i] = oi
				hess[i][i] = (fp - 2*f0 + fm) / (h * h)
				continue
			}
			oi, oj := xp[i], xp[j]
			xp[i], xp[j] = oi+h, oj+h
			fpp := f(xp)
			xp[i], xp[j] = oi+h, oj-h
			fpm := f(xp)
			xp[i], xp[j] = oi-h, oj+h
			fmp := f(xp)
			xp[i], xp[j] = oi-h, oj-h
			fmm := f(xp)
			xp[i], xp[j] = oi, oj
			v := (fpp - fpm - fmp + fmm) / (4 * h * h)
			hess[i][j] = v
			hess[j][i] = v
		}
	}
	return hess
}

// NumericJacobian approximates the Jacobian of the vector-valued function g at
// x using central finite differences with step h. If h <= 0 the DefaultStep is
// used. The result is an m-by-n matrix whose (i, j) entry is d g_i / d x_j.
func NumericJacobian(g VectorFunc, x []float64, h float64) [][]float64 {
	if h <= 0 {
		h = DefaultStep
	}
	n := len(x)
	xp := VecCopy(x)
	var cols [][]float64
	for j := 0; j < n; j++ {
		orig := xp[j]
		xp[j] = orig + h
		fp := g(xp)
		xp[j] = orig - h
		fm := g(xp)
		xp[j] = orig
		col := make([]float64, len(fp))
		for i := range fp {
			col[i] = (fp[i] - fm[i]) / (2 * h)
		}
		cols = append(cols, col)
	}
	// Assemble m-by-n from the n columns.
	m := 0
	if len(cols) > 0 {
		m = len(cols[0])
	}
	jac := make([][]float64, m)
	for i := 0; i < m; i++ {
		jac[i] = make([]float64, n)
		for j := 0; j < n; j++ {
			jac[i][j] = cols[j][i]
		}
	}
	return jac
}

// DirectionalDerivative approximates the derivative of f at x along the unit-
// scaled direction dir using a central finite difference with step h. If h <= 0
// the DefaultStep is used. The direction need not be normalized; the returned
// value is grad(f) . dir.
func DirectionalDerivative(f MultiFunc, x, dir []float64, h float64) float64 {
	if h <= 0 {
		h = DefaultStep
	}
	fp := f(VecAxpy(h, dir, x))
	fm := f(VecAxpy(-h, dir, x))
	return (fp - fm) / (2 * h)
}

// ---------------------------------------------------------------------------
// One-dimensional minimizers
// ---------------------------------------------------------------------------

// GoldenSection locates a minimum of the unimodal function f on the interval
// [a, b] by golden-section search, which shrinks the bracket by the constant
// factor 1/phi each step. It stops when the bracket width falls below tol or
// after maxIter iterations.
func GoldenSection(f ScalarFunc, a, b, tol float64, maxIter int) ScalarResult {
	if a > b {
		a, b = b, a
	}
	if tol <= 0 {
		tol = DefaultTol
	}
	if maxIter <= 0 {
		maxIter = DefaultMaxIter
	}
	c := b - InvGoldenRatio*(b-a)
	d := a + InvGoldenRatio*(b-a)
	fc, fd := f(c), f(d)
	iter := 0
	for ; iter < maxIter; iter++ {
		if math.Abs(b-a) <= tol {
			break
		}
		if fc < fd {
			b, d, fd = d, c, fc
			c = b - InvGoldenRatio*(b-a)
			fc = f(c)
		} else {
			a, c, fc = c, d, fd
			d = a + InvGoldenRatio*(b-a)
			fd = f(d)
		}
	}
	xm := 0.5 * (a + b)
	return ScalarResult{X: xm, F: f(xm), Iterations: iter, Converged: math.Abs(b-a) <= tol}
}

// BrentParabolic locates a minimum of the function f on the interval [a, b]
// using Brent's method, which combines the reliability of golden-section search
// with the fast convergence of successive parabolic interpolation. It stops
// when the estimate is bracketed to within tol or after maxIter iterations.
func BrentParabolic(f ScalarFunc, a, b, tol float64, maxIter int) ScalarResult {
	const cgold = 0.3819660112501051
	const zeps = 1e-12
	if a > b {
		a, b = b, a
	}
	if tol <= 0 {
		tol = DefaultTol
	}
	if maxIter <= 0 {
		maxIter = DefaultMaxIter
	}
	var d, e float64
	x := a + cgold*(b-a)
	w, v := x, x
	fx := f(x)
	fw, fv := fx, fx
	iter := 0
	converged := false
	for ; iter < maxIter; iter++ {
		xm := 0.5 * (a + b)
		tol1 := tol*math.Abs(x) + zeps
		tol2 := 2 * tol1
		if math.Abs(x-xm) <= (tol2 - 0.5*(b-a)) {
			converged = true
			break
		}
		useGolden := true
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
			if math.Abs(p) < math.Abs(0.5*q*etemp) && p > q*(a-x) && p < q*(b-x) {
				d = p / q
				u := x + d
				if (u-a) < tol2 || (b-u) < tol2 {
					d = math.Copysign(tol1, xm-x)
				}
				useGolden = false
			}
		}
		if useGolden {
			if x >= xm {
				e = a - x
			} else {
				e = b - x
			}
			d = cgold * e
		}
		var u float64
		if math.Abs(d) >= tol1 {
			u = x + d
		} else {
			u = x + math.Copysign(tol1, d)
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
	return ScalarResult{X: x, F: fx, Iterations: iter, Converged: converged}
}

// MinimizeScalar locates a minimum of f on [a, b] using BrentParabolic with the
// DefaultTol and DefaultMaxIter. It is a convenience wrapper for the common
// case where the caller does not need to tune the stopping rule.
func MinimizeScalar(f ScalarFunc, a, b float64) ScalarResult {
	return BrentParabolic(f, a, b, DefaultTol, DefaultMaxIter)
}

// ---------------------------------------------------------------------------
// Line search
// ---------------------------------------------------------------------------

// BacktrackingLineSearch returns a step length alpha along the descent
// direction dir from x that satisfies the Armijo (sufficient decrease)
// condition f(x + alpha*dir) <= f(x) + c*alpha*(grad . dir). Starting from
// alpha0 it repeatedly multiplies the step by rho (0 < rho < 1) until the
// condition holds or maxIter reductions have been made. grad must be the
// gradient of f at x and dir must be a descent direction (grad . dir < 0) for
// the guarantee to hold.
func BacktrackingLineSearch(f MultiFunc, x, dir, grad []float64, alpha0, c, rho float64, maxIter int) float64 {
	if alpha0 <= 0 {
		alpha0 = 1
	}
	if c <= 0 {
		c = 1e-4
	}
	if rho <= 0 || rho >= 1 {
		rho = 0.5
	}
	if maxIter <= 0 {
		maxIter = 50
	}
	fx := f(x)
	gd := VecDot(grad, dir)
	alpha := alpha0
	for i := 0; i < maxIter; i++ {
		if f(VecAxpy(alpha, dir, x)) <= fx+c*alpha*gd {
			break
		}
		alpha *= rho
	}
	return alpha
}

// ArmijoStep is a convenience wrapper around BacktrackingLineSearch using the
// conventional parameters alpha0 = 1, c = 1e-4, rho = 0.5 and 50 reductions.
func ArmijoStep(f MultiFunc, x, dir, grad []float64) float64 {
	return BacktrackingLineSearch(f, x, dir, grad, 1, 1e-4, 0.5, 50)
}

// ---------------------------------------------------------------------------
// First-order multivariate minimizers
// ---------------------------------------------------------------------------

// GradientDescentMomentum minimizes f by gradient descent with classical
// (heavy-ball) momentum starting from x0. The velocity is updated as
// v <- momentum*v - rate*grad and the point as x <- x + v. If grad is nil the
// gradient is estimated by central finite differences. Iteration stops when the
// gradient norm falls below tol or after maxIter steps.
func GradientDescentMomentum(f MultiFunc, grad GradFunc, x0 []float64, rate, momentum, tol float64, maxIter int) Result {
	g := optimizegrad(f, grad)
	if rate <= 0 {
		rate = DefaultLearningRate
	}
	if tol <= 0 {
		tol = DefaultTol
	}
	if maxIter <= 0 {
		maxIter = DefaultMaxIter
	}
	x := VecCopy(x0)
	v := make([]float64, len(x))
	iter := 0
	converged := false
	for ; iter < maxIter; iter++ {
		gr := g(x)
		if VecNorm(gr) < tol {
			converged = true
			break
		}
		for i := range x {
			v[i] = momentum*v[i] - rate*gr[i]
			x[i] += v[i]
		}
	}
	return Result{X: x, F: f(x), Iterations: iter, Converged: converged}
}

// GradientDescentNesterov minimizes f by gradient descent with Nesterov's
// accelerated (look-ahead) momentum starting from x0. The gradient is evaluated
// at the look-ahead point x + momentum*v. If grad is nil it is estimated by
// central finite differences. Iteration stops when the gradient norm at x falls
// below tol or after maxIter steps.
func GradientDescentNesterov(f MultiFunc, grad GradFunc, x0 []float64, rate, momentum, tol float64, maxIter int) Result {
	g := optimizegrad(f, grad)
	if rate <= 0 {
		rate = DefaultLearningRate
	}
	if tol <= 0 {
		tol = DefaultTol
	}
	if maxIter <= 0 {
		maxIter = DefaultMaxIter
	}
	x := VecCopy(x0)
	v := make([]float64, len(x))
	iter := 0
	converged := false
	for ; iter < maxIter; iter++ {
		if VecNorm(g(x)) < tol {
			converged = true
			break
		}
		look := make([]float64, len(x))
		for i := range x {
			look[i] = x[i] + momentum*v[i]
		}
		gr := g(look)
		for i := range x {
			v[i] = momentum*v[i] - rate*gr[i]
			x[i] += v[i]
		}
	}
	return Result{X: x, F: f(x), Iterations: iter, Converged: converged}
}

// CoordinateDescent minimizes f by cyclically minimizing along each coordinate
// axis in turn. Each coordinate sub-problem is solved by BrentParabolic on the
// window [x_i - window, x_i + window] centred on the current value. Iteration
// stops when a full sweep moves the point by less than tol or after maxIter
// sweeps. It is well suited to separable or mildly coupled convex objectives.
func CoordinateDescent(f MultiFunc, x0 []float64, window, tol float64, maxIter int) Result {
	if window <= 0 {
		window = 1
	}
	if tol <= 0 {
		tol = DefaultTol
	}
	if maxIter <= 0 {
		maxIter = DefaultMaxIter
	}
	x := VecCopy(x0)
	iter := 0
	converged := false
	for ; iter < maxIter; iter++ {
		prev := VecCopy(x)
		for i := range x {
			idx := i
			line := func(t float64) float64 {
				saved := x[idx]
				x[idx] = t
				val := f(x)
				x[idx] = saved
				return val
			}
			res := BrentParabolic(ScalarFunc(line), x[idx]-window, x[idx]+window, tol, 200)
			x[idx] = res.X
		}
		if VecDistance(x, prev) < tol {
			converged = true
			break
		}
	}
	return Result{X: x, F: f(x), Iterations: iter, Converged: converged}
}

// ConjugateGradient minimizes f by the nonlinear conjugate-gradient method with
// Fletcher-Reeves updates and Armijo backtracking line searches, starting from
// x0. If grad is nil the gradient is estimated by central finite differences.
// The search direction is reset to steepest descent every len(x0) iterations to
// preserve convergence. Iteration stops when the gradient norm falls below tol
// or after maxIter steps.
func ConjugateGradient(f MultiFunc, grad GradFunc, x0 []float64, tol float64, maxIter int) Result {
	g := optimizegrad(f, grad)
	if tol <= 0 {
		tol = DefaultTol
	}
	if maxIter <= 0 {
		maxIter = DefaultMaxIter
	}
	n := len(x0)
	x := VecCopy(x0)
	gr := g(x)
	dir := VecNegate(gr)
	iter := 0
	converged := false
	for ; iter < maxIter; iter++ {
		if VecNorm(gr) < tol {
			converged = true
			break
		}
		alpha := BacktrackingLineSearch(f, x, dir, gr, 1, 1e-4, 0.5, 60)
		x = VecAxpy(alpha, dir, x)
		grNew := g(x)
		denom := VecDot(gr, gr)
		var beta float64
		if denom > 0 {
			beta = VecDot(grNew, grNew) / denom
		}
		if (iter+1)%n == 0 {
			beta = 0 // periodic restart
		}
		dir = VecLinComb(-1, grNew, beta, dir)
		gr = grNew
	}
	return Result{X: x, F: f(x), Iterations: iter, Converged: converged}
}

// ---------------------------------------------------------------------------
// Second-order and quasi-Newton minimizers
// ---------------------------------------------------------------------------

// NewtonMultivariate minimizes f by the damped Newton method starting from x0.
// At each step it solves H*p = -g for the Newton direction p (where g and H are
// the gradient and Hessian) and takes an Armijo-backtracked step along p. If
// grad or hess is nil the corresponding quantity is estimated by finite
// differences. If the Hessian is singular the method falls back to a steepest-
// descent step. Iteration stops when the gradient norm falls below tol or after
// maxIter steps.
func NewtonMultivariate(f MultiFunc, grad GradFunc, hess HessFunc, x0 []float64, tol float64, maxIter int) Result {
	g := optimizegrad(f, grad)
	h := hess
	if h == nil {
		h = func(x []float64) [][]float64 { return NumericHessian(f, x, 1e-4) }
	}
	if tol <= 0 {
		tol = DefaultTol
	}
	if maxIter <= 0 {
		maxIter = DefaultMaxIter
	}
	x := VecCopy(x0)
	iter := 0
	converged := false
	for ; iter < maxIter; iter++ {
		gr := g(x)
		if VecNorm(gr) < tol {
			converged = true
			break
		}
		var dir []float64
		p, err := SolveLinearSystem(h(x), VecNegate(gr))
		if err != nil || VecDot(p, gr) >= 0 {
			// Not a descent direction (or singular): fall back to gradient.
			dir = VecNegate(gr)
		} else {
			dir = p
		}
		alpha := BacktrackingLineSearch(f, x, dir, gr, 1, 1e-4, 0.5, 50)
		x = VecAxpy(alpha, dir, x)
	}
	return Result{X: x, F: f(x), Iterations: iter, Converged: converged}
}

// BFGS minimizes f by the BFGS quasi-Newton method starting from x0, maintaining
// a dense approximation to the inverse Hessian that is refined by the rank-two
// BFGS update after each Armijo-backtracked line search. If grad is nil the
// gradient is estimated by central finite differences. Iteration stops when the
// gradient norm falls below tol or after maxIter steps.
func BFGS(f MultiFunc, grad GradFunc, x0 []float64, tol float64, maxIter int) Result {
	g := optimizegrad(f, grad)
	if tol <= 0 {
		tol = DefaultTol
	}
	if maxIter <= 0 {
		maxIter = DefaultMaxIter
	}
	n := len(x0)
	x := VecCopy(x0)
	gr := g(x)
	h := MatIdentity(n)
	iter := 0
	converged := false
	for ; iter < maxIter; iter++ {
		if VecNorm(gr) < tol {
			converged = true
			break
		}
		dir := VecNegate(MatVec(h, gr))
		if VecDot(dir, gr) >= 0 {
			// Reset to steepest descent if curvature estimate degrades.
			h = MatIdentity(n)
			dir = VecNegate(gr)
		}
		alpha := BacktrackingLineSearch(f, x, dir, gr, 1, 1e-4, 0.5, 50)
		s := VecScale(dir, alpha)
		xNew := VecAdd(x, s)
		grNew := g(xNew)
		y := VecSub(grNew, gr)
		sy := VecDot(s, y)
		if sy > 1e-12 {
			optimizebfgsUpdate(h, s, y, sy)
		}
		x, gr = xNew, grNew
	}
	return Result{X: x, F: f(x), Iterations: iter, Converged: converged}
}

// ---------------------------------------------------------------------------
// Derivative-free minimizers
// ---------------------------------------------------------------------------

// NelderMead minimizes f by the Nelder-Mead downhill simplex method starting
// from x0. An initial simplex is built by perturbing each coordinate of x0 by
// step. The method uses the standard reflection, expansion, contraction and
// shrink operations and requires no derivatives. Iteration stops when the spread
// of objective values across the simplex falls below tol or after maxIter
// iterations.
func NelderMead(f MultiFunc, x0 []float64, step, tol float64, maxIter int) Result {
	const (
		alpha = 1.0 // reflection
		gamma = 2.0 // expansion
		rho   = 0.5 // contraction
		sigma = 0.5 // shrink
	)
	if step == 0 {
		step = 0.05
	}
	if tol <= 0 {
		tol = DefaultTol
	}
	if maxIter <= 0 {
		maxIter = DefaultMaxIter
	}
	n := len(x0)
	simplex := make([][]float64, n+1)
	fvals := make([]float64, n+1)
	simplex[0] = VecCopy(x0)
	for i := 0; i < n; i++ {
		p := VecCopy(x0)
		if p[i] != 0 {
			p[i] += step * p[i]
		} else {
			p[i] = step
		}
		simplex[i+1] = p
	}
	for i := range simplex {
		fvals[i] = f(simplex[i])
	}
	iter := 0
	converged := false
	for ; iter < maxIter; iter++ {
		sort.Sort(optimizesimplexSorter{simplex, fvals})
		if math.Abs(fvals[n]-fvals[0]) <= tol*(math.Abs(fvals[0])+tol) {
			converged = true
			break
		}
		cen := Centroid(simplex[:n])
		worst := simplex[n]
		// Reflection.
		xr := VecLinComb(1+alpha, cen, -alpha, worst)
		fr := f(xr)
		switch {
		case fr < fvals[0]:
			// Expansion.
			xe := VecLinComb(1+gamma, cen, -gamma, worst)
			if fe := f(xe); fe < fr {
				simplex[n], fvals[n] = xe, fe
			} else {
				simplex[n], fvals[n] = xr, fr
			}
		case fr < fvals[n-1]:
			simplex[n], fvals[n] = xr, fr
		default:
			if fr < fvals[n] {
				// Outside contraction.
				xc := VecLinComb(1+rho*alpha, cen, -rho*alpha, worst)
				if fc := f(xc); fc <= fr {
					simplex[n], fvals[n] = xc, fc
					continue
				}
			} else {
				// Inside contraction.
				xc := VecLinComb(1-rho, cen, rho, worst)
				if fc := f(xc); fc < fvals[n] {
					simplex[n], fvals[n] = xc, fc
					continue
				}
			}
			// Shrink toward the best vertex.
			best := simplex[0]
			for i := 1; i <= n; i++ {
				simplex[i] = VecLinComb(1-sigma, best, sigma, simplex[i])
				fvals[i] = f(simplex[i])
			}
		}
	}
	sort.Sort(optimizesimplexSorter{simplex, fvals})
	return Result{X: VecCopy(simplex[0]), F: fvals[0], Iterations: iter, Converged: converged}
}

// SimulatedAnnealing minimizes f by simulated annealing starting from x0.
// Neighbours are proposed by adding independent Gaussian perturbations scaled by
// opts.StepSize, and uphill moves are accepted with the Metropolis probability
// exp(-Δf/T). The temperature starts at opts.InitialTemp and is multiplied by
// opts.CoolingRate each step, floored at opts.MinTemp. All randomness is drawn
// from a generator seeded by seed, so the result is fully deterministic for a
// given seed. The best point ever visited is returned.
func SimulatedAnnealing(f MultiFunc, x0 []float64, opts SAOptions, seed int64) Result {
	if opts.InitialTemp <= 0 {
		opts.InitialTemp = 10
	}
	if opts.CoolingRate <= 0 || opts.CoolingRate >= 1 {
		opts.CoolingRate = 0.995
	}
	if opts.StepSize <= 0 {
		opts.StepSize = 1
	}
	if opts.MaxIter <= 0 {
		opts.MaxIter = 10000
	}
	rng := rand.New(rand.NewSource(seed))
	x := VecCopy(x0)
	fx := f(x)
	best := VecCopy(x)
	fbest := fx
	temp := opts.InitialTemp
	iter := 0
	for ; iter < opts.MaxIter; iter++ {
		cand := make([]float64, len(x))
		for i := range x {
			cand[i] = x[i] + opts.StepSize*rng.NormFloat64()
		}
		fc := f(cand)
		if fc < fx || rng.Float64() < math.Exp((fx-fc)/temp) {
			x, fx = cand, fc
			if fc < fbest {
				best, fbest = VecCopy(cand), fc
			}
		}
		temp *= opts.CoolingRate
		if temp < opts.MinTemp {
			temp = opts.MinTemp
		}
	}
	return Result{X: best, F: fbest, Iterations: iter, Converged: true}
}

// ---------------------------------------------------------------------------
// Standard test objectives
// ---------------------------------------------------------------------------

// Sphere is the separable convex test objective sum_i x_i^2, with its unique
// minimum value 0 at the origin.
func Sphere(x []float64) float64 { return VecDot(x, x) }

// SphereGrad returns the exact gradient 2x of Sphere at x.
func SphereGrad(x []float64) []float64 { return VecScale(x, 2) }

// Rosenbrock is the classic non-convex test objective
// sum_i [100*(x_{i+1} - x_i^2)^2 + (1 - x_i)^2], with its global minimum value
// 0 at the all-ones point. It requires at least two coordinates.
func Rosenbrock(x []float64) float64 {
	var s float64
	for i := 0; i+1 < len(x); i++ {
		a := x[i+1] - x[i]*x[i]
		b := 1 - x[i]
		s += 100*a*a + b*b
	}
	return s
}

// RosenbrockGrad returns the exact gradient of Rosenbrock at x.
func RosenbrockGrad(x []float64) []float64 {
	n := len(x)
	g := make([]float64, n)
	for i := 0; i+1 < n; i++ {
		a := x[i+1] - x[i]*x[i]
		g[i] += -400*x[i]*a - 2*(1-x[i])
		g[i+1] += 200 * a
	}
	return g
}

// Booth is the two-dimensional convex test objective
// (x0 + 2*x1 - 7)^2 + (2*x0 + x1 - 5)^2, with its unique minimum value 0 at
// (1, 3). It panics if x does not have length 2.
func Booth(x []float64) float64 {
	optimizemustMatch(len(x), 2)
	a := x[0] + 2*x[1] - 7
	b := 2*x[0] + x[1] - 5
	return a*a + b*b
}

// Himmelblau is the two-dimensional test objective
// (x0^2 + x1 - 11)^2 + (x0 + x1^2 - 7)^2, which has four equal local minima of
// value 0, one of which is (3, 2). It panics if x does not have length 2.
func Himmelblau(x []float64) float64 {
	optimizemustMatch(len(x), 2)
	a := x[0]*x[0] + x[1] - 11
	b := x[0] + x[1]*x[1] - 7
	return a*a + b*b
}

// ---------------------------------------------------------------------------
// Unexported helpers (all prefixed "optimize" to avoid collision with the
// sibling file in this package).
// ---------------------------------------------------------------------------

// optimizeclamp confines v to the closed range [lo, hi].
func optimizeclamp(v, lo, hi float64) float64 {
	if v < lo {
		return lo
	}
	if v > hi {
		return hi
	}
	return v
}

// optimizemustMatch panics with ErrDimensionMismatch if a != b.
func optimizemustMatch(a, b int) {
	if a != b {
		panic(ErrDimensionMismatch)
	}
}

// optimizegrad returns grad if non-nil, otherwise a central finite-difference
// gradient estimator for f.
func optimizegrad(f MultiFunc, grad GradFunc) GradFunc {
	if grad != nil {
		return grad
	}
	return func(x []float64) []float64 { return NumericGradientCentral(f, x, DefaultStep) }
}

// optimizebfgsUpdate applies the BFGS inverse-Hessian update to h in place
// using the step s, gradient change y and their inner product sy = s . y.
func optimizebfgsUpdate(h [][]float64, s, y []float64, sy float64) {
	hy := MatVec(h, y)
	yhy := VecDot(y, hy)
	coef := (sy + yhy) / (sy * sy)
	for i := range h {
		for j := range h[i] {
			h[i][j] += coef*s[i]*s[j] - (hy[i]*s[j]+s[i]*hy[j])/sy
		}
	}
}

// optimizesimplexSorter jointly sorts a set of simplex vertices and their
// objective values into ascending order of objective value.
type optimizesimplexSorter struct {
	pts [][]float64
	fs  []float64
}

func (s optimizesimplexSorter) Len() int           { return len(s.fs) }
func (s optimizesimplexSorter) Less(i, j int) bool { return s.fs[i] < s.fs[j] }
func (s optimizesimplexSorter) Swap(i, j int) {
	s.fs[i], s.fs[j] = s.fs[j], s.fs[i]
	s.pts[i], s.pts[j] = s.pts[j], s.pts[i]
}
