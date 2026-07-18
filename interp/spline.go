// Package interp implements interpolation and approximation of one- and
// multi-dimensional sampled data.
//
// The one-dimensional constructors take a strictly increasing set of abscissae
// x together with matching ordinates y and return an interpolant object whose
// Eval method reconstructs a continuous function passing through (or, for the
// approximating variants, near) the samples. The package covers:
//
//   - Piecewise-linear interpolation (LinearInterp), together with the raw
//     scalar helpers Lerp, InverseLerp and LinearAt.
//   - Nearest-neighbour (NearestInterp) and piecewise-constant step
//     interpolation (StepInterp).
//   - Cubic splines with natural, clamped and not-a-knot boundary conditions
//     (CubicSpline), each supporting evaluation, the first three derivatives
//     and definite integration in closed form.
//   - The Akima spline (AkimaSpline), which suppresses the spurious
//     oscillations a global cubic spline can exhibit near sharp features.
//   - Monotone piecewise cubic Hermite interpolation (PCHIP), which preserves
//     the shape and monotonicity of the underlying data.
//   - Bilinear and trilinear interpolation on rectilinear grids (BilinearGrid,
//     TrilinearGrid) with the scalar helpers BilinearAt and TrilinearAt.
//
// Every routine is implemented with the Go standard library only, is
// deterministic, and validates its results against closed-form references in
// the accompanying tests.
package interp

import (
	"errors"
	"math"
)

// Sentinel errors returned by the constructors in this package.
var (
	// ErrTooFewPoints is returned when fewer sample points are supplied than
	// the chosen interpolant requires.
	ErrTooFewPoints = errors.New("interp: too few sample points")
	// ErrLengthMismatch is returned when the coordinate slices passed to a
	// constructor do not all have matching lengths.
	ErrLengthMismatch = errors.New("interp: coordinate length mismatch")
	// ErrNotSorted is returned when an abscissa slice is not strictly
	// increasing.
	ErrNotSorted = errors.New("interp: sample abscissae must be strictly increasing")
	// ErrGridShape is returned when the value array of a grid interpolant does
	// not match the lengths of its axis slices.
	ErrGridShape = errors.New("interp: grid value shape does not match axes")
)

// --- unexported helpers (all prefixed "interp" to avoid sibling collisions) ---

// interpCopy returns an independent copy of s.
func interpCopy(s []float64) []float64 {
	out := make([]float64, len(s))
	copy(out, s)
	return out
}

// interpCheckXY validates that x and y have equal length of at least min and
// that x is strictly increasing.
func interpCheckXY(x, y []float64, min int) error {
	if len(x) != len(y) {
		return ErrLengthMismatch
	}
	if len(x) < min {
		return ErrTooFewPoints
	}
	if err := interpCheckSorted(x); err != nil {
		return err
	}
	return nil
}

// interpCheckSorted reports whether x is strictly increasing.
func interpCheckSorted(x []float64) error {
	for i := 1; i < len(x); i++ {
		if !(x[i] > x[i-1]) {
			return ErrNotSorted
		}
	}
	return nil
}

// interpSearch returns the index i of the interval [xs[i], xs[i+1]] that
// contains x, clamped to the range [0, len(xs)-2]. Values outside the sampled
// range map to the first or last interval, so callers extrapolate.
func interpSearch(xs []float64, x float64) int {
	n := len(xs)
	if x <= xs[0] {
		return 0
	}
	if x >= xs[n-1] {
		return n - 2
	}
	lo, hi := 0, n-1
	for hi-lo > 1 {
		mid := int(uint(lo+hi) >> 1)
		if xs[mid] <= x {
			lo = mid
		} else {
			hi = mid
		}
	}
	return lo
}

// interpFloorIndex returns the index of the largest node not exceeding x,
// clamped to [0, len(xs)-1].
func interpFloorIndex(xs []float64, x float64) int {
	n := len(xs)
	if x < xs[0] {
		return 0
	}
	if x >= xs[n-1] {
		return n - 1
	}
	lo, hi := 0, n-1
	for hi-lo > 1 {
		mid := int(uint(lo+hi) >> 1)
		if xs[mid] <= x {
			lo = mid
		} else {
			hi = mid
		}
	}
	return lo
}

// interpCeilIndex returns the index of the smallest node not less than x,
// clamped to [0, len(xs)-1].
func interpCeilIndex(xs []float64, x float64) int {
	n := len(xs)
	if x <= xs[0] {
		return 0
	}
	if x > xs[n-1] {
		return n - 1
	}
	lo, hi := 0, n-1
	for hi-lo > 1 {
		mid := int(uint(lo+hi) >> 1)
		if xs[mid] < x {
			lo = mid
		} else {
			hi = mid
		}
	}
	return hi
}

// interpSign returns -1, 0 or +1 according to the sign of v.
func interpSign(v float64) float64 {
	switch {
	case v > 0:
		return 1
	case v < 0:
		return -1
	default:
		return 0
	}
}

// interpThomas solves a tridiagonal system with sub-diagonal sub, diagonal
// diag, super-diagonal sup and right-hand side rhs using the Thomas algorithm.
// sub[0] and sup[n-1] are ignored. The inputs are not modified.
func interpThomas(sub, diag, sup, rhs []float64) []float64 {
	n := len(diag)
	cp := make([]float64, n)
	dp := make([]float64, n)
	cp[0] = sup[0] / diag[0]
	dp[0] = rhs[0] / diag[0]
	for i := 1; i < n; i++ {
		m := diag[i] - sub[i]*cp[i-1]
		if i < n-1 {
			cp[i] = sup[i] / m
		}
		dp[i] = (rhs[i] - sub[i]*dp[i-1]) / m
	}
	x := make([]float64, n)
	x[n-1] = dp[n-1]
	for i := n - 2; i >= 0; i-- {
		x[i] = dp[i] - cp[i]*x[i+1]
	}
	return x
}

// interpSolveDense solves the dense linear system A x = b by Gaussian
// elimination with partial pivoting. A and b are not modified.
func interpSolveDense(A [][]float64, b []float64) []float64 {
	n := len(b)
	M := make([][]float64, n)
	for i := range M {
		M[i] = interpCopy(A[i])
	}
	x := interpCopy(b)
	for col := 0; col < n; col++ {
		p := col
		for r := col + 1; r < n; r++ {
			if math.Abs(M[r][col]) > math.Abs(M[p][col]) {
				p = r
			}
		}
		M[col], M[p] = M[p], M[col]
		x[col], x[p] = x[p], x[col]
		pv := M[col][col]
		for r := col + 1; r < n; r++ {
			f := M[r][col] / pv
			if f == 0 {
				continue
			}
			for k := col; k < n; k++ {
				M[r][k] -= f * M[col][k]
			}
			x[r] -= f * x[col]
		}
	}
	for i := n - 1; i >= 0; i-- {
		s := x[i]
		for k := i + 1; k < n; k++ {
			s -= M[i][k] * x[k]
		}
		x[i] = s / M[i][i]
	}
	return x
}

// interpCubic is the shared piecewise-cubic representation used by the cubic
// spline, Akima and PCHIP interpolants as well as by the piecewise-linear
// interpolant. On interval i, s(x) = y[i] + b[i]t + c[i]t^2 + d[i]t^3 where
// t = x - x[i].
type interpCubic struct {
	x, y    []float64 // nodes, length n
	b, c, d []float64 // coefficients, length n-1
}

// eval returns the interpolated value at x.
func (p *interpCubic) eval(x float64) float64 {
	i := interpSearch(p.x, x)
	t := x - p.x[i]
	return ((p.d[i]*t+p.c[i])*t+p.b[i])*t + p.y[i]
}

// deriv returns the first derivative at x.
func (p *interpCubic) deriv(x float64) float64 {
	i := interpSearch(p.x, x)
	t := x - p.x[i]
	return (3*p.d[i]*t+2*p.c[i])*t + p.b[i]
}

// deriv2 returns the second derivative at x.
func (p *interpCubic) deriv2(x float64) float64 {
	i := interpSearch(p.x, x)
	t := x - p.x[i]
	return 6*p.d[i]*t + 2*p.c[i]
}

// deriv3 returns the third derivative at x.
func (p *interpCubic) deriv3(x float64) float64 {
	i := interpSearch(p.x, x)
	return 6 * p.d[i]
}

// integral returns the definite integral of the interpolant from a to b. The
// limits are clamped to the sampled domain.
func (p *interpCubic) integral(a, b float64) float64 {
	if a == b {
		return 0
	}
	sign := 1.0
	if a > b {
		a, b = b, a
		sign = -1
	}
	lo, hi := p.x[0], p.x[len(p.x)-1]
	if a < lo {
		a = lo
	}
	if b > hi {
		b = hi
	}
	if a >= b {
		return 0
	}
	sum := 0.0
	i := interpSearch(p.x, a)
	for i < len(p.b) {
		segLo, segHi := p.x[i], p.x[i+1]
		l := math.Max(a, segLo)
		r := math.Min(b, segHi)
		if r > l {
			sum += interpAntideriv(p.y[i], p.b[i], p.c[i], p.d[i], l-segLo, r-segLo)
		}
		if segHi >= b {
			break
		}
		i++
	}
	return sign * sum
}

// interpAntideriv returns the integral of a + b t + c t^2 + d t^3 from t0 to t1.
func interpAntideriv(a, b, c, d, t0, t1 float64) float64 {
	f := func(t float64) float64 {
		return ((((d/4)*t+c/3)*t+b/2)*t + a) * t
	}
	return f(t1) - f(t0)
}

// interpHermite builds piecewise-cubic Hermite coefficients from node values y
// and node slopes m over abscissae x.
func interpHermite(x, y, m []float64) (b, c, d []float64) {
	n := len(x)
	b = make([]float64, n-1)
	c = make([]float64, n-1)
	d = make([]float64, n-1)
	for i := 0; i < n-1; i++ {
		h := x[i+1] - x[i]
		delta := (y[i+1] - y[i]) / h
		b[i] = m[i]
		c[i] = (3*delta - 2*m[i] - m[i+1]) / h
		d[i] = (m[i] + m[i+1] - 2*delta) / (h * h)
	}
	return b, c, d
}

// ---------------------------------------------------------------------------
// Scalar helpers
// ---------------------------------------------------------------------------

// Lerp returns the linear interpolation between a and b at parameter t, that is
// a + t*(b-a). t need not lie in [0, 1]; values outside extrapolate.
func Lerp(a, b, t float64) float64 {
	return a + t*(b-a)
}

// InverseLerp returns the parameter t for which Lerp(a, b, t) == v, that is
// (v-a)/(b-a). It returns 0 when a == b.
func InverseLerp(a, b, v float64) float64 {
	if a == b {
		return 0
	}
	return (v - a) / (b - a)
}

// LinearAt returns the value at x of the straight line through the points
// (x0, y0) and (x1, y1).
func LinearAt(x0, y0, x1, y1, x float64) float64 {
	return y0 + (y1-y0)*(x-x0)/(x1-x0)
}

// Clamp constrains v to the closed interval [lo, hi].
func Clamp(v, lo, hi float64) float64 {
	if v < lo {
		return lo
	}
	if v > hi {
		return hi
	}
	return v
}

// SmoothStep returns the cubic Hermite smoothstep 3u^2 - 2u^3 where u is x
// mapped from [edge0, edge1] to [0, 1] and clamped. It is 0 for x <= edge0 and
// 1 for x >= edge1, with zero slope at both ends.
func SmoothStep(edge0, edge1, x float64) float64 {
	u := Clamp((x-edge0)/(edge1-edge0), 0, 1)
	return u * u * (3 - 2*u)
}

// SmootherStep returns Perlin's quintic smoothstep 6u^5 - 15u^4 + 10u^3 where u
// is x mapped from [edge0, edge1] to [0, 1] and clamped. It has zero first and
// second derivatives at both ends.
func SmootherStep(edge0, edge1, x float64) float64 {
	u := Clamp((x-edge0)/(edge1-edge0), 0, 1)
	return u * u * u * (u*(u*6-15) + 10)
}

// BilinearAt returns the bilinear interpolation at (x, y) of the four corner
// values of the rectangle [x0, x1] x [y0, y1]: q00 at (x0, y0), q10 at
// (x1, y0), q01 at (x0, y1) and q11 at (x1, y1).
func BilinearAt(x0, x1, y0, y1, q00, q10, q01, q11, x, y float64) float64 {
	tx := (x - x0) / (x1 - x0)
	ty := (y - y0) / (y1 - y0)
	return interpBilinearUnit(q00, q10, q01, q11, tx, ty)
}

// interpBilinearUnit interpolates four corner values on the unit square at
// parameters (tx, ty).
func interpBilinearUnit(q00, q10, q01, q11, tx, ty float64) float64 {
	a := q00 + tx*(q10-q00)
	b := q01 + tx*(q11-q01)
	return a + ty*(b-a)
}

// TrilinearAt returns the trilinear interpolation at (x, y, z) of the eight
// corner values of the box [x0, x1] x [y0, y1] x [z0, z1]. Corner cijk sits at
// (x0 or x1, y0 or y1, z0 or z1) as selected by i, j, k in that order.
func TrilinearAt(x0, x1, y0, y1, z0, z1, c000, c100, c010, c110, c001, c101, c011, c111, x, y, z float64) float64 {
	tx := (x - x0) / (x1 - x0)
	ty := (y - y0) / (y1 - y0)
	tz := (z - z0) / (z1 - z0)
	return interpTrilinearUnit(c000, c100, c010, c110, c001, c101, c011, c111, tx, ty, tz)
}

// interpTrilinearUnit interpolates eight corner values on the unit cube at
// parameters (tx, ty, tz).
func interpTrilinearUnit(c000, c100, c010, c110, c001, c101, c011, c111, tx, ty, tz float64) float64 {
	z0 := interpBilinearUnit(c000, c100, c010, c110, tx, ty)
	z1 := interpBilinearUnit(c001, c101, c011, c111, tx, ty)
	return z0 + tz*(z1-z0)
}

// ---------------------------------------------------------------------------
// LinearInterp
// ---------------------------------------------------------------------------

// LinearInterp is a piecewise-linear interpolant through a set of samples.
// Evaluation outside the sampled range extrapolates along the end segments.
type LinearInterp struct {
	p interpCubic
}

// NewLinearInterp builds a piecewise-linear interpolant through the samples
// (x[i], y[i]). x must be strictly increasing and hold at least two points.
func NewLinearInterp(x, y []float64) (*LinearInterp, error) {
	if err := interpCheckXY(x, y, 2); err != nil {
		return nil, err
	}
	xc, yc := interpCopy(x), interpCopy(y)
	n := len(xc)
	b := make([]float64, n-1)
	for i := 0; i < n-1; i++ {
		b[i] = (yc[i+1] - yc[i]) / (xc[i+1] - xc[i])
	}
	return &LinearInterp{p: interpCubic{x: xc, y: yc, b: b, c: make([]float64, n-1), d: make([]float64, n-1)}}, nil
}

// Eval returns the interpolated value at x.
func (l *LinearInterp) Eval(x float64) float64 { return l.p.eval(x) }

// EvalDerivative returns the slope of the segment containing x.
func (l *LinearInterp) EvalDerivative(x float64) float64 { return l.p.deriv(x) }

// Integral returns the definite integral from a to b, with limits clamped to
// the sampled domain.
func (l *LinearInterp) Integral(a, b float64) float64 { return l.p.integral(a, b) }

// Domain returns the first and last abscissae of the samples.
func (l *LinearInterp) Domain() (lo, hi float64) {
	return l.p.x[0], l.p.x[len(l.p.x)-1]
}

// Len returns the number of samples.
func (l *LinearInterp) Len() int { return len(l.p.x) }

// Xs returns a copy of the sample abscissae.
func (l *LinearInterp) Xs() []float64 { return interpCopy(l.p.x) }

// Ys returns a copy of the sample ordinates.
func (l *LinearInterp) Ys() []float64 { return interpCopy(l.p.y) }

// EvalSlice returns the interpolant evaluated at every point of xs.
func (l *LinearInterp) EvalSlice(xs []float64) []float64 {
	out := make([]float64, len(xs))
	for i, x := range xs {
		out[i] = l.p.eval(x)
	}
	return out
}

// ---------------------------------------------------------------------------
// NearestInterp
// ---------------------------------------------------------------------------

// NearestInterp is a nearest-neighbour interpolant: it returns the ordinate of
// whichever sample abscissa is closest to the query point, breaking ties toward
// the lower index.
type NearestInterp struct {
	x, y []float64
}

// NewNearestInterp builds a nearest-neighbour interpolant. x must be strictly
// increasing and hold at least two points.
func NewNearestInterp(x, y []float64) (*NearestInterp, error) {
	if err := interpCheckXY(x, y, 2); err != nil {
		return nil, err
	}
	return &NearestInterp{x: interpCopy(x), y: interpCopy(y)}, nil
}

// Eval returns the ordinate of the sample nearest to x.
func (n *NearestInterp) Eval(x float64) float64 {
	m := len(n.x)
	if x <= n.x[0] {
		return n.y[0]
	}
	if x >= n.x[m-1] {
		return n.y[m-1]
	}
	i := interpSearch(n.x, x)
	if x-n.x[i] <= n.x[i+1]-x {
		return n.y[i]
	}
	return n.y[i+1]
}

// Domain returns the first and last abscissae of the samples.
func (n *NearestInterp) Domain() (lo, hi float64) { return n.x[0], n.x[len(n.x)-1] }

// Len returns the number of samples.
func (n *NearestInterp) Len() int { return len(n.x) }

// Xs returns a copy of the sample abscissae.
func (n *NearestInterp) Xs() []float64 { return interpCopy(n.x) }

// Ys returns a copy of the sample ordinates.
func (n *NearestInterp) Ys() []float64 { return interpCopy(n.y) }

// ---------------------------------------------------------------------------
// StepInterp
// ---------------------------------------------------------------------------

// StepInterp is a piecewise-constant (zero-order-hold) interpolant. A previous
// step holds the value of the nearest sample at or below the query point; a
// next step holds the value of the nearest sample at or above it.
type StepInterp struct {
	x, y []float64
	next bool
}

// NewPreviousStepInterp builds a piecewise-constant interpolant that holds the
// value of the nearest sample at or below the query point (a causal
// zero-order hold). x must be strictly increasing and hold at least two points.
func NewPreviousStepInterp(x, y []float64) (*StepInterp, error) {
	if err := interpCheckXY(x, y, 2); err != nil {
		return nil, err
	}
	return &StepInterp{x: interpCopy(x), y: interpCopy(y), next: false}, nil
}

// NewNextStepInterp builds a piecewise-constant interpolant that holds the
// value of the nearest sample at or above the query point. x must be strictly
// increasing and hold at least two points.
func NewNextStepInterp(x, y []float64) (*StepInterp, error) {
	if err := interpCheckXY(x, y, 2); err != nil {
		return nil, err
	}
	return &StepInterp{x: interpCopy(x), y: interpCopy(y), next: true}, nil
}

// Eval returns the piecewise-constant value at x.
func (s *StepInterp) Eval(x float64) float64 {
	if s.next {
		return s.y[interpCeilIndex(s.x, x)]
	}
	return s.y[interpFloorIndex(s.x, x)]
}

// Domain returns the first and last abscissae of the samples.
func (s *StepInterp) Domain() (lo, hi float64) { return s.x[0], s.x[len(s.x)-1] }

// Len returns the number of samples.
func (s *StepInterp) Len() int { return len(s.x) }

// Xs returns a copy of the sample abscissae.
func (s *StepInterp) Xs() []float64 { return interpCopy(s.x) }

// Ys returns a copy of the sample ordinates.
func (s *StepInterp) Ys() []float64 { return interpCopy(s.y) }

// ---------------------------------------------------------------------------
// CubicSpline
// ---------------------------------------------------------------------------

// CubicSpline is a C2-continuous piecewise-cubic interpolant. It is produced by
// one of the constructors NewNaturalCubicSpline, NewClampedCubicSpline or
// NewNotAKnotCubicSpline, which differ only in the boundary conditions imposed
// at the two ends of the data.
type CubicSpline struct {
	p interpCubic
}

// interpCubicFromMoments builds the piecewise coefficients from the node
// second-derivative moments M.
func interpCubicFromMoments(x, y, M []float64) interpCubic {
	n := len(x)
	b := make([]float64, n-1)
	c := make([]float64, n-1)
	d := make([]float64, n-1)
	for i := 0; i < n-1; i++ {
		h := x[i+1] - x[i]
		c[i] = M[i] / 2
		d[i] = (M[i+1] - M[i]) / (6 * h)
		b[i] = (y[i+1]-y[i])/h - h*(2*M[i]+M[i+1])/6
	}
	return interpCubic{x: x, y: y, b: b, c: c, d: d}
}

// NewNaturalCubicSpline builds a cubic spline with natural boundary conditions,
// that is with zero second derivative at both ends. x must be strictly
// increasing and hold at least two points.
func NewNaturalCubicSpline(x, y []float64) (*CubicSpline, error) {
	if err := interpCheckXY(x, y, 2); err != nil {
		return nil, err
	}
	xc, yc := interpCopy(x), interpCopy(y)
	n := len(xc)
	h := make([]float64, n-1)
	for i := 0; i < n-1; i++ {
		h[i] = xc[i+1] - xc[i]
	}
	sub := make([]float64, n)
	diag := make([]float64, n)
	sup := make([]float64, n)
	rhs := make([]float64, n)
	diag[0], diag[n-1] = 1, 1
	for i := 1; i < n-1; i++ {
		sub[i] = h[i-1]
		diag[i] = 2 * (h[i-1] + h[i])
		sup[i] = h[i]
		rhs[i] = 6 * ((yc[i+1]-yc[i])/h[i] - (yc[i]-yc[i-1])/h[i-1])
	}
	M := interpThomas(sub, diag, sup, rhs)
	return &CubicSpline{p: interpCubicFromMoments(xc, yc, M)}, nil
}

// NewClampedCubicSpline builds a cubic spline whose first derivative equals
// dLeft at the first node and dRight at the last node. x must be strictly
// increasing and hold at least two points.
func NewClampedCubicSpline(x, y []float64, dLeft, dRight float64) (*CubicSpline, error) {
	if err := interpCheckXY(x, y, 2); err != nil {
		return nil, err
	}
	xc, yc := interpCopy(x), interpCopy(y)
	n := len(xc)
	h := make([]float64, n-1)
	for i := 0; i < n-1; i++ {
		h[i] = xc[i+1] - xc[i]
	}
	sub := make([]float64, n)
	diag := make([]float64, n)
	sup := make([]float64, n)
	rhs := make([]float64, n)
	diag[0] = 2 * h[0]
	sup[0] = h[0]
	rhs[0] = 6 * ((yc[1]-yc[0])/h[0] - dLeft)
	sub[n-1] = h[n-2]
	diag[n-1] = 2 * h[n-2]
	rhs[n-1] = 6 * (dRight - (yc[n-1]-yc[n-2])/h[n-2])
	for i := 1; i < n-1; i++ {
		sub[i] = h[i-1]
		diag[i] = 2 * (h[i-1] + h[i])
		sup[i] = h[i]
		rhs[i] = 6 * ((yc[i+1]-yc[i])/h[i] - (yc[i]-yc[i-1])/h[i-1])
	}
	M := interpThomas(sub, diag, sup, rhs)
	return &CubicSpline{p: interpCubicFromMoments(xc, yc, M)}, nil
}

// NewNotAKnotCubicSpline builds a cubic spline with not-a-knot boundary
// conditions, requiring the third derivative to be continuous across the first
// and last interior knots. The interpolant reproduces any cubic polynomial
// exactly. x must be strictly increasing and hold at least four points.
func NewNotAKnotCubicSpline(x, y []float64) (*CubicSpline, error) {
	if err := interpCheckXY(x, y, 4); err != nil {
		return nil, err
	}
	xc, yc := interpCopy(x), interpCopy(y)
	n := len(xc)
	h := make([]float64, n-1)
	for i := 0; i < n-1; i++ {
		h[i] = xc[i+1] - xc[i]
	}
	A := make([][]float64, n)
	for i := range A {
		A[i] = make([]float64, n)
	}
	rhs := make([]float64, n)
	for i := 1; i < n-1; i++ {
		A[i][i-1] = h[i-1]
		A[i][i] = 2 * (h[i-1] + h[i])
		A[i][i+1] = h[i]
		rhs[i] = 6 * ((yc[i+1]-yc[i])/h[i] - (yc[i]-yc[i-1])/h[i-1])
	}
	A[0][0] = h[1]
	A[0][1] = -(h[0] + h[1])
	A[0][2] = h[0]
	A[n-1][n-3] = h[n-2]
	A[n-1][n-2] = -(h[n-3] + h[n-2])
	A[n-1][n-1] = h[n-3]
	M := interpSolveDense(A, rhs)
	return &CubicSpline{p: interpCubicFromMoments(xc, yc, M)}, nil
}

// Eval returns the spline value at x. Queries outside the sampled range are
// extrapolated with the boundary cubic.
func (s *CubicSpline) Eval(x float64) float64 { return s.p.eval(x) }

// EvalClamped returns the spline value at x with x first clamped to the sampled
// domain, so no extrapolation occurs.
func (s *CubicSpline) EvalClamped(x float64) float64 {
	lo, hi := s.Domain()
	return s.p.eval(Clamp(x, lo, hi))
}

// EvalDerivative returns the first derivative of the spline at x.
func (s *CubicSpline) EvalDerivative(x float64) float64 { return s.p.deriv(x) }

// EvalSecondDerivative returns the second derivative of the spline at x.
func (s *CubicSpline) EvalSecondDerivative(x float64) float64 { return s.p.deriv2(x) }

// EvalThirdDerivative returns the third derivative of the spline at x, constant
// within each interval.
func (s *CubicSpline) EvalThirdDerivative(x float64) float64 { return s.p.deriv3(x) }

// Integral returns the definite integral of the spline from a to b, with the
// limits clamped to the sampled domain.
func (s *CubicSpline) Integral(a, b float64) float64 { return s.p.integral(a, b) }

// Domain returns the first and last abscissae of the samples.
func (s *CubicSpline) Domain() (lo, hi float64) { return s.p.x[0], s.p.x[len(s.p.x)-1] }

// Len returns the number of samples.
func (s *CubicSpline) Len() int { return len(s.p.x) }

// Xs returns a copy of the sample abscissae.
func (s *CubicSpline) Xs() []float64 { return interpCopy(s.p.x) }

// Ys returns a copy of the sample ordinates.
func (s *CubicSpline) Ys() []float64 { return interpCopy(s.p.y) }

// EvalSlice returns the spline evaluated at every point of xs.
func (s *CubicSpline) EvalSlice(xs []float64) []float64 {
	out := make([]float64, len(xs))
	for i, x := range xs {
		out[i] = s.p.eval(x)
	}
	return out
}

// ---------------------------------------------------------------------------
// AkimaSpline
// ---------------------------------------------------------------------------

// AkimaSpline is a piecewise-cubic Hermite interpolant using Akima's local
// slope estimates. Because each slope depends only on nearby samples, the
// spline resists the spurious oscillations that a global cubic spline can
// produce near abrupt changes in the data.
type AkimaSpline struct {
	p interpCubic
}

// NewAkimaSpline builds an Akima spline through the samples (x[i], y[i]). x must
// be strictly increasing and hold at least three points.
func NewAkimaSpline(x, y []float64) (*AkimaSpline, error) {
	if err := interpCheckXY(x, y, 3); err != nil {
		return nil, err
	}
	xc, yc := interpCopy(x), interpCopy(y)
	n := len(xc)

	// Segment slopes with two phantom slopes prepended and appended.
	D := make([]float64, n+3)
	for i := 0; i < n-1; i++ {
		D[i+2] = (yc[i+1] - yc[i]) / (xc[i+1] - xc[i])
	}
	D[1] = 2*D[2] - D[3]
	D[0] = 2*D[1] - D[2]
	D[n+1] = 2*D[n] - D[n-1]
	D[n+2] = 2*D[n+1] - D[n]

	m := make([]float64, n)
	for i := 0; i < n; i++ {
		w1 := math.Abs(D[i+3] - D[i+2])
		w2 := math.Abs(D[i+1] - D[i])
		if w1+w2 == 0 {
			m[i] = (D[i+1] + D[i+2]) / 2
		} else {
			m[i] = (w1*D[i+1] + w2*D[i+2]) / (w1 + w2)
		}
	}
	b, c, d := interpHermite(xc, yc, m)
	return &AkimaSpline{p: interpCubic{x: xc, y: yc, b: b, c: c, d: d}}, nil
}

// Eval returns the spline value at x.
func (s *AkimaSpline) Eval(x float64) float64 { return s.p.eval(x) }

// EvalDerivative returns the first derivative of the spline at x.
func (s *AkimaSpline) EvalDerivative(x float64) float64 { return s.p.deriv(x) }

// EvalSecondDerivative returns the second derivative of the spline at x.
func (s *AkimaSpline) EvalSecondDerivative(x float64) float64 { return s.p.deriv2(x) }

// Integral returns the definite integral of the spline from a to b, with limits
// clamped to the sampled domain.
func (s *AkimaSpline) Integral(a, b float64) float64 { return s.p.integral(a, b) }

// Domain returns the first and last abscissae of the samples.
func (s *AkimaSpline) Domain() (lo, hi float64) { return s.p.x[0], s.p.x[len(s.p.x)-1] }

// Len returns the number of samples.
func (s *AkimaSpline) Len() int { return len(s.p.x) }

// Xs returns a copy of the sample abscissae.
func (s *AkimaSpline) Xs() []float64 { return interpCopy(s.p.x) }

// Ys returns a copy of the sample ordinates.
func (s *AkimaSpline) Ys() []float64 { return interpCopy(s.p.y) }

// EvalSlice returns the spline evaluated at every point of xs.
func (s *AkimaSpline) EvalSlice(xs []float64) []float64 {
	out := make([]float64, len(xs))
	for i, x := range xs {
		out[i] = s.p.eval(x)
	}
	return out
}

// ---------------------------------------------------------------------------
// PCHIP
// ---------------------------------------------------------------------------

// PCHIP is a monotone piecewise cubic Hermite interpolant using the
// Fritsch-Carlson slope limiter. It preserves the monotonicity of the data
// between samples and never overshoots, at the cost of only C1 continuity.
type PCHIP struct {
	p interpCubic
}

// interpPchipEnd computes a one-sided, shape-preserving end slope. h0 and d0 are
// the width and slope of the boundary interval; h1 and d1 those of the adjacent
// interval.
func interpPchipEnd(h0, h1, d0, d1 float64) float64 {
	m := ((2*h0+h1)*d0 - h0*d1) / (h0 + h1)
	if interpSign(m) != interpSign(d0) {
		return 0
	}
	if interpSign(d0) != interpSign(d1) && math.Abs(m) > 3*math.Abs(d0) {
		return 3 * d0
	}
	return m
}

// NewPCHIP builds a monotone piecewise cubic Hermite interpolant through the
// samples (x[i], y[i]). x must be strictly increasing and hold at least two
// points.
func NewPCHIP(x, y []float64) (*PCHIP, error) {
	if err := interpCheckXY(x, y, 2); err != nil {
		return nil, err
	}
	xc, yc := interpCopy(x), interpCopy(y)
	n := len(xc)
	h := make([]float64, n-1)
	delta := make([]float64, n-1)
	for i := 0; i < n-1; i++ {
		h[i] = xc[i+1] - xc[i]
		delta[i] = (yc[i+1] - yc[i]) / h[i]
	}
	m := make([]float64, n)
	if n == 2 {
		m[0], m[1] = delta[0], delta[0]
	} else {
		for i := 1; i < n-1; i++ {
			if delta[i-1]*delta[i] <= 0 {
				m[i] = 0
			} else {
				w1 := 2*h[i] + h[i-1]
				w2 := h[i] + 2*h[i-1]
				m[i] = (w1 + w2) / (w1/delta[i-1] + w2/delta[i])
			}
		}
		m[0] = interpPchipEnd(h[0], h[1], delta[0], delta[1])
		m[n-1] = interpPchipEnd(h[n-2], h[n-3], delta[n-2], delta[n-3])
	}
	b, c, d := interpHermite(xc, yc, m)
	return &PCHIP{p: interpCubic{x: xc, y: yc, b: b, c: c, d: d}}, nil
}

// Eval returns the interpolated value at x.
func (s *PCHIP) Eval(x float64) float64 { return s.p.eval(x) }

// EvalDerivative returns the first derivative of the interpolant at x.
func (s *PCHIP) EvalDerivative(x float64) float64 { return s.p.deriv(x) }

// EvalSecondDerivative returns the second derivative of the interpolant at x.
// Note that a PCHIP is only C1, so this is generally discontinuous at the nodes.
func (s *PCHIP) EvalSecondDerivative(x float64) float64 { return s.p.deriv2(x) }

// Integral returns the definite integral of the interpolant from a to b, with
// limits clamped to the sampled domain.
func (s *PCHIP) Integral(a, b float64) float64 { return s.p.integral(a, b) }

// IsMonotone reports whether the underlying sample ordinates are monotone
// (entirely non-decreasing or entirely non-increasing), the condition under
// which a PCHIP is itself monotone.
func (s *PCHIP) IsMonotone() bool {
	y := s.p.y
	inc, dec := true, true
	for i := 1; i < len(y); i++ {
		if y[i] < y[i-1] {
			inc = false
		}
		if y[i] > y[i-1] {
			dec = false
		}
	}
	return inc || dec
}

// Domain returns the first and last abscissae of the samples.
func (s *PCHIP) Domain() (lo, hi float64) { return s.p.x[0], s.p.x[len(s.p.x)-1] }

// Len returns the number of samples.
func (s *PCHIP) Len() int { return len(s.p.x) }

// Xs returns a copy of the sample abscissae.
func (s *PCHIP) Xs() []float64 { return interpCopy(s.p.x) }

// Ys returns a copy of the sample ordinates.
func (s *PCHIP) Ys() []float64 { return interpCopy(s.p.y) }

// EvalSlice returns the interpolant evaluated at every point of xs.
func (s *PCHIP) EvalSlice(xs []float64) []float64 {
	out := make([]float64, len(xs))
	for i, x := range xs {
		out[i] = s.p.eval(x)
	}
	return out
}

// ---------------------------------------------------------------------------
// BilinearGrid
// ---------------------------------------------------------------------------

// BilinearGrid is a bilinear interpolant on a rectilinear grid. The value at
// grid point (x[i], y[j]) is z[i][j]. Queries outside the grid extrapolate
// along the boundary cells unless EvalClamped is used.
type BilinearGrid struct {
	x, y []float64
	z    [][]float64
}

// NewBilinearGrid builds a bilinear grid interpolant. x and y must each be
// strictly increasing with at least two entries, and z must have shape
// len(x) by len(y) with z[i][j] the value at (x[i], y[j]).
func NewBilinearGrid(x, y []float64, z [][]float64) (*BilinearGrid, error) {
	if len(x) < 2 || len(y) < 2 {
		return nil, ErrTooFewPoints
	}
	if err := interpCheckSorted(x); err != nil {
		return nil, err
	}
	if err := interpCheckSorted(y); err != nil {
		return nil, err
	}
	if len(z) != len(x) {
		return nil, ErrGridShape
	}
	zc := make([][]float64, len(z))
	for i := range z {
		if len(z[i]) != len(y) {
			return nil, ErrGridShape
		}
		zc[i] = interpCopy(z[i])
	}
	return &BilinearGrid{x: interpCopy(x), y: interpCopy(y), z: zc}, nil
}

// Eval returns the bilinearly interpolated value at (x, y).
func (g *BilinearGrid) Eval(x, y float64) float64 {
	i := interpSearch(g.x, x)
	j := interpSearch(g.y, y)
	tx := (x - g.x[i]) / (g.x[i+1] - g.x[i])
	ty := (y - g.y[j]) / (g.y[j+1] - g.y[j])
	return interpBilinearUnit(g.z[i][j], g.z[i+1][j], g.z[i][j+1], g.z[i+1][j+1], tx, ty)
}

// EvalClamped returns the interpolated value at (x, y) with the query point
// first clamped to the grid extent, so no extrapolation occurs.
func (g *BilinearGrid) EvalClamped(x, y float64) float64 {
	x0, x1, y0, y1 := g.Domain()
	return g.Eval(Clamp(x, x0, x1), Clamp(y, y0, y1))
}

// Domain returns the grid extent as (xMin, xMax, yMin, yMax).
func (g *BilinearGrid) Domain() (xMin, xMax, yMin, yMax float64) {
	return g.x[0], g.x[len(g.x)-1], g.y[0], g.y[len(g.y)-1]
}

// Dims returns the number of grid points along the x and y axes.
func (g *BilinearGrid) Dims() (nx, ny int) { return len(g.x), len(g.y) }

// ---------------------------------------------------------------------------
// TrilinearGrid
// ---------------------------------------------------------------------------

// TrilinearGrid is a trilinear interpolant on a rectilinear grid. The value at
// grid point (x[i], y[j], z[k]) is v[i][j][k]. Queries outside the grid
// extrapolate along the boundary cells unless EvalClamped is used.
type TrilinearGrid struct {
	x, y, z []float64
	v       [][][]float64
}

// NewTrilinearGrid builds a trilinear grid interpolant. x, y and z must each be
// strictly increasing with at least two entries, and v must have shape
// len(x) by len(y) by len(z) with v[i][j][k] the value at (x[i], y[j], z[k]).
func NewTrilinearGrid(x, y, z []float64, v [][][]float64) (*TrilinearGrid, error) {
	if len(x) < 2 || len(y) < 2 || len(z) < 2 {
		return nil, ErrTooFewPoints
	}
	if err := interpCheckSorted(x); err != nil {
		return nil, err
	}
	if err := interpCheckSorted(y); err != nil {
		return nil, err
	}
	if err := interpCheckSorted(z); err != nil {
		return nil, err
	}
	if len(v) != len(x) {
		return nil, ErrGridShape
	}
	vc := make([][][]float64, len(v))
	for i := range v {
		if len(v[i]) != len(y) {
			return nil, ErrGridShape
		}
		vc[i] = make([][]float64, len(v[i]))
		for j := range v[i] {
			if len(v[i][j]) != len(z) {
				return nil, ErrGridShape
			}
			vc[i][j] = interpCopy(v[i][j])
		}
	}
	return &TrilinearGrid{x: interpCopy(x), y: interpCopy(y), z: interpCopy(z), v: vc}, nil
}

// Eval returns the trilinearly interpolated value at (x, y, z).
func (g *TrilinearGrid) Eval(x, y, z float64) float64 {
	i := interpSearch(g.x, x)
	j := interpSearch(g.y, y)
	k := interpSearch(g.z, z)
	tx := (x - g.x[i]) / (g.x[i+1] - g.x[i])
	ty := (y - g.y[j]) / (g.y[j+1] - g.y[j])
	tz := (z - g.z[k]) / (g.z[k+1] - g.z[k])
	return interpTrilinearUnit(
		g.v[i][j][k], g.v[i+1][j][k], g.v[i][j+1][k], g.v[i+1][j+1][k],
		g.v[i][j][k+1], g.v[i+1][j][k+1], g.v[i][j+1][k+1], g.v[i+1][j+1][k+1],
		tx, ty, tz)
}

// EvalClamped returns the interpolated value at (x, y, z) with the query point
// first clamped to the grid extent, so no extrapolation occurs.
func (g *TrilinearGrid) EvalClamped(x, y, z float64) float64 {
	x0, x1, y0, y1, z0, z1 := g.Domain()
	return g.Eval(Clamp(x, x0, x1), Clamp(y, y0, y1), Clamp(z, z0, z1))
}

// Domain returns the grid extent as (xMin, xMax, yMin, yMax, zMin, zMax).
func (g *TrilinearGrid) Domain() (xMin, xMax, yMin, yMax, zMin, zMax float64) {
	return g.x[0], g.x[len(g.x)-1], g.y[0], g.y[len(g.y)-1], g.z[0], g.z[len(g.z)-1]
}

// Dims returns the number of grid points along the x, y and z axes.
func (g *TrilinearGrid) Dims() (nx, ny, nz int) { return len(g.x), len(g.y), len(g.z) }

// ---------------------------------------------------------------------------
// One-shot convenience functions
// ---------------------------------------------------------------------------

// LinearInterpolate returns the piecewise-linear interpolation of the samples
// (x[i], y[i]) at the single query point xq. It is a convenience wrapper around
// NewLinearInterp for one-off evaluations.
func LinearInterpolate(x, y []float64, xq float64) (float64, error) {
	l, err := NewLinearInterp(x, y)
	if err != nil {
		return 0, err
	}
	return l.Eval(xq), nil
}

// NearestInterpolate returns the nearest-neighbour interpolation of the samples
// (x[i], y[i]) at the single query point xq.
func NearestInterpolate(x, y []float64, xq float64) (float64, error) {
	n, err := NewNearestInterp(x, y)
	if err != nil {
		return 0, err
	}
	return n.Eval(xq), nil
}

// SearchInterval returns the index i such that xs[i] <= x <= xs[i+1], clamped so
// that i lies in [0, len(xs)-2]. It requires len(xs) >= 2 and reports the
// containing interval for values outside the sampled range as the nearest
// boundary interval.
func SearchInterval(xs []float64, x float64) int {
	return interpSearch(xs, x)
}
