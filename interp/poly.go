package interp

import (
	"errors"
	"math"
)

// Sentinel errors returned by the interpolation constructors in this file.
// They use file-local unexported names so that this file never collides with
// the identifiers declared by the sibling files of package interp.
var (
	errPolyTooFew      = errors.New("interp: too few sample points")
	errPolyLenMismatch = errors.New("interp: coordinate length mismatch")
	errPolyDistinct    = errors.New("interp: sample abscissae must be distinct")
	errPolyDegree      = errors.New("interp: invalid polynomial degree")
	errPolyDim         = errors.New("interp: control points must share a common dimension")
	errPolyKnots       = errors.New("interp: invalid knot vector")
	errPolyPole        = errors.New("interp: rational interpolant has a pole at the query point")
)

// --- file-local helpers (all prefixed "interppoly" to avoid collisions) ---

// interppolyCopy returns an independent copy of s.
func interppolyCopy(s []float64) []float64 {
	out := make([]float64, len(s))
	copy(out, s)
	return out
}

// interppolyCopyMatrix returns an independent deep copy of m.
func interppolyCopyMatrix(m [][]float64) [][]float64 {
	out := make([][]float64, len(m))
	for i := range m {
		out[i] = interppolyCopy(m[i])
	}
	return out
}

// interppolyDistinct reports whether every element of x is unique.
func interppolyDistinct(x []float64) bool {
	for i := 0; i < len(x); i++ {
		for j := i + 1; j < len(x); j++ {
			if x[i] == x[j] {
				return false
			}
		}
	}
	return true
}

// interppolyBinom returns the binomial coefficient C(n, k) as a float64.
func interppolyBinom(n, k int) float64 {
	if k < 0 || k > n {
		return 0
	}
	if k > n-k {
		k = n - k
	}
	res := 1.0
	for i := 0; i < k; i++ {
		res = res * float64(n-i) / float64(i+1)
	}
	return res
}

// interppolyHornerNewton evaluates the Newton form with nodes x and
// coefficients c at xq, returning both the value and the first derivative.
func interppolyHornerNewton(x, c []float64, xq float64) (val, deriv float64) {
	n := len(c)
	if n == 0 {
		return 0, 0
	}
	val = c[n-1]
	deriv = 0
	for k := n - 2; k >= 0; k-- {
		deriv = deriv*(xq-x[k]) + val
		val = val*(xq-x[k]) + c[k]
	}
	return val, deriv
}

// interppolySolve solves the dense linear system A z = b in place using
// Gaussian elimination with partial pivoting. A is square and is overwritten.
// It returns errPolyDistinct-flavoured degeneracy as a nil result with false.
func interppolySolve(A [][]float64, b []float64) ([]float64, bool) {
	n := len(b)
	for i := 0; i < n; i++ {
		// Partial pivot.
		p := i
		max := math.Abs(A[i][i])
		for r := i + 1; r < n; r++ {
			if v := math.Abs(A[r][i]); v > max {
				max, p = v, r
			}
		}
		if max == 0 {
			return nil, false
		}
		if p != i {
			A[i], A[p] = A[p], A[i]
			b[i], b[p] = b[p], b[i]
		}
		for r := i + 1; r < n; r++ {
			f := A[r][i] / A[i][i]
			if f == 0 {
				continue
			}
			for c := i; c < n; c++ {
				A[r][c] -= f * A[i][c]
			}
			b[r] -= f * b[i]
		}
	}
	z := make([]float64, n)
	for i := n - 1; i >= 0; i-- {
		s := b[i]
		for c := i + 1; c < n; c++ {
			s -= A[i][c] * z[c]
		}
		z[i] = s / A[i][i]
	}
	return z, true
}

// interppolyClosestIndex returns the index of the sample abscissa nearest xq.
func interppolyClosestIndex(x []float64, xq float64) int {
	ns := 0
	best := math.Abs(xq - x[0])
	for i := 1; i < len(x); i++ {
		if d := math.Abs(xq - x[i]); d < best {
			best, ns = d, i
		}
	}
	return ns
}

// interppolyFindSpan returns the knot-span index i such that
// U[i] <= t < U[i+1], for a B-spline of the given degree with n+1 control
// points (so n = numCtrl-1) and knot vector U. The value is clamped to the
// half-open valid parameter range.
func interppolyFindSpan(n, degree int, t float64, U []float64) int {
	if t >= U[n+1] {
		return n
	}
	if t <= U[degree] {
		return degree
	}
	low, high := degree, n+1
	mid := (low + high) / 2
	for t < U[mid] || t >= U[mid+1] {
		if t < U[mid] {
			high = mid
		} else {
			low = mid
		}
		mid = (low + high) / 2
	}
	return mid
}

// ============================ Lagrange ============================

// LagrangeBasis returns the value at xq of the i-th Lagrange basis polynomial
// L_i built from the abscissae x. L_i equals one at x[i] and zero at every
// other node. The abscissae must be distinct.
func LagrangeBasis(x []float64, i int, xq float64) float64 {
	prod := 1.0
	for j := range x {
		if j == i {
			continue
		}
		prod *= (xq - x[j]) / (x[i] - x[j])
	}
	return prod
}

// LagrangeEval evaluates, at xq, the unique polynomial of degree < len(x) that
// passes through the samples (x[i], y[i]) using the Lagrange form. The two
// slices must have equal length and x must contain distinct abscissae.
func LagrangeEval(x, y []float64, xq float64) float64 {
	var sum float64
	for i := range x {
		sum += y[i] * LagrangeBasis(x, i, xq)
	}
	return sum
}

// LagrangeInterpolator is a reusable interpolant in Lagrange form built from a
// fixed set of samples.
type LagrangeInterpolator struct {
	x []float64
	y []float64
}

// NewLagrangeInterpolator builds a LagrangeInterpolator from the samples
// (x[i], y[i]). It requires at least one point, matching slice lengths and
// distinct abscissae.
func NewLagrangeInterpolator(x, y []float64) (*LagrangeInterpolator, error) {
	if len(x) != len(y) {
		return nil, errPolyLenMismatch
	}
	if len(x) < 1 {
		return nil, errPolyTooFew
	}
	if !interppolyDistinct(x) {
		return nil, errPolyDistinct
	}
	return &LagrangeInterpolator{x: interppolyCopy(x), y: interppolyCopy(y)}, nil
}

// Eval returns the interpolating polynomial evaluated at xq.
func (l *LagrangeInterpolator) Eval(xq float64) float64 {
	return LagrangeEval(l.x, l.y, xq)
}

// EvalSlice evaluates the interpolant at every element of xs.
func (l *LagrangeInterpolator) EvalSlice(xs []float64) []float64 {
	out := make([]float64, len(xs))
	for i, xq := range xs {
		out[i] = l.Eval(xq)
	}
	return out
}

// Len reports the number of samples used to build the interpolant.
func (l *LagrangeInterpolator) Len() int { return len(l.x) }

// Degree reports the degree of the interpolating polynomial, which is one less
// than the number of samples.
func (l *LagrangeInterpolator) Degree() int { return len(l.x) - 1 }

// ===================== Newton divided differences =====================

// DividedDifferences returns the Newton coefficients c of the interpolating
// polynomial written as c[0] + c[1](x-x0) + c[2](x-x0)(x-x1) + ... The inputs
// must have equal length and distinct abscissae.
func DividedDifferences(x, y []float64) []float64 {
	n := len(y)
	c := interppolyCopy(y)
	for j := 1; j < n; j++ {
		for i := n - 1; i >= j; i-- {
			c[i] = (c[i] - c[i-1]) / (x[i] - x[i-j])
		}
	}
	return c
}

// DividedDifferenceTable returns the full lower-triangular table of divided
// differences. Entry table[i][j] holds f[x_i, ..., x_{i+j}] and is defined for
// i+j < len(x); unused upper entries are left at zero. The leading column
// table[i][0] equals y[i] and the top row's j-th entry table[0][j] is the
// Newton coefficient c[j].
func DividedDifferenceTable(x, y []float64) [][]float64 {
	n := len(y)
	table := make([][]float64, n)
	for i := range table {
		table[i] = make([]float64, n)
		table[i][0] = y[i]
	}
	for j := 1; j < n; j++ {
		for i := 0; i < n-j; i++ {
			table[i][j] = (table[i+1][j-1] - table[i][j-1]) / (x[i+j] - x[i])
		}
	}
	return table
}

// NewtonEval evaluates the Newton form with nodes x and coefficients coeffs at
// xq using Horner's scheme. len(coeffs) must not exceed len(x).
func NewtonEval(x, coeffs []float64, xq float64) float64 {
	val, _ := interppolyHornerNewton(x, coeffs, xq)
	return val
}

// NewtonInterpolator is a reusable interpolant in Newton form that also
// supports incremental extension via AddPoint.
type NewtonInterpolator struct {
	x     []float64
	coeff []float64
}

// NewNewtonInterpolator builds a NewtonInterpolator from the samples
// (x[i], y[i]). It requires at least one point, matching lengths and distinct
// abscissae.
func NewNewtonInterpolator(x, y []float64) (*NewtonInterpolator, error) {
	if len(x) != len(y) {
		return nil, errPolyLenMismatch
	}
	if len(x) < 1 {
		return nil, errPolyTooFew
	}
	if !interppolyDistinct(x) {
		return nil, errPolyDistinct
	}
	return &NewtonInterpolator{
		x:     interppolyCopy(x),
		coeff: DividedDifferences(x, y),
	}, nil
}

// Eval returns the interpolating polynomial evaluated at xq.
func (n *NewtonInterpolator) Eval(xq float64) float64 {
	return NewtonEval(n.x, n.coeff, xq)
}

// EvalSlice evaluates the interpolant at every element of xs.
func (n *NewtonInterpolator) EvalSlice(xs []float64) []float64 {
	out := make([]float64, len(xs))
	for i, xq := range xs {
		out[i] = n.Eval(xq)
	}
	return out
}

// Coefficients returns a copy of the Newton coefficients.
func (n *NewtonInterpolator) Coefficients() []float64 { return interppolyCopy(n.coeff) }

// AddPoint extends the interpolant with one additional sample (xi, yi) in
// O(len) time, updating the Newton coefficients in place. The new abscissa must
// differ from every existing one.
func (n *NewtonInterpolator) AddPoint(xi, yi float64) error {
	for _, xv := range n.x {
		if xv == xi {
			return errPolyDistinct
		}
	}
	// The new top coefficient is the divided difference obtained by folding yi
	// through the existing nodes.
	acc := yi
	prodDen := 1.0
	// Evaluate the current Newton polynomial at xi and build up the product of
	// (xi - x[k]); the extra coefficient makes the polynomial pass through the
	// new point.
	pv := NewtonEval(n.x, n.coeff, xi)
	for _, xv := range n.x {
		prodDen *= (xi - xv)
	}
	acc = (yi - pv) / prodDen
	n.x = append(n.x, xi)
	n.coeff = append(n.coeff, acc)
	return nil
}

// Len reports the number of samples currently held by the interpolant.
func (n *NewtonInterpolator) Len() int { return len(n.x) }

// ============================ Barycentric ============================

// BarycentricWeights returns the barycentric weights w[j] = 1 / prod_{k!=j}
// (x[j] - x[k]) for the distinct abscissae x. These weights make the
// barycentric interpolation formula an O(n) evaluation.
func BarycentricWeights(x []float64) []float64 {
	n := len(x)
	w := make([]float64, n)
	for j := 0; j < n; j++ {
		prod := 1.0
		for k := 0; k < n; k++ {
			if k == j {
				continue
			}
			prod *= x[j] - x[k]
		}
		w[j] = 1 / prod
	}
	return w
}

// BarycentricEval evaluates the interpolating polynomial through (x[i], y[i])
// at xq using the second (true) barycentric formula with the supplied weights.
// If xq coincides with a node, the corresponding y value is returned exactly.
func BarycentricEval(x, y, w []float64, xq float64) float64 {
	var num, den float64
	for j := range x {
		d := xq - x[j]
		if d == 0 {
			return y[j]
		}
		t := w[j] / d
		num += t * y[j]
		den += t
	}
	return num / den
}

// BarycentricInterpolator is a reusable interpolant that evaluates in linear
// time via precomputed barycentric weights.
type BarycentricInterpolator struct {
	x []float64
	y []float64
	w []float64
}

// NewBarycentricInterpolator builds a BarycentricInterpolator from the samples
// (x[i], y[i]). It requires at least one point, matching lengths and distinct
// abscissae.
func NewBarycentricInterpolator(x, y []float64) (*BarycentricInterpolator, error) {
	if len(x) != len(y) {
		return nil, errPolyLenMismatch
	}
	if len(x) < 1 {
		return nil, errPolyTooFew
	}
	if !interppolyDistinct(x) {
		return nil, errPolyDistinct
	}
	return &BarycentricInterpolator{
		x: interppolyCopy(x),
		y: interppolyCopy(y),
		w: BarycentricWeights(x),
	}, nil
}

// Eval returns the interpolating polynomial evaluated at xq.
func (b *BarycentricInterpolator) Eval(xq float64) float64 {
	return BarycentricEval(b.x, b.y, b.w, xq)
}

// EvalSlice evaluates the interpolant at every element of xs.
func (b *BarycentricInterpolator) EvalSlice(xs []float64) []float64 {
	out := make([]float64, len(xs))
	for i, xq := range xs {
		out[i] = b.Eval(xq)
	}
	return out
}

// Weights returns a copy of the barycentric weights.
func (b *BarycentricInterpolator) Weights() []float64 { return interppolyCopy(b.w) }

// ============================ Hermite ============================

// interppolyHermiteCoeffs builds the doubled nodes z and the Newton
// coefficients of the Hermite interpolant that matches values y and first
// derivatives dy at the nodes x.
func interppolyHermiteCoeffs(x, y, dy []float64) (z, coeff []float64) {
	n := len(x)
	m := 2 * n
	z = make([]float64, m)
	Q := make([][]float64, m)
	for i := range Q {
		Q[i] = make([]float64, m)
	}
	for i := 0; i < n; i++ {
		z[2*i] = x[i]
		z[2*i+1] = x[i]
		Q[2*i][0] = y[i]
		Q[2*i+1][0] = y[i]
		Q[2*i+1][1] = dy[i]
		if i != 0 {
			Q[2*i][1] = (Q[2*i][0] - Q[2*i-1][0]) / (z[2*i] - z[2*i-1])
		}
	}
	for j := 2; j < m; j++ {
		for i := j; i < m; i++ {
			Q[i][j] = (Q[i][j-1] - Q[i-1][j-1]) / (z[i] - z[i-j])
		}
	}
	coeff = make([]float64, m)
	for i := 0; i < m; i++ {
		coeff[i] = Q[i][i]
	}
	return z, coeff
}

// HermiteEval evaluates, at xq, the osculating (Hermite) polynomial that
// matches both the values y and the first derivatives dy at the abscissae x.
// The three slices must share a common length and x must be distinct.
func HermiteEval(x, y, dy []float64, xq float64) float64 {
	z, c := interppolyHermiteCoeffs(x, y, dy)
	val, _ := interppolyHornerNewton(z, c, xq)
	return val
}

// HermiteInterpolator is a reusable osculating polynomial matching values and
// first derivatives at each node.
type HermiteInterpolator struct {
	z     []float64
	coeff []float64
}

// NewHermiteInterpolator builds a HermiteInterpolator matching the values y and
// first derivatives dy at the abscissae x. All three slices must share a length
// of at least one, and x must contain distinct abscissae.
func NewHermiteInterpolator(x, y, dy []float64) (*HermiteInterpolator, error) {
	if len(x) != len(y) || len(x) != len(dy) {
		return nil, errPolyLenMismatch
	}
	if len(x) < 1 {
		return nil, errPolyTooFew
	}
	if !interppolyDistinct(x) {
		return nil, errPolyDistinct
	}
	z, c := interppolyHermiteCoeffs(x, y, dy)
	return &HermiteInterpolator{z: z, coeff: c}, nil
}

// Eval returns the osculating polynomial evaluated at xq.
func (h *HermiteInterpolator) Eval(xq float64) float64 {
	val, _ := interppolyHornerNewton(h.z, h.coeff, xq)
	return val
}

// EvalDerivative returns the first derivative of the osculating polynomial at
// xq.
func (h *HermiteInterpolator) EvalDerivative(xq float64) float64 {
	_, d := interppolyHornerNewton(h.z, h.coeff, xq)
	return d
}

// Coefficients returns a copy of the Newton coefficients over the doubled node
// sequence.
func (h *HermiteInterpolator) Coefficients() []float64 { return interppolyCopy(h.coeff) }

// CubicHermiteSegment evaluates a single cubic Hermite segment on the unit
// interval at parameter t in [0, 1], given endpoint values p0, p1 and endpoint
// tangents m0, m1. It uses the standard Hermite basis functions.
func CubicHermiteSegment(p0, m0, p1, m1, t float64) float64 {
	t2 := t * t
	t3 := t2 * t
	h00 := 2*t3 - 3*t2 + 1
	h10 := t3 - 2*t2 + t
	h01 := -2*t3 + 3*t2
	h11 := t3 - t2
	return h00*p0 + h10*m0 + h01*p1 + h11*m1
}

// ============================ Neville ============================

// NevilleEval evaluates, at xq, the interpolating polynomial through the
// samples (x[i], y[i]) using Neville's iterated linear interpolation. The
// slices must share a length and x must be distinct.
func NevilleEval(x, y []float64, xq float64) float64 {
	n := len(x)
	p := interppolyCopy(y)
	for j := 1; j < n; j++ {
		for i := 0; i < n-j; i++ {
			p[i] = ((xq-x[i+j])*p[i] + (x[i]-xq)*p[i+1]) / (x[i] - x[i+j])
		}
	}
	return p[0]
}

// NevilleWithError evaluates the interpolating polynomial at xq and also
// returns the magnitude of the last correction, a heuristic estimate of the
// interpolation error. It follows the diagonal path of the Neville tableau
// nearest xq.
func NevilleWithError(x, y []float64, xq float64) (value, errEst float64) {
	n := len(x)
	c := interppolyCopy(y)
	d := interppolyCopy(y)
	ns := interppolyClosestIndex(x, xq)
	value = y[ns]
	ns--
	var dy float64
	for m := 1; m < n; m++ {
		for i := 0; i < n-m; i++ {
			ho := x[i] - xq
			hp := x[i+m] - xq
			w := c[i+1] - d[i]
			den := (ho - hp)
			den = w / den
			d[i] = hp * den
			c[i] = ho * den
		}
		if 2*(ns+1) < n-m {
			dy = c[ns+1]
		} else {
			dy = d[ns]
			ns--
		}
		value += dy
	}
	return value, math.Abs(dy)
}

// NevilleTable returns the full Neville tableau evaluated at xq. Row j (for
// 0 <= j < len(x)) has length len(x)-j; entry table[j][i] is the value at xq of
// the polynomial through the points i..i+j. Thus table[0] equals y and the
// single entry table[len(x)-1][0] is the final interpolated value.
func NevilleTable(x, y []float64, xq float64) [][]float64 {
	n := len(x)
	table := make([][]float64, n)
	table[0] = interppolyCopy(y)
	for j := 1; j < n; j++ {
		row := make([]float64, n-j)
		prev := table[j-1]
		for i := 0; i < n-j; i++ {
			row[i] = ((xq-x[i+j])*prev[i] + (x[i]-xq)*prev[i+1]) / (x[i] - x[i+j])
		}
		table[j] = row
	}
	return table
}

// ============================ Chebyshev ============================

// ChebyshevNodes returns the n roots of the degree-n Chebyshev polynomial of
// the first kind on [-1, 1], namely cos((2k+1)pi / (2n)) for k = 0..n-1. The
// nodes are returned in increasing order.
func ChebyshevNodes(n int) []float64 {
	out := make([]float64, n)
	for k := 0; k < n; k++ {
		out[n-1-k] = math.Cos(math.Pi * (float64(2*k + 1)) / float64(2*n))
	}
	return out
}

// ChebyshevNodesInterval returns the n Chebyshev roots mapped affinely from
// [-1, 1] onto [a, b], in increasing order.
func ChebyshevNodesInterval(n int, a, b float64) []float64 {
	base := ChebyshevNodes(n)
	out := make([]float64, n)
	for i, t := range base {
		out[i] = 0.5*(a+b) + 0.5*(b-a)*t
	}
	return out
}

// ChebyshevNodesSecondKind returns the n Chebyshev-Gauss-Lobatto points
// cos(k*pi / (n-1)) for k = 0..n-1 on [-1, 1] (the extrema of T_{n-1}),
// returned in increasing order. It requires n >= 2.
func ChebyshevNodesSecondKind(n int) []float64 {
	out := make([]float64, n)
	for k := 0; k < n; k++ {
		out[n-1-k] = math.Cos(math.Pi * float64(k) / float64(n-1))
	}
	return out
}

// ChebyshevT returns the value of the Chebyshev polynomial of the first kind of
// degree n at x, computed by the stable three-term recurrence.
func ChebyshevT(n int, x float64) float64 {
	if n == 0 {
		return 1
	}
	if n == 1 {
		return x
	}
	tm2, tm1 := 1.0, x
	var t float64
	for k := 2; k <= n; k++ {
		t = 2*x*tm1 - tm2
		tm2, tm1 = tm1, t
	}
	return t
}

// ChebyshevCoefficients samples f at the n Chebyshev roots on [a, b] and
// returns the n coefficients c of the Chebyshev series approximation
// f(x) ~= sum_{j=0}^{n-1} c[j] T_j(t) - c[0]/2, where t maps [a, b] to [-1, 1].
func ChebyshevCoefficients(f func(float64) float64, n int, a, b float64) []float64 {
	fk := make([]float64, n)
	bma := 0.5 * (b - a)
	bpa := 0.5 * (b + a)
	for k := 0; k < n; k++ {
		yk := math.Cos(math.Pi * (float64(k) + 0.5) / float64(n))
		fk[k] = f(yk*bma + bpa)
	}
	c := make([]float64, n)
	fac := 2.0 / float64(n)
	for j := 0; j < n; j++ {
		var sum float64
		for k := 0; k < n; k++ {
			sum += fk[k] * math.Cos(math.Pi*float64(j)*(float64(k)+0.5)/float64(n))
		}
		c[j] = fac * sum
	}
	return c
}

// ClenshawEval evaluates the Chebyshev series with coefficients c at the point
// y in [-1, 1] using Clenshaw recurrence, following the convention
// sum_{j} c[j] T_j(y) - c[0]/2.
func ClenshawEval(c []float64, y float64) float64 {
	var d, dd float64
	for j := len(c) - 1; j >= 1; j-- {
		d, dd = 2*y*d-dd+c[j], d
	}
	return y*d - dd + 0.5*c[0]
}

// ChebyshevInterpolator approximates a function on an interval [a, b] by
// truncated Chebyshev series obtained from samples at the Chebyshev roots.
type ChebyshevInterpolator struct {
	c    []float64
	a, b float64
}

// NewChebyshevInterpolator constructs a ChebyshevInterpolator of order n for f
// on [a, b] by sampling f at the n Chebyshev roots. It requires n >= 1 and
// a != b.
func NewChebyshevInterpolator(f func(float64) float64, n int, a, b float64) (*ChebyshevInterpolator, error) {
	if n < 1 {
		return nil, errPolyTooFew
	}
	if a == b {
		return nil, errPolyDistinct
	}
	return &ChebyshevInterpolator{c: ChebyshevCoefficients(f, n, a, b), a: a, b: b}, nil
}

// Eval returns the Chebyshev approximation at xq. Values of xq are mapped from
// [a, b] to [-1, 1] before the series is summed.
func (ci *ChebyshevInterpolator) Eval(xq float64) float64 {
	y := (2*xq - ci.a - ci.b) / (ci.b - ci.a)
	return ClenshawEval(ci.c, y)
}

// EvalSlice evaluates the approximation at every element of xs.
func (ci *ChebyshevInterpolator) EvalSlice(xs []float64) []float64 {
	out := make([]float64, len(xs))
	for i, xq := range xs {
		out[i] = ci.Eval(xq)
	}
	return out
}

// Coefficients returns a copy of the Chebyshev coefficients.
func (ci *ChebyshevInterpolator) Coefficients() []float64 { return interppolyCopy(ci.c) }

// ============================ Bezier ============================

// BernsteinBasis returns the value at t of the i-th Bernstein basis polynomial
// of degree n, C(n, i) t^i (1-t)^{n-i}.
func BernsteinBasis(n, i int, t float64) float64 {
	if i < 0 || i > n {
		return 0
	}
	return interppolyBinom(n, i) * math.Pow(t, float64(i)) * math.Pow(1-t, float64(n-i))
}

// DeCasteljau evaluates the scalar Bezier curve with the given control ordinates
// at parameter t in [0, 1] using de Casteljau's algorithm.
func DeCasteljau(points []float64, t float64) float64 {
	b := interppolyCopy(points)
	n := len(b)
	for r := 1; r < n; r++ {
		for i := 0; i < n-r; i++ {
			b[i] = (1-t)*b[i] + t*b[i+1]
		}
	}
	return b[0]
}

// BezierPoint evaluates the vector-valued Bezier curve with the given control
// points at parameter t in [0, 1] using de Casteljau's algorithm. Each control
// point is a coordinate slice and all must share the same dimension; the
// returned point has that dimension.
func BezierPoint(ctrl [][]float64, t float64) []float64 {
	n := len(ctrl)
	if n == 0 {
		return nil
	}
	dim := len(ctrl[0])
	b := interppolyCopyMatrix(ctrl)
	for r := 1; r < n; r++ {
		for i := 0; i < n-r; i++ {
			for k := 0; k < dim; k++ {
				b[i][k] = (1-t)*b[i][k] + t*b[i+1][k]
			}
		}
	}
	return b[0]
}

// BezierCurve evaluates the vector-valued Bezier curve at every parameter in ts
// and returns the resulting points in order.
func BezierCurve(ctrl [][]float64, ts []float64) [][]float64 {
	out := make([][]float64, len(ts))
	for i, t := range ts {
		out[i] = BezierPoint(ctrl, t)
	}
	return out
}

// BezierDerivative returns the derivative vector dC/dt of the Bezier curve at
// parameter t. The derivative is itself a degree n-1 Bezier curve whose control
// points are n (P_{i+1} - P_i).
func BezierDerivative(ctrl [][]float64, t float64) []float64 {
	n := len(ctrl)
	if n < 2 {
		if n == 1 {
			return make([]float64, len(ctrl[0]))
		}
		return nil
	}
	dim := len(ctrl[0])
	deg := n - 1
	q := make([][]float64, deg)
	for i := 0; i < deg; i++ {
		q[i] = make([]float64, dim)
		for k := 0; k < dim; k++ {
			q[i][k] = float64(deg) * (ctrl[i+1][k] - ctrl[i][k])
		}
	}
	return BezierPoint(q, t)
}

// BezierElevate returns the control points of the degree-elevated Bezier curve.
// The returned curve has one extra control point but describes exactly the same
// curve as the input.
func BezierElevate(ctrl [][]float64) [][]float64 {
	n := len(ctrl)
	if n == 0 {
		return nil
	}
	dim := len(ctrl[0])
	deg := n - 1
	out := make([][]float64, n+1)
	out[0] = interppolyCopy(ctrl[0])
	out[n] = interppolyCopy(ctrl[deg])
	for i := 1; i <= deg; i++ {
		alpha := float64(i) / float64(deg+1)
		p := make([]float64, dim)
		for k := 0; k < dim; k++ {
			p[k] = alpha*ctrl[i-1][k] + (1-alpha)*ctrl[i][k]
		}
		out[i] = p
	}
	return out
}

// Bezier is a reusable vector-valued Bezier curve of arbitrary degree.
type Bezier struct {
	ctrl [][]float64
}

// NewBezier builds a Bezier curve from the given control points. It requires at
// least one control point and a consistent dimension across all of them.
func NewBezier(ctrl [][]float64) (*Bezier, error) {
	if len(ctrl) < 1 {
		return nil, errPolyTooFew
	}
	dim := len(ctrl[0])
	for _, p := range ctrl {
		if len(p) != dim {
			return nil, errPolyDim
		}
	}
	return &Bezier{ctrl: interppolyCopyMatrix(ctrl)}, nil
}

// Eval returns the point on the curve at parameter t in [0, 1].
func (b *Bezier) Eval(t float64) []float64 { return BezierPoint(b.ctrl, t) }

// Derivative returns the derivative vector of the curve at parameter t.
func (b *Bezier) Derivative(t float64) []float64 { return BezierDerivative(b.ctrl, t) }

// Degree reports the polynomial degree of the curve, one less than the number
// of control points.
func (b *Bezier) Degree() int { return len(b.ctrl) - 1 }

// ============================ B-spline ============================

// BSplineBasis returns the value at t of the i-th B-spline basis function of
// the given degree over the knot vector knots, computed with the Cox-de Boor
// recurrence.
func BSplineBasis(i, degree int, knots []float64, t float64) float64 {
	if degree == 0 {
		if knots[i] <= t && t < knots[i+1] {
			return 1
		}
		// Include the right endpoint of the final span so the basis is
		// right-continuous at the domain end.
		if t == knots[len(knots)-1] && t == knots[i+1] && i+1 == len(knots)-1 {
			return 1
		}
		return 0
	}
	var left, right float64
	dl := knots[i+degree] - knots[i]
	if dl != 0 {
		left = (t - knots[i]) / dl * BSplineBasis(i, degree-1, knots, t)
	}
	dr := knots[i+degree+1] - knots[i+1]
	if dr != 0 {
		right = (knots[i+degree+1] - t) / dr * BSplineBasis(i+1, degree-1, knots, t)
	}
	return left + right
}

// UniformKnots returns a uniform (periodic-style) knot vector 0, 1, ..., m for a
// B-spline of the given degree with numCtrl control points. Its length is
// numCtrl + degree + 1 and its valid parameter range is [degree, numCtrl].
func UniformKnots(numCtrl, degree int) []float64 {
	m := numCtrl + degree + 1
	knots := make([]float64, m)
	for i := 0; i < m; i++ {
		knots[i] = float64(i)
	}
	return knots
}

// OpenUniformKnots returns a clamped (open uniform) knot vector on [0, 1] for a
// B-spline of the given degree with numCtrl control points, so that the curve
// interpolates its first and last control points. It requires
// numCtrl >= degree + 1.
func OpenUniformKnots(numCtrl, degree int) []float64 {
	m := numCtrl + degree + 1
	knots := make([]float64, m)
	interior := numCtrl - degree - 1
	for i := 0; i <= degree; i++ {
		knots[i] = 0
	}
	for i := 1; i <= interior; i++ {
		knots[degree+i] = float64(i) / float64(interior+1)
	}
	for i := m - degree - 1; i < m; i++ {
		knots[i] = 1
	}
	return knots
}

// DeBoor evaluates the vector-valued B-spline curve of the given degree with
// knot vector knots and control points ctrl at parameter t using de Boor's
// algorithm. t is clamped to the valid parameter range.
func DeBoor(degree int, knots []float64, ctrl [][]float64, t float64) []float64 {
	n := len(ctrl) - 1
	dim := len(ctrl[0])
	lo, hi := knots[degree], knots[n+1]
	if t < lo {
		t = lo
	}
	if t > hi {
		t = hi
	}
	s := interppolyFindSpan(n, degree, t, knots)
	d := make([][]float64, degree+1)
	for j := 0; j <= degree; j++ {
		d[j] = interppolyCopy(ctrl[j+s-degree])
	}
	for r := 1; r <= degree; r++ {
		for j := degree; j >= r; j-- {
			i := j + s - degree
			den := knots[i+degree-r+1] - knots[i]
			var alpha float64
			if den != 0 {
				alpha = (t - knots[i]) / den
			}
			for k := 0; k < dim; k++ {
				d[j][k] = (1-alpha)*d[j-1][k] + alpha*d[j][k]
			}
		}
	}
	return d[degree]
}

// BSpline is a reusable vector-valued B-spline curve.
type BSpline struct {
	degree int
	knots  []float64
	ctrl   [][]float64
}

// NewBSpline builds a B-spline of the given degree from the knot vector and
// control points. The knot vector length must equal len(ctrl) + degree + 1, the
// degree must be at least one, there must be at least degree+1 control points,
// and all control points must share a dimension.
func NewBSpline(degree int, knots []float64, ctrl [][]float64) (*BSpline, error) {
	if degree < 1 {
		return nil, errPolyDegree
	}
	if len(ctrl) < degree+1 {
		return nil, errPolyTooFew
	}
	if len(knots) != len(ctrl)+degree+1 {
		return nil, errPolyKnots
	}
	dim := len(ctrl[0])
	for _, p := range ctrl {
		if len(p) != dim {
			return nil, errPolyDim
		}
	}
	for i := 1; i < len(knots); i++ {
		if knots[i] < knots[i-1] {
			return nil, errPolyKnots
		}
	}
	return &BSpline{degree: degree, knots: interppolyCopy(knots), ctrl: interppolyCopyMatrix(ctrl)}, nil
}

// Eval returns the point on the curve at parameter t, clamped to the domain.
func (s *BSpline) Eval(t float64) []float64 {
	return DeBoor(s.degree, s.knots, s.ctrl, t)
}

// Domain returns the valid parameter range [lo, hi] of the curve.
func (s *BSpline) Domain() (lo, hi float64) {
	return s.knots[s.degree], s.knots[len(s.ctrl)]
}

// Degree reports the polynomial degree of the B-spline.
func (s *BSpline) Degree() int { return s.degree }

// ======================= Thiele rational =======================

// ThieleReciprocalDifferences returns the coefficients a of Thiele's continued
// fraction for the samples (x[i], y[i]):
//
//	R(x) = a[0] + (x-x0)/(a[1] + (x-x1)/(a[2] + ...)).
//
// The abscissae must be distinct and the reciprocal differences must be finite.
func ThieleReciprocalDifferences(x, y []float64) []float64 {
	n := len(x)
	rho := make([][]float64, n)
	for i := range rho {
		rho[i] = make([]float64, n)
		rho[i][0] = y[i]
	}
	for i := 0; i < n-1; i++ {
		rho[i][1] = (x[i] - x[i+1]) / (rho[i][0] - rho[i+1][0])
	}
	for j := 2; j < n; j++ {
		for i := 0; i < n-j; i++ {
			rho[i][j] = rho[i+1][j-2] + (x[i]-x[i+j])/(rho[i][j-1]-rho[i+1][j-1])
		}
	}
	a := make([]float64, n)
	a[0] = rho[0][0]
	if n > 1 {
		a[1] = rho[0][1]
	}
	for k := 2; k < n; k++ {
		a[k] = rho[0][k] - rho[0][k-2]
	}
	return a
}

// ThieleEval evaluates Thiele's continued fraction with nodes x and
// coefficients a (as returned by ThieleReciprocalDifferences) at xq.
func ThieleEval(x, a []float64, xq float64) float64 {
	n := len(a)
	if n == 0 {
		return 0
	}
	val := a[n-1]
	for k := n - 2; k >= 0; k-- {
		val = a[k] + (xq-x[k])/val
	}
	return val
}

// ThieleInterpolator is a reusable rational interpolant in Thiele
// continued-fraction form.
type ThieleInterpolator struct {
	x []float64
	a []float64
}

// NewThieleInterpolator builds a ThieleInterpolator from the samples
// (x[i], y[i]). It requires at least one point, matching lengths and distinct
// abscissae.
func NewThieleInterpolator(x, y []float64) (*ThieleInterpolator, error) {
	if len(x) != len(y) {
		return nil, errPolyLenMismatch
	}
	if len(x) < 1 {
		return nil, errPolyTooFew
	}
	if !interppolyDistinct(x) {
		return nil, errPolyDistinct
	}
	return &ThieleInterpolator{x: interppolyCopy(x), a: ThieleReciprocalDifferences(x, y)}, nil
}

// Eval returns the rational interpolant evaluated at xq.
func (t *ThieleInterpolator) Eval(xq float64) float64 {
	return ThieleEval(t.x, t.a, xq)
}

// EvalSlice evaluates the interpolant at every element of xs.
func (t *ThieleInterpolator) EvalSlice(xs []float64) []float64 {
	out := make([]float64, len(xs))
	for i, xq := range xs {
		out[i] = t.Eval(xq)
	}
	return out
}

// RationalInterpEval performs Bulirsch-Stoer diagonal rational function
// interpolation of the samples (x[i], y[i]) at xq, returning the interpolated
// value. It returns errPolyPole if the query point sits at a pole of the
// interpolant.
func RationalInterpEval(x, y []float64, xq float64) (float64, error) {
	const tiny = 1e-25
	n := len(x)
	c := interppolyCopy(y)
	d := make([]float64, n)
	for i := range d {
		d[i] = y[i] + tiny
	}
	ns := interppolyClosestIndex(x, xq)
	if x[ns] == xq {
		return y[ns], nil
	}
	value := y[ns]
	ns--
	for m := 1; m < n; m++ {
		for i := 0; i < n-m; i++ {
			w := c[i+1] - d[i]
			h := x[i+m] - xq
			tt := (x[i] - xq) * d[i] / h
			dd := tt - c[i+1]
			if dd == 0 {
				return 0, errPolyPole
			}
			dd = w / dd
			d[i] = c[i+1] * dd
			c[i] = tt * dd
		}
		var dy float64
		if 2*(ns+1) < n-m {
			dy = c[ns+1]
		} else {
			dy = d[ns]
			ns--
		}
		value += dy
	}
	return value, nil
}

// ===================== Least-squares polynomial fit =====================

// PolyVal evaluates the polynomial with ascending coefficients coeffs (so
// coeffs[0] is the constant term) at x using Horner's method.
func PolyVal(coeffs []float64, x float64) float64 {
	n := len(coeffs)
	if n == 0 {
		return 0
	}
	res := coeffs[n-1]
	for i := n - 2; i >= 0; i-- {
		res = res*x + coeffs[i]
	}
	return res
}

// PolyDerivativeCoeffs returns the ascending coefficients of the derivative of
// the polynomial given by the ascending coefficients coeffs.
func PolyDerivativeCoeffs(coeffs []float64) []float64 {
	if len(coeffs) <= 1 {
		return []float64{}
	}
	out := make([]float64, len(coeffs)-1)
	for i := 1; i < len(coeffs); i++ {
		out[i-1] = float64(i) * coeffs[i]
	}
	return out
}

// PolyIntegralCoeffs returns the ascending coefficients of an antiderivative of
// the polynomial given by coeffs, using constant of integration c.
func PolyIntegralCoeffs(coeffs []float64, c float64) []float64 {
	out := make([]float64, len(coeffs)+1)
	out[0] = c
	for i := 0; i < len(coeffs); i++ {
		out[i+1] = coeffs[i] / float64(i+1)
	}
	return out
}

// interppolyNormalFit solves the weighted least-squares normal equations for a
// polynomial of the given degree and returns its ascending coefficients.
func interppolyNormalFit(x, y, w []float64, degree int) ([]float64, error) {
	if len(x) != len(y) {
		return nil, errPolyLenMismatch
	}
	if degree < 0 {
		return nil, errPolyDegree
	}
	m := degree + 1
	if len(x) < m {
		return nil, errPolyTooFew
	}
	// Precompute weighted power sums to fill the symmetric normal matrix.
	sums := make([]float64, 2*degree+1)
	rhs := make([]float64, m)
	for k := range x {
		wk := 1.0
		if w != nil {
			wk = w[k]
		}
		xp := 1.0
		for p := 0; p < len(sums); p++ {
			sums[p] += wk * xp
			xp *= x[k]
		}
		xp = 1.0
		for p := 0; p < m; p++ {
			rhs[p] += wk * y[k] * xp
			xp *= x[k]
		}
	}
	A := make([][]float64, m)
	for i := 0; i < m; i++ {
		A[i] = make([]float64, m)
		for j := 0; j < m; j++ {
			A[i][j] = sums[i+j]
		}
	}
	sol, ok := interppolySolve(A, rhs)
	if !ok {
		return nil, errPolyDistinct
	}
	return sol, nil
}

// PolyFit returns the ascending coefficients of the degree-`degree`
// least-squares polynomial that best fits the samples (x[i], y[i]) in the
// ordinary (unweighted) sense. It requires at least degree+1 samples.
func PolyFit(x, y []float64, degree int) ([]float64, error) {
	return interppolyNormalFit(x, y, nil, degree)
}

// PolyFitWeighted returns the ascending coefficients of the degree-`degree`
// weighted least-squares polynomial fit, where w[i] is the weight of sample i.
func PolyFitWeighted(x, y, w []float64, degree int) ([]float64, error) {
	if len(w) != len(x) {
		return nil, errPolyLenMismatch
	}
	return interppolyNormalFit(x, y, w, degree)
}

// VandermondeSolve returns the ascending coefficients of the unique polynomial
// of degree len(x)-1 that passes exactly through the samples (x[i], y[i]),
// obtained by solving the Vandermonde system. The abscissae must be distinct.
func VandermondeSolve(x, y []float64) ([]float64, error) {
	if len(x) != len(y) {
		return nil, errPolyLenMismatch
	}
	n := len(x)
	if n < 1 {
		return nil, errPolyTooFew
	}
	A := make([][]float64, n)
	for i := 0; i < n; i++ {
		A[i] = make([]float64, n)
		xp := 1.0
		for j := 0; j < n; j++ {
			A[i][j] = xp
			xp *= x[i]
		}
	}
	sol, ok := interppolySolve(A, interppolyCopy(y))
	if !ok {
		return nil, errPolyDistinct
	}
	return sol, nil
}

// LeastSquaresLine returns the slope and intercept of the ordinary
// least-squares line y = slope*x + intercept fitting the samples. It requires
// at least two points with non-degenerate abscissae.
func LeastSquaresLine(x, y []float64) (slope, intercept float64, err error) {
	if len(x) != len(y) {
		return 0, 0, errPolyLenMismatch
	}
	n := len(x)
	if n < 2 {
		return 0, 0, errPolyTooFew
	}
	var sx, sy, sxx, sxy float64
	for i := 0; i < n; i++ {
		sx += x[i]
		sy += y[i]
		sxx += x[i] * x[i]
		sxy += x[i] * y[i]
	}
	fn := float64(n)
	den := fn*sxx - sx*sx
	if den == 0 {
		return 0, 0, errPolyDistinct
	}
	slope = (fn*sxy - sx*sy) / den
	intercept = (sy - slope*sx) / fn
	return slope, intercept, nil
}

// RSquared returns the coefficient of determination R^2 of the polynomial with
// ascending coefficients coeffs against the samples (x[i], y[i]).
func RSquared(x, y, coeffs []float64) float64 {
	var mean float64
	for _, v := range y {
		mean += v
	}
	mean /= float64(len(y))
	var ssRes, ssTot float64
	for i := range x {
		r := y[i] - PolyVal(coeffs, x[i])
		ssRes += r * r
		d := y[i] - mean
		ssTot += d * d
	}
	if ssTot == 0 {
		return 1
	}
	return 1 - ssRes/ssTot
}

// RMSE returns the root-mean-square error of the polynomial with ascending
// coefficients coeffs against the samples (x[i], y[i]).
func RMSE(x, y, coeffs []float64) float64 {
	var ss float64
	for i := range x {
		r := y[i] - PolyVal(coeffs, x[i])
		ss += r * r
	}
	return math.Sqrt(ss / float64(len(x)))
}

// PolynomialFit bundles the ascending coefficients of a least-squares fit with
// convenient evaluation.
type PolynomialFit struct {
	coeffs []float64
	degree int
}

// NewPolynomialFit builds a PolynomialFit of the given degree from the samples
// (x[i], y[i]) by ordinary least squares.
func NewPolynomialFit(x, y []float64, degree int) (*PolynomialFit, error) {
	c, err := PolyFit(x, y, degree)
	if err != nil {
		return nil, err
	}
	return &PolynomialFit{coeffs: c, degree: degree}, nil
}

// Eval evaluates the fitted polynomial at x.
func (p *PolynomialFit) Eval(x float64) float64 { return PolyVal(p.coeffs, x) }

// Coefficients returns a copy of the ascending fit coefficients.
func (p *PolynomialFit) Coefficients() []float64 { return interppolyCopy(p.coeffs) }

// Degree reports the degree of the fitted polynomial.
func (p *PolynomialFit) Degree() int { return p.degree }

// ===================== Trigonometric interpolation =====================

// interppolyTrigCoeffs returns the real Fourier cosine and sine coefficients
// a[k] = (2/N) sum_j y_j cos(k x_j) and b[k] = (2/N) sum_j y_j sin(k x_j) for
// k = 0..N/2, where the samples y_j are taken at x_j = 2*pi*j/N.
func interppolyTrigCoeffs(y []float64) (a, b []float64) {
	N := len(y)
	M := N / 2
	a = make([]float64, M+1)
	b = make([]float64, M+1)
	fac := 2.0 / float64(N)
	for k := 0; k <= M; k++ {
		var ca, cb float64
		for j := 0; j < N; j++ {
			ang := 2 * math.Pi * float64(k) * float64(j) / float64(N)
			ca += y[j] * math.Cos(ang)
			cb += y[j] * math.Sin(ang)
		}
		a[k] = fac * ca
		b[k] = fac * cb
	}
	return a, b
}

// TrigCoefficients returns the real trigonometric interpolation coefficients
// a and b (each of length N/2+1) for samples y taken at the equally spaced
// nodes x_j = 2*pi*j/N on [0, 2*pi). See TrigInterpEval for the reconstruction
// convention.
func TrigCoefficients(y []float64) (a, b []float64) {
	return interppolyTrigCoeffs(y)
}

// TrigInterpEval evaluates, at xq, the trigonometric polynomial that
// interpolates the samples y taken at the equally spaced nodes
// x_j = 2*pi*j/len(y) on [0, 2*pi). The interpolant is real and, for an even
// number of samples, uses the standard half-weight on the Nyquist term.
func TrigInterpEval(y []float64, xq float64) float64 {
	N := len(y)
	if N == 0 {
		return 0
	}
	a, b := interppolyTrigCoeffs(y)
	M := N / 2
	sum := a[0] / 2
	last := M
	if N%2 == 0 {
		last = M - 1
	}
	for k := 1; k <= last; k++ {
		sum += a[k]*math.Cos(float64(k)*xq) + b[k]*math.Sin(float64(k)*xq)
	}
	if N%2 == 0 {
		sum += 0.5 * a[M] * math.Cos(float64(M)*xq)
	}
	return sum
}

// TrigInterpolator is a reusable trigonometric interpolant of equally spaced
// periodic samples.
type TrigInterpolator struct {
	y []float64
}

// NewTrigInterpolator builds a TrigInterpolator from samples y taken at the
// nodes x_j = 2*pi*j/len(y). It requires at least one sample.
func NewTrigInterpolator(y []float64) (*TrigInterpolator, error) {
	if len(y) < 1 {
		return nil, errPolyTooFew
	}
	return &TrigInterpolator{y: interppolyCopy(y)}, nil
}

// Eval returns the trigonometric interpolant evaluated at xq.
func (t *TrigInterpolator) Eval(xq float64) float64 { return TrigInterpEval(t.y, xq) }

// DFTReal returns the discrete Fourier transform of the real input y, with
// re[k] and im[k] holding the real and imaginary parts of
// X_k = sum_j y_j exp(-2*pi*i*j*k/N).
func DFTReal(y []float64) (re, im []float64) {
	N := len(y)
	re = make([]float64, N)
	im = make([]float64, N)
	for k := 0; k < N; k++ {
		var sr, si float64
		for j := 0; j < N; j++ {
			ang := -2 * math.Pi * float64(k) * float64(j) / float64(N)
			sr += y[j] * math.Cos(ang)
			si += y[j] * math.Sin(ang)
		}
		re[k] = sr
		im[k] = si
	}
	return re, im
}
