package splines

// BoundaryType selects the end conditions used when constructing a cubic
// spline.
type BoundaryType int

const (
	// Natural sets the second derivative to zero at both ends.
	Natural BoundaryType = iota
	// Clamped fixes the first derivative at both ends to caller-supplied
	// values.
	Clamped
	// NotAKnot forces the third derivative to be continuous across the first
	// and last interior knots, so those knots are "not knots".
	NotAKnot
	// Periodic makes the spline and its first two derivatives match at the two
	// ends; it requires y[0] == y[len-1].
	Periodic
)

// CubicSpline is a piecewise cubic interpolant of one-dimensional data. It is
// produced by one of the constructors ([NaturalCubicSpline], [ClampedCubicSpline],
// [NotAKnotCubicSpline], [PeriodicCubicSpline] or [NewCubicSpline]) and is safe
// for concurrent evaluation once built.
type CubicSpline struct {
	x    []float64    // knot abscissae, strictly increasing, length n+1
	coef [][4]float64 // per-segment coefficients [c0,c1,c2,c3] in local var s=q-x[i]
	bt   BoundaryType
}

// NewCubicSpline builds a cubic spline through the samples (x[i], y[i]) with the
// requested boundary type. For [Clamped] the arguments dyLeft and dyRight give
// the imposed first derivatives at the two ends; they are ignored for the other
// boundary types. The abscissae must be strictly increasing.
func NewCubicSpline(x, y []float64, bt BoundaryType, dyLeft, dyRight float64) (*CubicSpline, error) {
	if len(x) != len(y) {
		return nil, ErrLenMismatch
	}
	if len(x) < 2 {
		return nil, ErrTooFewPoints
	}
	if !strictlyIncreasing(x) {
		return nil, ErrNotIncreasing
	}
	if bt == NotAKnot && len(x) < 4 {
		return nil, ErrTooFewPoints
	}
	if bt == Periodic && len(x) < 3 {
		return nil, ErrTooFewPoints
	}
	M, err := secondDerivatives(x, y, bt, dyLeft, dyRight)
	if err != nil {
		return nil, err
	}
	n := len(x) - 1
	coef := make([][4]float64, n)
	for i := 0; i < n; i++ {
		h := x[i+1] - x[i]
		c0 := y[i]
		c1 := (y[i+1]-y[i])/h - h/6*(2*M[i]+M[i+1])
		c2 := M[i] / 2
		c3 := (M[i+1] - M[i]) / (6 * h)
		coef[i] = [4]float64{c0, c1, c2, c3}
	}
	return &CubicSpline{x: append([]float64(nil), x...), coef: coef, bt: bt}, nil
}

// NaturalCubicSpline builds a natural cubic spline (zero second derivative at
// both ends) through (x[i], y[i]).
func NaturalCubicSpline(x, y []float64) (*CubicSpline, error) {
	return NewCubicSpline(x, y, Natural, 0, 0)
}

// ClampedCubicSpline builds a cubic spline whose first derivative equals dyLeft
// at the first knot and dyRight at the last knot.
func ClampedCubicSpline(x, y []float64, dyLeft, dyRight float64) (*CubicSpline, error) {
	return NewCubicSpline(x, y, Clamped, dyLeft, dyRight)
}

// NotAKnotCubicSpline builds a cubic spline with not-a-knot end conditions,
// the default used by many numerical packages. It requires at least four
// points.
func NotAKnotCubicSpline(x, y []float64) (*CubicSpline, error) {
	return NewCubicSpline(x, y, NotAKnot, 0, 0)
}

// PeriodicCubicSpline builds a periodic cubic spline; the first and last
// samples must have equal ordinates (y[0] == y[len-1]).
func PeriodicCubicSpline(x, y []float64) (*CubicSpline, error) {
	return NewCubicSpline(x, y, Periodic, 0, 0)
}

// secondDerivatives returns the vector of nodal second derivatives M[i] for the
// requested boundary type.
func secondDerivatives(x, y []float64, bt BoundaryType, dyL, dyR float64) ([]float64, error) {
	N := len(x)
	h := make([]float64, N-1)
	for i := 0; i < N-1; i++ {
		h[i] = x[i+1] - x[i]
	}
	if bt == Periodic {
		return periodicSecondDerivatives(x, y, h)
	}
	if bt == NotAKnot {
		return notAKnotSecondDerivatives(x, y, h)
	}
	a := make([]float64, N)
	b := make([]float64, N)
	c := make([]float64, N)
	r := make([]float64, N)
	for i := 1; i < N-1; i++ {
		a[i] = h[i-1]
		b[i] = 2 * (h[i-1] + h[i])
		c[i] = h[i]
		r[i] = 6 * ((y[i+1]-y[i])/h[i] - (y[i]-y[i-1])/h[i-1])
	}
	L := N - 1
	switch bt {
	case Natural:
		b[0], c[0], r[0] = 1, 0, 0
		a[L], b[L], r[L] = 0, 1, 0
	case Clamped:
		b[0], c[0] = 2*h[0], h[0]
		r[0] = 6 * ((y[1]-y[0])/h[0] - dyL)
		a[L], b[L] = h[L-1], 2*h[L-1]
		r[L] = 6 * (dyR - (y[L]-y[L-1])/h[L-1])
	default:
		return nil, ErrParam
	}
	return SolveTridiagonal(a, b, c, r)
}

// notAKnotSecondDerivatives solves for the nodal second derivatives under
// not-a-knot end conditions. The two boundary equations force the third
// derivative to be continuous across the first and last interior knots; each is
// used to express the endpoint second derivative (M0 and M_L) in terms of its
// neighbours, reducing the system to a tridiagonal one in the interior
// unknowns M1..M_{L-1}.
func notAKnotSecondDerivatives(x, y, h []float64) ([]float64, error) {
	N := len(x)
	L := N - 1
	m := N - 2 // interior unknowns M1..M_{L-1}
	interiorRHS := func(i int) float64 {
		return 6 * ((y[i+1]-y[i])/h[i] - (y[i]-y[i-1])/h[i-1])
	}
	a := make([]float64, m)
	b := make([]float64, m)
	c := make([]float64, m)
	r := make([]float64, m)
	// First reduced row (interior equation i=1 with M0 substituted).
	h0, h1 := h[0], h[1]
	b[0] = h0 + h0*h0/h1 + 2*(h0+h1)
	c[0] = h1 - h0*h0/h1
	r[0] = interiorRHS(1)
	// Interior rows i=2..L-2 (local index j=i-1).
	for i := 2; i <= L-2; i++ {
		j := i - 1
		a[j] = h[i-1]
		b[j] = 2 * (h[i-1] + h[i])
		c[j] = h[i]
		r[j] = interiorRHS(i)
	}
	// Last reduced row (interior equation i=L-1 with M_L substituted).
	hl2, hl1 := h[L-2], h[L-1]
	a[m-1] = hl2 - hl1*hl1/hl2
	b[m-1] = 2*(hl2+hl1) + hl1 + hl1*hl1/hl2
	r[m-1] = interiorRHS(L - 1)

	sol, err := SolveTridiagonal(a, b, c, r)
	if err != nil {
		return nil, err
	}
	M := make([]float64, N)
	for j := 0; j < m; j++ {
		M[j+1] = sol[j]
	}
	// Recover the endpoint second derivatives from the not-a-knot relations.
	M[0] = M[1]*(1+h0/h1) - (h0/h1)*M[2]
	M[L] = M[L-1]*(1+hl1/hl2) - (hl1/hl2)*M[L-2]
	return M, nil
}

// periodicSecondDerivatives solves the cyclic system for a periodic spline.
func periodicSecondDerivatives(x, y, h []float64) ([]float64, error) {
	N := len(x)
	if y[0] != y[N-1] {
		return nil, ErrParam
	}
	m := N - 1 // number of intervals and independent unknowns
	slope := make([]float64, m)
	for j := 0; j < m; j++ {
		slope[j] = (y[j+1] - y[j]) / h[j]
	}
	a := make([]float64, m)
	b := make([]float64, m)
	c := make([]float64, m)
	r := make([]float64, m)
	for i := 0; i < m; i++ {
		prev := (i - 1 + m) % m
		a[i] = h[prev]
		b[i] = 2 * (h[prev] + h[i])
		c[i] = h[i]
		r[i] = 6 * (slope[i] - slope[prev])
	}
	sol, err := SolveCyclicTridiagonal(a, b, c, r, a[0], c[m-1])
	if err != nil {
		return nil, err
	}
	M := make([]float64, N)
	copy(M, sol)
	M[N-1] = M[0]
	return M, nil
}

// Eval returns the spline value at q. Queries outside the knot range are
// evaluated by extrapolating the nearest end segment.
func (cs *CubicSpline) Eval(q float64) float64 {
	i := searchInterval(cs.x, q)
	s := q - cs.x[i]
	c := cs.coef[i]
	return c[0] + s*(c[1]+s*(c[2]+s*c[3]))
}

// EvalDerivative returns the first derivative of the spline at q.
func (cs *CubicSpline) EvalDerivative(q float64) float64 {
	i := searchInterval(cs.x, q)
	s := q - cs.x[i]
	c := cs.coef[i]
	return c[1] + s*(2*c[2]+3*s*c[3])
}

// EvalSecondDerivative returns the second derivative of the spline at q.
func (cs *CubicSpline) EvalSecondDerivative(q float64) float64 {
	i := searchInterval(cs.x, q)
	s := q - cs.x[i]
	c := cs.coef[i]
	return 2*c[2] + 6*s*c[3]
}

// EvalThirdDerivative returns the (piecewise constant) third derivative of the
// spline at q.
func (cs *CubicSpline) EvalThirdDerivative(q float64) float64 {
	i := searchInterval(cs.x, q)
	return 6 * cs.coef[i][3]
}

// Integrate returns the definite integral of the spline from a to b. The order
// a > b is allowed and yields the negated result.
func (cs *CubicSpline) Integrate(a, b float64) float64 {
	if a == b {
		return 0
	}
	sign := 1.0
	if a > b {
		a, b = b, a
		sign = -1
	}
	return sign * (cs.antideriv(b) - cs.antideriv(a))
}

// antideriv returns the value at q of the antiderivative that is zero at the
// left end of q's containing segment start accumulation, measured from x[0].
func (cs *CubicSpline) antideriv(q float64) float64 {
	total := 0.0
	seg := searchInterval(cs.x, q)
	for i := 0; i < seg; i++ {
		h := cs.x[i+1] - cs.x[i]
		total += segIntegral(cs.coef[i], h)
	}
	s := q - cs.x[seg]
	total += segIntegral(cs.coef[seg], s)
	return total
}

// segIntegral returns the integral of a local segment cubic from 0 to s.
func segIntegral(c [4]float64, s float64) float64 {
	return s * (c[0] + s*(c[1]/2+s*(c[2]/3+s*c[3]/4)))
}

// Domain returns the first and last knot abscissae.
func (cs *CubicSpline) Domain() (lo, hi float64) {
	return cs.x[0], cs.x[len(cs.x)-1]
}

// Knots returns a copy of the spline's knot abscissae.
func (cs *CubicSpline) Knots() []float64 {
	return append([]float64(nil), cs.x...)
}

// Len returns the number of polynomial segments in the spline.
func (cs *CubicSpline) Len() int { return len(cs.coef) }

// Boundary returns the boundary type used to build the spline.
func (cs *CubicSpline) Boundary() BoundaryType { return cs.bt }

// SegmentCoeffs returns the local monomial coefficients [c0,c1,c2,c3] of
// segment i, so that on that segment the spline equals
// c0 + c1*s + c2*s^2 + c3*s^3 with s = q - Knots()[i].
func (cs *CubicSpline) SegmentCoeffs(i int) [4]float64 { return cs.coef[i] }
