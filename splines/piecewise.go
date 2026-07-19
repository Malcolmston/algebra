package splines

import "math"

// PiecewiseCubic is a one-dimensional interpolant made of independent cubic
// pieces, each expressed in a local monomial basis. It backs the Hermite,
// Catmull-Rom, cardinal, monotone (PCHIP) and Akima constructors.
type PiecewiseCubic struct {
	x    []float64
	coef [][4]float64
}

// buildFromTangents assembles a PiecewiseCubic from nodal values y and nodal
// first derivatives (tangents) m using cubic Hermite interpolation on each
// segment.
func buildFromTangents(x, y, m []float64) *PiecewiseCubic {
	n := len(x) - 1
	coef := make([][4]float64, n)
	for i := 0; i < n; i++ {
		h := x[i+1] - x[i]
		dy := y[i+1] - y[i]
		c0 := y[i]
		c1 := m[i]
		c2 := (3*dy - h*(2*m[i]+m[i+1])) / (h * h)
		c3 := (-2*dy + h*(m[i]+m[i+1])) / (h * h * h)
		coef[i] = [4]float64{c0, c1, c2, c3}
	}
	return &PiecewiseCubic{x: append([]float64(nil), x...), coef: coef}
}

// HermiteSpline builds a piecewise cubic Hermite interpolant from nodal values
// y and nodal first derivatives m at the abscissae x. All three slices must
// share a length of at least two and x must be strictly increasing.
func HermiteSpline(x, y, m []float64) (*PiecewiseCubic, error) {
	if len(x) != len(y) || len(x) != len(m) {
		return nil, ErrLenMismatch
	}
	if len(x) < 2 {
		return nil, ErrTooFewPoints
	}
	if !strictlyIncreasing(x) {
		return nil, ErrNotIncreasing
	}
	return buildFromTangents(x, y, m), nil
}

// secants returns the divided differences d[k] = (y[k+1]-y[k])/(x[k+1]-x[k]).
func secants(x, y []float64) []float64 {
	d := make([]float64, len(x)-1)
	for k := range d {
		d[k] = (y[k+1] - y[k]) / (x[k+1] - x[k])
	}
	return d
}

// CardinalSpline builds a cardinal spline through (x[i], y[i]) with the given
// tension in [0,1]; tension 0 reproduces the Catmull-Rom spline while tension 1
// produces zero tangents (a piecewise-linear-like flat interpolant). The
// interior tangents use the non-uniform central difference and the ends use
// one-sided differences.
func CardinalSpline(x, y []float64, tension float64) (*PiecewiseCubic, error) {
	if len(x) != len(y) {
		return nil, ErrLenMismatch
	}
	if len(x) < 2 {
		return nil, ErrTooFewPoints
	}
	if !strictlyIncreasing(x) {
		return nil, ErrNotIncreasing
	}
	n := len(x)
	m := make([]float64, n)
	scale := 1 - tension
	m[0] = scale * (y[1] - y[0]) / (x[1] - x[0])
	m[n-1] = scale * (y[n-1] - y[n-2]) / (x[n-1] - x[n-2])
	for i := 1; i < n-1; i++ {
		m[i] = scale * (y[i+1] - y[i-1]) / (x[i+1] - x[i-1])
	}
	return buildFromTangents(x, y, m), nil
}

// CatmullRomSpline builds the Catmull-Rom interpolant through (x[i], y[i]); it
// is the cardinal spline with zero tension.
func CatmullRomSpline(x, y []float64) (*PiecewiseCubic, error) {
	return CardinalSpline(x, y, 0)
}

// MonotoneCubicSpline builds a shape-preserving monotone cubic interpolant
// using the Fritsch-Carlson method (the PCHIP algorithm). Where the sample data
// is monotone the interpolant is monotone too, avoiding the overshoot of an
// unconstrained cubic spline.
func MonotoneCubicSpline(x, y []float64) (*PiecewiseCubic, error) {
	if len(x) != len(y) {
		return nil, ErrLenMismatch
	}
	if len(x) < 2 {
		return nil, ErrTooFewPoints
	}
	if !strictlyIncreasing(x) {
		return nil, ErrNotIncreasing
	}
	n := len(x)
	if n == 2 {
		s := (y[1] - y[0]) / (x[1] - x[0])
		return buildFromTangents(x, y, []float64{s, s}), nil
	}
	h := make([]float64, n-1)
	d := make([]float64, n-1)
	for k := 0; k < n-1; k++ {
		h[k] = x[k+1] - x[k]
		d[k] = (y[k+1] - y[k]) / h[k]
	}
	m := make([]float64, n)
	// Interior tangents: weighted harmonic mean, zeroed at extrema.
	for i := 1; i < n-1; i++ {
		if d[i-1]*d[i] <= 0 {
			m[i] = 0
			continue
		}
		w1 := 2*h[i] + h[i-1]
		w2 := h[i] + 2*h[i-1]
		m[i] = (w1 + w2) / (w1/d[i-1] + w2/d[i])
	}
	// Endpoint tangents with the shape-preserving limiter.
	m[0] = pchipEnd(h[0], h[1], d[0], d[1])
	m[n-1] = pchipEnd(h[n-2], h[n-3], d[n-2], d[n-3])
	return buildFromTangents(x, y, m), nil
}

// pchipEnd computes a one-sided endpoint tangent with the PCHIP limiter, given
// the two adjacent interval widths (h0 nearest the end) and secant slopes.
func pchipEnd(h0, h1, d0, d1 float64) float64 {
	m := ((2*h0+h1)*d0 - h0*d1) / (h0 + h1)
	if sign(m) != sign(d0) {
		return 0
	}
	if sign(d0) != sign(d1) && math.Abs(m) > 3*math.Abs(d0) {
		return 3 * d0
	}
	return m
}

// AkimaSpline builds an Akima cubic interpolant through (x[i], y[i]). Akima's
// method uses locally weighted slopes, making it robust against wiggles caused
// by outliers while remaining smooth (C1).
func AkimaSpline(x, y []float64) (*PiecewiseCubic, error) {
	if len(x) != len(y) {
		return nil, ErrLenMismatch
	}
	if len(x) < 3 {
		return nil, ErrTooFewPoints
	}
	if !strictlyIncreasing(x) {
		return nil, ErrNotIncreasing
	}
	n := len(x)
	d := secants(x, y) // length n-1
	// Extended slope array e with two synthetic slopes on each side.
	e := make([]float64, n+3) // indices 0..n+2, e[k] = d[k-2]
	for k := 0; k < n-1; k++ {
		e[k+2] = d[k]
	}
	e[1] = 2*e[2] - e[3]
	e[0] = 2*e[1] - e[2]
	e[n] = 2*e[n-1] - e[n-2]
	e[n+1] = 2*e[n] - e[n-1]
	e[n+2] = 2*e[n+1] - e[n]
	m := make([]float64, n)
	for i := 0; i < n; i++ {
		w1 := math.Abs(e[i+3] - e[i+2])
		w2 := math.Abs(e[i+1] - e[i])
		if w1+w2 == 0 {
			m[i] = (e[i+1] + e[i+2]) / 2
		} else {
			m[i] = (w1*e[i+1] + w2*e[i+2]) / (w1 + w2)
		}
	}
	return buildFromTangents(x, y, m), nil
}

// Eval returns the interpolant value at q, extrapolating with the nearest end
// piece for queries outside the knot range.
func (p *PiecewiseCubic) Eval(q float64) float64 {
	i := searchInterval(p.x, q)
	s := q - p.x[i]
	c := p.coef[i]
	return c[0] + s*(c[1]+s*(c[2]+s*c[3]))
}

// EvalDerivative returns the first derivative of the interpolant at q.
func (p *PiecewiseCubic) EvalDerivative(q float64) float64 {
	i := searchInterval(p.x, q)
	s := q - p.x[i]
	c := p.coef[i]
	return c[1] + s*(2*c[2]+3*s*c[3])
}

// EvalSecondDerivative returns the second derivative of the interpolant at q.
func (p *PiecewiseCubic) EvalSecondDerivative(q float64) float64 {
	i := searchInterval(p.x, q)
	s := q - p.x[i]
	c := p.coef[i]
	return 2*c[2] + 6*s*c[3]
}

// Integrate returns the definite integral of the interpolant from a to b.
func (p *PiecewiseCubic) Integrate(a, b float64) float64 {
	if a == b {
		return 0
	}
	sign := 1.0
	if a > b {
		a, b = b, a
		sign = -1
	}
	return sign * (p.antideriv(b) - p.antideriv(a))
}

func (p *PiecewiseCubic) antideriv(q float64) float64 {
	total := 0.0
	seg := searchInterval(p.x, q)
	for i := 0; i < seg; i++ {
		h := p.x[i+1] - p.x[i]
		total += segIntegral(p.coef[i], h)
	}
	total += segIntegral(p.coef[seg], q-p.x[seg])
	return total
}

// Domain returns the first and last abscissae of the interpolant.
func (p *PiecewiseCubic) Domain() (lo, hi float64) { return p.x[0], p.x[len(p.x)-1] }

// Knots returns a copy of the interpolant's abscissae.
func (p *PiecewiseCubic) Knots() []float64 { return append([]float64(nil), p.x...) }

// Len returns the number of cubic pieces.
func (p *PiecewiseCubic) Len() int { return len(p.coef) }

// SegmentCoeffs returns the local monomial coefficients of piece i.
func (p *PiecewiseCubic) SegmentCoeffs(i int) [4]float64 { return p.coef[i] }

// sign returns -1, 0 or +1 following the sign of v.
func sign(v float64) float64 {
	switch {
	case v > 0:
		return 1
	case v < 0:
		return -1
	default:
		return 0
	}
}

// HermiteBasis00 is the Hermite basis function h00(t)=2t^3-3t^2+1, the weight
// of the start value in a unit cubic Hermite segment.
func HermiteBasis00(t float64) float64 { return t*t*(2*t-3) + 1 }

// HermiteBasis10 is the Hermite basis function h10(t)=t^3-2t^2+t, the weight of
// the (scaled) start tangent.
func HermiteBasis10(t float64) float64 { return t * (t*(t-2) + 1) }

// HermiteBasis01 is the Hermite basis function h01(t)=-2t^3+3t^2, the weight of
// the end value.
func HermiteBasis01(t float64) float64 { return t * t * (3 - 2*t) }

// HermiteBasis11 is the Hermite basis function h11(t)=t^3-t^2, the weight of the
// (scaled) end tangent.
func HermiteBasis11(t float64) float64 { return t * t * (t - 1) }

// CubicHermite evaluates a single cubic Hermite segment with endpoints y0, y1,
// endpoint tangents m0, m1, interval width h and local parameter t in [0,1].
func CubicHermite(y0, y1, m0, m1, h, t float64) float64 {
	return HermiteBasis00(t)*y0 + HermiteBasis10(t)*h*m0 +
		HermiteBasis01(t)*y1 + HermiteBasis11(t)*h*m1
}

// HermiteCurvePoint evaluates a vector-valued cubic Hermite segment with
// endpoints p0, p1 and endpoint tangents (derivatives) t0, t1 at local
// parameter t in [0,1]. The tangents are given directly (already scaled to the
// unit interval).
func HermiteCurvePoint(p0, p1, t0, t1 Vec, t float64) Vec {
	h00 := HermiteBasis00(t)
	h10 := HermiteBasis10(t)
	h01 := HermiteBasis01(t)
	h11 := HermiteBasis11(t)
	d := p0.Dim()
	out := make(Vec, d)
	for i := 0; i < d; i++ {
		out[i] = h00*p0[i] + h10*t0[i] + h01*p1[i] + h11*t1[i]
	}
	return out
}
