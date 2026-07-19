package splines

import (
	"math"
	"sort"
)

// Curve is the common interface implemented by the parametric curve types in
// this package. Domain reports the valid parameter interval, Eval returns the
// point at a parameter, and EvalDerivative returns the first derivative vector.
type Curve interface {
	Eval(t float64) Vec
	EvalDerivative(t float64) Vec
	Domain() (lo, hi float64)
}

// Domain returns the parameter domain of a Bezier curve, always [0,1].
func (b *BezierCurve) Domain() (lo, hi float64) { return 0, 1 }

// GaussLegendreRule returns the nodes and weights of the n-point Gauss-Legendre
// quadrature rule on the interval [-1,1]. The nodes are the roots of the
// degree-n Legendre polynomial found by Newton's iteration. It requires n >= 1.
func GaussLegendreRule(n int) (nodes, weights []float64) {
	nodes = make([]float64, n)
	weights = make([]float64, n)
	m := (n + 1) / 2
	for i := 0; i < m; i++ {
		// Initial guess for the i-th root.
		z := math.Cos(math.Pi * (float64(i) + 0.75) / (float64(n) + 0.5))
		var pp float64
		for iter := 0; iter < 100; iter++ {
			p1, p2 := 1.0, 0.0
			for j := 0; j < n; j++ {
				p3 := p2
				p2 = p1
				p1 = ((2*float64(j)+1)*z*p2 - float64(j)*p3) / float64(j+1)
			}
			pp = float64(n) * (z*p1 - p2) / (z*z - 1)
			z1 := z
			z = z1 - p1/pp
			if math.Abs(z-z1) < 1e-15 {
				break
			}
		}
		nodes[i] = -z
		nodes[n-1-i] = z
		w := 2 / ((1 - z*z) * pp * pp)
		weights[i] = w
		weights[n-1-i] = w
	}
	return nodes, weights
}

// GaussLegendreIntegrate approximates the integral of f over [a,b] with the
// n-point Gauss-Legendre rule.
func GaussLegendreIntegrate(f func(float64) float64, a, b float64, n int) float64 {
	nodes, weights := GaussLegendreRule(n)
	c1 := (b - a) / 2
	c2 := (b + a) / 2
	var s float64
	for i := range nodes {
		s += weights[i] * f(c1*nodes[i]+c2)
	}
	return c1 * s
}

// adaptiveGaussLength integrates the non-negative speed function f over [a,b]
// with recursive bisection driven by a 10-point Gauss-Legendre rule until the
// halves agree with the whole to within tol or the recursion depth is reached.
func adaptiveGaussLength(f func(float64) float64, a, b, tol float64, maxDepth int) float64 {
	whole := GaussLegendreIntegrate(f, a, b, 10)
	return adaptiveGaussRec(f, a, b, tol, whole, maxDepth)
}

func adaptiveGaussRec(f func(float64) float64, a, b, tol, whole float64, depth int) float64 {
	mid := (a + b) / 2
	left := GaussLegendreIntegrate(f, a, mid, 10)
	right := GaussLegendreIntegrate(f, mid, b, 10)
	if depth <= 0 || math.Abs(left+right-whole) <= tol {
		return left + right
	}
	return adaptiveGaussRec(f, a, mid, tol/2, left, depth-1) +
		adaptiveGaussRec(f, mid, b, tol/2, right, depth-1)
}

// CurveLength returns the total arc length of a curve over its full parameter
// domain, computed by adaptive Gauss-Legendre quadrature of the speed |C'(t)|.
func CurveLength(c Curve) float64 {
	lo, hi := c.Domain()
	return CurveLengthBetween(c, lo, hi)
}

// CurveLengthBetween returns the arc length of a curve between parameters a and
// b (a <= b assumed within the domain).
func CurveLengthBetween(c Curve, a, b float64) float64 {
	speed := func(t float64) float64 { return c.EvalDerivative(t).Norm() }
	return adaptiveGaussLength(speed, a, b, 1e-11, 30)
}

// ArcLengthParam provides an arc-length (unit-speed) reparameterisation of a
// curve. It stores a monotone table mapping parameter t to cumulative arc
// length s and inverts it to answer both directions.
type ArcLengthParam struct {
	c  Curve
	ts []float64 // sample parameters, increasing
	ss []float64 // cumulative arc length at each ts, increasing from 0
}

// NewArcLengthParam builds an arc-length parameterisation of c using n sample
// segments (n >= 1) across the curve's domain. More samples give a more
// accurate length table and faster inversion.
func NewArcLengthParam(c Curve, n int) (*ArcLengthParam, error) {
	if n < 1 {
		return nil, ErrTooFewPoints
	}
	lo, hi := c.Domain()
	ts := make([]float64, n+1)
	ss := make([]float64, n+1)
	ts[0] = lo
	ss[0] = 0
	for i := 1; i <= n; i++ {
		ts[i] = lo + (hi-lo)*float64(i)/float64(n)
		ss[i] = ss[i-1] + CurveLengthBetween(c, ts[i-1], ts[i])
	}
	return &ArcLengthParam{c: c, ts: ts, ss: ss}, nil
}

// TotalLength returns the total arc length recorded by the parameterisation.
func (ap *ArcLengthParam) TotalLength() float64 { return ap.ss[len(ap.ss)-1] }

// SForT returns the arc length from the start of the curve to parameter t using
// linear interpolation within the length table.
func (ap *ArcLengthParam) SForT(t float64) float64 {
	if t <= ap.ts[0] {
		return 0
	}
	if t >= ap.ts[len(ap.ts)-1] {
		return ap.TotalLength()
	}
	i := searchInterval(ap.ts, t)
	frac := (t - ap.ts[i]) / (ap.ts[i+1] - ap.ts[i])
	return ap.ss[i] + frac*(ap.ss[i+1]-ap.ss[i])
}

// TForS returns the curve parameter t at which the cumulative arc length equals
// s. The table is inverted by binary search and refined by a Newton step using
// the local speed, giving near machine-accurate unit-speed sampling.
func (ap *ArcLengthParam) TForS(s float64) float64 {
	total := ap.TotalLength()
	if s <= 0 {
		return ap.ts[0]
	}
	if s >= total {
		return ap.ts[len(ap.ts)-1]
	}
	i := sort.SearchFloat64s(ap.ss, s)
	if i > 0 {
		i--
	}
	ds := ap.ss[i+1] - ap.ss[i]
	frac := 0.0
	if ds > 0 {
		frac = (s - ap.ss[i]) / ds
	}
	t := ap.ts[i] + frac*(ap.ts[i+1]-ap.ts[i])
	// Newton refinement: solve S(t) - s = 0 with S'(t) = speed(t).
	for iter := 0; iter < 8; iter++ {
		st := ap.ss[i] + CurveLengthBetween(ap.c, ap.ts[i], t)
		sp := ap.c.EvalDerivative(t).Norm()
		if sp == 0 {
			break
		}
		dt := (st - s) / sp
		t -= dt
		if t < ap.ts[i] {
			t = ap.ts[i]
		}
		if t > ap.ts[i+1] {
			t = ap.ts[i+1]
		}
		if math.Abs(dt) < 1e-12 {
			break
		}
	}
	return t
}

// PointAtArcLength returns the point on the curve at cumulative arc length s,
// i.e. the unit-speed sample. It is shorthand for Eval(TForS(s)).
func (ap *ArcLengthParam) PointAtArcLength(s float64) Vec {
	return ap.c.Eval(ap.TForS(s))
}

// PointAtFraction returns the point on the curve at fraction f in [0,1] of the
// total arc length; f=0 is the start and f=1 is the end.
func (ap *ArcLengthParam) PointAtFraction(f float64) Vec {
	return ap.PointAtArcLength(f * ap.TotalLength())
}
