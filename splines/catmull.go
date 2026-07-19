package splines

import "math"

// CatmullRomPoint evaluates the uniform Catmull-Rom segment defined by four
// successive control points at local parameter t in [0,1]; the returned point
// lies on the curve between p1 and p2 and interpolates them at t=0 and t=1.
func CatmullRomPoint(p0, p1, p2, p3 Vec, t float64) Vec {
	t2 := t * t
	t3 := t2 * t
	d := p0.Dim()
	out := make(Vec, d)
	for i := 0; i < d; i++ {
		out[i] = 0.5 * ((2 * p1[i]) +
			(-p0[i]+p2[i])*t +
			(2*p0[i]-5*p1[i]+4*p2[i]-p3[i])*t2 +
			(-p0[i]+3*p1[i]-3*p2[i]+p3[i])*t3)
	}
	return out
}

// CatmullRomSegment evaluates a Catmull-Rom segment between p1 and p2 using the
// Barry-Goldman non-uniform formulation with parameterisation exponent alpha
// (0 = uniform, 0.5 = centripetal, 1 = chordal) at local parameter t in [0,1].
// The centripetal variant (alpha = 0.5) avoids the cusps and self-intersections
// that the uniform form can produce on unevenly spaced data.
func CatmullRomSegment(p0, p1, p2, p3 Vec, alpha, t float64) Vec {
	knot := func(ti float64, a, b Vec) float64 {
		return ti + math.Pow(a.Dist(b), alpha)
	}
	t0 := 0.0
	t1 := knot(t0, p0, p1)
	t2 := knot(t1, p1, p2)
	t3 := knot(t2, p2, p3)
	// Guard against coincident points collapsing a knot interval.
	if t1 == t0 {
		t1 = t0 + 1
	}
	if t2 == t1 {
		t2 = t1 + 1
	}
	if t3 == t2 {
		t3 = t2 + 1
	}
	tt := t1 + t*(t2-t1)
	a1 := lerpParam(p0, p1, t0, t1, tt)
	a2 := lerpParam(p1, p2, t1, t2, tt)
	a3 := lerpParam(p2, p3, t2, t3, tt)
	b1 := lerpParam(a1, a2, t0, t2, tt)
	b2 := lerpParam(a2, a3, t1, t3, tt)
	return lerpParam(b1, b2, t1, t2, tt)
}

// lerpParam linearly interpolates between a and b as the parameter moves from
// ta to tb, evaluated at t.
func lerpParam(a, b Vec, ta, tb, t float64) Vec {
	if tb == ta {
		return a.Clone()
	}
	return a.Lerp(b, (t-ta)/(tb-ta))
}

// CatmullRomCurve is a parametric Catmull-Rom spline through a sequence of
// points. Its parameter domain is [0, n-1] for n points, with integer values
// interpolating the corresponding control points. The alpha exponent selects
// uniform, centripetal or chordal knot spacing.
type CatmullRomCurve struct {
	pts   []Vec
	alpha float64
}

// NewCatmullRomCurve builds a Catmull-Rom curve through pts using the given
// alpha (0 uniform, 0.5 centripetal, 1 chordal). It needs at least two points;
// the ends are handled by reflecting the boundary points to form phantom
// neighbours.
func NewCatmullRomCurve(pts []Vec, alpha float64) (*CatmullRomCurve, error) {
	if len(pts) < 2 {
		return nil, ErrTooFewPoints
	}
	if _, err := commonDim(pts); err != nil {
		return nil, err
	}
	return &CatmullRomCurve{pts: CloneVecs(pts), alpha: alpha}, nil
}

// Domain returns the parameter interval [0, n-1] of the curve.
func (c *CatmullRomCurve) Domain() (lo, hi float64) {
	return 0, float64(len(c.pts) - 1)
}

// neighbours returns the four control points p0..p3 for the segment starting at
// control index k, reflecting at the boundaries to synthesise phantom points.
func (c *CatmullRomCurve) neighbours(k int) (p0, p1, p2, p3 Vec) {
	n := len(c.pts)
	p1 = c.pts[k]
	p2 = c.pts[k+1]
	if k-1 >= 0 {
		p0 = c.pts[k-1]
	} else {
		p0 = p1.Scale(2).Sub(p2) // reflect
	}
	if k+2 < n {
		p3 = c.pts[k+2]
	} else {
		p3 = p2.Scale(2).Sub(p1)
	}
	return p0, p1, p2, p3
}

// Eval returns the point on the curve at parameter u in [0, n-1].
func (c *CatmullRomCurve) Eval(u float64) Vec {
	n := len(c.pts)
	if u <= 0 {
		return c.pts[0].Clone()
	}
	if u >= float64(n-1) {
		return c.pts[n-1].Clone()
	}
	k := int(math.Floor(u))
	t := u - float64(k)
	p0, p1, p2, p3 := c.neighbours(k)
	return CatmullRomSegment(p0, p1, p2, p3, c.alpha, t)
}

// EvalDerivative returns the first derivative of the curve at parameter u,
// approximated by a symmetric finite difference (exact to rounding for the
// uniform case). It lets CatmullRomCurve satisfy the [Curve] interface.
func (c *CatmullRomCurve) EvalDerivative(u float64) Vec {
	lo, hi := c.Domain()
	h := 1e-6
	a := u - h
	b := u + h
	if a < lo {
		a = lo
	}
	if b > hi {
		b = hi
	}
	return c.Eval(b).Sub(c.Eval(a)).Scale(1 / (b - a))
}

// Point returns the point on segment k (between control points k and k+1) at
// local parameter t in [0,1].
func (c *CatmullRomCurve) Point(k int, t float64) Vec {
	p0, p1, p2, p3 := c.neighbours(k)
	return CatmullRomSegment(p0, p1, p2, p3, c.alpha, t)
}
