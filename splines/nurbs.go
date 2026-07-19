package splines

import "math"

// NURBSCurve is a Non-Uniform Rational B-Spline curve: a B-spline whose control
// points carry positive weights, letting it represent conic sections such as
// circles exactly.
type NURBSCurve struct {
	ctrl    []Vec
	weights []float64
	knots   []float64
	degree  int
}

// NewNURBSCurve builds a NURBS curve from control points, matching positive
// weights, a non-decreasing knot vector and a degree. The validation mirrors
// [NewBSplineCurve] and additionally requires each weight to be strictly
// positive.
func NewNURBSCurve(ctrl []Vec, weights, knots []float64, degree int) (*NURBSCurve, error) {
	if degree < 1 {
		return nil, ErrDegree
	}
	if len(ctrl) < degree+1 {
		return nil, ErrTooFewPoints
	}
	if len(weights) != len(ctrl) {
		return nil, ErrLenMismatch
	}
	if _, err := commonDim(ctrl); err != nil {
		return nil, err
	}
	if len(knots) != len(ctrl)+degree+1 {
		return nil, ErrKnots
	}
	if !nondecreasing(knots) {
		return nil, ErrKnots
	}
	for _, w := range weights {
		if w <= 0 {
			return nil, ErrWeights
		}
	}
	return &NURBSCurve{
		ctrl:    CloneVecs(ctrl),
		weights: append([]float64(nil), weights...),
		knots:   append([]float64(nil), knots...),
		degree:  degree,
	}, nil
}

// Degree returns the degree of the curve.
func (c *NURBSCurve) Degree() int { return c.degree }

// Dim returns the spatial dimension of the control points.
func (c *NURBSCurve) Dim() int { return c.ctrl[0].Dim() }

// ControlPoints returns an independent copy of the control points.
func (c *NURBSCurve) ControlPoints() []Vec { return CloneVecs(c.ctrl) }

// Weights returns an independent copy of the control-point weights.
func (c *NURBSCurve) Weights() []float64 { return append([]float64(nil), c.weights...) }

// Knots returns an independent copy of the knot vector.
func (c *NURBSCurve) Knots() []float64 { return append([]float64(nil), c.knots...) }

// Domain returns the valid parameter interval of the curve.
func (c *NURBSCurve) Domain() (lo, hi float64) {
	n := len(c.ctrl) - 1
	return c.knots[c.degree], c.knots[n+1]
}

// homogeneous returns the control points lifted to homogeneous coordinates,
// (w*x0, ..., w*x_{d-1}, w), for use with the non-rational B-spline machinery.
func (c *NURBSCurve) homogeneous() []Vec {
	d := c.Dim()
	pw := make([]Vec, len(c.ctrl))
	for i := range c.ctrl {
		h := make(Vec, d+1)
		for j := 0; j < d; j++ {
			h[j] = c.ctrl[i][j] * c.weights[i]
		}
		h[d] = c.weights[i]
		pw[i] = h
	}
	return pw
}

// Eval returns the point on the NURBS curve at parameter u by evaluating the
// weighted basis and dividing by the weight sum.
func (c *NURBSCurve) Eval(u float64) Vec {
	p := c.degree
	n := len(c.ctrl) - 1
	span := FindSpan(n, p, u, c.knots)
	N := BasisFuns(span, u, p, c.knots)
	d := c.Dim()
	num := make(Vec, d)
	var den float64
	for j := 0; j <= p; j++ {
		idx := span - p + j
		wN := N[j] * c.weights[idx]
		for k := 0; k < d; k++ {
			num[k] += wN * c.ctrl[idx][k]
		}
		den += wN
	}
	return num.Scale(1 / den)
}

// EvalDerivatives returns the point on the curve and its derivatives up to order
// d as a slice of length d+1, applying the rational-curve quotient rule of
// Algorithm A4.2 to the homogeneous B-spline derivatives.
func (c *NURBSCurve) EvalDerivatives(u float64, d int) []Vec {
	dim := c.Dim()
	bh := &BSplineCurve{ctrl: c.homogeneous(), knots: c.knots, degree: c.degree}
	hd := bh.EvalDerivatives(u, d)
	// Split into vector part Aders and scalar weight part wders.
	aders := make([]Vec, d+1)
	wders := make([]float64, d+1)
	for k := 0; k <= d; k++ {
		aders[k] = hd[k][:dim].Clone()
		wders[k] = hd[k][dim]
	}
	ck := make([]Vec, d+1)
	for k := 0; k <= d; k++ {
		v := aders[k].Clone()
		for i := 1; i <= k; i++ {
			bin := Binomial(k, i)
			v = v.AddScaled(-bin*wders[i], ck[k-i])
		}
		ck[k] = v.Scale(1 / wders[0])
	}
	return ck
}

// EvalDerivative returns the first derivative vector of the NURBS curve at
// parameter u.
func (c *NURBSCurve) EvalDerivative(u float64) Vec {
	return c.EvalDerivatives(u, 1)[1]
}

// NURBSCircle returns a full circle of the given radius centred at center,
// represented exactly as a quadratic NURBS curve with nine control points in the
// plane spanned by the first two coordinates. center may have any dimension
// >= 2; coordinates beyond the first two are copied from center unchanged.
func NURBSCircle(center Vec, radius float64) (*NURBSCurve, error) {
	if center.Dim() < 2 {
		return nil, ErrDim
	}
	if radius <= 0 {
		return nil, ErrParam
	}
	d := center.Dim()
	// Offsets of the nine control points on the unit circle's bounding square.
	offs := [9][2]float64{
		{1, 0}, {1, 1}, {0, 1}, {-1, 1}, {-1, 0},
		{-1, -1}, {0, -1}, {1, -1}, {1, 0},
	}
	ctrl := make([]Vec, 9)
	for i, o := range offs {
		p := center.Clone()
		p[0] = center[0] + radius*o[0]
		p[1] = center[1] + radius*o[1]
		_ = d
		ctrl[i] = p
	}
	s := math.Sqrt2 / 2
	weights := []float64{1, s, 1, s, 1, s, 1, s, 1}
	knots := []float64{0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1}
	return NewNURBSCurve(ctrl, weights, knots, 2)
}
