package splines

import "math"

// Binomial returns the binomial coefficient C(n, k) computed with an integer
// multiplicative recurrence that avoids intermediate overflow for moderate n.
// It returns 0 when k is out of the range [0, n].
func Binomial(n, k int) float64 {
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
	return math.Round(res)
}

// Bernstein returns the value of the Bernstein basis polynomial
// B(i,n)(t) = C(n,i) * t^i * (1-t)^(n-i).
func Bernstein(n, i int, t float64) float64 {
	if i < 0 || i > n {
		return 0
	}
	return Binomial(n, i) * math.Pow(t, float64(i)) * math.Pow(1-t, float64(n-i))
}

// BernsteinAll returns the full vector of degree-n Bernstein basis values at t,
// computed with a stable triangular recurrence (no explicit powers). The result
// has length n+1 and sums to one.
func BernsteinAll(n int, t float64) []float64 {
	b := make([]float64, n+1)
	b[0] = 1
	u := 1 - t
	for j := 1; j <= n; j++ {
		saved := 0.0
		for k := 0; k < j; k++ {
			tmp := b[k]
			b[k] = saved + u*tmp
			saved = t * tmp
		}
		b[j] = saved
	}
	return b
}

// BezierCurve is a polynomial Bezier curve defined by an ordered list of control
// points of a common dimension. Its degree is one less than the number of
// control points.
type BezierCurve struct {
	ctrl []Vec
}

// NewBezierCurve returns a Bezier curve with the given control points. It
// requires at least one control point and a uniform dimension.
func NewBezierCurve(ctrl ...Vec) (*BezierCurve, error) {
	if len(ctrl) == 0 {
		return nil, ErrTooFewPoints
	}
	if _, err := commonDim(ctrl); err != nil {
		return nil, err
	}
	return &BezierCurve{ctrl: CloneVecs(ctrl)}, nil
}

// Degree returns the polynomial degree of the curve (number of control points
// minus one).
func (b *BezierCurve) Degree() int { return len(b.ctrl) - 1 }

// Dim returns the spatial dimension of the curve's control points.
func (b *BezierCurve) Dim() int { return b.ctrl[0].Dim() }

// ControlPoints returns an independent copy of the control points.
func (b *BezierCurve) ControlPoints() []Vec { return CloneVecs(b.ctrl) }

// DeCasteljau evaluates the Bezier curve with the given control points at
// parameter t in [0,1] using de Casteljau's algorithm, the numerically stable
// standard method. It returns the point on the curve.
func DeCasteljau(ctrl []Vec, t float64) Vec {
	n := len(ctrl)
	tmp := CloneVecs(ctrl)
	for r := 1; r < n; r++ {
		for i := 0; i < n-r; i++ {
			tmp[i] = tmp[i].Lerp(tmp[i+1], t)
		}
	}
	return tmp[0]
}

// Eval returns the point on the curve at parameter t in [0,1].
func (b *BezierCurve) Eval(t float64) Vec { return DeCasteljau(b.ctrl, t) }

// Hodograph returns the derivative curve of b as a Bezier curve of degree one
// lower, whose control points are Q(i) = n*(P(i+1)-P(i)). Evaluating it gives
// the derivative of b with respect to its parameter.
func (b *BezierCurve) Hodograph() *BezierCurve {
	n := b.Degree()
	if n == 0 {
		return &BezierCurve{ctrl: []Vec{ZeroVec(b.Dim())}}
	}
	q := make([]Vec, n)
	for i := 0; i < n; i++ {
		q[i] = b.ctrl[i+1].Sub(b.ctrl[i]).Scale(float64(n))
	}
	return &BezierCurve{ctrl: q}
}

// EvalDerivative returns the first derivative vector of the curve at parameter
// t.
func (b *BezierCurve) EvalDerivative(t float64) Vec { return b.Hodograph().Eval(t) }

// EvalDerivativeK returns the k-th derivative vector of the curve at parameter
// t (k >= 0). A zero-vector is returned once k exceeds the degree.
func (b *BezierCurve) EvalDerivativeK(k int, t float64) Vec {
	cur := b
	for j := 0; j < k; j++ {
		cur = cur.Hodograph()
	}
	return cur.Eval(t)
}

// Tangent returns the unit tangent vector of the curve at parameter t. If the
// derivative vanishes it returns the raw (zero) derivative.
func (b *BezierCurve) Tangent(t float64) Vec { return b.EvalDerivative(t).Unit() }

// Elevate returns an equivalent Bezier curve of degree one higher. Degree
// elevation adds a control point without changing the curve's shape and is used
// to make two curves compatible for blending.
func (b *BezierCurve) Elevate() *BezierCurve {
	n := b.Degree()
	out := make([]Vec, n+2)
	out[0] = b.ctrl[0].Clone()
	out[n+1] = b.ctrl[n].Clone()
	for i := 1; i <= n; i++ {
		a := float64(i) / float64(n+1)
		out[i] = b.ctrl[i-1].Scale(a).Add(b.ctrl[i].Scale(1 - a))
	}
	return &BezierCurve{ctrl: out}
}

// Split subdivides the curve at parameter t and returns the two Bezier curves
// that together reproduce it: left covers the parameter range [0,t] and right
// covers [t,1], each reparameterised to [0,1].
func (b *BezierCurve) Split(t float64) (left, right *BezierCurve) {
	n := len(b.ctrl)
	tmp := CloneVecs(b.ctrl)
	lp := make([]Vec, n)
	rp := make([]Vec, n)
	lp[0] = tmp[0].Clone()
	rp[n-1] = tmp[n-1].Clone()
	for r := 1; r < n; r++ {
		for i := 0; i < n-r; i++ {
			tmp[i] = tmp[i].Lerp(tmp[i+1], t)
		}
		lp[r] = tmp[0].Clone()
		rp[n-1-r] = tmp[n-1-r].Clone()
	}
	return &BezierCurve{ctrl: lp}, &BezierCurve{ctrl: rp}
}

// Reverse returns the curve traversed in the opposite direction (control points
// reversed), so Reverse().Eval(t) == Eval(1-t).
func (b *BezierCurve) Reverse() *BezierCurve {
	n := len(b.ctrl)
	out := make([]Vec, n)
	for i := 0; i < n; i++ {
		out[i] = b.ctrl[n-1-i].Clone()
	}
	return &BezierCurve{ctrl: out}
}

// Polyline samples the curve at n+1 evenly spaced parameter values in [0,1] and
// returns the resulting points, a convenient coarse polygonal approximation.
func (b *BezierCurve) Polyline(n int) []Vec {
	if n < 1 {
		n = 1
	}
	pts := make([]Vec, n+1)
	for i := 0; i <= n; i++ {
		pts[i] = b.Eval(float64(i) / float64(n))
	}
	return pts
}

// Length returns an arc-length estimate of the curve obtained by adaptive
// Gauss-Legendre quadrature of the speed |C'(t)| over [0,1].
func (b *BezierCurve) Length() float64 {
	speed := func(t float64) float64 { return b.EvalDerivative(t).Norm() }
	return adaptiveGaussLength(speed, 0, 1, 1e-10, 24)
}
