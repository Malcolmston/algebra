package splines

// FindSpan returns the knot span index i such that U[i] <= u < U[i+1], for a
// B-spline of degree p with n+1 control points (so n = len(ctrl)-1) and knot
// vector U. It follows Algorithm A2.1 of The NURBS Book and clamps u to the
// last span at the right end of the domain.
func FindSpan(n, p int, u float64, U []float64) int {
	if u >= U[n+1] {
		return n
	}
	if u <= U[p] {
		return p
	}
	low, high := p, n+1
	mid := (low + high) / 2
	for u < U[mid] || u >= U[mid+1] {
		if u < U[mid] {
			high = mid
		} else {
			low = mid
		}
		mid = (low + high) / 2
	}
	return mid
}

// BasisFuns returns the p+1 non-zero B-spline basis function values at u for the
// span i (as returned by [FindSpan]). It uses the stable triangular recurrence
// of Algorithm A2.2, avoiding division by zero.
func BasisFuns(i int, u float64, p int, U []float64) []float64 {
	N := make([]float64, p+1)
	left := make([]float64, p+1)
	right := make([]float64, p+1)
	N[0] = 1
	for j := 1; j <= p; j++ {
		left[j] = u - U[i+1-j]
		right[j] = U[i+j] - u
		saved := 0.0
		for r := 0; r < j; r++ {
			temp := N[r] / (right[r+1] + left[j-r])
			N[r] = saved + right[r+1]*temp
			saved = left[j-r] * temp
		}
		N[j] = saved
	}
	return N
}

// DersBasisFuns returns the non-zero basis functions and their derivatives up to
// order d, evaluated at u for span i. The result ders[k][j] is the k-th
// derivative (k = 0..d) of the j-th non-zero basis function (j = 0..p). It
// implements Algorithm A2.3.
func DersBasisFuns(i int, u float64, p, d int, U []float64) [][]float64 {
	ndu := make([][]float64, p+1)
	for r := range ndu {
		ndu[r] = make([]float64, p+1)
	}
	left := make([]float64, p+1)
	right := make([]float64, p+1)
	ndu[0][0] = 1
	for j := 1; j <= p; j++ {
		left[j] = u - U[i+1-j]
		right[j] = U[i+j] - u
		saved := 0.0
		for r := 0; r < j; r++ {
			ndu[j][r] = right[r+1] + left[j-r]
			temp := ndu[r][j-1] / ndu[j][r]
			ndu[r][j] = saved + right[r+1]*temp
			saved = left[j-r] * temp
		}
		ndu[j][j] = saved
	}
	ders := make([][]float64, d+1)
	for k := range ders {
		ders[k] = make([]float64, p+1)
	}
	for j := 0; j <= p; j++ {
		ders[0][j] = ndu[j][p]
	}
	a := [2][]float64{make([]float64, p+1), make([]float64, p+1)}
	for r := 0; r <= p; r++ {
		s1, s2 := 0, 1
		a[0][0] = 1
		for k := 1; k <= d; k++ {
			dd := 0.0
			rk := r - k
			pk := p - k
			if r >= k {
				a[s2][0] = a[s1][0] / ndu[pk+1][rk]
				dd = a[s2][0] * ndu[rk][pk]
			}
			var j1, j2 int
			if rk >= -1 {
				j1 = 1
			} else {
				j1 = -rk
			}
			if r-1 <= pk {
				j2 = k - 1
			} else {
				j2 = p - r
			}
			for j := j1; j <= j2; j++ {
				a[s2][j] = (a[s1][j] - a[s1][j-1]) / ndu[pk+1][rk+j]
				dd += a[s2][j] * ndu[rk+j][pk]
			}
			if r <= pk {
				a[s2][k] = -a[s1][k-1] / ndu[pk+1][r]
				dd += a[s2][k] * ndu[r][pk]
			}
			ders[k][r] = dd
			s1, s2 = s2, s1
		}
	}
	rr := p
	for k := 1; k <= d; k++ {
		for j := 0; j <= p; j++ {
			ders[k][j] *= float64(rr)
		}
		rr *= p - k
	}
	return ders
}

// OneBasisFun returns the value of the single B-spline basis function N(i,p)(u)
// for knot vector U, implementing Algorithm A2.4. It is useful when only one
// basis function is required rather than the whole non-zero set.
func OneBasisFun(p int, U []float64, i int, u float64) float64 {
	m := len(U) - 1
	if (i == 0 && u == U[0]) || (i == m-p-1 && u == U[m]) {
		return 1
	}
	if u < U[i] || u >= U[i+p+1] {
		return 0
	}
	N := make([]float64, p+1)
	for j := 0; j <= p; j++ {
		if u >= U[i+j] && u < U[i+j+1] {
			N[j] = 1
		}
	}
	for k := 1; k <= p; k++ {
		var saved float64
		if N[0] != 0 {
			saved = (u - U[i]) * N[0] / (U[i+k] - U[i])
		}
		for j := 0; j < p-k+1; j++ {
			uLeft := U[i+j+1]
			uRight := U[i+j+k+1]
			if N[j+1] == 0 {
				N[j] = saved
				saved = 0
			} else {
				temp := N[j+1] / (uRight - uLeft)
				N[j] = saved + (uRight-u)*temp
				saved = (u - uLeft) * temp
			}
		}
	}
	return N[0]
}

// DeBoor evaluates a B-spline curve of degree p at parameter u directly with de
// Boor's algorithm, given the control points and knot vector. It is an
// allocation-light alternative to building a [BSplineCurve].
func DeBoor(p int, U []float64, ctrl []Vec, u float64) Vec {
	n := len(ctrl) - 1
	span := FindSpan(n, p, u, U)
	d := make([]Vec, p+1)
	for j := 0; j <= p; j++ {
		d[j] = ctrl[span-p+j].Clone()
	}
	for r := 1; r <= p; r++ {
		for j := p; j >= r; j-- {
			i := span - p + j
			denom := U[i+p-r+1] - U[i]
			var alpha float64
			if denom != 0 {
				alpha = (u - U[i]) / denom
			}
			d[j] = d[j-1].Scale(1 - alpha).Add(d[j].Scale(alpha))
		}
	}
	return d[p]
}

// UniformKnots returns a uniform (non-clamped) knot vector of the correct length
// len(U) = numCtrl + degree + 1 with unit spacing starting at zero.
func UniformKnots(numCtrl, degree int) []float64 {
	m := numCtrl + degree + 1
	U := make([]float64, m)
	for i := range U {
		U[i] = float64(i)
	}
	return U
}

// ClampedKnots returns an open uniform (clamped) knot vector so that the curve
// interpolates its first and last control points. The interior knots are evenly
// spaced on [0,1] and the ends have multiplicity degree+1.
func ClampedKnots(numCtrl, degree int) []float64 {
	m := numCtrl + degree + 1
	U := make([]float64, m)
	interior := numCtrl - degree - 1 // number of distinct interior knots
	for i := 0; i <= degree; i++ {
		U[i] = 0
		U[m-1-i] = 1
	}
	for i := 1; i <= interior; i++ {
		U[degree+i] = float64(i) / float64(interior+1)
	}
	return U
}

// BSplineCurve is a (non-rational) B-spline curve defined by control points, a
// knot vector and a degree.
type BSplineCurve struct {
	ctrl   []Vec
	knots  []float64
	degree int
}

// NewBSplineCurve builds a B-spline curve from control points, a non-decreasing
// knot vector and degree p. It validates that len(knots) == len(ctrl)+p+1, that
// there are at least p+1 control points and that the knot vector is
// non-decreasing.
func NewBSplineCurve(ctrl []Vec, knots []float64, degree int) (*BSplineCurve, error) {
	if degree < 1 {
		return nil, ErrDegree
	}
	if len(ctrl) < degree+1 {
		return nil, ErrTooFewPoints
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
	return &BSplineCurve{
		ctrl:   CloneVecs(ctrl),
		knots:  append([]float64(nil), knots...),
		degree: degree,
	}, nil
}

// Degree returns the degree of the curve.
func (b *BSplineCurve) Degree() int { return b.degree }

// Dim returns the spatial dimension of the control points.
func (b *BSplineCurve) Dim() int { return b.ctrl[0].Dim() }

// ControlPoints returns an independent copy of the control points.
func (b *BSplineCurve) ControlPoints() []Vec { return CloneVecs(b.ctrl) }

// Knots returns an independent copy of the knot vector.
func (b *BSplineCurve) Knots() []float64 { return append([]float64(nil), b.knots...) }

// Domain returns the valid parameter interval [U[p], U[n+1]] of the curve.
func (b *BSplineCurve) Domain() (lo, hi float64) {
	n := len(b.ctrl) - 1
	return b.knots[b.degree], b.knots[n+1]
}

// Eval returns the point on the curve at parameter u using the basis functions.
func (b *BSplineCurve) Eval(u float64) Vec {
	p := b.degree
	n := len(b.ctrl) - 1
	span := FindSpan(n, p, u, b.knots)
	N := BasisFuns(span, u, p, b.knots)
	out := ZeroVec(b.Dim())
	for j := 0; j <= p; j++ {
		out = out.AddScaled(N[j], b.ctrl[span-p+j])
	}
	return out
}

// EvalDerivatives returns the point on the curve together with its derivatives
// up to order d as a slice of length d+1 (index 0 is the point). It follows
// Algorithm A3.2.
func (b *BSplineCurve) EvalDerivatives(u float64, d int) []Vec {
	p := b.degree
	n := len(b.ctrl) - 1
	du := d
	if du > p {
		du = p
	}
	out := make([]Vec, d+1)
	for k := range out {
		out[k] = ZeroVec(b.Dim())
	}
	span := FindSpan(n, p, u, b.knots)
	nders := DersBasisFuns(span, u, p, du, b.knots)
	for k := 0; k <= du; k++ {
		acc := ZeroVec(b.Dim())
		for j := 0; j <= p; j++ {
			acc = acc.AddScaled(nders[k][j], b.ctrl[span-p+j])
		}
		out[k] = acc
	}
	return out
}

// EvalDerivative returns the first derivative vector of the curve at parameter
// u.
func (b *BSplineCurve) EvalDerivative(u float64) Vec {
	return b.EvalDerivatives(u, 1)[1]
}

// multiplicity returns how many times value v appears in the knot vector.
func multiplicity(U []float64, v float64) int {
	s := 0
	for _, k := range U {
		if k == v {
			s++
		}
	}
	return s
}

// InsertKnot returns a new curve with the knot u inserted r times using Boehm's
// algorithm (Algorithm A5.1). The geometry is unchanged; only the control-point
// representation is refined. It requires r >= 1 and r + (existing multiplicity)
// <= degree.
func (b *BSplineCurve) InsertKnot(u float64, r int) (*BSplineCurve, error) {
	if r < 1 {
		return nil, ErrParam
	}
	p := b.degree
	UP := b.knots
	Pw := b.ctrl
	np := len(Pw) - 1
	lo, hi := b.Domain()
	if u < lo || u > hi {
		return nil, ErrParam
	}
	s := multiplicity(UP, u)
	if s+r > p {
		return nil, ErrParam
	}
	k := FindSpan(np, p, u, UP)
	mp := np + p + 1
	nq := np + r
	UQ := make([]float64, mp+r+1)
	for i := 0; i <= k; i++ {
		UQ[i] = UP[i]
	}
	for i := 1; i <= r; i++ {
		UQ[k+i] = u
	}
	for i := k + 1; i <= mp; i++ {
		UQ[i+r] = UP[i]
	}
	Qw := make([]Vec, nq+1)
	for i := 0; i <= k-p; i++ {
		Qw[i] = Pw[i].Clone()
	}
	for i := k - s; i <= np; i++ {
		Qw[i+r] = Pw[i].Clone()
	}
	Rw := make([]Vec, p+1)
	for i := 0; i <= p-s; i++ {
		Rw[i] = Pw[k-p+i].Clone()
	}
	var L int
	for j := 1; j <= r; j++ {
		L = k - p + j
		for i := 0; i <= p-j-s; i++ {
			alpha := (u - UP[L+i]) / (UP[i+k+1] - UP[L+i])
			Rw[i] = Rw[i].Scale(1 - alpha).Add(Rw[i+1].Scale(alpha))
		}
		Qw[L] = Rw[0].Clone()
		Qw[k+r-j-s] = Rw[p-j-s].Clone()
	}
	for i := L + 1; i < k-s; i++ {
		Qw[i] = Rw[i-L].Clone()
	}
	return &BSplineCurve{ctrl: Qw, knots: UQ, degree: p}, nil
}
