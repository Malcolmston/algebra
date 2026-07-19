package splines

// checkNet validates that ctrl is a non-empty rectangular control net whose
// points share a dimension and returns the row count, column count and
// dimension.
func checkNet(ctrl [][]Vec) (rows, cols, dim int, err error) {
	rows = len(ctrl)
	if rows == 0 {
		return 0, 0, 0, ErrEmpty
	}
	cols = len(ctrl[0])
	if cols == 0 {
		return 0, 0, 0, ErrEmpty
	}
	dim = ctrl[0][0].Dim()
	if dim == 0 {
		return 0, 0, 0, ErrDim
	}
	for i := range ctrl {
		if len(ctrl[i]) != cols {
			return 0, 0, 0, ErrLenMismatch
		}
		for j := range ctrl[i] {
			if ctrl[i][j].Dim() != dim {
				return 0, 0, 0, ErrDim
			}
		}
	}
	return rows, cols, dim, nil
}

func cloneNet(ctrl [][]Vec) [][]Vec {
	out := make([][]Vec, len(ctrl))
	for i := range ctrl {
		out[i] = CloneVecs(ctrl[i])
	}
	return out
}

// BezierSurface is a tensor-product Bezier surface defined by a rectangular
// control net. Its degrees are one less than the numbers of rows and columns.
type BezierSurface struct {
	ctrl [][]Vec // ctrl[i][j], i indexes u, j indexes v
	dim  int
}

// NewBezierSurface builds a tensor-product Bezier surface from a rectangular
// control net ctrl[i][j], where i runs over the u direction and j over the v
// direction.
func NewBezierSurface(ctrl [][]Vec) (*BezierSurface, error) {
	_, _, dim, err := checkNet(ctrl)
	if err != nil {
		return nil, err
	}
	return &BezierSurface{ctrl: cloneNet(ctrl), dim: dim}, nil
}

// DegreeU returns the surface degree in the u direction.
func (s *BezierSurface) DegreeU() int { return len(s.ctrl) - 1 }

// DegreeV returns the surface degree in the v direction.
func (s *BezierSurface) DegreeV() int { return len(s.ctrl[0]) - 1 }

// Dim returns the spatial dimension of the control net.
func (s *BezierSurface) Dim() int { return s.dim }

// Eval returns the surface point at parameters (u,v) in [0,1]x[0,1] using the
// Bernstein basis in both directions.
func (s *BezierSurface) Eval(u, v float64) Vec {
	m := s.DegreeU()
	n := s.DegreeV()
	bu := BernsteinAll(m, u)
	bv := BernsteinAll(n, v)
	out := make(Vec, s.dim)
	for i := 0; i <= m; i++ {
		for j := 0; j <= n; j++ {
			c := bu[i] * bv[j]
			for k := 0; k < s.dim; k++ {
				out[k] += c * s.ctrl[i][j][k]
			}
		}
	}
	return out
}

// derivNetU returns the control net of the u-partial derivative surface, of
// degree (m-1, n).
func (s *BezierSurface) derivNetU() [][]Vec {
	m := s.DegreeU()
	n := s.DegreeV()
	out := make([][]Vec, m)
	for i := 0; i < m; i++ {
		out[i] = make([]Vec, n+1)
		for j := 0; j <= n; j++ {
			out[i][j] = s.ctrl[i+1][j].Sub(s.ctrl[i][j]).Scale(float64(m))
		}
	}
	return out
}

// derivNetV returns the control net of the v-partial derivative surface, of
// degree (m, n-1).
func (s *BezierSurface) derivNetV() [][]Vec {
	m := s.DegreeU()
	n := s.DegreeV()
	out := make([][]Vec, m+1)
	for i := 0; i <= m; i++ {
		out[i] = make([]Vec, n)
		for j := 0; j < n; j++ {
			out[i][j] = s.ctrl[i][j+1].Sub(s.ctrl[i][j]).Scale(float64(n))
		}
	}
	return out
}

// EvalPartialU returns the partial derivative of the surface with respect to u
// at (u,v).
func (s *BezierSurface) EvalPartialU(u, v float64) Vec {
	if s.DegreeU() == 0 {
		return make(Vec, s.dim)
	}
	ds := &BezierSurface{ctrl: s.derivNetU(), dim: s.dim}
	return ds.Eval(u, v)
}

// EvalPartialV returns the partial derivative of the surface with respect to v
// at (u,v).
func (s *BezierSurface) EvalPartialV(u, v float64) Vec {
	if s.DegreeV() == 0 {
		return make(Vec, s.dim)
	}
	ds := &BezierSurface{ctrl: s.derivNetV(), dim: s.dim}
	return ds.Eval(u, v)
}

// Normal returns the (unnormalised) surface normal at (u,v) as the cross product
// of the u- and v-partial derivatives. It is defined only for three-dimensional
// surfaces.
func (s *BezierSurface) Normal(u, v float64) (Vec, error) {
	if s.dim != 3 {
		return nil, ErrDim
	}
	return Cross(s.EvalPartialU(u, v), s.EvalPartialV(u, v)), nil
}

// Cross returns the 3-D cross product a x b. It requires both vectors to be
// three-dimensional.
func Cross(a, b Vec) Vec {
	return Vec{
		a[1]*b[2] - a[2]*b[1],
		a[2]*b[0] - a[0]*b[2],
		a[0]*b[1] - a[1]*b[0],
	}
}

// BSplineSurface is a tensor-product B-spline surface with independent degrees
// and knot vectors in the u and v directions.
type BSplineSurface struct {
	ctrl   [][]Vec // ctrl[i][j]
	knotsU []float64
	knotsV []float64
	degU   int
	degV   int
	dim    int
}

// NewBSplineSurface builds a tensor-product B-spline surface. The control net
// must be rectangular with len(ctrl) control rows in u and len(ctrl[0]) columns
// in v; the knot vectors must satisfy len(knotsU) == rows+degU+1 and
// len(knotsV) == cols+degV+1.
func NewBSplineSurface(ctrl [][]Vec, knotsU, knotsV []float64, degU, degV int) (*BSplineSurface, error) {
	rows, cols, dim, err := checkNet(ctrl)
	if err != nil {
		return nil, err
	}
	if degU < 1 || degV < 1 {
		return nil, ErrDegree
	}
	if rows < degU+1 || cols < degV+1 {
		return nil, ErrTooFewPoints
	}
	if len(knotsU) != rows+degU+1 || len(knotsV) != cols+degV+1 {
		return nil, ErrKnots
	}
	if !nondecreasing(knotsU) || !nondecreasing(knotsV) {
		return nil, ErrKnots
	}
	return &BSplineSurface{
		ctrl:   cloneNet(ctrl),
		knotsU: append([]float64(nil), knotsU...),
		knotsV: append([]float64(nil), knotsV...),
		degU:   degU,
		degV:   degV,
		dim:    dim,
	}, nil
}

// DegreeU returns the surface degree in the u direction.
func (s *BSplineSurface) DegreeU() int { return s.degU }

// DegreeV returns the surface degree in the v direction.
func (s *BSplineSurface) DegreeV() int { return s.degV }

// Dim returns the spatial dimension of the control net.
func (s *BSplineSurface) Dim() int { return s.dim }

// Domain returns the valid parameter rectangle [uLo,uHi] x [vLo,vHi].
func (s *BSplineSurface) Domain() (uLo, uHi, vLo, vHi float64) {
	nu := len(s.ctrl) - 1
	nv := len(s.ctrl[0]) - 1
	return s.knotsU[s.degU], s.knotsU[nu+1], s.knotsV[s.degV], s.knotsV[nv+1]
}

// Eval returns the surface point at (u,v) using the tensor product of the
// univariate B-spline basis functions in each direction.
func (s *BSplineSurface) Eval(u, v float64) Vec {
	nu := len(s.ctrl) - 1
	nv := len(s.ctrl[0]) - 1
	spanU := FindSpan(nu, s.degU, u, s.knotsU)
	spanV := FindSpan(nv, s.degV, v, s.knotsV)
	Nu := BasisFuns(spanU, u, s.degU, s.knotsU)
	Nv := BasisFuns(spanV, v, s.degV, s.knotsV)
	out := make(Vec, s.dim)
	for i := 0; i <= s.degU; i++ {
		ii := spanU - s.degU + i
		// Accumulate the v-direction combination for this control row first.
		tmp := make(Vec, s.dim)
		for j := 0; j <= s.degV; j++ {
			jj := spanV - s.degV + j
			for k := 0; k < s.dim; k++ {
				tmp[k] += Nv[j] * s.ctrl[ii][jj][k]
			}
		}
		for k := 0; k < s.dim; k++ {
			out[k] += Nu[i] * tmp[k]
		}
	}
	return out
}

// NURBSSurface is a tensor-product rational B-spline surface, adding a positive
// weight to every control point of a [BSplineSurface].
type NURBSSurface struct {
	ctrl    [][]Vec
	weights [][]float64
	knotsU  []float64
	knotsV  []float64
	degU    int
	degV    int
	dim     int
}

// NewNURBSSurface builds a tensor-product NURBS surface. The weight net must
// match the control net in shape and every weight must be strictly positive.
func NewNURBSSurface(ctrl [][]Vec, weights [][]float64, knotsU, knotsV []float64, degU, degV int) (*NURBSSurface, error) {
	rows, cols, dim, err := checkNet(ctrl)
	if err != nil {
		return nil, err
	}
	if len(weights) != rows {
		return nil, ErrLenMismatch
	}
	for i := range weights {
		if len(weights[i]) != cols {
			return nil, ErrLenMismatch
		}
		for _, w := range weights[i] {
			if w <= 0 {
				return nil, ErrWeights
			}
		}
	}
	if degU < 1 || degV < 1 {
		return nil, ErrDegree
	}
	if rows < degU+1 || cols < degV+1 {
		return nil, ErrTooFewPoints
	}
	if len(knotsU) != rows+degU+1 || len(knotsV) != cols+degV+1 {
		return nil, ErrKnots
	}
	if !nondecreasing(knotsU) || !nondecreasing(knotsV) {
		return nil, ErrKnots
	}
	w := make([][]float64, rows)
	for i := range weights {
		w[i] = append([]float64(nil), weights[i]...)
	}
	return &NURBSSurface{
		ctrl:    cloneNet(ctrl),
		weights: w,
		knotsU:  append([]float64(nil), knotsU...),
		knotsV:  append([]float64(nil), knotsV...),
		degU:    degU,
		degV:    degV,
		dim:     dim,
	}, nil
}

// DegreeU returns the surface degree in the u direction.
func (s *NURBSSurface) DegreeU() int { return s.degU }

// DegreeV returns the surface degree in the v direction.
func (s *NURBSSurface) DegreeV() int { return s.degV }

// Dim returns the spatial dimension of the control net.
func (s *NURBSSurface) Dim() int { return s.dim }

// Domain returns the valid parameter rectangle of the surface.
func (s *NURBSSurface) Domain() (uLo, uHi, vLo, vHi float64) {
	nu := len(s.ctrl) - 1
	nv := len(s.ctrl[0]) - 1
	return s.knotsU[s.degU], s.knotsU[nu+1], s.knotsV[s.degV], s.knotsV[nv+1]
}

// Eval returns the surface point at (u,v) using the rational tensor-product
// basis (weighted sum divided by the weight sum).
func (s *NURBSSurface) Eval(u, v float64) Vec {
	nu := len(s.ctrl) - 1
	nv := len(s.ctrl[0]) - 1
	spanU := FindSpan(nu, s.degU, u, s.knotsU)
	spanV := FindSpan(nv, s.degV, v, s.knotsV)
	Nu := BasisFuns(spanU, u, s.degU, s.knotsU)
	Nv := BasisFuns(spanV, v, s.degV, s.knotsV)
	num := make(Vec, s.dim)
	var den float64
	for i := 0; i <= s.degU; i++ {
		ii := spanU - s.degU + i
		for j := 0; j <= s.degV; j++ {
			jj := spanV - s.degV + j
			wN := Nu[i] * Nv[j] * s.weights[ii][jj]
			for k := 0; k < s.dim; k++ {
				num[k] += wN * s.ctrl[ii][jj][k]
			}
			den += wN
		}
	}
	return num.Scale(1 / den)
}
