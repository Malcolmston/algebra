package tropical

// PrincipalSolution returns the principal solution of the tropical linear
// system A (*) x = b obtained by residuation. For max-plus it is the greatest x
// with A (*) x <= b (componentwise); for min-plus it is the least x with
// A (*) x >= b. Its j-th entry is the dual tropical sum over i of the residual
// b[i] (/) A[i][j]. It returns an error if the dimensions or semirings do not
// match.
func (m Matrix) PrincipalSolution(b Vector) (Vector, error) {
	if m.sr.kind != b.sr.kind || m.rows != len(b.data) {
		return Vector{}, ErrDim
	}
	dual := m.sr.Dual()
	out := make([]float64, m.cols)
	for j := 0; j < m.cols; j++ {
		acc := dual.Zero()
		for i := 0; i < m.rows; i++ {
			acc = dual.Add(acc, m.sr.Div(b.data[i], m.data[i][j]))
		}
		out[j] = acc
	}
	return Vector{data: out, sr: m.sr}, nil
}

// SolveExact returns the principal solution x of A (*) x = b together with a
// boolean reporting whether it solves the system exactly to within tol. When
// the boolean is true, x is the greatest exact solution for max-plus and the
// least exact solution for min-plus. It returns an error if the dimensions or
// semirings do not match.
func (m Matrix) SolveExact(b Vector, tol float64) (Vector, bool, error) {
	x, err := m.PrincipalSolution(b)
	if err != nil {
		return Vector{}, false, err
	}
	ax, err := m.MulVec(x)
	if err != nil {
		return Vector{}, false, err
	}
	return x, ax.EqualTol(b, tol), nil
}

// GreatestSubSolution returns the greatest x with A (*) x <= b componentwise.
// This is the residuation-based principal solution and is defined for both
// semirings (for min-plus the inequality is the natural one on the reversed
// order). It returns an error if the dimensions or semirings do not match.
func (m Matrix) GreatestSubSolution(b Vector) (Vector, error) {
	return m.PrincipalSolution(b)
}

// SolveAffine returns the least solution of the affine tropical fixed-point
// equation x = A (*) x (+) b, which equals A* (*) b where A* is the Kleene
// star. It returns ErrDivergent when the star diverges and an error if the
// dimensions or semirings do not match.
func (m Matrix) SolveAffine(b Vector) (Vector, error) {
	if m.sr.kind != b.sr.kind || m.rows != len(b.data) {
		return Vector{}, ErrDim
	}
	if !m.IsSquare() {
		return Vector{}, ErrNotSquare
	}
	star, err := m.Star()
	if err != nil {
		return Vector{}, err
	}
	return star.MulVec(b)
}

// LeftResidual returns the greatest matrix X with A (*) X <= B, written A \ B.
// The (k,j) entry is the dual tropical sum over i of B[i][j] (/) A[i][k]. It
// requires A and B to have the same number of rows and semiring and returns an
// error otherwise.
func (m Matrix) LeftResidual(b Matrix) (Matrix, error) {
	if m.sr.kind != b.sr.kind || m.rows != b.rows {
		return Matrix{}, ErrDim
	}
	dual := m.sr.Dual()
	out := Zeros(m.sr, m.cols, b.cols)
	for k := 0; k < m.cols; k++ {
		for j := 0; j < b.cols; j++ {
			acc := dual.Zero()
			for i := 0; i < m.rows; i++ {
				acc = dual.Add(acc, m.sr.Div(b.data[i][j], m.data[i][k]))
			}
			out.data[k][j] = acc
		}
	}
	return out, nil
}

// RightResidual returns the greatest matrix X with X (*) A <= B, written B / A.
// The (i,k) entry is the dual tropical sum over j of B[i][j] (/) A[k][j]. It
// requires A and B to have the same number of columns and semiring and returns
// an error otherwise.
func (m Matrix) RightResidual(b Matrix) (Matrix, error) {
	if m.sr.kind != b.sr.kind || m.cols != b.cols {
		return Matrix{}, ErrDim
	}
	dual := m.sr.Dual()
	out := Zeros(m.sr, b.rows, m.rows)
	for i := 0; i < b.rows; i++ {
		for k := 0; k < m.rows; k++ {
			acc := dual.Zero()
			for j := 0; j < m.cols; j++ {
				acc = dual.Add(acc, m.sr.Div(b.data[i][j], m.data[k][j]))
			}
			out.data[i][k] = acc
		}
	}
	return out, nil
}
