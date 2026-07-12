package matrix

import "github.com/malcolmston/algebra"

// Add returns the entrywise sum m+n, each entry simplified. It returns
// [ErrDimension] if the shapes differ.
func (m *Matrix) Add(n *Matrix) (*Matrix, error) {
	if m.rows != n.rows || m.cols != n.cols {
		return nil, ErrDimension
	}
	out := New(m.rows, m.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			out.data[i][j] = simp(algebra.Add(m.data[i][j], n.data[i][j]))
		}
	}
	return out, nil
}

// Sub returns the entrywise difference m-n, each entry simplified. It returns
// [ErrDimension] if the shapes differ.
func (m *Matrix) Sub(n *Matrix) (*Matrix, error) {
	if m.rows != n.rows || m.cols != n.cols {
		return nil, ErrDimension
	}
	out := New(m.rows, m.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			out.data[i][j] = simp(algebra.Add(m.data[i][j], algebra.Mul(algebra.Int(-1), n.data[i][j])))
		}
	}
	return out, nil
}

// Neg returns -m (every entry negated).
func (m *Matrix) Neg() *Matrix { return m.ScalarMul(algebra.Int(-1)) }

// ScalarMul returns m with every entry multiplied by the expression s.
func (m *Matrix) ScalarMul(s algebra.Expr) *Matrix {
	out := New(m.rows, m.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			out.data[i][j] = simp(algebra.Mul(s, m.data[i][j]))
		}
	}
	return out
}

// Mul returns the matrix product m·n, each entry simplified. It returns
// [ErrDimension] if the inner dimensions disagree (m.Cols() != n.Rows()).
func (m *Matrix) Mul(n *Matrix) (*Matrix, error) {
	if m.cols != n.rows {
		return nil, ErrDimension
	}
	out := New(m.rows, n.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < n.cols; j++ {
			terms := make([]algebra.Expr, m.cols)
			for k := 0; k < m.cols; k++ {
				terms[k] = algebra.Mul(m.data[i][k], n.data[k][j])
			}
			out.data[i][j] = simp(algebra.Add(terms...))
		}
	}
	return out, nil
}

// MulVec returns the matrix-vector product m·v as a [Vector], treating v as a
// column vector. It returns [ErrDimension] if m.Cols() != v.Len().
func (m *Matrix) MulVec(v *Vector) (*Vector, error) {
	if m.cols != len(v.data) {
		return nil, ErrDimension
	}
	out := make([]algebra.Expr, m.rows)
	for i := 0; i < m.rows; i++ {
		terms := make([]algebra.Expr, m.cols)
		for k := 0; k < m.cols; k++ {
			terms[k] = algebra.Mul(m.data[i][k], v.data[k])
		}
		out[i] = simp(algebra.Add(terms...))
	}
	return &Vector{data: out}, nil
}

// Transpose returns the transpose of m.
func (m *Matrix) Transpose() *Matrix {
	out := New(m.cols, m.rows)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			out.data[j][i] = m.data[i][j]
		}
	}
	return out
}

// Pow returns m raised to the non-negative integer power p. Pow(0) is the
// identity of the same size. It returns [ErrNotSquare] for a non-square matrix
// and [ErrUnsupported] for a negative exponent (use [Matrix.Inverse] and then
// Pow for negative powers).
func (m *Matrix) Pow(p int) (*Matrix, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	if p < 0 {
		return nil, ErrUnsupported
	}
	result := Identity(m.rows)
	base := m.Clone()
	// Exponentiation by squaring.
	for p > 0 {
		if p&1 == 1 {
			var err error
			result, err = result.Mul(base)
			if err != nil {
				return nil, err
			}
		}
		p >>= 1
		if p > 0 {
			var err error
			base, err = base.Mul(base)
			if err != nil {
				return nil, err
			}
		}
	}
	return result, nil
}

// Kron returns the Kronecker product m⊗n, a block matrix of shape
// (m.Rows()·n.Rows())×(m.Cols()·n.Cols()) whose (i,j) block is m[i,j]·n.
func (m *Matrix) Kron(n *Matrix) *Matrix {
	out := New(m.rows*n.rows, m.cols*n.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			for p := 0; p < n.rows; p++ {
				for q := 0; q < n.cols; q++ {
					out.data[i*n.rows+p][j*n.cols+q] = simp(algebra.Mul(m.data[i][j], n.data[p][q]))
				}
			}
		}
	}
	return out
}
