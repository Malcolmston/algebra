package matrix

import "github.com/malcolmston/algebra"

// isZeroExpr reports whether e simplifies to the integer 0.
func isZeroExpr(e algebra.Expr) bool { return simp(e).Equal(zero()) }

// Trace returns the sum of the main-diagonal entries. It returns [ErrNotSquare]
// for a non-square matrix.
func (m *Matrix) Trace() (algebra.Expr, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	terms := make([]algebra.Expr, m.rows)
	for i := 0; i < m.rows; i++ {
		terms[i] = m.data[i][i]
	}
	return simp(algebra.Add(terms...)), nil
}

// Minor returns the (i,j) minor: the determinant of the submatrix obtained by
// deleting row i and column j. It returns [ErrNotSquare] for a non-square
// matrix and panics if the indices are out of range.
func (m *Matrix) Minor(i, j int) (algebra.Expr, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	m.checkIndex(i, j)
	return m.submatrix(i, j).Det()
}

// Cofactor returns the (i,j) cofactor: (-1)^(i+j) times the (i,j) minor.
func (m *Matrix) Cofactor(i, j int) (algebra.Expr, error) {
	minor, err := m.Minor(i, j)
	if err != nil {
		return nil, err
	}
	if (i+j)%2 == 1 {
		return simp(algebra.Mul(algebra.Int(-1), minor)), nil
	}
	return minor, nil
}

// submatrix returns the matrix with row r and column c removed.
func (m *Matrix) submatrix(r, c int) *Matrix {
	out := New(m.rows-1, m.cols-1)
	ri := 0
	for i := 0; i < m.rows; i++ {
		if i == r {
			continue
		}
		cj := 0
		for j := 0; j < m.cols; j++ {
			if j == c {
				continue
			}
			out.data[ri][cj] = m.data[i][j]
			cj++
		}
		ri++
	}
	return out
}

// Det returns the exact determinant computed by cofactor (Laplace) expansion,
// which performs no division and therefore stays exact for symbolic entries. It
// returns [ErrNotSquare] for a non-square matrix. The 0×0 determinant is 1.
func (m *Matrix) Det() (algebra.Expr, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	return simp(m.det()), nil
}

// det computes the determinant without the outer simplify wrapper, with fast
// paths for the common small sizes.
func (m *Matrix) det() algebra.Expr {
	n := m.rows
	d := m.data
	switch n {
	case 0:
		return one()
	case 1:
		return d[0][0]
	case 2:
		return algebra.Add(
			algebra.Mul(d[0][0], d[1][1]),
			algebra.Mul(algebra.Int(-1), d[0][1], d[1][0]),
		)
	case 3:
		return algebra.Add(
			algebra.Mul(d[0][0], d[1][1], d[2][2]),
			algebra.Mul(d[0][1], d[1][2], d[2][0]),
			algebra.Mul(d[0][2], d[1][0], d[2][1]),
			algebra.Mul(algebra.Int(-1), d[0][2], d[1][1], d[2][0]),
			algebra.Mul(algebra.Int(-1), d[0][0], d[1][2], d[2][1]),
			algebra.Mul(algebra.Int(-1), d[0][1], d[1][0], d[2][2]),
		)
	}
	// Laplace expansion along the first row, skipping zero cofactor entries.
	terms := make([]algebra.Expr, 0, n)
	for j := 0; j < n; j++ {
		if isZeroExpr(d[0][j]) {
			continue
		}
		minor := m.submatrix(0, j).det()
		term := algebra.Mul(d[0][j], minor)
		if j%2 == 1 {
			term = algebra.Mul(algebra.Int(-1), term)
		}
		terms = append(terms, term)
	}
	return algebra.Add(terms...)
}

// Adjugate returns the adjugate (classical adjoint): the transpose of the
// cofactor matrix. It returns [ErrNotSquare] for a non-square matrix.
func (m *Matrix) Adjugate() (*Matrix, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	out := New(m.rows, m.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			cof, err := m.Cofactor(i, j)
			if err != nil {
				return nil, err
			}
			// Transpose: cofactor (i,j) goes to position (j,i).
			out.data[j][i] = cof
		}
	}
	return out, nil
}

// Inverse returns the exact inverse computed as adjugate(m)/det(m). It returns
// [ErrNotSquare] for a non-square matrix and [ErrSingular] if the determinant
// simplifies to 0.
func (m *Matrix) Inverse() (*Matrix, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	det, err := m.Det()
	if err != nil {
		return nil, err
	}
	if isZeroExpr(det) {
		return nil, ErrSingular
	}
	adj, err := m.Adjugate()
	if err != nil {
		return nil, err
	}
	invDet := algebra.Pow(det, algebra.Int(-1))
	return adj.ScalarMul(invDet), nil
}

// RREF returns the reduced row-echelon form of m together with the rank (the
// number of nonzero pivot rows). Pivots are exact: a column entry is treated as
// a pivot only when it does not simplify to 0, and each pivot row is scaled by
// the inverse of its pivot. Symbolic pivots produce rational-function entries.
func (m *Matrix) RREF() (*Matrix, int) {
	out := m.Clone()
	lead := 0
	rank := 0
	for r := 0; r < out.rows; r++ {
		if lead >= out.cols {
			break
		}
		// Find a row at or below r whose entry in column lead is nonzero.
		i := r
		for isZeroExpr(out.data[i][lead]) {
			i++
			if i == out.rows {
				i = r
				lead++
				if lead == out.cols {
					return out, rank
				}
			}
		}
		out.data[i], out.data[r] = out.data[r], out.data[i]
		// Scale pivot row so the pivot becomes 1.
		pivInv := algebra.Pow(out.data[r][lead], algebra.Int(-1))
		for j := 0; j < out.cols; j++ {
			out.data[r][j] = simp(algebra.Mul(pivInv, out.data[r][j]))
		}
		// Eliminate the lead column from every other row.
		for k := 0; k < out.rows; k++ {
			if k == r {
				continue
			}
			factor := out.data[k][lead]
			if isZeroExpr(factor) {
				continue
			}
			for j := 0; j < out.cols; j++ {
				out.data[k][j] = simp(algebra.Add(
					out.data[k][j],
					algebra.Mul(algebra.Int(-1), factor, out.data[r][j]),
				))
			}
		}
		rank++
		lead++
	}
	return out, rank
}

// Rank returns the rank of m, the number of linearly independent rows,
// determined from its reduced row-echelon form.
func (m *Matrix) Rank() int {
	_, r := m.RREF()
	return r
}

// Floats converts every entry to a float64, returning [ErrUnsupported] wrapping
// the underlying evaluation error if any entry contains a free symbol or cannot
// be evaluated numerically. This is the numeric fast-path escape hatch.
func (m *Matrix) Floats() ([][]float64, error) {
	out := make([][]float64, m.rows)
	for i := 0; i < m.rows; i++ {
		out[i] = make([]float64, m.cols)
		for j := 0; j < m.cols; j++ {
			f, err := algebra.Evalf(m.data[i][j])
			if err != nil {
				return nil, err
			}
			out[i][j] = f
		}
	}
	return out, nil
}
