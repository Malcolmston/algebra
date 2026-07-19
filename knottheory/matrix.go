package knottheory

// LaurentMatrix is a dense matrix whose entries are Laurent polynomials. It is
// used for the Burau representation of a braid and for the Alexander matrix of a
// diagram. Entries are stored in row-major order.
type LaurentMatrix struct {
	rows, cols int
	data       [][]Laurent
}

// NewLaurentMatrix returns an r-by-c matrix with every entry the zero
// polynomial.
func NewLaurentMatrix(r, c int) *LaurentMatrix {
	d := make([][]Laurent, r)
	for i := range d {
		d[i] = make([]Laurent, c)
		for j := range d[i] {
			d[i][j] = ZeroLaurent()
		}
	}
	return &LaurentMatrix{rows: r, cols: c, data: d}
}

// IdentityLaurentMatrix returns the n-by-n identity matrix over the Laurent
// ring.
func IdentityLaurentMatrix(n int) *LaurentMatrix {
	m := NewLaurentMatrix(n, n)
	for i := 0; i < n; i++ {
		m.data[i][i] = OneLaurent()
	}
	return m
}

// Rows returns the number of rows of the matrix.
func (m *LaurentMatrix) Rows() int { return m.rows }

// Cols returns the number of columns of the matrix.
func (m *LaurentMatrix) Cols() int { return m.cols }

// At returns the entry in row i, column j.
func (m *LaurentMatrix) At(i, j int) Laurent { return m.data[i][j] }

// Set stores v in row i, column j.
func (m *LaurentMatrix) Set(i, j int, v Laurent) { m.data[i][j] = v }

// Clone returns an independent copy of the matrix.
func (m *LaurentMatrix) Clone() *LaurentMatrix {
	c := NewLaurentMatrix(m.rows, m.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			c.data[i][j] = m.data[i][j].Clone()
		}
	}
	return c
}

// Mul returns the matrix product m*other. It panics if the shapes are
// incompatible.
func (m *LaurentMatrix) Mul(other *LaurentMatrix) *LaurentMatrix {
	if m.cols != other.rows {
		panic("knottheory: LaurentMatrix.Mul shape mismatch")
	}
	out := NewLaurentMatrix(m.rows, other.cols)
	for i := 0; i < m.rows; i++ {
		for k := 0; k < m.cols; k++ {
			a := m.data[i][k]
			if a.IsZero() {
				continue
			}
			for j := 0; j < other.cols; j++ {
				out.data[i][j] = out.data[i][j].Add(a.Mul(other.data[k][j]))
			}
		}
	}
	return out
}

// Sub returns the entrywise difference m-other.
func (m *LaurentMatrix) Sub(other *LaurentMatrix) *LaurentMatrix {
	out := NewLaurentMatrix(m.rows, m.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			out.data[i][j] = m.data[i][j].Sub(other.data[i][j])
		}
	}
	return out
}

// Determinant returns the determinant of a square matrix computed by Laplace
// (cofactor) expansion over the Laurent ring. It panics for a non-square
// matrix. The method is intended for the small matrices that arise from braids
// and diagrams of modest crossing number.
func (m *LaurentMatrix) Determinant() Laurent {
	if m.rows != m.cols {
		panic("knottheory: Determinant of a non-square matrix")
	}
	return laurentDet(m.data)
}

// laurentDet computes the determinant of the square slice matrix a by cofactor
// expansion along the first row.
func laurentDet(a [][]Laurent) Laurent {
	n := len(a)
	switch n {
	case 0:
		return OneLaurent()
	case 1:
		return a[0][0].Clone()
	case 2:
		return a[0][0].Mul(a[1][1]).Sub(a[0][1].Mul(a[1][0]))
	}
	det := ZeroLaurent()
	for j := 0; j < n; j++ {
		if a[0][j].IsZero() {
			continue
		}
		minor := make([][]Laurent, n-1)
		for i := 1; i < n; i++ {
			row := make([]Laurent, 0, n-1)
			for k := 0; k < n; k++ {
				if k == j {
					continue
				}
				row = append(row, a[i][k])
			}
			minor[i-1] = row
		}
		term := a[0][j].Mul(laurentDet(minor))
		if j%2 == 0 {
			det = det.Add(term)
		} else {
			det = det.Sub(term)
		}
	}
	return det
}

// minorDet returns the determinant of the submatrix of a obtained by deleting
// the given row and column.
func minorDet(a [][]Laurent, delRow, delCol int) Laurent {
	n := len(a)
	sub := make([][]Laurent, 0, n-1)
	for i := 0; i < n; i++ {
		if i == delRow {
			continue
		}
		row := make([]Laurent, 0, n-1)
		for j := 0; j < n; j++ {
			if j == delCol {
				continue
			}
			row = append(row, a[i][j])
		}
		sub = append(sub, row)
	}
	return laurentDet(sub)
}
