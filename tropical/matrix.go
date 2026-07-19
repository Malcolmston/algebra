package tropical

import (
	"errors"
	"strings"
)

// Matrix is a dense tropical matrix stored in row-major order together with the
// semiring under which its entries are interpreted.
type Matrix struct {
	rows, cols int
	data       [][]float64
	sr         Semiring
}

// ErrDim is returned when two matrices or vectors have incompatible shapes.
var ErrDim = errors.New("tropical: incompatible dimensions")

// ErrNotSquare is returned by operations that require a square matrix.
var ErrNotSquare = errors.New("tropical: matrix must be square")

// ErrDivergent is returned when a tropical closure diverges because of a
// negative cycle (min-plus) or a positive cycle (max-plus).
var ErrDivergent = errors.New("tropical: closure diverges (bad cycle)")

// NewMatrix returns a Matrix over sr holding a deep copy of the rectangular
// slice data. It panics if data is not rectangular.
func NewMatrix(sr Semiring, data [][]float64) Matrix {
	r := len(data)
	c := 0
	if r > 0 {
		c = len(data[0])
	}
	d := make([][]float64, r)
	for i := range data {
		if len(data[i]) != c {
			panic("tropical: NewMatrix requires a rectangular slice")
		}
		row := make([]float64, c)
		copy(row, data[i])
		d[i] = row
	}
	return Matrix{rows: r, cols: c, data: d, sr: sr}
}

// MinPlusMatrix returns a min-plus Matrix holding a copy of data.
func MinPlusMatrix(data [][]float64) Matrix { return NewMatrix(MinPlusSemiring(), data) }

// MaxPlusMatrix returns a max-plus Matrix holding a copy of data.
func MaxPlusMatrix(data [][]float64) Matrix { return NewMatrix(MaxPlusSemiring(), data) }

// Zeros returns an r-by-c matrix whose entries are all the tropical zero of the
// semiring. It panics if r or c is negative.
func Zeros(sr Semiring, r, c int) Matrix {
	if r < 0 || c < 0 {
		panic("tropical: Zeros requires non-negative dimensions")
	}
	z := sr.Zero()
	d := make([][]float64, r)
	for i := range d {
		row := make([]float64, c)
		for j := range row {
			row[j] = z
		}
		d[i] = row
	}
	return Matrix{rows: r, cols: c, data: d, sr: sr}
}

// Constant returns an r-by-c matrix with every entry equal to value.
func Constant(sr Semiring, r, c int, value float64) Matrix {
	d := make([][]float64, r)
	for i := range d {
		row := make([]float64, c)
		for j := range row {
			row[j] = value
		}
		d[i] = row
	}
	return Matrix{rows: r, cols: c, data: d, sr: sr}
}

// Identity returns the n-by-n tropical identity matrix: the tropical one (0) on
// the diagonal and the tropical zero elsewhere.
func Identity(sr Semiring, n int) Matrix {
	m := Zeros(sr, n, n)
	for i := 0; i < n; i++ {
		m.data[i][i] = 0
	}
	return m
}

// Diag returns a square matrix with the given values on the diagonal and the
// tropical zero elsewhere.
func Diag(sr Semiring, values []float64) Matrix {
	n := len(values)
	m := Zeros(sr, n, n)
	for i := 0; i < n; i++ {
		m.data[i][i] = values[i]
	}
	return m
}

// Rows returns the number of rows.
func (m Matrix) Rows() int { return m.rows }

// Cols returns the number of columns.
func (m Matrix) Cols() int { return m.cols }

// Semiring returns the semiring under which the matrix is interpreted.
func (m Matrix) Semiring() Semiring { return m.sr }

// At returns the entry in row i and column j. It panics if the indices are out
// of range.
func (m Matrix) At(i, j int) float64 { return m.data[i][j] }

// Set stores value in row i and column j. It panics if the indices are out of
// range.
func (m Matrix) Set(i, j int, value float64) { m.data[i][j] = value }

// Clone returns an independent deep copy of the matrix.
func (m Matrix) Clone() Matrix {
	return NewMatrix(m.sr, m.data)
}

// Raw returns a fresh deep copy of the entries as a [][]float64.
func (m Matrix) Raw() [][]float64 {
	d := make([][]float64, m.rows)
	for i := range d {
		row := make([]float64, m.cols)
		copy(row, m.data[i])
		d[i] = row
	}
	return d
}

// Row returns a fresh copy of row i.
func (m Matrix) Row(i int) []float64 {
	c := make([]float64, m.cols)
	copy(c, m.data[i])
	return c
}

// Col returns a fresh copy of column j.
func (m Matrix) Col(j int) []float64 {
	c := make([]float64, m.rows)
	for i := 0; i < m.rows; i++ {
		c[i] = m.data[i][j]
	}
	return c
}

// RowVector returns row i as a Vector.
func (m Matrix) RowVector(i int) Vector { return Vector{data: m.Row(i), sr: m.sr} }

// ColVector returns column j as a Vector.
func (m Matrix) ColVector(j int) Vector { return Vector{data: m.Col(j), sr: m.sr} }

// Diagonal returns a fresh copy of the main diagonal. For a non-square matrix
// it returns the leading square diagonal of length min(rows, cols).
func (m Matrix) Diagonal() []float64 {
	n := m.rows
	if m.cols < n {
		n = m.cols
	}
	d := make([]float64, n)
	for i := 0; i < n; i++ {
		d[i] = m.data[i][i]
	}
	return d
}

// IsSquare reports whether the matrix has equally many rows and columns.
func (m Matrix) IsSquare() bool { return m.rows == m.cols }

// Equal reports whether m and n have the same semiring, shape and identical
// entries.
func (m Matrix) Equal(n Matrix) bool {
	if m.sr.kind != n.sr.kind || m.rows != n.rows || m.cols != n.cols {
		return false
	}
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			if m.data[i][j] != n.data[i][j] {
				return false
			}
		}
	}
	return true
}

// EqualTol reports whether m and n have the same semiring and shape and every
// pair of entries agrees to within tol (infinities must match exactly).
func (m Matrix) EqualTol(n Matrix, tol float64) bool {
	if m.sr.kind != n.sr.kind || m.rows != n.rows || m.cols != n.cols {
		return false
	}
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			if !closeScalar(m.data[i][j], n.data[i][j], tol) {
				return false
			}
		}
	}
	return true
}

// Transpose returns the transpose of the matrix.
func (m Matrix) Transpose() Matrix {
	out := Zeros(m.sr, m.cols, m.rows)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			out.data[j][i] = m.data[i][j]
		}
	}
	return out
}

// SubMatrix returns the submatrix spanning rows [r0, r1) and columns [c0, c1).
// It panics if the ranges are invalid.
func (m Matrix) SubMatrix(r0, r1, c0, c1 int) Matrix {
	if r0 < 0 || c0 < 0 || r1 > m.rows || c1 > m.cols || r0 > r1 || c0 > c1 {
		panic("tropical: SubMatrix range out of bounds")
	}
	out := Zeros(m.sr, r1-r0, c1-c0)
	for i := r0; i < r1; i++ {
		for j := c0; j < c1; j++ {
			out.data[i-r0][j-c0] = m.data[i][j]
		}
	}
	return out
}

// Minor returns the matrix with row p and column q removed. It panics if the
// indices are out of range.
func (m Matrix) Minor(p, q int) Matrix {
	if p < 0 || p >= m.rows || q < 0 || q >= m.cols {
		panic("tropical: Minor index out of range")
	}
	out := Zeros(m.sr, m.rows-1, m.cols-1)
	ri := 0
	for i := 0; i < m.rows; i++ {
		if i == p {
			continue
		}
		cj := 0
		for j := 0; j < m.cols; j++ {
			if j == q {
				continue
			}
			out.data[ri][cj] = m.data[i][j]
			cj++
		}
		ri++
	}
	return out
}

// Add returns the elementwise tropical sum of m and n. It returns an error if
// the shapes or semirings differ.
func (m Matrix) Add(n Matrix) (Matrix, error) {
	if m.sr.kind != n.sr.kind || m.rows != n.rows || m.cols != n.cols {
		return Matrix{}, ErrDim
	}
	out := Zeros(m.sr, m.rows, m.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			out.data[i][j] = m.sr.Add(m.data[i][j], n.data[i][j])
		}
	}
	return out, nil
}

// Mul returns the tropical matrix product m (*) n, where the (i,j) entry is the
// tropical sum over k of m[i][k] (*) n[k][j]. Over min-plus this is the
// shortest-path relaxation. It returns an error if the inner dimensions or the
// semirings do not match.
func (m Matrix) Mul(n Matrix) (Matrix, error) {
	if m.sr.kind != n.sr.kind || m.cols != n.rows {
		return Matrix{}, ErrDim
	}
	out := Zeros(m.sr, m.rows, n.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < n.cols; j++ {
			acc := m.sr.Zero()
			for k := 0; k < m.cols; k++ {
				acc = m.sr.Add(acc, m.sr.Mul(m.data[i][k], n.data[k][j]))
			}
			out.data[i][j] = acc
		}
	}
	return out, nil
}

// ScalarMul returns the matrix obtained by tropically multiplying every entry
// by c.
func (m Matrix) ScalarMul(c float64) Matrix {
	out := Zeros(m.sr, m.rows, m.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			out.data[i][j] = m.sr.Mul(m.data[i][j], c)
		}
	}
	return out
}

// ScalarAdd returns the matrix obtained by tropically adding c to every entry.
func (m Matrix) ScalarAdd(c float64) Matrix {
	out := Zeros(m.sr, m.rows, m.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			out.data[i][j] = m.sr.Add(m.data[i][j], c)
		}
	}
	return out
}

// MulVec returns the tropical matrix-vector product m (*) v. It returns an
// error if the dimensions or semirings do not match.
func (m Matrix) MulVec(v Vector) (Vector, error) {
	if m.sr.kind != v.sr.kind || m.cols != len(v.data) {
		return Vector{}, ErrDim
	}
	out := make([]float64, m.rows)
	for i := 0; i < m.rows; i++ {
		acc := m.sr.Zero()
		for k := 0; k < m.cols; k++ {
			acc = m.sr.Add(acc, m.sr.Mul(m.data[i][k], v.data[k]))
		}
		out[i] = acc
	}
	return Vector{data: out, sr: m.sr}, nil
}

// VecMul returns the tropical row-vector-matrix product v^T (*) m. It returns an
// error if the dimensions or semirings do not match.
func (m Matrix) VecMul(v Vector) (Vector, error) {
	if m.sr.kind != v.sr.kind || m.rows != len(v.data) {
		return Vector{}, ErrDim
	}
	out := make([]float64, m.cols)
	for j := 0; j < m.cols; j++ {
		acc := m.sr.Zero()
		for k := 0; k < m.rows; k++ {
			acc = m.sr.Add(acc, m.sr.Mul(v.data[k], m.data[k][j]))
		}
		out[j] = acc
	}
	return Vector{data: out, sr: m.sr}, nil
}

// Pow returns the tropical matrix power m^n by repeated squaring. The zeroth
// power is the identity. It returns ErrNotSquare for a non-square matrix and an
// error for negative n.
func (m Matrix) Pow(n int) (Matrix, error) {
	if !m.IsSquare() {
		return Matrix{}, ErrNotSquare
	}
	if n < 0 {
		return Matrix{}, errors.New("tropical: Pow requires n >= 0")
	}
	result := Identity(m.sr, m.rows)
	base := m.Clone()
	for n > 0 {
		if n&1 == 1 {
			result, _ = result.Mul(base)
		}
		n >>= 1
		if n > 0 {
			base, _ = base.Mul(base)
		}
	}
	return result, nil
}

// Trace returns the tropical trace, the tropical sum of the diagonal entries.
// It returns ErrNotSquare for a non-square matrix.
func (m Matrix) Trace() (float64, error) {
	if !m.IsSquare() {
		return 0, ErrNotSquare
	}
	acc := m.sr.Zero()
	for i := 0; i < m.rows; i++ {
		acc = m.sr.Add(acc, m.data[i][i])
	}
	return acc, nil
}

// IsIdentity reports whether the matrix is square and equals the tropical
// identity to within tol.
func (m Matrix) IsIdentity(tol float64) bool {
	if !m.IsSquare() {
		return false
	}
	z := m.sr.Zero()
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			if i == j {
				if !closeScalar(m.data[i][j], 0, tol) {
					return false
				}
			} else if m.data[i][j] != z {
				return false
			}
		}
	}
	return true
}

// IsDiagonal reports whether every off-diagonal entry is the tropical zero.
func (m Matrix) IsDiagonal() bool {
	z := m.sr.Zero()
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			if i != j && m.data[i][j] != z {
				return false
			}
		}
	}
	return true
}

// IsZero reports whether every entry is the tropical zero.
func (m Matrix) IsZero() bool {
	z := m.sr.Zero()
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			if m.data[i][j] != z {
				return false
			}
		}
	}
	return true
}

// String renders the matrix with one row per line and space-separated entries.
func (m Matrix) String() string {
	var b strings.Builder
	for i := 0; i < m.rows; i++ {
		parts := make([]string, m.cols)
		for j := 0; j < m.cols; j++ {
			parts[j] = m.sr.FormatScalar(m.data[i][j])
		}
		b.WriteString("[" + strings.Join(parts, " ") + "]")
		if i != m.rows-1 {
			b.WriteByte('\n')
		}
	}
	return b.String()
}
