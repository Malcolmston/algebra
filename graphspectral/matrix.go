package graphspectral

import (
	"fmt"
	"math"
	"strings"
)

// Matrix is a dense matrix of float64 values stored in row-major order.
type Matrix struct {
	rows, cols int
	data       []float64
}

// NewMatrix returns a new rows-by-cols matrix with every entry zero. It panics
// if rows or cols is negative.
func NewMatrix(rows, cols int) *Matrix {
	if rows < 0 || cols < 0 {
		panic("graphspectral: negative matrix dimension")
	}
	return &Matrix{rows: rows, cols: cols, data: make([]float64, rows*cols)}
}

// NewMatrixFromRows builds a matrix from a slice of equal-length rows. It returns
// ErrDimensionMismatch if the rows have differing lengths and ErrEmpty if no
// rows are given.
func NewMatrixFromRows(rows [][]float64) (*Matrix, error) {
	if len(rows) == 0 {
		return nil, ErrEmpty
	}
	c := len(rows[0])
	m := NewMatrix(len(rows), c)
	for i, r := range rows {
		if len(r) != c {
			return nil, ErrDimensionMismatch
		}
		copy(m.data[i*c:(i+1)*c], r)
	}
	return m, nil
}

// IdentityMatrix returns the n-by-n identity matrix.
func IdentityMatrix(n int) *Matrix {
	m := NewMatrix(n, n)
	for i := 0; i < n; i++ {
		m.data[i*n+i] = 1
	}
	return m
}

// DiagMatrix returns a square diagonal matrix whose diagonal is d.
func DiagMatrix(d []float64) *Matrix {
	n := len(d)
	m := NewMatrix(n, n)
	for i := 0; i < n; i++ {
		m.data[i*n+i] = d[i]
	}
	return m
}

// Rows returns the number of rows.
func (m *Matrix) Rows() int { return m.rows }

// Cols returns the number of columns.
func (m *Matrix) Cols() int { return m.cols }

// IsSquare reports whether the matrix has equal numbers of rows and columns.
func (m *Matrix) IsSquare() bool { return m.rows == m.cols }

// At returns the entry in row i, column j. It panics if the index is out of
// range.
func (m *Matrix) At(i, j int) float64 {
	if i < 0 || i >= m.rows || j < 0 || j >= m.cols {
		panic("graphspectral: matrix index out of range")
	}
	return m.data[i*m.cols+j]
}

// Set assigns v to the entry in row i, column j. It panics if the index is out
// of range.
func (m *Matrix) Set(i, j int, v float64) {
	if i < 0 || i >= m.rows || j < 0 || j >= m.cols {
		panic("graphspectral: matrix index out of range")
	}
	m.data[i*m.cols+j] = v
}

// Add increments the entry in row i, column j by v.
func (m *Matrix) Add(i, j int, v float64) {
	m.Set(i, j, m.At(i, j)+v)
}

// Clone returns an independent copy of the matrix.
func (m *Matrix) Clone() *Matrix {
	c := &Matrix{rows: m.rows, cols: m.cols, data: make([]float64, len(m.data))}
	copy(c.data, m.data)
	return c
}

// Row returns a copy of row i.
func (m *Matrix) Row(i int) []float64 {
	out := make([]float64, m.cols)
	copy(out, m.data[i*m.cols:(i+1)*m.cols])
	return out
}

// Col returns a copy of column j.
func (m *Matrix) Col(j int) []float64 {
	out := make([]float64, m.rows)
	for i := 0; i < m.rows; i++ {
		out[i] = m.data[i*m.cols+j]
	}
	return out
}

// SetRow overwrites row i with the values in r. It returns ErrDimensionMismatch
// if len(r) differs from the number of columns.
func (m *Matrix) SetRow(i int, r []float64) error {
	if len(r) != m.cols {
		return ErrDimensionMismatch
	}
	copy(m.data[i*m.cols:(i+1)*m.cols], r)
	return nil
}

// Diagonal returns a copy of the main diagonal.
func (m *Matrix) Diagonal() []float64 {
	n := m.rows
	if m.cols < n {
		n = m.cols
	}
	out := make([]float64, n)
	for i := 0; i < n; i++ {
		out[i] = m.data[i*m.cols+i]
	}
	return out
}

// Trace returns the sum of the diagonal entries. It returns 0 for a non-square
// matrix's shorter diagonal is used.
func (m *Matrix) Trace() float64 {
	var s float64
	for _, x := range m.Diagonal() {
		s += x
	}
	return s
}

// Equals reports whether m and o have the same shape and identical entries.
func (m *Matrix) Equals(o *Matrix) bool {
	if m.rows != o.rows || m.cols != o.cols {
		return false
	}
	for i := range m.data {
		if m.data[i] != o.data[i] {
			return false
		}
	}
	return true
}

// ApproxEqual reports whether m and o have the same shape and agree entry-wise
// to within absolute tolerance tol.
func (m *Matrix) ApproxEqual(o *Matrix, tol float64) bool {
	if m.rows != o.rows || m.cols != o.cols {
		return false
	}
	for i := range m.data {
		if math.Abs(m.data[i]-o.data[i]) > tol {
			return false
		}
	}
	return true
}

// Plus returns the sum m+o. It returns ErrDimensionMismatch on a shape mismatch.
func (m *Matrix) Plus(o *Matrix) (*Matrix, error) {
	if m.rows != o.rows || m.cols != o.cols {
		return nil, ErrDimensionMismatch
	}
	r := m.Clone()
	for i := range r.data {
		r.data[i] += o.data[i]
	}
	return r, nil
}

// Minus returns the difference m-o. It returns ErrDimensionMismatch on a shape
// mismatch.
func (m *Matrix) Minus(o *Matrix) (*Matrix, error) {
	if m.rows != o.rows || m.cols != o.cols {
		return nil, ErrDimensionMismatch
	}
	r := m.Clone()
	for i := range r.data {
		r.data[i] -= o.data[i]
	}
	return r, nil
}

// Scale returns a new matrix equal to s*m.
func (m *Matrix) Scale(s float64) *Matrix {
	r := m.Clone()
	for i := range r.data {
		r.data[i] *= s
	}
	return r
}

// Mul returns the matrix product m*o. It returns ErrDimensionMismatch when the
// inner dimensions disagree.
func (m *Matrix) Mul(o *Matrix) (*Matrix, error) {
	if m.cols != o.rows {
		return nil, ErrDimensionMismatch
	}
	r := NewMatrix(m.rows, o.cols)
	for i := 0; i < m.rows; i++ {
		for k := 0; k < m.cols; k++ {
			a := m.data[i*m.cols+k]
			if a == 0 {
				continue
			}
			for j := 0; j < o.cols; j++ {
				r.data[i*r.cols+j] += a * o.data[k*o.cols+j]
			}
		}
	}
	return r, nil
}

// MulVec returns the matrix-vector product m*v. It returns ErrDimensionMismatch
// if len(v) differs from the number of columns.
func (m *Matrix) MulVec(v []float64) ([]float64, error) {
	if len(v) != m.cols {
		return nil, ErrDimensionMismatch
	}
	out := make([]float64, m.rows)
	for i := 0; i < m.rows; i++ {
		var s float64
		for j := 0; j < m.cols; j++ {
			s += m.data[i*m.cols+j] * v[j]
		}
		out[i] = s
	}
	return out, nil
}

// Transpose returns the transpose of the matrix.
func (m *Matrix) Transpose() *Matrix {
	r := NewMatrix(m.cols, m.rows)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			r.data[j*r.cols+i] = m.data[i*m.cols+j]
		}
	}
	return r
}

// IsSymmetric reports whether the matrix is square and symmetric to within
// absolute tolerance tol.
func (m *Matrix) IsSymmetric(tol float64) bool {
	if !m.IsSquare() {
		return false
	}
	for i := 0; i < m.rows; i++ {
		for j := i + 1; j < m.cols; j++ {
			if math.Abs(m.data[i*m.cols+j]-m.data[j*m.cols+i]) > tol {
				return false
			}
		}
	}
	return true
}

// Symmetrize returns (m + mᵀ)/2, the symmetric part of a square matrix.
func (m *Matrix) Symmetrize() *Matrix {
	r := NewMatrix(m.rows, m.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			r.data[i*r.cols+j] = 0.5 * (m.data[i*m.cols+j] + m.data[j*m.cols+i])
		}
	}
	return r
}

// RowSums returns the sum of each row.
func (m *Matrix) RowSums() []float64 {
	out := make([]float64, m.rows)
	for i := 0; i < m.rows; i++ {
		var s float64
		for j := 0; j < m.cols; j++ {
			s += m.data[i*m.cols+j]
		}
		out[i] = s
	}
	return out
}

// ColSums returns the sum of each column.
func (m *Matrix) ColSums() []float64 {
	out := make([]float64, m.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			out[j] += m.data[i*m.cols+j]
		}
	}
	return out
}

// FrobeniusNorm returns the square root of the sum of squared entries.
func (m *Matrix) FrobeniusNorm() float64 {
	var s float64
	for _, x := range m.data {
		s += x * x
	}
	return math.Sqrt(s)
}

// MaxAbs returns the largest absolute entry.
func (m *Matrix) MaxAbs() float64 {
	var mx float64
	for _, x := range m.data {
		if a := math.Abs(x); a > mx {
			mx = a
		}
	}
	return mx
}

// SubMatrix returns the matrix obtained by deleting the given set of row indices
// and column indices. Indices out of range are ignored.
func (m *Matrix) SubMatrix(dropRows, dropCols map[int]bool) *Matrix {
	var rk, ck []int
	for i := 0; i < m.rows; i++ {
		if !dropRows[i] {
			rk = append(rk, i)
		}
	}
	for j := 0; j < m.cols; j++ {
		if !dropCols[j] {
			ck = append(ck, j)
		}
	}
	r := NewMatrix(len(rk), len(ck))
	for a, i := range rk {
		for b, j := range ck {
			r.data[a*r.cols+b] = m.data[i*m.cols+j]
		}
	}
	return r
}

// Apply returns a new matrix whose entries are f applied element-wise.
func (m *Matrix) Apply(f func(float64) float64) *Matrix {
	r := m.Clone()
	for i := range r.data {
		r.data[i] = f(r.data[i])
	}
	return r
}

// Pow returns m raised to the non-negative integer power p by repeated squaring.
// It returns ErrNotSquare for a non-square matrix and ErrInvalidArgument for a
// negative exponent.
func (m *Matrix) Pow(p int) (*Matrix, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	if p < 0 {
		return nil, ErrInvalidArgument
	}
	result := IdentityMatrix(m.rows)
	base := m.Clone()
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

// String renders the matrix with each entry formatted using %g.
func (m *Matrix) String() string {
	var b strings.Builder
	for i := 0; i < m.rows; i++ {
		b.WriteByte('[')
		for j := 0; j < m.cols; j++ {
			if j > 0 {
				b.WriteByte(' ')
			}
			fmt.Fprintf(&b, "%g", m.data[i*m.cols+j])
		}
		b.WriteByte(']')
		if i < m.rows-1 {
			b.WriteByte('\n')
		}
	}
	return b.String()
}

// ToRows returns the matrix as a fresh slice of row slices.
func (m *Matrix) ToRows() [][]float64 {
	out := make([][]float64, m.rows)
	for i := 0; i < m.rows; i++ {
		out[i] = m.Row(i)
	}
	return out
}

// OuterProduct returns the rank-one matrix a·bᵀ.
func OuterProduct(a, b []float64) *Matrix {
	r := NewMatrix(len(a), len(b))
	for i := range a {
		for j := range b {
			r.data[i*r.cols+j] = a[i] * b[j]
		}
	}
	return r
}
