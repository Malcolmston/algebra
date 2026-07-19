package randommatrix

import (
	"errors"
	"fmt"
	"math"
	"strings"
)

// ErrDimensionMismatch is returned when two matrices have incompatible shapes
// for the requested operation.
var ErrDimensionMismatch = errors.New("randommatrix: dimension mismatch")

// ErrNotSquare is returned when a square matrix is required but a rectangular
// one was supplied.
var ErrNotSquare = errors.New("randommatrix: matrix is not square")

// ErrNotSymmetric is returned when a symmetric matrix is required.
var ErrNotSymmetric = errors.New("randommatrix: matrix is not symmetric")

// Matrix is a dense, row-major matrix of float64 values.
type Matrix struct {
	rows, cols int
	data       []float64
}

// NewMatrix returns a rows-by-cols matrix whose entries are all zero. It panics
// if either dimension is negative.
func NewMatrix(rows, cols int) *Matrix {
	if rows < 0 || cols < 0 {
		panic("randommatrix: negative matrix dimension")
	}
	return &Matrix{rows: rows, cols: cols, data: make([]float64, rows*cols)}
}

// NewMatrixFrom builds a rows-by-cols matrix from a flat, row-major slice. The
// slice is copied. It panics if len(data) != rows*cols.
func NewMatrixFrom(rows, cols int, data []float64) *Matrix {
	if len(data) != rows*cols {
		panic("randommatrix: data length does not match dimensions")
	}
	cp := make([]float64, len(data))
	copy(cp, data)
	return &Matrix{rows: rows, cols: cols, data: cp}
}

// NewMatrixFromRows builds a matrix from a slice of equal-length rows.
func NewMatrixFromRows(rows [][]float64) *Matrix {
	if len(rows) == 0 {
		return NewMatrix(0, 0)
	}
	c := len(rows[0])
	m := NewMatrix(len(rows), c)
	for i, r := range rows {
		if len(r) != c {
			panic("randommatrix: ragged rows")
		}
		copy(m.data[i*c:(i+1)*c], r)
	}
	return m
}

// Zeros returns a rows-by-cols zero matrix; it is an alias for NewMatrix.
func Zeros(rows, cols int) *Matrix { return NewMatrix(rows, cols) }

// Identity returns the n-by-n identity matrix.
func Identity(n int) *Matrix {
	m := NewMatrix(n, n)
	for i := 0; i < n; i++ {
		m.data[i*n+i] = 1
	}
	return m
}

// Diag returns a square diagonal matrix whose diagonal is d.
func Diag(d []float64) *Matrix {
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

// IsSquare reports whether the matrix is square.
func (m *Matrix) IsSquare() bool { return m.rows == m.cols }

// At returns the entry at row i, column j.
func (m *Matrix) At(i, j int) float64 { return m.data[i*m.cols+j] }

// Set assigns v to the entry at row i, column j.
func (m *Matrix) Set(i, j int, v float64) { m.data[i*m.cols+j] = v }

// Add increments the entry at row i, column j by v.
func (m *Matrix) Add(i, j int, v float64) { m.data[i*m.cols+j] += v }

// Data returns a copy of the underlying row-major slice.
func (m *Matrix) Data() []float64 {
	cp := make([]float64, len(m.data))
	copy(cp, m.data)
	return cp
}

// Row returns a copy of row i.
func (m *Matrix) Row(i int) []float64 {
	r := make([]float64, m.cols)
	copy(r, m.data[i*m.cols:(i+1)*m.cols])
	return r
}

// Col returns a copy of column j.
func (m *Matrix) Col(j int) []float64 {
	c := make([]float64, m.rows)
	for i := 0; i < m.rows; i++ {
		c[i] = m.data[i*m.cols+j]
	}
	return c
}

// Diagonal returns a copy of the main diagonal.
func (m *Matrix) Diagonal() []float64 {
	n := m.rows
	if m.cols < n {
		n = m.cols
	}
	d := make([]float64, n)
	for i := 0; i < n; i++ {
		d[i] = m.data[i*m.cols+i]
	}
	return d
}

// Clone returns an independent deep copy of the matrix.
func (m *Matrix) Clone() *Matrix {
	return &Matrix{rows: m.rows, cols: m.cols, data: m.Data()}
}

// Equals reports whether every entry of m and other agree to within tol.
func (m *Matrix) Equals(other *Matrix, tol float64) bool {
	if m.rows != other.rows || m.cols != other.cols {
		return false
	}
	for i := range m.data {
		if math.Abs(m.data[i]-other.data[i]) > tol {
			return false
		}
	}
	return true
}

// Transpose returns the transpose of the matrix.
func (m *Matrix) Transpose() *Matrix {
	t := NewMatrix(m.cols, m.rows)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			t.data[j*m.rows+i] = m.data[i*m.cols+j]
		}
	}
	return t
}

// Plus returns the entrywise sum m + other. It returns an error on shape
// mismatch.
func (m *Matrix) Plus(other *Matrix) (*Matrix, error) {
	if m.rows != other.rows || m.cols != other.cols {
		return nil, ErrDimensionMismatch
	}
	r := m.Clone()
	for i := range r.data {
		r.data[i] += other.data[i]
	}
	return r, nil
}

// Minus returns the entrywise difference m - other.
func (m *Matrix) Minus(other *Matrix) (*Matrix, error) {
	if m.rows != other.rows || m.cols != other.cols {
		return nil, ErrDimensionMismatch
	}
	r := m.Clone()
	for i := range r.data {
		r.data[i] -= other.data[i]
	}
	return r, nil
}

// Scale returns a new matrix equal to m with every entry multiplied by s.
func (m *Matrix) Scale(s float64) *Matrix {
	r := m.Clone()
	for i := range r.data {
		r.data[i] *= s
	}
	return r
}

// Mul returns the matrix product m * other.
func (m *Matrix) Mul(other *Matrix) (*Matrix, error) {
	if m.cols != other.rows {
		return nil, ErrDimensionMismatch
	}
	r := NewMatrix(m.rows, other.cols)
	for i := 0; i < m.rows; i++ {
		for k := 0; k < m.cols; k++ {
			a := m.data[i*m.cols+k]
			if a == 0 {
				continue
			}
			for j := 0; j < other.cols; j++ {
				r.data[i*other.cols+j] += a * other.data[k*other.cols+j]
			}
		}
	}
	return r, nil
}

// MulVec returns the matrix-vector product m * v.
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

// Trace returns the sum of the diagonal entries. It panics for a non-square
// matrix.
func (m *Matrix) Trace() float64 {
	if !m.IsSquare() {
		panic(ErrNotSquare)
	}
	var s float64
	for i := 0; i < m.rows; i++ {
		s += m.data[i*m.cols+i]
	}
	return s
}

// FrobeniusNorm returns the Frobenius norm, sqrt of the sum of squared entries.
func (m *Matrix) FrobeniusNorm() float64 {
	var s float64
	for _, v := range m.data {
		s += v * v
	}
	return math.Sqrt(s)
}

// MaxAbs returns the largest absolute value of any entry.
func (m *Matrix) MaxAbs() float64 {
	var mx float64
	for _, v := range m.data {
		if a := math.Abs(v); a > mx {
			mx = a
		}
	}
	return mx
}

// IsSymmetric reports whether m is square and symmetric to within tol.
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
	if !m.IsSquare() {
		panic(ErrNotSquare)
	}
	n := m.rows
	r := NewMatrix(n, n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			r.data[i*n+j] = 0.5 * (m.data[i*n+j] + m.data[j*n+i])
		}
	}
	return r
}

// SubMatrix returns the r0:r1, c0:c1 block (half-open ranges).
func (m *Matrix) SubMatrix(r0, r1, c0, c1 int) *Matrix {
	out := NewMatrix(r1-r0, c1-c0)
	for i := r0; i < r1; i++ {
		for j := c0; j < c1; j++ {
			out.data[(i-r0)*out.cols+(j-c0)] = m.data[i*m.cols+j]
		}
	}
	return out
}

// String renders the matrix with entries formatted to three decimals.
func (m *Matrix) String() string {
	var b strings.Builder
	for i := 0; i < m.rows; i++ {
		b.WriteByte('[')
		for j := 0; j < m.cols; j++ {
			if j > 0 {
				b.WriteByte(' ')
			}
			fmt.Fprintf(&b, "%8.3f", m.data[i*m.cols+j])
		}
		b.WriteString("]\n")
	}
	return b.String()
}
