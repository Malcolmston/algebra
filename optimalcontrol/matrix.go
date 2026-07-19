package optimalcontrol

import (
	"errors"
	"fmt"
	"math"
)

// Matrix is a dense real matrix stored in row-major order. It is the basic
// linear-algebra container used throughout the optimalcontrol package. The zero
// value is not usable; construct matrices with NewMatrix, Zeros and friends.
type Matrix struct {
	rows, cols int
	data       []float64
}

// ErrDim is returned when matrix or vector dimensions are incompatible.
var ErrDim = errors.New("optimalcontrol: incompatible dimensions")

// ErrSingular is returned when a matrix that must be invertible is singular.
var ErrSingular = errors.New("optimalcontrol: matrix is singular")

// ErrNotConverged is returned by iterative solvers that fail to converge within
// the allotted number of iterations.
var ErrNotConverged = errors.New("optimalcontrol: iteration did not converge")

// NewMatrix builds an r×c matrix from the supplied row-major data. The data
// slice is copied. It panics if len(data) != r*c.
func NewMatrix(r, c int, data []float64) *Matrix {
	if r < 0 || c < 0 {
		panic("optimalcontrol: negative dimension")
	}
	if len(data) != r*c {
		panic("optimalcontrol: data length does not match dimensions")
	}
	cp := make([]float64, len(data))
	copy(cp, data)
	return &Matrix{rows: r, cols: c, data: cp}
}

// Zeros returns an r×c matrix of zeros.
func Zeros(r, c int) *Matrix {
	return &Matrix{rows: r, cols: c, data: make([]float64, r*c)}
}

// Ones returns an r×c matrix whose entries are all one.
func Ones(r, c int) *Matrix {
	m := Zeros(r, c)
	for i := range m.data {
		m.data[i] = 1
	}
	return m
}

// Eye returns the n×n identity matrix.
func Eye(n int) *Matrix {
	m := Zeros(n, n)
	for i := 0; i < n; i++ {
		m.data[i*n+i] = 1
	}
	return m
}

// Diag returns a square diagonal matrix whose diagonal is v.
func Diag(v []float64) *Matrix {
	n := len(v)
	m := Zeros(n, n)
	for i := 0; i < n; i++ {
		m.data[i*n+i] = v[i]
	}
	return m
}

// FromRows builds a matrix from a slice of equal-length rows.
func FromRows(rows [][]float64) *Matrix {
	r := len(rows)
	if r == 0 {
		return Zeros(0, 0)
	}
	c := len(rows[0])
	m := Zeros(r, c)
	for i := 0; i < r; i++ {
		if len(rows[i]) != c {
			panic("optimalcontrol: ragged rows")
		}
		copy(m.data[i*c:(i+1)*c], rows[i])
	}
	return m
}

// ColumnVector builds an n×1 matrix from v.
func ColumnVector(v []float64) *Matrix {
	return NewMatrix(len(v), 1, v)
}

// RowVector builds a 1×n matrix from v.
func RowVector(v []float64) *Matrix {
	return NewMatrix(1, len(v), v)
}

// Rows returns the number of rows.
func (m *Matrix) Rows() int { return m.rows }

// Cols returns the number of columns.
func (m *Matrix) Cols() int { return m.cols }

// IsSquare reports whether the matrix is square.
func (m *Matrix) IsSquare() bool { return m.rows == m.cols }

// At returns the element at row i, column j.
func (m *Matrix) At(i, j int) float64 { return m.data[i*m.cols+j] }

// Set assigns v to the element at row i, column j.
func (m *Matrix) Set(i, j int, v float64) { m.data[i*m.cols+j] = v }

// Add accumulates v into the element at row i, column j.
func (m *Matrix) Add(i, j int, v float64) { m.data[i*m.cols+j] += v }

// Data returns a copy of the underlying row-major data.
func (m *Matrix) Data() []float64 {
	cp := make([]float64, len(m.data))
	copy(cp, m.data)
	return cp
}

// Clone returns a deep copy of the matrix.
func (m *Matrix) Clone() *Matrix {
	return NewMatrix(m.rows, m.cols, m.data)
}

// Row returns a copy of row i as a slice.
func (m *Matrix) Row(i int) []float64 {
	out := make([]float64, m.cols)
	copy(out, m.data[i*m.cols:(i+1)*m.cols])
	return out
}

// Col returns a copy of column j as a slice.
func (m *Matrix) Col(j int) []float64 {
	out := make([]float64, m.rows)
	for i := 0; i < m.rows; i++ {
		out[i] = m.data[i*m.cols+j]
	}
	return out
}

// SetRow overwrites row i with v.
func (m *Matrix) SetRow(i int, v []float64) {
	if len(v) != m.cols {
		panic("optimalcontrol: row length mismatch")
	}
	copy(m.data[i*m.cols:(i+1)*m.cols], v)
}

// SetCol overwrites column j with v.
func (m *Matrix) SetCol(j int, v []float64) {
	if len(v) != m.rows {
		panic("optimalcontrol: column length mismatch")
	}
	for i := 0; i < m.rows; i++ {
		m.data[i*m.cols+j] = v[i]
	}
}

// Equal reports whether a and b have the same shape and identical entries.
func (m *Matrix) Equal(b *Matrix) bool {
	if m.rows != b.rows || m.cols != b.cols {
		return false
	}
	for i := range m.data {
		if m.data[i] != b.data[i] {
			return false
		}
	}
	return true
}

// ApproxEqual reports whether a and b have the same shape and entries within
// absolute tolerance tol.
func (m *Matrix) ApproxEqual(b *Matrix, tol float64) bool {
	if m.rows != b.rows || m.cols != b.cols {
		return false
	}
	for i := range m.data {
		if math.Abs(m.data[i]-b.data[i]) > tol {
			return false
		}
	}
	return true
}

// Plus returns the sum a+b.
func (m *Matrix) Plus(b *Matrix) *Matrix {
	if m.rows != b.rows || m.cols != b.cols {
		panic(ErrDim)
	}
	out := Zeros(m.rows, m.cols)
	for i := range m.data {
		out.data[i] = m.data[i] + b.data[i]
	}
	return out
}

// Minus returns the difference a-b.
func (m *Matrix) Minus(b *Matrix) *Matrix {
	if m.rows != b.rows || m.cols != b.cols {
		panic(ErrDim)
	}
	out := Zeros(m.rows, m.cols)
	for i := range m.data {
		out.data[i] = m.data[i] - b.data[i]
	}
	return out
}

// Scale returns the matrix scaled by s.
func (m *Matrix) Scale(s float64) *Matrix {
	out := Zeros(m.rows, m.cols)
	for i := range m.data {
		out.data[i] = m.data[i] * s
	}
	return out
}

// Neg returns the additive inverse of the matrix.
func (m *Matrix) Neg() *Matrix { return m.Scale(-1) }

// Mul returns the matrix product a·b.
func (m *Matrix) Mul(b *Matrix) *Matrix {
	if m.cols != b.rows {
		panic(ErrDim)
	}
	out := Zeros(m.rows, b.cols)
	for i := 0; i < m.rows; i++ {
		for k := 0; k < m.cols; k++ {
			aik := m.data[i*m.cols+k]
			if aik == 0 {
				continue
			}
			for j := 0; j < b.cols; j++ {
				out.data[i*b.cols+j] += aik * b.data[k*b.cols+j]
			}
		}
	}
	return out
}

// MulVec returns the matrix-vector product a·v.
func (m *Matrix) MulVec(v []float64) []float64 {
	if m.cols != len(v) {
		panic(ErrDim)
	}
	out := make([]float64, m.rows)
	for i := 0; i < m.rows; i++ {
		var s float64
		for j := 0; j < m.cols; j++ {
			s += m.data[i*m.cols+j] * v[j]
		}
		out[i] = s
	}
	return out
}

// Transpose returns the transpose of the matrix.
func (m *Matrix) Transpose() *Matrix {
	out := Zeros(m.cols, m.rows)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			out.data[j*m.rows+i] = m.data[i*m.cols+j]
		}
	}
	return out
}

// Trace returns the sum of the diagonal entries of a square matrix.
func (m *Matrix) Trace() float64 {
	if !m.IsSquare() {
		panic(ErrDim)
	}
	var s float64
	for i := 0; i < m.rows; i++ {
		s += m.data[i*m.cols+i]
	}
	return s
}

// Symmetrize returns (A + Aᵀ)/2, the symmetric part of a square matrix.
func (m *Matrix) Symmetrize() *Matrix {
	if !m.IsSquare() {
		panic(ErrDim)
	}
	out := Zeros(m.rows, m.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			out.data[i*m.cols+j] = 0.5 * (m.data[i*m.cols+j] + m.data[j*m.cols+i])
		}
	}
	return out
}

// IsSymmetric reports whether the matrix is square and symmetric within tol.
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

// FrobeniusNorm returns the Frobenius norm of the matrix.
func (m *Matrix) FrobeniusNorm() float64 {
	var s float64
	for _, x := range m.data {
		s += x * x
	}
	return math.Sqrt(s)
}

// MaxAbs returns the largest absolute entry of the matrix.
func (m *Matrix) MaxAbs() float64 {
	var mx float64
	for _, x := range m.data {
		if a := math.Abs(x); a > mx {
			mx = a
		}
	}
	return mx
}

// OneNorm returns the maximum absolute column sum of the matrix.
func (m *Matrix) OneNorm() float64 {
	var mx float64
	for j := 0; j < m.cols; j++ {
		var s float64
		for i := 0; i < m.rows; i++ {
			s += math.Abs(m.data[i*m.cols+j])
		}
		if s > mx {
			mx = s
		}
	}
	return mx
}

// InfNorm returns the maximum absolute row sum of the matrix.
func (m *Matrix) InfNorm() float64 {
	var mx float64
	for i := 0; i < m.rows; i++ {
		var s float64
		for j := 0; j < m.cols; j++ {
			s += math.Abs(m.data[i*m.cols+j])
		}
		if s > mx {
			mx = s
		}
	}
	return mx
}

// Submatrix returns the r0..r1-1 × c0..c1-1 block of the matrix.
func (m *Matrix) Submatrix(r0, r1, c0, c1 int) *Matrix {
	out := Zeros(r1-r0, c1-c0)
	for i := r0; i < r1; i++ {
		for j := c0; j < c1; j++ {
			out.data[(i-r0)*out.cols+(j-c0)] = m.data[i*m.cols+j]
		}
	}
	return out
}

// SetBlock copies b into the block of m with top-left corner at (r0, c0).
func (m *Matrix) SetBlock(r0, c0 int, b *Matrix) {
	for i := 0; i < b.rows; i++ {
		for j := 0; j < b.cols; j++ {
			m.data[(r0+i)*m.cols+(c0+j)] = b.data[i*b.cols+j]
		}
	}
}

// HStack returns the horizontal concatenation [a b].
func HStack(a, b *Matrix) *Matrix {
	if a.rows != b.rows {
		panic(ErrDim)
	}
	out := Zeros(a.rows, a.cols+b.cols)
	out.SetBlock(0, 0, a)
	out.SetBlock(0, a.cols, b)
	return out
}

// VStack returns the vertical concatenation [a; b].
func VStack(a, b *Matrix) *Matrix {
	if a.cols != b.cols {
		panic(ErrDim)
	}
	out := Zeros(a.rows+b.rows, a.cols)
	out.SetBlock(0, 0, a)
	out.SetBlock(a.rows, 0, b)
	return out
}

// BlockMatrix assembles the 2×2 block matrix [[a b];[c d]].
func BlockMatrix(a, b, c, d *Matrix) *Matrix {
	return VStack(HStack(a, b), HStack(c, d))
}

// Kron returns the Kronecker product a⊗b.
func Kron(a, b *Matrix) *Matrix {
	out := Zeros(a.rows*b.rows, a.cols*b.cols)
	for i := 0; i < a.rows; i++ {
		for j := 0; j < a.cols; j++ {
			aij := a.data[i*a.cols+j]
			for p := 0; p < b.rows; p++ {
				for q := 0; q < b.cols; q++ {
					out.data[(i*b.rows+p)*out.cols+(j*b.cols+q)] = aij * b.data[p*b.cols+q]
				}
			}
		}
	}
	return out
}

// Vec returns the column-major vectorization of the matrix (columns stacked).
func (m *Matrix) Vec() []float64 {
	out := make([]float64, m.rows*m.cols)
	k := 0
	for j := 0; j < m.cols; j++ {
		for i := 0; i < m.rows; i++ {
			out[k] = m.data[i*m.cols+j]
			k++
		}
	}
	return out
}

// Unvec reshapes a column-major vector into an r×c matrix.
func Unvec(v []float64, r, c int) *Matrix {
	if len(v) != r*c {
		panic(ErrDim)
	}
	m := Zeros(r, c)
	k := 0
	for j := 0; j < c; j++ {
		for i := 0; i < r; i++ {
			m.data[i*c+j] = v[k]
			k++
		}
	}
	return m
}

// String renders the matrix for debugging.
func (m *Matrix) String() string {
	s := fmt.Sprintf("Matrix %dx%d:\n", m.rows, m.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			s += fmt.Sprintf("% .6g ", m.data[i*m.cols+j])
		}
		s += "\n"
	}
	return s
}
