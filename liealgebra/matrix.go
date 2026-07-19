package liealgebra

import (
	"math"
)

// Matrix is a dense real matrix stored in row-major order. The zero value is
// not usable; construct matrices with [NewMatrix] and the related helpers.
type Matrix struct {
	Rows int
	Cols int
	Data []float64 // len == Rows*Cols, element (i,j) at Data[i*Cols+j]
}

// NewMatrix returns a rows-by-cols zero matrix.
func NewMatrix(rows, cols int) *Matrix {
	if rows < 0 || cols < 0 {
		rows, cols = 0, 0
	}
	return &Matrix{Rows: rows, Cols: cols, Data: make([]float64, rows*cols)}
}

// NewMatrixFromRows builds a matrix from a slice of equal-length rows.
// It returns [ErrDim] if the rows are ragged.
func NewMatrixFromRows(rows [][]float64) (*Matrix, error) {
	r := len(rows)
	if r == 0 {
		return NewMatrix(0, 0), nil
	}
	c := len(rows[0])
	m := NewMatrix(r, c)
	for i := 0; i < r; i++ {
		if len(rows[i]) != c {
			return nil, ErrDim
		}
		copy(m.Data[i*c:(i+1)*c], rows[i])
	}
	return m, nil
}

// ZeroMatrix returns an n-by-n zero matrix.
func ZeroMatrix(n int) *Matrix { return NewMatrix(n, n) }

// IdentityMatrix returns the n-by-n identity matrix.
func IdentityMatrix(n int) *Matrix {
	m := NewMatrix(n, n)
	for i := 0; i < n; i++ {
		m.Data[i*n+i] = 1
	}
	return m
}

// DiagMatrix returns the square matrix with the given diagonal entries.
func DiagMatrix(d []float64) *Matrix {
	n := len(d)
	m := NewMatrix(n, n)
	for i := 0; i < n; i++ {
		m.Data[i*n+i] = d[i]
	}
	return m
}

// At returns element (i,j). It panics if the index is out of range.
func (m *Matrix) At(i, j int) float64 {
	if i < 0 || i >= m.Rows || j < 0 || j >= m.Cols {
		panic("liealgebra: Matrix.At index out of range")
	}
	return m.Data[i*m.Cols+j]
}

// Set assigns v to element (i,j). It panics if the index is out of range.
func (m *Matrix) Set(i, j int, v float64) {
	if i < 0 || i >= m.Rows || j < 0 || j >= m.Cols {
		panic("liealgebra: Matrix.Set index out of range")
	}
	m.Data[i*m.Cols+j] = v
}

// Dims returns the number of rows and columns.
func (m *Matrix) Dims() (int, int) { return m.Rows, m.Cols }

// IsSquare reports whether the matrix has equal row and column counts.
func (m *Matrix) IsSquare() bool { return m.Rows == m.Cols }

// Clone returns a deep copy of the matrix.
func (m *Matrix) Clone() *Matrix {
	c := NewMatrix(m.Rows, m.Cols)
	copy(c.Data, m.Data)
	return c
}

// Row returns a copy of row i.
func (m *Matrix) Row(i int) []float64 {
	out := make([]float64, m.Cols)
	copy(out, m.Data[i*m.Cols:(i+1)*m.Cols])
	return out
}

// Col returns a copy of column j.
func (m *Matrix) Col(j int) []float64 {
	out := make([]float64, m.Rows)
	for i := 0; i < m.Rows; i++ {
		out[i] = m.Data[i*m.Cols+j]
	}
	return out
}

// Equal reports whether a and b have the same shape and identical entries.
func (m *Matrix) Equal(b *Matrix) bool {
	if m.Rows != b.Rows || m.Cols != b.Cols {
		return false
	}
	for i := range m.Data {
		if m.Data[i] != b.Data[i] {
			return false
		}
	}
	return true
}

// ApproxEqual reports whether a and b have the same shape and every entry
// agrees to within absolute tolerance tol.
func (m *Matrix) ApproxEqual(b *Matrix, tol float64) bool {
	if m.Rows != b.Rows || m.Cols != b.Cols {
		return false
	}
	for i := range m.Data {
		if math.Abs(m.Data[i]-b.Data[i]) > tol {
			return false
		}
	}
	return true
}

// Add returns m+b or [ErrDim] on a shape mismatch.
func (m *Matrix) Add(b *Matrix) (*Matrix, error) {
	if m.Rows != b.Rows || m.Cols != b.Cols {
		return nil, ErrDim
	}
	out := NewMatrix(m.Rows, m.Cols)
	for i := range m.Data {
		out.Data[i] = m.Data[i] + b.Data[i]
	}
	return out, nil
}

// Sub returns m-b or [ErrDim] on a shape mismatch.
func (m *Matrix) Sub(b *Matrix) (*Matrix, error) {
	if m.Rows != b.Rows || m.Cols != b.Cols {
		return nil, ErrDim
	}
	out := NewMatrix(m.Rows, m.Cols)
	for i := range m.Data {
		out.Data[i] = m.Data[i] - b.Data[i]
	}
	return out, nil
}

// Scale returns the matrix multiplied by the scalar s.
func (m *Matrix) Scale(s float64) *Matrix {
	out := NewMatrix(m.Rows, m.Cols)
	for i := range m.Data {
		out.Data[i] = s * m.Data[i]
	}
	return out
}

// Neg returns the additive inverse of the matrix.
func (m *Matrix) Neg() *Matrix { return m.Scale(-1) }

// Mul returns the matrix product m*b or [ErrDim] if the inner dimensions
// disagree.
func (m *Matrix) Mul(b *Matrix) (*Matrix, error) {
	if m.Cols != b.Rows {
		return nil, ErrDim
	}
	out := NewMatrix(m.Rows, b.Cols)
	for i := 0; i < m.Rows; i++ {
		for k := 0; k < m.Cols; k++ {
			a := m.Data[i*m.Cols+k]
			if a == 0 {
				continue
			}
			for j := 0; j < b.Cols; j++ {
				out.Data[i*b.Cols+j] += a * b.Data[k*b.Cols+j]
			}
		}
	}
	return out, nil
}

// MatVec returns the matrix-vector product m*x or [ErrDim] on a mismatch.
func (m *Matrix) MatVec(x []float64) ([]float64, error) {
	if m.Cols != len(x) {
		return nil, ErrDim
	}
	out := make([]float64, m.Rows)
	for i := 0; i < m.Rows; i++ {
		s := 0.0
		for j := 0; j < m.Cols; j++ {
			s += m.Data[i*m.Cols+j] * x[j]
		}
		out[i] = s
	}
	return out, nil
}

// Transpose returns the transpose of the matrix.
func (m *Matrix) Transpose() *Matrix {
	out := NewMatrix(m.Cols, m.Rows)
	for i := 0; i < m.Rows; i++ {
		for j := 0; j < m.Cols; j++ {
			out.Data[j*m.Rows+i] = m.Data[i*m.Cols+j]
		}
	}
	return out
}

// Trace returns the sum of the diagonal entries. It returns [ErrNotSquare] for
// a non-square matrix.
func (m *Matrix) Trace() (float64, error) {
	if !m.IsSquare() {
		return 0, ErrNotSquare
	}
	s := 0.0
	for i := 0; i < m.Rows; i++ {
		s += m.Data[i*m.Cols+i]
	}
	return s, nil
}

// FrobeniusNorm returns the square root of the sum of squared entries.
func (m *Matrix) FrobeniusNorm() float64 {
	s := 0.0
	for _, v := range m.Data {
		s += v * v
	}
	return math.Sqrt(s)
}

// MaxAbs returns the largest absolute entry, or 0 for an empty matrix.
func (m *Matrix) MaxAbs() float64 {
	max := 0.0
	for _, v := range m.Data {
		if a := math.Abs(v); a > max {
			max = a
		}
	}
	return max
}

// IsSymmetric reports whether m equals its transpose within tolerance tol.
func (m *Matrix) IsSymmetric(tol float64) bool {
	if !m.IsSquare() {
		return false
	}
	n := m.Rows
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if math.Abs(m.Data[i*n+j]-m.Data[j*n+i]) > tol {
				return false
			}
		}
	}
	return true
}

// IsAntisymmetric reports whether m equals the negative of its transpose within
// tolerance tol (a real antisymmetric/skew matrix).
func (m *Matrix) IsAntisymmetric(tol float64) bool {
	if !m.IsSquare() {
		return false
	}
	n := m.Rows
	for i := 0; i < n; i++ {
		if math.Abs(m.Data[i*n+i]) > tol {
			return false
		}
		for j := i + 1; j < n; j++ {
			if math.Abs(m.Data[i*n+j]+m.Data[j*n+i]) > tol {
				return false
			}
		}
	}
	return true
}

// IsDiagonal reports whether all off-diagonal entries vanish within tol.
func (m *Matrix) IsDiagonal(tol float64) bool {
	if !m.IsSquare() {
		return false
	}
	n := m.Rows
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if i != j && math.Abs(m.Data[i*n+j]) > tol {
				return false
			}
		}
	}
	return true
}

// IsTraceless reports whether the trace is zero within tolerance tol.
func (m *Matrix) IsTraceless(tol float64) bool {
	t, err := m.Trace()
	if err != nil {
		return false
	}
	return math.Abs(t) <= tol
}

// Kronecker returns the Kronecker (tensor) product m ⊗ b.
func (m *Matrix) Kronecker(b *Matrix) *Matrix {
	out := NewMatrix(m.Rows*b.Rows, m.Cols*b.Cols)
	oc := out.Cols
	for i := 0; i < m.Rows; i++ {
		for j := 0; j < m.Cols; j++ {
			a := m.Data[i*m.Cols+j]
			for p := 0; p < b.Rows; p++ {
				for q := 0; q < b.Cols; q++ {
					ri := i*b.Rows + p
					cj := j*b.Cols + q
					out.Data[ri*oc+cj] = a * b.Data[p*b.Cols+q]
				}
			}
		}
	}
	return out
}

// Symmetrize returns (m+mᵀ)/2, the symmetric part of a square matrix.
func (m *Matrix) Symmetrize() (*Matrix, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	t := m.Transpose()
	s, _ := m.Add(t)
	return s.Scale(0.5), nil
}

// AntiSymmetrize returns (m-mᵀ)/2, the antisymmetric part of a square matrix.
func (m *Matrix) AntiSymmetrize() (*Matrix, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	t := m.Transpose()
	s, _ := m.Sub(t)
	return s.Scale(0.5), nil
}

// Vec flattens the matrix into a single row-major slice (a fresh copy).
func (m *Matrix) Vec() []float64 {
	out := make([]float64, len(m.Data))
	copy(out, m.Data)
	return out
}
