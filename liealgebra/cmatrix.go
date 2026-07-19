package liealgebra

import (
	"math"
	"math/cmplx"
)

// CMatrix is a dense complex matrix stored in row-major order. Complex matrices
// are needed for the unitary Lie algebras su(n), whose generators are naturally
// complex. The zero value is not usable; construct with [NewCMatrix].
type CMatrix struct {
	Rows int
	Cols int
	Data []complex128
}

// NewCMatrix returns a rows-by-cols zero complex matrix.
func NewCMatrix(rows, cols int) *CMatrix {
	if rows < 0 || cols < 0 {
		rows, cols = 0, 0
	}
	return &CMatrix{Rows: rows, Cols: cols, Data: make([]complex128, rows*cols)}
}

// NewCMatrixFromRows builds a complex matrix from equal-length rows, returning
// [ErrDim] if the rows are ragged.
func NewCMatrixFromRows(rows [][]complex128) (*CMatrix, error) {
	r := len(rows)
	if r == 0 {
		return NewCMatrix(0, 0), nil
	}
	c := len(rows[0])
	m := NewCMatrix(r, c)
	for i := 0; i < r; i++ {
		if len(rows[i]) != c {
			return nil, ErrDim
		}
		copy(m.Data[i*c:(i+1)*c], rows[i])
	}
	return m, nil
}

// IdentityCMatrix returns the n-by-n complex identity matrix.
func IdentityCMatrix(n int) *CMatrix {
	m := NewCMatrix(n, n)
	for i := 0; i < n; i++ {
		m.Data[i*n+i] = 1
	}
	return m
}

// DiagCMatrix returns the square complex matrix with the given diagonal.
func DiagCMatrix(d []complex128) *CMatrix {
	n := len(d)
	m := NewCMatrix(n, n)
	for i := 0; i < n; i++ {
		m.Data[i*n+i] = d[i]
	}
	return m
}

// RealToComplex promotes a real matrix to a complex matrix.
func RealToComplex(m *Matrix) *CMatrix {
	c := NewCMatrix(m.Rows, m.Cols)
	for i, v := range m.Data {
		c.Data[i] = complex(v, 0)
	}
	return c
}

// At returns element (i,j), panicking on an out-of-range index.
func (m *CMatrix) At(i, j int) complex128 {
	if i < 0 || i >= m.Rows || j < 0 || j >= m.Cols {
		panic("liealgebra: CMatrix.At index out of range")
	}
	return m.Data[i*m.Cols+j]
}

// Set assigns v to element (i,j), panicking on an out-of-range index.
func (m *CMatrix) Set(i, j int, v complex128) {
	if i < 0 || i >= m.Rows || j < 0 || j >= m.Cols {
		panic("liealgebra: CMatrix.Set index out of range")
	}
	m.Data[i*m.Cols+j] = v
}

// Dims returns the number of rows and columns.
func (m *CMatrix) Dims() (int, int) { return m.Rows, m.Cols }

// IsSquare reports whether the matrix is square.
func (m *CMatrix) IsSquare() bool { return m.Rows == m.Cols }

// Clone returns a deep copy.
func (m *CMatrix) Clone() *CMatrix {
	c := NewCMatrix(m.Rows, m.Cols)
	copy(c.Data, m.Data)
	return c
}

// Equal reports exact equality of shape and entries.
func (m *CMatrix) Equal(b *CMatrix) bool {
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

// ApproxEqual reports shape equality with every entry agreeing to within tol in
// complex modulus.
func (m *CMatrix) ApproxEqual(b *CMatrix, tol float64) bool {
	if m.Rows != b.Rows || m.Cols != b.Cols {
		return false
	}
	for i := range m.Data {
		if cmplx.Abs(m.Data[i]-b.Data[i]) > tol {
			return false
		}
	}
	return true
}

// Add returns m+b or [ErrDim].
func (m *CMatrix) Add(b *CMatrix) (*CMatrix, error) {
	if m.Rows != b.Rows || m.Cols != b.Cols {
		return nil, ErrDim
	}
	out := NewCMatrix(m.Rows, m.Cols)
	for i := range m.Data {
		out.Data[i] = m.Data[i] + b.Data[i]
	}
	return out, nil
}

// Sub returns m-b or [ErrDim].
func (m *CMatrix) Sub(b *CMatrix) (*CMatrix, error) {
	if m.Rows != b.Rows || m.Cols != b.Cols {
		return nil, ErrDim
	}
	out := NewCMatrix(m.Rows, m.Cols)
	for i := range m.Data {
		out.Data[i] = m.Data[i] - b.Data[i]
	}
	return out, nil
}

// Scale returns the matrix times the complex scalar s.
func (m *CMatrix) Scale(s complex128) *CMatrix {
	out := NewCMatrix(m.Rows, m.Cols)
	for i := range m.Data {
		out.Data[i] = s * m.Data[i]
	}
	return out
}

// Neg returns the additive inverse.
func (m *CMatrix) Neg() *CMatrix { return m.Scale(-1) }

// Mul returns the matrix product m*b or [ErrDim].
func (m *CMatrix) Mul(b *CMatrix) (*CMatrix, error) {
	if m.Cols != b.Rows {
		return nil, ErrDim
	}
	out := NewCMatrix(m.Rows, b.Cols)
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

// Transpose returns the (unconjugated) transpose.
func (m *CMatrix) Transpose() *CMatrix {
	out := NewCMatrix(m.Cols, m.Rows)
	for i := 0; i < m.Rows; i++ {
		for j := 0; j < m.Cols; j++ {
			out.Data[j*m.Rows+i] = m.Data[i*m.Cols+j]
		}
	}
	return out
}

// Conjugate returns the entrywise complex conjugate.
func (m *CMatrix) Conjugate() *CMatrix {
	out := NewCMatrix(m.Rows, m.Cols)
	for i, v := range m.Data {
		out.Data[i] = cmplx.Conj(v)
	}
	return out
}

// ConjugateTranspose returns the Hermitian conjugate (dagger) mᴴ.
func (m *CMatrix) ConjugateTranspose() *CMatrix {
	out := NewCMatrix(m.Cols, m.Rows)
	for i := 0; i < m.Rows; i++ {
		for j := 0; j < m.Cols; j++ {
			out.Data[j*m.Rows+i] = cmplx.Conj(m.Data[i*m.Cols+j])
		}
	}
	return out
}

// Dagger is an alias for [CMatrix.ConjugateTranspose].
func (m *CMatrix) Dagger() *CMatrix { return m.ConjugateTranspose() }

// Trace returns the sum of diagonal entries, or [ErrNotSquare].
func (m *CMatrix) Trace() (complex128, error) {
	if !m.IsSquare() {
		return 0, ErrNotSquare
	}
	var s complex128
	for i := 0; i < m.Rows; i++ {
		s += m.Data[i*m.Cols+i]
	}
	return s, nil
}

// Real returns the real part of each entry as a real matrix.
func (m *CMatrix) Real() *Matrix {
	out := NewMatrix(m.Rows, m.Cols)
	for i, v := range m.Data {
		out.Data[i] = real(v)
	}
	return out
}

// Imag returns the imaginary part of each entry as a real matrix.
func (m *CMatrix) Imag() *Matrix {
	out := NewMatrix(m.Rows, m.Cols)
	for i, v := range m.Data {
		out.Data[i] = imag(v)
	}
	return out
}

// FrobeniusNorm returns the square root of the sum of squared moduli.
func (m *CMatrix) FrobeniusNorm() float64 {
	s := 0.0
	for _, v := range m.Data {
		a := cmplx.Abs(v)
		s += a * a
	}
	return math.Sqrt(s)
}

// MaxAbs returns the largest entry modulus.
func (m *CMatrix) MaxAbs() float64 {
	max := 0.0
	for _, v := range m.Data {
		if a := cmplx.Abs(v); a > max {
			max = a
		}
	}
	return max
}

// IsHermitian reports whether m equals its conjugate transpose within tol.
func (m *CMatrix) IsHermitian(tol float64) bool {
	if !m.IsSquare() {
		return false
	}
	n := m.Rows
	for i := 0; i < n; i++ {
		for j := i; j < n; j++ {
			if cmplx.Abs(m.Data[i*n+j]-cmplx.Conj(m.Data[j*n+i])) > tol {
				return false
			}
		}
	}
	return true
}

// IsAntiHermitian reports whether m equals the negative of its conjugate
// transpose within tol (a skew-Hermitian matrix, as in su(n)).
func (m *CMatrix) IsAntiHermitian(tol float64) bool {
	if !m.IsSquare() {
		return false
	}
	n := m.Rows
	for i := 0; i < n; i++ {
		for j := i; j < n; j++ {
			if cmplx.Abs(m.Data[i*n+j]+cmplx.Conj(m.Data[j*n+i])) > tol {
				return false
			}
		}
	}
	return true
}

// IsUnitary reports whether mᴴm is the identity within tol.
func (m *CMatrix) IsUnitary(tol float64) bool {
	if !m.IsSquare() {
		return false
	}
	p, err := m.Dagger().Mul(m)
	if err != nil {
		return false
	}
	n := m.Rows
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			want := complex128(0)
			if i == j {
				want = 1
			}
			if cmplx.Abs(p.Data[i*n+j]-want) > tol {
				return false
			}
		}
	}
	return true
}

// IsTraceless reports whether the trace is zero within tol.
func (m *CMatrix) IsTraceless(tol float64) bool {
	t, err := m.Trace()
	if err != nil {
		return false
	}
	return cmplx.Abs(t) <= tol
}
