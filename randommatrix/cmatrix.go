package randommatrix

import (
	"math"
	"math/cmplx"
)

// CMatrix is a dense, row-major matrix of complex128 values.
type CMatrix struct {
	rows, cols int
	data       []complex128
}

// NewCMatrix returns a rows-by-cols complex matrix of zeros.
func NewCMatrix(rows, cols int) *CMatrix {
	if rows < 0 || cols < 0 {
		panic("randommatrix: negative matrix dimension")
	}
	return &CMatrix{rows: rows, cols: cols, data: make([]complex128, rows*cols)}
}

// NewCMatrixFrom builds a complex matrix from a flat row-major slice, which is
// copied.
func NewCMatrixFrom(rows, cols int, data []complex128) *CMatrix {
	if len(data) != rows*cols {
		panic("randommatrix: data length does not match dimensions")
	}
	cp := make([]complex128, len(data))
	copy(cp, data)
	return &CMatrix{rows: rows, cols: cols, data: cp}
}

// CIdentity returns the n-by-n complex identity matrix.
func CIdentity(n int) *CMatrix {
	m := NewCMatrix(n, n)
	for i := 0; i < n; i++ {
		m.data[i*n+i] = 1
	}
	return m
}

// Rows returns the number of rows.
func (m *CMatrix) Rows() int { return m.rows }

// Cols returns the number of columns.
func (m *CMatrix) Cols() int { return m.cols }

// IsSquare reports whether the matrix is square.
func (m *CMatrix) IsSquare() bool { return m.rows == m.cols }

// At returns the entry at row i, column j.
func (m *CMatrix) At(i, j int) complex128 { return m.data[i*m.cols+j] }

// Set assigns v to the entry at row i, column j.
func (m *CMatrix) Set(i, j int, v complex128) { m.data[i*m.cols+j] = v }

// Clone returns an independent deep copy.
func (m *CMatrix) Clone() *CMatrix {
	cp := make([]complex128, len(m.data))
	copy(cp, m.data)
	return &CMatrix{rows: m.rows, cols: m.cols, data: cp}
}

// ConjugateTranspose returns the conjugate (Hermitian) transpose m†.
func (m *CMatrix) ConjugateTranspose() *CMatrix {
	t := NewCMatrix(m.cols, m.rows)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			t.data[j*m.rows+i] = cmplx.Conj(m.data[i*m.cols+j])
		}
	}
	return t
}

// Plus returns the entrywise sum m + other.
func (m *CMatrix) Plus(other *CMatrix) (*CMatrix, error) {
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
func (m *CMatrix) Minus(other *CMatrix) (*CMatrix, error) {
	if m.rows != other.rows || m.cols != other.cols {
		return nil, ErrDimensionMismatch
	}
	r := m.Clone()
	for i := range r.data {
		r.data[i] -= other.data[i]
	}
	return r, nil
}

// Scale multiplies every entry by the complex scalar s.
func (m *CMatrix) Scale(s complex128) *CMatrix {
	r := m.Clone()
	for i := range r.data {
		r.data[i] *= s
	}
	return r
}

// Mul returns the matrix product m * other.
func (m *CMatrix) Mul(other *CMatrix) (*CMatrix, error) {
	if m.cols != other.rows {
		return nil, ErrDimensionMismatch
	}
	r := NewCMatrix(m.rows, other.cols)
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

// Trace returns the sum of the diagonal entries.
func (m *CMatrix) Trace() complex128 {
	if !m.IsSquare() {
		panic(ErrNotSquare)
	}
	var s complex128
	for i := 0; i < m.rows; i++ {
		s += m.data[i*m.cols+i]
	}
	return s
}

// FrobeniusNorm returns the Frobenius norm of the complex matrix.
func (m *CMatrix) FrobeniusNorm() float64 {
	var s float64
	for _, v := range m.data {
		s += real(v)*real(v) + imag(v)*imag(v)
	}
	return math.Sqrt(s)
}

// IsHermitian reports whether m is square and Hermitian to within tol.
func (m *CMatrix) IsHermitian(tol float64) bool {
	if !m.IsSquare() {
		return false
	}
	for i := 0; i < m.rows; i++ {
		if math.Abs(imag(m.data[i*m.cols+i])) > tol {
			return false
		}
		for j := i + 1; j < m.cols; j++ {
			d := m.data[i*m.cols+j] - cmplx.Conj(m.data[j*m.cols+i])
			if cmplx.Abs(d) > tol {
				return false
			}
		}
	}
	return true
}

// Hermitianize returns (m + m†)/2, the Hermitian part of a square matrix.
func (m *CMatrix) Hermitianize() *CMatrix {
	if !m.IsSquare() {
		panic(ErrNotSquare)
	}
	n := m.rows
	r := NewCMatrix(n, n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			r.data[i*n+j] = 0.5 * (m.data[i*n+j] + cmplx.Conj(m.data[j*n+i]))
		}
	}
	return r
}

// RealPart returns the matrix of real parts.
func (m *CMatrix) RealPart() *Matrix {
	r := NewMatrix(m.rows, m.cols)
	for i := range m.data {
		r.data[i] = real(m.data[i])
	}
	return r
}

// ImagPart returns the matrix of imaginary parts.
func (m *CMatrix) ImagPart() *Matrix {
	r := NewMatrix(m.rows, m.cols)
	for i := range m.data {
		r.data[i] = imag(m.data[i])
	}
	return r
}
