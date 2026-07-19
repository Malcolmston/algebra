package operatortheory

import (
	"fmt"
	"math"
	"math/cmplx"
	"math/rand"
	"strings"
)

// defaultTol is the tolerance used by predicates and iterative methods when the
// caller does not supply one explicitly.
const defaultTol = 1e-9

// Matrix is a dense complex matrix stored in row-major order. It is interpreted
// throughout the package as a bounded linear operator on C^n (when square).
type Matrix struct {
	rows, cols int
	data       []complex128
}

// NewMatrix returns a rows-by-cols zero matrix. It panics if either dimension
// is negative.
func NewMatrix(rows, cols int) *Matrix {
	if rows < 0 || cols < 0 {
		panic("operatortheory: negative matrix dimension")
	}
	return &Matrix{rows: rows, cols: cols, data: make([]complex128, rows*cols)}
}

// NewMatrixFromData returns a rows-by-cols matrix filled from data in row-major
// order. It returns ErrDimensionMismatch if len(data) != rows*cols.
func NewMatrixFromData(rows, cols int, data []complex128) (*Matrix, error) {
	if rows < 0 || cols < 0 {
		return nil, ErrInvalidArgument
	}
	if len(data) != rows*cols {
		return nil, ErrDimensionMismatch
	}
	m := NewMatrix(rows, cols)
	copy(m.data, data)
	return m, nil
}

// FromRows builds a matrix from a slice of rows. All rows must have equal
// length. It returns ErrDimensionMismatch on ragged input and ErrEmpty when no
// rows are given.
func FromRows(rows [][]complex128) (*Matrix, error) {
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

// FromReal returns a complex matrix whose entries are the given real numbers
// (in row-major order) with zero imaginary part.
func FromReal(rows, cols int, data []float64) (*Matrix, error) {
	if len(data) != rows*cols {
		return nil, ErrDimensionMismatch
	}
	m := NewMatrix(rows, cols)
	for i, x := range data {
		m.data[i] = complex(x, 0)
	}
	return m, nil
}

// Identity returns the n-by-n identity operator.
func Identity(n int) *Matrix {
	m := NewMatrix(n, n)
	for i := 0; i < n; i++ {
		m.data[i*n+i] = 1
	}
	return m
}

// Zero returns the rows-by-cols zero matrix. It is a synonym for NewMatrix.
func Zero(rows, cols int) *Matrix { return NewMatrix(rows, cols) }

// Diagonal returns the square diagonal matrix with the given diagonal entries.
func Diagonal(d []complex128) *Matrix {
	n := len(d)
	m := NewMatrix(n, n)
	for i := 0; i < n; i++ {
		m.data[i*n+i] = d[i]
	}
	return m
}

// RealDiagonal returns the square diagonal matrix with the given real diagonal
// entries.
func RealDiagonal(d []float64) *Matrix {
	n := len(d)
	m := NewMatrix(n, n)
	for i := 0; i < n; i++ {
		m.data[i*n+i] = complex(d[i], 0)
	}
	return m
}

// Companion returns the companion matrix of the monic polynomial whose
// coefficients (excluding the leading 1) are given from the constant term up to
// the degree-(n-1) term: p(x) = x^n + c[n-1] x^(n-1) + ... + c[1] x + c[0].
// Its eigenvalues are the roots of p. It returns ErrEmpty for no coefficients.
func Companion(coeffs []complex128) (*Matrix, error) {
	n := len(coeffs)
	if n == 0 {
		return nil, ErrEmpty
	}
	m := NewMatrix(n, n)
	for i := 1; i < n; i++ {
		m.data[i*n+(i-1)] = 1
	}
	for i := 0; i < n; i++ {
		m.data[i*n+(n-1)] = -coeffs[i]
	}
	return m, nil
}

// JordanBlock returns the n-by-n Jordan block with eigenvalue lambda: lambda on
// the diagonal and 1 on the first superdiagonal.
func JordanBlock(lambda complex128, n int) *Matrix {
	m := NewMatrix(n, n)
	for i := 0; i < n; i++ {
		m.data[i*n+i] = lambda
		if i+1 < n {
			m.data[i*n+i+1] = 1
		}
	}
	return m
}

// Rows returns the number of rows.
func (m *Matrix) Rows() int { return m.rows }

// Cols returns the number of columns.
func (m *Matrix) Cols() int { return m.cols }

// Dims returns the number of rows and columns.
func (m *Matrix) Dims() (int, int) { return m.rows, m.cols }

// IsSquare reports whether the matrix has equal numbers of rows and columns.
func (m *Matrix) IsSquare() bool { return m.rows == m.cols }

// At returns the entry in row i and column j. It panics if the indices are out
// of range.
func (m *Matrix) At(i, j int) complex128 {
	if i < 0 || i >= m.rows || j < 0 || j >= m.cols {
		panic("operatortheory: index out of range")
	}
	return m.data[i*m.cols+j]
}

// Set stores v in row i and column j. It panics if the indices are out of
// range.
func (m *Matrix) Set(i, j int, v complex128) {
	if i < 0 || i >= m.rows || j < 0 || j >= m.cols {
		panic("operatortheory: index out of range")
	}
	m.data[i*m.cols+j] = v
}

// Clone returns an independent copy of the matrix.
func (m *Matrix) Clone() *Matrix {
	c := NewMatrix(m.rows, m.cols)
	copy(c.data, m.data)
	return c
}

// Row returns a copy of row i as a Vector.
func (m *Matrix) Row(i int) Vector {
	v := make(Vector, m.cols)
	copy(v, m.data[i*m.cols:(i+1)*m.cols])
	return v
}

// Col returns a copy of column j as a Vector.
func (m *Matrix) Col(j int) Vector {
	v := make(Vector, m.rows)
	for i := 0; i < m.rows; i++ {
		v[i] = m.data[i*m.cols+j]
	}
	return v
}

// SetRow overwrites row i with the entries of v. It panics on a length
// mismatch.
func (m *Matrix) SetRow(i int, v Vector) {
	if len(v) != m.cols {
		panic("operatortheory: row length mismatch")
	}
	copy(m.data[i*m.cols:(i+1)*m.cols], v)
}

// SetCol overwrites column j with the entries of v. It panics on a length
// mismatch.
func (m *Matrix) SetCol(j int, v Vector) {
	if len(v) != m.rows {
		panic("operatortheory: column length mismatch")
	}
	for i := 0; i < m.rows; i++ {
		m.data[i*m.cols+j] = v[i]
	}
}

// Diag returns the main diagonal of the matrix as a Vector of length
// min(rows,cols).
func (m *Matrix) Diag() Vector {
	n := m.rows
	if m.cols < n {
		n = m.cols
	}
	v := make(Vector, n)
	for i := 0; i < n; i++ {
		v[i] = m.data[i*m.cols+i]
	}
	return v
}

// Equal reports whether m and other have the same shape and agree entrywise to
// within tol.
func (m *Matrix) Equal(other *Matrix, tol float64) bool {
	if m.rows != other.rows || m.cols != other.cols {
		return false
	}
	for i := range m.data {
		if cmplx.Abs(m.data[i]-other.data[i]) > tol {
			return false
		}
	}
	return true
}

// Add returns m + b. It returns ErrDimensionMismatch on a shape mismatch.
func (m *Matrix) Add(b *Matrix) (*Matrix, error) {
	if m.rows != b.rows || m.cols != b.cols {
		return nil, ErrDimensionMismatch
	}
	r := NewMatrix(m.rows, m.cols)
	for i := range m.data {
		r.data[i] = m.data[i] + b.data[i]
	}
	return r, nil
}

// Sub returns m - b. It returns ErrDimensionMismatch on a shape mismatch.
func (m *Matrix) Sub(b *Matrix) (*Matrix, error) {
	if m.rows != b.rows || m.cols != b.cols {
		return nil, ErrDimensionMismatch
	}
	r := NewMatrix(m.rows, m.cols)
	for i := range m.data {
		r.data[i] = m.data[i] - b.data[i]
	}
	return r, nil
}

// Neg returns -m.
func (m *Matrix) Neg() *Matrix { return m.Scale(-1) }

// Scale returns the matrix s*m.
func (m *Matrix) Scale(s complex128) *Matrix {
	r := NewMatrix(m.rows, m.cols)
	for i := range m.data {
		r.data[i] = s * m.data[i]
	}
	return r
}

// ScaleReal returns the matrix s*m for a real scalar s.
func (m *Matrix) ScaleReal(s float64) *Matrix { return m.Scale(complex(s, 0)) }

// Mul returns the matrix product m*b. It returns ErrDimensionMismatch when the
// inner dimensions do not agree.
func (m *Matrix) Mul(b *Matrix) (*Matrix, error) {
	if m.cols != b.rows {
		return nil, ErrDimensionMismatch
	}
	r := NewMatrix(m.rows, b.cols)
	for i := 0; i < m.rows; i++ {
		for k := 0; k < m.cols; k++ {
			a := m.data[i*m.cols+k]
			if a == 0 {
				continue
			}
			for j := 0; j < b.cols; j++ {
				r.data[i*r.cols+j] += a * b.data[k*b.cols+j]
			}
		}
	}
	return r, nil
}

// MulVec returns the matrix-vector product m*v. It returns ErrDimensionMismatch
// when len(v) != cols.
func (m *Matrix) MulVec(v Vector) (Vector, error) {
	if len(v) != m.cols {
		return nil, ErrDimensionMismatch
	}
	r := make(Vector, m.rows)
	for i := 0; i < m.rows; i++ {
		var s complex128
		for j := 0; j < m.cols; j++ {
			s += m.data[i*m.cols+j] * v[j]
		}
		r[i] = s
	}
	return r, nil
}

// Transpose returns the transpose of m (no conjugation).
func (m *Matrix) Transpose() *Matrix {
	r := NewMatrix(m.cols, m.rows)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			r.data[j*r.cols+i] = m.data[i*m.cols+j]
		}
	}
	return r
}

// Conjugate returns the entrywise complex conjugate of m.
func (m *Matrix) Conjugate() *Matrix {
	r := NewMatrix(m.rows, m.cols)
	for i := range m.data {
		r.data[i] = cmplx.Conj(m.data[i])
	}
	return r
}

// Adjoint returns the conjugate transpose (Hermitian adjoint) m^H.
func (m *Matrix) Adjoint() *Matrix {
	r := NewMatrix(m.cols, m.rows)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			r.data[j*r.cols+i] = cmplx.Conj(m.data[i*m.cols+j])
		}
	}
	return r
}

// Trace returns the sum of the diagonal entries. It panics if the matrix is not
// square.
func (m *Matrix) Trace() complex128 {
	if !m.IsSquare() {
		panic("operatortheory: trace of non-square matrix")
	}
	var s complex128
	for i := 0; i < m.rows; i++ {
		s += m.data[i*m.cols+i]
	}
	return s
}

// HadamardProduct returns the entrywise (Schur) product of m and b.
func (m *Matrix) HadamardProduct(b *Matrix) (*Matrix, error) {
	if m.rows != b.rows || m.cols != b.cols {
		return nil, ErrDimensionMismatch
	}
	r := NewMatrix(m.rows, m.cols)
	for i := range m.data {
		r.data[i] = m.data[i] * b.data[i]
	}
	return r, nil
}

// Kron returns the Kronecker product m (x) b.
func (m *Matrix) Kron(b *Matrix) *Matrix {
	r := NewMatrix(m.rows*b.rows, m.cols*b.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			a := m.data[i*m.cols+j]
			for p := 0; p < b.rows; p++ {
				for q := 0; q < b.cols; q++ {
					ri := i*b.rows + p
					rj := j*b.cols + q
					r.data[ri*r.cols+rj] = a * b.data[p*b.cols+q]
				}
			}
		}
	}
	return r
}

// DirectSum returns the block-diagonal matrix diag(m, b).
func (m *Matrix) DirectSum(b *Matrix) *Matrix {
	r := NewMatrix(m.rows+b.rows, m.cols+b.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			r.data[i*r.cols+j] = m.data[i*m.cols+j]
		}
	}
	for i := 0; i < b.rows; i++ {
		for j := 0; j < b.cols; j++ {
			r.data[(m.rows+i)*r.cols+(m.cols+j)] = b.data[i*b.cols+j]
		}
	}
	return r
}

// Commutator returns the commutator [m,b] = m*b - b*m. It returns an error if
// the products are not defined.
func (m *Matrix) Commutator(b *Matrix) (*Matrix, error) {
	mb, err := m.Mul(b)
	if err != nil {
		return nil, err
	}
	bm, err := b.Mul(m)
	if err != nil {
		return nil, err
	}
	return mb.Sub(bm)
}

// AntiCommutator returns the anticommutator {m,b} = m*b + b*m.
func (m *Matrix) AntiCommutator(b *Matrix) (*Matrix, error) {
	mb, err := m.Mul(b)
	if err != nil {
		return nil, err
	}
	bm, err := b.Mul(m)
	if err != nil {
		return nil, err
	}
	return mb.Add(bm)
}

// Power returns m raised to the non-negative integer power k using binary
// exponentiation. Power(0) is the identity. It panics if m is not square or k
// is negative.
func (m *Matrix) Power(k int) *Matrix {
	if !m.IsSquare() {
		panic("operatortheory: power of non-square matrix")
	}
	if k < 0 {
		panic("operatortheory: negative integer power")
	}
	result := Identity(m.rows)
	base := m.Clone()
	for k > 0 {
		if k&1 == 1 {
			result, _ = result.Mul(base)
		}
		k >>= 1
		if k > 0 {
			base, _ = base.Mul(base)
		}
	}
	return result
}

// HermitianPart returns (m + m^H)/2, the Hermitian part of a square matrix.
func (m *Matrix) HermitianPart() *Matrix {
	a := m.Adjoint()
	s, _ := m.Add(a)
	return s.Scale(0.5)
}

// SkewHermitianPart returns (m - m^H)/2, the skew-Hermitian part of a square
// matrix.
func (m *Matrix) SkewHermitianPart() *Matrix {
	a := m.Adjoint()
	s, _ := m.Sub(a)
	return s.Scale(0.5)
}

// RealPart returns the entrywise real part of m as a complex matrix with zero
// imaginary part.
func (m *Matrix) RealPart() *Matrix {
	r := NewMatrix(m.rows, m.cols)
	for i := range m.data {
		r.data[i] = complex(real(m.data[i]), 0)
	}
	return r
}

// ImagPart returns the entrywise imaginary part of m as a real-valued complex
// matrix.
func (m *Matrix) ImagPart() *Matrix {
	r := NewMatrix(m.rows, m.cols)
	for i := range m.data {
		r.data[i] = complex(imag(m.data[i]), 0)
	}
	return r
}

// Submatrix returns the r1..r2-1 by c1..c2-1 block of m. It returns
// ErrOutOfRange for an invalid range.
func (m *Matrix) Submatrix(r1, r2, c1, c2 int) (*Matrix, error) {
	if r1 < 0 || c1 < 0 || r2 > m.rows || c2 > m.cols || r1 > r2 || c1 > c2 {
		return nil, ErrOutOfRange
	}
	r := NewMatrix(r2-r1, c2-c1)
	for i := r1; i < r2; i++ {
		for j := c1; j < c2; j++ {
			r.data[(i-r1)*r.cols+(j-c1)] = m.data[i*m.cols+j]
		}
	}
	return r, nil
}

// Apply returns the image m*v of the vector v under the operator, panicking on
// a dimension mismatch. It is a convenience wrapper around MulVec.
func (m *Matrix) Apply(v Vector) Vector {
	r, err := m.MulVec(v)
	if err != nil {
		panic("operatortheory: dimension mismatch in Apply")
	}
	return r
}

// String renders the matrix with each entry formatted to three decimal places.
func (m *Matrix) String() string {
	var b strings.Builder
	for i := 0; i < m.rows; i++ {
		b.WriteByte('[')
		for j := 0; j < m.cols; j++ {
			if j > 0 {
				b.WriteByte(' ')
			}
			z := m.data[i*m.cols+j]
			fmt.Fprintf(&b, "%+.3f%+.3fi", real(z), imag(z))
		}
		b.WriteByte(']')
		if i+1 < m.rows {
			b.WriteByte('\n')
		}
	}
	return b.String()
}

// RandomMatrix returns an m-by-n matrix whose entries have real and imaginary
// parts drawn independently from the standard normal distribution, using a
// deterministic generator seeded by seed.
func RandomMatrix(rows, cols int, seed int64) *Matrix {
	rng := rand.New(rand.NewSource(seed))
	m := NewMatrix(rows, cols)
	for i := range m.data {
		m.data[i] = complex(rng.NormFloat64(), rng.NormFloat64())
	}
	return m
}

// RandomHermitian returns a random n-by-n Hermitian matrix built as
// (A + A^H)/2 from a random matrix A seeded by seed.
func RandomHermitian(n int, seed int64) *Matrix {
	a := RandomMatrix(n, n, seed)
	return a.HermitianPart()
}

// RandomUnitary returns a random n-by-n unitary matrix obtained from the QR
// factorisation of a random complex matrix seeded by seed. The construction
// follows the standard recipe that yields a Haar-like distribution.
func RandomUnitary(n int, seed int64) *Matrix {
	a := RandomMatrix(n, n, seed)
	q, r := a.qrRaw()
	// Adjust the signs so that the result is a well-defined function of a.
	for j := 0; j < n; j++ {
		d := r.data[j*r.cols+j]
		var phase complex128 = 1
		if a := cmplx.Abs(d); a > 0 {
			phase = d / complex(a, 0)
		}
		for i := 0; i < n; i++ {
			q.data[i*q.cols+j] *= phase
		}
	}
	return q
}

// MaxAbs returns the largest modulus among the entries of m.
func (m *Matrix) MaxAbs() float64 {
	var mx float64
	for _, z := range m.data {
		if a := cmplx.Abs(z); a > mx {
			mx = a
		}
	}
	return mx
}

// isFinite reports whether all entries are finite.
func (m *Matrix) isFinite() bool {
	for _, z := range m.data {
		if math.IsNaN(real(z)) || math.IsInf(real(z), 0) ||
			math.IsNaN(imag(z)) || math.IsInf(imag(z), 0) {
			return false
		}
	}
	return true
}
