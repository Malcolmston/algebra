package grouprep

import (
	"errors"
	"fmt"
	"math"
	"math/cmplx"
	"strings"
)

// Matrix is a dense complex matrix stored row-major: m[i][j] is the entry in
// row i and column j. The zero value is not usable; construct matrices with
// [NewMatrix], [MatrixFromRows] and friends.
type Matrix [][]complex128

// NewMatrix returns an r×c matrix with every entry equal to zero. It panics if
// r or c is negative.
func NewMatrix(r, c int) Matrix {
	if r < 0 || c < 0 {
		panic("grouprep: NewMatrix requires non-negative dimensions")
	}
	m := make(Matrix, r)
	for i := range m {
		m[i] = make([]complex128, c)
	}
	return m
}

// ZeroMatrix is an alias for [NewMatrix] returning the r×c zero matrix.
func ZeroMatrix(r, c int) Matrix {
	return NewMatrix(r, c)
}

// MatrixFromRows builds a matrix from an explicit slice of rows. The input is
// deep-copied. It panics if the rows are ragged.
func MatrixFromRows(rows [][]complex128) Matrix {
	if len(rows) == 0 {
		return Matrix{}
	}
	c := len(rows[0])
	m := make(Matrix, len(rows))
	for i, row := range rows {
		if len(row) != c {
			panic("grouprep: MatrixFromRows requires rectangular input")
		}
		m[i] = append([]complex128(nil), row...)
	}
	return m
}

// MatrixFromReal builds a complex matrix whose entries are the given real
// numbers (with zero imaginary part). It panics if the input is ragged.
func MatrixFromReal(rows [][]float64) Matrix {
	if len(rows) == 0 {
		return Matrix{}
	}
	c := len(rows[0])
	m := make(Matrix, len(rows))
	for i, row := range rows {
		if len(row) != c {
			panic("grouprep: MatrixFromReal requires rectangular input")
		}
		m[i] = make([]complex128, c)
		for j, v := range row {
			m[i][j] = complex(v, 0)
		}
	}
	return m
}

// IdentityMatrix returns the n×n identity matrix. It panics if n is negative.
func IdentityMatrix(n int) Matrix {
	m := NewMatrix(n, n)
	for i := 0; i < n; i++ {
		m[i][i] = 1
	}
	return m
}

// ScalarMatrix returns the n×n diagonal matrix with z on the diagonal, i.e.
// z times the identity.
func ScalarMatrix(n int, z complex128) Matrix {
	m := NewMatrix(n, n)
	for i := 0; i < n; i++ {
		m[i][i] = z
	}
	return m
}

// DiagMatrix returns the square diagonal matrix whose diagonal entries are the
// elements of diag.
func DiagMatrix(diag []complex128) Matrix {
	n := len(diag)
	m := NewMatrix(n, n)
	for i := 0; i < n; i++ {
		m[i][i] = diag[i]
	}
	return m
}

// PermutationMatrix returns the n×n permutation matrix P of the permutation p
// (in one-line notation, p[i] the image of i), defined by P[p[i]][i] = 1. Thus
// P acts on a standard basis vector e_i by sending it to e_{p[i]}. It panics if
// p is not a permutation of {0,...,n-1}.
func PermutationMatrix(p []int) Matrix {
	n := len(p)
	seen := make([]bool, n)
	for _, v := range p {
		if v < 0 || v >= n || seen[v] {
			panic("grouprep: PermutationMatrix requires a valid permutation")
		}
		seen[v] = true
	}
	m := NewMatrix(n, n)
	for i, v := range p {
		m[v][i] = 1
	}
	return m
}

// RotationMatrix returns the 2×2 real rotation matrix through angle theta
// (counter-clockwise).
func RotationMatrix(theta float64) Matrix {
	c, s := math.Cos(theta), math.Sin(theta)
	return Matrix{
		{complex(c, 0), complex(-s, 0)},
		{complex(s, 0), complex(c, 0)},
	}
}

// ReflectionMatrix returns the 2×2 real matrix reflecting the plane across the
// line through the origin at angle theta to the x-axis.
func ReflectionMatrix(theta float64) Matrix {
	c, s := math.Cos(2*theta), math.Sin(2*theta)
	return Matrix{
		{complex(c, 0), complex(s, 0)},
		{complex(s, 0), complex(-c, 0)},
	}
}

// Rows returns the number of rows of m.
func (m Matrix) Rows() int { return len(m) }

// Cols returns the number of columns of m (zero for an empty matrix).
func (m Matrix) Cols() int {
	if len(m) == 0 {
		return 0
	}
	return len(m[0])
}

// At returns the entry m[i][j].
func (m Matrix) At(i, j int) complex128 { return m[i][j] }

// Set assigns z to entry m[i][j].
func (m Matrix) Set(i, j int, z complex128) { m[i][j] = z }

// IsSquare reports whether m has equal row and column counts.
func (m Matrix) IsSquare() bool { return m.Rows() == m.Cols() }

// Clone returns an independent deep copy of m.
func (m Matrix) Clone() Matrix {
	c := make(Matrix, len(m))
	for i := range m {
		c[i] = append([]complex128(nil), m[i]...)
	}
	return c
}

// Equal reports whether m and n have the same shape and identical entries.
func (m Matrix) Equal(n Matrix) bool {
	if m.Rows() != n.Rows() || m.Cols() != n.Cols() {
		return false
	}
	for i := range m {
		for j := range m[i] {
			if m[i][j] != n[i][j] {
				return false
			}
		}
	}
	return true
}

// ApproxEqual reports whether m and n have the same shape and entries agreeing
// to within tol in modulus.
func (m Matrix) ApproxEqual(n Matrix, tol float64) bool {
	if m.Rows() != n.Rows() || m.Cols() != n.Cols() {
		return false
	}
	for i := range m {
		for j := range m[i] {
			if cmplx.Abs(m[i][j]-n[i][j]) > tol {
				return false
			}
		}
	}
	return true
}

// Row returns a copy of row i.
func (m Matrix) Row(i int) []complex128 {
	return append([]complex128(nil), m[i]...)
}

// Col returns a copy of column j.
func (m Matrix) Col(j int) []complex128 {
	out := make([]complex128, m.Rows())
	for i := range m {
		out[i] = m[i][j]
	}
	return out
}

// Trace returns the sum of the diagonal entries of the square matrix m. It
// panics if m is not square.
func (m Matrix) Trace() complex128 {
	if !m.IsSquare() {
		panic("grouprep: Trace requires a square matrix")
	}
	var s complex128
	for i := range m {
		s += m[i][i]
	}
	return s
}

// Transpose returns the transpose of m.
func (m Matrix) Transpose() Matrix {
	t := NewMatrix(m.Cols(), m.Rows())
	for i := range m {
		for j := range m[i] {
			t[j][i] = m[i][j]
		}
	}
	return t
}

// Conjugate returns the entrywise complex conjugate of m.
func (m Matrix) Conjugate() Matrix {
	c := NewMatrix(m.Rows(), m.Cols())
	for i := range m {
		for j := range m[i] {
			c[i][j] = cmplx.Conj(m[i][j])
		}
	}
	return c
}

// ConjugateTranspose returns the conjugate transpose (Hermitian adjoint) m†.
func (m Matrix) ConjugateTranspose() Matrix {
	t := NewMatrix(m.Cols(), m.Rows())
	for i := range m {
		for j := range m[i] {
			t[j][i] = cmplx.Conj(m[i][j])
		}
	}
	return t
}

// Add returns m+n. It returns an error if the shapes differ.
func (m Matrix) Add(n Matrix) (Matrix, error) {
	if m.Rows() != n.Rows() || m.Cols() != n.Cols() {
		return nil, errors.New("grouprep: Add requires equal shapes")
	}
	out := NewMatrix(m.Rows(), m.Cols())
	for i := range m {
		for j := range m[i] {
			out[i][j] = m[i][j] + n[i][j]
		}
	}
	return out, nil
}

// Sub returns m-n. It returns an error if the shapes differ.
func (m Matrix) Sub(n Matrix) (Matrix, error) {
	if m.Rows() != n.Rows() || m.Cols() != n.Cols() {
		return nil, errors.New("grouprep: Sub requires equal shapes")
	}
	out := NewMatrix(m.Rows(), m.Cols())
	for i := range m {
		for j := range m[i] {
			out[i][j] = m[i][j] - n[i][j]
		}
	}
	return out, nil
}

// Scale returns z·m, the matrix with every entry multiplied by z.
func (m Matrix) Scale(z complex128) Matrix {
	out := NewMatrix(m.Rows(), m.Cols())
	for i := range m {
		for j := range m[i] {
			out[i][j] = z * m[i][j]
		}
	}
	return out
}

// Mul returns the matrix product m·n. It returns an error if the inner
// dimensions do not match.
func (m Matrix) Mul(n Matrix) (Matrix, error) {
	if m.Cols() != n.Rows() {
		return nil, fmt.Errorf("grouprep: Mul dimension mismatch %dx%d · %dx%d",
			m.Rows(), m.Cols(), n.Rows(), n.Cols())
	}
	out := NewMatrix(m.Rows(), n.Cols())
	for i := range m {
		for k := 0; k < m.Cols(); k++ {
			a := m[i][k]
			if a == 0 {
				continue
			}
			for j := 0; j < n.Cols(); j++ {
				out[i][j] += a * n[k][j]
			}
		}
	}
	return out, nil
}

// MulVec returns the matrix-vector product m·v. It returns an error if the
// dimensions are incompatible.
func (m Matrix) MulVec(v []complex128) ([]complex128, error) {
	if m.Cols() != len(v) {
		return nil, errors.New("grouprep: MulVec dimension mismatch")
	}
	out := make([]complex128, m.Rows())
	for i := range m {
		var s complex128
		for j := range m[i] {
			s += m[i][j] * v[j]
		}
		out[i] = s
	}
	return out, nil
}

// Pow returns m raised to the non-negative integer power k. Pow(m, 0) is the
// identity. It returns an error if m is not square.
func (m Matrix) Pow(k int) (Matrix, error) {
	if !m.IsSquare() {
		return nil, errors.New("grouprep: Pow requires a square matrix")
	}
	if k < 0 {
		return nil, errors.New("grouprep: Pow requires k >= 0")
	}
	result := IdentityMatrix(m.Rows())
	base := m.Clone()
	for k > 0 {
		if k&1 == 1 {
			result = mustMul(result, base)
		}
		k >>= 1
		if k > 0 {
			base = mustMul(base, base)
		}
	}
	return result, nil
}

// Kronecker returns the Kronecker (tensor) product m⊗n, an
// (r_m·r_n)×(c_m·c_n) matrix. It underlies the tensor product of
// representations.
func (m Matrix) Kronecker(n Matrix) Matrix {
	rm, cm := m.Rows(), m.Cols()
	rn, cn := n.Rows(), n.Cols()
	out := NewMatrix(rm*rn, cm*cn)
	for i := 0; i < rm; i++ {
		for j := 0; j < cm; j++ {
			a := m[i][j]
			for p := 0; p < rn; p++ {
				for q := 0; q < cn; q++ {
					out[i*rn+p][j*cn+q] = a * n[p][q]
				}
			}
		}
	}
	return out
}

// DirectSum returns the block-diagonal matrix diag(m, n), the direct sum m⊕n.
func (m Matrix) DirectSum(n Matrix) Matrix {
	out := NewMatrix(m.Rows()+n.Rows(), m.Cols()+n.Cols())
	for i := range m {
		for j := range m[i] {
			out[i][j] = m[i][j]
		}
	}
	ro, co := m.Rows(), m.Cols()
	for i := range n {
		for j := range n[i] {
			out[ro+i][co+j] = n[i][j]
		}
	}
	return out
}

// Submatrix returns the matrix formed by selecting the given rows and columns
// in the given order.
func (m Matrix) Submatrix(rows, cols []int) Matrix {
	out := NewMatrix(len(rows), len(cols))
	for i, r := range rows {
		for j, c := range cols {
			out[i][j] = m[r][c]
		}
	}
	return out
}

// Det returns the determinant of the square matrix m, computed by Gaussian
// elimination with partial pivoting. It panics if m is not square.
func (m Matrix) Det() complex128 {
	if !m.IsSquare() {
		panic("grouprep: Det requires a square matrix")
	}
	n := m.Rows()
	a := m.Clone()
	det := complex128(1)
	for col := 0; col < n; col++ {
		piv := -1
		best := 0.0
		for r := col; r < n; r++ {
			if v := cmplx.Abs(a[r][col]); v > best {
				best = v
				piv = r
			}
		}
		if piv == -1 {
			return 0
		}
		if piv != col {
			a[piv], a[col] = a[col], a[piv]
			det = -det
		}
		det *= a[col][col]
		inv := 1 / a[col][col]
		for r := col + 1; r < n; r++ {
			factor := a[r][col] * inv
			if factor == 0 {
				continue
			}
			for c := col; c < n; c++ {
				a[r][c] -= factor * a[col][c]
			}
		}
	}
	return det
}

// Inverse returns the inverse of the square matrix m via Gauss-Jordan
// elimination. It returns an error if m is not square or is singular.
func (m Matrix) Inverse() (Matrix, error) {
	if !m.IsSquare() {
		return nil, errors.New("grouprep: Inverse requires a square matrix")
	}
	n := m.Rows()
	a := m.Clone()
	inv := IdentityMatrix(n)
	for col := 0; col < n; col++ {
		piv := -1
		best := 0.0
		for r := col; r < n; r++ {
			if v := cmplx.Abs(a[r][col]); v > best {
				best = v
				piv = r
			}
		}
		if piv == -1 || best == 0 {
			return nil, errors.New("grouprep: Inverse of a singular matrix")
		}
		if piv != col {
			a[piv], a[col] = a[col], a[piv]
			inv[piv], inv[col] = inv[col], inv[piv]
		}
		p := a[col][col]
		for c := 0; c < n; c++ {
			a[col][c] /= p
			inv[col][c] /= p
		}
		for r := 0; r < n; r++ {
			if r == col {
				continue
			}
			f := a[r][col]
			if f == 0 {
				continue
			}
			for c := 0; c < n; c++ {
				a[r][c] -= f * a[col][c]
				inv[r][c] -= f * inv[col][c]
			}
		}
	}
	return inv, nil
}

// FrobeniusNorm returns the Frobenius norm of m, the square root of the sum of
// the squared moduli of all entries.
func (m Matrix) FrobeniusNorm() float64 {
	var s float64
	for i := range m {
		for j := range m[i] {
			a := cmplx.Abs(m[i][j])
			s += a * a
		}
	}
	return math.Sqrt(s)
}

// MaxAbsDiff returns the largest modulus |m[i][j]-n[i][j]| over all entries. It
// panics if the shapes differ.
func (m Matrix) MaxAbsDiff(n Matrix) float64 {
	if m.Rows() != n.Rows() || m.Cols() != n.Cols() {
		panic("grouprep: MaxAbsDiff requires equal shapes")
	}
	var mx float64
	for i := range m {
		for j := range m[i] {
			if v := cmplx.Abs(m[i][j] - n[i][j]); v > mx {
				mx = v
			}
		}
	}
	return mx
}

// IsIdentity reports whether m is the identity matrix to within tol.
func (m Matrix) IsIdentity(tol float64) bool {
	if !m.IsSquare() {
		return false
	}
	for i := range m {
		for j := range m[i] {
			want := complex128(0)
			if i == j {
				want = 1
			}
			if cmplx.Abs(m[i][j]-want) > tol {
				return false
			}
		}
	}
	return true
}

// IsZero reports whether every entry of m is zero to within tol.
func (m Matrix) IsZero(tol float64) bool {
	for i := range m {
		for j := range m[i] {
			if cmplx.Abs(m[i][j]) > tol {
				return false
			}
		}
	}
	return true
}

// IsDiagonal reports whether all off-diagonal entries of the square matrix m
// vanish to within tol.
func (m Matrix) IsDiagonal(tol float64) bool {
	if !m.IsSquare() {
		return false
	}
	for i := range m {
		for j := range m[i] {
			if i != j && cmplx.Abs(m[i][j]) > tol {
				return false
			}
		}
	}
	return true
}

// IsHermitian reports whether m equals its conjugate transpose to within tol.
func (m Matrix) IsHermitian(tol float64) bool {
	if !m.IsSquare() {
		return false
	}
	for i := range m {
		for j := range m[i] {
			if cmplx.Abs(m[i][j]-cmplx.Conj(m[j][i])) > tol {
				return false
			}
		}
	}
	return true
}

// IsUnitary reports whether m is unitary, i.e. m†m is the identity to within
// tol.
func (m Matrix) IsUnitary(tol float64) bool {
	if !m.IsSquare() {
		return false
	}
	prod := mustMul(m.ConjugateTranspose(), m)
	return prod.IsIdentity(tol)
}

// Round returns a copy of m with every entry rounded to the given number of
// decimal places.
func (m Matrix) Round(decimals int) Matrix {
	out := NewMatrix(m.Rows(), m.Cols())
	for i := range m {
		for j := range m[i] {
			out[i][j] = RoundC(m[i][j], decimals)
		}
	}
	return out
}

// Rank returns the numerical rank of m, the number of pivots found by Gaussian
// elimination with entries below tol treated as zero.
func (m Matrix) Rank(tol float64) int {
	a := m.Clone()
	rows, cols := a.Rows(), a.Cols()
	rank := 0
	for col := 0; col < cols && rank < rows; col++ {
		piv := -1
		best := tol
		for r := rank; r < rows; r++ {
			if v := cmplx.Abs(a[r][col]); v > best {
				best = v
				piv = r
			}
		}
		if piv == -1 {
			continue
		}
		a[piv], a[rank] = a[rank], a[piv]
		p := a[rank][col]
		for r := 0; r < rows; r++ {
			if r == rank {
				continue
			}
			f := a[r][col] / p
			for c := col; c < cols; c++ {
				a[r][c] -= f * a[rank][c]
			}
		}
		rank++
	}
	return rank
}

// String renders m with entries rounded to three decimals, one row per line.
func (m Matrix) String() string {
	var b strings.Builder
	for i := range m {
		b.WriteString("[")
		for j := range m[i] {
			if j > 0 {
				b.WriteString(" ")
			}
			z := RoundC(m[i][j], 3)
			b.WriteString(fmt.Sprintf("%g%+gi", real(z), imag(z)))
		}
		b.WriteString("]\n")
	}
	return b.String()
}

// MatMul returns the product a·b. It returns an error on a dimension mismatch;
// it is the free-function form of [Matrix.Mul].
func MatMul(a, b Matrix) (Matrix, error) { return a.Mul(b) }

// MatAdd returns a+b, erroring on a shape mismatch.
func MatAdd(a, b Matrix) (Matrix, error) { return a.Add(b) }

// MatSub returns a-b, erroring on a shape mismatch.
func MatSub(a, b Matrix) (Matrix, error) { return a.Sub(b) }

// MatKronecker returns the Kronecker product a⊗b.
func MatKronecker(a, b Matrix) Matrix { return a.Kronecker(b) }

// MatDirectSum returns the block-diagonal direct sum a⊕b.
func MatDirectSum(a, b Matrix) Matrix { return a.DirectSum(b) }

// Commutator returns the matrix commutator ab-ba. It panics if the products
// are undefined.
func Commutator(a, b Matrix) Matrix {
	ab := mustMul(a, b)
	ba := mustMul(b, a)
	out, err := ab.Sub(ba)
	if err != nil {
		panic("grouprep: Commutator requires square matrices of equal size")
	}
	return out
}

// mustMul multiplies two matrices, panicking on a dimension mismatch. It is an
// internal helper for contexts where compatibility is guaranteed.
func mustMul(a, b Matrix) Matrix {
	out, err := a.Mul(b)
	if err != nil {
		panic(err)
	}
	return out
}
