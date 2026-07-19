package fem

import (
	"errors"
	"math"
)

// Vector is a dense real vector stored as a slice of float64 values.
type Vector []float64

// NewVector returns a zero Vector of length n. It panics if n is negative.
func NewVector(n int) Vector {
	if n < 0 {
		panic("fem: negative vector length")
	}
	return make(Vector, n)
}

// VectorOf returns a Vector containing the supplied values.
func VectorOf(vals ...float64) Vector {
	v := make(Vector, len(vals))
	copy(v, vals)
	return v
}

// Len returns the number of entries in the vector.
func (v Vector) Len() int { return len(v) }

// Clone returns an independent copy of the vector.
func (v Vector) Clone() Vector {
	c := make(Vector, len(v))
	copy(c, v)
	return c
}

// Dot returns the Euclidean inner product of v and w. It panics on a length
// mismatch.
func (v Vector) Dot(w Vector) float64 {
	if len(v) != len(w) {
		panic("fem: vector length mismatch in Dot")
	}
	var s float64
	for i := range v {
		s += v[i] * w[i]
	}
	return s
}

// Norm2 returns the Euclidean (L2) norm of the vector.
func (v Vector) Norm2() float64 { return math.Sqrt(v.Dot(v)) }

// NormInf returns the maximum absolute value of the vector entries.
func (v Vector) NormInf() float64 {
	var m float64
	for _, x := range v {
		if a := math.Abs(x); a > m {
			m = a
		}
	}
	return m
}

// Sum returns the sum of the vector entries.
func (v Vector) Sum() float64 {
	var s float64
	for _, x := range v {
		s += x
	}
	return s
}

// Scale returns a new vector equal to a*v.
func (v Vector) Scale(a float64) Vector {
	out := make(Vector, len(v))
	for i := range v {
		out[i] = a * v[i]
	}
	return out
}

// Add returns the elementwise sum v+w. It panics on a length mismatch.
func (v Vector) Add(w Vector) Vector {
	if len(v) != len(w) {
		panic("fem: vector length mismatch in Add")
	}
	out := make(Vector, len(v))
	for i := range v {
		out[i] = v[i] + w[i]
	}
	return out
}

// Sub returns the elementwise difference v-w. It panics on a length mismatch.
func (v Vector) Sub(w Vector) Vector {
	if len(v) != len(w) {
		panic("fem: vector length mismatch in Sub")
	}
	out := make(Vector, len(v))
	for i := range v {
		out[i] = v[i] - w[i]
	}
	return out
}

// AXPY returns a*x+y, the classic BLAS operation. It panics on a length
// mismatch.
func AXPY(a float64, x, y Vector) Vector {
	if len(x) != len(y) {
		panic("fem: vector length mismatch in AXPY")
	}
	out := make(Vector, len(x))
	for i := range x {
		out[i] = a*x[i] + y[i]
	}
	return out
}

// Matrix is a dense row-major real matrix.
type Matrix struct {
	rows, cols int
	data       []float64
}

// NewMatrix returns a zero r×c matrix. It panics if either dimension is
// negative.
func NewMatrix(r, c int) *Matrix {
	if r < 0 || c < 0 {
		panic("fem: negative matrix dimension")
	}
	return &Matrix{rows: r, cols: c, data: make([]float64, r*c)}
}

// Zeros is an alias for NewMatrix returning an r×c zero matrix.
func Zeros(r, c int) *Matrix { return NewMatrix(r, c) }

// Identity returns the n×n identity matrix.
func Identity(n int) *Matrix {
	m := NewMatrix(n, n)
	for i := 0; i < n; i++ {
		m.data[i*n+i] = 1
	}
	return m
}

// MatrixFromRows builds a matrix from a slice of equal-length rows.
func MatrixFromRows(rows [][]float64) *Matrix {
	r := len(rows)
	if r == 0 {
		return NewMatrix(0, 0)
	}
	c := len(rows[0])
	m := NewMatrix(r, c)
	for i := 0; i < r; i++ {
		if len(rows[i]) != c {
			panic("fem: ragged rows in MatrixFromRows")
		}
		copy(m.data[i*c:(i+1)*c], rows[i])
	}
	return m
}

// Rows returns the number of rows.
func (m *Matrix) Rows() int { return m.rows }

// Cols returns the number of columns.
func (m *Matrix) Cols() int { return m.cols }

// At returns the entry at row i, column j.
func (m *Matrix) At(i, j int) float64 { return m.data[i*m.cols+j] }

// Set assigns v to the entry at row i, column j.
func (m *Matrix) Set(i, j int, v float64) { m.data[i*m.cols+j] = v }

// Add adds v to the entry at row i, column j (scatter-add used in assembly).
func (m *Matrix) Add(i, j int, v float64) { m.data[i*m.cols+j] += v }

// Clone returns an independent copy of the matrix.
func (m *Matrix) Clone() *Matrix {
	c := NewMatrix(m.rows, m.cols)
	copy(c.data, m.data)
	return c
}

// Row returns a copy of row i as a Vector.
func (m *Matrix) Row(i int) Vector {
	out := make(Vector, m.cols)
	copy(out, m.data[i*m.cols:(i+1)*m.cols])
	return out
}

// Col returns a copy of column j as a Vector.
func (m *Matrix) Col(j int) Vector {
	out := make(Vector, m.rows)
	for i := 0; i < m.rows; i++ {
		out[i] = m.data[i*m.cols+j]
	}
	return out
}

// Scale returns a*m as a new matrix.
func (m *Matrix) Scale(a float64) *Matrix {
	out := NewMatrix(m.rows, m.cols)
	for i, x := range m.data {
		out.data[i] = a * x
	}
	return out
}

// AddMatrix returns the sum m+b. It panics on a dimension mismatch.
func (m *Matrix) AddMatrix(b *Matrix) *Matrix {
	if m.rows != b.rows || m.cols != b.cols {
		panic("fem: dimension mismatch in AddMatrix")
	}
	out := NewMatrix(m.rows, m.cols)
	for i := range m.data {
		out.data[i] = m.data[i] + b.data[i]
	}
	return out
}

// MulVec returns the matrix–vector product m*x. It panics on a length
// mismatch.
func (m *Matrix) MulVec(x Vector) Vector {
	if len(x) != m.cols {
		panic("fem: dimension mismatch in MulVec")
	}
	out := make(Vector, m.rows)
	for i := 0; i < m.rows; i++ {
		var s float64
		base := i * m.cols
		for j := 0; j < m.cols; j++ {
			s += m.data[base+j] * x[j]
		}
		out[i] = s
	}
	return out
}

// Mul returns the matrix product m*b. It panics on a dimension mismatch.
func (m *Matrix) Mul(b *Matrix) *Matrix {
	if m.cols != b.rows {
		panic("fem: dimension mismatch in Mul")
	}
	out := NewMatrix(m.rows, b.cols)
	for i := 0; i < m.rows; i++ {
		for k := 0; k < m.cols; k++ {
			a := m.data[i*m.cols+k]
			if a == 0 {
				continue
			}
			for j := 0; j < b.cols; j++ {
				out.data[i*b.cols+j] += a * b.data[k*b.cols+j]
			}
		}
	}
	return out
}

// Transpose returns the transpose of the matrix.
func (m *Matrix) Transpose() *Matrix {
	out := NewMatrix(m.cols, m.rows)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			out.data[j*m.rows+i] = m.data[i*m.cols+j]
		}
	}
	return out
}

// IsSymmetric reports whether the matrix is square and symmetric within tol.
func (m *Matrix) IsSymmetric(tol float64) bool {
	if m.rows != m.cols {
		return false
	}
	for i := 0; i < m.rows; i++ {
		for j := i + 1; j < m.cols; j++ {
			if math.Abs(m.At(i, j)-m.At(j, i)) > tol {
				return false
			}
		}
	}
	return true
}

// ErrSingular is returned when a matrix is (numerically) singular.
var ErrSingular = errors.New("fem: matrix is singular")

// ErrNotConverged is returned when an iterative solver fails to converge.
var ErrNotConverged = errors.New("fem: iterative solver did not converge")

// LU holds an in-place LU factorisation with partial pivoting.
type LU struct {
	lu   *Matrix
	piv  []int
	sign float64
}

// LUDecompose computes the LU factorisation of a square matrix using partial
// pivoting. It returns ErrSingular if a zero pivot is encountered.
func LUDecompose(a *Matrix) (*LU, error) {
	if a.rows != a.cols {
		panic("fem: LUDecompose requires a square matrix")
	}
	n := a.rows
	lu := a.Clone()
	piv := make([]int, n)
	for i := range piv {
		piv[i] = i
	}
	sign := 1.0
	for k := 0; k < n; k++ {
		p := k
		max := math.Abs(lu.data[k*n+k])
		for i := k + 1; i < n; i++ {
			if v := math.Abs(lu.data[i*n+k]); v > max {
				max = v
				p = i
			}
		}
		if max == 0 {
			return nil, ErrSingular
		}
		if p != k {
			for j := 0; j < n; j++ {
				lu.data[k*n+j], lu.data[p*n+j] = lu.data[p*n+j], lu.data[k*n+j]
			}
			piv[k], piv[p] = piv[p], piv[k]
			sign = -sign
		}
		pivot := lu.data[k*n+k]
		for i := k + 1; i < n; i++ {
			f := lu.data[i*n+k] / pivot
			lu.data[i*n+k] = f
			for j := k + 1; j < n; j++ {
				lu.data[i*n+j] -= f * lu.data[k*n+j]
			}
		}
	}
	return &LU{lu: lu, piv: piv, sign: sign}, nil
}

// Solve solves A x = b using the stored factorisation.
func (f *LU) Solve(b Vector) Vector {
	n := f.lu.rows
	if len(b) != n {
		panic("fem: dimension mismatch in LU.Solve")
	}
	x := make(Vector, n)
	for i := 0; i < n; i++ {
		x[i] = b[f.piv[i]]
	}
	// Forward substitution (unit lower triangular).
	for i := 0; i < n; i++ {
		for j := 0; j < i; j++ {
			x[i] -= f.lu.data[i*n+j] * x[j]
		}
	}
	// Back substitution (upper triangular).
	for i := n - 1; i >= 0; i-- {
		for j := i + 1; j < n; j++ {
			x[i] -= f.lu.data[i*n+j] * x[j]
		}
		x[i] /= f.lu.data[i*n+i]
	}
	return x
}

// Det returns the determinant of the factorised matrix.
func (f *LU) Det() float64 {
	n := f.lu.rows
	d := f.sign
	for i := 0; i < n; i++ {
		d *= f.lu.data[i*n+i]
	}
	return d
}

// SolveDense solves the dense linear system A x = b by LU factorisation.
func SolveDense(a *Matrix, b Vector) (Vector, error) {
	f, err := LUDecompose(a)
	if err != nil {
		return nil, err
	}
	return f.Solve(b), nil
}

// Determinant returns the determinant of a square matrix, or 0 if it is
// singular.
func Determinant(a *Matrix) float64 {
	f, err := LUDecompose(a)
	if err != nil {
		return 0
	}
	return f.Det()
}

// Inverse returns the inverse of a square matrix, or an error if it is
// singular.
func Inverse(a *Matrix) (*Matrix, error) {
	f, err := LUDecompose(a)
	if err != nil {
		return nil, err
	}
	n := a.rows
	inv := NewMatrix(n, n)
	for j := 0; j < n; j++ {
		e := make(Vector, n)
		e[j] = 1
		col := f.Solve(e)
		for i := 0; i < n; i++ {
			inv.data[i*n+j] = col[i]
		}
	}
	return inv, nil
}
