package spectralpde

import (
	"errors"
	"math"
)

// ErrDimensionMismatch is returned when operands have incompatible shapes.
var ErrDimensionMismatch = errors.New("spectralpde: dimension mismatch")

// ErrSingularMatrix is returned when a matrix is (numerically) singular.
var ErrSingularMatrix = errors.New("spectralpde: singular matrix")

// ErrInvalidArgument is returned for out-of-range or malformed arguments.
var ErrInvalidArgument = errors.New("spectralpde: invalid argument")

// Matrix is a dense, row-major real matrix used by the linear-algebra helpers.
type Matrix struct {
	Rows int
	Cols int
	Data [][]float64
}

// NewMatrix allocates an r-by-c zero matrix.
func NewMatrix(r, c int) *Matrix {
	d := make([][]float64, r)
	for i := range d {
		d[i] = make([]float64, c)
	}
	return &Matrix{Rows: r, Cols: c, Data: d}
}

// NewMatrixFrom wraps an existing row-major slice as a Matrix. The slice is
// referenced, not copied.
func NewMatrixFrom(data [][]float64) *Matrix {
	r := len(data)
	c := 0
	if r > 0 {
		c = len(data[0])
	}
	return &Matrix{Rows: r, Cols: c, Data: data}
}

// At returns the entry at row i, column j.
func (m *Matrix) At(i, j int) float64 { return m.Data[i][j] }

// Set assigns v to the entry at row i, column j.
func (m *Matrix) Set(i, j int, v float64) { m.Data[i][j] = v }

// Clone returns a deep copy of the matrix.
func (m *Matrix) Clone() *Matrix {
	out := NewMatrix(m.Rows, m.Cols)
	for i := 0; i < m.Rows; i++ {
		copy(out.Data[i], m.Data[i])
	}
	return out
}

// Zeros returns an r-by-c matrix of zeros as a plain [][]float64.
func Zeros(r, c int) [][]float64 {
	d := make([][]float64, r)
	for i := range d {
		d[i] = make([]float64, c)
	}
	return d
}

// Identity returns the n-by-n identity matrix as [][]float64.
func Identity(n int) [][]float64 {
	d := Zeros(n, n)
	for i := 0; i < n; i++ {
		d[i][i] = 1
	}
	return d
}

// Diag returns a square matrix with the given vector on its main diagonal.
func Diag(v []float64) [][]float64 {
	n := len(v)
	d := Zeros(n, n)
	for i := 0; i < n; i++ {
		d[i][i] = v[i]
	}
	return d
}

// CopyMatrix returns a deep copy of a [][]float64.
func CopyMatrix(a [][]float64) [][]float64 {
	out := make([][]float64, len(a))
	for i := range a {
		out[i] = make([]float64, len(a[i]))
		copy(out[i], a[i])
	}
	return out
}

// Transpose returns the transpose of a.
func Transpose(a [][]float64) [][]float64 {
	r := len(a)
	if r == 0 {
		return [][]float64{}
	}
	c := len(a[0])
	out := Zeros(c, r)
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			out[j][i] = a[i][j]
		}
	}
	return out
}

// MatVec returns the matrix-vector product a*x.
func MatVec(a [][]float64, x []float64) []float64 {
	r := len(a)
	out := make([]float64, r)
	for i := 0; i < r; i++ {
		var s float64
		row := a[i]
		for j := 0; j < len(row); j++ {
			s += row[j] * x[j]
		}
		out[i] = s
	}
	return out
}

// MatMul returns the matrix product a*b.
func MatMul(a, b [][]float64) [][]float64 {
	r := len(a)
	inner := len(b)
	c := 0
	if inner > 0 {
		c = len(b[0])
	}
	out := Zeros(r, c)
	for i := 0; i < r; i++ {
		for k := 0; k < inner; k++ {
			aik := a[i][k]
			if aik == 0 {
				continue
			}
			brow := b[k]
			orow := out[i]
			for j := 0; j < c; j++ {
				orow[j] += aik * brow[j]
			}
		}
	}
	return out
}

// MatAdd returns the elementwise sum a+b.
func MatAdd(a, b [][]float64) [][]float64 {
	r := len(a)
	out := make([][]float64, r)
	for i := 0; i < r; i++ {
		out[i] = make([]float64, len(a[i]))
		for j := range a[i] {
			out[i][j] = a[i][j] + b[i][j]
		}
	}
	return out
}

// MatSub returns the elementwise difference a-b.
func MatSub(a, b [][]float64) [][]float64 {
	r := len(a)
	out := make([][]float64, r)
	for i := 0; i < r; i++ {
		out[i] = make([]float64, len(a[i]))
		for j := range a[i] {
			out[i][j] = a[i][j] - b[i][j]
		}
	}
	return out
}

// MatScale returns s*a.
func MatScale(a [][]float64, s float64) [][]float64 {
	r := len(a)
	out := make([][]float64, r)
	for i := 0; i < r; i++ {
		out[i] = make([]float64, len(a[i]))
		for j := range a[i] {
			out[i][j] = s * a[i][j]
		}
	}
	return out
}

// MatPow returns a raised to the non-negative integer power p (a must be
// square). MatPow(a,0) is the identity.
func MatPow(a [][]float64, p int) [][]float64 {
	if p < 0 {
		panic(ErrInvalidArgument)
	}
	n := len(a)
	result := Identity(n)
	base := CopyMatrix(a)
	for p > 0 {
		if p&1 == 1 {
			result = MatMul(result, base)
		}
		p >>= 1
		if p > 0 {
			base = MatMul(base, base)
		}
	}
	return result
}

// Kron returns the Kronecker product of a and b.
func Kron(a, b [][]float64) [][]float64 {
	ra, ca := len(a), 0
	if ra > 0 {
		ca = len(a[0])
	}
	rb, cb := len(b), 0
	if rb > 0 {
		cb = len(b[0])
	}
	out := Zeros(ra*rb, ca*cb)
	for i := 0; i < ra; i++ {
		for j := 0; j < ca; j++ {
			aij := a[i][j]
			for k := 0; k < rb; k++ {
				for l := 0; l < cb; l++ {
					out[i*rb+k][j*cb+l] = aij * b[k][l]
				}
			}
		}
	}
	return out
}

// Trace returns the sum of the diagonal entries of a square matrix.
func Trace(a [][]float64) float64 {
	var s float64
	for i := 0; i < len(a); i++ {
		s += a[i][i]
	}
	return s
}

// FrobeniusNorm returns the Frobenius norm of a.
func FrobeniusNorm(a [][]float64) float64 {
	var s float64
	for i := range a {
		for j := range a[i] {
			s += a[i][j] * a[i][j]
		}
	}
	return math.Sqrt(s)
}

// MatMaxAbs returns the largest absolute entry of a.
func MatMaxAbs(a [][]float64) float64 {
	var m float64
	for i := range a {
		for j := range a[i] {
			if v := math.Abs(a[i][j]); v > m {
				m = v
			}
		}
	}
	return m
}

// LUDecompose computes an in-place-style LU factorization with partial
// pivoting. It returns the combined LU matrix, the pivot permutation and the
// sign of the permutation. The input is not modified.
func LUDecompose(a [][]float64) (lu [][]float64, piv []int, sign float64, err error) {
	n := len(a)
	if n == 0 || len(a[0]) != n {
		return nil, nil, 0, ErrInvalidArgument
	}
	lu = CopyMatrix(a)
	piv = make([]int, n)
	for i := range piv {
		piv[i] = i
	}
	sign = 1
	for k := 0; k < n; k++ {
		p := k
		max := math.Abs(lu[k][k])
		for i := k + 1; i < n; i++ {
			if v := math.Abs(lu[i][k]); v > max {
				max = v
				p = i
			}
		}
		if max == 0 {
			return nil, nil, 0, ErrSingularMatrix
		}
		if p != k {
			lu[k], lu[p] = lu[p], lu[k]
			piv[k], piv[p] = piv[p], piv[k]
			sign = -sign
		}
		pivot := lu[k][k]
		for i := k + 1; i < n; i++ {
			lu[i][k] /= pivot
			f := lu[i][k]
			for j := k + 1; j < n; j++ {
				lu[i][j] -= f * lu[k][j]
			}
		}
	}
	return lu, piv, sign, nil
}

// LUSolve solves A*x = b given the LU factorization and pivot from
// LUDecompose.
func LUSolve(lu [][]float64, piv []int, b []float64) []float64 {
	n := len(lu)
	x := make([]float64, n)
	for i := 0; i < n; i++ {
		x[i] = b[piv[i]]
	}
	// Forward substitution (unit lower triangular).
	for i := 0; i < n; i++ {
		s := x[i]
		for j := 0; j < i; j++ {
			s -= lu[i][j] * x[j]
		}
		x[i] = s
	}
	// Back substitution.
	for i := n - 1; i >= 0; i-- {
		s := x[i]
		for j := i + 1; j < n; j++ {
			s -= lu[i][j] * x[j]
		}
		x[i] = s / lu[i][i]
	}
	return x
}

// SolveLinearSystem solves A*x = b for a square, non-singular A.
func SolveLinearSystem(a [][]float64, b []float64) ([]float64, error) {
	if len(a) != len(b) {
		return nil, ErrDimensionMismatch
	}
	lu, piv, _, err := LUDecompose(a)
	if err != nil {
		return nil, err
	}
	return LUSolve(lu, piv, b), nil
}

// Inverse returns the inverse of a square, non-singular matrix.
func Inverse(a [][]float64) ([][]float64, error) {
	n := len(a)
	lu, piv, _, err := LUDecompose(a)
	if err != nil {
		return nil, err
	}
	inv := Zeros(n, n)
	col := make([]float64, n)
	for c := 0; c < n; c++ {
		for i := range col {
			col[i] = 0
		}
		col[c] = 1
		x := LUSolve(lu, piv, col)
		for i := 0; i < n; i++ {
			inv[i][c] = x[i]
		}
	}
	return inv, nil
}

// Determinant returns the determinant of a square matrix.
func Determinant(a [][]float64) (float64, error) {
	lu, _, sign, err := LUDecompose(a)
	if err != nil {
		if errors.Is(err, ErrSingularMatrix) {
			return 0, nil
		}
		return 0, err
	}
	det := sign
	for i := 0; i < len(lu); i++ {
		det *= lu[i][i]
	}
	return det, nil
}

// SolveMatrix solves A*X = B for X, where B (and X) have multiple columns.
func SolveMatrix(a, b [][]float64) ([][]float64, error) {
	n := len(a)
	if len(b) != n {
		return nil, ErrDimensionMismatch
	}
	lu, piv, _, err := LUDecompose(a)
	if err != nil {
		return nil, err
	}
	cols := 0
	if n > 0 {
		cols = len(b[0])
	}
	x := Zeros(n, cols)
	rhs := make([]float64, n)
	for c := 0; c < cols; c++ {
		for i := 0; i < n; i++ {
			rhs[i] = b[i][c]
		}
		sol := LUSolve(lu, piv, rhs)
		for i := 0; i < n; i++ {
			x[i][c] = sol[i]
		}
	}
	return x, nil
}

// JacobiEigenSymmetric computes all eigenvalues and eigenvectors of a real
// symmetric matrix using the cyclic Jacobi method. The returned eigenvalues
// are sorted ascending and the columns of the returned matrix are the
// corresponding orthonormal eigenvectors.
func JacobiEigenSymmetric(a [][]float64) (vals []float64, vecs [][]float64, err error) {
	n := len(a)
	if n == 0 || len(a[0]) != n {
		return nil, nil, ErrInvalidArgument
	}
	m := CopyMatrix(a)
	v := Identity(n)
	const maxSweeps = 100
	for sweep := 0; sweep < maxSweeps; sweep++ {
		off := 0.0
		for p := 0; p < n; p++ {
			for q := p + 1; q < n; q++ {
				off += m[p][q] * m[p][q]
			}
		}
		if off < 1e-30 {
			break
		}
		for p := 0; p < n; p++ {
			for q := p + 1; q < n; q++ {
				if math.Abs(m[p][q]) < 1e-300 {
					continue
				}
				theta := (m[q][q] - m[p][p]) / (2 * m[p][q])
				t := math.Copysign(1, theta) / (math.Abs(theta) + math.Sqrt(theta*theta+1))
				if theta == 0 {
					t = 1
				}
				c := 1 / math.Sqrt(t*t+1)
				s := t * c
				for k := 0; k < n; k++ {
					mkp := m[k][p]
					mkq := m[k][q]
					m[k][p] = c*mkp - s*mkq
					m[k][q] = s*mkp + c*mkq
				}
				for k := 0; k < n; k++ {
					mpk := m[p][k]
					mqk := m[q][k]
					m[p][k] = c*mpk - s*mqk
					m[q][k] = s*mpk + c*mqk
				}
				for k := 0; k < n; k++ {
					vkp := v[k][p]
					vkq := v[k][q]
					v[k][p] = c*vkp - s*vkq
					v[k][q] = s*vkp + c*vkq
				}
			}
		}
	}
	vals = make([]float64, n)
	for i := 0; i < n; i++ {
		vals[i] = m[i][i]
	}
	idx := make([]int, n)
	for i := range idx {
		idx[i] = i
	}
	// Simple selection sort on eigenvalues, permuting eigenvectors along.
	for i := 0; i < n; i++ {
		min := i
		for j := i + 1; j < n; j++ {
			if vals[idx[j]] < vals[idx[min]] {
				min = j
			}
		}
		idx[i], idx[min] = idx[min], idx[i]
	}
	sortedVals := make([]float64, n)
	vecs = Zeros(n, n)
	for c := 0; c < n; c++ {
		sortedVals[c] = vals[idx[c]]
		for r := 0; r < n; r++ {
			vecs[r][c] = v[r][idx[c]]
		}
	}
	return sortedVals, vecs, nil
}
