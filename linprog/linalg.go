package linprog

import (
	"errors"
	"math"
)

// ErrSingular is returned by the linear-system solvers when the coefficient
// matrix is singular (or numerically singular) so that no unique solution
// exists.
var ErrSingular = errors.New("linprog: singular matrix")

// ErrDimension is returned when the dimensions of the supplied vectors and
// matrices are inconsistent with the requested operation.
var ErrDimension = errors.New("linprog: dimension mismatch")

// Dot returns the inner product of two equal-length vectors.
// It panics if the lengths differ.
func Dot(a, b []float64) float64 {
	if len(a) != len(b) {
		panic(ErrDimension)
	}
	var s float64
	for i := range a {
		s += a[i] * b[i]
	}
	return s
}

// MatVec returns the matrix-vector product A*x where A is an m-by-n row-major
// matrix and x has length n. The result has length m. It panics on a
// dimension mismatch.
func MatVec(a [][]float64, x []float64) []float64 {
	out := make([]float64, len(a))
	for i, row := range a {
		if len(row) != len(x) {
			panic(ErrDimension)
		}
		var s float64
		for j, v := range row {
			s += v * x[j]
		}
		out[i] = s
	}
	return out
}

// MatTVec returns the product A^T*x, i.e. the matrix-vector product of the
// transpose of the m-by-n matrix A with the length-m vector x. The result has
// length n. It panics on a dimension mismatch.
func MatTVec(a [][]float64, x []float64) []float64 {
	m := len(a)
	if m == 0 {
		return nil
	}
	if len(x) != m {
		panic(ErrDimension)
	}
	n := len(a[0])
	out := make([]float64, n)
	for i := 0; i < m; i++ {
		if len(a[i]) != n {
			panic(ErrDimension)
		}
		xi := x[i]
		for j := 0; j < n; j++ {
			out[j] += a[i][j] * xi
		}
	}
	return out
}

// Transpose returns the transpose of the m-by-n row-major matrix A as an
// n-by-m matrix. A rectangular but ragged input panics.
func Transpose(a [][]float64) [][]float64 {
	m := len(a)
	if m == 0 {
		return nil
	}
	n := len(a[0])
	out := make([][]float64, n)
	for j := 0; j < n; j++ {
		out[j] = make([]float64, m)
	}
	for i := 0; i < m; i++ {
		if len(a[i]) != n {
			panic(ErrDimension)
		}
		for j := 0; j < n; j++ {
			out[j][i] = a[i][j]
		}
	}
	return out
}

// CloneMatrix returns an independent deep copy of the row-major matrix A.
func CloneMatrix(a [][]float64) [][]float64 {
	out := make([][]float64, len(a))
	for i, row := range a {
		out[i] = append([]float64(nil), row...)
	}
	return out
}

// InfNorm returns the infinity norm (maximum absolute value) of a vector.
// It returns 0 for an empty vector.
func InfNorm(v []float64) float64 {
	var m float64
	for _, x := range v {
		if a := math.Abs(x); a > m {
			m = a
		}
	}
	return m
}

// SolveLinear solves the square linear system A*x = b for x using Gaussian
// elimination with partial pivoting. A must be an n-by-n matrix and b must
// have length n. It returns [ErrSingular] if A is (numerically) singular and
// [ErrDimension] if the shapes are inconsistent. A and b are not modified.
func SolveLinear(a [][]float64, b []float64) ([]float64, error) {
	n := len(a)
	if n == 0 {
		return nil, nil
	}
	if len(b) != n {
		return nil, ErrDimension
	}
	// Build an augmented working copy.
	m := make([][]float64, n)
	for i := 0; i < n; i++ {
		if len(a[i]) != n {
			return nil, ErrDimension
		}
		m[i] = make([]float64, n+1)
		copy(m[i], a[i])
		m[i][n] = b[i]
	}
	for col := 0; col < n; col++ {
		// Partial pivot: largest magnitude in the column.
		piv := col
		best := math.Abs(m[col][col])
		for r := col + 1; r < n; r++ {
			if v := math.Abs(m[r][col]); v > best {
				best = v
				piv = r
			}
		}
		if best < 1e-14 {
			return nil, ErrSingular
		}
		m[col], m[piv] = m[piv], m[col]
		// Eliminate below.
		for r := 0; r < n; r++ {
			if r == col {
				continue
			}
			f := m[r][col] / m[col][col]
			if f == 0 {
				continue
			}
			for c := col; c <= n; c++ {
				m[r][c] -= f * m[col][c]
			}
		}
	}
	x := make([]float64, n)
	for i := 0; i < n; i++ {
		x[i] = m[i][n] / m[i][i]
	}
	return x, nil
}
