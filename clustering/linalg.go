package clustering

import (
	"errors"
	"math"
)

// ErrSingularMatrix is returned when a matrix operation requires a nonsingular
// matrix but the supplied matrix is (numerically) singular.
var ErrSingularMatrix = errors.New("clustering: matrix is singular")

// ErrDimensionMismatch is returned when the dimensions of the inputs to a
// linear-algebra or clustering routine are incompatible.
var ErrDimensionMismatch = errors.New("clustering: dimension mismatch")

// ErrEmptyData is returned when a routine that requires at least one sample is
// given an empty dataset.
var ErrEmptyData = errors.New("clustering: empty dataset")

// ErrInvalidK is returned when the requested number of clusters is not valid
// for the supplied data (for example k <= 0 or k greater than the number of
// samples).
var ErrInvalidK = errors.New("clustering: invalid number of clusters")

// CloneMatrix returns a deep copy of the matrix m.
func CloneMatrix(m [][]float64) [][]float64 {
	out := make([][]float64, len(m))
	for i, row := range m {
		out[i] = append([]float64(nil), row...)
	}
	return out
}

// CloneVector returns a deep copy of the vector v.
func CloneVector(v []float64) []float64 {
	return append([]float64(nil), v...)
}

// Identity returns the n x n identity matrix.
func Identity(n int) [][]float64 {
	m := make([][]float64, n)
	for i := range m {
		m[i] = make([]float64, n)
		m[i][i] = 1
	}
	return m
}

// Zeros returns a rows x cols matrix filled with zeros.
func Zeros(rows, cols int) [][]float64 {
	m := make([][]float64, rows)
	for i := range m {
		m[i] = make([]float64, cols)
	}
	return m
}

// Transpose returns the transpose of the matrix m.
func Transpose(m [][]float64) [][]float64 {
	if len(m) == 0 {
		return [][]float64{}
	}
	rows, cols := len(m), len(m[0])
	t := make([][]float64, cols)
	for j := 0; j < cols; j++ {
		t[j] = make([]float64, rows)
		for i := 0; i < rows; i++ {
			t[j][i] = m[i][j]
		}
	}
	return t
}

// MatMul returns the matrix product a*b. It returns ErrDimensionMismatch if the
// inner dimensions do not agree.
func MatMul(a, b [][]float64) ([][]float64, error) {
	if len(a) == 0 || len(b) == 0 {
		return nil, ErrDimensionMismatch
	}
	n, m, p := len(a), len(a[0]), len(b[0])
	if len(b) != m {
		return nil, ErrDimensionMismatch
	}
	out := make([][]float64, n)
	for i := 0; i < n; i++ {
		out[i] = make([]float64, p)
		for k := 0; k < m; k++ {
			aik := a[i][k]
			if aik == 0 {
				continue
			}
			brow := b[k]
			for j := 0; j < p; j++ {
				out[i][j] += aik * brow[j]
			}
		}
	}
	return out, nil
}

// MatVecMul returns the product of matrix m and column vector v.
func MatVecMul(m [][]float64, v []float64) ([]float64, error) {
	if len(m) == 0 {
		return []float64{}, nil
	}
	if len(m[0]) != len(v) {
		return nil, ErrDimensionMismatch
	}
	out := make([]float64, len(m))
	for i, row := range m {
		var s float64
		for j, x := range row {
			s += x * v[j]
		}
		out[i] = s
	}
	return out, nil
}

// luDecompose computes an LU decomposition with partial pivoting of the square
// matrix a. It returns the combined LU factors, the pivot permutation and the
// sign of the permutation (used for the determinant).
func luDecompose(a [][]float64) (lu [][]float64, piv []int, sign float64, err error) {
	n := len(a)
	if n == 0 || len(a[0]) != n {
		return nil, nil, 0, ErrDimensionMismatch
	}
	lu = CloneMatrix(a)
	piv = make([]int, n)
	for i := range piv {
		piv[i] = i
	}
	sign = 1
	for k := 0; k < n; k++ {
		// Find pivot.
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
			lu[p], lu[k] = lu[k], lu[p]
			piv[p], piv[k] = piv[k], piv[p]
			sign = -sign
		}
		for i := k + 1; i < n; i++ {
			lu[i][k] /= lu[k][k]
			factor := lu[i][k]
			for j := k + 1; j < n; j++ {
				lu[i][j] -= factor * lu[k][j]
			}
		}
	}
	return lu, piv, sign, nil
}

// Determinant returns the determinant of the square matrix a.
func Determinant(a [][]float64) (float64, error) {
	lu, _, sign, err := luDecompose(a)
	if err != nil {
		if errors.Is(err, ErrSingularMatrix) {
			return 0, nil
		}
		return 0, err
	}
	det := sign
	for i := range lu {
		det *= lu[i][i]
	}
	return det, nil
}

// Inverse returns the inverse of the square matrix a, or ErrSingularMatrix if a
// is singular.
func Inverse(a [][]float64) ([][]float64, error) {
	n := len(a)
	lu, piv, _, err := luDecompose(a)
	if err != nil {
		return nil, err
	}
	inv := make([][]float64, n)
	for i := range inv {
		inv[i] = make([]float64, n)
	}
	// Solve for each column of the identity.
	for col := 0; col < n; col++ {
		b := make([]float64, n)
		for i := 0; i < n; i++ {
			if piv[i] == col {
				b[i] = 1
			} else {
				b[i] = 0
			}
		}
		// Forward substitution (Ly = b).
		for i := 0; i < n; i++ {
			for j := 0; j < i; j++ {
				b[i] -= lu[i][j] * b[j]
			}
		}
		// Back substitution (Ux = y).
		for i := n - 1; i >= 0; i-- {
			for j := i + 1; j < n; j++ {
				b[i] -= lu[i][j] * b[j]
			}
			b[i] /= lu[i][i]
		}
		for i := 0; i < n; i++ {
			inv[i][col] = b[i]
		}
	}
	return inv, nil
}

// SolveLinearSystem solves a*x = b for x, where a is a square nonsingular
// matrix.
func SolveLinearSystem(a [][]float64, b []float64) ([]float64, error) {
	n := len(a)
	if len(b) != n {
		return nil, ErrDimensionMismatch
	}
	lu, piv, _, err := luDecompose(a)
	if err != nil {
		return nil, err
	}
	x := make([]float64, n)
	for i := 0; i < n; i++ {
		x[i] = b[piv[i]]
	}
	for i := 0; i < n; i++ {
		for j := 0; j < i; j++ {
			x[i] -= lu[i][j] * x[j]
		}
	}
	for i := n - 1; i >= 0; i-- {
		for j := i + 1; j < n; j++ {
			x[i] -= lu[i][j] * x[j]
		}
		x[i] /= lu[i][i]
	}
	return x, nil
}

// Cholesky returns the lower-triangular Cholesky factor L of a symmetric
// positive-definite matrix a, such that a = L*L^T. It returns an error if a is
// not positive definite.
func Cholesky(a [][]float64) ([][]float64, error) {
	n := len(a)
	if n == 0 || len(a[0]) != n {
		return nil, ErrDimensionMismatch
	}
	l := make([][]float64, n)
	for i := range l {
		l[i] = make([]float64, n)
	}
	for i := 0; i < n; i++ {
		for j := 0; j <= i; j++ {
			sum := a[i][j]
			for k := 0; k < j; k++ {
				sum -= l[i][k] * l[j][k]
			}
			if i == j {
				if sum <= 0 {
					return nil, errors.New("clustering: matrix is not positive definite")
				}
				l[i][j] = math.Sqrt(sum)
			} else {
				l[i][j] = sum / l[j][j]
			}
		}
	}
	return l, nil
}

// Trace returns the sum of the diagonal entries of the square matrix a.
func Trace(a [][]float64) float64 {
	var s float64
	n := len(a)
	for i := 0; i < n && i < len(a[i]); i++ {
		s += a[i][i]
	}
	return s
}
