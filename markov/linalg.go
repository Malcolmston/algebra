package markov

import (
	"errors"
	"math"
)

// Errors returned by the linear-algebra and chain routines.
var (
	// ErrNotSquare indicates that a matrix required to be square was not.
	ErrNotSquare = errors.New("markov: matrix is not square")
	// ErrDimMismatch indicates incompatible matrix or vector dimensions.
	ErrDimMismatch = errors.New("markov: dimension mismatch")
	// ErrSingular indicates that a matrix was singular (or numerically so).
	ErrSingular = errors.New("markov: matrix is singular")
	// ErrEmpty indicates that an empty matrix or vector was supplied.
	ErrEmpty = errors.New("markov: empty input")
	// ErrNotStochastic indicates that a matrix was not row-stochastic.
	ErrNotStochastic = errors.New("markov: matrix is not row-stochastic")
	// ErrNoConvergence indicates an iterative method failed to converge.
	ErrNoConvergence = errors.New("markov: iteration did not converge")
	// ErrNotAbsorbing indicates a chain has no absorbing structure required by
	// the requested operation.
	ErrNotAbsorbing = errors.New("markov: chain is not absorbing")
	// ErrNotErgodic indicates a chain is not ergodic where ergodicity is required.
	ErrNotErgodic = errors.New("markov: chain is not ergodic")
)

// Identity returns the n×n identity matrix. It returns nil for n <= 0.
func Identity(n int) [][]float64 {
	if n <= 0 {
		return nil
	}
	m := make([][]float64, n)
	for i := range m {
		m[i] = make([]float64, n)
		m[i][i] = 1
	}
	return m
}

// Zeros returns a rows×cols matrix filled with zeros. It returns nil if either
// dimension is non-positive.
func Zeros(rows, cols int) [][]float64 {
	if rows <= 0 || cols <= 0 {
		return nil
	}
	m := make([][]float64, rows)
	for i := range m {
		m[i] = make([]float64, cols)
	}
	return m
}

// CopyMatrix returns a deep copy of a.
func CopyMatrix(a [][]float64) [][]float64 {
	if a == nil {
		return nil
	}
	b := make([][]float64, len(a))
	for i := range a {
		b[i] = make([]float64, len(a[i]))
		copy(b[i], a[i])
	}
	return b
}

// CopyVector returns a copy of v.
func CopyVector(v []float64) []float64 {
	if v == nil {
		return nil
	}
	w := make([]float64, len(v))
	copy(w, v)
	return w
}

// IsSquare reports whether a is a non-empty square matrix.
func IsSquare(a [][]float64) bool {
	n := len(a)
	if n == 0 {
		return false
	}
	for _, row := range a {
		if len(row) != n {
			return false
		}
	}
	return true
}

// Rows returns the number of rows of a.
func Rows(a [][]float64) int { return len(a) }

// Cols returns the number of columns of a (the length of its first row, or 0).
func Cols(a [][]float64) int {
	if len(a) == 0 {
		return 0
	}
	return len(a[0])
}

// Transpose returns the transpose of the (possibly rectangular) matrix a.
func Transpose(a [][]float64) [][]float64 {
	if len(a) == 0 {
		return nil
	}
	c := len(a[0])
	t := make([][]float64, c)
	for j := 0; j < c; j++ {
		t[j] = make([]float64, len(a))
		for i := range a {
			t[j][i] = a[i][j]
		}
	}
	return t
}

// MatAdd returns a+b. It returns nil if the shapes differ.
func MatAdd(a, b [][]float64) [][]float64 {
	if len(a) != len(b) {
		return nil
	}
	c := make([][]float64, len(a))
	for i := range a {
		if len(a[i]) != len(b[i]) {
			return nil
		}
		c[i] = make([]float64, len(a[i]))
		for j := range a[i] {
			c[i][j] = a[i][j] + b[i][j]
		}
	}
	return c
}

// MatSub returns a-b. It returns nil if the shapes differ.
func MatSub(a, b [][]float64) [][]float64 {
	if len(a) != len(b) {
		return nil
	}
	c := make([][]float64, len(a))
	for i := range a {
		if len(a[i]) != len(b[i]) {
			return nil
		}
		c[i] = make([]float64, len(a[i]))
		for j := range a[i] {
			c[i][j] = a[i][j] - b[i][j]
		}
	}
	return c
}

// MatScale returns s·a.
func MatScale(a [][]float64, s float64) [][]float64 {
	c := make([][]float64, len(a))
	for i := range a {
		c[i] = make([]float64, len(a[i]))
		for j := range a[i] {
			c[i][j] = s * a[i][j]
		}
	}
	return c
}

// MatMul returns the matrix product a·b. It returns nil if the inner
// dimensions do not agree.
func MatMul(a, b [][]float64) [][]float64 {
	if len(a) == 0 || len(b) == 0 {
		return nil
	}
	n, k, m := len(a), len(b), len(b[0])
	if len(a[0]) != k {
		return nil
	}
	c := make([][]float64, n)
	for i := 0; i < n; i++ {
		c[i] = make([]float64, m)
		for l := 0; l < k; l++ {
			ail := a[i][l]
			if ail == 0 {
				continue
			}
			brow := b[l]
			for j := 0; j < m; j++ {
				c[i][j] += ail * brow[j]
			}
		}
	}
	return c
}

// MatVec returns the matrix-vector product a·v (v as a column).
func MatVec(a [][]float64, v []float64) []float64 {
	if len(a) == 0 || len(a[0]) != len(v) {
		return nil
	}
	out := make([]float64, len(a))
	for i := range a {
		var s float64
		for j := range v {
			s += a[i][j] * v[j]
		}
		out[i] = s
	}
	return out
}

// VecMat returns the row-vector–matrix product v·a.
func VecMat(v []float64, a [][]float64) []float64 {
	if len(a) == 0 || len(a) != len(v) {
		return nil
	}
	m := len(a[0])
	out := make([]float64, m)
	for i := range v {
		vi := v[i]
		if vi == 0 {
			continue
		}
		for j := 0; j < m; j++ {
			out[j] += vi * a[i][j]
		}
	}
	return out
}

// MatPow returns a raised to the non-negative integer power k, computed by
// binary exponentiation. MatPow(a, 0) is the identity. It returns nil if a is
// not square or k is negative.
func MatPow(a [][]float64, k int) [][]float64 {
	if !IsSquare(a) || k < 0 {
		return nil
	}
	n := len(a)
	result := Identity(n)
	base := CopyMatrix(a)
	for k > 0 {
		if k&1 == 1 {
			result = MatMul(result, base)
		}
		k >>= 1
		if k > 0 {
			base = MatMul(base, base)
		}
	}
	return result
}

// Trace returns the sum of the diagonal entries of the square matrix a.
func Trace(a [][]float64) float64 {
	var s float64
	for i := range a {
		if i < len(a[i]) {
			s += a[i][i]
		}
	}
	return s
}

// MatEqual reports whether a and b have the same shape and all corresponding
// entries differ by at most tol in absolute value.
func MatEqual(a, b [][]float64, tol float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if len(a[i]) != len(b[i]) {
			return false
		}
		for j := range a[i] {
			if math.Abs(a[i][j]-b[i][j]) > tol {
				return false
			}
		}
	}
	return true
}

// VecEqual reports whether a and b have equal length and every corresponding
// entry differs by at most tol.
func VecEqual(a, b []float64, tol float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if math.Abs(a[i]-b[i]) > tol {
			return false
		}
	}
	return true
}

// luDecompose computes an in-place LU decomposition with partial pivoting of a
// copy of the square matrix a. It returns the combined LU matrix, the pivot
// permutation, and the sign of the permutation (±1). It errors if a is not
// square or is singular.
func luDecompose(a [][]float64) (lu [][]float64, piv []int, sign float64, err error) {
	if !IsSquare(a) {
		return nil, nil, 0, ErrNotSquare
	}
	n := len(a)
	lu = CopyMatrix(a)
	piv = make([]int, n)
	for i := range piv {
		piv[i] = i
	}
	sign = 1
	for col := 0; col < n; col++ {
		// Partial pivot: find largest magnitude in this column.
		p := col
		max := math.Abs(lu[col][col])
		for r := col + 1; r < n; r++ {
			if v := math.Abs(lu[r][col]); v > max {
				max = v
				p = r
			}
		}
		if max == 0 {
			return nil, nil, 0, ErrSingular
		}
		if p != col {
			lu[p], lu[col] = lu[col], lu[p]
			piv[p], piv[col] = piv[col], piv[p]
			sign = -sign
		}
		pivot := lu[col][col]
		for r := col + 1; r < n; r++ {
			f := lu[r][col] / pivot
			lu[r][col] = f
			for c := col + 1; c < n; c++ {
				lu[r][c] -= f * lu[col][c]
			}
		}
	}
	return lu, piv, sign, nil
}

// luSolve solves LU x = P b for a single right-hand side b given the LU
// factorization and pivot vector produced by luDecompose.
func luSolve(lu [][]float64, piv []int, b []float64) []float64 {
	n := len(lu)
	x := make([]float64, n)
	for i := 0; i < n; i++ {
		x[i] = b[piv[i]]
	}
	// Forward substitution (unit lower triangular).
	for i := 0; i < n; i++ {
		for j := 0; j < i; j++ {
			x[i] -= lu[i][j] * x[j]
		}
	}
	// Back substitution (upper triangular).
	for i := n - 1; i >= 0; i-- {
		for j := i + 1; j < n; j++ {
			x[i] -= lu[i][j] * x[j]
		}
		x[i] /= lu[i][i]
	}
	return x
}

// LUDecompose returns the LU decomposition of the square matrix a with partial
// pivoting. The returned matrix stores the unit-lower-triangular factor L
// (below the diagonal, the implied unit diagonal omitted) and the upper factor
// U (on and above the diagonal). piv is the row permutation applied to a and
// sign is the determinant sign of that permutation.
func LUDecompose(a [][]float64) (lu [][]float64, piv []int, sign float64, err error) {
	return luDecompose(a)
}

// SolveLinear solves the linear system a·x = b for a square matrix a. It
// returns ErrSingular if a is singular.
func SolveLinear(a [][]float64, b []float64) ([]float64, error) {
	if !IsSquare(a) {
		return nil, ErrNotSquare
	}
	if len(b) != len(a) {
		return nil, ErrDimMismatch
	}
	lu, piv, _, err := luDecompose(a)
	if err != nil {
		return nil, err
	}
	return luSolve(lu, piv, b), nil
}

// SolveMatrix solves a·X = B for a square matrix a and a matrix right-hand side
// B, returning X with the same number of columns as B.
func SolveMatrix(a, b [][]float64) ([][]float64, error) {
	if !IsSquare(a) {
		return nil, ErrNotSquare
	}
	if len(b) != len(a) {
		return nil, ErrDimMismatch
	}
	lu, piv, _, err := luDecompose(a)
	if err != nil {
		return nil, err
	}
	m := len(b[0])
	// Solve column by column.
	cols := make([][]float64, m)
	for j := 0; j < m; j++ {
		rhs := make([]float64, len(b))
		for i := range b {
			rhs[i] = b[i][j]
		}
		cols[j] = luSolve(lu, piv, rhs)
	}
	x := make([][]float64, len(a))
	for i := range x {
		x[i] = make([]float64, m)
		for j := 0; j < m; j++ {
			x[i][j] = cols[j][i]
		}
	}
	return x, nil
}

// MatInverse returns the inverse of the square matrix a, or ErrSingular if a is
// singular.
func MatInverse(a [][]float64) ([][]float64, error) {
	if !IsSquare(a) {
		return nil, ErrNotSquare
	}
	n := len(a)
	lu, piv, _, err := luDecompose(a)
	if err != nil {
		return nil, err
	}
	inv := make([][]float64, n)
	for i := range inv {
		inv[i] = make([]float64, n)
	}
	for j := 0; j < n; j++ {
		e := make([]float64, n)
		e[j] = 1
		col := luSolve(lu, piv, e)
		for i := 0; i < n; i++ {
			inv[i][j] = col[i]
		}
	}
	return inv, nil
}

// Determinant returns the determinant of the square matrix a. A singular
// matrix yields 0 (with a nil error).
func Determinant(a [][]float64) (float64, error) {
	if !IsSquare(a) {
		return 0, ErrNotSquare
	}
	lu, _, sign, err := luDecompose(a)
	if err != nil {
		if errors.Is(err, ErrSingular) {
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

// Dot returns the inner product of vectors a and b. It returns 0 if the lengths
// differ.
func Dot(a, b []float64) float64 {
	if len(a) != len(b) {
		return 0
	}
	var s float64
	for i := range a {
		s += a[i] * b[i]
	}
	return s
}

// VecAdd returns a+b, or nil if the lengths differ.
func VecAdd(a, b []float64) []float64 {
	if len(a) != len(b) {
		return nil
	}
	c := make([]float64, len(a))
	for i := range a {
		c[i] = a[i] + b[i]
	}
	return c
}

// VecSub returns a-b, or nil if the lengths differ.
func VecSub(a, b []float64) []float64 {
	if len(a) != len(b) {
		return nil
	}
	c := make([]float64, len(a))
	for i := range a {
		c[i] = a[i] - b[i]
	}
	return c
}

// VecScale returns s·v.
func VecScale(v []float64, s float64) []float64 {
	c := make([]float64, len(v))
	for i := range v {
		c[i] = s * v[i]
	}
	return c
}

// VecSum returns the sum of the entries of v.
func VecSum(v []float64) float64 {
	var s float64
	for _, x := range v {
		s += x
	}
	return s
}

// VecNorm1 returns the L1 norm (sum of absolute values) of v.
func VecNorm1(v []float64) float64 {
	var s float64
	for _, x := range v {
		s += math.Abs(x)
	}
	return s
}

// VecNorm2 returns the Euclidean (L2) norm of v.
func VecNorm2(v []float64) float64 {
	var s float64
	for _, x := range v {
		s += x * x
	}
	return math.Sqrt(s)
}

// VecNormInf returns the maximum absolute value among the entries of v.
func VecNormInf(v []float64) float64 {
	var m float64
	for _, x := range v {
		if a := math.Abs(x); a > m {
			m = a
		}
	}
	return m
}

// MatNorm1 returns the induced 1-norm (maximum absolute column sum) of a.
func MatNorm1(a [][]float64) float64 {
	if len(a) == 0 {
		return 0
	}
	c := len(a[0])
	var max float64
	for j := 0; j < c; j++ {
		var s float64
		for i := range a {
			s += math.Abs(a[i][j])
		}
		if s > max {
			max = s
		}
	}
	return max
}

// MatNormInf returns the induced infinity-norm (maximum absolute row sum) of a.
func MatNormInf(a [][]float64) float64 {
	var max float64
	for i := range a {
		var s float64
		for j := range a[i] {
			s += math.Abs(a[i][j])
		}
		if s > max {
			max = s
		}
	}
	return max
}

// MatNormFro returns the Frobenius norm of a.
func MatNormFro(a [][]float64) float64 {
	var s float64
	for i := range a {
		for j := range a[i] {
			s += a[i][j] * a[i][j]
		}
	}
	return math.Sqrt(s)
}
