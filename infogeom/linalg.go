package infogeom

import (
	"errors"
	"math"
)

// Sentinel errors returned by the package. They can be compared with
// errors.Is by callers that need to distinguish failure modes.
var (
	// ErrDim is returned when vectors or matrices have incompatible or empty
	// dimensions.
	ErrDim = errors.New("infogeom: incompatible dimensions")
	// ErrNotProb is returned when a slice that is required to be a probability
	// distribution has a negative entry or does not sum to one.
	ErrNotProb = errors.New("infogeom: not a probability distribution")
	// ErrDomain is returned when a parameter lies outside the domain on which
	// a quantity is defined (for example a non-positive variance).
	ErrDomain = errors.New("infogeom: argument out of domain")
	// ErrSingular is returned when a matrix that must be invertible or
	// positive definite is singular.
	ErrSingular = errors.New("infogeom: matrix is singular")
)

// probTol is the tolerance used when checking that a distribution sums to one.
const probTol = 1e-9

// NewVector returns a new zero vector of length n. It returns ErrDim when n is
// not positive.
func NewVector(n int) ([]float64, error) {
	if n <= 0 {
		return nil, ErrDim
	}
	return make([]float64, n), nil
}

// CloneVector returns an independent copy of v.
func CloneVector(v []float64) []float64 {
	out := make([]float64, len(v))
	copy(out, v)
	return out
}

// Dot returns the Euclidean inner product of x and y. It returns ErrDim when
// the vectors have different lengths.
func Dot(x, y []float64) (float64, error) {
	if len(x) != len(y) || len(x) == 0 {
		return 0, ErrDim
	}
	var s float64
	for i := range x {
		s += x[i] * y[i]
	}
	return s, nil
}

// Norm2 returns the Euclidean (L2) norm of x.
func Norm2(x []float64) float64 {
	var s float64
	for _, xi := range x {
		s += xi * xi
	}
	return math.Sqrt(s)
}

// Norm1 returns the L1 (taxicab) norm of x.
func Norm1(x []float64) float64 {
	var s float64
	for _, xi := range x {
		s += math.Abs(xi)
	}
	return s
}

// NormInf returns the L-infinity (maximum absolute value) norm of x.
func NormInf(x []float64) float64 {
	var m float64
	for _, xi := range x {
		if a := math.Abs(xi); a > m {
			m = a
		}
	}
	return m
}

// AddVectors returns the element-wise sum x+y. It returns ErrDim when the
// vectors have different lengths.
func AddVectors(x, y []float64) ([]float64, error) {
	if len(x) != len(y) {
		return nil, ErrDim
	}
	out := make([]float64, len(x))
	for i := range x {
		out[i] = x[i] + y[i]
	}
	return out, nil
}

// SubVectors returns the element-wise difference x-y. It returns ErrDim when
// the vectors have different lengths.
func SubVectors(x, y []float64) ([]float64, error) {
	if len(x) != len(y) {
		return nil, ErrDim
	}
	out := make([]float64, len(x))
	for i := range x {
		out[i] = x[i] - y[i]
	}
	return out, nil
}

// ScaleVector returns the vector a*x.
func ScaleVector(a float64, x []float64) []float64 {
	out := make([]float64, len(x))
	for i := range x {
		out[i] = a * x[i]
	}
	return out
}

// AxpyVector returns a*x + y (the BLAS "axpy" operation). It returns ErrDim
// when the vectors have different lengths.
func AxpyVector(a float64, x, y []float64) ([]float64, error) {
	if len(x) != len(y) {
		return nil, ErrDim
	}
	out := make([]float64, len(x))
	for i := range x {
		out[i] = a*x[i] + y[i]
	}
	return out, nil
}

// Sum returns the sum of the entries of x.
func Sum(x []float64) float64 {
	var s float64
	for _, xi := range x {
		s += xi
	}
	return s
}

// NewMatrix returns a new r-by-c zero matrix. It returns ErrDim when either
// dimension is not positive.
func NewMatrix(r, c int) ([][]float64, error) {
	if r <= 0 || c <= 0 {
		return nil, ErrDim
	}
	m := make([][]float64, r)
	for i := range m {
		m[i] = make([]float64, c)
	}
	return m, nil
}

// CloneMatrix returns an independent deep copy of a.
func CloneMatrix(a [][]float64) [][]float64 {
	out := make([][]float64, len(a))
	for i := range a {
		out[i] = make([]float64, len(a[i]))
		copy(out[i], a[i])
	}
	return out
}

// Identity returns the n-by-n identity matrix. It returns ErrDim when n is not
// positive.
func Identity(n int) ([][]float64, error) {
	m, err := NewMatrix(n, n)
	if err != nil {
		return nil, err
	}
	for i := 0; i < n; i++ {
		m[i][i] = 1
	}
	return m, nil
}

// Diagonal returns the square matrix whose main diagonal is d and whose
// off-diagonal entries are zero.
func Diagonal(d []float64) [][]float64 {
	n := len(d)
	m := make([][]float64, n)
	for i := 0; i < n; i++ {
		m[i] = make([]float64, n)
		m[i][i] = d[i]
	}
	return m
}

// isRectangular reports whether a is a non-empty rectangular matrix and
// returns its dimensions.
func isRectangular(a [][]float64) (int, int, bool) {
	r := len(a)
	if r == 0 {
		return 0, 0, false
	}
	c := len(a[0])
	if c == 0 {
		return 0, 0, false
	}
	for i := 1; i < r; i++ {
		if len(a[i]) != c {
			return 0, 0, false
		}
	}
	return r, c, true
}

// Dims returns the number of rows and columns of a. It returns ErrDim when a
// is empty or ragged.
func Dims(a [][]float64) (int, int, error) {
	r, c, ok := isRectangular(a)
	if !ok {
		return 0, 0, ErrDim
	}
	return r, c, nil
}

// Transpose returns the transpose of a. It returns ErrDim when a is empty or
// ragged.
func Transpose(a [][]float64) ([][]float64, error) {
	r, c, ok := isRectangular(a)
	if !ok {
		return nil, ErrDim
	}
	out := make([][]float64, c)
	for j := 0; j < c; j++ {
		out[j] = make([]float64, r)
		for i := 0; i < r; i++ {
			out[j][i] = a[i][j]
		}
	}
	return out, nil
}

// MatVec returns the matrix-vector product a*x. It returns ErrDim on a shape
// mismatch.
func MatVec(a [][]float64, x []float64) ([]float64, error) {
	r, c, ok := isRectangular(a)
	if !ok || c != len(x) {
		return nil, ErrDim
	}
	out := make([]float64, r)
	for i := 0; i < r; i++ {
		var s float64
		for j := 0; j < c; j++ {
			s += a[i][j] * x[j]
		}
		out[i] = s
	}
	return out, nil
}

// VecMat returns the row-vector-matrix product x*a. It returns ErrDim on a
// shape mismatch.
func VecMat(x []float64, a [][]float64) ([]float64, error) {
	r, c, ok := isRectangular(a)
	if !ok || r != len(x) {
		return nil, ErrDim
	}
	out := make([]float64, c)
	for j := 0; j < c; j++ {
		var s float64
		for i := 0; i < r; i++ {
			s += x[i] * a[i][j]
		}
		out[j] = s
	}
	return out, nil
}

// MatMul returns the matrix product a*b. It returns ErrDim on a shape
// mismatch.
func MatMul(a, b [][]float64) ([][]float64, error) {
	ra, ca, ok := isRectangular(a)
	if !ok {
		return nil, ErrDim
	}
	rb, cb, ok := isRectangular(b)
	if !ok || ca != rb {
		return nil, ErrDim
	}
	out := make([][]float64, ra)
	for i := 0; i < ra; i++ {
		out[i] = make([]float64, cb)
		for k := 0; k < ca; k++ {
			aik := a[i][k]
			if aik == 0 {
				continue
			}
			for j := 0; j < cb; j++ {
				out[i][j] += aik * b[k][j]
			}
		}
	}
	return out, nil
}

// QuadraticForm returns the scalar x^T A x. It returns ErrDim when A is not
// square of order len(x).
func QuadraticForm(a [][]float64, x []float64) (float64, error) {
	ax, err := MatVec(a, x)
	if err != nil {
		return 0, err
	}
	return Dot(ax, x)
}

// BilinearForm returns the scalar x^T A y. It returns ErrDim on a shape
// mismatch.
func BilinearForm(a [][]float64, x, y []float64) (float64, error) {
	ay, err := MatVec(a, y)
	if err != nil {
		return 0, err
	}
	return Dot(x, ay)
}

// Trace returns the sum of the diagonal entries of the square matrix a. It
// returns ErrDim when a is not square.
func Trace(a [][]float64) (float64, error) {
	r, c, ok := isRectangular(a)
	if !ok || r != c {
		return 0, ErrDim
	}
	var s float64
	for i := 0; i < r; i++ {
		s += a[i][i]
	}
	return s, nil
}

// AddMatrices returns the element-wise sum a+b. It returns ErrDim on a shape
// mismatch.
func AddMatrices(a, b [][]float64) ([][]float64, error) {
	ra, ca, ok := isRectangular(a)
	if !ok {
		return nil, ErrDim
	}
	rb, cb, ok := isRectangular(b)
	if !ok || ra != rb || ca != cb {
		return nil, ErrDim
	}
	out := make([][]float64, ra)
	for i := 0; i < ra; i++ {
		out[i] = make([]float64, ca)
		for j := 0; j < ca; j++ {
			out[i][j] = a[i][j] + b[i][j]
		}
	}
	return out, nil
}

// ScaleMatrix returns the matrix s*a.
func ScaleMatrix(s float64, a [][]float64) [][]float64 {
	out := make([][]float64, len(a))
	for i := range a {
		out[i] = make([]float64, len(a[i]))
		for j := range a[i] {
			out[i][j] = s * a[i][j]
		}
	}
	return out
}

// luDecompose computes an in-place LU decomposition with partial pivoting of a
// copy of a. It returns the factored matrix, the pivot permutation and the
// sign of the permutation, or ErrSingular / ErrDim on failure.
func luDecompose(a [][]float64) ([][]float64, []int, float64, error) {
	n, c, ok := isRectangular(a)
	if !ok || n != c {
		return nil, nil, 0, ErrDim
	}
	lu := CloneMatrix(a)
	piv := make([]int, n)
	for i := range piv {
		piv[i] = i
	}
	sign := 1.0
	for k := 0; k < n; k++ {
		// pivot: find the largest magnitude entry in column k.
		p := k
		max := math.Abs(lu[k][k])
		for i := k + 1; i < n; i++ {
			if v := math.Abs(lu[i][k]); v > max {
				max = v
				p = i
			}
		}
		if max == 0 {
			return nil, nil, 0, ErrSingular
		}
		if p != k {
			lu[p], lu[k] = lu[k], lu[p]
			piv[p], piv[k] = piv[k], piv[p]
			sign = -sign
		}
		for i := k + 1; i < n; i++ {
			lu[i][k] /= lu[k][k]
			f := lu[i][k]
			for j := k + 1; j < n; j++ {
				lu[i][j] -= f * lu[k][j]
			}
		}
	}
	return lu, piv, sign, nil
}

// Determinant returns the determinant of the square matrix a. It returns
// ErrDim when a is not square; a singular matrix yields determinant zero.
func Determinant(a [][]float64) (float64, error) {
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

// luSolve solves A x = b given the LU factorisation (lu, piv) of A.
func luSolve(lu [][]float64, piv []int, b []float64) []float64 {
	n := len(lu)
	x := make([]float64, n)
	for i := 0; i < n; i++ {
		x[i] = b[piv[i]]
	}
	// forward substitution (unit lower triangular).
	for i := 0; i < n; i++ {
		for j := 0; j < i; j++ {
			x[i] -= lu[i][j] * x[j]
		}
	}
	// back substitution (upper triangular).
	for i := n - 1; i >= 0; i-- {
		for j := i + 1; j < n; j++ {
			x[i] -= lu[i][j] * x[j]
		}
		x[i] /= lu[i][i]
	}
	return x
}

// Solve returns the solution x of the linear system A x = b. It returns ErrDim
// on a shape mismatch and ErrSingular when A is singular.
func Solve(a [][]float64, b []float64) ([]float64, error) {
	n, c, ok := isRectangular(a)
	if !ok || n != c || n != len(b) {
		return nil, ErrDim
	}
	lu, piv, _, err := luDecompose(a)
	if err != nil {
		return nil, err
	}
	return luSolve(lu, piv, b), nil
}

// Inverse returns the inverse of the square matrix a. It returns ErrDim when a
// is not square and ErrSingular when a is singular.
func Inverse(a [][]float64) ([][]float64, error) {
	n, c, ok := isRectangular(a)
	if !ok || n != c {
		return nil, ErrDim
	}
	lu, piv, _, err := luDecompose(a)
	if err != nil {
		return nil, err
	}
	inv := make([][]float64, n)
	e := make([]float64, n)
	for j := 0; j < n; j++ {
		for i := range e {
			e[i] = 0
		}
		e[j] = 1
		col := luSolve(lu, piv, e)
		for i := 0; i < n; i++ {
			if inv[i] == nil {
				inv[i] = make([]float64, n)
			}
			inv[i][j] = col[i]
		}
	}
	return inv, nil
}

// IsSymmetric reports whether a is a square matrix that is symmetric to within
// the given absolute tolerance.
func IsSymmetric(a [][]float64, tol float64) bool {
	n, c, ok := isRectangular(a)
	if !ok || n != c {
		return false
	}
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if math.Abs(a[i][j]-a[j][i]) > tol {
				return false
			}
		}
	}
	return true
}

// Cholesky returns the lower-triangular Cholesky factor L of a symmetric
// positive-definite matrix a, so that a = L L^T. It returns ErrDim when a is
// not square and ErrSingular when a is not positive definite.
func Cholesky(a [][]float64) ([][]float64, error) {
	n, c, ok := isRectangular(a)
	if !ok || n != c {
		return nil, ErrDim
	}
	l := make([][]float64, n)
	for i := range l {
		l[i] = make([]float64, n)
	}
	for i := 0; i < n; i++ {
		for j := 0; j <= i; j++ {
			s := a[i][j]
			for k := 0; k < j; k++ {
				s -= l[i][k] * l[j][k]
			}
			if i == j {
				if s <= 0 {
					return nil, ErrSingular
				}
				l[i][j] = math.Sqrt(s)
			} else {
				l[i][j] = s / l[j][j]
			}
		}
	}
	return l, nil
}

// IsPositiveDefinite reports whether the square matrix a is symmetric positive
// definite.
func IsPositiveDefinite(a [][]float64) bool {
	if !IsSymmetric(a, 1e-12) {
		return false
	}
	_, err := Cholesky(a)
	return err == nil
}

// LogDet returns the natural logarithm of the determinant of a symmetric
// positive-definite matrix a, computed stably from its Cholesky factor. It
// returns ErrSingular when a is not positive definite.
func LogDet(a [][]float64) (float64, error) {
	l, err := Cholesky(a)
	if err != nil {
		return 0, err
	}
	var s float64
	for i := range l {
		s += math.Log(l[i][i])
	}
	return 2 * s, nil
}
