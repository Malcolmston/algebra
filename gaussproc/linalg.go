package gaussproc

import (
	"errors"
	"math"
)

// ErrDimensionMismatch is returned when two vectors or matrices have
// incompatible shapes for the requested operation.
var ErrDimensionMismatch = errors.New("gaussproc: dimension mismatch")

// ErrNotPositiveDefinite is returned by [Cholesky] and dependent routines when
// the supplied matrix is not symmetric positive definite.
var ErrNotPositiveDefinite = errors.New("gaussproc: matrix is not positive definite")

// ErrSingular is returned when a triangular system has a zero pivot.
var ErrSingular = errors.New("gaussproc: singular triangular matrix")

// ErrEmpty is returned when an operation requires at least one data point but
// receives none.
var ErrEmpty = errors.New("gaussproc: empty input")

// Matrix is a dense, row-major real matrix stored as a slice of rows. All rows
// are expected to have the same length.
type Matrix [][]float64

// Dot returns the Euclidean inner product of vectors a and b. It panics if the
// vectors have different lengths.
func Dot(a, b []float64) float64 {
	if len(a) != len(b) {
		panic(ErrDimensionMismatch)
	}
	var s float64
	for i := range a {
		s += a[i] * b[i]
	}
	return s
}

// Norm returns the Euclidean (L2) norm of vector a.
func Norm(a []float64) float64 {
	return math.Sqrt(SquaredNorm(a))
}

// SquaredNorm returns the squared Euclidean norm of vector a.
func SquaredNorm(a []float64) float64 {
	var s float64
	for _, v := range a {
		s += v * v
	}
	return s
}

// SquaredDistance returns the squared Euclidean distance between a and b. It
// panics if the vectors have different lengths.
func SquaredDistance(a, b []float64) float64 {
	if len(a) != len(b) {
		panic(ErrDimensionMismatch)
	}
	var s float64
	for i := range a {
		d := a[i] - b[i]
		s += d * d
	}
	return s
}

// Distance returns the Euclidean distance between vectors a and b.
func Distance(a, b []float64) float64 {
	return math.Sqrt(SquaredDistance(a, b))
}

// L1Distance returns the Manhattan (L1) distance between vectors a and b.
func L1Distance(a, b []float64) float64 {
	if len(a) != len(b) {
		panic(ErrDimensionMismatch)
	}
	var s float64
	for i := range a {
		s += math.Abs(a[i] - b[i])
	}
	return s
}

// VecAdd returns the element-wise sum a+b as a new slice.
func VecAdd(a, b []float64) []float64 {
	if len(a) != len(b) {
		panic(ErrDimensionMismatch)
	}
	out := make([]float64, len(a))
	for i := range a {
		out[i] = a[i] + b[i]
	}
	return out
}

// VecSub returns the element-wise difference a-b as a new slice.
func VecSub(a, b []float64) []float64 {
	if len(a) != len(b) {
		panic(ErrDimensionMismatch)
	}
	out := make([]float64, len(a))
	for i := range a {
		out[i] = a[i] - b[i]
	}
	return out
}

// VecScale returns the vector a scaled by the scalar s as a new slice.
func VecScale(s float64, a []float64) []float64 {
	out := make([]float64, len(a))
	for i := range a {
		out[i] = s * a[i]
	}
	return out
}

// VecCopy returns an independent copy of vector a.
func VecCopy(a []float64) []float64 {
	out := make([]float64, len(a))
	copy(out, a)
	return out
}

// VecShift returns a new slice with the scalar c added to every entry of a.
func VecShift(c float64, a []float64) []float64 {
	out := make([]float64, len(a))
	for i := range a {
		out[i] = a[i] + c
	}
	return out
}

// VectorsEqual reports whether a and b have equal length and identical entries.
func VectorsEqual(a, b []float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// MaxAbsDiff returns the largest absolute difference between corresponding
// entries of a and b. It panics if the vectors have different lengths.
func MaxAbsDiff(a, b []float64) float64 {
	if len(a) != len(b) {
		panic(ErrDimensionMismatch)
	}
	var m float64
	for i := range a {
		if d := math.Abs(a[i] - b[i]); d > m {
			m = d
		}
	}
	return m
}

// NewMatrix returns a rows-by-cols matrix with all entries set to zero.
func NewMatrix(rows, cols int) Matrix {
	m := make(Matrix, rows)
	for i := range m {
		m[i] = make([]float64, cols)
	}
	return m
}

// Identity returns the n-by-n identity matrix.
func Identity(n int) Matrix {
	m := NewMatrix(n, n)
	for i := 0; i < n; i++ {
		m[i][i] = 1
	}
	return m
}

// Diag returns a square matrix with the entries of d on its main diagonal and
// zeros elsewhere.
func Diag(d []float64) Matrix {
	n := len(d)
	m := NewMatrix(n, n)
	for i := 0; i < n; i++ {
		m[i][i] = d[i]
	}
	return m
}

// Rows returns the number of rows of m.
func (m Matrix) Rows() int { return len(m) }

// Cols returns the number of columns of m, or zero if m has no rows.
func (m Matrix) Cols() int {
	if len(m) == 0 {
		return 0
	}
	return len(m[0])
}

// At returns the entry in row i and column j of m.
func (m Matrix) At(i, j int) float64 { return m[i][j] }

// Set assigns v to the entry in row i and column j of m.
func (m Matrix) Set(i, j int, v float64) { m[i][j] = v }

// Clone returns an independent deep copy of m.
func (m Matrix) Clone() Matrix {
	out := make(Matrix, len(m))
	for i := range m {
		out[i] = make([]float64, len(m[i]))
		copy(out[i], m[i])
	}
	return out
}

// Transpose returns a new matrix that is the transpose of m.
func (m Matrix) Transpose() Matrix {
	r, c := m.Rows(), m.Cols()
	t := NewMatrix(c, r)
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			t[j][i] = m[i][j]
		}
	}
	return t
}

// Diagonal returns the main diagonal of m as a new slice. Its length is the
// smaller of the row and column counts.
func (m Matrix) Diagonal() []float64 {
	n := m.Rows()
	if c := m.Cols(); c < n {
		n = c
	}
	d := make([]float64, n)
	for i := 0; i < n; i++ {
		d[i] = m[i][i]
	}
	return d
}

// Trace returns the sum of the main-diagonal entries of square matrix m.
func (m Matrix) Trace() float64 {
	var s float64
	for _, v := range m.Diagonal() {
		s += v
	}
	return s
}

// IsSquare reports whether m has equal row and column counts.
func (m Matrix) IsSquare() bool { return m.Rows() == m.Cols() }

// IsSymmetric reports whether m is square and symmetric within the absolute
// tolerance tol.
func (m Matrix) IsSymmetric(tol float64) bool {
	if !m.IsSquare() {
		return false
	}
	n := m.Rows()
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if math.Abs(m[i][j]-m[j][i]) > tol {
				return false
			}
		}
	}
	return true
}

// MaxAbsMatrixDiff returns the largest absolute difference between
// corresponding entries of a and b. It panics if the shapes differ.
func MaxAbsMatrixDiff(a, b Matrix) float64 {
	if a.Rows() != b.Rows() || a.Cols() != b.Cols() {
		panic(ErrDimensionMismatch)
	}
	var m float64
	for i := range a {
		for j := range a[i] {
			if d := math.Abs(a[i][j] - b[i][j]); d > m {
				m = d
			}
		}
	}
	return m
}

// MatVec returns the matrix-vector product A·x. It panics if the number of
// columns of A differs from the length of x.
func MatVec(a Matrix, x []float64) []float64 {
	if a.Cols() != len(x) {
		panic(ErrDimensionMismatch)
	}
	out := make([]float64, a.Rows())
	for i := range a {
		out[i] = Dot(a[i], x)
	}
	return out
}

// MatMul returns the matrix product A·B. It panics if the number of columns of
// A differs from the number of rows of B.
func MatMul(a, b Matrix) Matrix {
	if a.Cols() != b.Rows() {
		panic(ErrDimensionMismatch)
	}
	r, inner, c := a.Rows(), b.Rows(), b.Cols()
	out := NewMatrix(r, c)
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

// OuterProduct returns the rank-one matrix a·bᵀ.
func OuterProduct(a, b []float64) Matrix {
	out := NewMatrix(len(a), len(b))
	for i := range a {
		for j := range b {
			out[i][j] = a[i] * b[j]
		}
	}
	return out
}

// MatAdd returns the element-wise sum a+b. It panics if the shapes differ.
func MatAdd(a, b Matrix) Matrix {
	if a.Rows() != b.Rows() || a.Cols() != b.Cols() {
		panic(ErrDimensionMismatch)
	}
	out := NewMatrix(a.Rows(), a.Cols())
	for i := range a {
		for j := range a[i] {
			out[i][j] = a[i][j] + b[i][j]
		}
	}
	return out
}

// MatSub returns the element-wise difference a-b. It panics if the shapes
// differ.
func MatSub(a, b Matrix) Matrix {
	if a.Rows() != b.Rows() || a.Cols() != b.Cols() {
		panic(ErrDimensionMismatch)
	}
	out := NewMatrix(a.Rows(), a.Cols())
	for i := range a {
		for j := range a[i] {
			out[i][j] = a[i][j] - b[i][j]
		}
	}
	return out
}

// MatScale returns matrix m with every entry multiplied by the scalar s.
func MatScale(s float64, m Matrix) Matrix {
	out := NewMatrix(m.Rows(), m.Cols())
	for i := range m {
		for j := range m[i] {
			out[i][j] = s * m[i][j]
		}
	}
	return out
}

// AddToDiagonal returns a copy of square matrix m with the scalar c added to
// every main-diagonal entry.
func AddToDiagonal(m Matrix, c float64) Matrix {
	out := m.Clone()
	n := out.Rows()
	for i := 0; i < n; i++ {
		out[i][i] += c
	}
	return out
}

// Cholesky computes the lower-triangular Cholesky factor L of a symmetric
// positive-definite matrix a, so that a = L·Lᵀ. It returns
// [ErrNotPositiveDefinite] if a non-positive pivot is encountered.
func Cholesky(a Matrix) (Matrix, error) {
	n := a.Rows()
	if n == 0 {
		return nil, ErrEmpty
	}
	if a.Cols() != n {
		return nil, ErrDimensionMismatch
	}
	l := NewMatrix(n, n)
	for j := 0; j < n; j++ {
		sum := a[j][j]
		for k := 0; k < j; k++ {
			sum -= l[j][k] * l[j][k]
		}
		if sum <= 0 || math.IsNaN(sum) {
			return nil, ErrNotPositiveDefinite
		}
		ljj := math.Sqrt(sum)
		l[j][j] = ljj
		for i := j + 1; i < n; i++ {
			s := a[i][j]
			for k := 0; k < j; k++ {
				s -= l[i][k] * l[j][k]
			}
			l[i][j] = s / ljj
		}
	}
	return l, nil
}

// ForwardSubstitution solves the lower-triangular system L·x = b and returns x.
// It returns [ErrSingular] if a diagonal entry of L is zero.
func ForwardSubstitution(l Matrix, b []float64) ([]float64, error) {
	n := l.Rows()
	if len(b) != n {
		return nil, ErrDimensionMismatch
	}
	x := make([]float64, n)
	for i := 0; i < n; i++ {
		s := b[i]
		for k := 0; k < i; k++ {
			s -= l[i][k] * x[k]
		}
		if l[i][i] == 0 {
			return nil, ErrSingular
		}
		x[i] = s / l[i][i]
	}
	return x, nil
}

// BackSubstitution solves the upper-triangular system U·x = b and returns x. It
// returns [ErrSingular] if a diagonal entry of U is zero.
func BackSubstitution(u Matrix, b []float64) ([]float64, error) {
	n := u.Rows()
	if len(b) != n {
		return nil, ErrDimensionMismatch
	}
	x := make([]float64, n)
	for i := n - 1; i >= 0; i-- {
		s := b[i]
		for k := i + 1; k < n; k++ {
			s -= u[i][k] * x[k]
		}
		if u[i][i] == 0 {
			return nil, ErrSingular
		}
		x[i] = s / u[i][i]
	}
	return x, nil
}

// solveLowerTranspose solves Lᵀ·x = b for a lower-triangular L, without
// forming the transpose explicitly.
func solveLowerTranspose(l Matrix, b []float64) ([]float64, error) {
	n := l.Rows()
	if len(b) != n {
		return nil, ErrDimensionMismatch
	}
	x := make([]float64, n)
	for i := n - 1; i >= 0; i-- {
		s := b[i]
		for k := i + 1; k < n; k++ {
			s -= l[k][i] * x[k]
		}
		if l[i][i] == 0 {
			return nil, ErrSingular
		}
		x[i] = s / l[i][i]
	}
	return x, nil
}

// CholeskySolve solves the symmetric positive-definite system A·x = b given the
// lower-triangular Cholesky factor L of A (as returned by [Cholesky]).
func CholeskySolve(l Matrix, b []float64) ([]float64, error) {
	y, err := ForwardSubstitution(l, b)
	if err != nil {
		return nil, err
	}
	return solveLowerTranspose(l, y)
}

// CholeskySolveMatrix solves A·X = B column by column given the Cholesky factor
// L of A and a right-hand-side matrix B, returning X.
func CholeskySolveMatrix(l Matrix, b Matrix) (Matrix, error) {
	n := l.Rows()
	if b.Rows() != n {
		return nil, ErrDimensionMismatch
	}
	cols := b.Cols()
	x := NewMatrix(n, cols)
	col := make([]float64, n)
	for j := 0; j < cols; j++ {
		for i := 0; i < n; i++ {
			col[i] = b[i][j]
		}
		sol, err := CholeskySolve(l, col)
		if err != nil {
			return nil, err
		}
		for i := 0; i < n; i++ {
			x[i][j] = sol[i]
		}
	}
	return x, nil
}

// SolveSPD solves the symmetric positive-definite system A·x = b by Cholesky
// factorisation.
func SolveSPD(a Matrix, b []float64) ([]float64, error) {
	l, err := Cholesky(a)
	if err != nil {
		return nil, err
	}
	return CholeskySolve(l, b)
}

// InvertSPD returns the inverse of a symmetric positive-definite matrix a using
// its Cholesky factorisation.
func InvertSPD(a Matrix) (Matrix, error) {
	l, err := Cholesky(a)
	if err != nil {
		return nil, err
	}
	return CholeskySolveMatrix(l, Identity(a.Rows()))
}

// LogDetFromCholesky returns log(det A) computed from the lower-triangular
// Cholesky factor L of A, as 2·Σ log L[i][i].
func LogDetFromCholesky(l Matrix) float64 {
	var s float64
	for i := 0; i < l.Rows(); i++ {
		s += math.Log(l[i][i])
	}
	return 2 * s
}

// DeterminantFromCholesky returns det A computed from the Cholesky factor L of
// A, as the square of the product of the diagonal of L.
func DeterminantFromCholesky(l Matrix) float64 {
	p := 1.0
	for i := 0; i < l.Rows(); i++ {
		p *= l[i][i]
	}
	return p * p
}

// LogDetSPD returns log(det A) for a symmetric positive-definite matrix a.
func LogDetSPD(a Matrix) (float64, error) {
	l, err := Cholesky(a)
	if err != nil {
		return 0, err
	}
	return LogDetFromCholesky(l), nil
}

// QuadraticForm returns the scalar xᵀ·A·x.
func QuadraticForm(a Matrix, x []float64) float64 {
	return Dot(x, MatVec(a, x))
}
