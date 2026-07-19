package graphspectral

import (
	"math"
	"sort"
)

// EigenResult holds the result of a symmetric eigendecomposition: the
// eigenvalues and their corresponding eigenvectors, stored as the columns of a
// matrix. The pair (Values[i], Vectors.Col(i)) satisfies A·v = λ·v.
type EigenResult struct {
	// Values holds the eigenvalues.
	Values []float64
	// Vectors holds the eigenvectors as columns; column i corresponds to
	// Values[i]. The eigenvectors are orthonormal.
	Vectors *Matrix
}

// Len returns the number of eigenpairs.
func (e *EigenResult) Len() int { return len(e.Values) }

// Vector returns the i-th eigenvector (column i of Vectors).
func (e *EigenResult) Vector(i int) []float64 { return e.Vectors.Col(i) }

// Pair returns the i-th eigenvalue together with its eigenvector.
func (e *EigenResult) Pair(i int) (float64, []float64) {
	return e.Values[i], e.Vectors.Col(i)
}

// SortAscending reorders the eigenpairs so the eigenvalues increase.
func (e *EigenResult) SortAscending() { e.sortBy(true) }

// SortDescending reorders the eigenpairs so the eigenvalues decrease.
func (e *EigenResult) SortDescending() { e.sortBy(false) }

func (e *EigenResult) sortBy(asc bool) {
	n := len(e.Values)
	idx := make([]int, n)
	for i := range idx {
		idx[i] = i
	}
	sort.SliceStable(idx, func(a, b int) bool {
		if asc {
			return e.Values[idx[a]] < e.Values[idx[b]]
		}
		return e.Values[idx[a]] > e.Values[idx[b]]
	})
	newVals := make([]float64, n)
	newVec := NewMatrix(e.Vectors.rows, n)
	for newCol, old := range idx {
		newVals[newCol] = e.Values[old]
		for r := 0; r < e.Vectors.rows; r++ {
			newVec.Set(r, newCol, e.Vectors.At(r, old))
		}
	}
	e.Values = newVals
	e.Vectors = newVec
}

// Smallest returns the smallest eigenvalue and its eigenvector.
func (e *EigenResult) Smallest() (float64, []float64) {
	i := ArgMin(e.Values)
	return e.Values[i], e.Vectors.Col(i)
}

// Largest returns the largest eigenvalue and its eigenvector.
func (e *EigenResult) Largest() (float64, []float64) {
	i := ArgMax(e.Values)
	return e.Values[i], e.Vectors.Col(i)
}

// EigenSymmetric computes all eigenvalues and eigenvectors of a real symmetric
// matrix using the cyclic Jacobi method. It returns ErrNotSquare if the matrix
// is not square and ErrNotSymmetric if it is not symmetric to within 1e-9. The
// returned eigenvectors are orthonormal; the pairs are not sorted.
func EigenSymmetric(m *Matrix) (*EigenResult, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	if !m.IsSymmetric(1e-9) {
		return nil, ErrNotSymmetric
	}
	n := m.rows
	if n == 0 {
		return &EigenResult{Values: nil, Vectors: NewMatrix(0, 0)}, nil
	}
	a := m.ToRows()
	v := make([][]float64, n)
	for i := range v {
		v[i] = make([]float64, n)
		v[i][i] = 1
	}
	rotate := func(mat [][]float64, i, j, k, l int, s, tau float64) {
		g := mat[i][j]
		h := mat[k][l]
		mat[i][j] = g - s*(h+g*tau)
		mat[k][l] = h + s*(g-h*tau)
	}
	const maxSweeps = 100
	for sweep := 0; sweep < maxSweeps; sweep++ {
		var sm float64
		for p := 0; p < n-1; p++ {
			for q := p + 1; q < n; q++ {
				sm += math.Abs(a[p][q])
			}
		}
		if sm == 0 {
			break
		}
		var thresh float64
		if sweep < 3 {
			thresh = 0.2 * sm / float64(n*n)
		}
		for p := 0; p < n-1; p++ {
			for q := p + 1; q < n; q++ {
				g := 100 * math.Abs(a[p][q])
				if sweep > 3 &&
					math.Abs(a[p][p])+g == math.Abs(a[p][p]) &&
					math.Abs(a[q][q])+g == math.Abs(a[q][q]) {
					a[p][q] = 0
					continue
				}
				if math.Abs(a[p][q]) <= thresh {
					continue
				}
				h := a[q][q] - a[p][p]
				var t float64
				if math.Abs(h)+g == math.Abs(h) {
					t = a[p][q] / h
				} else {
					theta := 0.5 * h / a[p][q]
					t = 1 / (math.Abs(theta) + math.Sqrt(1+theta*theta))
					if theta < 0 {
						t = -t
					}
				}
				c := 1 / math.Sqrt(1+t*t)
				s := t * c
				tau := s / (1 + c)
				hh := t * a[p][q]
				a[p][p] -= hh
				a[q][q] += hh
				a[p][q] = 0
				for j := 0; j < p; j++ {
					rotate(a, j, p, j, q, s, tau)
				}
				for j := p + 1; j < q; j++ {
					rotate(a, p, j, j, q, s, tau)
				}
				for j := q + 1; j < n; j++ {
					rotate(a, p, j, q, j, s, tau)
				}
				for j := 0; j < n; j++ {
					rotate(v, j, p, j, q, s, tau)
				}
			}
		}
	}
	vals := make([]float64, n)
	for i := 0; i < n; i++ {
		vals[i] = a[i][i]
	}
	vecs, err := NewMatrixFromRows(v)
	if err != nil {
		return nil, err
	}
	return &EigenResult{Values: vals, Vectors: vecs}, nil
}

// Eigenvalues returns the eigenvalues of a symmetric matrix in ascending order.
func Eigenvalues(m *Matrix) ([]float64, error) {
	e, err := EigenSymmetric(m)
	if err != nil {
		return nil, err
	}
	e.SortAscending()
	return e.Values, nil
}

// SpectralRadius returns the largest absolute eigenvalue of a symmetric matrix.
func SpectralRadius(m *Matrix) (float64, error) {
	vals, err := Eigenvalues(m)
	if err != nil {
		return 0, err
	}
	var r float64
	for _, x := range vals {
		if a := math.Abs(x); a > r {
			r = a
		}
	}
	return r, nil
}

// RayleighQuotient returns xᵀAx / xᵀx, an estimate of an eigenvalue associated
// with the direction x. It returns 0 for a zero vector or a dimension mismatch.
func RayleighQuotient(a *Matrix, x []float64) float64 {
	if len(x) != a.cols {
		return 0
	}
	ax, err := a.MulVec(x)
	if err != nil {
		return 0
	}
	den := Dot(x, x)
	if den == 0 {
		return 0
	}
	return Dot(x, ax) / den
}

// PowerIteration estimates the dominant eigenpair (largest magnitude eigenvalue)
// of a square matrix by the power method, starting from a uniform vector. It
// returns the eigenvalue, its unit eigenvector, and ErrNoConvergence if the
// Rayleigh quotient has not settled within maxIter iterations to tolerance tol.
func PowerIteration(a *Matrix, maxIter int, tol float64) (float64, []float64, error) {
	if !a.IsSquare() {
		return 0, nil, ErrNotSquare
	}
	n := a.rows
	if n == 0 {
		return 0, nil, ErrEmpty
	}
	if maxIter <= 0 {
		maxIter = 1000
	}
	x := Ones(n)
	x = Normalize(x)
	prev := math.NaN()
	for it := 0; it < maxIter; it++ {
		y, err := a.MulVec(x)
		if err != nil {
			return 0, nil, err
		}
		nrm := Norm2(y)
		if nrm == 0 {
			return 0, x, nil
		}
		for i := range y {
			y[i] /= nrm
		}
		lambda := RayleighQuotient(a, y)
		x = y
		if !math.IsNaN(prev) && math.Abs(lambda-prev) <= tol {
			// fix sign so the largest-magnitude component is positive
			return lambda, canonicalSign(x), nil
		}
		prev = lambda
	}
	return RayleighQuotient(a, x), canonicalSign(x), ErrNoConvergence
}

// canonicalSign flips a vector so that its largest-magnitude entry is positive.
func canonicalSign(v []float64) []float64 {
	idx := 0
	for i := range v {
		if math.Abs(v[i]) > math.Abs(v[idx]) {
			idx = i
		}
	}
	if len(v) > 0 && v[idx] < 0 {
		return VecScale(v, -1)
	}
	return v
}
