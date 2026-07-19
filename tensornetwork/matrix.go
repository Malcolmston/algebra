package tensornetwork

import (
	"errors"
	"fmt"
	"math"
	"math/rand"
	"sort"
)

// Matrix is a dense, row-major matrix of float64 values. It is the numerical
// linear-algebra backend used by the tensor decompositions in this package.
type Matrix struct {
	rows int
	cols int
	data []float64
}

// NewMatrix returns a new rows x cols matrix with all entries set to zero. It
// panics if rows or cols is negative.
func NewMatrix(rows, cols int) *Matrix {
	if rows < 0 || cols < 0 {
		panic("tensornetwork: negative matrix dimension")
	}
	return &Matrix{rows: rows, cols: cols, data: make([]float64, rows*cols)}
}

// NewMatrixData wraps the given row-major data as a rows x cols matrix. The
// slice is copied. It returns an error if len(data) != rows*cols.
func NewMatrixData(rows, cols int, data []float64) (*Matrix, error) {
	if rows < 0 || cols < 0 {
		return nil, errors.New("tensornetwork: negative matrix dimension")
	}
	if len(data) != rows*cols {
		return nil, fmt.Errorf("tensornetwork: data length %d != %d", len(data), rows*cols)
	}
	cp := make([]float64, len(data))
	copy(cp, data)
	return &Matrix{rows: rows, cols: cols, data: cp}, nil
}

// MatrixFromRows builds a matrix from a slice of equal-length rows. It returns
// an error if the rows are ragged.
func MatrixFromRows(rows [][]float64) (*Matrix, error) {
	if len(rows) == 0 {
		return &Matrix{}, nil
	}
	c := len(rows[0])
	m := NewMatrix(len(rows), c)
	for i, r := range rows {
		if len(r) != c {
			return nil, errors.New("tensornetwork: ragged rows")
		}
		copy(m.data[i*c:(i+1)*c], r)
	}
	return m, nil
}

// IdentityMatrix returns the n x n identity matrix.
func IdentityMatrix(n int) *Matrix {
	m := NewMatrix(n, n)
	for i := 0; i < n; i++ {
		m.data[i*n+i] = 1
	}
	return m
}

// DiagMatrix returns a square diagonal matrix whose diagonal is v.
func DiagMatrix(v []float64) *Matrix {
	n := len(v)
	m := NewMatrix(n, n)
	for i := 0; i < n; i++ {
		m.data[i*n+i] = v[i]
	}
	return m
}

// RandMatrix returns a rows x cols matrix filled with independent samples from
// the standard normal distribution, drawn from a deterministic source seeded by
// seed.
func RandMatrix(seed int64, rows, cols int) *Matrix {
	r := rand.New(rand.NewSource(seed))
	m := NewMatrix(rows, cols)
	for i := range m.data {
		m.data[i] = r.NormFloat64()
	}
	return m
}

// Rows returns the number of rows.
func (m *Matrix) Rows() int { return m.rows }

// Cols returns the number of columns.
func (m *Matrix) Cols() int { return m.cols }

// At returns the entry at row i, column j.
func (m *Matrix) At(i, j int) float64 { return m.data[i*m.cols+j] }

// Set stores v at row i, column j.
func (m *Matrix) Set(i, j int, v float64) { m.data[i*m.cols+j] = v }

// Data returns the underlying row-major slice. Callers must not assume they own
// it; use [Matrix.Clone] for an independent copy.
func (m *Matrix) Data() []float64 { return m.data }

// Clone returns an independent deep copy of m.
func (m *Matrix) Clone() *Matrix {
	cp := make([]float64, len(m.data))
	copy(cp, m.data)
	return &Matrix{rows: m.rows, cols: m.cols, data: cp}
}

// Row returns a copy of row i.
func (m *Matrix) Row(i int) []float64 {
	out := make([]float64, m.cols)
	copy(out, m.data[i*m.cols:(i+1)*m.cols])
	return out
}

// Col returns a copy of column j.
func (m *Matrix) Col(j int) []float64 {
	out := make([]float64, m.rows)
	for i := 0; i < m.rows; i++ {
		out[i] = m.data[i*m.cols+j]
	}
	return out
}

// Equal reports whether m and n have the same shape and all corresponding
// entries differ by at most tol in absolute value.
func (m *Matrix) Equal(n *Matrix, tol float64) bool {
	if m.rows != n.rows || m.cols != n.cols {
		return false
	}
	for i := range m.data {
		if math.Abs(m.data[i]-n.data[i]) > tol {
			return false
		}
	}
	return true
}

// Transpose returns the transpose of m.
func (m *Matrix) Transpose() *Matrix {
	out := NewMatrix(m.cols, m.rows)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			out.data[j*m.rows+i] = m.data[i*m.cols+j]
		}
	}
	return out
}

// Mul returns the matrix product m*n. It returns an error if the inner
// dimensions do not agree.
func (m *Matrix) Mul(n *Matrix) (*Matrix, error) {
	if m.cols != n.rows {
		return nil, fmt.Errorf("tensornetwork: shape mismatch %dx%d * %dx%d", m.rows, m.cols, n.rows, n.cols)
	}
	out := NewMatrix(m.rows, n.cols)
	for i := 0; i < m.rows; i++ {
		for k := 0; k < m.cols; k++ {
			a := m.data[i*m.cols+k]
			if a == 0 {
				continue
			}
			for j := 0; j < n.cols; j++ {
				out.data[i*n.cols+j] += a * n.data[k*n.cols+j]
			}
		}
	}
	return out, nil
}

// MatMul returns a*b or panics on a dimension mismatch. It is a convenience
// wrapper around [Matrix.Mul] for callers that treat a mismatch as a bug.
func MatMul(a, b *Matrix) *Matrix {
	c, err := a.Mul(b)
	if err != nil {
		panic(err)
	}
	return c
}

// MulVec returns the matrix-vector product m*v. It returns an error if
// len(v) != m.Cols().
func (m *Matrix) MulVec(v []float64) ([]float64, error) {
	if len(v) != m.cols {
		return nil, fmt.Errorf("tensornetwork: length %d != cols %d", len(v), m.cols)
	}
	out := make([]float64, m.rows)
	for i := 0; i < m.rows; i++ {
		s := 0.0
		for j := 0; j < m.cols; j++ {
			s += m.data[i*m.cols+j] * v[j]
		}
		out[i] = s
	}
	return out, nil
}

// Add returns the entrywise sum m+n. It returns an error on a shape mismatch.
func (m *Matrix) Add(n *Matrix) (*Matrix, error) {
	if m.rows != n.rows || m.cols != n.cols {
		return nil, errors.New("tensornetwork: shape mismatch in Add")
	}
	out := m.Clone()
	for i := range out.data {
		out.data[i] += n.data[i]
	}
	return out, nil
}

// Sub returns the entrywise difference m-n. It returns an error on a shape
// mismatch.
func (m *Matrix) Sub(n *Matrix) (*Matrix, error) {
	if m.rows != n.rows || m.cols != n.cols {
		return nil, errors.New("tensornetwork: shape mismatch in Sub")
	}
	out := m.Clone()
	for i := range out.data {
		out.data[i] -= n.data[i]
	}
	return out, nil
}

// Scale returns m with every entry multiplied by s.
func (m *Matrix) Scale(s float64) *Matrix {
	out := m.Clone()
	for i := range out.data {
		out.data[i] *= s
	}
	return out
}

// Hadamard returns the entrywise (Hadamard) product m∘n. It returns an error on
// a shape mismatch.
func (m *Matrix) Hadamard(n *Matrix) (*Matrix, error) {
	if m.rows != n.rows || m.cols != n.cols {
		return nil, errors.New("tensornetwork: shape mismatch in Hadamard")
	}
	out := m.Clone()
	for i := range out.data {
		out.data[i] *= n.data[i]
	}
	return out, nil
}

// Diag returns the main diagonal of m as a slice of length min(rows, cols).
func (m *Matrix) Diag() []float64 {
	n := m.rows
	if m.cols < n {
		n = m.cols
	}
	out := make([]float64, n)
	for i := 0; i < n; i++ {
		out[i] = m.data[i*m.cols+i]
	}
	return out
}

// Trace returns the sum of the diagonal entries of m.
func (m *Matrix) Trace() float64 {
	s := 0.0
	for _, d := range m.Diag() {
		s += d
	}
	return s
}

// FrobeniusNorm returns the Frobenius norm of m, sqrt of the sum of squares of
// all entries.
func (m *Matrix) FrobeniusNorm() float64 {
	s := 0.0
	for _, v := range m.data {
		s += v * v
	}
	return math.Sqrt(s)
}

// SubMatrix returns the r0:r1 x c0:c1 block of m (half-open ranges) as a fresh
// matrix.
func (m *Matrix) SubMatrix(r0, r1, c0, c1 int) *Matrix {
	out := NewMatrix(r1-r0, c1-c0)
	for i := r0; i < r1; i++ {
		for j := c0; j < c1; j++ {
			out.data[(i-r0)*out.cols+(j-c0)] = m.data[i*m.cols+j]
		}
	}
	return out
}

// KroneckerMatrix returns the Kronecker product a⊗b, an (a.rows*b.rows) x
// (a.cols*b.cols) matrix.
func KroneckerMatrix(a, b *Matrix) *Matrix {
	out := NewMatrix(a.rows*b.rows, a.cols*b.cols)
	for ai := 0; ai < a.rows; ai++ {
		for aj := 0; aj < a.cols; aj++ {
			av := a.data[ai*a.cols+aj]
			for bi := 0; bi < b.rows; bi++ {
				for bj := 0; bj < b.cols; bj++ {
					ri := ai*b.rows + bi
					rj := aj*b.cols + bj
					out.data[ri*out.cols+rj] = av * b.data[bi*b.cols+bj]
				}
			}
		}
	}
	return out
}

// KhatriRaoMatrix returns the column-wise Khatri-Rao product a⊙b. The two
// matrices must have the same number of columns; the result has
// a.rows*b.rows rows and that many columns, with the first factor varying
// slowest within each column. It returns an error on a column mismatch.
func KhatriRaoMatrix(a, b *Matrix) (*Matrix, error) {
	if a.cols != b.cols {
		return nil, fmt.Errorf("tensornetwork: Khatri-Rao column mismatch %d != %d", a.cols, b.cols)
	}
	out := NewMatrix(a.rows*b.rows, a.cols)
	for c := 0; c < a.cols; c++ {
		for ai := 0; ai < a.rows; ai++ {
			av := a.data[ai*a.cols+c]
			for bi := 0; bi < b.rows; bi++ {
				out.data[(ai*b.rows+bi)*out.cols+c] = av * b.data[bi*b.cols+c]
			}
		}
	}
	return out, nil
}

// KhatriRao returns the Khatri-Rao product of a non-empty list of matrices, all
// with the same number of columns, combined left to right so that the first
// matrix varies slowest. It returns an error if the list is empty or the column
// counts disagree.
func KhatriRao(mats ...*Matrix) (*Matrix, error) {
	if len(mats) == 0 {
		return nil, errors.New("tensornetwork: KhatriRao needs at least one matrix")
	}
	acc := mats[0]
	for _, m := range mats[1:] {
		var err error
		acc, err = KhatriRaoMatrix(acc, m)
		if err != nil {
			return nil, err
		}
	}
	return acc, nil
}

// EigSym computes the eigenvalues and eigenvectors of a real symmetric matrix
// using the cyclic Jacobi method. The eigenvalues are returned in descending
// order and the columns of the returned matrix are the corresponding
// orthonormal eigenvectors. It returns an error if m is not square.
func (m *Matrix) EigSym() (values []float64, vectors *Matrix, err error) {
	if m.rows != m.cols {
		return nil, nil, errors.New("tensornetwork: EigSym requires a square matrix")
	}
	n := m.rows
	a := m.Clone()
	v := IdentityMatrix(n)
	const maxSweeps = 100
	for sweep := 0; sweep < maxSweeps; sweep++ {
		off := 0.0
		for p := 0; p < n; p++ {
			for q := p + 1; q < n; q++ {
				off += a.data[p*n+q] * a.data[p*n+q]
			}
		}
		if off < 1e-30 {
			break
		}
		for p := 0; p < n; p++ {
			for q := p + 1; q < n; q++ {
				apq := a.data[p*n+q]
				if math.Abs(apq) < 1e-300 {
					continue
				}
				app := a.data[p*n+p]
				aqq := a.data[q*n+q]
				theta := (aqq - app) / (2 * apq)
				t := math.Copysign(1, theta) / (math.Abs(theta) + math.Sqrt(theta*theta+1))
				if theta == 0 {
					t = 1
				}
				c := 1 / math.Sqrt(t*t+1)
				s := t * c
				for k := 0; k < n; k++ {
					akp := a.data[k*n+p]
					akq := a.data[k*n+q]
					a.data[k*n+p] = c*akp - s*akq
					a.data[k*n+q] = s*akp + c*akq
				}
				for k := 0; k < n; k++ {
					apk := a.data[p*n+k]
					aqk := a.data[q*n+k]
					a.data[p*n+k] = c*apk - s*aqk
					a.data[q*n+k] = s*apk + c*aqk
				}
				for k := 0; k < n; k++ {
					vkp := v.data[k*n+p]
					vkq := v.data[k*n+q]
					v.data[k*n+p] = c*vkp - s*vkq
					v.data[k*n+q] = s*vkp + c*vkq
				}
			}
		}
	}
	vals := make([]float64, n)
	for i := 0; i < n; i++ {
		vals[i] = a.data[i*n+i]
	}
	idx := make([]int, n)
	for i := range idx {
		idx[i] = i
	}
	sort.Slice(idx, func(i, j int) bool { return vals[idx[i]] > vals[idx[j]] })
	outVals := make([]float64, n)
	outVec := NewMatrix(n, n)
	for newc, oldc := range idx {
		outVals[newc] = vals[oldc]
		for r := 0; r < n; r++ {
			outVec.data[r*n+newc] = v.data[r*n+oldc]
		}
	}
	return outVals, outVec, nil
}

// SVD computes a thin singular value decomposition of m using one-sided Jacobi
// rotations. It returns matrices U (rows x k) and V (cols x k) with orthonormal
// columns and singular values s (length k, descending) such that
// m = U * diag(s) * Vᵀ, where k = min(rows, cols).
func (m *Matrix) SVD() (U *Matrix, s []float64, V *Matrix, err error) {
	if m.rows >= m.cols {
		return svdTall(m)
	}
	// Compute SVD of the transpose and swap the roles of U and V.
	Ut, st, Vt, e := svdTall(m.Transpose())
	if e != nil {
		return nil, nil, nil, e
	}
	return Vt, st, Ut, nil
}

// svdTall computes the SVD of a matrix with rows >= cols by one-sided Jacobi.
func svdTall(m *Matrix) (*Matrix, []float64, *Matrix, error) {
	rows, cols := m.rows, m.cols
	u := m.Clone() // becomes U*Sigma
	v := IdentityMatrix(cols)
	const maxSweeps = 60
	for sweep := 0; sweep < maxSweeps; sweep++ {
		changed := false
		for p := 0; p < cols-1; p++ {
			for q := p + 1; q < cols; q++ {
				var alpha, beta, gamma float64
				for i := 0; i < rows; i++ {
					up := u.data[i*cols+p]
					uq := u.data[i*cols+q]
					alpha += up * up
					beta += uq * uq
					gamma += up * uq
				}
				if math.Abs(gamma) <= 1e-15*math.Sqrt(alpha*beta) || gamma == 0 {
					continue
				}
				changed = true
				zeta := (beta - alpha) / (2 * gamma)
				t := math.Copysign(1, zeta) / (math.Abs(zeta) + math.Sqrt(1+zeta*zeta))
				c := 1 / math.Sqrt(1+t*t)
				s := c * t
				for i := 0; i < rows; i++ {
					up := u.data[i*cols+p]
					uq := u.data[i*cols+q]
					u.data[i*cols+p] = c*up - s*uq
					u.data[i*cols+q] = s*up + c*uq
				}
				for i := 0; i < cols; i++ {
					vp := v.data[i*cols+p]
					vq := v.data[i*cols+q]
					v.data[i*cols+p] = c*vp - s*vq
					v.data[i*cols+q] = s*vp + c*vq
				}
			}
		}
		if !changed {
			break
		}
	}
	sigma := make([]float64, cols)
	for j := 0; j < cols; j++ {
		var nrm float64
		for i := 0; i < rows; i++ {
			nrm += u.data[i*cols+j] * u.data[i*cols+j]
		}
		sigma[j] = math.Sqrt(nrm)
	}
	idx := make([]int, cols)
	for i := range idx {
		idx[i] = i
	}
	sort.Slice(idx, func(i, j int) bool { return sigma[idx[i]] > sigma[idx[j]] })
	U := NewMatrix(rows, cols)
	V := NewMatrix(cols, cols)
	s := make([]float64, cols)
	for newc, oldc := range idx {
		s[newc] = sigma[oldc]
		if sigma[oldc] > 1e-300 {
			for i := 0; i < rows; i++ {
				U.data[i*cols+newc] = u.data[i*cols+oldc] / sigma[oldc]
			}
		} else {
			// Degenerate column: leave as a zero column (fine for thin SVD use).
			for i := 0; i < rows; i++ {
				U.data[i*cols+newc] = 0
			}
		}
		for i := 0; i < cols; i++ {
			V.data[i*cols+newc] = v.data[i*cols+oldc]
		}
	}
	return U, s, V, nil
}

// QR computes a reduced QR factorization of m using Householder reflections.
// For an m x n matrix with rows >= cols it returns Q (rows x cols) with
// orthonormal columns and R (cols x cols) upper triangular such that m = Q*R.
// It returns an error if rows < cols.
func (m *Matrix) QR() (Q *Matrix, R *Matrix, err error) {
	if m.rows < m.cols {
		return nil, nil, errors.New("tensornetwork: QR requires rows >= cols")
	}
	rows, cols := m.rows, m.cols
	r := m.Clone()
	q := IdentityMatrix(rows)
	for k := 0; k < cols; k++ {
		// Householder vector for column k below the diagonal.
		var normx float64
		for i := k; i < rows; i++ {
			normx += r.data[i*cols+k] * r.data[i*cols+k]
		}
		normx = math.Sqrt(normx)
		if normx == 0 {
			continue
		}
		alpha := -math.Copysign(normx, r.data[k*cols+k])
		v := make([]float64, rows)
		v[k] = r.data[k*cols+k] - alpha
		for i := k + 1; i < rows; i++ {
			v[i] = r.data[i*cols+k]
		}
		var vnorm2 float64
		for i := k; i < rows; i++ {
			vnorm2 += v[i] * v[i]
		}
		if vnorm2 < 1e-300 {
			continue
		}
		// Apply H = I - 2 v vᵀ / (vᵀv) to R (columns k..cols).
		for j := k; j < cols; j++ {
			var dot float64
			for i := k; i < rows; i++ {
				dot += v[i] * r.data[i*cols+j]
			}
			f := 2 * dot / vnorm2
			for i := k; i < rows; i++ {
				r.data[i*cols+j] -= f * v[i]
			}
		}
		// Accumulate Q = Q * H (H symmetric).
		for i := 0; i < rows; i++ {
			var dot float64
			for jj := k; jj < rows; jj++ {
				dot += q.data[i*rows+jj] * v[jj]
			}
			f := 2 * dot / vnorm2
			for jj := k; jj < rows; jj++ {
				q.data[i*rows+jj] -= f * v[jj]
			}
		}
	}
	Q = NewMatrix(rows, cols)
	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			Q.data[i*cols+j] = q.data[i*rows+j]
		}
	}
	R = NewMatrix(cols, cols)
	for i := 0; i < cols; i++ {
		for j := i; j < cols; j++ {
			R.data[i*cols+j] = r.data[i*cols+j]
		}
	}
	return Q, R, nil
}

// Rank returns the numerical rank of m, the number of singular values that
// exceed tol times the largest singular value. If tol is not positive a default
// relative tolerance is used.
func (m *Matrix) Rank(tol float64) int {
	_, s, _, err := m.SVD()
	if err != nil || len(s) == 0 {
		return 0
	}
	if tol <= 0 {
		tol = 1e-12
	}
	thr := s[0] * tol
	r := 0
	for _, sv := range s {
		if sv > thr {
			r++
		}
	}
	return r
}

// SpectralNorm returns the largest singular value of m.
func (m *Matrix) SpectralNorm() float64 {
	_, s, _, err := m.SVD()
	if err != nil || len(s) == 0 {
		return 0
	}
	return s[0]
}

// Cond returns the 2-norm condition number of m, the ratio of its largest to
// its smallest singular value. It returns +Inf if the smallest singular value
// is zero.
func (m *Matrix) Cond() float64 {
	_, s, _, err := m.SVD()
	if err != nil || len(s) == 0 {
		return math.Inf(1)
	}
	last := s[len(s)-1]
	if last <= 0 {
		return math.Inf(1)
	}
	return s[0] / last
}

// Pinv returns the Moore-Penrose pseudoinverse of m computed from its SVD.
// Singular values below tol times the largest are treated as zero; if tol is
// not positive a default is used.
func (m *Matrix) Pinv(tol float64) *Matrix {
	U, s, V, err := m.SVD()
	if err != nil {
		return NewMatrix(m.cols, m.rows)
	}
	if tol <= 0 {
		tol = 1e-12
	}
	var smax float64
	if len(s) > 0 {
		smax = s[0]
	}
	thr := smax * tol
	// pinv = V * diag(1/s) * Uᵀ.
	k := len(s)
	sv := NewMatrix(V.rows, k)
	for j := 0; j < k; j++ {
		inv := 0.0
		if s[j] > thr {
			inv = 1 / s[j]
		}
		for i := 0; i < V.rows; i++ {
			sv.data[i*k+j] = V.data[i*k+j] * inv
		}
	}
	res, _ := sv.Mul(U.Transpose())
	return res
}

// Solve returns the least-squares solution x minimizing ‖m·x − b‖ using the
// pseudoinverse. It returns an error if len(b) != m.Rows().
func (m *Matrix) Solve(b []float64) ([]float64, error) {
	if len(b) != m.rows {
		return nil, fmt.Errorf("tensornetwork: length %d != rows %d", len(b), m.rows)
	}
	return m.Pinv(0).MulVec(b)
}

// GramColumns returns the Gram matrix mᵀ·m of the columns of m.
func (m *Matrix) GramColumns() *Matrix {
	g, _ := m.Transpose().Mul(m)
	return g
}
