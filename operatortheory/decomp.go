package operatortheory

import (
	"math"
	"math/cmplx"
	"sort"
)

// Determinant returns the determinant of a square matrix. It returns
// ErrNotSquare for a non-square matrix.
func (m *Matrix) Determinant() (complex128, error) {
	if !m.IsSquare() {
		return 0, ErrNotSquare
	}
	if m.rows == 0 {
		return 1, nil
	}
	return luDecompose(m).det(), nil
}

// Inverse returns the inverse of a square matrix. It returns ErrNotSquare for a
// non-square matrix and ErrSingular when the matrix is (numerically) singular.
func (m *Matrix) Inverse() (*Matrix, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	n := m.rows
	f := luDecompose(m)
	if f.singular {
		return nil, ErrSingular
	}
	inv := NewMatrix(n, n)
	for j := 0; j < n; j++ {
		e := make(Vector, n)
		e[j] = 1
		x := f.solve(e)
		for i := 0; i < n; i++ {
			inv.data[i*n+j] = x[i]
		}
	}
	return inv, nil
}

// Solve returns the solution X of the linear system m*X = b for a matrix
// right-hand side. It returns ErrNotSquare, ErrDimensionMismatch or ErrSingular
// as appropriate.
func (m *Matrix) Solve(b *Matrix) (*Matrix, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	if b.rows != m.rows {
		return nil, ErrDimensionMismatch
	}
	f := luDecompose(m)
	if f.singular {
		return nil, ErrSingular
	}
	x := NewMatrix(m.rows, b.cols)
	for j := 0; j < b.cols; j++ {
		col := f.solve(b.Col(j))
		for i := 0; i < m.rows; i++ {
			x.data[i*x.cols+j] = col[i]
		}
	}
	return x, nil
}

// SolveVec returns the solution x of the linear system m*x = b for a vector
// right-hand side.
func (m *Matrix) SolveVec(b Vector) (Vector, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	if len(b) != m.rows {
		return nil, ErrDimensionMismatch
	}
	f := luDecompose(m)
	if f.singular {
		return nil, ErrSingular
	}
	return f.solve(b), nil
}

// QR returns a reduced QR factorisation m = Q*R computed with complex
// Householder reflections. Q has orthonormal columns and R is upper triangular.
func (m *Matrix) QR() (q, r *Matrix) {
	return m.qrRaw()
}

// Hessenberg returns the upper-Hessenberg reduction H together with the unitary
// Q such that m = Q*H*Q^H. It returns ErrNotSquare for a non-square matrix.
func (m *Matrix) Hessenberg() (h, q *Matrix, err error) {
	if !m.IsSquare() {
		return nil, nil, ErrNotSquare
	}
	h, q = m.hessenbergRaw()
	return h, q, nil
}

// HermitianEigen returns the real eigenvalues (in ascending order) and the
// orthonormal eigenvectors (as the columns of the returned matrix) of a
// Hermitian matrix. The matrix is not required to be exactly Hermitian; its
// Hermitian part is used, so callers should ensure the input truly is Hermitian
// for meaningful results. It returns ErrNotSquare for a non-square matrix.
func (m *Matrix) HermitianEigen() (values []float64, vectors *Matrix, err error) {
	if !m.IsSquare() {
		return nil, nil, ErrNotSquare
	}
	h := m.HermitianPart()
	vals, vecs := hermitianEigenRaw(h)
	return vals, vecs, nil
}

// Eigenvalues returns all eigenvalues of a square matrix, counted with
// multiplicity, computed with the shifted QR algorithm. The order is
// unspecified. It returns ErrNotSquare for a non-square matrix.
func (m *Matrix) Eigenvalues() ([]complex128, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	if m.rows == 0 {
		return nil, nil
	}
	return eigenvaluesQR(m), nil
}

// Eigenpair bundles an eigenvalue with a corresponding unit eigenvector.
type Eigenpair struct {
	// Value is the eigenvalue.
	Value complex128
	// Vector is a unit eigenvector for Value.
	Vector Vector
}

// Eigen returns eigenvalue/eigenvector pairs for a square matrix. Eigenvalues
// are found with the QR algorithm and eigenvectors by inverse iteration. The
// eigenvectors are reliable for matrices with well-separated eigenvalues; for
// Hermitian matrices prefer HermitianEigen. It returns ErrNotSquare for a
// non-square matrix.
func (m *Matrix) Eigen() ([]Eigenpair, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	vals := eigenvaluesQR(m)
	pairs := make([]Eigenpair, len(vals))
	for i, lam := range vals {
		v := m.eigenvector(lam)
		pairs[i] = Eigenpair{Value: lam, Vector: v}
	}
	return pairs, nil
}

// eigenvector computes a unit eigenvector for the (approximate) eigenvalue lam
// by inverse iteration on m - (lam + eps) I.
func (m *Matrix) eigenvector(lam complex128) Vector {
	n := m.rows
	scale := m.FrobeniusNorm()
	if scale == 0 {
		scale = 1
	}
	eps := complex(1e-10*scale, 1e-10*scale)
	shifted := m.Clone()
	for i := 0; i < n; i++ {
		shifted.data[i*n+i] -= lam + eps
	}
	f := luDecompose(shifted)
	// Start from a fixed spread-out vector.
	x := make(Vector, n)
	for i := range x {
		x[i] = complex(1+0.1*float64(i), 0.1*float64(i))
	}
	x, _ = x.Normalize()
	if f.singular {
		return x
	}
	for iter := 0; iter < 8; iter++ {
		y := f.solve(x)
		ny := y.Norm()
		if ny == 0 {
			break
		}
		x = y.Scale(complex(1/ny, 0))
	}
	// Fix a canonical global phase: make the largest-modulus entry real
	// positive.
	best := 0
	for i := 1; i < n; i++ {
		if cmplx.Abs(x[i]) > cmplx.Abs(x[best]) {
			best = i
		}
	}
	if a := cmplx.Abs(x[best]); a > 0 {
		phase := x[best] / complex(a, 0)
		x = x.Scale(cmplx.Conj(phase))
	}
	return x
}

// SVD returns a reduced singular value decomposition m = U * diag(s) * V^H,
// where U has orthonormal columns (m-by-k), s holds the k = min(rows,cols)
// singular values in descending order, and V has orthonormal columns
// (cols-by-k).
func (m *Matrix) SVD() (u *Matrix, s []float64, v *Matrix) {
	rows, cols := m.rows, m.cols
	k := cols
	if rows < k {
		k = rows
	}
	if rows >= cols {
		// Right vectors from A^H A.
		ah := m.Adjoint()
		g, _ := ah.Mul(m) // cols-by-cols Hermitian
		vals, vecs := hermitianEigenRaw(g)
		order := descendingOrder(vals)
		s = make([]float64, k)
		v = NewMatrix(cols, k)
		uCols := make([]Vector, k)
		for c := 0; c < k; c++ {
			idx := order[c]
			sig := math.Sqrt(math.Max(0, vals[idx]))
			s[c] = sig
			vc := colOf(vecs, idx)
			for i := 0; i < cols; i++ {
				v.data[i*k+c] = vc[i]
			}
			av, _ := m.MulVec(vc)
			if sig > 1e-12*(1+s[0]) {
				uCols[c] = av.Scale(complex(1/sig, 0))
			} else {
				uCols[c] = nil
			}
		}
		fillMissing(uCols, rows)
		u = columnsToMatrix(uCols, rows)
		return u, s, v
	}
	// rows < cols: left vectors from A A^H.
	ah := m.Adjoint()
	g, _ := m.Mul(ah) // rows-by-rows Hermitian
	vals, vecs := hermitianEigenRaw(g)
	order := descendingOrder(vals)
	s = make([]float64, k)
	u = NewMatrix(rows, k)
	vCols := make([]Vector, k)
	for c := 0; c < k; c++ {
		idx := order[c]
		sig := math.Sqrt(math.Max(0, vals[idx]))
		s[c] = sig
		uc := colOf(vecs, idx)
		for i := 0; i < rows; i++ {
			u.data[i*k+c] = uc[i]
		}
		av, _ := ah.MulVec(uc)
		if sig > 1e-12*(1+s[0]) {
			vCols[c] = av.Scale(complex(1/sig, 0))
		} else {
			vCols[c] = nil
		}
	}
	fillMissing(vCols, cols)
	v = columnsToMatrix(vCols, cols)
	return u, s, v
}

// SingularValues returns the singular values of the matrix in descending order.
func (m *Matrix) SingularValues() []float64 {
	_, s, _ := m.SVD()
	return s
}

// Rank returns the numerical rank of the matrix, the number of singular values
// exceeding tol * (largest singular value). If tol <= 0 a default relative
// tolerance is used.
func (m *Matrix) Rank(tol float64) int {
	s := m.SingularValues()
	if len(s) == 0 {
		return 0
	}
	if tol <= 0 {
		tol = 1e-12
	}
	thresh := tol * s[0] * float64(maxInt(m.rows, m.cols))
	r := 0
	for _, sv := range s {
		if sv > thresh {
			r++
		}
	}
	return r
}

// ConditionNumber returns the 2-norm condition number, the ratio of the largest
// to the smallest singular value. It returns +Inf for a singular matrix.
func (m *Matrix) ConditionNumber() float64 {
	s := m.SingularValues()
	if len(s) == 0 {
		return 0
	}
	smin := s[len(s)-1]
	if smin == 0 {
		return math.Inf(1)
	}
	return s[0] / smin
}

// PolarDecomposition returns the polar factors of a square matrix, U unitary
// (or a partial isometry when singular) and P Hermitian positive semidefinite,
// such that m = U*P. It returns ErrNotSquare for a non-square matrix.
func (m *Matrix) PolarDecomposition() (unitary, positive *Matrix, err error) {
	if !m.IsSquare() {
		return nil, nil, ErrNotSquare
	}
	u, s, v := m.SVD()
	vh := v.Adjoint()
	// P = V diag(s) V^H.
	sv := v.Clone()
	n := m.rows
	for c := 0; c < n; c++ {
		for i := 0; i < n; i++ {
			sv.data[i*n+c] *= complex(s[c], 0)
		}
	}
	positive, _ = sv.Mul(vh)
	// W = U V^H.
	unitary, _ = u.Mul(vh)
	return unitary, positive, nil
}

// descendingOrder returns indices ordering vals from largest to smallest.
func descendingOrder(vals []float64) []int {
	idx := make([]int, len(vals))
	for i := range idx {
		idx[i] = i
	}
	sort.Slice(idx, func(a, b int) bool { return vals[idx[a]] > vals[idx[b]] })
	return idx
}

// colOf extracts column j of the square matrix as a Vector.
func colOf(m *Matrix, j int) Vector {
	v := make(Vector, m.rows)
	for i := 0; i < m.rows; i++ {
		v[i] = m.data[i*m.cols+j]
	}
	return v
}

// fillMissing replaces nil columns by orthonormal completions (via
// Gram-Schmidt against the present columns and the standard basis).
func fillMissing(cols []Vector, dim int) {
	var present []Vector
	for _, c := range cols {
		if c != nil {
			present = append(present, c)
		}
	}
	next := 0
	for i := range cols {
		if cols[i] != nil {
			continue
		}
		for next < dim {
			cand := BasisVector(dim, next)
			next++
			for _, q := range present {
				cand = cand.Sub(q.Scale(q.Dot(cand)))
			}
			nrm := cand.Norm()
			if nrm > 1e-9 {
				cand = cand.Scale(complex(1/nrm, 0))
				cols[i] = cand
				present = append(present, cand)
				break
			}
		}
		if cols[i] == nil {
			cols[i] = NewVector(dim)
		}
	}
}

// columnsToMatrix assembles the given columns into a dim-by-len matrix.
func columnsToMatrix(cols []Vector, dim int) *Matrix {
	k := len(cols)
	m := NewMatrix(dim, k)
	for c := 0; c < k; c++ {
		for i := 0; i < dim; i++ {
			m.data[i*k+c] = cols[c][i]
		}
	}
	return m
}

func maxInt(a, b int) int {
	if a > b {
		return a
	}
	return b
}
