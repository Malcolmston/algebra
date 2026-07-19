package fem

import "math"

// SparseMatrix is a square sparse matrix in dictionary-of-keys form. It is the
// natural target for finite element assembly, where element contributions are
// scattered into a small number of global entries.
type SparseMatrix struct {
	n    int
	rows []map[int]float64
}

// NewSparseMatrix returns an n×n zero sparse matrix.
func NewSparseMatrix(n int) *SparseMatrix {
	if n < 0 {
		panic("fem: negative dimension in NewSparseMatrix")
	}
	rows := make([]map[int]float64, n)
	for i := range rows {
		rows[i] = make(map[int]float64)
	}
	return &SparseMatrix{n: n, rows: rows}
}

// Dim returns the dimension n of the n×n matrix.
func (s *SparseMatrix) Dim() int { return s.n }

// At returns the entry at row i, column j (zero if absent).
func (s *SparseMatrix) At(i, j int) float64 { return s.rows[i][j] }

// Set assigns v to entry (i,j), removing the entry if v is exactly zero.
func (s *SparseMatrix) Set(i, j int, v float64) {
	if v == 0 {
		delete(s.rows[i], j)
		return
	}
	s.rows[i][j] = v
}

// AddEntry scatters v into entry (i,j), the fundamental assembly operation.
func (s *SparseMatrix) AddEntry(i, j int, v float64) {
	if v == 0 {
		return
	}
	s.rows[i][j] += v
}

// NNZ returns the number of stored (structurally non-zero) entries.
func (s *SparseMatrix) NNZ() int {
	c := 0
	for _, r := range s.rows {
		c += len(r)
	}
	return c
}

// MulVec returns the sparse matrix–vector product s*x.
func (s *SparseMatrix) MulVec(x Vector) Vector {
	if len(x) != s.n {
		panic("fem: dimension mismatch in SparseMatrix.MulVec")
	}
	out := make(Vector, s.n)
	for i := 0; i < s.n; i++ {
		var sum float64
		for j, v := range s.rows[i] {
			sum += v * x[j]
		}
		out[i] = sum
	}
	return out
}

// ToDense returns the dense representation of the sparse matrix.
func (s *SparseMatrix) ToDense() *Matrix {
	m := NewMatrix(s.n, s.n)
	for i := 0; i < s.n; i++ {
		for j, v := range s.rows[i] {
			m.Set(i, j, v)
		}
	}
	return m
}

// ScaleRow multiplies every stored entry of row i by a.
func (s *SparseMatrix) ScaleRow(i int, a float64) {
	for j := range s.rows[i] {
		s.rows[i][j] *= a
	}
}

// ClearRow removes every entry of row i.
func (s *SparseMatrix) ClearRow(i int) {
	s.rows[i] = make(map[int]float64)
}

// IsSymmetric reports whether the sparse matrix is symmetric within tol.
func (s *SparseMatrix) IsSymmetric(tol float64) bool {
	for i := 0; i < s.n; i++ {
		for j, v := range s.rows[i] {
			if math.Abs(v-s.At(j, i)) > tol {
				return false
			}
		}
	}
	return true
}

// Diagonal returns the main diagonal of the sparse matrix as a Vector.
func (s *SparseMatrix) Diagonal() Vector {
	d := make(Vector, s.n)
	for i := 0; i < s.n; i++ {
		d[i] = s.rows[i][i]
	}
	return d
}

// ConjugateGradient solves the symmetric positive-definite system A x = b using
// the (unpreconditioned) conjugate gradient method. It returns the solution and
// the number of iterations used, or ErrNotConverged.
func ConjugateGradient(a *SparseMatrix, b Vector, tol float64, maxIter int) (Vector, int, error) {
	n := a.n
	if len(b) != n {
		panic("fem: dimension mismatch in ConjugateGradient")
	}
	if maxIter <= 0 {
		maxIter = 10 * (n + 1)
	}
	x := make(Vector, n)
	r := b.Clone()
	p := r.Clone()
	rsold := r.Dot(r)
	bnorm := b.Norm2()
	if bnorm == 0 {
		bnorm = 1
	}
	if math.Sqrt(rsold)/bnorm <= tol {
		return x, 0, nil
	}
	for k := 1; k <= maxIter; k++ {
		ap := a.MulVec(p)
		pap := p.Dot(ap)
		if pap == 0 {
			return x, k, ErrNotConverged
		}
		alpha := rsold / pap
		for i := range x {
			x[i] += alpha * p[i]
			r[i] -= alpha * ap[i]
		}
		rsnew := r.Dot(r)
		if math.Sqrt(rsnew)/bnorm <= tol {
			return x, k, nil
		}
		beta := rsnew / rsold
		for i := range p {
			p[i] = r[i] + beta*p[i]
		}
		rsold = rsnew
	}
	return x, maxIter, ErrNotConverged
}

// PCGJacobi solves A x = b with Jacobi (diagonal) preconditioned conjugate
// gradients. It is suited to symmetric positive-definite systems.
func PCGJacobi(a *SparseMatrix, b Vector, tol float64, maxIter int) (Vector, int, error) {
	n := a.n
	if len(b) != n {
		panic("fem: dimension mismatch in PCGJacobi")
	}
	if maxIter <= 0 {
		maxIter = 10 * (n + 1)
	}
	diag := a.Diagonal()
	inv := make(Vector, n)
	for i, d := range diag {
		if d == 0 {
			inv[i] = 1
		} else {
			inv[i] = 1 / d
		}
	}
	x := make(Vector, n)
	r := b.Clone()
	z := make(Vector, n)
	for i := range z {
		z[i] = inv[i] * r[i]
	}
	p := z.Clone()
	rz := r.Dot(z)
	bnorm := b.Norm2()
	if bnorm == 0 {
		bnorm = 1
	}
	for k := 1; k <= maxIter; k++ {
		if r.Norm2()/bnorm <= tol {
			return x, k - 1, nil
		}
		ap := a.MulVec(p)
		pap := p.Dot(ap)
		if pap == 0 {
			return x, k, ErrNotConverged
		}
		alpha := rz / pap
		for i := range x {
			x[i] += alpha * p[i]
			r[i] -= alpha * ap[i]
		}
		for i := range z {
			z[i] = inv[i] * r[i]
		}
		rznew := r.Dot(z)
		beta := rznew / rz
		for i := range p {
			p[i] = z[i] + beta*p[i]
		}
		rz = rznew
	}
	if r.Norm2()/bnorm <= tol {
		return x, maxIter, nil
	}
	return x, maxIter, ErrNotConverged
}

// SolveSPD solves a symmetric positive-definite sparse system, preferring the
// dense LU path for small problems and Jacobi-preconditioned CG otherwise.
func SolveSPD(a *SparseMatrix, b Vector) (Vector, error) {
	if a.n <= 64 {
		return SolveDense(a.ToDense(), b)
	}
	x, _, err := PCGJacobi(a, b, 1e-12, 20*(a.n+1))
	if err != nil {
		return SolveDense(a.ToDense(), b)
	}
	return x, nil
}
