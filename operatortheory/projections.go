package operatortheory

import (
	"math"
	"sort"
)

// OrthogonalProjector returns the orthogonal projection operator onto the
// subspace spanned by the given vectors of common length n. Linearly dependent
// vectors are handled correctly. It returns ErrEmpty when no vectors are given
// and ErrDimensionMismatch on unequal lengths.
func OrthogonalProjector(vectors []Vector) (*Matrix, error) {
	if len(vectors) == 0 {
		return nil, ErrEmpty
	}
	n := len(vectors[0])
	for _, v := range vectors {
		if len(v) != n {
			return nil, ErrDimensionMismatch
		}
	}
	basis := GramSchmidt(vectors, 1e-12)
	p := NewMatrix(n, n)
	for _, q := range basis {
		add, _ := p.Add(q.OuterProduct(q))
		p = add
	}
	return p, nil
}

// RangeProjector returns the orthogonal projection onto the column space (range)
// of m.
func (m *Matrix) RangeProjector() *Matrix {
	var cols []Vector
	for j := 0; j < m.cols; j++ {
		cols = append(cols, m.Col(j))
	}
	p, err := OrthogonalProjector(cols)
	if err != nil {
		return NewMatrix(m.rows, m.rows)
	}
	return p
}

// KernelBasis returns an orthonormal basis for the null space (kernel) of a
// square matrix, obtained from the right singular vectors with negligible
// singular values. tol sets the relative threshold below which a singular value
// is treated as zero.
func (m *Matrix) KernelBasis(tol float64) []Vector {
	if tol <= 0 {
		tol = 1e-9
	}
	_, s, v := m.SVD()
	if len(s) == 0 {
		return nil
	}
	thresh := tol * s[0] * float64(maxInt(m.rows, m.cols))
	var basis []Vector
	for j := 0; j < len(s); j++ {
		if s[j] <= thresh {
			basis = append(basis, colOf(v, j))
		}
	}
	// Columns beyond len(s) (when cols > rows) are also in the kernel; SVD is
	// reduced, so account for them via the range projector complement.
	return basis
}

// SpectralComponent bundles a (real) eigenvalue of a Hermitian operator with
// the orthogonal projector onto its eigenspace.
type SpectralComponent struct {
	// Eigenvalue is the (real) eigenvalue.
	Eigenvalue float64
	// Multiplicity is the dimension of the eigenspace.
	Multiplicity int
	// Projector is the orthogonal projection onto the eigenspace.
	Projector *Matrix
}

// SpectralDecomposition returns the spectral resolution of a Hermitian matrix:
// one SpectralComponent per distinct eigenvalue (eigenvalues within tol are
// merged), sorted in ascending order. The projectors sum to the identity and
// m = sum lambda_k P_k. It returns ErrNotSquare for a non-square matrix.
func (m *Matrix) SpectralDecomposition(tol float64) ([]SpectralComponent, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	tol = orDefault(tol)
	vals, vecs := hermitianEigenRaw(m.HermitianPart())
	n := m.rows
	type group struct {
		val  float64
		cols []Vector
	}
	var groups []group
	// vals are ascending from hermitianEigenRaw.
	order := make([]int, len(vals))
	for i := range order {
		order[i] = i
	}
	sort.Slice(order, func(a, b int) bool { return vals[order[a]] < vals[order[b]] })
	for _, idx := range order {
		v := vals[idx]
		col := colOf(vecs, idx)
		if len(groups) > 0 && math.Abs(groups[len(groups)-1].val-v) <= tol {
			g := &groups[len(groups)-1]
			g.cols = append(g.cols, col)
		} else {
			groups = append(groups, group{val: v, cols: []Vector{col}})
		}
	}
	comps := make([]SpectralComponent, len(groups))
	for i, g := range groups {
		p := NewMatrix(n, n)
		for _, q := range g.cols {
			add, _ := p.Add(q.OuterProduct(q))
			p = add
		}
		comps[i] = SpectralComponent{
			Eigenvalue:   g.val,
			Multiplicity: len(g.cols),
			Projector:    p,
		}
	}
	return comps, nil
}

// DiagonalizeHermitian returns a unitary U and a real diagonal matrix D such
// that m = U*D*U^H, for a Hermitian matrix. It returns ErrNotSquare for a
// non-square matrix.
func (m *Matrix) DiagonalizeHermitian() (u, d *Matrix, err error) {
	if !m.IsSquare() {
		return nil, nil, ErrNotSquare
	}
	vals, vecs := hermitianEigenRaw(m.HermitianPart())
	dd := NewMatrix(m.rows, m.rows)
	for i, v := range vals {
		dd.data[i*m.rows+i] = complex(v, 0)
	}
	return vecs, dd, nil
}

// SpectralProjector returns the orthogonal projector onto the eigenspace of a
// Hermitian matrix associated with the eigenvalue closest to target, provided
// it lies within tol; otherwise it returns the zero matrix. It is the operator
// obtained from functional calculus with the indicator of {target}.
func (m *Matrix) SpectralProjector(target, tol float64) *Matrix {
	comps, err := m.SpectralDecomposition(tol)
	if err != nil {
		return NewMatrix(m.rows, m.rows)
	}
	best := -1
	bestDist := math.Inf(1)
	for i, c := range comps {
		d := math.Abs(c.Eigenvalue - target)
		if d < bestDist {
			bestDist = d
			best = i
		}
	}
	if best < 0 || bestDist > orDefault(tol) {
		return NewMatrix(m.rows, m.rows)
	}
	return comps[best].Projector
}

// Reflection returns the Householder reflection I - 2 P, where P is the
// orthogonal projector onto the line spanned by the nonzero unit-normalisable
// vector v. The result is a unitary involution.
func Reflection(v Vector) (*Matrix, error) {
	u, n := v.Normalize()
	if n == 0 {
		return nil, ErrInvalidArgument
	}
	dim := len(v)
	p := u.OuterProduct(u)
	r := Identity(dim)
	sub, _ := r.Sub(p.Scale(2))
	return sub, nil
}
