package lattices

import (
	"math"
	"math/big"
)

// Basis is a lattice basis: an ordered list of row vectors that generate the
// lattice by integer linear combination. The vectors should be linearly
// independent; the number of vectors is the lattice rank and the length of each
// vector is the ambient dimension.
type Basis []Vec

// NewBasis returns a Basis containing copies of the given row vectors. It
// panics if the rows have differing lengths.
func NewBasis(rows ...Vec) Basis {
	b := make(Basis, len(rows))
	dim := -1
	for i, r := range rows {
		if dim == -1 {
			dim = len(r)
		} else if len(r) != dim {
			panic("lattices: basis rows of unequal length")
		}
		b[i] = r.Clone()
	}
	return b
}

// BasisFromRows builds a Basis from a slice of float64 row slices.
func BasisFromRows(rows [][]float64) Basis {
	b := make(Basis, len(rows))
	for i, r := range rows {
		b[i] = NewVec(r...)
	}
	return b
}

// BasisFromMatrix builds a Basis whose vectors are the rows of m.
func BasisFromMatrix(m Matrix) Basis {
	b := make(Basis, m.rows)
	for i := 0; i < m.rows; i++ {
		b[i] = m.Row(i)
	}
	return b
}

// IdentityBasis returns the standard basis of Z^n (the n rows of the identity).
func IdentityBasis(n int) Basis {
	b := make(Basis, n)
	for i := 0; i < n; i++ {
		b[i] = UnitVec(n, i)
	}
	return b
}

// Rank returns the number of vectors in the basis.
func (b Basis) Rank() int { return len(b) }

// Dim returns the ambient dimension (length of each vector), or 0 for an empty
// basis.
func (b Basis) Dim() int {
	if len(b) == 0 {
		return 0
	}
	return len(b[0])
}

// Clone returns an independent deep copy of the basis.
func (b Basis) Clone() Basis {
	c := make(Basis, len(b))
	for i := range b {
		c[i] = b[i].Clone()
	}
	return c
}

// Row returns the i-th basis vector (the stored slice, not a copy).
func (b Basis) Row(i int) Vec { return b[i] }

// Swap exchanges basis vectors i and j in place.
func (b Basis) Swap(i, j int) { b[i], b[j] = b[j], b[i] }

// Matrix returns the basis as a Matrix whose rows are the basis vectors.
func (b Basis) Matrix() Matrix {
	rows := make([][]float64, len(b))
	for i := range b {
		rows[i] = b[i]
	}
	return NewMatrix(rows)
}

// Gram returns the Gram matrix G with G[i][j] = <b_i, b_j>.
func (b Basis) Gram() Matrix {
	n := len(b)
	g := ZeroMatrix(n, n)
	for i := 0; i < n; i++ {
		for j := i; j < n; j++ {
			d := b[i].Dot(b[j])
			g.data[i][j] = d
			g.data[j][i] = d
		}
	}
	return g
}

// GramRat returns the exact rational Gram matrix of the basis.
func (b Basis) GramRat() RatMatrix {
	n := len(b)
	rb := make([]RatVec, n)
	for i := range b {
		rb[i] = RatVecFromFloats(b[i]...)
	}
	g := ZeroRatMatrix(n, n)
	for i := 0; i < n; i++ {
		for j := i; j < n; j++ {
			d := rb[i].Dot(rb[j])
			g.data[i][j].Set(d)
			g.data[j][i].Set(d)
		}
	}
	return g
}

// GramDeterminant returns det(G), the determinant of the Gram matrix. For a
// full-rank basis this equals the square of the lattice covolume.
func (b Basis) GramDeterminant() float64 {
	d, _ := b.Gram().Det()
	return d
}

// GramDeterminantRat returns the exact determinant of the Gram matrix.
func (b Basis) GramDeterminantRat() *big.Rat {
	d, err := b.GramRat().Det()
	if err != nil || d == nil {
		return new(big.Rat)
	}
	return d
}

// Determinant returns the lattice determinant (covolume) sqrt(det G). It is the
// volume of a fundamental parallelepiped and is invariant under change of
// basis.
func (b Basis) Determinant() float64 {
	d := b.GramDeterminant()
	if d < 0 {
		d = 0
	}
	return math.Sqrt(d)
}

// Volume is an alias for Determinant: the covolume of the lattice.
func (b Basis) Volume() float64 { return b.Determinant() }

// IsFullRank reports whether the basis vectors are linearly independent, i.e.
// the Gram determinant is nonzero within the tolerance tol.
func (b Basis) IsFullRank(tol float64) bool {
	return math.Abs(b.GramDeterminant()) > tol
}

// IsSquare reports whether the number of vectors equals the ambient dimension.
func (b Basis) IsSquare() bool { return b.Rank() == b.Dim() }

// Point returns the lattice point sum_i coeffs[i]*b_i for integer coefficients.
// It panics if len(coeffs) != Rank.
func (b Basis) Point(coeffs []int64) Vec {
	if len(coeffs) != len(b) {
		panic(ErrDimMismatch)
	}
	r := ZeroVec(b.Dim())
	for i, c := range coeffs {
		r = r.AddScaled(b[i], float64(c))
	}
	return r
}

// PointFloat returns sum_i coeffs[i]*b_i for real coefficients (a point in the
// real span of the lattice). It panics if len(coeffs) != Rank.
func (b Basis) PointFloat(coeffs []float64) Vec {
	if len(coeffs) != len(b) {
		panic(ErrDimMismatch)
	}
	r := ZeroVec(b.Dim())
	for i, c := range coeffs {
		r = r.AddScaled(b[i], c)
	}
	return r
}

// Coordinates returns the coefficient vector x such that b represents v, that
// is sum_i x[i]*b_i = v, by solving the (possibly overdetermined) system in a
// least-squares sense via the normal equations G x = B v. It returns
// ErrNotFullRank when the Gram matrix is singular.
func (b Basis) Coordinates(v Vec) (Vec, error) {
	if len(b) == 0 {
		return nil, ErrEmpty
	}
	if len(v) != b.Dim() {
		return nil, ErrDimMismatch
	}
	g := b.Gram()
	rhs := make(Vec, len(b))
	for i := range b {
		rhs[i] = b[i].Dot(v)
	}
	x, err := g.Solve(rhs)
	if err != nil {
		return nil, ErrNotFullRank
	}
	return x, nil
}

// OrthogonalityDefect returns the ratio prod_i |b_i| / covolume, a measure of
// how far the basis is from orthogonal. It equals 1 exactly for an orthogonal
// basis and is at least 1 otherwise (Hadamard's inequality).
func (b Basis) OrthogonalityDefect() float64 {
	det := b.Determinant()
	if det == 0 {
		return math.Inf(1)
	}
	prod := 1.0
	for _, v := range b {
		prod *= v.Norm()
	}
	return prod / det
}

// HadamardRatio returns covolume / prod_i |b_i|, a value in (0, 1]. Values near
// 1 indicate a nearly orthogonal ("good") basis.
func (b Basis) HadamardRatio() float64 {
	prod := 1.0
	for _, v := range b {
		prod *= v.Norm()
	}
	if prod == 0 {
		return 0
	}
	r := b.Determinant() / prod
	return math.Pow(r, 1.0/float64(len(b)))
}

// ShortestBasisVector returns the index and Euclidean norm of the shortest
// nonzero vector currently in the basis (not the lattice minimum). It returns
// -1 and 0 for an empty basis.
func (b Basis) ShortestBasisVector() (int, float64) {
	best := -1
	bestNorm := math.Inf(1)
	for i, v := range b {
		if v.IsZero() {
			continue
		}
		if n := v.Norm(); n < bestNorm {
			bestNorm = n
			best = i
		}
	}
	if best == -1 {
		return -1, 0
	}
	return best, bestNorm
}

// MaxNorm returns the largest Euclidean norm among the basis vectors.
func (b Basis) MaxNorm() float64 {
	var m float64
	for _, v := range b {
		if n := v.Norm(); n > m {
			m = n
		}
	}
	return m
}
