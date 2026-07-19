package lattices

import "math/big"

// Dual returns the dual (reciprocal) lattice basis d_0, ..., d_{n-1}, the unique
// vectors in the span of the basis satisfying <d_i, b_j> = 1 if i == j and 0
// otherwise. It is computed as G^{-1} B, where G is the Gram matrix and B the
// basis rows. It returns ErrEmpty for an empty basis and ErrNotFullRank when the
// Gram matrix is singular.
func (b Basis) Dual() (Basis, error) {
	if len(b) == 0 {
		return nil, ErrEmpty
	}
	g := b.Gram()
	gi, err := g.Inverse()
	if err != nil {
		return nil, ErrNotFullRank
	}
	n := len(b)
	dual := make(Basis, n)
	for i := 0; i < n; i++ {
		d := ZeroVec(b.Dim())
		for k := 0; k < n; k++ {
			d = d.AddScaled(b[k], gi.At(i, k))
		}
		dual[i] = d
	}
	return dual, nil
}

// DualRat returns the exact rational dual basis, computed with exact rational
// arithmetic. It returns ErrEmpty or ErrNotFullRank as appropriate.
func (b Basis) DualRat() ([]RatVec, error) {
	if len(b) == 0 {
		return nil, ErrEmpty
	}
	g := b.GramRat()
	gi, err := g.Inverse()
	if err != nil {
		return nil, ErrNotFullRank
	}
	n := len(b)
	rows := make([]RatVec, n)
	for i := range b {
		rows[i] = RatVecFromFloats(b[i]...)
	}
	dual := make([]RatVec, n)
	for i := 0; i < n; i++ {
		d := NewRatVec(b.Dim())
		for k := 0; k < n; k++ {
			d = d.AddScaled(rows[k], gi.data[i][k])
		}
		dual[i] = d
	}
	return dual, nil
}

// DualDeterminant returns the covolume of the dual lattice, which equals the
// reciprocal of the covolume of the primal lattice. It returns 0 for a
// degenerate lattice.
func (b Basis) DualDeterminant() float64 {
	d := b.Determinant()
	if d == 0 {
		return 0
	}
	return 1 / d
}

// Biorthogonality returns the matrix M[i][j] = <d_i, b_j> where d is the dual
// basis; it should equal the identity matrix up to numerical error. It returns
// ErrEmpty or ErrNotFullRank as appropriate.
func (b Basis) Biorthogonality() (Matrix, error) {
	dual, err := b.Dual()
	if err != nil {
		return Matrix{}, err
	}
	n := len(b)
	m := ZeroMatrix(n, n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			m.data[i][j] = dual[i].Dot(b[j])
		}
	}
	return m, nil
}

// IsDualTo reports whether other is (numerically, within tol) the dual basis of
// b, i.e. whether <other_i, b_j> is the identity matrix within tol.
func (b Basis) IsDualTo(other Basis, tol float64) bool {
	if len(other) != len(b) || other.Dim() != b.Dim() {
		return false
	}
	n := len(b)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			want := 0.0
			if i == j {
				want = 1
			}
			if v := other[i].Dot(b[j]); v > want+tol || v < want-tol {
				return false
			}
		}
	}
	return true
}

// GramInverseRat returns the exact inverse of the rational Gram matrix, the
// Gram matrix of the dual basis. It returns ErrEmpty or ErrNotFullRank as
// appropriate.
func (b Basis) GramInverseRat() (RatMatrix, error) {
	if len(b) == 0 {
		return RatMatrix{}, ErrEmpty
	}
	gi, err := b.GramRat().Inverse()
	if err != nil {
		return RatMatrix{}, ErrNotFullRank
	}
	return gi, nil
}

// DualGramDeterminantRat returns the exact determinant of the dual Gram matrix,
// which equals the reciprocal of the primal Gram determinant.
func (b Basis) DualGramDeterminantRat() (*big.Rat, error) {
	d := b.GramDeterminantRat()
	if d.Sign() == 0 {
		return nil, ErrNotFullRank
	}
	return new(big.Rat).Inv(d), nil
}
