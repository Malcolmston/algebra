package lattices

import "math"

// GramSchmidt holds the result of Gram-Schmidt orthogonalization of a basis:
// the orthogonal vectors Star (the b*_i), the strictly lower-triangular matrix
// of Gram-Schmidt coefficients Mu (Mu[i][j] = <b_i, b*_j>/<b*_j, b*_j> for
// j < i), and the squared norms Norm2 of the orthogonal vectors.
type GramSchmidt struct {
	Star  []Vec
	Mu    [][]float64
	Norm2 []float64
}

// Orthogonalize computes the Gram-Schmidt orthogonalization of the basis
// without normalizing the resulting vectors. The i-th orthogonal vector is
// b_i minus its projection onto the span of the earlier orthogonal vectors.
func (b Basis) Orthogonalize() GramSchmidt {
	n := len(b)
	star := make([]Vec, n)
	mu := make([][]float64, n)
	norm2 := make([]float64, n)
	for i := 0; i < n; i++ {
		mu[i] = make([]float64, n)
		star[i] = b[i].Clone()
		for j := 0; j < i; j++ {
			if norm2[j] == 0 {
				mu[i][j] = 0
				continue
			}
			m := b[i].Dot(star[j]) / norm2[j]
			mu[i][j] = m
			star[i] = star[i].AddScaled(star[j], -m)
		}
		mu[i][i] = 1
		norm2[i] = star[i].Norm2()
	}
	return GramSchmidt{Star: star, Mu: mu, Norm2: norm2}
}

// StarBasis returns the orthogonalized vectors b*_i as a Basis.
func (b Basis) StarBasis() Basis {
	gs := b.Orthogonalize()
	return Basis(gs.Star)
}

// MuMatrix returns the strictly lower-triangular matrix of Gram-Schmidt
// coefficients (with ones on the diagonal) as a Matrix.
func (b Basis) MuMatrix() Matrix {
	gs := b.Orthogonalize()
	return NewMatrix(gs.Mu)
}

// GramSchmidtNorms returns the squared norms of the Gram-Schmidt vectors b*_i.
func (b Basis) GramSchmidtNorms() []float64 {
	gs := b.Orthogonalize()
	out := make([]float64, len(gs.Norm2))
	copy(out, gs.Norm2)
	return out
}

// OrthogonalNormal returns the Gram-Schmidt vectors normalized to unit length,
// forming an orthonormal basis of the span. Zero vectors are left unchanged.
func (b Basis) OrthogonalNormal() Basis {
	gs := b.Orthogonalize()
	out := make(Basis, len(gs.Star))
	for i, v := range gs.Star {
		out[i] = v.Normalize()
	}
	return out
}

// IsSizeReduced reports whether every Gram-Schmidt coefficient Mu[i][j] with
// j < i has absolute value at most 1/2 + tol, the size-reduced condition.
func (b Basis) IsSizeReduced(tol float64) bool {
	gs := b.Orthogonalize()
	for i := range gs.Mu {
		for j := 0; j < i; j++ {
			if math.Abs(gs.Mu[i][j]) > 0.5+tol {
				return false
			}
		}
	}
	return true
}

// ProductOfStarNorms returns prod_i |b*_i|, which for a full-rank basis equals
// the lattice covolume exactly (Gram-Schmidt preserves volume).
func (b Basis) ProductOfStarNorms() float64 {
	gs := b.Orthogonalize()
	p := 1.0
	for _, n2 := range gs.Norm2 {
		p *= math.Sqrt(n2)
	}
	return p
}
