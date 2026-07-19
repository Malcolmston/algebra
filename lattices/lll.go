package lattices

import "math"

// DefaultDelta is the standard Lovasz parameter used by LLL reduction. Values
// in the open interval (1/4, 1) guarantee termination; 3/4 is the classical
// choice.
const DefaultDelta = 0.75

// reduceRow subtracts round(mu[k][l]) times basis vector l from vector k. It
// reports whether the row changed.
func reduceRow(b Basis, k, l int, mu [][]float64) bool {
	q := math.Round(mu[k][l])
	if q == 0 {
		return false
	}
	b[k] = b[k].AddScaled(b[l], -q)
	return true
}

// SizeReduced returns a copy of the basis that is size reduced: every
// Gram-Schmidt coefficient mu[i][j] with j < i is made to satisfy
// |mu[i][j]| <= 1/2 by subtracting integer multiples of earlier vectors. Size
// reduction does not reorder the basis and leaves the lattice unchanged.
func (b Basis) SizeReduced() Basis {
	c := b.Clone()
	n := len(c)
	for i := 1; i < n; i++ {
		for j := i - 1; j >= 0; j-- {
			gs := c.Orthogonalize()
			reduceRow(c, i, j, gs.Mu)
		}
	}
	return c
}

// LLL returns an LLL-reduced basis for the same lattice using the given Lovasz
// parameter delta, which must lie in (1/4, 1]. It applies the classical
// Lenstra-Lenstra-Lovasz algorithm: repeated size reduction interleaved with
// swaps whenever the Lovasz condition fails. It returns ErrBadParameter for an
// out-of-range delta and ErrEmpty for an empty basis.
func (b Basis) LLL(delta float64) (Basis, error) {
	if delta <= 0.25 || delta > 1 {
		return nil, ErrBadParameter
	}
	if len(b) == 0 {
		return nil, ErrEmpty
	}
	c := b.Clone()
	n := len(c)
	if n < 2 {
		return c, nil
	}
	k := 1
	for k < n {
		gs := c.Orthogonalize()
		reduceRow(c, k, k-1, gs.Mu)
		gs = c.Orthogonalize()
		if gs.Norm2[k] >= (delta-gs.Mu[k][k-1]*gs.Mu[k][k-1])*gs.Norm2[k-1] {
			for l := k - 2; l >= 0; l-- {
				g2 := c.Orthogonalize()
				reduceRow(c, k, l, g2.Mu)
			}
			k++
		} else {
			c.Swap(k, k-1)
			if k > 1 {
				k--
			}
		}
	}
	return c, nil
}

// LLLDefault returns an LLL-reduced basis using DefaultDelta (3/4).
func (b Basis) LLLDefault() Basis {
	r, _ := b.LLL(DefaultDelta)
	return r
}

// IsLLLReduced reports whether the basis is LLL reduced for parameter delta
// within tolerance tol: it must be size reduced and satisfy the Lovasz
// condition |b*_k|^2 >= (delta - mu[k][k-1]^2) |b*_{k-1}|^2 for every k.
func (b Basis) IsLLLReduced(delta, tol float64) bool {
	gs := b.Orthogonalize()
	n := len(b)
	for i := 0; i < n; i++ {
		for j := 0; j < i; j++ {
			if math.Abs(gs.Mu[i][j]) > 0.5+tol {
				return false
			}
		}
	}
	for k := 1; k < n; k++ {
		lhs := gs.Norm2[k]
		rhs := (delta - gs.Mu[k][k-1]*gs.Mu[k][k-1]) * gs.Norm2[k-1]
		if lhs < rhs-tol {
			return false
		}
	}
	return true
}

// LLLGramDefect returns the orthogonality defect of the LLL-reduced basis, a
// convenient quality measure after reduction.
func (b Basis) LLLGramDefect() float64 {
	return b.LLLDefault().OrthogonalityDefect()
}
