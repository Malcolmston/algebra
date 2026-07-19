package simplicial

import (
	"math/big"
	"strings"
)

// BettiNumber returns the k-th Betti number b_k of the complex, the rank of the
// k-th homology group with rational coefficients. It equals dim ker ∂_k −
// rank ∂_{k+1} = (n_k − rank ∂_k) − rank ∂_{k+1}, where n_k is the number of
// k-simplices. Negative k yields 0.
func (c *Complex) BettiNumber(k int) int {
	if k < 0 {
		return 0
	}
	nk := c.ChainRank(k)
	rk := c.BoundaryRankQ(k)      // rank ∂_k
	rk1 := c.BoundaryRankQ(k + 1) // rank ∂_{k+1}
	b := nk - rk - rk1
	if b < 0 {
		return 0
	}
	return b
}

// BettiNumberGF2 returns the k-th Betti number computed with GF(2) coefficients.
// It can differ from the rational Betti number in the presence of 2-torsion —
// for example the real projective plane has GF(2) Betti numbers (1,1,1) but
// rational Betti numbers (1,0,0).
func (c *Complex) BettiNumberGF2(k int) int {
	if k < 0 {
		return 0
	}
	nk := c.ChainRank(k)
	rk := c.BoundaryRankGF2(k)
	rk1 := c.BoundaryRankGF2(k + 1)
	b := nk - rk - rk1
	if b < 0 {
		return 0
	}
	return b
}

// BettiNumbers returns the rational Betti numbers b_0,…,b_d of the complex,
// where d is its dimension. The empty complex yields an empty slice.
func (c *Complex) BettiNumbers() []int {
	d := c.Dimension()
	if d < 0 {
		return nil
	}
	out := make([]int, d+1)
	for k := 0; k <= d; k++ {
		out[k] = c.BettiNumber(k)
	}
	return out
}

// BettiNumbersGF2 returns the GF(2) Betti numbers b_0,…,b_d of the complex.
func (c *Complex) BettiNumbersGF2() []int {
	d := c.Dimension()
	if d < 0 {
		return nil
	}
	out := make([]int, d+1)
	for k := 0; k <= d; k++ {
		out[k] = c.BettiNumberGF2(k)
	}
	return out
}

// HomologyRank is an alias for [Complex.BettiNumber]: the free rank of H_k.
func (c *Complex) HomologyRank(k int) int { return c.BettiNumber(k) }

// EulerCharacteristicFromBetti returns the Euler characteristic computed as the
// alternating sum of the Betti numbers, Σ_k (−1)^k b_k. By the Euler–Poincaré
// theorem it equals [Complex.EulerCharacteristic].
func (c *Complex) EulerCharacteristicFromBetti() int {
	chi := 0
	sign := 1
	for _, b := range c.BettiNumbers() {
		chi += sign * b
		sign = -sign
	}
	return chi
}

// TorsionCoefficients returns the torsion coefficients of the k-th integral
// homology group H_k(c; Z): the elementary divisors greater than one of the
// boundary map ∂_{k+1}, obtained from its Smith normal form. An empty result
// means H_k is free. For example the real projective plane returns [2] for
// k = 1.
func (c *Complex) TorsionCoefficients(k int) []*big.Int {
	if k < 0 {
		return nil
	}
	b := c.BoundaryMatrixZ(k + 1)
	if b.cols == 0 || b.rows == 0 {
		return nil
	}
	return b.SmithNormalForm().TorsionFactors()
}

// HomologyGroup describes a finitely generated abelian group H_k = Z^r ⊕
// (Z/t_1 ⊕ … ⊕ Z/t_m) by its free rank r and its torsion coefficients
// t_1 | t_2 | … | t_m, each greater than one.
type HomologyGroup struct {
	// FreeRank is the rank of the free part (the Betti number).
	FreeRank int
	// Torsion lists the torsion coefficients, each greater than one and each
	// dividing the next.
	Torsion []*big.Int
}

// Betti returns the free rank of the group.
func (h HomologyGroup) Betti() int { return h.FreeRank }

// IsFree reports whether the group is torsion-free.
func (h HomologyGroup) IsFree() bool { return len(h.Torsion) == 0 }

// IsTrivial reports whether the group is the zero group.
func (h HomologyGroup) IsTrivial() bool { return h.FreeRank == 0 && len(h.Torsion) == 0 }

// String renders the group in the usual notation, e.g. "Z^2 + Z/2" or "0".
func (h HomologyGroup) String() string {
	var parts []string
	switch {
	case h.FreeRank == 1:
		parts = append(parts, "Z")
	case h.FreeRank > 1:
		parts = append(parts, "Z^"+itoa(h.FreeRank))
	}
	for _, t := range h.Torsion {
		parts = append(parts, "Z/"+t.String())
	}
	if len(parts) == 0 {
		return "0"
	}
	return strings.Join(parts, " + ")
}

// HomologyZ returns the k-th integral homology group H_k(c; Z) as a
// [HomologyGroup], combining the Betti number (free rank) with the torsion
// coefficients from the Smith normal form of ∂_{k+1}.
func (c *Complex) HomologyZ(k int) HomologyGroup {
	return HomologyGroup{
		FreeRank: c.BettiNumber(k),
		Torsion:  c.TorsionCoefficients(k),
	}
}

// Homology returns the integral homology groups H_0,…,H_d of the complex.
func (c *Complex) Homology() []HomologyGroup {
	d := c.Dimension()
	if d < 0 {
		return nil
	}
	out := make([]HomologyGroup, d+1)
	for k := 0; k <= d; k++ {
		out[k] = c.HomologyZ(k)
	}
	return out
}

// itoa is a tiny helper avoiding an fmt dependency in hot paths.
func itoa(n int) string {
	return big.NewInt(int64(n)).String()
}
