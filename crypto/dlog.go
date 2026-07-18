package crypto

import (
	"errors"
	"math/big"
)

// ErrNoDiscreteLog is returned when a discrete logarithm does not exist for the
// given base, target and modulus.
var ErrNoDiscreteLog = errors.New("crypto: no discrete logarithm exists")

// BabyStepGiantStep solves the discrete logarithm g^x ≡ h (mod p) for x in
// [0, order), using Shanks' baby-step giant-step algorithm in O(sqrt(order))
// time and space. order must be a positive multiple of the multiplicative order
// of g (for a prime modulus p, passing order = p-1 always works). It returns
// the least such non-negative x, or ErrNoDiscreteLog if none exists. The
// modulus p must be positive.
func BabyStepGiantStep(g, h, p, order *big.Int) (*big.Int, error) {
	if p.Sign() <= 0 {
		panic("crypto: BabyStepGiantStep requires modulus p > 0")
	}
	if order.Sign() <= 0 {
		panic("crypto: BabyStepGiantStep requires order > 0")
	}
	gg := new(big.Int).Mod(g, p)
	hh := new(big.Int).Mod(h, p)

	// m = ceil(sqrt(order)).
	m := new(big.Int).Sqrt(order)
	if new(big.Int).Mul(m, m).Cmp(order) < 0 {
		m.Add(m, cryptoOne)
	}

	// Baby steps: table maps g^j -> j for j in [0, m).
	table := make(map[string]int64)
	e := big.NewInt(1)
	mInt := m.Int64()
	for j := int64(0); j < mInt; j++ {
		key := e.String()
		if _, seen := table[key]; !seen {
			table[key] = j
		}
		e.Mul(e, gg)
		e.Mod(e, p)
	}

	// factor = g^(-m) mod p.
	gm := ModExp(gg, m, p)
	factor := new(big.Int).ModInverse(gm, p)
	if factor == nil {
		return nil, ErrNoDiscreteLog
	}

	gamma := new(big.Int).Set(hh)
	for i := int64(0); i <= mInt; i++ {
		if j, ok := table[gamma.String()]; ok {
			x := new(big.Int).SetInt64(i)
			x.Mul(x, m)
			x.Add(x, big.NewInt(j))
			x.Mod(x, order)
			return x, nil
		}
		gamma.Mul(gamma, factor)
		gamma.Mod(gamma, p)
	}
	return nil, ErrNoDiscreteLog
}

// DiscreteLog solves g^x ≡ h (mod p) for a prime modulus p, searching the full
// range [0, p-1). It is a convenience wrapper around BabyStepGiantStep with
// order = p-1. It returns ErrNoDiscreteLog if no solution exists.
func DiscreteLog(g, h, p *big.Int) (*big.Int, error) {
	order := new(big.Int).Sub(p, cryptoOne)
	return BabyStepGiantStep(g, h, p, order)
}
