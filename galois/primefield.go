package galois

import (
	"errors"
	"fmt"
	"math/big"
)

// Mod returns a reduced into the canonical range [0, p) for the modulus p > 0.
func Mod(a, p *big.Int) *big.Int {
	return new(big.Int).Mod(a, p)
}

// AddMod returns (a + b) mod p in the range [0, p).
func AddMod(a, b, p *big.Int) *big.Int {
	r := new(big.Int).Add(a, b)
	return r.Mod(r, p)
}

// SubMod returns (a - b) mod p in the range [0, p).
func SubMod(a, b, p *big.Int) *big.Int {
	r := new(big.Int).Sub(a, b)
	return r.Mod(r, p)
}

// MulMod returns (a * b) mod p in the range [0, p).
func MulMod(a, b, p *big.Int) *big.Int {
	r := new(big.Int).Mul(a, b)
	return r.Mod(r, p)
}

// NegMod returns (-a) mod p in the range [0, p).
func NegMod(a, p *big.Int) *big.Int {
	r := new(big.Int).Neg(a)
	return r.Mod(r, p)
}

// InvMod returns the multiplicative inverse of a modulo p, or an error when a
// is not invertible (gcd(a, p) != 1).
func InvMod(a, p *big.Int) (*big.Int, error) {
	r := new(big.Int).ModInverse(new(big.Int).Mod(a, p), p)
	if r == nil {
		return nil, fmt.Errorf("galois: %s has no inverse mod %s", a.String(), p.String())
	}
	return r, nil
}

// DivMod returns (a / b) mod p, that is a times the inverse of b, or an error
// when b is not invertible.
func DivMod(a, b, p *big.Int) (*big.Int, error) {
	inv, err := InvMod(b, p)
	if err != nil {
		return nil, err
	}
	return MulMod(a, inv, p), nil
}

// PowMod returns a raised to the power e modulo p. Negative exponents invert a
// first and therefore require a to be invertible.
func PowMod(a, e, p *big.Int) (*big.Int, error) {
	if e.Sign() < 0 {
		inv, err := InvMod(a, p)
		if err != nil {
			return nil, err
		}
		return new(big.Int).Exp(inv, new(big.Int).Neg(e), p), nil
	}
	return new(big.Int).Exp(a, e, p), nil
}

// Jacobi returns the Jacobi symbol (a/n) for odd n > 0, one of -1, 0, +1.
func Jacobi(a, n *big.Int) int {
	return big.Jacobi(a, n)
}

// Legendre returns the Legendre symbol (a/p) for an odd prime p: it is 0 when
// p divides a, +1 when a is a non-zero quadratic residue, and -1 otherwise.
func Legendre(a, p *big.Int) int {
	return big.Jacobi(new(big.Int).Mod(a, p), p)
}

// IsQuadraticResidue reports whether a is a square modulo the prime p. Zero is
// treated as a residue.
func IsQuadraticResidue(a, p *big.Int) bool {
	am := new(big.Int).Mod(a, p)
	if am.Sign() == 0 {
		return true
	}
	return big.Jacobi(am, p) == 1
}

// SqrtMod returns a square root r of a modulo the prime p, so that
// r*r ≡ a (mod p). It uses the Tonelli–Shanks algorithm and returns an error
// when a is a non-residue. When a root exists, p-r is the other root.
func SqrtMod(a, p *big.Int) (*big.Int, error) {
	am := new(big.Int).Mod(a, p)
	if am.Sign() == 0 {
		return big.NewInt(0), nil
	}
	if p.Cmp(big2) == 0 {
		return am, nil
	}
	if big.Jacobi(am, p) != 1 {
		return nil, errors.New("galois: argument is not a quadratic residue")
	}
	// p ≡ 3 (mod 4): r = a^((p+1)/4).
	if new(big.Int).And(p, big.NewInt(3)).Cmp(big.NewInt(3)) == 0 {
		e := new(big.Int).Add(p, big1)
		e.Rsh(e, 2)
		return new(big.Int).Exp(am, e, p), nil
	}
	// Write p-1 = q * 2^s with q odd.
	q := new(big.Int).Sub(p, big1)
	s := 0
	for q.Bit(0) == 0 {
		q.Rsh(q, 1)
		s++
	}
	// Find a non-residue z.
	z := big.NewInt(2)
	for big.Jacobi(z, p) != -1 {
		z.Add(z, big1)
	}
	m := s
	c := new(big.Int).Exp(z, q, p)
	t := new(big.Int).Exp(am, q, p)
	rExp := new(big.Int).Add(q, big1)
	rExp.Rsh(rExp, 1)
	r := new(big.Int).Exp(am, rExp, p)
	for t.Cmp(big1) != 0 {
		// find least i, 0 < i < m, with t^(2^i) == 1.
		i := 0
		tt := clone(t)
		for tt.Cmp(big1) != 0 {
			tt.Mul(tt, tt)
			tt.Mod(tt, p)
			i++
			if i == m {
				return nil, errors.New("galois: SqrtMod failed to converge")
			}
		}
		// b = c^(2^(m-i-1)).
		b := clone(c)
		for j := 0; j < m-i-1; j++ {
			b.Mul(b, b)
			b.Mod(b, p)
		}
		m = i
		c = new(big.Int).Mul(b, b)
		c.Mod(c, p)
		t.Mul(t, c)
		t.Mod(t, p)
		r.Mul(r, b)
		r.Mod(r, p)
	}
	return r, nil
}

// MultiplicativeOrder returns the multiplicative order of a in the group of
// units modulo the prime p: the smallest positive k with a^k ≡ 1 (mod p). It
// returns an error when a is 0 modulo p. The modulus is assumed prime, so the
// group order is p-1.
func MultiplicativeOrder(a, p *big.Int) (*big.Int, error) {
	am := new(big.Int).Mod(a, p)
	if am.Sign() == 0 {
		return nil, errors.New("galois: 0 has no multiplicative order")
	}
	order := new(big.Int).Sub(p, big1)
	for _, q := range PrimeFactors(order) {
		for new(big.Int).Mod(order, q).Sign() == 0 {
			cand := new(big.Int).Div(order, q)
			if new(big.Int).Exp(am, cand, p).Cmp(big1) == 0 {
				order = cand
			} else {
				break
			}
		}
	}
	return order, nil
}

// IsPrimitiveRoot reports whether a is a primitive root modulo the prime p,
// i.e. a generator of the cyclic group of units of order p-1.
func IsPrimitiveRoot(a, p *big.Int) bool {
	am := new(big.Int).Mod(a, p)
	if am.Sign() == 0 {
		return false
	}
	phi := new(big.Int).Sub(p, big1)
	for _, q := range PrimeFactors(phi) {
		e := new(big.Int).Div(phi, q)
		if new(big.Int).Exp(am, e, p).Cmp(big1) == 0 {
			return false
		}
	}
	return true
}

// PrimitiveRoot returns the smallest primitive root modulo the prime p.
func PrimitiveRoot(p *big.Int) (*big.Int, error) {
	if p.Cmp(big2) == 0 {
		return big.NewInt(1), nil
	}
	if !IsPrimeInt(p) {
		return nil, errors.New("galois: PrimitiveRoot requires a prime modulus")
	}
	g := big.NewInt(2)
	for g.Cmp(p) < 0 {
		if IsPrimitiveRoot(g, p) {
			return clone(g), nil
		}
		g.Add(g, big1)
	}
	return nil, errors.New("galois: no primitive root found")
}
