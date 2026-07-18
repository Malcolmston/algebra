package crypto

import (
	"errors"
	"math/big"
)

// Jacobi returns the Jacobi symbol (a/n) for an odd positive integer n. The
// value is one of -1, 0 or +1. When n is an odd prime the Jacobi symbol equals
// the Legendre symbol. Jacobi panics if n is even or non-positive.
func Jacobi(a, n *big.Int) int {
	if n.Sign() <= 0 || n.Bit(0) == 0 {
		panic("crypto: Jacobi requires an odd positive modulus n")
	}
	x := new(big.Int).Mod(a, n)
	y := new(big.Int).Set(n)
	result := 1
	for x.Sign() != 0 {
		// Factor out powers of two from x, tracking the sign contribution.
		for x.Bit(0) == 0 {
			x.Rsh(x, 1)
			m8 := y.Bit(0) + y.Bit(1)*2 + y.Bit(2)*4 // y mod 8
			if m8 == 3 || m8 == 5 {                  // y ≡ 3 or 5 (mod 8)
				result = -result
			}
		}
		x, y = y, x
		// Quadratic reciprocity flip.
		if x.Bit(0) == 1 && y.Bit(0) == 1 && x.Bit(1) == 1 && y.Bit(1) == 1 {
			// both ≡ 3 (mod 4)
			result = -result
		}
		x.Mod(x, y)
	}
	if y.Cmp(cryptoOne) == 0 {
		return result
	}
	return 0
}

// Legendre returns the Legendre symbol (a/p) for an odd prime p: +1 if a is a
// non-zero quadratic residue mod p, -1 if it is a non-residue, and 0 if p
// divides a. The result is computed via Euler's criterion a^((p-1)/2) mod p.
// Legendre panics if p is not an odd prime candidate (even or <= 2); it does
// not itself verify primality, so the result is only meaningful for prime p.
func Legendre(a, p *big.Int) int {
	if p.Cmp(cryptoTwo) <= 0 || p.Bit(0) == 0 {
		panic("crypto: Legendre requires an odd prime p > 2")
	}
	e := new(big.Int).Sub(p, cryptoOne)
	e.Rsh(e, 1)
	r := ModExp(a, e, p)
	switch {
	case r.Sign() == 0:
		return 0
	case r.Cmp(cryptoOne) == 0:
		return 1
	default:
		return -1
	}
}

// IsQuadraticResidue reports whether a is a quadratic residue modulo the odd
// prime p, i.e. whether there exists x with x^2 ≡ a (mod p). Zero is treated as
// a residue (0 = 0^2). It panics under the same conditions as Legendre.
func IsQuadraticResidue(a, p *big.Int) bool {
	return Legendre(a, p) >= 0
}

// ModSqrt returns a square root r of a modulo the odd prime p (so that r^2 ≡ a
// (mod p)), using the Tonelli-Shanks algorithm. When a solution exists the
// smaller of the two roots r and p-r is returned. It returns an error if a is a
// quadratic non-residue modulo p. It panics if p is even or <= 2.
func ModSqrt(a, p *big.Int) (*big.Int, error) {
	if p.Cmp(cryptoTwo) <= 0 || p.Bit(0) == 0 {
		panic("crypto: ModSqrt requires an odd prime p > 2")
	}
	aa := new(big.Int).Mod(a, p)
	if aa.Sign() == 0 {
		return big.NewInt(0), nil
	}
	if Legendre(aa, p) == -1 {
		return nil, errors.New("crypto: ModSqrt argument is a quadratic non-residue")
	}
	// Case p ≡ 3 (mod 4): direct formula r = a^((p+1)/4).
	if new(big.Int).And(p, cryptoThree).Cmp(cryptoThree) == 0 {
		e := new(big.Int).Add(p, cryptoOne)
		e.Rsh(e, 2)
		r := ModExp(aa, e, p)
		return cryptoSmallerRoot(r, p), nil
	}
	// Write p-1 = q * 2^s with q odd.
	q := new(big.Int).Sub(p, cryptoOne)
	s := 0
	for q.Bit(0) == 0 {
		q.Rsh(q, 1)
		s++
	}
	// Find a quadratic non-residue z.
	z := big.NewInt(2)
	for Legendre(z, p) != -1 {
		z.Add(z, cryptoOne)
	}
	c := ModExp(z, q, p)
	e := new(big.Int).Add(q, cryptoOne)
	e.Rsh(e, 1)
	r := ModExp(aa, e, p)
	t := ModExp(aa, q, p)
	m := s
	for t.Cmp(cryptoOne) != 0 {
		// Find least i, 0 < i < m, with t^(2^i) ≡ 1.
		i := 0
		t2 := new(big.Int).Set(t)
		for t2.Cmp(cryptoOne) != 0 {
			t2.Mul(t2, t2)
			t2.Mod(t2, p)
			i++
			if i == m {
				return nil, errors.New("crypto: ModSqrt failed to converge")
			}
		}
		b := new(big.Int).Set(c)
		for j := 0; j < m-i-1; j++ {
			b.Mul(b, b)
			b.Mod(b, p)
		}
		r.Mul(r, b)
		r.Mod(r, p)
		c.Mul(b, b)
		c.Mod(c, p)
		t.Mul(t, c)
		t.Mod(t, p)
		m = i
	}
	return cryptoSmallerRoot(r, p), nil
}

// cryptoSmallerRoot returns whichever of r and p-r is numerically smaller.
func cryptoSmallerRoot(r, p *big.Int) *big.Int {
	other := new(big.Int).Sub(p, r)
	if other.Cmp(r) < 0 {
		return other
	}
	return new(big.Int).Set(r)
}

// MultiplicativeOrder returns the multiplicative order of a modulo n, that is
// the smallest positive integer k with a^k ≡ 1 (mod n). It requires gcd(a, n)
// = 1 and n > 1, returning an error otherwise. The order is found by reducing
// the Carmichael function λ(n) (a known multiple of every element's order) over
// its prime factors.
func MultiplicativeOrder(a, n *big.Int) (*big.Int, error) {
	if n.Cmp(cryptoOne) <= 0 {
		return nil, errors.New("crypto: MultiplicativeOrder requires n > 1")
	}
	if GCD(a, n).Cmp(cryptoOne) != 0 {
		return nil, errors.New("crypto: MultiplicativeOrder requires gcd(a, n) == 1")
	}
	order := CarmichaelLambda(n)
	factors := Factorization(order)
	for _, f := range factors {
		for i := 0; i < f.Exponent; i++ {
			cand := new(big.Int).Div(order, f.Prime)
			if ModExp(a, cand, n).Cmp(cryptoOne) == 0 {
				order = cand
			} else {
				break
			}
		}
	}
	return order, nil
}

// PrimitiveRoot returns the smallest primitive root modulo the odd prime p,
// i.e. the least generator g of the cyclic multiplicative group (Z/pZ)* whose
// multiplicative order equals p-1. It panics if p <= 2. The result is only
// meaningful when p is prime.
func PrimitiveRoot(p *big.Int) *big.Int {
	if p.Cmp(cryptoTwo) <= 0 {
		panic("crypto: PrimitiveRoot requires a prime p > 2")
	}
	phi := new(big.Int).Sub(p, cryptoOne)
	factors := Factorization(phi)
	g := big.NewInt(2)
	for g.Cmp(p) < 0 {
		if cryptoIsGenerator(g, phi, p, factors) {
			return new(big.Int).Set(g)
		}
		g.Add(g, cryptoOne)
	}
	panic("crypto: PrimitiveRoot found no generator (p is not prime)")
}

// IsPrimitiveRoot reports whether g is a primitive root modulo the odd prime p,
// i.e. whether its multiplicative order is exactly p-1. It panics if p <= 2 and
// is only meaningful for prime p.
func IsPrimitiveRoot(g, p *big.Int) bool {
	if p.Cmp(cryptoTwo) <= 0 {
		panic("crypto: IsPrimitiveRoot requires a prime p > 2")
	}
	gg := new(big.Int).Mod(g, p)
	if gg.Sign() == 0 {
		return false
	}
	phi := new(big.Int).Sub(p, cryptoOne)
	return cryptoIsGenerator(gg, phi, p, Factorization(phi))
}

// cryptoIsGenerator reports whether g generates a group of order phi modulo p,
// given the prime factorization of phi. g is a generator exactly when
// g^(phi/q) != 1 for every distinct prime q dividing phi.
func cryptoIsGenerator(g, phi, p *big.Int, factors []Factor) bool {
	for _, f := range factors {
		e := new(big.Int).Div(phi, f.Prime)
		if ModExp(g, e, p).Cmp(cryptoOne) == 0 {
			return false
		}
	}
	return true
}
