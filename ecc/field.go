package ecc

import "math/big"

// eccMod returns a mod m normalized into the range [0, m) for a positive
// modulus m. It never mutates its arguments.
func eccMod(a, m *big.Int) *big.Int {
	r := new(big.Int).Mod(a, m)
	if r.Sign() < 0 {
		r.Add(r, m)
	}
	return r
}

// Mod returns the representative of a modulo p in the range [0, p). The modulus
// p must be positive. Unlike the raw big.Int remainder, the result is always
// non-negative.
func Mod(a, p *big.Int) *big.Int {
	if p.Sign() <= 0 {
		panic("ecc: Mod requires a positive modulus")
	}
	return eccMod(a, p)
}

// ModAdd returns (a + b) mod p in the range [0, p) for a positive modulus p.
func ModAdd(a, b, p *big.Int) *big.Int {
	return eccMod(new(big.Int).Add(a, b), p)
}

// ModSub returns (a - b) mod p in the range [0, p) for a positive modulus p.
func ModSub(a, b, p *big.Int) *big.Int {
	return eccMod(new(big.Int).Sub(a, b), p)
}

// ModMul returns (a * b) mod p in the range [0, p) for a positive modulus p.
func ModMul(a, b, p *big.Int) *big.Int {
	return eccMod(new(big.Int).Mul(a, b), p)
}

// ModNeg returns (-a) mod p in the range [0, p) for a positive modulus p.
func ModNeg(a, p *big.Int) *big.Int {
	return eccMod(new(big.Int).Neg(a), p)
}

// ModExp returns (base ^ exp) mod p for a positive modulus p. Negative
// exponents are supported when base is invertible modulo p; in that case the
// modular inverse of base is raised to the magnitude of exp.
func ModExp(base, exp, p *big.Int) *big.Int {
	if exp.Sign() < 0 {
		inv, ok := ModInverse(base, p)
		if !ok {
			panic("ecc: ModExp of non-invertible base with negative exponent")
		}
		return new(big.Int).Exp(inv, new(big.Int).Neg(exp), p)
	}
	return new(big.Int).Exp(base, exp, p)
}

// ModInverse returns the multiplicative inverse of a modulo p together with a
// boolean reporting whether the inverse exists. The inverse exists exactly when
// gcd(a, p) = 1. When it does not exist the returned integer is nil.
func ModInverse(a, p *big.Int) (*big.Int, bool) {
	inv := new(big.Int).ModInverse(eccMod(a, p), p)
	if inv == nil {
		return nil, false
	}
	return inv, true
}

// ModDiv returns (a / b) mod p, that is a times the modular inverse of b, along
// with a boolean reporting whether b is invertible modulo p. When b is not
// invertible the returned integer is nil.
func ModDiv(a, b, p *big.Int) (*big.Int, bool) {
	inv, ok := ModInverse(b, p)
	if !ok {
		return nil, false
	}
	return ModMul(a, inv, p), true
}

// ModSqrt returns a square root of a modulo an odd prime p together with a
// boolean reporting whether a is a quadratic residue. When a has no square root
// the returned integer is nil. The smaller of the two roots r and p-r is
// returned, so the result is canonical. This routine assumes p is prime.
func ModSqrt(a, p *big.Int) (*big.Int, bool) {
	am := eccMod(a, p)
	r := new(big.Int).ModSqrt(am, p)
	if r == nil {
		return nil, false
	}
	other := new(big.Int).Sub(p, r)
	if other.Cmp(r) < 0 {
		r = other
	}
	return r, true
}

// Legendre returns the Legendre symbol (a / p) for an odd prime p, one of the
// values -1, 0 or +1: 0 when p divides a, +1 when a is a non-zero quadratic
// residue modulo p, and -1 when a is a quadratic non-residue. It is implemented
// via the Jacobi symbol, which coincides with the Legendre symbol for prime p.
func Legendre(a, p *big.Int) int {
	return big.Jacobi(eccMod(a, p), p)
}

// IsQuadraticResidue reports whether a is a non-zero quadratic residue modulo
// the odd prime p, i.e. whether there exists x with x^2 = a (mod p) and
// a is not congruent to 0.
func IsQuadraticResidue(a, p *big.Int) bool {
	return Legendre(a, p) == 1
}

// eccISqrt returns the integer square root floor(sqrt(n)) for n >= 0.
func eccISqrt(n *big.Int) *big.Int {
	if n.Sign() < 0 {
		panic("ecc: integer square root of a negative number")
	}
	return new(big.Int).Sqrt(n)
}
