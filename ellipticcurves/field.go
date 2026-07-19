package ellipticcurves

import (
	"errors"
	"math/big"
	"math/rand"
)

// Common errors returned by the field-arithmetic helpers.
var (
	// ErrNotPrime indicates that a modulus expected to be an odd prime was not.
	ErrNotPrime = errors.New("ellipticcurves: modulus is not an odd prime")
	// ErrNoSquareRoot indicates that an element is a quadratic non-residue.
	ErrNoSquareRoot = errors.New("ellipticcurves: element has no square root")
	// ErrNotInvertible indicates that an element is not invertible modulo p.
	ErrNotInvertible = errors.New("ellipticcurves: element is not invertible")
	// ErrZeroModulus indicates that a zero or negative modulus was supplied.
	ErrZeroModulus = errors.New("ellipticcurves: modulus must be positive")
)

var (
	bigZero  = big.NewInt(0)
	bigOne   = big.NewInt(1)
	bigTwo   = big.NewInt(2)
	bigThree = big.NewInt(3)
	bigFour  = big.NewInt(4)
)

// Mod returns the canonical representative of a modulo m in the range [0, m).
// The result is always non-negative regardless of the sign of a. It does not
// mutate its arguments.
func Mod(a, m *big.Int) *big.Int {
	r := new(big.Int).Mod(a, m)
	if r.Sign() < 0 {
		r.Add(r, m)
	}
	return r
}

// ModAdd returns (a + b) mod m reduced into [0, m).
func ModAdd(a, b, m *big.Int) *big.Int {
	return Mod(new(big.Int).Add(a, b), m)
}

// ModSub returns (a - b) mod m reduced into [0, m).
func ModSub(a, b, m *big.Int) *big.Int {
	return Mod(new(big.Int).Sub(a, b), m)
}

// ModMul returns (a * b) mod m reduced into [0, m).
func ModMul(a, b, m *big.Int) *big.Int {
	return Mod(new(big.Int).Mul(a, b), m)
}

// ModNeg returns (-a) mod m reduced into [0, m).
func ModNeg(a, m *big.Int) *big.Int {
	return Mod(new(big.Int).Neg(a), m)
}

// ModDouble returns (2*a) mod m reduced into [0, m).
func ModDouble(a, m *big.Int) *big.Int {
	return Mod(new(big.Int).Lsh(a, 1), m)
}

// ModSquare returns (a*a) mod m reduced into [0, m).
func ModSquare(a, m *big.Int) *big.Int {
	return Mod(new(big.Int).Mul(a, a), m)
}

// ModExp returns a^e mod m using big.Int.Exp. Negative exponents are supported
// when a is invertible modulo m; otherwise ModExp returns ErrNotInvertible.
func ModExp(a, e, m *big.Int) (*big.Int, error) {
	if e.Sign() >= 0 {
		return new(big.Int).Exp(a, e, m), nil
	}
	inv, err := ModInverse(a, m)
	if err != nil {
		return nil, err
	}
	pe := new(big.Int).Neg(e)
	return new(big.Int).Exp(inv, pe, m), nil
}

// ModInverse returns the multiplicative inverse of a modulo m, or
// ErrNotInvertible when gcd(a, m) != 1.
func ModInverse(a, m *big.Int) (*big.Int, error) {
	inv := new(big.Int).ModInverse(Mod(a, m), m)
	if inv == nil {
		return nil, ErrNotInvertible
	}
	return inv, nil
}

// ModDiv returns (a * b^-1) mod m, or ErrNotInvertible when b is not
// invertible modulo m.
func ModDiv(a, b, m *big.Int) (*big.Int, error) {
	inv, err := ModInverse(b, m)
	if err != nil {
		return nil, err
	}
	return ModMul(a, inv, m), nil
}

// LegendreSymbol returns the Legendre symbol (a/p) for an odd prime p, one of
// -1, 0 or 1. It computes a^((p-1)/2) mod p. The caller is responsible for
// supplying a prime p; use JacobiSymbol for general odd moduli.
func LegendreSymbol(a, p *big.Int) int {
	am := Mod(a, p)
	if am.Sign() == 0 {
		return 0
	}
	e := new(big.Int).Rsh(new(big.Int).Sub(p, bigOne), 1)
	r := new(big.Int).Exp(am, e, p)
	if r.Cmp(bigOne) == 0 {
		return 1
	}
	return -1
}

// JacobiSymbol returns the Jacobi symbol (a/n) for a positive odd integer n.
// For prime n it coincides with the Legendre symbol. It returns 0 when
// gcd(a, n) != 1. The behaviour for even or non-positive n is undefined and
// reported as 0.
func JacobiSymbol(a, n *big.Int) int {
	if n.Sign() <= 0 || n.Bit(0) == 0 {
		return 0
	}
	x := Mod(a, n)
	y := new(big.Int).Set(n)
	result := 1
	for x.Sign() != 0 {
		for x.Bit(0) == 0 {
			x.Rsh(x, 1)
			r := new(big.Int).Mod(y, big.NewInt(8)).Int64()
			if r == 3 || r == 5 {
				result = -result
			}
		}
		x, y = y, x
		if x.Bit(0) == 1 && y.Bit(0) == 1 {
			xm := new(big.Int).Mod(x, bigFour).Int64()
			ym := new(big.Int).Mod(y, bigFour).Int64()
			if xm == 3 && ym == 3 {
				result = -result
			}
		}
		x.Mod(x, y)
	}
	if y.Cmp(bigOne) == 0 {
		return result
	}
	return 0
}

// IsQuadraticResidue reports whether a is a non-zero quadratic residue modulo
// the odd prime p. Zero is not considered a residue by this predicate.
func IsQuadraticResidue(a, p *big.Int) bool {
	return LegendreSymbol(a, p) == 1
}

// ModSqrt returns a square root of a modulo the odd prime p, or ErrNoSquareRoot
// when a is a quadratic non-residue. When a root r is returned, p-r is the
// other root. It delegates to big.Int.ModSqrt.
func ModSqrt(a, p *big.Int) (*big.Int, error) {
	am := Mod(a, p)
	r := new(big.Int).ModSqrt(am, p)
	if r == nil {
		return nil, ErrNoSquareRoot
	}
	return r, nil
}

// ModSqrtTonelli returns a square root of a modulo the odd prime p using an
// explicit Tonelli-Shanks implementation. It returns ErrNoSquareRoot for a
// non-residue. The result matches ModSqrt up to the choice of sign.
func ModSqrtTonelli(a, p *big.Int) (*big.Int, error) {
	am := Mod(a, p)
	if am.Sign() == 0 {
		return big.NewInt(0), nil
	}
	if LegendreSymbol(am, p) != 1 {
		return nil, ErrNoSquareRoot
	}
	// p = 3 mod 4 fast path.
	if new(big.Int).Mod(p, bigFour).Int64() == 3 {
		e := new(big.Int).Rsh(new(big.Int).Add(p, bigOne), 2)
		return new(big.Int).Exp(am, e, p), nil
	}
	// Write p-1 = q * 2^s with q odd.
	q := new(big.Int).Sub(p, bigOne)
	s := 0
	for q.Bit(0) == 0 {
		q.Rsh(q, 1)
		s++
	}
	// Find a non-residue z.
	z := big.NewInt(2)
	for LegendreSymbol(z, p) != -1 {
		z.Add(z, bigOne)
	}
	m := s
	c := new(big.Int).Exp(z, q, p)
	t := new(big.Int).Exp(am, q, p)
	r := new(big.Int).Exp(am, new(big.Int).Rsh(new(big.Int).Add(q, bigOne), 1), p)
	for t.Cmp(bigOne) != 0 {
		// Find least i, 0 < i < m, with t^(2^i) = 1.
		i := 0
		tt := new(big.Int).Set(t)
		for tt.Cmp(bigOne) != 0 {
			tt = ModSquare(tt, p)
			i++
			if i == m {
				return nil, ErrNoSquareRoot
			}
		}
		b := new(big.Int).Set(c)
		for j := 0; j < m-i-1; j++ {
			b = ModSquare(b, p)
		}
		m = i
		c = ModSquare(b, p)
		t = ModMul(t, c, p)
		r = ModMul(r, b, p)
	}
	return r, nil
}

// SqrtCandidates returns both square roots of a modulo the odd prime p as
// (r, p-r). It returns ErrNoSquareRoot when a is a non-residue. When a == 0 the
// two roots coincide at 0.
func SqrtCandidates(a, p *big.Int) (*big.Int, *big.Int, error) {
	r, err := ModSqrt(a, p)
	if err != nil {
		return nil, nil, err
	}
	other := ModNeg(r, p)
	return r, other, nil
}

// MultiplicativeOrder returns the multiplicative order of a modulo the odd
// prime p, that is the least k > 0 with a^k = 1 mod p. It returns 0 when a is
// not invertible modulo p. The order divides p-1; the routine factors p-1 and
// removes prime factors to find the exact order.
func MultiplicativeOrder(a, p *big.Int) *big.Int {
	am := Mod(a, p)
	if am.Sign() == 0 {
		return big.NewInt(0)
	}
	order := new(big.Int).Sub(p, bigOne)
	factors := Factorize(order)
	for prime, exp := range factors {
		for e := 0; e < exp; e++ {
			cand := new(big.Int).Div(order, prime)
			if new(big.Int).Exp(am, cand, p).Cmp(bigOne) == 0 {
				order = cand
			} else {
				break
			}
		}
	}
	return order
}

// PrimitiveRoot returns a generator of the multiplicative group modulo the odd
// prime p, i.e. an element of order p-1. It searches upward from 2.
func PrimitiveRoot(p *big.Int) *big.Int {
	target := new(big.Int).Sub(p, bigOne)
	g := big.NewInt(2)
	for g.Cmp(p) < 0 {
		if MultiplicativeOrder(g, p).Cmp(target) == 0 {
			return new(big.Int).Set(g)
		}
		g.Add(g, bigOne)
	}
	return nil
}

// IsProbablePrime reports whether n is prime according to the Miller-Rabin and
// Baillie-PSW tests in big.Int.ProbablyPrime with the given number of rounds.
func IsProbablePrime(n *big.Int, rounds int) bool {
	if n.Sign() <= 0 {
		return false
	}
	return n.ProbablyPrime(rounds)
}

// NextPrime returns the least prime strictly greater than n.
func NextPrime(n *big.Int) *big.Int {
	c := new(big.Int).Add(n, bigOne)
	if c.Cmp(bigTwo) <= 0 {
		return big.NewInt(2)
	}
	if c.Bit(0) == 0 {
		c.Add(c, bigOne)
	}
	for !c.ProbablyPrime(20) {
		c.Add(c, bigTwo)
	}
	return c
}

// RandomFieldElement returns a uniformly random element of F_p using the
// supplied source. The result lies in [0, p).
func RandomFieldElement(p *big.Int, rng *rand.Rand) *big.Int {
	return new(big.Int).Rand(rng, p)
}

// RandomNonzeroFieldElement returns a uniformly random element of F_p^* using
// the supplied source. The result lies in [1, p).
func RandomNonzeroFieldElement(p *big.Int, rng *rand.Rand) *big.Int {
	if p.Cmp(bigTwo) <= 0 {
		return big.NewInt(1)
	}
	m := new(big.Int).Sub(p, bigOne)
	r := new(big.Int).Rand(rng, m)
	return r.Add(r, bigOne)
}
