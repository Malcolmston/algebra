package padic

import (
	"errors"
	"math/big"
	"math/rand"
)

// Common package errors.
var (
	// ErrNotPrime is returned when a modulus expected to be prime is not.
	ErrNotPrime = errors.New("padic: modulus is not prime")
	// ErrZeroDivision is returned by division or inversion of the p-adic zero.
	ErrZeroDivision = errors.New("padic: division by zero")
	// ErrNotInvertible is returned when a modular inverse does not exist.
	ErrNotInvertible = errors.New("padic: element not invertible")
	// ErrPrimeMismatch is returned when combining p-adics for different primes.
	ErrPrimeMismatch = errors.New("padic: mismatched primes")
	// ErrNoRoot is returned when a requested root does not exist.
	ErrNoRoot = errors.New("padic: root does not exist")
	// ErrDomain is returned when an argument is outside a function's domain.
	ErrDomain = errors.New("padic: argument out of domain")
	// ErrPrecision is returned when the requested precision is non-positive.
	ErrPrecision = errors.New("padic: precision must be positive")
)

var (
	bigZero = big.NewInt(0)
	bigOne  = big.NewInt(1)
	bigTwo  = big.NewInt(2)
)

// PPow returns p raised to the non-negative power k as a fresh big.Int. It
// panics only if k is negative, which is a programming error.
func PPow(p *big.Int, k int) *big.Int {
	if k < 0 {
		k = 0
	}
	return new(big.Int).Exp(p, big.NewInt(int64(k)), nil)
}

// PowMod returns base^exp mod m using fast modular exponentiation. A negative
// exp is interpreted as (base^-1)^|exp| when base is invertible mod m; if it
// is not invertible, the result is nil.
func PowMod(base, exp, m *big.Int) *big.Int {
	if exp.Sign() >= 0 {
		return new(big.Int).Exp(base, exp, m)
	}
	inv := new(big.Int).ModInverse(base, m)
	if inv == nil {
		return nil
	}
	posExp := new(big.Int).Neg(exp)
	return new(big.Int).Exp(inv, posExp, m)
}

// InvMod returns the multiplicative inverse of a modulo m, in the range
// [0, m), or an error if a is not invertible modulo m.
func InvMod(a, m *big.Int) (*big.Int, error) {
	inv := new(big.Int).ModInverse(a, m)
	if inv == nil {
		return nil, ErrNotInvertible
	}
	return inv, nil
}

// Mod returns a reduced into the canonical range [0, m).
func Mod(a, m *big.Int) *big.Int {
	return new(big.Int).Mod(a, m)
}

// GCD returns the greatest common divisor of a and b as a fresh big.Int.
func GCD(a, b *big.Int) *big.Int {
	return new(big.Int).GCD(nil, nil, new(big.Int).Abs(a), new(big.Int).Abs(b))
}

// Coprime reports whether a and b share no common factor greater than 1.
func Coprime(a, b *big.Int) bool {
	return GCD(a, b).Cmp(bigOne) == 0
}

// IsPrime reports whether n is a prime, using a strong probabilistic test
// (Baillie-PSW plus Miller-Rabin rounds) that is deterministic in practice
// for the sizes used here.
func IsPrime(n *big.Int) bool {
	if n == nil || n.Sign() <= 0 {
		return false
	}
	return n.ProbablyPrime(20)
}

// MillerRabin reports whether n is probably prime using the given number of
// Miller-Rabin rounds with witnesses drawn from rng. A nil rng disables the
// random rounds and falls back to the standard library test.
func MillerRabin(n *big.Int, rounds int, rng *rand.Rand) bool {
	if n == nil || n.Sign() <= 0 {
		return false
	}
	if n.Cmp(bigTwo) < 0 {
		return false
	}
	if n.Bit(0) == 0 {
		return n.Cmp(bigTwo) == 0
	}
	if rng == nil {
		return n.ProbablyPrime(rounds)
	}
	// n-1 = 2^s * d with d odd.
	nm1 := new(big.Int).Sub(n, bigOne)
	d := new(big.Int).Set(nm1)
	s := 0
	for d.Bit(0) == 0 {
		d.Rsh(d, 1)
		s++
	}
	nm3 := new(big.Int).Sub(n, big.NewInt(3))
	for i := 0; i < rounds; i++ {
		// random witness a in [2, n-2]
		a := new(big.Int).Rand(rng, nm3)
		a.Add(a, bigTwo)
		x := new(big.Int).Exp(a, d, n)
		if x.Cmp(bigOne) == 0 || x.Cmp(nm1) == 0 {
			continue
		}
		composite := true
		for j := 0; j < s-1; j++ {
			x.Mul(x, x).Mod(x, n)
			if x.Cmp(nm1) == 0 {
				composite = false
				break
			}
		}
		if composite {
			return false
		}
	}
	return true
}

// NextPrime returns the smallest prime strictly greater than n.
func NextPrime(n *big.Int) *big.Int {
	c := new(big.Int).Add(n, bigOne)
	if c.Cmp(bigTwo) < 0 {
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

// PrevPrime returns the largest prime strictly less than n, or nil if there is
// none (n <= 2).
func PrevPrime(n *big.Int) *big.Int {
	if n.Cmp(big.NewInt(3)) <= 0 {
		if n.Cmp(big.NewInt(3)) == 0 {
			return big.NewInt(2)
		}
		return nil
	}
	c := new(big.Int).Sub(n, bigOne)
	if c.Bit(0) == 0 {
		c.Sub(c, bigOne)
	}
	for c.Cmp(bigTwo) > 0 && !c.ProbablyPrime(20) {
		c.Sub(c, bigTwo)
	}
	if c.ProbablyPrime(20) {
		return c
	}
	return big.NewInt(2)
}

// LegendreSymbol returns the Legendre symbol (a/p) for an odd prime p, one of
// -1, 0, or +1. It does not verify that p is prime.
func LegendreSymbol(a, p *big.Int) int {
	am := new(big.Int).Mod(a, p)
	if am.Sign() == 0 {
		return 0
	}
	e := new(big.Int).Rsh(new(big.Int).Sub(p, bigOne), 1) // (p-1)/2
	r := new(big.Int).Exp(am, e, p)
	if r.Cmp(bigOne) == 0 {
		return 1
	}
	return -1
}

// JacobiSymbol returns the Jacobi symbol (a/n) for an odd positive n, one of
// -1, 0, or +1.
func JacobiSymbol(a, n *big.Int) int {
	if n.Sign() <= 0 || n.Bit(0) == 0 {
		return 0
	}
	x := new(big.Int).Mod(a, n)
	y := new(big.Int).Set(n)
	result := 1
	for x.Sign() != 0 {
		for x.Bit(0) == 0 {
			x.Rsh(x, 1)
			m8 := new(big.Int).And(y, big.NewInt(7)).Int64()
			if m8 == 3 || m8 == 5 {
				result = -result
			}
		}
		x, y = y, x
		if new(big.Int).And(x, big.NewInt(3)).Int64() == 3 &&
			new(big.Int).And(y, big.NewInt(3)).Int64() == 3 {
			result = -result
		}
		x.Mod(x, y)
	}
	if y.Cmp(bigOne) == 0 {
		return result
	}
	return 0
}

// IsQuadraticResidue reports whether a is a non-zero quadratic residue modulo
// the odd prime p.
func IsQuadraticResidue(a, p *big.Int) bool {
	return LegendreSymbol(a, p) == 1
}

// SqrtModP returns a square root of a modulo the odd prime p using the
// Tonelli-Shanks algorithm. rng supplies the randomness used to find a
// quadratic non-residue; pass a seeded *rand.Rand for reproducibility. It
// returns ErrNoRoot when a is a non-residue.
func SqrtModP(a, p *big.Int, rng *rand.Rand) (*big.Int, error) {
	am := new(big.Int).Mod(a, p)
	if am.Sign() == 0 {
		return big.NewInt(0), nil
	}
	if p.Cmp(bigTwo) == 0 {
		return new(big.Int).Set(am), nil
	}
	if LegendreSymbol(am, p) != 1 {
		return nil, ErrNoRoot
	}
	// p ≡ 3 mod 4 fast path.
	if new(big.Int).And(p, big.NewInt(3)).Int64() == 3 {
		e := new(big.Int).Rsh(new(big.Int).Add(p, bigOne), 2) // (p+1)/4
		return new(big.Int).Exp(am, e, p), nil
	}
	// Write p-1 = q * 2^s with q odd.
	q := new(big.Int).Sub(p, bigOne)
	s := 0
	for q.Bit(0) == 0 {
		q.Rsh(q, 1)
		s++
	}
	// Find a non-residue z. With rng == nil, scan deterministically.
	z := big.NewInt(2)
	if rng == nil {
		for LegendreSymbol(z, p) != -1 {
			z.Add(z, bigOne)
		}
	} else {
		pm2 := new(big.Int).Sub(p, bigTwo)
		for LegendreSymbol(z, p) != -1 {
			z.Rand(rng, pm2)
			z.Add(z, bigTwo)
		}
	}
	m := s
	c := new(big.Int).Exp(z, q, p)
	t := new(big.Int).Exp(am, q, p)
	r := new(big.Int).Exp(am, new(big.Int).Rsh(new(big.Int).Add(q, bigOne), 1), p)
	for t.Cmp(bigOne) != 0 {
		// Find least i, 0<i<m, with t^(2^i) = 1.
		i := 0
		tt := new(big.Int).Set(t)
		for i = 1; i < m; i++ {
			tt.Mul(tt, tt).Mod(tt, p)
			if tt.Cmp(bigOne) == 0 {
				break
			}
		}
		b := new(big.Int).Set(c)
		for j := 0; j < m-i-1; j++ {
			b.Mul(b, b).Mod(b, p)
		}
		m = i
		c.Mul(b, b).Mod(c, p)
		t.Mul(t, c).Mod(t, p)
		r.Mul(r, b).Mod(r, p)
	}
	return r, nil
}

// abs returns |x| as a fresh big.Int.
func absInt(x *big.Int) *big.Int { return new(big.Int).Abs(x) }

// minInt returns the smaller of a and b.
func minInt(a, b int) int {
	if a < b {
		return a
	}
	return b
}

// maxInt returns the larger of a and b.
func maxInt(a, b int) int {
	if a > b {
		return a
	}
	return b
}

// ValuationInt64 returns the p-adic valuation of the machine integer n for the
// machine prime p (p >= 2). It returns -1 for n == 0.
func ValuationInt64(p, n int64) int {
	if n == 0 {
		return -1
	}
	if n < 0 {
		n = -n
	}
	v := 0
	for n%p == 0 {
		n /= p
		v++
	}
	return v
}

// Factorial returns n! as a big.Int for n >= 0.
func Factorial(n int) *big.Int {
	r := big.NewInt(1)
	for i := int64(2); i <= int64(n); i++ {
		r.Mul(r, big.NewInt(i))
	}
	return r
}

// ValuationFactorial returns the p-adic valuation of n! using Legendre's
// formula, v_p(n!) = sum_{i>=1} floor(n / p^i) = (n - s_p(n)) / (p - 1).
func ValuationFactorial(p *big.Int, n int) int {
	if n <= 0 {
		return 0
	}
	total := 0
	pk := new(big.Int).Set(p)
	nb := big.NewInt(int64(n))
	for pk.Cmp(nb) <= 0 {
		q := new(big.Int).Div(nb, pk)
		total += int(q.Int64())
		pk.Mul(pk, p)
	}
	return total
}

// Binomial returns the binomial coefficient C(n, k) as a big.Int.
func Binomial(n, k int) *big.Int {
	if k < 0 || k > n {
		return big.NewInt(0)
	}
	return new(big.Int).Binomial(int64(n), int64(k))
}

// MultiplicativeOrderModP returns the multiplicative order of a modulo the
// prime p (the least e > 0 with a^e == 1 mod p), or 0 if a is divisible by p.
func MultiplicativeOrderModP(a, p *big.Int) int {
	am := new(big.Int).Mod(a, p)
	if am.Sign() == 0 {
		return 0
	}
	cur := new(big.Int).Set(am)
	for e := 1; ; e++ {
		if cur.Cmp(bigOne) == 0 {
			return e
		}
		cur.Mul(cur, am).Mod(cur, p)
		if e > int(new(big.Int).Sub(p, bigOne).Int64())+1 {
			return 0
		}
	}
}

// CRTPair returns the unique residue x modulo m1*m2 with x == r1 mod m1 and
// x == r2 mod m2, for coprime moduli m1, m2. It returns an error if the moduli
// are not coprime.
func CRTPair(r1, m1, r2, m2 *big.Int) (*big.Int, *big.Int, error) {
	if !Coprime(m1, m2) {
		return nil, nil, ErrNotInvertible
	}
	m := new(big.Int).Mul(m1, m2)
	m1inv := new(big.Int).ModInverse(new(big.Int).Mod(m1, m2), m2)
	// x = r1 + m1 * ((r2 - r1) * m1inv mod m2)
	diff := new(big.Int).Sub(r2, r1)
	t := new(big.Int).Mul(diff, m1inv)
	t.Mod(t, m2)
	x := new(big.Int).Add(r1, new(big.Int).Mul(m1, t))
	x.Mod(x, m)
	return x, m, nil
}

// IsPrimePower reports whether n is a prime power p^k (k >= 1); if so it also
// returns p and k. Otherwise ok is false.
func IsPrimePower(n *big.Int) (p *big.Int, k int, ok bool) {
	if n.Cmp(bigTwo) < 0 {
		return nil, 0, false
	}
	// find the smallest prime factor by trial division
	d := big.NewInt(2)
	r := new(big.Int)
	for new(big.Int).Mul(d, d).Cmp(n) <= 0 {
		if new(big.Int).Mod(n, d).Sign() == 0 {
			m := new(big.Int).Set(n)
			e := 0
			for {
				q := new(big.Int)
				q.QuoRem(m, d, r)
				if r.Sign() != 0 {
					break
				}
				m.Set(q)
				e++
			}
			if m.Cmp(bigOne) == 0 {
				return new(big.Int).Set(d), e, true
			}
			return nil, 0, false
		}
		d.Add(d, bigOne)
	}
	// n itself is prime
	return new(big.Int).Set(n), 1, true
}
