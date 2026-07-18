package crypto

import (
	"math/big"
	"sort"
)

// Factor represents a prime power in an integer's factorization: Prime raised
// to the given positive Exponent.
type Factor struct {
	Prime    *big.Int
	Exponent int
}

// PollardRho returns a non-trivial factor of the composite integer n using
// Pollard's rho algorithm with the iteration x -> x^2 + c. It returns 2 for
// even n. If n is prime a copy of n itself is returned, since it has no proper
// factor. It panics if n <= 1. The polynomial constant c is advanced
// deterministically on failure, so the routine needs no randomness and always
// returns the same factor for the same input.
func PollardRho(n *big.Int) *big.Int {
	if n.Cmp(cryptoOne) <= 0 {
		panic("crypto: PollardRho requires n > 1")
	}
	if n.Bit(0) == 0 {
		return big.NewInt(2)
	}
	if IsPrime(n) {
		return new(big.Int).Set(n)
	}
	x := new(big.Int)
	y := new(big.Int)
	d := new(big.Int)
	diff := new(big.Int)
	for c := int64(1); ; c++ {
		bc := big.NewInt(c)
		x.SetInt64(2)
		y.SetInt64(2)
		d.SetInt64(1)
		for d.Cmp(cryptoOne) == 0 {
			cryptoRhoStep(x, bc, n)
			cryptoRhoStep(y, bc, n)
			cryptoRhoStep(y, bc, n)
			diff.Sub(x, y)
			diff.Abs(diff)
			if diff.Sign() == 0 {
				break
			}
			d.GCD(nil, nil, diff, n)
		}
		if d.Cmp(cryptoOne) != 0 && d.Cmp(n) != 0 {
			return new(big.Int).Set(d)
		}
		// Otherwise retry with the next constant c.
	}
}

// cryptoRhoStep advances x in place to (x^2 + c) mod n.
func cryptoRhoStep(x, c, n *big.Int) {
	x.Mul(x, x)
	x.Add(x, c)
	x.Mod(x, n)
}

// PollardRhoBrent returns a non-trivial factor of the composite integer n using
// Brent's improvement to Pollard's rho, which reduces the number of gcd
// computations. It behaves like PollardRho with respect to even and prime
// inputs and is deterministic. It panics if n <= 1.
func PollardRhoBrent(n *big.Int) *big.Int {
	if n.Cmp(cryptoOne) <= 0 {
		panic("crypto: PollardRhoBrent requires n > 1")
	}
	if n.Bit(0) == 0 {
		return big.NewInt(2)
	}
	if IsPrime(n) {
		return new(big.Int).Set(n)
	}
	for c := int64(1); ; c++ {
		bc := big.NewInt(c)
		if g := cryptoBrentTry(n, bc); g != nil {
			return g
		}
	}
}

// cryptoBrentTry runs one attempt of Brent's cycle-finding rho with constant c,
// returning a non-trivial factor of n or nil if this attempt degenerated.
func cryptoBrentTry(n, c *big.Int) *big.Int {
	y := big.NewInt(2)
	m := big.NewInt(128)
	g := big.NewInt(1)
	r := big.NewInt(1)
	q := big.NewInt(1)
	var x, ys *big.Int
	tmp := new(big.Int)
	for g.Cmp(cryptoOne) == 0 {
		x = new(big.Int).Set(y)
		for i := new(big.Int).Set(r); i.Sign() > 0; i.Sub(i, cryptoOne) {
			cryptoRhoStep(y, c, n)
		}
		k := big.NewInt(0)
		for k.Cmp(r) < 0 && g.Cmp(cryptoOne) == 0 {
			ys = new(big.Int).Set(y)
			bound := new(big.Int).Sub(r, k)
			if bound.Cmp(m) > 0 {
				bound.Set(m)
			}
			for i := new(big.Int).Set(bound); i.Sign() > 0; i.Sub(i, cryptoOne) {
				cryptoRhoStep(y, c, n)
				tmp.Sub(x, y)
				tmp.Abs(tmp)
				q.Mul(q, tmp)
				q.Mod(q, n)
			}
			g.GCD(nil, nil, q, n)
			k.Add(k, m)
		}
		r.Lsh(r, 1)
	}
	if g.Cmp(n) == 0 {
		// Back-track to find the factor one step at a time.
		for {
			cryptoRhoStep(ys, c, n)
			tmp.Sub(x, ys)
			tmp.Abs(tmp)
			if tmp.Sign() == 0 {
				return nil
			}
			g.GCD(nil, nil, tmp, n)
			if g.Cmp(cryptoOne) > 0 {
				break
			}
		}
	}
	if g.Cmp(cryptoOne) != 0 && g.Cmp(n) != 0 {
		return g
	}
	return nil
}

// PrimeFactors returns the distinct-with-multiplicity prime factors of n in
// ascending order (each prime appears as many times as it divides n). For
// example PrimeFactors(360) yields [2 2 2 3 3 5]. It panics if n < 1;
// PrimeFactors(1) is the empty slice.
func PrimeFactors(n *big.Int) []*big.Int {
	if n.Sign() < 1 {
		panic("crypto: PrimeFactors requires n >= 1")
	}
	var factors []*big.Int
	m := new(big.Int).Set(n)
	// Strip small primes first for speed.
	for _, p := range cryptoSmallPrimes {
		bp := big.NewInt(p)
		for new(big.Int).Mod(m, bp).Sign() == 0 {
			factors = append(factors, big.NewInt(p))
			m.Div(m, bp)
		}
	}
	cryptoFactorRec(m, &factors)
	sort.Slice(factors, func(i, j int) bool {
		return factors[i].Cmp(factors[j]) < 0
	})
	return factors
}

// cryptoFactorRec appends the prime factors of m (with multiplicity) to *out,
// recursing through Pollard's rho for composite parts.
func cryptoFactorRec(m *big.Int, out *[]*big.Int) {
	if m.Cmp(cryptoOne) <= 0 {
		return
	}
	if IsPrime(m) {
		*out = append(*out, new(big.Int).Set(m))
		return
	}
	d := PollardRho(m)
	cryptoFactorRec(d, out)
	cryptoFactorRec(new(big.Int).Div(m, d), out)
}

// Factorization returns the prime factorization of n as a slice of Factor
// (prime, exponent) pairs sorted by ascending prime. It panics if n < 1;
// Factorization(1) is the empty slice.
func Factorization(n *big.Int) []Factor {
	primes := PrimeFactors(n)
	var result []Factor
	for _, p := range primes {
		if len(result) > 0 && result[len(result)-1].Prime.Cmp(p) == 0 {
			result[len(result)-1].Exponent++
		} else {
			result = append(result, Factor{Prime: new(big.Int).Set(p), Exponent: 1})
		}
	}
	return result
}

// SmallestPrimeFactor returns the least prime factor of n. It panics if n < 2.
func SmallestPrimeFactor(n *big.Int) *big.Int {
	if n.Cmp(cryptoTwo) < 0 {
		panic("crypto: SmallestPrimeFactor requires n >= 2")
	}
	primes := PrimeFactors(n)
	return new(big.Int).Set(primes[0])
}

// EulerTotient returns Euler's totient φ(n): the count of integers in [1, n]
// that are coprime to n. It is computed from the prime factorization as
// n · ∏ (1 - 1/p). It panics if n < 1; φ(1) = 1.
func EulerTotient(n *big.Int) *big.Int {
	if n.Sign() < 1 {
		panic("crypto: EulerTotient requires n >= 1")
	}
	if n.Cmp(cryptoOne) == 0 {
		return big.NewInt(1)
	}
	result := new(big.Int)
	for _, f := range Factorization(n) {
		// p^(e-1) * (p-1)
		term := new(big.Int).Exp(f.Prime, big.NewInt(int64(f.Exponent-1)), nil)
		term.Mul(term, new(big.Int).Sub(f.Prime, cryptoOne))
		if result.Sign() == 0 {
			result.Set(term)
		} else {
			result.Mul(result, term)
		}
	}
	return result
}

// CarmichaelLambda returns the Carmichael function λ(n): the smallest positive
// integer m such that a^m ≡ 1 (mod n) for every a coprime to n. It is the least
// common multiple of λ over each prime-power factor. It panics if n < 1; λ(1) =
// 1.
func CarmichaelLambda(n *big.Int) *big.Int {
	if n.Sign() < 1 {
		panic("crypto: CarmichaelLambda requires n >= 1")
	}
	if n.Cmp(cryptoOne) == 0 {
		return big.NewInt(1)
	}
	result := big.NewInt(1)
	for _, f := range Factorization(n) {
		var lambda *big.Int
		if f.Prime.Cmp(cryptoTwo) == 0 && f.Exponent >= 3 {
			// λ(2^e) = 2^(e-2) for e >= 3.
			lambda = new(big.Int).Lsh(cryptoOne, uint(f.Exponent-2))
		} else {
			// λ(p^e) = p^(e-1) * (p-1) for odd p (and 2, 4).
			lambda = new(big.Int).Exp(f.Prime, big.NewInt(int64(f.Exponent-1)), nil)
			lambda.Mul(lambda, new(big.Int).Sub(f.Prime, cryptoOne))
		}
		result = LCM(result, lambda)
	}
	return result
}
