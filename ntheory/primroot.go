package ntheory

import "sort"

// ntheoryLCMU64 returns the least common multiple of two positive unsigned
// integers. It divides before multiplying so intermediate values stay as small
// as possible. ntheoryLCMU64(0, x) and ntheoryLCMU64(x, 0) return 0.
func ntheoryLCMU64(a, b uint64) uint64 {
	if a == 0 || b == 0 {
		return 0
	}
	return a / ntheoryGCDU64(a, b) * b
}

// ntheoryPhiU64 returns Euler's totient φ(n) computed from the prime
// factorization via [FactorizeU64], using the product form
// n · ∏(1 − 1/p) evaluated as n/p·(p−1) to avoid overflow and fractions.
// ntheoryPhiU64(0) == 0 and ntheoryPhiU64(1) == 1.
func ntheoryPhiU64(n uint64) uint64 {
	if n == 0 {
		return 0
	}
	if n == 1 {
		return 1
	}
	result := n
	for p := range FactorizeU64(n) {
		result = result / p * (p - 1)
	}
	return result
}

// ntheoryHasPrimitiveRootU64 reports whether the multiplicative group modulo m
// is cyclic, i.e. whether a primitive root modulo m exists. This holds exactly
// for m in {1, 2, 4}, for m = p^k with p an odd prime, and for m = 2·p^k.
func ntheoryHasPrimitiveRootU64(m uint64) bool {
	if m == 1 || m == 2 || m == 4 {
		return true
	}
	if m == 0 {
		return false
	}
	// Strip a single factor of 2; what remains must be an odd prime power.
	r := m
	if r%2 == 0 {
		r /= 2
		if r%2 == 0 {
			return false // divisible by 4 and > 4
		}
	}
	if r == 1 {
		return false // guards leftover even values
	}
	return len(FactorizeU64(r)) == 1
}

// Carmichael returns λ(n), the Carmichael function: the exponent of the
// multiplicative group modulo n, i.e. the smallest m > 0 with a^m ≡ 1 (mod n)
// for every a coprime to n.
//
// It is computed from the prime factorization ([FactorizeU64]) as the least
// common multiple of λ over the prime powers of n, where λ(2) = 1, λ(4) = 2,
// λ(2^k) = 2^(k−2) for k ≥ 3, and λ(p^k) = p^(k−1)(p−1) for odd primes p.
// By convention λ(1) = 1 (and λ(0) = 0). n is unsigned, so its sign is
// irrelevant.
func Carmichael(n uint64) uint64 {
	if n == 0 {
		return 0
	}
	if n == 1 {
		return 1
	}
	result := uint64(1)
	for p, k := range FactorizeU64(n) {
		var lambdaPK uint64
		if p == 2 {
			switch {
			case k == 1:
				lambdaPK = 1
			case k == 2:
				lambdaPK = 2
			default:
				lambdaPK = uint64(1) << (k - 2)
			}
		} else {
			// p^(k-1) * (p-1)
			lambdaPK = p - 1
			for i := 0; i < k-1; i++ {
				lambdaPK *= p
			}
		}
		result = ntheoryLCMU64(result, lambdaPK)
	}
	return result
}

// MultiplicativeOrder returns the multiplicative order of a modulo m — the
// smallest k > 0 with a^k ≡ 1 (mod m) — together with ok == true. It returns
// (0, false) when gcd(a, m) ≠ 1, since a then has no order. m must be positive.
//
// It is computed the fast way: the order necessarily divides λ(m), so it starts
// from ord = [Carmichael](m) and, for each prime q dividing λ(m), divides q out
// of ord as long as a^(ord/q) ≡ 1 (mod m) still holds ([ModPowU64]). This tests
// only the prime-quotient candidates rather than enumerating every divisor, so
// it is asymptotically far cheaper than the divisor-scanning [Order].
func MultiplicativeOrder(a, m uint64) (uint64, bool) {
	if m == 0 {
		panic("ntheory: MultiplicativeOrder requires a positive modulus")
	}
	if m == 1 {
		return 1, true
	}
	a %= m
	if ntheoryGCDU64(a, m) != 1 {
		return 0, false
	}
	ord := Carmichael(m)
	for q := range FactorizeU64(ord) {
		for ord%q == 0 && ModPowU64(a, ord/q, m) == 1 {
			ord /= q
		}
	}
	return ord, true
}

// PrimitiveRoot returns the smallest primitive root modulo m, together with
// ok == true. It returns (0, false) when no primitive root exists; primitive
// roots exist exactly for m in {1, 2, 4}, m = p^k with p an odd prime, and
// m = 2·p^k. By convention PrimitiveRoot(1) == (0, true). m must be positive.
//
// Existence is decided structurally, then candidates g are tested with the fast
// prime-quotient criterion of [IsPrimitiveRoot]: g is a primitive root iff
// g^(λ(m)/q) ≠ 1 (mod m) for every prime q dividing λ(m). Because a primitive
// root exists only when the group is cyclic (λ(m) = φ(m)), this reuses
// [FactorizeU64] of φ(m) and never scans all divisors.
func PrimitiveRoot(m uint64) (uint64, bool) {
	if m == 0 {
		panic("ntheory: PrimitiveRoot requires a positive modulus")
	}
	if m == 1 {
		return 0, true
	}
	if !ntheoryHasPrimitiveRootU64(m) {
		return 0, false
	}
	for g := uint64(1); g < m; g++ {
		if IsPrimitiveRoot(g, m) {
			return g, true
		}
	}
	return 0, false
}

// IsPrimitiveRoot reports whether g is a primitive root modulo m, i.e. whether
// the multiplicative order of g modulo m equals φ(m). g must be coprime to m for
// the result to be true. m must be positive; IsPrimitiveRoot(g, 1) == true.
//
// Rather than compute the order outright, it applies the fast criterion: with
// gcd(g, m) = 1, g is a primitive root iff g^(φ(m)/q) ≠ 1 (mod m) for every
// prime q dividing φ(m) ([FactorizeU64] and [ModPowU64]). This also correctly
// reports false when the group is not cyclic, since then no g attains order
// φ(m).
func IsPrimitiveRoot(g, m uint64) bool {
	if m == 0 {
		panic("ntheory: IsPrimitiveRoot requires a positive modulus")
	}
	if m == 1 {
		return true
	}
	g %= m
	if ntheoryGCDU64(g, m) != 1 {
		return false
	}
	phi := ntheoryPhiU64(m)
	for q := range FactorizeU64(phi) {
		if ModPowU64(g, phi/q, m) == 1 {
			return false
		}
	}
	return true
}

// PrimitiveRoots returns all primitive roots modulo m in ascending order, or nil
// when none exist. When the group is cyclic there are exactly φ(φ(m)) of them.
//
// It finds the smallest primitive root g with [PrimitiveRoot] and generates the
// rest as the powers g^k for which gcd(k, φ(m)) = 1, k ranging over [1, φ(m)];
// these powers are precisely the generators of the cyclic group. The results are
// sorted before being returned so the slice is deterministic and ascending.
func PrimitiveRoots(m uint64) []uint64 {
	g, ok := PrimitiveRoot(m)
	if !ok {
		return nil
	}
	phi := ntheoryPhiU64(m)
	roots := make([]uint64, 0)
	for k := uint64(1); k <= phi; k++ {
		if ntheoryGCDU64(k, phi) == 1 {
			roots = append(roots, ModPowU64(g, k, m))
		}
	}
	sort.Slice(roots, func(i, j int) bool { return roots[i] < roots[j] })
	return roots
}
