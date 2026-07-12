package ntheory

import "math/big"

// ModPow returns base**exp mod m computed with math/big internally so no
// intermediate value overflows. The result is normalized to the range [0, m).
//
// m must be positive. A negative exp is supported when base is invertible
// modulo m: the inverse of base is raised to |exp|; if base is not invertible
// modulo m, ModPow panics.
func ModPow(base, exp, m int64) int64 {
	if m <= 0 {
		panic("ntheory: ModPow requires a positive modulus")
	}
	if m == 1 {
		return 0
	}
	r := new(big.Int).Exp(big.NewInt(base), big.NewInt(exp), big.NewInt(m))
	if r == nil {
		panic("ntheory: ModPow with negative exponent requires base coprime to modulus")
	}
	return r.Int64()
}

// ModInverse returns the modular multiplicative inverse of a modulo m: the
// unique x in [0, m) with a*x ≡ 1 (mod m), together with ok == true. If a and m
// are not coprime no inverse exists and it returns (0, false). m must be
// positive.
func ModInverse(a, m int64) (inv int64, ok bool) {
	if m <= 0 {
		panic("ntheory: ModInverse requires a positive modulus")
	}
	if m == 1 {
		return 0, true
	}
	g, x, _ := ExtendedGCD(((a%m)+m)%m, m)
	if g != 1 {
		return 0, false
	}
	return ((x % m) + m) % m, true
}

// CRT solves a system of simultaneous congruences via the Chinese Remainder
// Theorem. Given residues r and moduli n of equal length, it returns the unique
// x in [0, N) with x ≡ r[i] (mod n[i]) for all i, where N is the least common
// multiple of the moduli, together with N and ok == true.
//
// The moduli need not be pairwise coprime; if the congruences are inconsistent
// (no solution exists) it returns (0, 0, false). Every modulus must be
// positive, and r and n must have the same non-zero length.
func CRT(r, n []int64) (x, mod int64, ok bool) {
	if len(r) != len(n) || len(r) == 0 {
		panic("ntheory: CRT requires equal, non-empty residue and modulus slices")
	}
	x = ((r[0] % n[0]) + n[0]) % n[0]
	mod = n[0]
	if mod <= 0 {
		panic("ntheory: CRT requires positive moduli")
	}
	for i := 1; i < len(n); i++ {
		ni := n[i]
		if ni <= 0 {
			panic("ntheory: CRT requires positive moduli")
		}
		ri := ((r[i] % ni) + ni) % ni
		g, p, _ := ExtendedGCD(mod, ni)
		diff := ri - x
		if diff%g != 0 {
			return 0, 0, false
		}
		lcm := mod / g * ni
		// x_new = x + mod * (diff/g * p mod (ni/g)), all done in big to be safe.
		bx := big.NewInt(x)
		bmod := big.NewInt(mod)
		step := big.NewInt(diff / g)
		step.Mul(step, big.NewInt(p))
		step.Mod(step, big.NewInt(ni/g))
		step.Mul(step, bmod)
		bx.Add(bx, step)
		blcm := big.NewInt(lcm)
		bx.Mod(bx, blcm)
		if bx.Sign() < 0 {
			bx.Add(bx, blcm)
		}
		x = bx.Int64()
		mod = lcm
	}
	return x, mod, true
}

// LegendreSymbol returns the Legendre symbol (a/p) for an odd prime p:
//
//	 0 if a ≡ 0 (mod p),
//	+1 if a is a non-zero quadratic residue modulo p,
//	-1 if a is a non-residue.
//
// p must be an odd prime; the behavior is otherwise undefined.
func LegendreSymbol(a, p int64) int {
	if p <= 2 {
		panic("ntheory: LegendreSymbol requires an odd prime p")
	}
	a = ((a % p) + p) % p
	if a == 0 {
		return 0
	}
	r := ModPow(a, (p-1)/2, p)
	if r == p-1 {
		return -1
	}
	return int(r) // 1
}

// JacobiSymbol returns the Jacobi symbol (a/n) for an odd positive integer n.
// It generalizes the Legendre symbol to composite n; note that a result of +1
// does not by itself prove a is a quadratic residue when n is composite. n must
// be odd and positive.
func JacobiSymbol(a, n int64) int {
	if n <= 0 || n%2 == 0 {
		panic("ntheory: JacobiSymbol requires an odd positive n")
	}
	a = ((a % n) + n) % n
	result := 1
	for a != 0 {
		for a%2 == 0 {
			a /= 2
			if r := n % 8; r == 3 || r == 5 {
				result = -result
			}
		}
		a, n = n, a
		if a%4 == 3 && n%4 == 3 {
			result = -result
		}
		a %= n
	}
	if n == 1 {
		return result
	}
	return 0
}

// IsQuadraticResidue reports whether a is a quadratic residue modulo the prime
// p, i.e. whether some x satisfies x**2 ≡ a (mod p). Zero counts as a residue
// (x == 0). p must be prime.
func IsQuadraticResidue(a, p int64) bool {
	if p <= 1 {
		panic("ntheory: IsQuadraticResidue requires a prime p")
	}
	a = ((a % p) + p) % p
	if p == 2 {
		return true // 0 and 1 are both residues mod 2
	}
	if a == 0 {
		return true
	}
	return LegendreSymbol(a, p) == 1
}

// DiscreteLog returns the smallest non-negative x with g**x ≡ h (mod m), found
// via the baby-step giant-step algorithm, together with ok == true. If no such
// x below m exists it returns (0, false). m must be positive and g should be a
// unit modulo m for a solution to be guaranteed searchable; the search space is
// bounded by m.
func DiscreteLog(g, h, m int64) (x int64, ok bool) {
	if m <= 0 {
		panic("ntheory: DiscreteLog requires a positive modulus")
	}
	g = ((g % m) + m) % m
	h = ((h % m) + m) % m
	if m == 1 {
		return 0, true
	}
	// n = ceil(sqrt(m)).
	n := int64(1)
	for n*n < m {
		n++
	}
	// Baby steps: table of g^j -> j for j in [0, n).
	table := make(map[int64]int64, n)
	cur := int64(1)
	for j := int64(0); j < n; j++ {
		if _, seen := table[cur]; !seen {
			table[cur] = j
		}
		cur = mulMod(cur, g, m)
	}
	// factor = g^(-n) mod m.
	gn := ModPow(g, n, m)
	factor, inv := ModInverse(gn, m)
	if !inv {
		// g not invertible; fall back to a bounded linear scan.
		val := int64(1)
		for i := int64(0); i < m; i++ {
			if val == h {
				return i, true
			}
			val = mulMod(val, g, m)
		}
		return 0, false
	}
	gamma := h
	for i := int64(0); i <= n; i++ {
		if j, found := table[gamma]; found {
			return i*n + j, true
		}
		gamma = mulMod(gamma, factor, m)
	}
	return 0, false
}

// Order returns the multiplicative order of a modulo m: the smallest positive
// integer k with a**k ≡ 1 (mod m), together with ok == true. If a and m are not
// coprime the order does not exist and it returns (0, false). m must be
// positive.
func Order(a, m int64) (ord int64, ok bool) {
	if m <= 0 {
		panic("ntheory: Order requires a positive modulus")
	}
	if m == 1 {
		return 1, true
	}
	a = ((a % m) + m) % m
	if GCD(a, m) != 1 {
		return 0, false
	}
	phi := EulerPhi(m)
	// The order divides phi(m); test divisors of phi in ascending order.
	divs := Divisors(phi)
	for _, d := range divs {
		if ModPow(a, d, m) == 1 {
			return d, true
		}
	}
	return phi, true
}
