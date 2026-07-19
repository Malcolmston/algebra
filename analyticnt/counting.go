package analyticnt

import "math"

// LegendrePhi returns Legendre's φ(x, a): the count of integers in [1, x] not
// divisible by any of the first a primes. It is evaluated by the standard
// recursion φ(x, a) = φ(x, a−1) − φ(⌊x/p_a⌋, a−1) with φ(x, 0) = ⌊x⌋.
func LegendrePhi(x int64, a int) int64 {
	if x < 0 {
		return 0
	}
	primes := firstPrimes(a)
	return legendrePhi(x, a, primes)
}

// firstPrimes returns the first a primes.
func firstPrimes(a int) []int64 {
	if a <= 0 {
		return nil
	}
	// Estimate an upper bound for the a-th prime.
	limit := 20
	if a >= 6 {
		fa := float64(a)
		limit = int(fa*(math.Log(fa)+math.Log(math.Log(fa)))) + 3
	}
	var ps []int64
	for {
		ps = PrimesUpTo(limit)
		if len(ps) >= a {
			return ps[:a]
		}
		limit *= 2
	}
}

// legendrePhi is the memo-free recursion driven by a precomputed prime prefix.
func legendrePhi(x int64, a int, primes []int64) int64 {
	if a == 0 {
		return x
	}
	if x == 0 {
		return 0
	}
	return legendrePhi(x, a-1, primes) - legendrePhi(x/primes[a-1], a-1, primes)
}

// PrimePiLegendre returns π(x) via Legendre's formula
// π(x) = φ(x, a) + a − 1, where a = π(⌊√x⌋). It is exact.
func PrimePiLegendre(x int64) int64 {
	if x < 2 {
		return 0
	}
	root := int64(math.Sqrt(float64(x)))
	for (root+1)*(root+1) <= x {
		root++
	}
	for root*root > x {
		root--
	}
	a := int(PrimePi(root))
	primes := firstPrimes(a)
	return legendrePhi(x, a, primes) + int64(a) - 1
}

// PrimePiMeissel returns π(x) via the Meissel formula, which reduces the work
// relative to Legendre's method. With a = π(x^{1/3}) and b = π(x^{1/2}),
// π(x) = φ(x, a) + (b + a − 2)(b − a + 1)/2 − Σ_{i=a+1}^{b} π(⌊x/p_i⌋). It is
// exact and agrees with PrimePi.
func PrimePiMeissel(x int64) int64 {
	if x < 2 {
		return 0
	}
	if x < 8 {
		return PrimePi(x)
	}
	cbrt := int64(math.Cbrt(float64(x)))
	for (cbrt+1)*(cbrt+1)*(cbrt+1) <= x {
		cbrt++
	}
	for cbrt*cbrt*cbrt > x {
		cbrt--
	}
	sqrt := int64(math.Sqrt(float64(x)))
	for (sqrt+1)*(sqrt+1) <= x {
		sqrt++
	}
	for sqrt*sqrt > x {
		sqrt--
	}
	a := int(PrimePi(cbrt))
	b := int(PrimePi(sqrt))
	primes := firstPrimes(b)
	phi := legendrePhi(x, a, primes)
	// Combinatorial correction.
	corr := int64(b+a-2) * int64(b-a+1) / 2
	var sum int64
	for i := a + 1; i <= b; i++ {
		sum += PrimePi(x / primes[i-1])
	}
	return phi + corr - sum
}

// PrimePiCheck cross-checks all three exact counters (sieve, Legendre, Meissel)
// at x and reports whether they agree, returning the common value and the
// agreement flag.
func PrimePiCheck(x int64) (value int64, agree bool) {
	a := PrimePi(x)
	b := PrimePiLegendre(x)
	c := PrimePiMeissel(x)
	return a, a == b && b == c
}
