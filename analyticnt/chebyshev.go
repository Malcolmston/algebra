package analyticnt

import "math"

// VonMangoldt returns the von Mangoldt function Λ(n): ln p if n = p^k for a
// prime p and integer k >= 1, and 0 otherwise. Λ(1) = 0.
func VonMangoldt(n int64) float64 {
	if n < 2 {
		return 0
	}
	ok, p, _ := IsPrimePower(n)
	if ok {
		return math.Log(float64(p))
	}
	return 0
}

// VonMangoldtLambda is a synonym for VonMangoldt.
func VonMangoldtLambda(n int64) float64 { return VonMangoldt(n) }

// ChebyshevTheta returns the first Chebyshev function θ(x) = Σ_{p ≤ x} ln p,
// summing the natural logarithm over primes p at most x.
func ChebyshevTheta(x float64) float64 {
	if x < 2 {
		return 0
	}
	n := int(math.Floor(x))
	s := Sieve(n)
	sum := 0.0
	for i := 2; i <= n; i++ {
		if s[i] {
			sum += math.Log(float64(i))
		}
	}
	return sum
}

// ChebyshevPsi returns the second Chebyshev function ψ(x) = Σ_{n ≤ x} Λ(n) =
// Σ_{p^k ≤ x} ln p. By the prime number theorem ψ(x) ~ x.
func ChebyshevPsi(x float64) float64 {
	if x < 2 {
		return 0
	}
	n := int(math.Floor(x))
	// ψ(x) = Σ_{p ≤ x} ⌊log_p x⌋ · ln p.
	s := Sieve(n)
	sum := 0.0
	for p := 2; p <= n; p++ {
		if !s[p] {
			continue
		}
		lp := math.Log(float64(p))
		pk := int64(p)
		for pk <= int64(n) {
			sum += lp
			if pk > int64(n)/int64(p) {
				break
			}
			pk *= int64(p)
		}
	}
	return sum
}

// SecondChebyshev is a synonym for ChebyshevPsi.
func SecondChebyshev(x float64) float64 { return ChebyshevPsi(x) }

// ChebyshevThetaExp returns e^{θ(x)} = Π_{p ≤ x} p, the primorial of x, as a
// float64. For large x this overflows to +Inf; use with care.
func ChebyshevThetaExp(x float64) float64 {
	return math.Exp(ChebyshevTheta(x))
}

// PrimeHarmonic returns Σ_{p ≤ x} 1/p, the sum of reciprocals of primes up to
// x. By Mertens' theorem it grows like ln ln x + M.
func PrimeHarmonic(x float64) float64 {
	if x < 2 {
		return 0
	}
	n := int(math.Floor(x))
	s := Sieve(n)
	sum := 0.0
	for i := 2; i <= n; i++ {
		if s[i] {
			sum += 1 / float64(i)
		}
	}
	return sum
}

// PrimeLogHarmonic returns Σ_{p ≤ x} (ln p)/p, which by Mertens' first theorem
// is asymptotic to ln x.
func PrimeLogHarmonic(x float64) float64 {
	if x < 2 {
		return 0
	}
	n := int(math.Floor(x))
	s := Sieve(n)
	sum := 0.0
	for i := 2; i <= n; i++ {
		if s[i] {
			sum += math.Log(float64(i)) / float64(i)
		}
	}
	return sum
}

// MertensProduct returns Π_{p ≤ x} (1 − 1/p). By Mertens' third theorem this is
// asymptotic to e^{−γ}/ln x.
func MertensProduct(x float64) float64 {
	if x < 2 {
		return 1
	}
	n := int(math.Floor(x))
	s := Sieve(n)
	prod := 1.0
	for i := 2; i <= n; i++ {
		if s[i] {
			prod *= 1 - 1/float64(i)
		}
	}
	return prod
}

// VonMangoldtSum returns Σ_{n ≤ N} Λ(n), which equals ψ(N). It is provided as an
// integer-argument companion to ChebyshevPsi.
func VonMangoldtSum(N int64) float64 {
	return ChebyshevPsi(float64(N))
}
