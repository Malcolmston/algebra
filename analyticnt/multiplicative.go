package analyticnt

import "math"

// MobiusMu returns the Möbius function μ(n): 1 if n is a square-free positive
// integer with an even number of prime factors, −1 if square-free with an odd
// number, and 0 if n has a squared prime factor. μ(1) = 1.
func MobiusMu(n int64) int {
	if n <= 0 {
		panic("analyticnt: MobiusMu requires n >= 1")
	}
	if n == 1 {
		return 1
	}
	factors := 0
	for p := int64(2); p*p <= n; p++ {
		if n%p == 0 {
			n /= p
			if n%p == 0 {
				return 0
			}
			factors++
		}
	}
	if n > 1 {
		factors++
	}
	if factors%2 == 0 {
		return 1
	}
	return -1
}

// MobiusSieve returns μ(k) for every k from 0 to n as a slice of length n+1,
// with entry 0 left as 0. It is far faster than repeated MobiusMu calls.
func MobiusSieve(n int) []int {
	if n < 0 {
		panic("analyticnt: MobiusSieve requires n >= 0")
	}
	mu := make([]int, n+1)
	primes := make([]int, 0)
	if n >= 1 {
		mu[1] = 1
	}
	isComposite := make([]bool, n+1)
	for i := 2; i <= n; i++ {
		if !isComposite[i] {
			primes = append(primes, i)
			mu[i] = -1
		}
		for _, p := range primes {
			if i*p > n {
				break
			}
			isComposite[i*p] = true
			if i%p == 0 {
				mu[i*p] = 0
				break
			}
			mu[i*p] = -mu[i]
		}
	}
	return mu
}

// MertensFunction returns the Mertens function M(n) = Σ_{k=1}^{n} μ(k). n must
// be >= 0; M(0) = 0.
func MertensFunction(n int) int64 {
	if n < 0 {
		panic("analyticnt: MertensFunction requires n >= 0")
	}
	mu := MobiusSieve(n)
	var s int64
	for k := 1; k <= n; k++ {
		s += int64(mu[k])
	}
	return s
}

// MertensSequence returns the running Mertens values M(0), M(1), …, M(n) as a
// slice of length n+1.
func MertensSequence(n int) []int64 {
	if n < 0 {
		panic("analyticnt: MertensSequence requires n >= 0")
	}
	mu := MobiusSieve(n)
	out := make([]int64, n+1)
	var s int64
	for k := 1; k <= n; k++ {
		s += int64(mu[k])
		out[k] = s
	}
	return out
}

// EulerPhi returns Euler's totient φ(n), the number of integers in [1, n]
// coprime to n. n must be >= 1.
func EulerPhi(n int64) int64 {
	if n <= 0 {
		panic("analyticnt: EulerPhi requires n >= 1")
	}
	result := n
	m := n
	for p := int64(2); p*p <= m; p++ {
		if m%p == 0 {
			for m%p == 0 {
				m /= p
			}
			result -= result / p
		}
	}
	if m > 1 {
		result -= result / m
	}
	return result
}

// Liouville returns the Liouville function λ(n) = (−1)^{Ω(n)}, where Ω counts
// prime factors with multiplicity. λ(1) = 1.
func Liouville(n int64) int {
	if BigOmega(n)%2 == 0 {
		return 1
	}
	return -1
}

// Omega returns ω(n), the number of distinct prime factors of n. ω(1) = 0.
func Omega(n int64) int {
	if n <= 1 {
		return 0
	}
	count := 0
	for p := int64(2); p*p <= n; p++ {
		if n%p == 0 {
			count++
			for n%p == 0 {
				n /= p
			}
		}
	}
	if n > 1 {
		count++
	}
	return count
}

// BigOmega returns Ω(n), the number of prime factors of n counted with
// multiplicity. Ω(1) = 0.
func BigOmega(n int64) int {
	if n <= 1 {
		return 0
	}
	count := 0
	for p := int64(2); p*p <= n; p++ {
		for n%p == 0 {
			n /= p
			count++
		}
	}
	if n > 1 {
		count++
	}
	return count
}

// Radical returns rad(n), the product of the distinct prime factors of n.
// rad(1) = 1.
func Radical(n int64) int64 {
	if n <= 1 {
		return 1
	}
	r := int64(1)
	for p := int64(2); p*p <= n; p++ {
		if n%p == 0 {
			r *= p
			for n%p == 0 {
				n /= p
			}
		}
	}
	if n > 1 {
		r *= n
	}
	return r
}

// IsSquareFree reports whether n has no repeated prime factor.
func IsSquareFree(n int64) bool {
	if n <= 0 {
		return false
	}
	return MobiusMu(n) != 0
}

// DivisorCount returns τ(n) = d(n), the number of positive divisors of n.
func DivisorCount(n int64) int64 {
	if n <= 0 {
		panic("analyticnt: DivisorCount requires n >= 1")
	}
	count := int64(1)
	for p := int64(2); p*p <= n; p++ {
		if n%p == 0 {
			e := int64(0)
			for n%p == 0 {
				n /= p
				e++
			}
			count *= e + 1
		}
	}
	if n > 1 {
		count *= 2
	}
	return count
}

// DivisorSigma returns σ_k(n), the sum of the k-th powers of the divisors of n.
// σ_0 = τ (divisor count) and σ_1 is the sum of divisors.
func DivisorSigma(k int, n int64) int64 {
	if n <= 0 {
		panic("analyticnt: DivisorSigma requires n >= 1")
	}
	if k == 0 {
		return DivisorCount(n)
	}
	result := int64(1)
	for p := int64(2); p*p <= n; p++ {
		if n%p == 0 {
			e := 0
			for n%p == 0 {
				n /= p
				e++
			}
			// (p^{k(e+1)} - 1)/(p^k - 1)
			pk := ipow(p, k)
			num := ipow(pk, e+1) - 1
			result *= num / (pk - 1)
		}
	}
	if n > 1 {
		pk := ipow(n, k)
		result *= pk + 1
	}
	return result
}

// ipow returns base^exp for non-negative exp using integer arithmetic.
func ipow(base int64, exp int) int64 {
	r := int64(1)
	for i := 0; i < exp; i++ {
		r *= base
	}
	return r
}

// MangoldtSummatory returns Σ_{n≤N} Λ(n)/n as a partial sum related to the
// logarithmic derivative of ζ; it is a small utility for PNT experiments.
func MangoldtSummatory(N int64) float64 {
	sum := 0.0
	for n := int64(2); n <= N; n++ {
		l := VonMangoldt(n)
		if l != 0 {
			sum += l / float64(n)
		}
	}
	return sum
}

// MobiusSummatoryReciprocal returns Σ_{n≤N} μ(n)/n, which tends to 0 and whose
// value is bounded by 1 in absolute value.
func MobiusSummatoryReciprocal(N int) float64 {
	mu := MobiusSieve(N)
	sum := 0.0
	for n := 1; n <= N; n++ {
		if mu[n] != 0 {
			sum += float64(mu[n]) / float64(n)
		}
	}
	return sum
}

// JordanTotient returns the Jordan totient J_k(n) = n^k Π_{p|n}(1 − p^{-k}),
// generalizing Euler's totient (k = 1).
func JordanTotient(k int, n int64) int64 {
	if n <= 0 || k < 1 {
		panic("analyticnt: JordanTotient requires n >= 1 and k >= 1")
	}
	result := ipow(n, k)
	m := n
	for p := int64(2); p*p <= m; p++ {
		if m%p == 0 {
			for m%p == 0 {
				m /= p
			}
			pk := ipow(p, k)
			result -= result / pk
		}
	}
	if m > 1 {
		pk := ipow(m, k)
		result -= result / pk
	}
	return result
}

// TotientSummatory returns Φ(n) = Σ_{k=1}^{n} φ(k), whose asymptotic is
// 3n²/π².
func TotientSummatory(n int) int64 {
	if n < 0 {
		panic("analyticnt: TotientSummatory requires n >= 0")
	}
	phi := make([]int64, n+1)
	for i := 0; i <= n; i++ {
		phi[i] = int64(i)
	}
	for i := 2; i <= n; i++ {
		if phi[i] == int64(i) { // i is prime
			for j := i; j <= n; j += i {
				phi[j] -= phi[j] / int64(i)
			}
		}
	}
	var s int64
	for k := 1; k <= n; k++ {
		s += phi[k]
	}
	return s
}

// squareFreeDensityEstimate returns the expected count of square-free numbers up
// to n, ~ 6n/π², used in tests as a sanity comparison.
func squareFreeDensityEstimate(n float64) float64 {
	return 6 * n / (math.Pi * math.Pi)
}
