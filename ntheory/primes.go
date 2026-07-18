package ntheory

import (
	"math/big"
	"sort"
)

// mulMod returns (a * b) mod m computed without intermediate overflow by using
// math/big. m must be positive.
func mulMod(a, b, m int64) int64 {
	if m == 1 {
		return 0
	}
	x := big.NewInt(a)
	x.Mul(x, big.NewInt(b))
	x.Mod(x, big.NewInt(m))
	return x.Int64()
}

// witnesses used by the deterministic Miller-Rabin test. Testing against this
// fixed set of bases is proven to be a correct primality test for every
// n < 3,317,044,064,679,887,385,961,981, which covers the entire uint64 range
// and therefore every non-negative int64.
var mrWitnesses = []int64{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}

// IsPrime reports whether n is prime using a deterministic Miller-Rabin test.
//
// The chosen set of witness bases makes the test provably correct for every
// value representable by an int64, so the answer is exact (never merely
// probable). Negative numbers, 0 and 1 are not prime. Carmichael numbers such
// as 561 are correctly reported as composite.
func IsPrime(n int64) bool {
	if n < 2 {
		return false
	}
	for _, p := range []int64{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37} {
		if n == p {
			return true
		}
		if n%p == 0 {
			return false
		}
	}
	// Write n-1 = d * 2^s with d odd.
	d := n - 1
	s := 0
	for d%2 == 0 {
		d /= 2
		s++
	}
	for _, a := range mrWitnesses {
		if a%n == 0 {
			continue
		}
		x := ModPow(a, d, n)
		if x == 1 || x == n-1 {
			continue
		}
		composite := true
		for i := 0; i < s-1; i++ {
			x = mulMod(x, x, n)
			if x == n-1 {
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

// IsPrimeBig reports whether n is (very probably) prime using
// math/big.Int.ProbablyPrime with 20 Miller-Rabin rounds plus a Baillie-PSW
// test. No composite number is known to pass, making it suitable for the
// arbitrary-precision case where a deterministic int64 test does not apply.
// Values below 2 are not prime.
func IsPrimeBig(n *big.Int) bool {
	if n.Sign() < 0 {
		return false
	}
	return n.ProbablyPrime(20)
}

// NextPrime returns the smallest prime strictly greater than n.
func NextPrime(n int64) int64 {
	if n < 2 {
		return 2
	}
	candidate := n + 1
	if candidate%2 == 0 {
		if candidate == 2 {
			return 2
		}
		candidate++
	}
	for !IsPrime(candidate) {
		candidate += 2
	}
	return candidate
}

// PrevPrime returns the greatest prime strictly less than n together with
// ok == true. When no prime is smaller than n (that is, n <= 2) it returns
// (0, false). It is the int64 counterpart of [NextPrime] and delegates to
// [PrevPrimeU64].
func PrevPrime(n int64) (prime int64, ok bool) {
	if n <= 2 {
		return 0, false
	}
	p, found := PrevPrimeU64(uint64(n))
	return int64(p), found
}

// PrimesUpTo returns all primes p with p <= n in ascending order, computed with
// the sieve of Eratosthenes. It returns nil for n < 2.
func PrimesUpTo(n int64) []int64 {
	if n < 2 {
		return nil
	}
	sieve := make([]bool, n+1) // sieve[i] == true means i is composite.
	var primes []int64
	for i := int64(2); i <= n; i++ {
		if sieve[i] {
			continue
		}
		primes = append(primes, i)
		for j := i * i; j <= n; j += i {
			sieve[j] = true
		}
	}
	return primes
}

// PrimePi returns π(n), the number of primes p with p <= n.
func PrimePi(n int64) int64 {
	return int64(len(PrimesUpTo(n)))
}

// Factorize returns the prime factorization of |n| as a map from each prime
// factor to its exponent. Factorize(0) and Factorize(±1) return an empty
// (non-nil) map. The sign of n is ignored.
//
// See [FactorList] for a deterministically ordered slice representation.
func Factorize(n int64) map[int64]int {
	n = abs64(n)
	factors := make(map[int64]int)
	if n < 2 {
		return factors
	}
	for n%2 == 0 {
		factors[2]++
		n /= 2
	}
	for d := int64(3); d*d <= n; d += 2 {
		for n%d == 0 {
			factors[d]++
			n /= d
		}
	}
	if n > 1 {
		factors[n]++
	}
	return factors
}

// PrimePower pairs a prime with its exponent in a factorization.
type PrimePower struct {
	Prime    int64 // Prime is the prime base.
	Exponent int   // Exponent is the power to which Prime is raised.
}

// FactorList returns the prime factorization of |n| as a slice of
// [PrimePower] values sorted by ascending prime. It is the ordered counterpart
// of [Factorize]. FactorList(0) and FactorList(±1) return nil.
func FactorList(n int64) []PrimePower {
	factors := Factorize(n)
	if len(factors) == 0 {
		return nil
	}
	list := make([]PrimePower, 0, len(factors))
	for p, e := range factors {
		list = append(list, PrimePower{Prime: p, Exponent: e})
	}
	sort.Slice(list, func(i, j int) bool { return list[i].Prime < list[j].Prime })
	return list
}

// EulerPhi returns Euler's totient φ(n): the count of integers in [1, n] that
// are coprime to n. The sign of n is ignored. By convention EulerPhi(0) == 0
// and EulerPhi(1) == 1.
func EulerPhi(n int64) int64 {
	n = abs64(n)
	if n == 0 {
		return 0
	}
	result := n
	for p := range Factorize(n) {
		result -= result / p
	}
	return result
}

// MobiusMu returns the Möbius function μ(n): 0 if n is divisible by a square
// greater than 1, otherwise (-1)^k where k is the number of distinct prime
// factors. The sign of n is ignored and μ(1) == 1. MobiusMu(0) == 0.
func MobiusMu(n int64) int {
	n = abs64(n)
	if n == 0 {
		return 0
	}
	if n == 1 {
		return 1
	}
	result := 1
	for _, e := range Factorize(n) {
		if e > 1 {
			return 0
		}
		result = -result
	}
	return result
}

// Radical returns the radical of n: the product of its distinct prime factors.
// The sign of n is ignored. Radical(0) == 0 and Radical(1) == 1.
func Radical(n int64) int64 {
	n = abs64(n)
	if n == 0 {
		return 0
	}
	result := int64(1)
	for p := range Factorize(n) {
		result *= p
	}
	return result
}
