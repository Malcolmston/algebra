package analyticnt

import (
	"errors"
	"math"
	"sort"
)

// ErrNonPositive is returned or panicked with when a routine that requires a
// positive argument receives a non-positive one.
var ErrNonPositive = errors.New("analyticnt: argument must be positive")

// Sieve returns a boolean slice s of length n+1 where s[i] reports whether i is
// prime, for 0 <= i <= n. It uses the sieve of Eratosthenes. n must be
// non-negative.
func Sieve(n int) []bool {
	if n < 0 {
		panic("analyticnt: Sieve requires n >= 0")
	}
	s := make([]bool, n+1)
	for i := 2; i <= n; i++ {
		s[i] = true
	}
	for p := 2; p*p <= n; p++ {
		if s[p] {
			for m := p * p; m <= n; m += p {
				s[m] = false
			}
		}
	}
	return s
}

// PrimesUpTo returns all primes p with p <= n in increasing order. n may be
// negative, in which case the empty slice is returned.
func PrimesUpTo(n int) []int64 {
	if n < 2 {
		return []int64{}
	}
	s := Sieve(n)
	primes := make([]int64, 0, int(float64(n)/(math.Log(float64(n))-1))+16)
	for i := 2; i <= n; i++ {
		if s[i] {
			primes = append(primes, int64(i))
		}
	}
	return primes
}

// PrimesInRange returns all primes p with lo <= p <= hi in increasing order.
func PrimesInRange(lo, hi int) []int64 {
	if hi < 2 || hi < lo {
		return []int64{}
	}
	if lo < 2 {
		lo = 2
	}
	s := Sieve(hi)
	out := []int64{}
	for i := lo; i <= hi; i++ {
		if s[i] {
			out = append(out, int64(i))
		}
	}
	return out
}

// IsPrime reports whether n is prime using deterministic trial division backed
// by a 6k±1 wheel. It is correct for all int64 inputs.
func IsPrime(n int64) bool {
	if n < 2 {
		return false
	}
	if n%2 == 0 {
		return n == 2
	}
	if n%3 == 0 {
		return n == 3
	}
	for i := int64(5); i*i <= n; i += 6 {
		if n%i == 0 || n%(i+2) == 0 {
			return false
		}
	}
	return true
}

// IsPrimePower reports whether n = p^k for some prime p and k >= 1, and if so
// returns the base prime p and exponent k. When n is not a prime power it
// returns (false, 0, 0).
func IsPrimePower(n int64) (ok bool, p int64, k int) {
	if n < 2 {
		return false, 0, 0
	}
	for d := int64(2); d*d <= n; d++ {
		if n%d == 0 {
			m := n
			e := 0
			for m%d == 0 {
				m /= d
				e++
			}
			if m == 1 {
				return true, d, e
			}
			return false, 0, 0
		}
	}
	// n is prime.
	return true, n, 1
}

// NthPrime returns the n-th prime (1-indexed): NthPrime(1) == 2, NthPrime(2)
// == 3, and so on. n must be >= 1.
func NthPrime(n int) int64 {
	if n < 1 {
		panic("analyticnt: NthPrime requires n >= 1")
	}
	// Upper bound for the n-th prime (Rosser). Valid for n >= 6.
	var limit int
	if n < 6 {
		limit = 15
	} else {
		fn := float64(n)
		limit = int(fn*(math.Log(fn)+math.Log(math.Log(fn)))) + 3
	}
	for {
		s := Sieve(limit)
		count := 0
		for i := 2; i <= limit; i++ {
			if s[i] {
				count++
				if count == n {
					return int64(i)
				}
			}
		}
		limit *= 2
	}
}

// NextPrime returns the smallest prime strictly greater than n.
func NextPrime(n int64) int64 {
	if n < 2 {
		return 2
	}
	c := n + 1
	for !IsPrime(c) {
		c++
	}
	return c
}

// PrevPrime returns the largest prime strictly less than n, or 0 if there is
// none (n <= 2).
func PrevPrime(n int64) int64 {
	if n <= 2 {
		return 0
	}
	c := n - 1
	for c >= 2 && !IsPrime(c) {
		c--
	}
	if c < 2 {
		return 0
	}
	return c
}

// PrimePi returns the exact prime-counting function pi(x): the number of primes
// p with p <= x. It sieves up to floor(x) and is intended for moderate x.
func PrimePi(x int64) int64 {
	if x < 2 {
		return 0
	}
	s := Sieve(int(x))
	var count int64
	for i := 2; i <= int(x); i++ {
		if s[i] {
			count++
		}
	}
	return count
}

// PrimePiSieve returns pi(k) for every k from 0 to n as a cumulative slice of
// length n+1, where result[k] is the number of primes <= k.
func PrimePiSieve(n int) []int64 {
	if n < 0 {
		panic("analyticnt: PrimePiSieve requires n >= 0")
	}
	s := Sieve(n)
	out := make([]int64, n+1)
	var c int64
	for i := 0; i <= n; i++ {
		if i >= 2 && s[i] {
			c++
		}
		out[i] = c
	}
	return out
}

// PrimeCountBetween returns the number of primes p with lo < p <= hi.
func PrimeCountBetween(lo, hi int64) int64 {
	return PrimePi(hi) - PrimePi(lo)
}

// PrimorialPrimes returns the first k primes whose product is the k-th
// primorial. It is a convenience wrapper returning the prime list only.
func PrimorialPrimes(k int) []int64 {
	if k < 0 {
		panic("analyticnt: PrimorialPrimes requires k >= 0")
	}
	out := make([]int64, 0, k)
	p := int64(1)
	for len(out) < k {
		p = NextPrime(p)
		out = append(out, p)
	}
	return out
}

// PrimeIndex returns the 1-based index of prime p in the sequence of primes
// (so PrimeIndex(2) == 1). If p is not prime it returns 0.
func PrimeIndex(p int64) int64 {
	if !IsPrime(p) {
		return 0
	}
	return PrimePi(p)
}

// primeCache lazily caches a prefix of the primes for internal counting
// routines. It returns primes up to at least limit.
func primesAtLeast(limit int) []int64 {
	return PrimesUpTo(limit)
}

// sortedSearchInt64 returns the number of elements in the sorted slice a that
// are <= v.
func sortedSearchInt64(a []int64, v int64) int {
	return sort.Search(len(a), func(i int) bool { return a[i] > v })
}
