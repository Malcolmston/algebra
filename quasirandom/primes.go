package quasirandom

import (
	"errors"
	"sort"
)

// ErrNonPositive is returned when a strictly positive argument was required but
// a value less than one was supplied.
var ErrNonPositive = errors.New("quasirandom: argument must be positive")

// ErrBadBase is returned when a radix argument is smaller than two.
var ErrBadBase = errors.New("quasirandom: base must be at least 2")

// ErrDimension is returned when a dimension argument is out of range.
var ErrDimension = errors.New("quasirandom: dimension out of range")

// IsPrime reports whether n is a prime number. Values below two are not prime.
func IsPrime(n int) bool {
	if n < 2 {
		return false
	}
	if n%2 == 0 {
		return n == 2
	}
	if n%3 == 0 {
		return n == 3
	}
	for i := 5; i*i <= n; i += 6 {
		if n%i == 0 || n%(i+2) == 0 {
			return false
		}
	}
	return true
}

// NextPrime returns the smallest prime strictly greater than n.
func NextPrime(n int) int {
	if n < 2 {
		return 2
	}
	c := n + 1
	for !IsPrime(c) {
		c++
	}
	return c
}

// PrevPrime returns the largest prime strictly smaller than n, or zero when no
// such prime exists (n <= 2).
func PrevPrime(n int) int {
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

// NextPrimeGE returns the smallest prime greater than or equal to n.
func NextPrimeGE(n int) int {
	if n <= 2 {
		return 2
	}
	if IsPrime(n) {
		return n
	}
	return NextPrime(n)
}

// PrimeSieve returns all primes less than or equal to limit in increasing
// order using the sieve of Eratosthenes. A limit below two yields an empty
// slice.
func PrimeSieve(limit int) []int {
	if limit < 2 {
		return []int{}
	}
	composite := make([]bool, limit+1)
	primes := []int{}
	for i := 2; i <= limit; i++ {
		if !composite[i] {
			primes = append(primes, i)
			for j := i * i; j <= limit; j += i {
				composite[j] = true
			}
		}
	}
	return primes
}

// Primes returns the first n prime numbers (2, 3, 5, ...) in increasing order.
// It returns an error when n is negative.
func Primes(n int) ([]int, error) {
	if n < 0 {
		return nil, ErrNonPositive
	}
	out := make([]int, 0, n)
	c := 2
	for len(out) < n {
		if IsPrime(c) {
			out = append(out, c)
		}
		c++
	}
	return out, nil
}

// Prime returns the i-th prime number using one-based indexing, so Prime(1)==2,
// Prime(2)==3 and so on. It returns an error when i is not positive.
func Prime(i int) (int, error) {
	if i < 1 {
		return 0, ErrNonPositive
	}
	ps, err := Primes(i)
	if err != nil {
		return 0, err
	}
	return ps[i-1], nil
}

// PrimeIndex returns the one-based index of the prime p, or zero when p is not
// prime.
func PrimeIndex(p int) int {
	if !IsPrime(p) {
		return 0
	}
	idx := 0
	for c := 2; c <= p; c++ {
		if IsPrime(c) {
			idx++
		}
	}
	return idx
}

// PrimeBases returns the first dim prime numbers, the canonical choice of
// pairwise-coprime bases for a Halton sequence of the given dimension.
func PrimeBases(dim int) ([]int, error) {
	return Primes(dim)
}

// CountPrimesBelow returns the number of primes strictly less than n
// (the prime-counting function pi evaluated just below n).
func CountPrimesBelow(n int) int {
	return len(PrimeSieve(n - 1))
}

// AreCoprimeBases reports whether every pair of the supplied bases is coprime,
// the condition under which a multi-dimensional Halton construction is
// well distributed. An empty or singleton list is trivially coprime.
func AreCoprimeBases(bases []int) bool {
	for i := 0; i < len(bases); i++ {
		for j := i + 1; j < len(bases); j++ {
			if gcd(bases[i], bases[j]) != 1 {
				return false
			}
		}
	}
	return true
}

// gcd returns the greatest common divisor of a and b.
func gcd(a, b int) int {
	if a < 0 {
		a = -a
	}
	if b < 0 {
		b = -b
	}
	for b != 0 {
		a, b = b, a%b
	}
	return a
}

// GCD returns the greatest common divisor of a and b, with GCD(0,0)==0.
func GCD(a, b int) int { return gcd(a, b) }

// sortedUnique returns the sorted, de-duplicated copy of xs.
func sortedUnique(xs []float64) []float64 {
	out := append([]float64(nil), xs...)
	sort.Float64s(out)
	w := 0
	for i, v := range out {
		if i == 0 || v != out[w-1] {
			out[w] = v
			w++
		}
	}
	return out[:w]
}
