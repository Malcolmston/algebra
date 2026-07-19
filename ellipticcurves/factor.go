package ellipticcurves

import (
	"math/big"
	"sort"
)

// Gcd returns the greatest common divisor of a and b as a non-negative integer.
func Gcd(a, b *big.Int) *big.Int {
	return new(big.Int).GCD(nil, nil, new(big.Int).Abs(a), new(big.Int).Abs(b))
}

// Lcm returns the least common multiple of a and b as a non-negative integer.
// Lcm(0, x) is defined to be 0.
func Lcm(a, b *big.Int) *big.Int {
	if a.Sign() == 0 || b.Sign() == 0 {
		return big.NewInt(0)
	}
	g := Gcd(a, b)
	l := new(big.Int).Div(new(big.Int).Abs(a), g)
	l.Mul(l, new(big.Int).Abs(b))
	return l
}

// LcmSlice returns the least common multiple of a slice of integers, or 1 for
// an empty slice.
func LcmSlice(xs []*big.Int) *big.Int {
	acc := big.NewInt(1)
	for _, x := range xs {
		acc = Lcm(acc, x)
	}
	return acc
}

// IntSqrt returns the integer square root floor(sqrt(n)) for n >= 0. It returns
// 0 for negative input.
func IntSqrt(n *big.Int) *big.Int {
	if n.Sign() < 0 {
		return big.NewInt(0)
	}
	return new(big.Int).Sqrt(n)
}

// IsPerfectSquare reports whether n is a perfect square, returning the square
// root when it is.
func IsPerfectSquare(n *big.Int) (*big.Int, bool) {
	if n.Sign() < 0 {
		return nil, false
	}
	r := new(big.Int).Sqrt(n)
	if new(big.Int).Mul(r, r).Cmp(n) == 0 {
		return r, true
	}
	return nil, false
}

// TrialDivision returns the smallest prime factor of n found by trial division
// up to the given limit, or nil if none is found in range. It is used as the
// first stage of Factorize.
func TrialDivision(n, limit *big.Int) *big.Int {
	if n.Sign() < 0 {
		n = new(big.Int).Abs(n)
	}
	if n.Bit(0) == 0 {
		return big.NewInt(2)
	}
	d := big.NewInt(3)
	for d.Cmp(limit) <= 0 {
		if new(big.Int).Mod(n, d).Sign() == 0 {
			return new(big.Int).Set(d)
		}
		d.Add(d, bigTwo)
	}
	return nil
}

// PollardRho returns a non-trivial factor of the composite odd integer n using
// Pollard's rho algorithm, or nil on failure. It is deterministic given its
// internal polynomial and does not consult any random source.
func PollardRho(n *big.Int) *big.Int {
	if n.Bit(0) == 0 {
		return big.NewInt(2)
	}
	for c := int64(1); c < 100; c++ {
		x := big.NewInt(2)
		y := big.NewInt(2)
		d := big.NewInt(1)
		cc := big.NewInt(c)
		f := func(v *big.Int) *big.Int {
			t := new(big.Int).Mul(v, v)
			t.Add(t, cc)
			return t.Mod(t, n)
		}
		for d.Cmp(bigOne) == 0 {
			x = f(x)
			y = f(f(y))
			diff := new(big.Int).Sub(x, y)
			diff.Abs(diff)
			if diff.Sign() == 0 {
				break
			}
			d = Gcd(diff, n)
		}
		if d.Cmp(bigOne) != 0 && d.Cmp(n) != 0 {
			return d
		}
	}
	return nil
}

// Factorize returns the prime factorization of |n| as a map from prime to
// exponent. Factorize(0) and Factorize(1) return an empty map. It combines
// trial division with Pollard's rho and a primality test.
func Factorize(n *big.Int) map[*big.Int]int {
	result := make(map[*big.Int]int)
	m := new(big.Int).Abs(n)
	if m.Cmp(bigOne) <= 0 {
		return result
	}
	addFactor := func(p *big.Int) {
		for key := range result {
			if key.Cmp(p) == 0 {
				result[key]++
				return
			}
		}
		result[new(big.Int).Set(p)] = 1
	}
	// Strip small primes.
	small := []int64{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}
	for _, s := range small {
		sp := big.NewInt(s)
		for new(big.Int).Mod(m, sp).Sign() == 0 {
			m.Div(m, sp)
			addFactor(sp)
		}
	}
	if m.Cmp(bigOne) == 0 {
		return result
	}
	var recurse func(v *big.Int)
	recurse = func(v *big.Int) {
		if v.Cmp(bigOne) == 0 {
			return
		}
		if v.ProbablyPrime(20) {
			addFactor(v)
			return
		}
		if r, ok := IsPerfectSquare(v); ok {
			recurse(r)
			recurse(r)
			return
		}
		var factor *big.Int
		limit := big.NewInt(1000000)
		factor = TrialDivision(v, limit)
		if factor == nil {
			factor = PollardRho(v)
		}
		if factor == nil {
			// Give up and treat v as prime (best effort).
			addFactor(v)
			return
		}
		q := new(big.Int).Div(v, factor)
		recurse(factor)
		recurse(q)
	}
	recurse(m)
	return result
}

// PrimeFactors returns the distinct prime factors of |n| in ascending order.
func PrimeFactors(n *big.Int) []*big.Int {
	fs := Factorize(n)
	out := make([]*big.Int, 0, len(fs))
	for p := range fs {
		out = append(out, p)
	}
	sort.Slice(out, func(i, j int) bool { return out[i].Cmp(out[j]) < 0 })
	return out
}

// Divisors returns all positive divisors of |n| in ascending order. Divisors(0)
// returns an empty slice.
func Divisors(n *big.Int) []*big.Int {
	m := new(big.Int).Abs(n)
	if m.Sign() == 0 {
		return nil
	}
	fs := Factorize(m)
	divs := []*big.Int{big.NewInt(1)}
	for prime, exp := range fs {
		var next []*big.Int
		pe := big.NewInt(1)
		for e := 0; e <= exp; e++ {
			for _, d := range divs {
				next = append(next, new(big.Int).Mul(d, pe))
			}
			pe = new(big.Int).Mul(pe, prime)
		}
		divs = next
	}
	sort.Slice(divs, func(i, j int) bool { return divs[i].Cmp(divs[j]) < 0 })
	return divs
}

// NumDivisors returns the number of positive divisors of |n|.
func NumDivisors(n *big.Int) int {
	fs := Factorize(n)
	count := 1
	for _, exp := range fs {
		count *= exp + 1
	}
	return count
}

// EulerTotient returns Euler's totient phi(|n|), the count of integers in
// [1, n] coprime to n. EulerTotient(0) is 0 and EulerTotient(1) is 1.
func EulerTotient(n *big.Int) *big.Int {
	m := new(big.Int).Abs(n)
	if m.Sign() == 0 {
		return big.NewInt(0)
	}
	if m.Cmp(bigOne) == 0 {
		return big.NewInt(1)
	}
	result := new(big.Int).Set(m)
	for _, p := range PrimeFactors(m) {
		// result = result * (1 - 1/p) = result / p * (p-1).
		result.Div(result, p)
		result.Mul(result, new(big.Int).Sub(p, bigOne))
	}
	return result
}
