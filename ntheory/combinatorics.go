package ntheory

import "math/big"

// Factorial returns n! as an arbitrary-precision integer. n must be
// non-negative; Factorial(0) == 1.
func Factorial(n int64) *big.Int {
	if n < 0 {
		panic("ntheory: Factorial requires n >= 0")
	}
	result := big.NewInt(1)
	for i := int64(2); i <= n; i++ {
		result.Mul(result, big.NewInt(i))
	}
	return result
}

// DoubleFactorial returns n!! = n·(n-2)·(n-4)·…, the product of every second
// integer down to 1 or 2. By convention DoubleFactorial(0) ==
// DoubleFactorial(-1) == 1. n must be >= -1.
func DoubleFactorial(n int64) *big.Int {
	if n < -1 {
		panic("ntheory: DoubleFactorial requires n >= -1")
	}
	result := big.NewInt(1)
	for i := n; i > 1; i -= 2 {
		result.Mul(result, big.NewInt(i))
	}
	return result
}

// Binomial returns the binomial coefficient C(n, k) = n! / (k!·(n-k)!) as an
// arbitrary-precision integer. It returns 0 when k < 0 or k > n. n must be
// non-negative.
func Binomial(n, k int64) *big.Int {
	if n < 0 {
		panic("ntheory: Binomial requires n >= 0")
	}
	if k < 0 || k > n {
		return big.NewInt(0)
	}
	return new(big.Int).Binomial(n, k)
}

// Multinomial returns the multinomial coefficient (Σkᵢ)! / (k₀!·k₁!·…) as an
// arbitrary-precision integer. Every kᵢ must be non-negative. Multinomial with
// no arguments returns 1.
func Multinomial(ks ...int64) *big.Int {
	var total int64
	for _, k := range ks {
		if k < 0 {
			panic("ntheory: Multinomial requires non-negative arguments")
		}
		total += k
	}
	// (total)! / Π kᵢ! computed as a running product of binomials, which keeps
	// every intermediate value integral: C(total, k0)·C(total-k0, k1)·…
	result := big.NewInt(1)
	remaining := total
	for _, k := range ks {
		result.Mul(result, new(big.Int).Binomial(remaining, k))
		remaining -= k
	}
	return result
}

// Permutations returns the number of ordered k-permutations of n objects,
// nPr = n! / (n-k)!, as an arbitrary-precision integer. It returns 0 when k < 0
// or k > n. n must be non-negative.
func Permutations(n, k int64) *big.Int {
	if n < 0 {
		panic("ntheory: Permutations requires n >= 0")
	}
	if k < 0 || k > n {
		return big.NewInt(0)
	}
	result := big.NewInt(1)
	for i := n - k + 1; i <= n; i++ {
		result.Mul(result, big.NewInt(i))
	}
	return result
}

// CatalanNumber returns the n-th Catalan number Cₙ = C(2n, n) / (n+1) as an
// arbitrary-precision integer. n must be non-negative; CatalanNumber(0) == 1.
func CatalanNumber(n int64) *big.Int {
	if n < 0 {
		panic("ntheory: CatalanNumber requires n >= 0")
	}
	c := new(big.Int).Binomial(2*n, n)
	return c.Div(c, big.NewInt(n+1))
}

// StirlingSecond returns the Stirling number of the second kind S(n, k): the
// number of ways to partition a set of n labeled elements into k non-empty
// unlabeled subsets, as an arbitrary-precision integer. n and k must be
// non-negative. S(0, 0) == 1 and S(n, 0) == 0 for n > 0.
func StirlingSecond(n, k int64) *big.Int {
	if n < 0 || k < 0 {
		panic("ntheory: StirlingSecond requires n >= 0 and k >= 0")
	}
	if k > n {
		return big.NewInt(0)
	}
	// Row-by-row dynamic programming on S(i, j) = j·S(i-1, j) + S(i-1, j-1).
	prev := make([]*big.Int, k+1)
	for j := range prev {
		prev[j] = big.NewInt(0)
	}
	prev[0] = big.NewInt(1) // S(0, 0) = 1
	for i := int64(1); i <= n; i++ {
		curr := make([]*big.Int, k+1)
		curr[0] = big.NewInt(0)
		for j := int64(1); j <= k; j++ {
			term := new(big.Int).Mul(big.NewInt(j), prev[j])
			term.Add(term, prev[j-1])
			curr[j] = term
		}
		prev = curr
	}
	return prev[k]
}

// Partition returns p(n), the number of distinct ways to write the
// non-negative integer n as a sum of positive integers where order does not
// matter, as an arbitrary-precision integer. Partition(0) == 1.
//
// It uses Euler's pentagonal number theorem recurrence, which is exact.
func Partition(n int64) *big.Int {
	if n < 0 {
		panic("ntheory: Partition requires n >= 0")
	}
	p := make([]*big.Int, n+1)
	p[0] = big.NewInt(1)
	for m := int64(1); m <= n; m++ {
		sum := big.NewInt(0)
		// Generalized pentagonal numbers g = k(3k-1)/2 for k = 1, -1, 2, -2, …
		for k := int64(1); ; k++ {
			g1 := k * (3*k - 1) / 2
			g2 := k * (3*k + 1) / 2
			if g1 > m && g2 > m {
				break
			}
			sign := 1
			if k%2 == 0 {
				sign = -1
			}
			if g1 <= m {
				if sign > 0 {
					sum.Add(sum, p[m-g1])
				} else {
					sum.Sub(sum, p[m-g1])
				}
			}
			if g2 <= m {
				if sign > 0 {
					sum.Add(sum, p[m-g2])
				} else {
					sum.Sub(sum, p[m-g2])
				}
			}
		}
		p[m] = sum
	}
	return p[n]
}
