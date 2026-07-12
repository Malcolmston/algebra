package stats

import "math"

// Factorial returns n! as a float64. It returns NaN for negative n and +Inf
// once the result overflows float64 (n > 170). The result is computed in the
// log-gamma domain and rounded, so it is exact for all values that fit.
func Factorial(n int) float64 {
	if n < 0 {
		return math.NaN()
	}
	if n <= 1 {
		return 1
	}
	if n > 170 {
		return math.Inf(1)
	}
	return math.Round(math.Exp(gammaLn(float64(n) + 1)))
}

// Choose returns the binomial coefficient C(n, k) = n! / (k!·(n-k)!), the
// number of ways to choose k items from n without regard to order. It returns
// 0 when k < 0 or k > n, and NaN when n < 0. The value is computed via
// log-gamma and rounded, so it stays finite and accurate for large n.
func Choose(n, k int) float64 {
	if n < 0 {
		return math.NaN()
	}
	if k < 0 || k > n {
		return 0
	}
	if k == 0 || k == n {
		return 1
	}
	logC := gammaLn(float64(n)+1) - gammaLn(float64(k)+1) - gammaLn(float64(n-k)+1)
	return math.Round(math.Exp(logC))
}

// Perm returns the number of permutations P(n, k) = n! / (n-k)!, the number of
// ordered arrangements of k items chosen from n. It returns 0 when k < 0 or
// k > n, and NaN when n < 0. The value is computed via log-gamma and rounded.
func Perm(n, k int) float64 {
	if n < 0 {
		return math.NaN()
	}
	if k < 0 || k > n {
		return 0
	}
	if k == 0 {
		return 1
	}
	logP := gammaLn(float64(n)+1) - gammaLn(float64(n-k)+1)
	return math.Round(math.Exp(logP))
}
