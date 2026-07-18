package stats

import "math"

// This file adds discrete distributions that mirror the Binomial/Poisson API
// (PMF, CDF, Quantile, Mean, Variance) but use the "count of failures"
// support convention.
//
// The exported names deliberately differ from the Geometric and
// NegativeBinomial types declared elsewhere in this package: those use the
// "number of trials" (k >= 1) and real-valued R parameterizations, whereas the
// types here count failures before a success (support k >= 0) and, for the
// negative binomial, take an integer number of successes R. Keeping distinct
// names lets both parameterizations coexist in the same package.

// GeometricFailures is a geometric distribution parameterized by the number of
// failures before the first success. Each independent Bernoulli trial succeeds
// with probability P in (0, 1]. Its support is the integers k = 0, 1, 2, …
// (the "failures" parameterization used by SciPy's geom with loc = -1, i.e. the
// number of failures preceding the first success). This convention is what
// makes its mean (1-P)/P rather than 1/P.
type GeometricFailures struct {
	P float64 // P is the per-trial success probability in (0, 1].
}

// valid reports whether the parameter P lies in the admissible range (0, 1].
func (g GeometricFailures) valid() bool {
	return g.P > 0 && g.P <= 1
}

// PMF returns the probability mass P(X = k), the probability of exactly k
// failures before the first success: (1-P)^k · P for k >= 0, and 0 for k < 0.
// It returns NaN if P is outside (0, 1].
func (g GeometricFailures) PMF(k int) float64 {
	if !g.valid() {
		return math.NaN()
	}
	if k < 0 {
		return 0
	}
	return math.Pow(1-g.P, float64(k)) * g.P
}

// CDF returns the cumulative probability P(X <= k) = 1 - (1-P)^(k+1). It is 0
// for k < 0 and returns NaN if P is outside (0, 1].
func (g GeometricFailures) CDF(k int) float64 {
	if !g.valid() {
		return math.NaN()
	}
	if k < 0 {
		return 0
	}
	return 1 - math.Pow(1-g.P, float64(k+1))
}

// Quantile returns the smallest integer k such that CDF(k) >= p, for p in
// [0, 1]. It uses the closed form ceil(ln(1-p)/ln(1-P)) - 1 clamped at 0, which
// is exact and avoids an iterative search. It returns the sentinel 0 for
// invalid P or for p outside [0, 1], and math.MaxInt for p = 1 when P < 1,
// where no finite k attains a cumulative probability of exactly 1.
func (g GeometricFailures) Quantile(p float64) int {
	if !g.valid() || math.IsNaN(p) || p < 0 || p > 1 {
		return 0
	}
	if p == 0 {
		return 0
	}
	if g.P == 1 {
		// All probability mass sits at k = 0, so CDF(0) = 1 >= p.
		return 0
	}
	if p == 1 {
		return math.MaxInt
	}
	q := math.Ceil(math.Log(1-p)/math.Log(1-g.P)) - 1
	if q < 0 {
		return 0
	}
	return int(q)
}

// Mean returns the mean of the distribution, (1-P)/P. It returns NaN if P is
// outside (0, 1].
func (g GeometricFailures) Mean() float64 {
	if !g.valid() {
		return math.NaN()
	}
	return (1 - g.P) / g.P
}

// Variance returns the variance of the distribution, (1-P)/P². It returns NaN
// if P is outside (0, 1].
func (g GeometricFailures) Variance() float64 {
	if !g.valid() {
		return math.NaN()
	}
	return (1 - g.P) / (g.P * g.P)
}

// NegativeBinomialInt is a negative binomial distribution parameterized by an
// integer number of successes R (>= 1) and success probability P in (0, 1]. It
// models the number of failures k = 0, 1, 2, … before the R-th success in a
// sequence of independent Bernoulli trials (the "failures" parameterization
// used by SciPy's nbinom). With R = 1 it reduces to GeometricFailures.
type NegativeBinomialInt struct {
	R int     // R is the target number of successes (>= 1).
	P float64 // P is the per-trial success probability in (0, 1].
}

// valid reports whether the parameters satisfy R >= 1 and P in (0, 1].
func (nb NegativeBinomialInt) valid() bool {
	return nb.R >= 1 && nb.P > 0 && nb.P <= 1
}

// PMF returns the probability mass P(X = k), the probability of exactly k
// failures before the R-th success. It is evaluated in the log domain via
// gammaLn to stay finite for large arguments:
//
//	exp(gammaLn(k+R) - gammaLn(k+1) - gammaLn(R) + R·ln P + k·ln(1-P)).
//
// It is 0 for k < 0 and returns NaN for invalid parameters.
func (nb NegativeBinomialInt) PMF(k int) float64 {
	if !nb.valid() {
		return math.NaN()
	}
	if k < 0 {
		return 0
	}
	if nb.P >= 1 {
		// P == 1: every trial succeeds, so there are certainly 0 failures.
		if k == 0 {
			return 1
		}
		return 0
	}
	fk := float64(k)
	r := float64(nb.R)
	logP := gammaLn(fk+r) - gammaLn(fk+1) - gammaLn(r) + r*math.Log(nb.P) + fk*math.Log(1-nb.P)
	return math.Exp(logP)
}

// CDF returns the cumulative probability P(X <= k). As a performance measure it
// is evaluated in closed form with the regularized incomplete beta function,
//
//	I_P(R, k+1),
//
// which yields the exact cumulative probability in O(1) rather than summing k+1
// PMF terms. It is 0 for k < 0 and returns NaN for invalid parameters.
func (nb NegativeBinomialInt) CDF(k int) float64 {
	if !nb.valid() {
		return math.NaN()
	}
	if k < 0 {
		return 0
	}
	return regularizedIncompleteBeta(float64(nb.R), float64(k+1), nb.P)
}

// Quantile returns the smallest integer k such that CDF(k) >= p, for p in
// [0, 1], found by scanning k upward from the lower support bound 0 while
// accumulating a running cumulative sum of the PMF. It returns the sentinel 0
// for invalid parameters or for p outside [0, 1], and math.MaxInt for p = 1
// when P < 1, where no finite k attains a cumulative probability of exactly 1.
func (nb NegativeBinomialInt) Quantile(p float64) int {
	if !nb.valid() || math.IsNaN(p) || p < 0 || p > 1 {
		return 0
	}
	if p == 1 {
		// With P == 1 all mass is at 0, so the quantile is 0; otherwise the
		// support is unbounded and CDF(k) < 1 for every finite k, so p = 1 is
		// unattainable. Handle this explicitly because the running cumulative
		// sum below rounds up to exactly 1.0 at a large but finite k.
		if nb.P >= 1 {
			return 0
		}
		return math.MaxInt
	}
	cum := 0.0
	for k := 0; ; k++ {
		pmf := nb.PMF(k)
		cum += pmf
		if cum >= p {
			return k
		}
		if pmf == 0 && k > 0 {
			// The tail has underflowed to zero without reaching p; the
			// requested quantile is unattainable at any finite k.
			return math.MaxInt
		}
	}
}

// Mean returns the mean of the distribution, R·(1-P)/P. It returns NaN for
// invalid parameters.
func (nb NegativeBinomialInt) Mean() float64 {
	if !nb.valid() {
		return math.NaN()
	}
	return float64(nb.R) * (1 - nb.P) / nb.P
}

// Variance returns the variance of the distribution, R·(1-P)/P². It returns NaN
// for invalid parameters.
func (nb NegativeBinomialInt) Variance() float64 {
	if !nb.valid() {
		return math.NaN()
	}
	return float64(nb.R) * (1 - nb.P) / (nb.P * nb.P)
}

// Hypergeometric is a hypergeometric distribution: the number of successes k in
// Draws draws made without replacement from a finite population of N items
// containing K successes. The parameters must satisfy 0 <= K <= N and
// 0 <= Draws <= N. Its support is the integers max(0, Draws-(N-K)) <= k <=
// min(K, Draws).
type Hypergeometric struct {
	N     int // N is the population size.
	K     int // K is the number of successes in the population.
	Draws int // Draws is the number of items drawn without replacement.
}

// valid reports whether the parameters satisfy N >= 1, 0 <= K <= N and
// 0 <= Draws <= N.
func (h Hypergeometric) valid() bool {
	return h.N >= 1 && h.K >= 0 && h.K <= h.N && h.Draws >= 0 && h.Draws <= h.N
}

// support returns the inclusive lower and upper bounds of the support,
// max(0, Draws-(N-K)) and min(K, Draws).
func (h Hypergeometric) support() (int, int) {
	lo := 0
	if v := h.Draws - (h.N - h.K); v > lo {
		lo = v
	}
	hi := h.K
	if h.Draws < hi {
		hi = h.Draws
	}
	return lo, hi
}

// PMF returns the probability mass P(X = k) = C(K, k)·C(N-K, Draws-k)/C(N, Draws)
// using the exact integer binomial coefficient Choose. It is 0 for k outside
// the support and returns NaN for invalid parameters.
func (h Hypergeometric) PMF(k int) float64 {
	if !h.valid() {
		return math.NaN()
	}
	lo, hi := h.support()
	if k < lo || k > hi {
		return 0
	}
	return Choose(h.K, k) * Choose(h.N-h.K, h.Draws-k) / Choose(h.N, h.Draws)
}

// CDF returns the cumulative probability P(X <= k) by accumulating PMF terms
// from the lower support bound up to k. It is 0 below the support, 1 at or
// above the upper support bound, and returns NaN for invalid parameters.
func (h Hypergeometric) CDF(k int) float64 {
	if !h.valid() {
		return math.NaN()
	}
	lo, hi := h.support()
	if k < lo {
		return 0
	}
	if k >= hi {
		return 1
	}
	sum := 0.0
	for i := lo; i <= k; i++ {
		sum += h.PMF(i)
	}
	return sum
}

// Quantile returns the smallest integer k such that CDF(k) >= p, for p in
// [0, 1], found by scanning k upward from the lower support bound while
// carrying a running cumulative sum of the PMF (so the CDF is not recomputed
// from scratch at each step). It returns the sentinel 0 for invalid parameters
// or for p outside [0, 1]. Because the support is finite, p = 1 returns the
// upper support bound.
func (h Hypergeometric) Quantile(p float64) int {
	if !h.valid() || math.IsNaN(p) || p < 0 || p > 1 {
		return 0
	}
	lo, hi := h.support()
	cum := 0.0
	for k := lo; k < hi; k++ {
		cum += h.PMF(k)
		if cum >= p {
			return k
		}
	}
	return hi
}

// Mean returns the mean of the distribution, Draws·K/N. It returns NaN for
// invalid parameters.
func (h Hypergeometric) Mean() float64 {
	if !h.valid() {
		return math.NaN()
	}
	return float64(h.Draws) * float64(h.K) / float64(h.N)
}

// Variance returns the variance of the distribution,
// Draws·(K/N)·(1-K/N)·(N-Draws)/(N-1). A single-item population (N = 1) has no
// variability and yields 0. It returns NaN for invalid parameters.
func (h Hypergeometric) Variance() float64 {
	if !h.valid() {
		return math.NaN()
	}
	if h.N == 1 {
		return 0
	}
	N := float64(h.N)
	K := float64(h.K)
	n := float64(h.Draws)
	return n * (K / N) * (1 - K/N) * (N - n) / (N - 1)
}
