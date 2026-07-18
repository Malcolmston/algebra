package stats

import "math"

// This file adds the two discrete distributions in the "count of trials"
// parameterizations that complement the discrete distributions in
// discretedist.go and the continuous families in contdist.go.

// Geometric is a geometric distribution: the number of Bernoulli trials up to
// and including the first success, each trial succeeding with probability P in
// (0, 1]. Its support is the integers k = 1, 2, 3, … (the "shifted" or trials
// parameterization used by SciPy's geom).
type Geometric struct {
	P float64
}

// PMF returns the probability mass P(X = k), the probability that the first
// success occurs on trial k. It is 0 for k < 1.
func (g Geometric) PMF(k int) float64 {
	if k < 1 {
		return 0
	}
	return math.Pow(1-g.P, float64(k-1)) * g.P
}

// CDF returns the cumulative probability P(X <= k). It is 0 for k < 1.
func (g Geometric) CDF(k int) float64 {
	if k < 1 {
		return 0
	}
	return 1 - math.Pow(1-g.P, float64(k))
}

// Quantile returns the smallest integer k such that CDF(k) >= p, for p in
// [0, 1], as a float64. It returns NaN if p is outside [0, 1].
func (g Geometric) Quantile(p float64) float64 {
	if p < 0 || p > 1 || math.IsNaN(p) {
		return math.NaN()
	}
	if p == 0 {
		return 1
	}
	if p >= 1 {
		return math.Inf(1)
	}
	return math.Ceil(math.Log(1-p) / math.Log(1-g.P))
}

// Mean returns the mean of the distribution, 1/P.
func (g Geometric) Mean() float64 { return 1 / g.P }

// Variance returns the variance of the distribution, (1-P)/P².
func (g Geometric) Variance() float64 { return (1 - g.P) / (g.P * g.P) }

// NegativeBinomial is a negative binomial distribution: the number of failures
// before the R-th success in a sequence of independent Bernoulli trials, each
// succeeding with probability P. Its support is the integers k = 0, 1, 2, …
// (the failures parameterization used by SciPy's nbinom). R may be any positive
// real number.
type NegativeBinomial struct {
	R float64
	P float64
}

// PMF returns the probability mass P(X = k), the probability of exactly k
// failures before the R-th success. It is 0 for k < 0.
func (nb NegativeBinomial) PMF(k int) float64 {
	if k < 0 {
		return 0
	}
	fk := float64(k)
	logC := gammaLn(fk+nb.R) - gammaLn(fk+1) - gammaLn(nb.R)
	logP := logC + nb.R*math.Log(nb.P) + fk*math.Log(1-nb.P)
	return math.Exp(logP)
}

// CDF returns the cumulative probability P(X <= k), evaluated with the
// regularized incomplete beta function. It is 0 for k < 0.
func (nb NegativeBinomial) CDF(k int) float64 {
	if k < 0 {
		return 0
	}
	return regularizedIncompleteBeta(nb.R, float64(k)+1, nb.P)
}

// Mean returns the mean of the distribution, R·(1-P)/P.
func (nb NegativeBinomial) Mean() float64 { return nb.R * (1 - nb.P) / nb.P }

// Variance returns the variance of the distribution, R·(1-P)/P².
func (nb NegativeBinomial) Variance() float64 { return nb.R * (1 - nb.P) / (nb.P * nb.P) }
