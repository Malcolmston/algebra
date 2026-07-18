package probability

import (
	"math"
	"sort"
)

// Distribution is a finitely-supported discrete probability distribution over
// real-valued outcomes. Outcomes[i] occurs with probability Probs[i]. The two
// slices always have equal length.
//
// Distributions produced by the constructors and by the transform and
// convolution methods maintain a canonical form: Outcomes is sorted in strictly
// ascending order (duplicate outcomes are merged by summing their
// probabilities) and every probability is non-negative and sums to one within
// [probabilityTol]. Methods assume this invariant; mutating the fields directly
// may invalidate cumulative and quantile queries.
type Distribution struct {
	// Outcomes holds the distinct outcome values in ascending order.
	Outcomes []float64
	// Probs holds the probability mass of the corresponding outcome.
	Probs []float64
}

// NewDistribution builds a Distribution from parallel outcome and probability
// slices. Duplicate outcomes are merged and the support is sorted. It returns an
// error if the slices differ in length, are empty, contain a negative or
// non-finite probability, or if the probabilities do not sum to one within
// [probabilityTol].
func NewDistribution(outcomes, probs []float64) (Distribution, error) {
	if len(outcomes) != len(probs) {
		return Distribution{}, probabilityErrorf("NewDistribution: outcomes/probs length mismatch %d != %d", len(outcomes), len(probs))
	}
	if len(outcomes) == 0 {
		return Distribution{}, probabilityErrorf("NewDistribution: empty support")
	}
	for i, p := range probs {
		if math.IsNaN(p) || math.IsInf(p, 0) {
			return Distribution{}, probabilityErrorf("NewDistribution: non-finite probability at index %d", i)
		}
		if p < 0 {
			return Distribution{}, probabilityErrorf("NewDistribution: negative probability %g at index %d", p, i)
		}
		if math.IsNaN(outcomes[i]) || math.IsInf(outcomes[i], 0) {
			return Distribution{}, probabilityErrorf("NewDistribution: non-finite outcome at index %d", i)
		}
	}
	if s := probabilitySum(probs); probabilityAbs(s-1) > probabilityTol {
		return Distribution{}, probabilityErrorf("NewDistribution: probabilities sum to %g, not 1", s)
	}
	outs, ps := probabilityMerge(outcomes, probs)
	return Distribution{Outcomes: outs, Probs: ps}, nil
}

// Uniform returns the uniform distribution over the given distinct outcomes,
// assigning probability 1/n to each. It returns an error if outcomes is empty.
func Uniform(outcomes []float64) (Distribution, error) {
	if len(outcomes) == 0 {
		return Distribution{}, probabilityErrorf("Uniform: empty support")
	}
	p := 1.0 / float64(len(outcomes))
	probs := make([]float64, len(outcomes))
	for i := range probs {
		probs[i] = p
	}
	return NewDistribution(outcomes, probs)
}

// DiscreteUniform returns the uniform distribution over the consecutive integers
// a, a+1, …, b (inclusive), each with probability 1/(b-a+1). It returns an error
// if b < a.
func DiscreteUniform(a, b int) (Distribution, error) {
	if b < a {
		return Distribution{}, probabilityErrorf("DiscreteUniform: empty range [%d,%d]", a, b)
	}
	n := b - a + 1
	outs := make([]float64, n)
	for i := 0; i < n; i++ {
		outs[i] = float64(a + i)
	}
	return Uniform(outs)
}

// Bernoulli returns the Bernoulli distribution with success probability p,
// supported on {0, 1} with P(1) = p and P(0) = 1-p. It returns an error if p is
// outside [0, 1].
func Bernoulli(p float64) (Distribution, error) {
	if p < 0 || p > 1 || math.IsNaN(p) {
		return Distribution{}, probabilityErrorf("Bernoulli: p=%g out of range [0,1]", p)
	}
	return NewDistribution([]float64{0, 1}, []float64{1 - p, p})
}

// Binomial returns the binomial distribution for n independent Bernoulli trials
// each with success probability p, supported on {0, 1, …, n}. It returns an
// error if n is negative or p is outside [0, 1].
func Binomial(n int, p float64) (Distribution, error) {
	if n < 0 {
		return Distribution{}, probabilityErrorf("Binomial: negative n=%d", n)
	}
	if p < 0 || p > 1 || math.IsNaN(p) {
		return Distribution{}, probabilityErrorf("Binomial: p=%g out of range [0,1]", p)
	}
	outs := make([]float64, n+1)
	probs := make([]float64, n+1)
	for k := 0; k <= n; k++ {
		outs[k] = float64(k)
		probs[k] = probabilityBinomPMF(n, k, p)
	}
	return NewDistribution(outs, probs)
}

// Poisson returns the Poisson distribution with rate lambda truncated to the
// support {0, 1, …, kmax} and renormalized so its probabilities sum to one. For
// a kmax comfortably larger than lambda the truncation error is negligible. It
// returns an error if lambda is negative or kmax is negative.
func Poisson(lambda float64, kmax int) (Distribution, error) {
	if lambda < 0 || math.IsNaN(lambda) {
		return Distribution{}, probabilityErrorf("Poisson: negative lambda=%g", lambda)
	}
	if kmax < 0 {
		return Distribution{}, probabilityErrorf("Poisson: negative kmax=%d", kmax)
	}
	outs := make([]float64, kmax+1)
	probs := make([]float64, kmax+1)
	// P(k) = e^{-lambda} lambda^k / k!, built via the recurrence P(k) = P(k-1)*lambda/k.
	term := math.Exp(-lambda)
	sum := 0.0
	for k := 0; k <= kmax; k++ {
		if k > 0 {
			term *= lambda / float64(k)
		}
		outs[k] = float64(k)
		probs[k] = term
		sum += term
	}
	for k := range probs {
		probs[k] /= sum
	}
	return NewDistribution(outs, probs)
}

// Geometric returns the geometric distribution for the number of trials up to
// and including the first success (support {1, 2, …, kmax}) with per-trial
// success probability p, truncated to kmax and renormalized. It returns an error
// if p is outside (0, 1] or kmax is less than one.
func Geometric(p float64, kmax int) (Distribution, error) {
	if p <= 0 || p > 1 || math.IsNaN(p) {
		return Distribution{}, probabilityErrorf("Geometric: p=%g out of range (0,1]", p)
	}
	if kmax < 1 {
		return Distribution{}, probabilityErrorf("Geometric: kmax=%d must be >= 1", kmax)
	}
	outs := make([]float64, kmax)
	probs := make([]float64, kmax)
	sum := 0.0
	term := p // P(1) = p; P(k) = P(k-1)*(1-p).
	for k := 1; k <= kmax; k++ {
		outs[k-1] = float64(k)
		probs[k-1] = term
		sum += term
		term *= (1 - p)
	}
	for i := range probs {
		probs[i] /= sum
	}
	return NewDistribution(outs, probs)
}

// probabilityBinomPMF returns the binomial probability mass C(n,k) p^k (1-p)^{n-k}
// computed via lgamma for numerical stability across the full range of n.
func probabilityBinomPMF(n, k int, p float64) float64 {
	if k < 0 || k > n {
		return 0
	}
	if p == 0 {
		if k == 0 {
			return 1
		}
		return 0
	}
	if p == 1 {
		if k == n {
			return 1
		}
		return 0
	}
	logC, _ := math.Lgamma(float64(n + 1))
	lg1, _ := math.Lgamma(float64(k + 1))
	lg2, _ := math.Lgamma(float64(n - k + 1))
	logC -= lg1 + lg2
	logP := logC + float64(k)*math.Log(p) + float64(n-k)*math.Log(1-p)
	return math.Exp(logP)
}

// Len returns the number of distinct outcomes in the support.
func (d Distribution) Len() int { return len(d.Outcomes) }

// Support returns a copy of the outcome values in ascending order.
func (d Distribution) Support() []float64 {
	out := make([]float64, len(d.Outcomes))
	copy(out, d.Outcomes)
	return out
}

// Min returns the smallest outcome in the support.
func (d Distribution) Min() float64 { return d.Outcomes[0] }

// Max returns the largest outcome in the support.
func (d Distribution) Max() float64 { return d.Outcomes[len(d.Outcomes)-1] }

// Validate reports whether the distribution is well-formed: equal-length
// non-empty slices, non-negative finite probabilities, and a total mass of one
// within [probabilityTol]. It returns nil when the distribution is valid.
func (d Distribution) Validate() error {
	if len(d.Outcomes) != len(d.Probs) {
		return probabilityErrorf("Validate: length mismatch %d != %d", len(d.Outcomes), len(d.Probs))
	}
	if len(d.Outcomes) == 0 {
		return probabilityErrorf("Validate: empty support")
	}
	for i, p := range d.Probs {
		if math.IsNaN(p) || p < 0 {
			return probabilityErrorf("Validate: invalid probability %g at index %d", p, i)
		}
	}
	if s := probabilitySum(d.Probs); probabilityAbs(s-1) > probabilityTol {
		return probabilityErrorf("Validate: probabilities sum to %g, not 1", s)
	}
	return nil
}

// Normalize returns a copy of the distribution whose probabilities are rescaled
// to sum to exactly one. It is useful after manually assembling unnormalized
// weights. It returns an error if the total weight is not strictly positive.
func (d Distribution) Normalize() (Distribution, error) {
	s := probabilitySum(d.Probs)
	if s <= 0 || math.IsNaN(s) {
		return Distribution{}, probabilityErrorf("Normalize: non-positive total weight %g", s)
	}
	outs := make([]float64, len(d.Outcomes))
	probs := make([]float64, len(d.Probs))
	copy(outs, d.Outcomes)
	for i, p := range d.Probs {
		probs[i] = p / s
	}
	return Distribution{Outcomes: outs, Probs: probs}, nil
}

// PMF returns the probability mass P(X = x), i.e. the probability assigned to
// the outcome exactly equal to x, or zero if x is not in the support.
func (d Distribution) PMF(x float64) float64 {
	i := sort.SearchFloat64s(d.Outcomes, x)
	if i < len(d.Outcomes) && d.Outcomes[i] == x {
		return d.Probs[i]
	}
	return 0
}

// CDF returns the cumulative distribution function P(X <= x): the total
// probability of all outcomes less than or equal to x.
func (d Distribution) CDF(x float64) float64 {
	sum := 0.0
	for i, o := range d.Outcomes {
		if o <= x {
			sum += d.Probs[i]
		} else {
			break
		}
	}
	return sum
}

// Quantile returns the smallest outcome x such that P(X <= x) >= q, i.e. the
// generalized inverse of the CDF. q is clamped to [0, 1]; a q of zero returns
// the minimum outcome and a q of one returns the maximum.
func (d Distribution) Quantile(q float64) float64 {
	if q <= 0 {
		return d.Outcomes[0]
	}
	if q >= 1 {
		return d.Outcomes[len(d.Outcomes)-1]
	}
	cum := 0.0
	for i, o := range d.Outcomes {
		cum += d.Probs[i]
		if cum >= q-probabilityTol {
			return o
		}
	}
	return d.Outcomes[len(d.Outcomes)-1]
}

// Median returns the 0.5 quantile of the distribution.
func (d Distribution) Median() float64 { return d.Quantile(0.5) }

// Mode returns the outcome carrying the greatest probability mass. When several
// outcomes tie, the smallest such outcome is returned.
func (d Distribution) Mode() float64 {
	best := 0
	for i := 1; i < len(d.Probs); i++ {
		if d.Probs[i] > d.Probs[best] {
			best = i
		}
	}
	return d.Outcomes[best]
}
