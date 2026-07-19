package bayesian

import (
	"math"
)

// ------------------------------------------------------------------
// Standard non-informative Beta priors
// ------------------------------------------------------------------

// UniformBetaPrior returns Beta(1,1), the uniform prior on a probability.
func UniformBetaPrior() Beta { return Beta{Alpha: 1, Beta: 1} }

// JeffreysBetaPrior returns Beta(½,½), the Jeffreys prior for a Bernoulli
// success probability.
func JeffreysBetaPrior() Beta { return Beta{Alpha: 0.5, Beta: 0.5} }

// HaldaneBetaPrior returns Beta(ε,ε) with a tiny ε, an approximation of the
// improper Haldane prior Beta(0,0) whose posterior mean equals the sample
// proportion.
func HaldaneBetaPrior() Beta { return Beta{Alpha: 1e-9, Beta: 1e-9} }

// ------------------------------------------------------------------
// Beta-Bernoulli / Beta-Binomial conjugate update
// ------------------------------------------------------------------

// BetaBernoulliPosterior returns the posterior Beta(α+successes, β+failures)
// after observing the given counts of Bernoulli successes and failures.
func BetaBernoulliPosterior(prior Beta, successes, failures int) Beta {
	return Beta{Alpha: prior.Alpha + float64(successes), Beta: prior.Beta + float64(failures)}
}

// BetaBinomialPosterior returns the posterior after observing k successes in n
// Binomial trials. It is identical to a Beta-Bernoulli update with k successes
// and n−k failures.
func BetaBinomialPosterior(prior Beta, k, n int) Beta {
	return Beta{Alpha: prior.Alpha + float64(k), Beta: prior.Beta + float64(n-k)}
}

// BetaBernoulliSequentialPosterior applies a sequence of Bernoulli outcomes
// (true = success) to the prior and returns the resulting posterior.
func BetaBernoulliSequentialPosterior(prior Beta, outcomes []bool) Beta {
	post := prior
	for _, o := range outcomes {
		if o {
			post.Alpha++
		} else {
			post.Beta++
		}
	}
	return post
}

// BetaBernoulliPredictiveProb returns the probability that the next Bernoulli
// trial is a success under the posterior, which equals the posterior mean.
func BetaBernoulliPredictiveProb(post Beta) float64 {
	return post.Mean()
}

// ------------------------------------------------------------------
// Gamma-Poisson conjugate update
// ------------------------------------------------------------------

// GammaPoissonPosterior returns the posterior Gamma(shape+Σcounts, rate+n)
// after observing n Poisson counts whose total is sumCounts.
func GammaPoissonPosterior(prior Gamma, sumCounts, n int) Gamma {
	return Gamma{Shape: prior.Shape + float64(sumCounts), Rate: prior.Rate + float64(n)}
}

// GammaPoissonPosteriorData returns the Gamma-Poisson posterior for a slice of
// observed non-negative integer counts.
func GammaPoissonPosteriorData(prior Gamma, counts []int) Gamma {
	sum := 0
	for _, c := range counts {
		sum += c
	}
	return GammaPoissonPosterior(prior, sum, len(counts))
}

// GammaExponentialPosterior returns the posterior for the rate of an
// exponential likelihood with a Gamma prior, after observing n samples whose
// sum is sumX: Gamma(shape+n, rate+sumX).
func GammaExponentialPosterior(prior Gamma, sumX float64, n int) Gamma {
	return Gamma{Shape: prior.Shape + float64(n), Rate: prior.Rate + sumX}
}

// ------------------------------------------------------------------
// Normal-Normal update with known sampling variance
// ------------------------------------------------------------------

// NormalKnownVariancePosterior returns the posterior over the mean of a normal
// likelihood with known sampling variance sigma², given a normal prior and the
// summary statistics n (sample size) and xbar (sample mean).
func NormalKnownVariancePosterior(prior Normal, sigma2 float64, n int, xbar float64) Normal {
	priorPrec := 1 / (prior.Sigma * prior.Sigma)
	dataPrec := float64(n) / sigma2
	postPrec := priorPrec + dataPrec
	postVar := 1 / postPrec
	postMean := postVar * (priorPrec*prior.Mu + dataPrec*xbar)
	return Normal{Mu: postMean, Sigma: math.Sqrt(postVar)}
}

// NormalKnownVariancePosteriorData is like NormalKnownVariancePosterior but
// takes the raw data slice and known sampling variance sigma².
func NormalKnownVariancePosteriorData(prior Normal, sigma2 float64, data []float64) Normal {
	n := len(data)
	if n == 0 {
		return prior
	}
	var sum float64
	for _, x := range data {
		sum += x
	}
	return NormalKnownVariancePosterior(prior, sigma2, n, sum/float64(n))
}

// NormalKnownVariancePredictive returns the posterior predictive distribution
// for a single new observation with known sampling variance sigma², given the
// posterior over the mean. Its variance is the posterior variance plus sigma².
func NormalKnownVariancePredictive(post Normal, sigma2 float64) Normal {
	return Normal{Mu: post.Mu, Sigma: math.Sqrt(post.Sigma*post.Sigma + sigma2)}
}

// ------------------------------------------------------------------
// Inverse-Gamma update for the variance with a known mean
// ------------------------------------------------------------------

// InverseGammaVariancePosterior returns the posterior over the variance of a
// normal likelihood with known mean mu, given an Inverse-Gamma prior and the
// data. The update adds n/2 to the shape and ½Σ(xᵢ−mu)² to the scale.
func InverseGammaVariancePosterior(prior InverseGamma, mu float64, data []float64) InverseGamma {
	n := len(data)
	var ss float64
	for _, x := range data {
		d := x - mu
		ss += d * d
	}
	return InverseGamma{Shape: prior.Shape + float64(n)/2, Scale: prior.Scale + 0.5*ss}
}

// ------------------------------------------------------------------
// Normal-Inverse-Gamma: normal likelihood, unknown mean and variance
// ------------------------------------------------------------------

// NormalInverseGamma is the conjugate prior for a normal likelihood with both
// mean and variance unknown. Its hierarchy is σ² ~ InverseGamma(Alpha, Beta)
// and μ | σ² ~ Normal(Mu, σ²/Kappa).
type NormalInverseGamma struct {
	Mu    float64 // prior mean location
	Kappa float64 // prior pseudo-count on the mean
	Alpha float64 // Inverse-Gamma shape
	Beta  float64 // Inverse-Gamma scale
}

// NewNormalInverseGamma constructs a Normal-Inverse-Gamma prior, returning
// ErrParam if Kappa, Alpha or Beta is non-positive.
func NewNormalInverseGamma(mu, kappa, alpha, beta float64) (NormalInverseGamma, error) {
	if kappa <= 0 || alpha <= 0 || beta <= 0 {
		return NormalInverseGamma{}, ErrParam
	}
	return NormalInverseGamma{Mu: mu, Kappa: kappa, Alpha: alpha, Beta: beta}, nil
}

// NormalInverseGammaPosterior updates a Normal-Inverse-Gamma prior with the
// observed data and returns the posterior parameters.
func NormalInverseGammaPosterior(prior NormalInverseGamma, data []float64) NormalInverseGamma {
	n := float64(len(data))
	if n == 0 {
		return prior
	}
	var sum float64
	for _, x := range data {
		sum += x
	}
	xbar := sum / n
	var ss float64
	for _, x := range data {
		d := x - xbar
		ss += d * d
	}
	kappaN := prior.Kappa + n
	muN := (prior.Kappa*prior.Mu + sum) / kappaN
	alphaN := prior.Alpha + n/2
	betaN := prior.Beta + 0.5*ss + 0.5*prior.Kappa*n*(xbar-prior.Mu)*(xbar-prior.Mu)/kappaN
	return NormalInverseGamma{Mu: muN, Kappa: kappaN, Alpha: alphaN, Beta: betaN}
}

// MarginalVariance returns the marginal posterior over the variance, an
// Inverse-Gamma(Alpha, Beta) distribution.
func (p NormalInverseGamma) MarginalVariance() InverseGamma {
	return InverseGamma{Shape: p.Alpha, Scale: p.Beta}
}

// MarginalMean returns the marginal posterior over the mean, a location-scale
// Student-t with 2·Alpha degrees of freedom.
func (p NormalInverseGamma) MarginalMean() StudentT {
	scale := math.Sqrt(p.Beta / (p.Alpha * p.Kappa))
	return StudentT{Nu: 2 * p.Alpha, Loc: p.Mu, Scale: scale}
}

// Predictive returns the posterior predictive distribution for a single new
// observation, a location-scale Student-t with 2·Alpha degrees of freedom.
func (p NormalInverseGamma) Predictive() StudentT {
	scale := math.Sqrt(p.Beta * (p.Kappa + 1) / (p.Alpha * p.Kappa))
	return StudentT{Nu: 2 * p.Alpha, Loc: p.Mu, Scale: scale}
}

// ExpectedVariance returns E[σ²] under the marginal posterior, Beta/(Alpha−1)
// for Alpha > 1.
func (p NormalInverseGamma) ExpectedVariance() float64 {
	return p.MarginalVariance().Mean()
}

// ------------------------------------------------------------------
// Dirichlet-Multinomial conjugate update
// ------------------------------------------------------------------

// DirichletMultinomialPosterior returns the posterior Dirichlet after observing
// category counts. The counts slice must have the same length as the prior
// concentration vector; a length mismatch returns the prior unchanged.
func DirichletMultinomialPosterior(prior Dirichlet, counts []int) Dirichlet {
	if len(counts) != len(prior.Alpha) {
		return prior
	}
	out := make([]float64, len(prior.Alpha))
	for i := range out {
		out[i] = prior.Alpha[i] + float64(counts[i])
	}
	return Dirichlet{Alpha: out}
}

// DirichletCategoricalPosterior updates a Dirichlet prior with a sequence of
// observed category indices in [0,K) and returns the posterior.
func DirichletCategoricalPosterior(prior Dirichlet, categories []int) Dirichlet {
	out := make([]float64, len(prior.Alpha))
	copy(out, prior.Alpha)
	for _, c := range categories {
		if c >= 0 && c < len(out) {
			out[c]++
		}
	}
	return Dirichlet{Alpha: out}
}

// DirichletMultinomialPredictive returns the predictive category probabilities
// for the next single draw under the posterior, i.e. the posterior mean αᵢ/α₀.
func DirichletMultinomialPredictive(post Dirichlet) []float64 {
	return post.Mean()
}
