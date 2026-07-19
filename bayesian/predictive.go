package bayesian

// BetaBinomialPredictive returns the posterior predictive distribution for the
// number of successes in m future Binomial trials, given a Beta posterior.
func BetaBinomialPredictive(post Beta, m int) BetaBinomial {
	return BetaBinomial{N: m, Alpha: post.Alpha, Beta: post.Beta}
}

// GammaPoissonPredictive returns the posterior predictive distribution for the
// count of a single future Poisson observation, a Negative-Binomial with R =
// shape and P = rate/(rate+1).
func GammaPoissonPredictive(post Gamma) NegativeBinomial {
	return NegativeBinomial{R: post.Shape, P: post.Rate / (post.Rate + 1)}
}

// GammaPoissonPredictiveExposure returns the posterior predictive for the count
// observed over a future exposure t, a Negative-Binomial with R = shape and
// P = rate/(rate+t).
func GammaPoissonPredictiveExposure(post Gamma, t float64) NegativeBinomial {
	return NegativeBinomial{R: post.Shape, P: post.Rate / (post.Rate + t)}
}

// NormalNormalPredictive returns the posterior predictive distribution for a
// single new observation with known sampling variance sigma², given the
// posterior over the mean. It is an alias for NormalKnownVariancePredictive.
func NormalNormalPredictive(post Normal, sigma2 float64) Normal {
	return NormalKnownVariancePredictive(post, sigma2)
}

// PosteriorMean returns the posterior mean of a distribution exposing a Mean
// method. It is provided for symmetry with PosteriorVariance.
func PosteriorMean(d interface{ Mean() float64 }) float64 { return d.Mean() }

// PosteriorVariance returns the posterior variance of a distribution exposing a
// Variance method.
func PosteriorVariance(d interface{ Variance() float64 }) float64 { return d.Variance() }
