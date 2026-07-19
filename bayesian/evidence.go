package bayesian

import (
	"math"
)

// LogMarginalLikelihoodBetaBinomial returns the log marginal likelihood (log
// evidence) of observing k successes in n Binomial trials under a Beta(a,b)
// prior on the success probability.
func LogMarginalLikelihoodBetaBinomial(prior Beta, k, n int) float64 {
	return LogChoose(n, k) +
		LogBeta(prior.Alpha+float64(k), prior.Beta+float64(n-k)) -
		LogBeta(prior.Alpha, prior.Beta)
}

// MarginalLikelihoodBetaBinomial returns the marginal likelihood (evidence) of
// k successes in n trials under a Beta prior.
func MarginalLikelihoodBetaBinomial(prior Beta, k, n int) float64 {
	return math.Exp(LogMarginalLikelihoodBetaBinomial(prior, k, n))
}

// LogMarginalLikelihoodGammaPoisson returns the log marginal likelihood of the
// observed Poisson counts under a Gamma(shape, rate) prior on the rate.
func LogMarginalLikelihoodGammaPoisson(prior Gamma, counts []int) float64 {
	n := len(counts)
	sum := 0
	var logFact float64
	for _, c := range counts {
		sum += c
		logFact += LogFactorial(c)
	}
	a, b := prior.Shape, prior.Rate
	fs := float64(sum)
	return -logFact +
		LogGamma(a+fs) - LogGamma(a) +
		a*math.Log(b) - (a+fs)*math.Log(b+float64(n))
}

// LogMarginalLikelihoodDirichletMultinomial returns the log marginal likelihood
// of the observed category counts under a Dirichlet prior.
func LogMarginalLikelihoodDirichletMultinomial(prior Dirichlet, counts []int) float64 {
	if len(counts) != len(prior.Alpha) {
		return math.NaN()
	}
	alpha0 := prior.Sum()
	total := 0
	logMult := 0.0
	var num, den float64
	for i, c := range counts {
		total += c
		logMult -= LogFactorial(c)
		num += LogGamma(prior.Alpha[i] + float64(c))
		den += LogGamma(prior.Alpha[i])
	}
	logMult += LogFactorial(total)
	return logMult + (num - den) - (LogGamma(alpha0+float64(total)) - LogGamma(alpha0))
}

// LogMarginalLikelihoodNormalKnownVariance returns the log marginal likelihood
// of the data under a normal likelihood with known sampling variance sigma²
// and a Normal prior on the mean. It is computed exactly via the sequential
// posterior-predictive factorization.
func LogMarginalLikelihoodNormalKnownVariance(prior Normal, sigma2 float64, data []float64) float64 {
	post := prior
	var ll float64
	for _, x := range data {
		pred := Normal{Mu: post.Mu, Sigma: math.Sqrt(post.Sigma*post.Sigma + sigma2)}
		ll += pred.LogPDF(x)
		post = NormalKnownVariancePosterior(post, sigma2, 1, x)
	}
	return ll
}

// LogMarginalLikelihoodNormalInverseGamma returns the log marginal likelihood
// of the data under a normal likelihood with unknown mean and variance and a
// Normal-Inverse-Gamma prior.
func LogMarginalLikelihoodNormalInverseGamma(prior NormalInverseGamma, data []float64) float64 {
	n := float64(len(data))
	if n == 0 {
		return 0
	}
	post := NormalInverseGammaPosterior(prior, data)
	// log p(D) = log Γ(αₙ) − log Γ(α₀) + α₀ log β₀ − αₙ log βₙ
	//            + ½(log κ₀ − log κₙ) − (n/2) log(2π)
	return LogGamma(post.Alpha) - LogGamma(prior.Alpha) +
		prior.Alpha*math.Log(prior.Beta) - post.Alpha*math.Log(post.Beta) +
		0.5*(math.Log(prior.Kappa)-math.Log(post.Kappa)) -
		(n/2)*math.Log(2*math.Pi)
}

// BayesFactor returns the Bayes factor B₁₂ = exp(logEvidence1 − logEvidence2)
// comparing model 1 against model 2 from their log marginal likelihoods.
func BayesFactor(logEvidence1, logEvidence2 float64) float64 {
	return math.Exp(logEvidence1 - logEvidence2)
}

// LogBayesFactor returns the log Bayes factor logEvidence1 − logEvidence2.
func LogBayesFactor(logEvidence1, logEvidence2 float64) float64 {
	return logEvidence1 - logEvidence2
}

// PosteriorOdds returns the posterior odds of model 1 versus model 2, the
// product of the Bayes factor and the prior odds.
func PosteriorOdds(bayesFactor, priorOdds float64) float64 {
	return bayesFactor * priorOdds
}

// PosteriorModelProbabilities returns the posterior probabilities of a set of
// models given their log marginal likelihoods and (unnormalized) prior
// probabilities. The two slices must have equal length; a mismatch returns nil.
func PosteriorModelProbabilities(logEvidences, priors []float64) []float64 {
	if len(logEvidences) != len(priors) || len(logEvidences) == 0 {
		return nil
	}
	logJoint := make([]float64, len(logEvidences))
	var priorSum float64
	for _, p := range priors {
		priorSum += p
	}
	for i := range logEvidences {
		logJoint[i] = logEvidences[i] + math.Log(priors[i]/priorSum)
	}
	norm := LogSumExp(logJoint)
	out := make([]float64, len(logEvidences))
	for i := range out {
		out[i] = math.Exp(logJoint[i] - norm)
	}
	return out
}

// BayesFactorInterpretation returns a qualitative label for the strength of
// evidence carried by a Bayes factor on the Jeffreys/Kass–Raftery scale
// ("negligible", "positive", "strong", or "decisive"). The magnitude is taken
// so that values below 1 are interpreted in favor of the alternative.
func BayesFactorInterpretation(bf float64) string {
	if bf < 1 {
		bf = 1 / bf
	}
	switch {
	case bf < 3.2:
		return "negligible"
	case bf < 10:
		return "positive"
	case bf < 100:
		return "strong"
	default:
		return "decisive"
	}
}
