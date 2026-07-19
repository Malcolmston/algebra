// Package bayesian is a dependency-free toolkit for Bayesian inference written
// entirely with the Go standard library.
//
// It brings together the workhorses of applied Bayesian statistics behind a
// small, value-oriented API. Nothing samples or requires random numbers; every
// routine is deterministic and analytic (or, where a closed form does not
// exist, evaluated by a well-conditioned numerical method).
//
// # Distributions
//
// The package provides the continuous conjugate families [Beta], [Gamma],
// [InverseGamma], [Normal] and the location-scale [StudentT], the multivariate
// [Dirichlet], and the discrete [Poisson], [Binomial], [BetaBinomial] and
// [NegativeBinomial] laws. Each distribution offers the usual density/mass,
// cumulative, quantile, mean, variance and (where meaningful) entropy methods.
//
// # Conjugate updates
//
// Closed-form posterior updates are provided for the standard conjugate models:
// Beta-Bernoulli/Binomial ([BetaBernoulliPosterior], [BetaBinomialPosterior]),
// Gamma-Poisson ([GammaPoissonPosterior]), Normal-Normal with known variance
// ([NormalKnownVariancePosterior]), the Normal-Inverse-Gamma model for unknown
// mean and variance ([NormalInverseGammaPosterior]), Inverse-Gamma updates for
// a variance with known mean ([InverseGammaVariancePosterior]), and
// Dirichlet-Multinomial ([DirichletMultinomialPosterior]).
//
// # Predictions and summaries
//
// Posterior predictive distributions ([BetaBinomialPredictive],
// [GammaPoissonPredictive], [NormalNormalPredictive] and the Student-t
// predictive of the Normal-Inverse-Gamma model), equal-tailed and
// highest-density credible intervals ([EqualTailedInterval],
// [HighestDensityInterval]), and tail/interval probabilities round out the
// summary tooling.
//
// # Model comparison
//
// Log marginal likelihoods are available for each conjugate model
// (for example [LogMarginalLikelihoodBetaBinomial]), together with
// [BayesFactor], [PosteriorOdds], [PosteriorModelProbabilities] and a
// qualitative [BayesFactorInterpretation].
//
// # Classifiers and graphical models
//
// The package includes Gaussian, multinomial and Bernoulli naive-Bayes
// classifiers ([GaussianNB], [MultinomialNB], [BernoulliNB]) and a discrete
// Bayesian network ([BayesianNetwork]) supporting joint, marginal and
// conditional inference by variable elimination over its [Factor] potentials.
//
// # Special functions
//
// Underlying everything is a small set of exported special functions —
// [LogGamma], [LogBeta], [Digamma], [Trigamma], [RegularizedIncompleteBeta],
// [RegularizedGammaP] and their inverses — that may be useful on their own.
package bayesian
