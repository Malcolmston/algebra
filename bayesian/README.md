# bayesian — dependency-free toolkit for Bayesian inference written entirely with the Go standard

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/bayesian.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/bayesian)

Package bayesian is a dependency-free toolkit for Bayesian inference written
entirely with the Go standard library.

It brings together the workhorses of applied Bayesian statistics behind a
small, value-oriented API. Nothing samples or requires random numbers; every
routine is deterministic and analytic (or, where a closed form does not exist,
evaluated by a well-conditioned numerical method).

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/bayesian"
```

## Usage

This is the package's own `ExampleBayesianNetwork`, so it compiles and its output is
asserted on every `go test ./bayesian/`.

```go
bn := bs.NewBayesianNetwork()
	_ = bn.AddVariable("Disease", 2)
	_ = bn.AddVariable("Test", 2)
	_ = bn.SetCPT("Disease", nil, []float64{0.999, 0.001})
	_ = bn.SetCPT("Test", []string{"Disease"}, []float64{0.99, 0.01, 0.01, 0.99})
	p, err := bn.ConditionalProbability("Disease", 1, map[string]int{"Test": 1})
	if err != nil {
		panic(err)
	}
	fmt.Printf("%.4f\n", p)
```

```
0.0902
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func Accuracy(clf interface { Predict([]float64) (int, error) }, X [][]float64, y []int) (float64, error)` | Accuracy returns the fraction of samples in X whose predicted label from the given classifier matches the true label in y. |
| `func BayesFactor(logEvidence1, logEvidence2 float64) float64` | BayesFactor returns the Bayes factor B₁₂ = exp(logEvidence1 − logEvidence2) comparing model 1 against model 2 from their log marginal… |
| `func BayesFactorInterpretation(bf float64) string` | BayesFactorInterpretation returns a qualitative label for the strength of evidence carried by a Bayes factor on the Jeffreys/Kass–Raftery scale… |
| `func BetaBernoulliPredictiveProb(post Beta) float64` | BetaBernoulliPredictiveProb returns the probability that the next Bernoulli trial is a success under the posterior, which equals the posterior mean. |
| `func BetaFunc(a, b float64) float64` | BetaFunc returns the Beta function B(a,b) = Γ(a)Γ(b)/Γ(a+b) for a,b > 0. |
| `func Digamma(x float64) float64` | Digamma returns the digamma function ψ(x) = d/dx ln Γ(x). |
| `func DirichletMultinomialPredictive(post Dirichlet) []float64` | DirichletMultinomialPredictive returns the predictive category probabilities for the next single draw under the posterior, i.e. |
| `func ErfInv(x float64) float64` | ErfInv returns the inverse error function, the value y such that erf(y) = x for x in (−1,1). |
| `func GammaFunc(x float64) float64` | GammaFunc returns the Gamma function Γ(x). |
| `func InverseRegularizedGammaP(p, a float64) float64` | InverseRegularizedGammaP returns the value x such that RegularizedGammaP(a,x) = p, for a > 0 and p in [0,1]. |
| `func InverseRegularizedIncompleteBeta(p, a, b float64) float64` | InverseRegularizedIncompleteBeta returns the value x in [0,1] such that RegularizedIncompleteBeta(x,a,b) = p, i.e. |
| `func JensenShannonCategorical(p, q []float64) float64` | JensenShannonCategorical returns the Jensen-Shannon divergence in nats between two discrete probability vectors. |
| `func KLDivergenceBeta(p, q Beta) float64` | KLDivergenceBeta returns the Kullback-Leibler divergence D(p ‖ q) in nats between two Beta distributions. |
| `func KLDivergenceCategorical(p, q []float64) float64` | KLDivergenceCategorical returns the Kullback-Leibler divergence D(p ‖ q) in nats between two discrete distributions given as probability vectors of… |
| `func KLDivergenceDirichlet(p, q Dirichlet) float64` | KLDivergenceDirichlet returns the Kullback-Leibler divergence D(p ‖ q) in nats between two Dirichlet distributions over the same simplex. |
| `func KLDivergenceGamma(p, q Gamma) float64` | KLDivergenceGamma returns the Kullback-Leibler divergence D(p ‖ q) in nats between two Gamma distributions in the shape/rate parameterization. |
| `func KLDivergenceNormal(p, q Normal) float64` | KLDivergenceNormal returns the Kullback-Leibler divergence D(p ‖ q) in nats between two univariate normal distributions. |
| `func LogBayesFactor(logEvidence1, logEvidence2 float64) float64` | LogBayesFactor returns the log Bayes factor logEvidence1 − logEvidence2. |
| `func LogBeta(a, b float64) float64` | LogBeta returns the natural logarithm of the Beta function, ln B(a,b) = ln Γ(a) + ln Γ(b) − ln Γ(a+b), for a,b > 0. |
| `func LogChoose(n, k int) float64` | LogChoose returns the natural logarithm of the binomial coefficient ln C(n,k). |
| `func LogFactorial(n int) float64` | LogFactorial returns ln(n!) computed in the log-gamma domain as ln Γ(n+1). |
| `func LogGamma(x float64) float64` | LogGamma returns the natural logarithm of the absolute value of the Gamma function, ln\|Γ(x)\|. |
| `func LogMarginalLikelihoodBetaBinomial(prior Beta, k, n int) float64` | LogMarginalLikelihoodBetaBinomial returns the log marginal likelihood (log evidence) of observing k successes in n Binomial trials under a Beta(a,b)… |
| `func LogMarginalLikelihoodDirichletMultinomial(prior Dirichlet, counts []int) float64` | LogMarginalLikelihoodDirichletMultinomial returns the log marginal likelihood of the observed category counts under a Dirichlet prior. |
| `func LogMarginalLikelihoodGammaPoisson(prior Gamma, counts []int) float64` | LogMarginalLikelihoodGammaPoisson returns the log marginal likelihood of the observed Poisson counts under a Gamma(shape, rate) prior on the rate. |
| `func LogMarginalLikelihoodNormalInverseGamma(prior NormalInverseGamma, data []float64) float64` | LogMarginalLikelihoodNormalInverseGamma returns the log marginal likelihood of the data under a normal likelihood with unknown mean and variance and… |
| `func LogMarginalLikelihoodNormalKnownVariance(prior Normal, sigma2 float64, data []float64) float64` | LogMarginalLikelihoodNormalKnownVariance returns the log marginal likelihood of the data under a normal likelihood with known sampling variance… |
| `func LogMultinomialCoefficient(ks []int) float64` | LogMultinomialCoefficient returns the natural logarithm of the multinomial coefficient (Σkᵢ)! / Πkᵢ! for the non-negative counts ks. |
| `func LogSumExp(xs []float64) float64` | LogSumExp returns ln(Σ exp(xᵢ)) computed in a numerically stable way by factoring out the maximum element. |
| `func MarginalLikelihoodBetaBinomial(prior Beta, k, n int) float64` | MarginalLikelihoodBetaBinomial returns the marginal likelihood (evidence) of k successes in n trials under a Beta prior. |
| `func PosteriorMean(d interface{ Mean() float64 }) float64` | PosteriorMean returns the posterior mean of a distribution exposing a Mean method. |
| `func PosteriorModelProbabilities(logEvidences, priors []float64) []float64` | PosteriorModelProbabilities returns the posterior probabilities of a set of models given their log marginal likelihoods and (unnormalized) prior… |
| `func PosteriorOdds(bayesFactor, priorOdds float64) float64` | PosteriorOdds returns the posterior odds of model 1 versus model 2, the product of the Bayes factor and the prior odds. |
| `func PosteriorVariance(d interface{ Variance() float64 }) float64` | PosteriorVariance returns the posterior variance of a distribution exposing a Variance method. |
| `func ProbabilityBetaExceedsBeta(a, b Beta) float64` | ProbabilityBetaExceedsBeta returns the probability P(X > Y) that a draw from Beta a exceeds an independent draw from Beta b, computed by numerical… |
| `func ProbabilityBetaGreater(d Beta, threshold float64) float64` | ProbabilityBetaGreater returns the posterior probability that a Beta random variable exceeds threshold, P(θ > t). |
| `func ProbabilityBetween(d interface{ CDF(float64) float64 }, lo, hi float64) float64` | ProbabilityBetween returns P(lo < X ≤ hi) for a distribution exposing a CDF. |
| `func ProbabilityGreater(d interface{ CDF(float64) float64 }, threshold float64) float64` | ProbabilityGreater returns P(X > threshold) for a distribution exposing a CDF method. |
| `func RegularizedGammaP(a, x float64) float64` | RegularizedGammaP returns the regularized lower incomplete gamma function P(a,x) = γ(a,x)/Γ(a) for a > 0 and x >= 0. |
| `func RegularizedGammaQ(a, x float64) float64` | RegularizedGammaQ returns the regularized upper incomplete gamma function Q(a,x) = 1 − P(a,x) = Γ(a,x)/Γ(a) for a > 0 and x >= 0. |
| `func RegularizedIncompleteBeta(x, a, b float64) float64` | RegularizedIncompleteBeta returns the regularized incomplete beta function Iₓ(a,b) for x in [0,1] and a,b > 0. |
| `func StdNormalCDF(z float64) float64` | StdNormalCDF returns the standard normal cumulative distribution function Φ(z) = ½(1 + erf(z/√2)). |
| `func StdNormalPDF(z float64) float64` | StdNormalPDF returns the standard normal probability density φ(z). |
| `func StdNormalQuantile(p float64) float64` | StdNormalQuantile returns the p-quantile of the standard normal distribution, the value z such that Φ(z) = p for p in (0,1). |
| `func Trigamma(x float64) float64` | Trigamma returns the trigamma function ψ₁(x) = d²/dx² ln Γ(x), the derivative of the digamma function. |

### Types

| Type | What it is |
| --- | --- |
| `BayesianNetwork` | BayesianNetwork is a discrete directed graphical model. |
| `BernoulliNB` | BernoulliNB is a Bernoulli naive-Bayes classifier for binary feature vectors; each feature contributes both when present and, via its complement,… |
| `Beta` | Beta is a Beta distribution with shape parameters Alpha, Beta > 0. |
| `BetaBinomial` | BetaBinomial is the Beta-Binomial distribution: the number of successes in N trials where the success probability is drawn from Beta(Alpha, Beta). |
| `Binomial` | Binomial is a Binomial distribution with N ≥ 0 trials and success probability P in [0,1]. |
| `Dirichlet` | Dirichlet is a Dirichlet distribution over the probability simplex with concentration vector Alpha, all entries strictly positive. |
| `Factor` | Factor is a discrete potential over a set of named categorical variables. |
| `Gamma` | Gamma is a Gamma distribution with Shape > 0 and Rate > 0 (rate = 1/scale). |
| `GaussianNB` | GaussianNB is a Gaussian naive-Bayes classifier: each continuous feature is modeled as a class-conditional normal distribution with features assumed… |
| `Interval` | Interval is a closed real interval [Lower, Upper]. |
| `InverseGamma` | InverseGamma is an Inverse-Gamma distribution with Shape > 0 and Scale > 0. |
| `MultinomialNB` | MultinomialNB is a multinomial naive-Bayes classifier suitable for count features such as word frequencies, using additive (Laplace/Lidstone)… |
| `NegativeBinomial` | NegativeBinomial is parameterized by a real number R > 0 of "failures" and a success probability P in (0,1); the mass is over the number of counts k… |
| `Normal` | Normal is a univariate normal (Gaussian) distribution with mean Mu and standard deviation Sigma > 0. |
| `NormalInverseGamma` | NormalInverseGamma is the conjugate prior for a normal likelihood with both mean and variance unknown. |
| `Poisson` | Poisson is a Poisson distribution with rate Lambda ≥ 0. |
| `StudentT` | StudentT is a location-scale Student's t distribution with Nu > 0 degrees of freedom, location Loc and scale Scale > 0. |

<details>
<summary><code>BayesianNetwork</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewBayesianNetwork() *BayesianNetwork` | NewBayesianNetwork returns an empty Bayesian network. |
| `func (bn *BayesianNetwork) AddVariable(name string, card int) error` | AddVariable declares a categorical variable with the given cardinality. |
| `func (bn *BayesianNetwork) CPT(name string) (Factor, bool)` | CPT returns the conditional probability factor for name and reports whether it has been set. |
| `func (bn *BayesianNetwork) Cardinality(name string) int` | Cardinality returns the cardinality of the named variable, or 0 if unknown. |
| `func (bn *BayesianNetwork) Conditional(query []string, evidence map[string]int) (Factor, error)` | Conditional returns the posterior distribution over the query variables given the evidence, as a normalized Factor. |
| `func (bn *BayesianNetwork) ConditionalProbability(variable string, value int, evidence map[string]int) (float64, error)` | ConditionalProbability returns P(variable = value \| evidence), a scalar posterior for a single variable and category. |
| `func (bn *BayesianNetwork) EvidenceProbability(evidence map[string]int) (float64, error)` | EvidenceProbability returns P(evidence), the probability of the observed assignment obtained by summing the joint distribution consistent with it. |
| `func (bn *BayesianNetwork) JointProbability(assign map[string]int) float64` | JointProbability returns P(assignment) for a full assignment of every variable, the product of each CPT evaluated at the assignment. |
| `func (bn *BayesianNetwork) Marginal(query ...string) (Factor, error)` | Marginal returns the marginal distribution over the query variables as a normalized Factor. |
| `func (bn *BayesianNetwork) MarginalProbability(variable string, value int) (float64, error)` | MarginalProbability returns P(variable = value), a scalar marginal for a single variable and category. |
| `func (bn *BayesianNetwork) Parents(name string) []string` | Parents returns the parent variables of name in the order they were declared. |
| `func (bn *BayesianNetwork) SetCPT(name string, parents []string, table []float64) error` | SetCPT installs the conditional probability table P(name \| parents). |
| `func (bn *BayesianNetwork) Variables() []string` | Variables returns the network's variable names in declaration order. |

</details>

<details>
<summary><code>BernoulliNB</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewBernoulliNB(alpha, threshold float64) *BernoulliNB` | NewBernoulliNB returns an unfitted Bernoulli naive-Bayes classifier. |
| `func (b *BernoulliNB) Classes() []int` | Classes returns a copy of the sorted class labels learned by Fit. |
| `func (b *BernoulliNB) Fit(X [][]float64, y []int) error` | Fit estimates class priors and smoothed feature-presence probabilities. |
| `func (b *BernoulliNB) JointLogLikelihood(x []float64) ([]float64, error)` | JointLogLikelihood returns the unnormalized log joint probabilities for the binary feature vector x under each class. |
| `func (b *BernoulliNB) Predict(x []float64) (int, error)` | Predict returns the most probable class label for x. |
| `func (b *BernoulliNB) PredictBatch(X [][]float64) ([]int, error)` | PredictBatch returns predictions for every row of X. |
| `func (b *BernoulliNB) PredictProba(x []float64) ([]float64, error)` | PredictProba returns the posterior class probabilities for x. |
| `func (b *BernoulliNB) Priors() []float64` | Priors returns a copy of the estimated class prior probabilities, aligned with Classes. |

</details>

<details>
<summary><code>Beta</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BetaBernoulliPosterior(prior Beta, successes, failures int) Beta` | BetaBernoulliPosterior returns the posterior Beta(α+successes, β+failures) after observing the given counts of Bernoulli successes and failures. |
| `func BetaBernoulliSequentialPosterior(prior Beta, outcomes []bool) Beta` | BetaBernoulliSequentialPosterior applies a sequence of Bernoulli outcomes (true = success) to the prior and returns the resulting posterior. |
| `func BetaBinomialPosterior(prior Beta, k, n int) Beta` | BetaBinomialPosterior returns the posterior after observing k successes in n Binomial trials. |
| `func BetaGeometricPosterior(prior Beta, failures []int) Beta` | BetaGeometricPosterior returns the posterior over the success probability of a geometric likelihood with a Beta prior. |
| `func BetaNegativeBinomialPosterior(prior Beta, r int, failures []int) Beta` | BetaNegativeBinomialPosterior returns the posterior over the success probability of a negative-binomial likelihood with a Beta prior, where each… |
| `func HaldaneBetaPrior() Beta` | HaldaneBetaPrior returns Beta(ε,ε) with a tiny ε, an approximation of the improper Haldane prior Beta(0,0) whose posterior mean equals the sample… |
| `func JeffreysBetaPrior() Beta` | JeffreysBetaPrior returns Beta(½,½), the Jeffreys prior for a Bernoulli success probability. |
| `func NewBeta(alpha, beta float64) (Beta, error)` | NewBeta constructs a Beta distribution, returning ErrParam if either shape is non-positive. |
| `func UniformBetaPrior() Beta` | UniformBetaPrior returns Beta(1,1), the uniform prior on a probability. |
| `func (d Beta) CDF(x float64) float64` | CDF returns the cumulative distribution function P(X ≤ x). |
| `func (d Beta) Entropy() float64` | Entropy returns the differential entropy of the distribution in nats. |
| `func (d Beta) LogPDF(x float64) float64` | LogPDF returns the natural logarithm of the density at x in (0,1). |
| `func (d Beta) Mean() float64` | Mean returns the expected value α/(α+β). |
| `func (d Beta) MeanLog() float64` | MeanLog returns E[ln X] = ψ(α) − ψ(α+β). |
| `func (d Beta) Mode() float64` | Mode returns the mode (α−1)/(α+β−2). |
| `func (d Beta) PDF(x float64) float64` | PDF returns the probability density at x in [0,1]. |
| `func (d Beta) Quantile(p float64) float64` | Quantile returns the p-quantile (inverse CDF) for p in [0,1]. |
| `func (d Beta) Skewness() float64` | Skewness returns the distribution's skewness. |
| `func (d Beta) StdDev() float64` | StdDev returns the standard deviation, the square root of the variance. |
| `func (d Beta) Variance() float64` | Variance returns αβ/((α+β)²(α+β+1)). |

</details>

<details>
<summary><code>BetaBinomial</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BetaBinomialPredictive(post Beta, m int) BetaBinomial` | BetaBinomialPredictive returns the posterior predictive distribution for the number of successes in m future Binomial trials, given a Beta posterior. |
| `func NewBetaBinomial(n int, alpha, beta float64) (BetaBinomial, error)` | NewBetaBinomial constructs a Beta-Binomial distribution, returning ErrParam for invalid parameters. |
| `func (d BetaBinomial) CDF(k int) float64` | CDF returns P(X ≤ k) computed by summation over the support. |
| `func (d BetaBinomial) LogPMF(k int) float64` | LogPMF returns the natural logarithm of the mass at k. |
| `func (d BetaBinomial) Mean() float64` | Mean returns N·α/(α+β). |
| `func (d BetaBinomial) PMF(k int) float64` | PMF returns the probability mass P(X = k) for k in [0,N]. |
| `func (d BetaBinomial) Variance() float64` | Variance returns the variance of the Beta-Binomial distribution. |

</details>

<details>
<summary><code>Binomial</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewBinomial(n int, p float64) (Binomial, error)` | NewBinomial constructs a Binomial distribution, returning ErrParam for an out-of-range probability or negative N. |
| `func (d Binomial) CDF(k int) float64` | CDF returns P(X ≤ k). |
| `func (d Binomial) LogPMF(k int) float64` | LogPMF returns the natural logarithm of the mass at k. |
| `func (d Binomial) Mean() float64` | Mean returns N·P. |
| `func (d Binomial) PMF(k int) float64` | PMF returns the probability mass P(X = k). |
| `func (d Binomial) Variance() float64` | Variance returns N·P·(1−P). |

</details>

<details>
<summary><code>Dirichlet</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DirichletCategoricalPosterior(prior Dirichlet, categories []int) Dirichlet` | DirichletCategoricalPosterior updates a Dirichlet prior with a sequence of observed category indices in [0,K) and returns the posterior. |
| `func DirichletMultinomialPosterior(prior Dirichlet, counts []int) Dirichlet` | DirichletMultinomialPosterior returns the posterior Dirichlet after observing category counts. |
| `func JeffreysDirichletPrior(k int) Dirichlet` | JeffreysDirichletPrior returns the Dirichlet(½, …, ½) prior over k categories, the Jeffreys prior for a categorical likelihood. |
| `func NewDirichlet(alpha []float64) (Dirichlet, error)` | NewDirichlet constructs a Dirichlet distribution from a concentration vector, returning ErrParam if it is empty or has a non-positive entry. |
| `func UniformDirichletPrior(k int) Dirichlet` | UniformDirichletPrior returns the Dirichlet(1, …, 1) prior over k categories, which is uniform on the simplex. |
| `func (d Dirichlet) Covariance() [][]float64` | Covariance returns the K×K covariance matrix of the distribution as a slice of rows. |
| `func (d Dirichlet) Dim() int` | Dim returns the number of categories. |
| `func (d Dirichlet) Entropy() float64` | Entropy returns the differential entropy of the distribution in nats. |
| `func (d Dirichlet) LogPDF(x []float64) float64` | LogPDF returns the natural logarithm of the density at the point x on the simplex. |
| `func (d Dirichlet) Marginal(i int) Beta` | Marginal returns the marginal distribution of component i, which is Beta(αᵢ, α₀−αᵢ). |
| `func (d Dirichlet) Mean() []float64` | Mean returns the mean vector αᵢ/α₀. |
| `func (d Dirichlet) MeanLog() []float64` | MeanLog returns the vector E[ln Xᵢ] = ψ(αᵢ) − ψ(α₀). |
| `func (d Dirichlet) Mode() []float64` | Mode returns the mode vector (αᵢ−1)/(α₀−K); it is defined only when every αᵢ > 1, otherwise the returned slice is nil. |
| `func (d Dirichlet) PDF(x []float64) float64` | PDF returns the density at the point x on the simplex. |
| `func (d Dirichlet) Sum() float64` | Sum returns the total concentration α₀ = Σαᵢ. |
| `func (d Dirichlet) Variance() []float64` | Variance returns the per-component marginal variances. |

</details>

<details>
<summary><code>Factor</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewFactor(vars []string, card []int, table []float64) (Factor, error)` | NewFactor constructs a Factor from parallel variable-name and cardinality slices and a value table. |
| `func (f Factor) Get(assign map[string]int) float64` | Get returns the factor value for a full or partial assignment given as a map from variable name to category index. |
| `func (f Factor) Marginalize(name string) Factor` | Marginalize returns the factor obtained by summing out the named variable. |
| `func (a Factor) Multiply(b Factor) Factor` | Multiply returns the factor product of a and b, whose scope is the union of their variables. |
| `func (f Factor) Normalize() Factor` | Normalize returns a copy of the factor scaled so its entries sum to one. |
| `func (f Factor) Reduce(evidence map[string]int) Factor` | Reduce returns the factor obtained by fixing the variables in evidence to the given category indices and dropping them from scope. |
| `func (f Factor) Scope() []string` | Scope returns a sorted copy of the factor's variable names. |
| `func (f Factor) Size() int` | Size returns the number of entries in the factor's table. |
| `func (f Factor) Sum() float64` | Sum returns the total mass of the factor. |

</details>

<details>
<summary><code>Gamma</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func GammaExponentialPosterior(prior Gamma, sumX float64, n int) Gamma` | GammaExponentialPosterior returns the posterior for the rate of an exponential likelihood with a Gamma prior, after observing n samples whose sum is… |
| `func GammaPoissonPosterior(prior Gamma, sumCounts, n int) Gamma` | GammaPoissonPosterior returns the posterior Gamma(shape+Σcounts, rate+n) after observing n Poisson counts whose total is sumCounts. |
| `func GammaPoissonPosteriorData(prior Gamma, counts []int) Gamma` | GammaPoissonPosteriorData returns the Gamma-Poisson posterior for a slice of observed non-negative integer counts. |
| `func GammaPoissonPosteriorExposure(prior Gamma, sumCounts int, exposure float64) Gamma` | GammaPoissonPosteriorExposure returns the Gamma-Poisson posterior when the counts are observed over a total exposure (time, area or population)… |
| `func GammaPrecisionPosterior(prior Gamma, mu float64, data []float64) Gamma` | GammaPrecisionPosterior returns the posterior over the precision (the reciprocal of the variance) of a normal likelihood with known mean mu, given a… |
| `func JeffreysGammaPrior() Gamma` | JeffreysGammaPrior returns Gamma(½, ε) with a tiny rate, an approximation of the improper Jeffreys prior λ^(-1/2) for a Poisson rate. |
| `func NewGamma(shape, rate float64) (Gamma, error)` | NewGamma constructs a Gamma distribution in the shape/rate parameterization, returning ErrParam if either parameter is non-positive. |
| `func (d Gamma) CDF(x float64) float64` | CDF returns the cumulative distribution function P(X ≤ x). |
| `func (d Gamma) Entropy() float64` | Entropy returns the differential entropy in nats. |
| `func (d Gamma) LogPDF(x float64) float64` | LogPDF returns the natural logarithm of the density at x > 0. |
| `func (d Gamma) Mean() float64` | Mean returns Shape/Rate. |
| `func (d Gamma) MeanLog() float64` | MeanLog returns E[ln X] = ψ(Shape) − ln(Rate). |
| `func (d Gamma) Mode() float64` | Mode returns (Shape−1)/Rate for Shape ≥ 1, and 0 otherwise. |
| `func (d Gamma) PDF(x float64) float64` | PDF returns the density at x ≥ 0. |
| `func (d Gamma) Quantile(p float64) float64` | Quantile returns the p-quantile for p in [0,1]. |
| `func (d Gamma) Scale() float64` | Scale returns the scale parameter 1/Rate. |
| `func (d Gamma) Skewness() float64` | Skewness returns 2/√Shape. |
| `func (d Gamma) StdDev() float64` | StdDev returns the standard deviation. |
| `func (d Gamma) Variance() float64` | Variance returns Shape/Rate². |

</details>

<details>
<summary><code>GaussianNB</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewGaussianNB() *GaussianNB` | NewGaussianNB returns an unfitted Gaussian naive-Bayes classifier. |
| `func (g *GaussianNB) Classes() []int` | Classes returns a copy of the sorted class labels learned by Fit. |
| `func (g *GaussianNB) Fit(X [][]float64, y []int) error` | Fit estimates class priors and per-feature class-conditional means and variances from the training matrix X (rows are samples) and labels y. |
| `func (g *GaussianNB) JointLogLikelihood(x []float64) ([]float64, error)` | JointLogLikelihood returns, for the sample x, the unnormalized log joint probabilities log P(class) + log P(x \| class) for each class. |
| `func (g *GaussianNB) Predict(x []float64) (int, error)` | Predict returns the most probable class label for the sample x. |
| `func (g *GaussianNB) PredictBatch(X [][]float64) ([]int, error)` | PredictBatch returns predictions for every row of X. |
| `func (g *GaussianNB) PredictProba(x []float64) ([]float64, error)` | PredictProba returns the posterior class probabilities for the sample x. |
| `func (g *GaussianNB) Priors() []float64` | Priors returns a copy of the estimated class prior probabilities, aligned with Classes. |

</details>

<details>
<summary><code>Interval</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BetaCredibleInterval(d Beta, level float64) Interval` | BetaCredibleInterval returns the equal-tailed credible interval of the given mass for a Beta distribution. |
| `func BetaHDI(d Beta, level float64) Interval` | BetaHDI returns the highest-density credible interval of the given mass for a Beta distribution. |
| `func EqualTailedInterval(d quantiler, level float64) Interval` | EqualTailedInterval returns the equal-tailed credible interval of mass level (for example 0.95) for any distribution with a Quantile method: it cuts… |
| `func GammaCredibleInterval(d Gamma, level float64) Interval` | GammaCredibleInterval returns the equal-tailed credible interval of the given mass for a Gamma distribution. |
| `func GammaHDI(d Gamma, level float64) Interval` | GammaHDI returns the highest-density credible interval of the given mass for a Gamma distribution. |
| `func HighestDensityInterval(d densityQuantiler, level float64) Interval` | HighestDensityInterval returns the highest-density credible interval of the given mass for a unimodal continuous distribution exposing PDF and… |
| `func InverseGammaCredibleInterval(d InverseGamma, level float64) Interval` | InverseGammaCredibleInterval returns the equal-tailed credible interval of the given mass for an Inverse-Gamma distribution. |
| `func InverseGammaHDI(d InverseGamma, level float64) Interval` | InverseGammaHDI returns the highest-density credible interval of the given mass for an Inverse-Gamma distribution. |
| `func NormalCredibleInterval(d Normal, level float64) Interval` | NormalCredibleInterval returns the equal-tailed credible interval of the given mass for a Normal distribution. |
| `func NormalHDI(d Normal, level float64) Interval` | NormalHDI returns the highest-density credible interval of the given mass for a Normal distribution. |
| `func StudentTCredibleInterval(d StudentT, level float64) Interval` | StudentTCredibleInterval returns the equal-tailed credible interval of the given mass for a Student-t distribution. |
| `func StudentTHDI(d StudentT, level float64) Interval` | StudentTHDI returns the highest-density credible interval of the given mass for a Student-t distribution. |
| `func (iv Interval) Contains(x float64) bool` | Contains reports whether x lies within the closed interval. |
| `func (iv Interval) Width() float64` | Width returns Upper − Lower. |

</details>

<details>
<summary><code>InverseGamma</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func InverseGammaVariancePosterior(prior InverseGamma, mu float64, data []float64) InverseGamma` | InverseGammaVariancePosterior returns the posterior over the variance of a normal likelihood with known mean mu, given an Inverse-Gamma prior and the… |
| `func NewInverseGamma(shape, scale float64) (InverseGamma, error)` | NewInverseGamma constructs an Inverse-Gamma distribution, returning ErrParam for non-positive parameters. |
| `func (d InverseGamma) CDF(x float64) float64` | CDF returns the cumulative distribution function P(X ≤ x). |
| `func (d InverseGamma) LogPDF(x float64) float64` | LogPDF returns the natural logarithm of the density at x > 0. |
| `func (d InverseGamma) Mean() float64` | Mean returns Scale/(Shape−1) for Shape > 1, else +Inf. |
| `func (d InverseGamma) Mode() float64` | Mode returns Scale/(Shape+1). |
| `func (d InverseGamma) PDF(x float64) float64` | PDF returns the density at x > 0. |
| `func (d InverseGamma) Quantile(p float64) float64` | Quantile returns the p-quantile for p in [0,1]. |
| `func (d InverseGamma) Variance() float64` | Variance returns the variance for Shape > 2, else +Inf. |

</details>

<details>
<summary><code>MultinomialNB</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewMultinomialNB(alpha float64) *MultinomialNB` | NewMultinomialNB returns an unfitted multinomial naive-Bayes classifier with the given additive smoothing parameter alpha (alpha = 1 is Laplace… |
| `func (m *MultinomialNB) Classes() []int` | Classes returns a copy of the sorted class labels learned by Fit. |
| `func (m *MultinomialNB) Fit(X [][]float64, y []int) error` | Fit estimates class priors and smoothed per-feature log probabilities from the count matrix X and labels y. |
| `func (m *MultinomialNB) JointLogLikelihood(x []float64) ([]float64, error)` | JointLogLikelihood returns the unnormalized log joint probabilities for the count vector x under each class. |
| `func (m *MultinomialNB) Predict(x []float64) (int, error)` | Predict returns the most probable class label for the count vector x. |
| `func (m *MultinomialNB) PredictBatch(X [][]float64) ([]int, error)` | PredictBatch returns predictions for every row of X. |
| `func (m *MultinomialNB) PredictProba(x []float64) ([]float64, error)` | PredictProba returns the posterior class probabilities for the count vector x. |
| `func (m *MultinomialNB) Priors() []float64` | Priors returns a copy of the estimated class prior probabilities, aligned with Classes. |

</details>

<details>
<summary><code>NegativeBinomial</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func GammaPoissonPredictive(post Gamma) NegativeBinomial` | GammaPoissonPredictive returns the posterior predictive distribution for the count of a single future Poisson observation, a Negative-Binomial with R… |
| `func GammaPoissonPredictiveExposure(post Gamma, t float64) NegativeBinomial` | GammaPoissonPredictiveExposure returns the posterior predictive for the count observed over a future exposure t, a Negative-Binomial with R = shape… |
| `func NewNegativeBinomial(r, p float64) (NegativeBinomial, error)` | NewNegativeBinomial constructs a Negative-Binomial distribution, returning ErrParam for invalid parameters. |
| `func (d NegativeBinomial) CDF(k int) float64` | CDF returns P(X ≤ k) via the regularized incomplete beta function. |
| `func (d NegativeBinomial) LogPMF(k int) float64` | LogPMF returns the natural logarithm of the mass at k ≥ 0 using the generalized (real-valued R) binomial coefficient. |
| `func (d NegativeBinomial) Mean() float64` | Mean returns R·(1−P)/P. |
| `func (d NegativeBinomial) PMF(k int) float64` | PMF returns the probability mass P(X = k) for k ≥ 0. |
| `func (d NegativeBinomial) Variance() float64` | Variance returns R·(1−P)/P². |

</details>

<details>
<summary><code>Normal</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewNormal(mu, sigma float64) (Normal, error)` | NewNormal constructs a Normal distribution, returning ErrParam if Sigma ≤ 0. |
| `func NormalKnownVariancePosterior(prior Normal, sigma2 float64, n int, xbar float64) Normal` | NormalKnownVariancePosterior returns the posterior over the mean of a normal likelihood with known sampling variance sigma², given a normal prior… |
| `func NormalKnownVariancePosteriorData(prior Normal, sigma2 float64, data []float64) Normal` | NormalKnownVariancePosteriorData is like NormalKnownVariancePosterior but takes the raw data slice and known sampling variance sigma². |
| `func NormalKnownVariancePredictive(post Normal, sigma2 float64) Normal` | NormalKnownVariancePredictive returns the posterior predictive distribution for a single new observation with known sampling variance sigma², given… |
| `func NormalNormalPredictive(post Normal, sigma2 float64) Normal` | NormalNormalPredictive returns the posterior predictive distribution for a single new observation with known sampling variance sigma², given the… |
| `func (d Normal) CDF(x float64) float64` | CDF returns the cumulative distribution function P(X ≤ x). |
| `func (d Normal) Entropy() float64` | Entropy returns the differential entropy in nats. |
| `func (d Normal) LogPDF(x float64) float64` | LogPDF returns the natural logarithm of the density at x. |
| `func (d Normal) Mean() float64` | Mean returns Mu. |
| `func (d Normal) PDF(x float64) float64` | PDF returns the density at x. |
| `func (d Normal) Quantile(p float64) float64` | Quantile returns the p-quantile for p in [0,1]. |
| `func (d Normal) StdDev() float64` | StdDev returns Sigma. |
| `func (d Normal) Variance() float64` | Variance returns Sigma². |

</details>

<details>
<summary><code>NormalInverseGamma</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewNormalInverseGamma(mu, kappa, alpha, beta float64) (NormalInverseGamma, error)` | NewNormalInverseGamma constructs a Normal-Inverse-Gamma prior, returning ErrParam if Kappa, Alpha or Beta is non-positive. |
| `func NormalInverseGammaPosterior(prior NormalInverseGamma, data []float64) NormalInverseGamma` | NormalInverseGammaPosterior updates a Normal-Inverse-Gamma prior with the observed data and returns the posterior parameters. |
| `func (p NormalInverseGamma) ExpectedVariance() float64` | ExpectedVariance returns E[σ²] under the marginal posterior, Beta/(Alpha−1) for Alpha > 1. |
| `func (p NormalInverseGamma) MarginalMean() StudentT` | MarginalMean returns the marginal posterior over the mean, a location-scale Student-t with 2·Alpha degrees of freedom. |
| `func (p NormalInverseGamma) MarginalVariance() InverseGamma` | MarginalVariance returns the marginal posterior over the variance, an Inverse-Gamma(Alpha, Beta) distribution. |
| `func (p NormalInverseGamma) Predictive() StudentT` | Predictive returns the posterior predictive distribution for a single new observation, a location-scale Student-t with 2·Alpha degrees of freedom. |

</details>

<details>
<summary><code>Poisson</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewPoisson(lambda float64) (Poisson, error)` | NewPoisson constructs a Poisson distribution, returning ErrParam for a negative rate. |
| `func (d Poisson) CDF(k int) float64` | CDF returns P(X ≤ k). |
| `func (d Poisson) LogPMF(k int) float64` | LogPMF returns the natural logarithm of the mass at k ≥ 0. |
| `func (d Poisson) Mean() float64` | Mean returns Lambda. |
| `func (d Poisson) PMF(k int) float64` | PMF returns the probability mass P(X = k) for integer k ≥ 0. |
| `func (d Poisson) Variance() float64` | Variance returns Lambda. |

</details>

<details>
<summary><code>StudentT</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewStudentT(nu, loc, scale float64) (StudentT, error)` | NewStudentT constructs a location-scale Student-t distribution, returning ErrParam for non-positive Nu or Scale. |
| `func (d StudentT) CDF(x float64) float64` | CDF returns the cumulative distribution function P(X ≤ x). |
| `func (d StudentT) LogPDF(x float64) float64` | LogPDF returns the natural logarithm of the density at x. |
| `func (d StudentT) Mean() float64` | Mean returns Loc for Nu > 1, else NaN. |
| `func (d StudentT) PDF(x float64) float64` | PDF returns the density at x. |
| `func (d StudentT) Quantile(p float64) float64` | Quantile returns the p-quantile for p in [0,1]. |
| `func (d StudentT) StdDev() float64` | StdDev returns the standard deviation. |
| `func (d StudentT) Variance() float64` | Variance returns the variance for Nu > 2, +Inf for 1 < Nu ≤ 2, else NaN. |

</details>

### Variables

`ErrFactor`, `ErrNetwork`, `ErrNotFitted`, `ErrParam`, `ErrShape`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/bayesian).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
