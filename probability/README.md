# probability — provides tools for working with discrete probability spaces, random variables

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/probability.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/probability)

Package probability provides tools for working with discrete probability
spaces, random variables, and finite-state Markov chains using only the Go
standard library.

The central type is `Distribution`, a finitely-supported discrete probability
distribution over real-valued outcomes. It offers the full complement of
summary statistics (expectation, variance, higher moments, skewness, excess
kurtosis, entropy), the standard generating functions (moment-generating,
cumulant-generating, probability-generating and the characteristic function),
cumulative and quantile queries, random-variable transforms (affine maps and
arbitrary functions), and convolution of independent sums.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/probability"
```

## Usage

This is the package's own `ExampleBayesPosterior`, so it compiles and its output is
asserted on every `go test ./probability/`.

```go
post, err := probability.BayesPosterior(
		[]float64{1.0 / 3, 1.0 / 3, 1.0 / 3},
		[]float64{0.1, 0.3, 0.6},
	)
	if err != nil {
		panic(err)
	}
	for _, p := range post {
		fmt.Printf("%.4f ", p)
	}
	fmt.Println()
```

```
0.1000 0.3000 0.6000
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AreIndependent(pA, pB, pAandB float64) bool` | AreIndependent reports whether events A and B are independent, i.e. |
| `func Bayes(pBgivenA, pA, pB float64) (float64, error)` | Bayes returns the posterior probability P(A \| B) via Bayes' theorem, P(B \| A)·P(A) / P(B). |
| `func BayesPosterior(priors, likelihoods []float64) ([]float64, error)` | BayesPosterior returns the full posterior distribution over a partition of hypotheses {A_i}: posterior`i` = priors`i`·likelihoods`i` / Σ_k… |
| `func Complement(p float64) float64` | Complement returns the probability of the complementary event, 1 - p. |
| `func ConditionalProbability(pAandB, pB float64) (float64, error)` | ConditionalProbability returns P(A \| B) = P(A ∩ B) / P(B). |
| `func JointProbabilityIndependent(pA, pB float64) float64` | JointProbabilityIndependent returns P(A ∩ B) = P(A)·P(B) under the assumption that A and B are independent. |
| `func OddsFromProbability(p float64) float64` | OddsFromProbability converts a probability p in [0, 1) to odds p / (1 - p). |
| `func ProbabilityFromOdds(o float64) float64` | ProbabilityFromOdds converts non-negative odds o to the probability o / (1 + o). |
| `func TotalProbability(priors, likelihoods []float64) (float64, error)` | TotalProbability returns P(B) = Σ_i priors`i`·likelihoods`i` via the law of total probability, where priors is a partition {A_i} with P(A_i) =… |
| `func UnionProbability(pA, pB, pAandB float64) float64` | UnionProbability returns P(A ∪ B) = P(A) + P(B) - P(A ∩ B) via the inclusion-exclusion principle. |

### Types

| Type | What it is |
| --- | --- |
| `Distribution` | Distribution is a finitely-supported discrete probability distribution over real-valued outcomes. |
| `JointDistribution` | JointDistribution is a bivariate discrete distribution over a grid of outcomes. |
| `MarkovChain` | MarkovChain is a finite-state discrete-time Markov chain represented by a row-stochastic transition matrix P, where P`i``j` is the probability of… |

<details>
<summary><code>Distribution</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Bernoulli(p float64) (Distribution, error)` | Bernoulli returns the Bernoulli distribution with success probability p, supported on {0, 1} with P(1) = p and P(0) = 1-p. |
| `func Binomial(n int, p float64) (Distribution, error)` | Binomial returns the binomial distribution for n independent Bernoulli trials each with success probability p, supported on {0, 1, …, n}. |
| `func DiscreteUniform(a, b int) (Distribution, error)` | DiscreteUniform returns the uniform distribution over the consecutive integers a, a+1, …, b (inclusive), each with probability 1/(b-a+1). |
| `func Geometric(p float64, kmax int) (Distribution, error)` | Geometric returns the geometric distribution for the number of trials up to and including the first success (support {1, 2, …, kmax}) with… |
| `func Mixture(weights []float64, components []Distribution) (Distribution, error)` | Mixture returns the finite mixture Σ_k weights`k`·components`k`, a single Distribution whose outcomes are the union of the components' supports… |
| `func NewDistribution(outcomes, probs []float64) (Distribution, error)` | NewDistribution builds a Distribution from parallel outcome and probability slices. |
| `func Poisson(lambda float64, kmax int) (Distribution, error)` | Poisson returns the Poisson distribution with rate lambda truncated to the support {0, 1, …, kmax} and renormalized so its probabilities sum to one. |
| `func Uniform(outcomes []float64) (Distribution, error)` | Uniform returns the uniform distribution over the given distinct outcomes, assigning probability 1/n to each. |
| `func (d Distribution) Affine(a, b float64) Distribution` | Affine returns the distribution of aX + b, applying a scale followed by a shift in a single pass. |
| `func (d Distribution) CDF(x float64) float64` | CDF returns the cumulative distribution function P(X <= x): the total probability of all outcomes less than or equal to x. |
| `func (d Distribution) CGF(t float64) float64` | CGF returns the cumulant-generating function K(t) = ln E[e^{tX}] = ln M(t) evaluated at t. |
| `func (d Distribution) CentralMoment(k int) float64` | CentralMoment returns the k-th central moment E[(X - E`X`)^k]. |
| `func (d Distribution) CharacteristicFunction(t float64) complex128` | CharacteristicFunction returns φ(t) = E[e^{itX}], the characteristic function of X evaluated at real argument t, as a complex value. |
| `func (d Distribution) Convolve(other Distribution) Distribution` | Convolve returns the distribution of the sum X + Y of two independent random variables X (the receiver) and Y (the argument): every pair of outcomes… |
| `func (d Distribution) ConvolvePower(n int) (Distribution, error)` | ConvolvePower returns the distribution of the sum of n independent copies of X (the receiver): X_1 + X_2 + … + X_n. |
| `func (d Distribution) Entropy() float64` | Entropy returns the Shannon entropy H(X) = -Σ_i Probs`i` ln Probs`i` measured in nats. |
| `func (d Distribution) EntropyBits() float64` | EntropyBits returns the Shannon entropy measured in bits, i.e. |
| `func (d Distribution) Expectation(g func(float64) float64) float64` | Expectation returns E[g(X)] = Σ_i g(Outcomes`i`) · Probs`i`, the expected value of the arbitrary function g applied to the random variable. |
| `func (d Distribution) Kurtosis() float64` | Kurtosis returns the excess kurtosis, E[(X - E`X`)^4] / σ^4 - 3, which is zero for a normal distribution. |
| `func (d Distribution) Len() int` | Len returns the number of distinct outcomes in the support. |
| `func (d Distribution) MGF(t float64) float64` | MGF returns the moment-generating function M(t) = E[e^{tX}] evaluated numerically at t. |
| `func (d Distribution) MGFDerivativeMoment(k int, h float64) float64` | MGFDerivativeMoment returns the k-th raw moment estimated by numerically differentiating the moment-generating function k times at t = 0 using a… |
| `func (d Distribution) Max() float64` | Max returns the largest outcome in the support. |
| `func (d Distribution) Mean() float64` | Mean returns the expected value E`X` of the distribution. |
| `func (d Distribution) Median() float64` | Median returns the 0.5 quantile of the distribution. |
| `func (d Distribution) Min() float64` | Min returns the smallest outcome in the support. |
| `func (d Distribution) Mode() float64` | Mode returns the outcome carrying the greatest probability mass. |
| `func (d Distribution) Moment(k int) float64` | Moment returns the k-th raw moment E[X^k]. |
| `func (d Distribution) Normalize() (Distribution, error)` | Normalize returns a copy of the distribution whose probabilities are rescaled to sum to exactly one. |
| `func (d Distribution) PGF(z float64) float64` | PGF returns the probability-generating function G(z) = E[z^X] = Σ_i Probs`i` z^{Outcomes`i`} evaluated at z. |
| `func (d Distribution) PMF(x float64) float64` | PMF returns the probability mass P(X = x), i.e. |
| `func (d Distribution) Quantile(q float64) float64` | Quantile returns the smallest outcome x such that P(X <= x) >= q, i.e. |
| `func (d Distribution) Scale(a float64) Distribution` | Scale returns the distribution of aX, multiplying every outcome by a. |
| `func (d Distribution) Shift(b float64) Distribution` | Shift returns the distribution of X + b, translating every outcome by b. |
| `func (d Distribution) Skewness() float64` | Skewness returns the standardized third central moment, E[(X - E`X`)^3] / σ^3, a dimensionless measure of asymmetry. |
| `func (d Distribution) Standardize() (Distribution, error)` | Standardize returns the distribution of (X - E`X`) / σ, the standardized random variable with mean zero and unit variance. |
| `func (d Distribution) StdDev() float64` | StdDev returns the standard deviation, the non-negative square root of the variance. |
| `func (d Distribution) Support() []float64` | Support returns a copy of the outcome values in ascending order. |
| `func (d Distribution) Transform(g func(float64) float64) Distribution` | Transform returns the distribution of the random variable Y = g(X), applying g to every outcome. |
| `func (d Distribution) Validate() error` | Validate reports whether the distribution is well-formed: equal-length non-empty slices, non-negative finite probabilities, and a total mass of one… |
| `func (d Distribution) Variance() float64` | Variance returns Var(X) = E[(X - E`X`)^2], the population variance of the distribution. |

</details>

<details>
<summary><code>JointDistribution</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func IndependentJoint(a, b Distribution) JointDistribution` | IndependentJoint builds the joint distribution of two independent random variables from their marginals, setting P`i``j` = a.Probs`i`·b.Probs`j`. |
| `func NewJointDistribution(x, y []float64, p [][]float64) (JointDistribution, error)` | NewJointDistribution builds a JointDistribution and validates it: X and Y must be non-empty, P must be a len(X)-by-len(Y) grid of non-negative finite… |
| `func (j JointDistribution) ConditionalXGivenY(k int) (Distribution, error)` | ConditionalXGivenY returns the conditional distribution of X given Y = Y`k`, i.e. |
| `func (j JointDistribution) ConditionalYGivenX(i int) (Distribution, error)` | ConditionalYGivenX returns the conditional distribution of Y given X = X`i`, i.e. |
| `func (j JointDistribution) Correlation() float64` | Correlation returns the Pearson correlation coefficient Cov(X, Y) / (σ_X · σ_Y), a value in [-1, 1]. |
| `func (j JointDistribution) Covariance() float64` | Covariance returns Cov(X, Y) = E`XY` - E`X`·E`Y`, the covariance of the two variables under the joint distribution. |
| `func (j JointDistribution) ExpectationXY(g func(x, y float64) float64) float64` | ExpectationXY returns E[g(X, Y)] = Σ_{i,k} g(X`i`, Y`k`) · P`i``k` for an arbitrary function g of the two variables. |
| `func (j JointDistribution) Independent() bool` | Independent reports whether X and Y are independent, i.e. |
| `func (j JointDistribution) MarginalX() Distribution` | MarginalX returns the marginal distribution of X obtained by summing the joint probabilities across the Y axis for each row. |
| `func (j JointDistribution) MarginalY() Distribution` | MarginalY returns the marginal distribution of Y obtained by summing the joint probabilities across the X axis for each column. |

</details>

<details>
<summary><code>MarkovChain</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewMarkovChain(p [][]float64) (MarkovChain, error)` | NewMarkovChain builds a MarkovChain from a transition matrix, validating that it is square, non-negative, and row-stochastic (each row sums to one… |
| `func (m MarkovChain) AbsorbingStates() []int` | AbsorbingStates returns the indices of the absorbing states, those i with P`i``i` equal to one (and hence no probability of leaving), in ascending… |
| `func (m MarkovChain) AbsorptionProbabilities() ([][]float64, error)` | AbsorptionProbabilities returns the matrix B = N·R of absorption probabilities, where N is the fundamental matrix and R is the… |
| `func (m MarkovChain) DistributionAfter(initial []float64, n int) ([]float64, error)` | DistributionAfter returns the state distribution after n steps starting from the initial distribution, i.e. |
| `func (m MarkovChain) ExpectedStepsToAbsorption() ([]float64, error)` | ExpectedStepsToAbsorption returns, for each state, the expected number of steps until the chain is absorbed. |
| `func (m MarkovChain) FundamentalMatrix() ([][]float64, error)` | FundamentalMatrix returns the fundamental matrix N = (I - Q)^{-1} of an absorbing chain, where Q is the transient-to-transient submatrix of P. Rows… |
| `func (m MarkovChain) IsAbsorbing() bool` | IsAbsorbing reports whether the chain is an absorbing Markov chain: it has at least one absorbing state and every state can reach an absorbing state. |
| `func (m MarkovChain) IsIrreducible() bool` | IsIrreducible reports whether the chain is irreducible, i.e. |
| `func (m MarkovChain) IsRegular() bool` | IsRegular reports whether the chain is regular (primitive): some power P^k has all strictly positive entries. |
| `func (m MarkovChain) MeanRecurrenceTimes() ([]float64, error)` | MeanRecurrenceTimes returns the vector of mean recurrence times of an irreducible chain, where entry i is 1/π_i and π is the stationary… |
| `func (m MarkovChain) NStep(n int) ([][]float64, error)` | NStep returns the n-step transition matrix P^n, whose (i, j) entry is the probability of moving from state i to state j in exactly n steps. |
| `func (m MarkovChain) Reachable(i, j int) bool` | Reachable reports whether state j is reachable from state i in zero or more steps (a state is always reachable from itself). |
| `func (m MarkovChain) Size() int` | Size returns the number of states in the chain. |
| `func (m MarkovChain) StationaryDistribution() ([]float64, error)` | StationaryDistribution returns a stationary distribution π satisfying πP = π and Σ π = 1. |
| `func (m MarkovChain) Step(dist []float64) ([]float64, error)` | Step advances a distribution over states by one transition, returning the row vector dist·P. dist must have length equal to the number of states. |
| `func (m MarkovChain) TransientStates() []int` | TransientStates returns the indices of the non-absorbing states in ascending order. |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/probability).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
