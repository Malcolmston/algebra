# stats — small, dependency-free statistics and probability toolkit written entirely with the Go

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/stats.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/stats)

Package stats is a small, dependency-free statistics and probability toolkit
written entirely with the Go standard library (only math and sort). It
operates on plain []float64 data and is completely standalone: it does not
import the parent algebra package.

The package is organised into four groups.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/stats"
```

## Usage

This is the package's own `ExampleBinomial`, so it compiles and its output is
asserted on every `go test ./stats/`.

```go
b := stats.Binomial{N: 5, P: 0.5}
	fmt.Printf("PMF(2)=%.4f mean=%.1f var=%.2f\n",
		b.PMF(2), b.Mean(), b.Variance())
```

```
PMF(2)=0.3125 mean=2.5 var=1.25
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func Choose(n, k int) float64` | Choose returns the binomial coefficient C(n, k) = n! / (k!·(n-k)!), the number of ways to choose k items from n without regard to order. |
| `func Correlation(xs, ys []float64) float64` | Correlation returns the Pearson product-moment correlation coefficient between xs and ys, a value in [-1, 1]. |
| `func CorrelationMatrix(cols [][]float64) [][]float64` | CorrelationMatrix returns the symmetric k×k Pearson correlation matrix for k variables, where cols`j` holds the samples of variable j. |
| `func Covariance(xs, ys []float64) float64` | Covariance returns the unbiased sample covariance between xs and ys using the n-1 denominator. |
| `func CovarianceMatrix(cols [][]float64) [][]float64` | CovarianceMatrix returns the symmetric k×k unbiased (n-1) sample covariance matrix for k variables, where cols`j` holds the samples of variable j… |
| `func Factorial(n int) float64` | Factorial returns n! as a float64. |
| `func GeometricMean(xs []float64) float64` | GeometricMean returns the geometric mean of xs, the n-th root of the product of the values. |
| `func HarmonicMean(xs []float64) float64` | HarmonicMean returns the harmonic mean of xs, the reciprocal of the mean of the reciprocals. |
| `func IQR(xs []float64) float64` | IQR returns the interquartile range of xs, the difference between the 75th and 25th percentiles. |
| `func Kurtosis(xs []float64) float64` | Kurtosis returns the population excess kurtosis of xs (the fourth standardized moment minus 3, so a normal distribution has kurtosis 0). |
| `func LinearRegression(xs, ys []float64) (slope, intercept, r2 float64)` | LinearRegression fits an ordinary-least-squares straight line y = slope·x + intercept to the paired data (xs, ys) and returns its slope, intercept… |
| `func Max(xs []float64) float64` | Max returns the largest value in xs, or NaN if xs is empty. |
| `func Mean(xs []float64) float64` | Mean returns the arithmetic mean of xs, or NaN if xs is empty. |
| `func MeanConfidenceInterval(xs []float64, level float64) (lo, hi float64)` | MeanConfidenceInterval returns the two-sided confidence interval [lo, hi] for the population mean of xs at the given confidence level (for example… |
| `func MeanConfidenceIntervalZ(xs []float64, sigma, level float64) (lo, hi float64)` | MeanConfidenceIntervalZ returns the two-sided confidence interval [lo, hi] for the population mean of xs at the given confidence level using the… |
| `func Median(xs []float64) float64` | Median returns the median of xs, or NaN if xs is empty. |
| `func Min(xs []float64) float64` | Min returns the smallest value in xs, or NaN if xs is empty. |
| `func Mode(xs []float64) []float64` | Mode returns the value or values that occur most frequently in xs, sorted in ascending order. |
| `func MultipleLinearRegression(X [][]float64, y []float64, intercept bool) (coeffs []float64, r2 float64)` | MultipleLinearRegression fits an ordinary-least-squares linear model predicting y from the predictors in X and returns the fitted coefficients… |
| `func NormalCDF(x, mu, sigma float64) float64` | NormalCDF returns the cumulative probability P(X <= x) for a normal distribution with the given mean mu and standard deviation sigma. |
| `func NormalPDF(x, mu, sigma float64) float64` | NormalPDF returns the density at x of a normal distribution with the given mean mu and standard deviation sigma. |
| `func Percentile(xs []float64, p float64) float64` | Percentile returns the p-th percentile of xs for p in [0, 100]. |
| `func Perm(n, k int) float64` | Perm returns the number of permutations P(n, k) = n! / (n-k)!, the number of ordered arrangements of k items chosen from n. |
| `func PopStdDev(xs []float64) float64` | PopStdDev returns the population standard deviation of xs (the square root of `PopVariance`). |
| `func PopVariance(xs []float64) float64` | PopVariance returns the population variance of xs using the n denominator. |
| `func Predict(coeffs []float64, x []float64, intercept bool) float64` | Predict returns the linear prediction for a single predictor row x given a coefficient vector, honoring the same intercept convention as the fitting… |
| `func Product(xs []float64) float64` | Product returns the product of xs. |
| `func ProportionConfidenceInterval(successes, n int, level float64) (lo, hi float64)` | ProportionConfidenceInterval returns the two-sided Wilson score confidence interval [lo, hi] for a binomial proportion given successes out of n… |
| `func Quantile(xs []float64, q float64) float64` | Quantile returns the q-quantile of xs for q in [0, 1] using linear interpolation between the closest ranks (the "type 7" method used by NumPy and R… |
| `func Range(xs []float64) float64` | Range returns the difference between the largest and smallest values in xs, or NaN if xs is empty. |
| `func Resample(xs []float64, r *rand.Rand) []float64` | Resample returns a bootstrap sample of xs: a new slice of the same length whose elements are drawn from xs uniformly with replacement using r. |
| `func RidgeRegression(X [][]float64, y []float64, lambda float64, intercept bool) []float64` | RidgeRegression fits an L2-penalized (ridge) linear model, solving (XᵀX + λI)·b = Xᵀy, and returns the fitted coefficients. |
| `func SampleIntN(d DiscreteSampler, n int, r *rand.Rand) []int` | SampleIntN draws n independent integer variates from d using r and returns them in a new slice. |
| `func SampleN(d ContinuousSampler, n int, r *rand.Rand) []float64` | SampleN draws n independent variates from d using r and returns them in a new slice. |
| `func Shuffle(xs []float64, r *rand.Rand)` | Shuffle permutes xs in place using r with the Fisher–Yates algorithm, giving each of the len(xs)! orderings equal probability. |
| `func Skewness(xs []float64) float64` | Skewness returns the population coefficient of skewness of xs (the third standardized moment). |
| `func StdDev(xs []float64) float64` | StdDev returns the sample standard deviation of xs (the square root of `Variance`). |
| `func Sum(xs []float64) float64` | Sum returns the sum of xs. |
| `func Variance(xs []float64) float64` | Variance returns the unbiased sample variance of xs using the n-1 (Bessel) denominator. |
| `func WeightedMean(xs, weights []float64) float64` | WeightedMean returns the weighted arithmetic mean of xs with the given weights: sum(w*x) / sum(w). |
| `func ZScore(x, mean, stdDev float64) float64` | ZScore returns the standard score of x given a distribution mean and standard deviation: (x - mean) / stdDev. |

### Types

| Type | What it is |
| --- | --- |
| `Accumulator` | Accumulator is a constant-memory, single-pass streaming replacement for the multi-pass descriptive functions (`Mean`, `Variance`, `Skewness`,… |
| `Beta` | Beta is a beta distribution on the interval [0, 1] with shape parameters Alpha (> 0) and Beta (> 0). |
| `Binomial` | Binomial is a binomial distribution: the number of successes in N independent trials, each succeeding with probability P. N must be >= 0 and P must… |
| `BootstrapResult` | BootstrapResult holds the outcome of a bootstrap resampling analysis: Estimate is the statistic evaluated on the original sample, StdErr is the… |
| `ChiSquared` | ChiSquared is a chi-squared distribution with K degrees of freedom (K > 0). |
| `ContinuousSampler` | ContinuousSampler is implemented by continuous distributions that can draw a single random variate. |
| `CovAccumulator` | CovAccumulator is a constant-memory, single-pass streaming replacement for the multi-pass bivariate functions `Covariance`, `Correlation` and… |
| `DiscreteSampler` | DiscreteSampler is implemented by discrete distributions that can draw a single integer variate. |
| `Exponential` | Exponential is an exponential distribution with rate parameter Lambda > 0, modelling the waiting time between events in a Poisson process. |
| `FDist` | FDist is Fisher's F distribution with degrees of freedom D1 and D2 (both > 0), the ratio of two independent chi-squared variates each divided by its… |
| `Gamma` | Gamma is a gamma distribution parameterized by Shape (k, > 0) and Scale (θ, > 0). |
| `Geometric` | Geometric is a geometric distribution: the number of Bernoulli trials up to and including the first success, each trial succeeding with probability P… |
| `GeometricFailures` | GeometricFailures is a geometric distribution parameterized by the number of failures before the first success. |
| `Hypergeometric` | Hypergeometric is a hypergeometric distribution: the number of successes k in Draws draws made without replacement from a finite population of N… |
| `LogNormal` | LogNormal is a log-normal distribution: a positive variable whose natural logarithm is normally distributed with mean Mu and standard deviation Sigma… |
| `NegativeBinomial` | NegativeBinomial is a negative binomial distribution: the number of failures before the R-th success in a sequence of independent Bernoulli trials,… |
| `NegativeBinomialInt` | NegativeBinomialInt is a negative binomial distribution parameterized by an integer number of successes R (>= 1) and success probability P in (0, 1]. |
| `Normal` | Normal is a normal (Gaussian) distribution with mean Mu and standard deviation Sigma. |
| `Poisson` | Poisson is a Poisson distribution with rate parameter Lambda > 0: the number of events in a fixed interval when events occur independently at a… |
| `StudentT` | StudentT is a Student's t distribution with Nu degrees of freedom (Nu > 0). |
| `TestResult` | TestResult is the common outcome of the classical significance tests in this file. |
| `Uniform` | Uniform is a continuous uniform distribution on the closed interval [A, B] with A < B. |
| `Weibull` | Weibull is a Weibull distribution with shape parameter Shape (k, > 0) and scale parameter Scale (λ, > 0). |

<details>
<summary><code>Accumulator</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (a *Accumulator) Count() int` | Count returns the number of observations pushed so far. |
| `func (a *Accumulator) Kurtosis() float64` | Kurtosis returns the population excess kurtosis (the fourth standardized moment minus 3, so a normal distribution has kurtosis 0), computed from the… |
| `func (a *Accumulator) Max() float64` | Max returns the largest observation pushed so far, or NaN if none have been pushed. |
| `func (a *Accumulator) Mean() float64` | Mean returns the arithmetic mean of the observations, or NaN if none have been pushed. |
| `func (a *Accumulator) Merge(o Accumulator)` | Merge folds the observations summarized by o into the receiver as if they had been pushed directly, using the pairwise-moment merge formulas of Chan,… |
| `func (a *Accumulator) Min() float64` | Min returns the smallest observation pushed so far, or NaN if none have been pushed. |
| `func (a *Accumulator) PopStdDev() float64` | PopStdDev returns the population standard deviation (the square root of `Accumulator.PopVariance`), or NaN if no observations have been pushed. |
| `func (a *Accumulator) PopVariance() float64` | PopVariance returns the population variance using the n denominator, or NaN if no observations have been pushed. |
| `func (a *Accumulator) Push(x float64)` | Push incorporates a single observation x into the accumulator, updating the count, mean, central moments M2, M3 and M4, and the running minimum and… |
| `func (a *Accumulator) PushAll(xs []float64)` | PushAll incorporates every observation in xs by calling `Accumulator.Push` on each in order. |
| `func (a *Accumulator) Reset()` | Reset returns the accumulator to its empty zero-value state so it can be reused without allocating. |
| `func (a *Accumulator) Skewness() float64` | Skewness returns the population coefficient of skewness (the third standardized moment), computed from the running moments as sqrt(n)*M3/M2^1.5. |
| `func (a *Accumulator) StdDev() float64` | StdDev returns the sample standard deviation (the square root of `Accumulator.Variance`), or NaN if fewer than two observations have been pushed. |
| `func (a *Accumulator) Sum() float64` | Sum returns the sum of all observations, computed as the mean times the count. |
| `func (a *Accumulator) Variance() float64` | Variance returns the unbiased sample variance using the n-1 (Bessel) denominator, or NaN if fewer than two observations have been pushed. |

</details>

<details>
<summary><code>Beta</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (b Beta) CDF(x float64) float64` | CDF returns the cumulative probability P(X <= x), the regularized incomplete beta function I_x(Alpha, Beta). |
| `func (b Beta) Mean() float64` | Mean returns the mean of the distribution, Alpha/(Alpha+Beta). |
| `func (b Beta) PDF(x float64) float64` | PDF returns the probability density at x. |
| `func (b Beta) Quantile(p float64) float64` | Quantile returns the inverse CDF: the value x such that CDF(x) = p, for p in [0, 1]. |
| `func (b Beta) Sample(r *rand.Rand) float64` | Sample draws a single variate from the beta distribution using r, as G1/(G1+G2) where G1 ~ Gamma(Alpha, 1) and G2 ~ Gamma(Beta, 1) are independent. |
| `func (b Beta) Variance() float64` | Variance returns the variance of the distribution, Alpha·Beta / ((Alpha+Beta)²·(Alpha+Beta+1)). |

</details>

<details>
<summary><code>Binomial</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (b Binomial) CDF(k int) float64` | CDF returns the cumulative probability P(X <= k). |
| `func (b Binomial) Mean() float64` | Mean returns the mean of the distribution, N·P. |
| `func (b Binomial) PMF(k int) float64` | PMF returns the probability mass P(X = k), the probability of exactly k successes. |
| `func (b Binomial) SampleInt(r *rand.Rand) int` | SampleInt draws a single count of successes from the binomial distribution using r. |
| `func (b Binomial) Variance() float64` | Variance returns the variance of the distribution, N·P·(1-P). |

</details>

<details>
<summary><code>BootstrapResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Bootstrap(xs []float64, stat func([]float64) float64, iters int, level float64, r *rand.Rand) BootstrapResult` | Bootstrap performs a nonparametric bootstrap of the statistic stat over xs. |
| `func BootstrapMean(xs []float64, iters int, level float64, r *rand.Rand) BootstrapResult` | BootstrapMean is a convenience wrapper around Bootstrap that resamples the sample mean of xs, returning its bootstrap estimate, standard error, and… |

</details>

<details>
<summary><code>ChiSquared</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (c ChiSquared) CDF(x float64) float64` | CDF returns the cumulative probability P(X <= x), evaluated with the regularized lower incomplete gamma function. |
| `func (c ChiSquared) Mean() float64` | Mean returns the mean of the distribution, K. |
| `func (c ChiSquared) PDF(x float64) float64` | PDF returns the probability density at x. |
| `func (c ChiSquared) Sample(r *rand.Rand) float64` | Sample draws a single variate from the chi-squared distribution using r, as a Gamma(K/2, 2) variate. |
| `func (c ChiSquared) Variance() float64` | Variance returns the variance of the distribution, 2K. |

</details>

<details>
<summary><code>CovAccumulator</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (a *CovAccumulator) Correlation() float64` | Correlation returns the Pearson product-moment correlation coefficient between x and y, a value in [-1, 1]. |
| `func (a *CovAccumulator) Count() int` | Count returns the number of paired observations pushed so far. |
| `func (a *CovAccumulator) Covariance() float64` | Covariance returns the unbiased sample covariance between x and y using the n-1 denominator, or NaN if fewer than two pairs have been pushed. |
| `func (a *CovAccumulator) Intercept() float64` | Intercept returns the intercept of the ordinary-least-squares line y = slope*x + intercept fitted to the pushed pairs, equal to the intercept… |
| `func (a *CovAccumulator) Merge(o CovAccumulator)` | Merge folds the paired observations summarized by o into the receiver as if they had been pushed directly, using the pairwise co-moment merge formula. |
| `func (a *CovAccumulator) Push(x, y float64)` | Push incorporates a single paired observation (x, y), updating the running means, the second moments of x and y, and the co-moment C in one pass… |
| `func (a *CovAccumulator) PushAll(xs, ys []float64)` | PushAll incorporates every paired observation in xs and ys by calling `CovAccumulator.Push` on each pair in order. |
| `func (a *CovAccumulator) Reset()` | Reset returns the accumulator to its empty zero-value state so it can be reused without allocating. |
| `func (a *CovAccumulator) Slope() float64` | Slope returns the slope of the ordinary-least-squares line y = slope*x + intercept fitted to the pushed pairs, equal to the slope reported by… |

</details>

<details>
<summary><code>Exponential</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (e Exponential) CDF(x float64) float64` | CDF returns the cumulative probability P(X <= x). |
| `func (e Exponential) Mean() float64` | Mean returns the mean of the distribution, 1/Lambda. |
| `func (e Exponential) PDF(x float64) float64` | PDF returns the probability density at x. |
| `func (e Exponential) Quantile(p float64) float64` | Quantile returns the inverse CDF for p in [0, 1]. |
| `func (e Exponential) Sample(r *rand.Rand) float64` | Sample draws a single variate from the exponential distribution using r via inverse-CDF sampling, reusing the distribution's Quantile method. |
| `func (e Exponential) Variance() float64` | Variance returns the variance of the distribution, 1/Lambda². |

</details>

<details>
<summary><code>FDist</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (f FDist) CDF(x float64) float64` | CDF returns the cumulative probability P(X <= x), the regularized incomplete beta function I_{D1·x/(D1·x+D2)}(D1/2, D2/2). |
| `func (f FDist) Mean() float64` | Mean returns the mean of the distribution, D2/(D2-2) for D2 > 2, and NaN otherwise. |
| `func (f FDist) PDF(x float64) float64` | PDF returns the probability density at x. |
| `func (f FDist) Quantile(p float64) float64` | Quantile returns the inverse CDF: the value x such that CDF(x) = p, for p in [0, 1]. |
| `func (f FDist) Sample(r *rand.Rand) float64` | Sample draws a single variate from the F distribution using r, as the ratio (X1/D1)/(X2/D2) of two independent chi-squared variates scaled by their… |
| `func (f FDist) Variance() float64` | Variance returns the variance of the distribution, 2·D2²·(D1+D2-2) / (D1·(D2-2)²·(D2-4)) for D2 > 4. |

</details>

<details>
<summary><code>Gamma</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (g Gamma) CDF(x float64) float64` | CDF returns the cumulative probability P(X <= x), evaluated with the regularized lower incomplete gamma function. |
| `func (g Gamma) Mean() float64` | Mean returns the mean of the distribution, Shape·Scale. |
| `func (g Gamma) PDF(x float64) float64` | PDF returns the probability density at x. |
| `func (g Gamma) Sample(r *rand.Rand) float64` | Sample draws a single variate from the gamma distribution using r via the Marsaglia–Tsang squeeze. |
| `func (g Gamma) Variance() float64` | Variance returns the variance of the distribution, Shape·Scale². |

</details>

<details>
<summary><code>Geometric</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (g Geometric) CDF(k int) float64` | CDF returns the cumulative probability P(X <= k). |
| `func (g Geometric) Mean() float64` | Mean returns the mean of the distribution, 1/P. |
| `func (g Geometric) PMF(k int) float64` | PMF returns the probability mass P(X = k), the probability that the first success occurs on trial k. |
| `func (g Geometric) Quantile(p float64) float64` | Quantile returns the smallest integer k such that CDF(k) >= p, for p in [0, 1], as a float64. |
| `func (g Geometric) SampleInt(r *rand.Rand) int` | SampleInt draws a single variate from the geometric distribution using r via inverse-transform sampling, floor(ln(1-U)/ln(1-P)). |
| `func (g Geometric) Variance() float64` | Variance returns the variance of the distribution, (1-P)/P². |

</details>

<details>
<summary><code>GeometricFailures</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (g GeometricFailures) CDF(k int) float64` | CDF returns the cumulative probability P(X <= k) = 1 - (1-P)^(k+1). |
| `func (g GeometricFailures) Mean() float64` | Mean returns the mean of the distribution, (1-P)/P. It returns NaN if P is outside (0, 1]. |
| `func (g GeometricFailures) PMF(k int) float64` | PMF returns the probability mass P(X = k), the probability of exactly k failures before the first success: (1-P)^k · P for k >= 0, and 0 for k < 0. |
| `func (g GeometricFailures) Quantile(p float64) int` | Quantile returns the smallest integer k such that CDF(k) >= p, for p in [0, 1]. |
| `func (g GeometricFailures) Variance() float64` | Variance returns the variance of the distribution, (1-P)/P². |

</details>

<details>
<summary><code>Hypergeometric</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (h Hypergeometric) CDF(k int) float64` | CDF returns the cumulative probability P(X <= k) by accumulating PMF terms from the lower support bound up to k. |
| `func (h Hypergeometric) Mean() float64` | Mean returns the mean of the distribution, Draws·K/N. It returns NaN for invalid parameters. |
| `func (h Hypergeometric) PMF(k int) float64` | PMF returns the probability mass P(X = k) = C(K, k)·C(N-K, Draws-k)/C(N, Draws) using the exact integer binomial coefficient Choose. |
| `func (h Hypergeometric) Quantile(p float64) int` | Quantile returns the smallest integer k such that CDF(k) >= p, for p in [0, 1], found by scanning k upward from the lower support bound while… |
| `func (h Hypergeometric) SampleInt(r *rand.Rand) int` | SampleInt draws a single count of successes from the hypergeometric distribution using r by simulating an urn drawn without replacement: Draws items… |
| `func (h Hypergeometric) Variance() float64` | Variance returns the variance of the distribution, Draws·(K/N)·(1-K/N)·(N-Draws)/(N-1). |

</details>

<details>
<summary><code>LogNormal</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (l LogNormal) CDF(x float64) float64` | CDF returns the cumulative probability P(X <= x), the underlying normal CDF evaluated at ln x. |
| `func (l LogNormal) Mean() float64` | Mean returns the mean of the distribution, exp(Mu + Sigma²/2). |
| `func (l LogNormal) PDF(x float64) float64` | PDF returns the probability density at x, computed as the underlying normal density at ln x divided by x. |
| `func (l LogNormal) Quantile(p float64) float64` | Quantile returns the inverse CDF for p in [0, 1], exp(Mu + Sigma·Φ⁻¹(p)) where Φ⁻¹ is the standard normal quantile. |
| `func (l LogNormal) Sample(r *rand.Rand) float64` | Sample draws a single variate from the log-normal distribution using r via inverse-CDF sampling, reusing the distribution's Quantile method. |
| `func (l LogNormal) Variance() float64` | Variance returns the variance of the distribution, (exp(Sigma²)-1)·exp(2·Mu + Sigma²). |

</details>

<details>
<summary><code>NegativeBinomial</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (nb NegativeBinomial) CDF(k int) float64` | CDF returns the cumulative probability P(X <= k), evaluated with the regularized incomplete beta function. |
| `func (nb NegativeBinomial) Mean() float64` | Mean returns the mean of the distribution, R·(1-P)/P. |
| `func (nb NegativeBinomial) PMF(k int) float64` | PMF returns the probability mass P(X = k), the probability of exactly k failures before the R-th success. |
| `func (nb NegativeBinomial) SampleInt(r *rand.Rand) int` | SampleInt draws a single count of failures from the negative binomial distribution using r via the Poisson–gamma mixture: a Gamma(R, (1-P)/P) rate… |
| `func (nb NegativeBinomial) Variance() float64` | Variance returns the variance of the distribution, R·(1-P)/P². |

</details>

<details>
<summary><code>NegativeBinomialInt</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (nb NegativeBinomialInt) CDF(k int) float64` | CDF returns the cumulative probability P(X <= k). |
| `func (nb NegativeBinomialInt) Mean() float64` | Mean returns the mean of the distribution, R·(1-P)/P. It returns NaN for invalid parameters. |
| `func (nb NegativeBinomialInt) PMF(k int) float64` | PMF returns the probability mass P(X = k), the probability of exactly k failures before the R-th success. |
| `func (nb NegativeBinomialInt) Quantile(p float64) int` | Quantile returns the smallest integer k such that CDF(k) >= p, for p in [0, 1], found by scanning k upward from the lower support bound 0 while… |
| `func (nb NegativeBinomialInt) Variance() float64` | Variance returns the variance of the distribution, R·(1-P)/P². |

</details>

<details>
<summary><code>Normal</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (n Normal) CDF(x float64) float64` | CDF returns the cumulative probability P(X <= x), evaluated with math.Erf. |
| `func (n Normal) Mean() float64` | Mean returns the mean of the distribution, Mu. |
| `func (n Normal) PDF(x float64) float64` | PDF returns the probability density at x. |
| `func (n Normal) Quantile(p float64) float64` | Quantile returns the inverse CDF: the value x such that CDF(x) = p, for p in [0, 1]. |
| `func (n Normal) Sample(r *rand.Rand) float64` | Sample draws a single variate from the normal distribution using r. |
| `func (n Normal) Variance() float64` | Variance returns the variance of the distribution, Sigma². |

</details>

<details>
<summary><code>Poisson</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (p Poisson) CDF(k int) float64` | CDF returns the cumulative probability P(X <= k). |
| `func (p Poisson) Mean() float64` | Mean returns the mean of the distribution, Lambda. |
| `func (p Poisson) PMF(k int) float64` | PMF returns the probability mass P(X = k). |
| `func (p Poisson) SampleInt(r *rand.Rand) int` | SampleInt draws a single count from the Poisson distribution using r. |
| `func (p Poisson) Variance() float64` | Variance returns the variance of the distribution, Lambda. |

</details>

<details>
<summary><code>StudentT</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (t StudentT) CDF(x float64) float64` | CDF returns the cumulative probability P(X <= x), evaluated with the regularized incomplete beta function. |
| `func (t StudentT) Mean() float64` | Mean returns the mean of the distribution, 0 for Nu > 1 and NaN otherwise. |
| `func (t StudentT) PDF(x float64) float64` | PDF returns the probability density at x. |
| `func (t StudentT) Quantile(p float64) float64` | Quantile returns the p-quantile (inverse CDF) of the Student's t distribution: the value x such that CDF(x) = p, for p in (0, 1). |
| `func (t StudentT) Sample(r *rand.Rand) float64` | Sample draws a single variate from Student's t distribution using r, as Z / sqrt(V/Nu) where Z is standard normal and V ~ ChiSquared(Nu) is… |
| `func (t StudentT) Variance() float64` | Variance returns the variance of the distribution: Nu/(Nu-2) for Nu > 2, +Inf for 1 < Nu <= 2, and NaN otherwise. |

</details>

<details>
<summary><code>TestResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ChiSquareGoodnessOfFit(observed, expected []float64) TestResult` | ChiSquareGoodnessOfFit performs Pearson's chi-squared goodness-of-fit test comparing observed counts against expected counts. |
| `func ChiSquareIndependence(table [][]float64) TestResult` | ChiSquareIndependence performs Pearson's chi-squared test of independence on a contingency table given as a slice of equal-length rows of counts. |
| `func MannWhitneyU(xs, ys []float64) TestResult` | MannWhitneyU performs the two-sided Mann–Whitney U test (equivalently the Wilcoxon rank-sum test) of the null hypothesis that xs and ys come from… |
| `func OneSampleTTest(xs []float64, mu0 float64) TestResult` | OneSampleTTest performs a two-sided one-sample Student's t-test of the null hypothesis that the population mean of xs equals mu0. |
| `func OneWayANOVA(groups ...[]float64) TestResult` | OneWayANOVA performs a one-way analysis of variance testing the null hypothesis that all supplied groups share a common mean. |
| `func PairedTTest(xs, ys []float64) TestResult` | PairedTTest performs a two-sided paired-sample t-test of the null hypothesis that the mean of the paired differences xs`i`-ys`i` is zero. |
| `func TwoSampleTTest(xs, ys []float64, equalVar bool) TestResult` | TwoSampleTTest performs a two-sided two-sample t-test of the null hypothesis that xs and ys share a common mean. |

</details>

<details>
<summary><code>Uniform</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (u Uniform) CDF(x float64) float64` | CDF returns the cumulative probability P(X <= x). |
| `func (u Uniform) Mean() float64` | Mean returns the mean of the distribution, (A+B)/2. |
| `func (u Uniform) PDF(x float64) float64` | PDF returns the probability density at x: 1/(B-A) inside [A, B] and 0 outside. |
| `func (u Uniform) Quantile(p float64) float64` | Quantile returns the inverse CDF for p in [0, 1]. |
| `func (u Uniform) Sample(r *rand.Rand) float64` | Sample draws a single variate from the uniform distribution on [A, B) using a single call to r.Float64: A + (B-A)·U. |
| `func (u Uniform) Variance() float64` | Variance returns the variance of the distribution, (B-A)²/12. |

</details>

<details>
<summary><code>Weibull</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (w Weibull) CDF(x float64) float64` | CDF returns the cumulative probability P(X <= x), 1 - exp(-(x/λ)^k). |
| `func (w Weibull) Mean() float64` | Mean returns the mean of the distribution, λ·Γ(1 + 1/k). |
| `func (w Weibull) PDF(x float64) float64` | PDF returns the probability density at x, (k/λ)·(x/λ)^{k-1}·exp(-(x/λ)^k). |
| `func (w Weibull) Quantile(p float64) float64` | Quantile returns the inverse CDF for p in [0, 1], λ·(-ln(1-p))^{1/k}. |
| `func (w Weibull) Sample(r *rand.Rand) float64` | Sample draws a single variate from the Weibull distribution using r via inverse-CDF sampling, reusing the distribution's Quantile method. |
| `func (w Weibull) Variance() float64` | Variance returns the variance of the distribution, λ²·(Γ(1 + 2/k) - Γ(1 + 1/k)²). |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/stats).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
