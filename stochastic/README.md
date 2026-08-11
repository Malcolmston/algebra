# stochastic — implements stochastic processes and stochastic differential equation (SDE) simulation

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/stochastic.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/stochastic)

Package stochastic implements stochastic processes and stochastic differential
equation (SDE) simulation in pure Go.

- Random variate generation: uniform, normal, exponential, gamma, beta,
Poisson, binomial, geometric and many other one-dimensional distributions,
plus categorical sampling. All draws come from a caller-supplied
*math/rand.Rand, so results are fully deterministic given a seed.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/stochastic"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AbsorbingRandomWalk(rng *rand.Rand, start, lo, hi int, p float64) (position, steps int)` | AbsorbingRandomWalk simulates a simple random walk that starts at start, steps +1 with probability p, and stops when it first reaches either boundary… |
| `func ArrivalsToInterArrivals(arrivals []float64) []float64` | ArrivalsToInterArrivals converts a sorted slice of arrival times into the gaps between successive arrivals (the first gap is measured from time 0). |
| `func Autocorrelation(x []float64, lag int) float64` | Autocorrelation returns the sample autocorrelation of x at the given lag, normalized by the lag-0 autocovariance. |
| `func Autocovariance(x []float64, lag int) float64` | Autocovariance returns the sample autocovariance of x at the given lag. |
| `func BernoulliSample(rng *rand.Rand, p float64) bool` | BernoulliSample returns true with probability p and false otherwise. |
| `func BetaSample(rng *rand.Rand, alpha, beta float64) float64` | BetaSample returns a sample from the beta distribution with shape parameters alpha and beta on the interval (0, 1). |
| `func BinomialSample(rng *rand.Rand, n int, p float64) int` | BinomialSample returns a sample from the binomial distribution Binomial(n, p), the number of successes in n independent Bernoulli(p) trials. |
| `func BoxMullerPair(rng *rand.Rand) (float64, float64)` | BoxMullerPair returns two independent standard normal variates generated from a single pair of uniforms using the Box-Muller transform. |
| `func BrownianBridgeCovariance(s, t, T float64) float64` | BrownianBridgeCovariance returns the covariance s*(T-t)/T (for s <= t) of a standard Brownian bridge on [0, T]. |
| `func BrownianCovariance(s, t float64) float64` | BrownianCovariance returns the covariance min(s, t) of standard Brownian motion at times s and t. |
| `func BrownianHittingProbabilityBeforeT(a, t, sigma float64) float64` | BrownianHittingProbabilityBeforeT returns the probability that a driftless Brownian motion with volatility sigma reaches level a (a != 0) at some… |
| `func BrownianIncrements(rng *rand.Rand, dt float64, n int) []float64` | BrownianIncrements returns n independent Brownian increments over time steps of size dt, each distributed as N(0, dt). |
| `func BrownianMaxCDF(a, t, sigma float64) float64` | BrownianMaxCDF returns the CDF P(max_{0<=s<=t} W(s) <= a) of the running maximum of a driftless Brownian motion with volatility sigma. |
| `func CategoricalSample(rng *rand.Rand, probs []float64) int` | CategoricalSample is a convenience wrapper around DiscreteSample for probability vectors that already sum to one. |
| `func CauchySample(rng *rand.Rand, location, scale float64) float64` | CauchySample returns a sample from the Cauchy distribution with the given location and scale. |
| `func ChiSquaredSample(rng *rand.Rand, k float64) float64` | ChiSquaredSample returns a sample from the chi-squared distribution with k degrees of freedom. |
| `func CompoundPoisson(rng *rand.Rand, rate, T float64, jump func(*rand.Rand) float64) float64` | CompoundPoisson returns the value X(T) = sum of jumps of a compound Poisson process with the given event rate on [0, T], where each jump is drawn… |
| `func CompoundPoissonMean(rate, T, jumpMean float64) float64` | CompoundPoissonMean returns the analytic mean rate*T*E`jump` of a compound Poisson process observed on [0, T], given the mean jump size. |
| `func CompoundPoissonVariance(rate, T, jumpMean, jumpVar float64) float64` | CompoundPoissonVariance returns the analytic variance rate*T*E[jump^2] of a compound Poisson process on [0, T], given the mean and variance of the… |
| `func ConfidenceInterval95(x []float64) (lo, hi float64)` | ConfidenceInterval95 returns the approximate 95% confidence interval for the mean of x using the normal quantile 1.959964. |
| `func DiscreteSample(rng *rand.Rand, weights []float64) int` | DiscreteSample returns an index in [0, len(weights)) chosen with probability proportional to the corresponding weight. |
| `func EnsembleFinalValues(paths []Path) []float64` | EnsembleFinalValues returns the terminal value of each path. |
| `func EnsembleMean(paths []Path) []float64` | EnsembleMean returns the pointwise mean across a set of paths sampled on the same grid. |
| `func EnsembleVariance(paths []Path) []float64` | EnsembleVariance returns the pointwise unbiased variance across a set of paths sampled on the same grid. |
| `func EstimateDrift(p Path) float64` | EstimateDrift estimates the constant drift of a path with equally spaced times as the mean increment divided by the time step. |
| `func EstimateGBMParams(p Path) (mu, sigma float64)` | EstimateGBMParams estimates the drift mu and volatility sigma of geometric Brownian motion from a positive-valued sample path with equally spaced… |
| `func EstimateOUParams(p Path) (theta, mu, sigma float64)` | EstimateOUParams estimates the mean-reversion speed theta, long-run mean mu and volatility sigma of an Ornstein-Uhlenbeck process from an equally… |
| `func EstimateVolatility(p Path) float64` | EstimateVolatility estimates the diffusion volatility of a path with equally spaced times as the standard deviation of increments divided by sqrt(dt). |
| `func EulerMaruyamaStep(sde SDE, t, x, dt, dW float64) float64` | EulerMaruyamaStep advances one Euler-Maruyama step from state x at time t over a step of size dt using the supplied Brownian increment dW. |
| `func ExponentialSample(rng *rand.Rand, rate float64) float64` | ExponentialSample returns a sample from the exponential distribution with the given rate (inverse mean). |
| `func ExponentialVector(rng *rand.Rand, n int, rate float64) []float64` | ExponentialVector returns a slice of n independent exponential variates with the given rate. |
| `func FirstPassageDriftCDF(a, drift, sigma, t float64) float64` | FirstPassageDriftCDF returns P(T_a <= t), the probability that a Brownian motion with drift and volatility sigma reaches level a > 0 by time t, using… |
| `func FirstPassageDriftMean(a, drift float64) float64` | FirstPassageDriftMean returns the mean first-passage time a/drift of a Brownian motion X(t) = drift*t + sigma*W(t) to a level a > 0, valid for drift… |
| `func FirstPassageTime(p Path, level float64) (float64, bool)` | FirstPassageTime returns the first time index and time at which a path reaches or crosses the given level, relative to its starting value. |
| `func GBMMean(s0, mu, t float64) float64` | GBMMean returns the analytic mean S0*exp(mu*t) of geometric Brownian motion at time t. |
| `func GBMVariance(s0, mu, sigma, t float64) float64` | GBMVariance returns the analytic variance of geometric Brownian motion at time t: S0^2*exp(2 mu t)*(exp(sigma^2 t) - 1). |
| `func GamblersRuinExpectedDuration(i, N int, p float64) float64` | GamblersRuinExpectedDuration returns the expected number of steps until a simple random walk starting at i, absorbing at 0 and N and stepping up with… |
| `func GamblersRuinProbability(i, N int, p float64) float64` | GamblersRuinProbability returns the probability that a simple random walk starting at i, absorbing at 0 and N, and stepping up with probability p,… |
| `func GamblersRuinWinProbability(i, N int, p float64) float64` | GamblersRuinWinProbability returns the probability that the walk of GamblersRuinProbability reaches N before 0. |
| `func GammaSample(rng *rand.Rand, shape, scale float64) float64` | GammaSample returns a sample from the gamma distribution with the given shape (k) and scale (theta) parameters, using the Marsaglia-Tsang method. |
| `func GeometricSample(rng *rand.Rand, p float64) int` | GeometricSample returns the number of Bernoulli(p) trials up to and including the first success. |
| `func InhomogeneousPoissonArrivals(rng *rand.Rand, lambda func(float64) float64, maxRate, T float64) []float64` | InhomogeneousPoissonArrivals returns the arrival times of an inhomogeneous Poisson process on [0, T] whose intensity is the function lambda(t), using… |
| `func InterArrivalTimes(rng *rand.Rand, rate float64, n int) []float64` | InterArrivalTimes returns n independent exponential inter-arrival times of a Poisson process with the given rate. |
| `func InverseGaussianCDF(x, mu, lambda float64) float64` | InverseGaussianCDF returns the cumulative distribution function of the inverse Gaussian distribution with mean mu and shape lambda at x > 0. |
| `func InverseGaussianMean(mu, lambda float64) float64` | InverseGaussianMean returns the mean mu of the inverse Gaussian distribution. |
| `func InverseGaussianPDF(x, mu, lambda float64) float64` | InverseGaussianPDF returns the density of the inverse Gaussian distribution with mean mu and shape lambda at x > 0. |
| `func InverseGaussianVariance(mu, lambda float64) float64` | InverseGaussianVariance returns the variance mu^3/lambda of the inverse Gaussian distribution. |
| `func LaplaceSample(rng *rand.Rand, location, scale float64) float64` | LaplaceSample returns a sample from the Laplace (double-exponential) distribution with the given location and scale (b). |
| `func LastPassageTime(p Path, level float64) (float64, bool)` | LastPassageTime returns the last time at which the path value equals or crosses level (detected as a sign change of value-level between consecutive… |
| `func LazyRandomWalk(rng *rand.Rand, n int, p, stay float64) []int` | LazyRandomWalk returns an n-step lazy random walk that stays put with probability stay, steps +1 with probability (1-stay)*p and -1 otherwise. |
| `func LogNormalSample(rng *rand.Rand, mu, sigma float64) float64` | LogNormalSample returns a sample from the log-normal distribution whose natural logarithm is N(mu, sigma^2). |
| `func LogReturns(p Path) []float64` | LogReturns returns the log-returns log(v`i`/v`i-1`) of the path values. |
| `func MassAction(k float64, reactants []int) func(state []int) float64` | MassAction returns a propensity function for a mass-action reaction with rate constant k and the given reactant multiplicities (reactants`i` copies… |
| `func MonteCarloExpectation(seed int64, sde SDE, x0, T float64, nSteps, paths int, payoff func(Path) float64) (mean, stderr float64)` | MonteCarloExpectation estimates E[payoff(X)] where X is the terminal-value Path of the SDE, by averaging over the given number of independent paths. |
| `func NewRNG(seed int64) *rand.Rand` | NewRNG returns a new deterministic pseudo-random source seeded with the given value. |
| `func NormalCDF(x, mu, sigma float64) float64` | NormalCDF returns the cumulative distribution function of N(mu, sigma^2) at x. |
| `func NormalPDF(x, mu, sigma float64) float64` | NormalPDF returns the probability density function of N(mu, sigma^2) at x. |
| `func NormalSample(rng *rand.Rand, mu, sigma float64) float64` | NormalSample returns a sample from the normal distribution N(mu, sigma^2) with mean mu and standard deviation sigma. |
| `func NormalVector(rng *rand.Rand, n int) []float64` | NormalVector returns a slice of n independent standard normal variates. |
| `func OUMean(x0, theta, mu, t float64) float64` | OUMean returns the analytic mean mu + (x0 - mu)*exp(-theta*t) of the Ornstein-Uhlenbeck process at time t. |
| `func OUStationaryVariance(theta, sigma float64) float64` | OUStationaryVariance returns the stationary variance sigma^2/(2 theta) of the Ornstein-Uhlenbeck process. |
| `func OUVariance(theta, sigma, t float64) float64` | OUVariance returns the analytic variance sigma^2*(1 - exp(-2 theta t))/(2 theta) of the Ornstein-Uhlenbeck process at time t. |
| `func OccupationTime(p Path, lo, hi float64) float64` | OccupationTime returns the total time the linearly interpolated path spends in the closed interval [lo, hi], computed by summing the durations of… |
| `func ParetoSample(rng *rand.Rand, xm, alpha float64) float64` | ParetoSample returns a sample from the Pareto distribution with minimum value xm and tail index alpha. |
| `func Permutation(rng *rand.Rand, n int) []int` | Permutation returns a uniformly random permutation of the integers [0, n) using the Fisher-Yates shuffle. |
| `func PoissonArrivals(rng *rand.Rand, rate, T float64) []float64` | PoissonArrivals returns the ordered arrival times of a homogeneous Poisson process of the given rate on the interval [0, T]. |
| `func PoissonArrivalsN(rng *rand.Rand, rate float64, n int) []float64` | PoissonArrivalsN returns the first n arrival times of a homogeneous Poisson process of the given rate, starting from time 0. |
| `func PoissonCDF(k int, lambda float64) float64` | PoissonCDF returns the cumulative probability P(N <= k) for a Poisson distribution with mean lambda. |
| `func PoissonPMF(k int, lambda float64) float64` | PoissonPMF returns the probability mass P(N = k) for a Poisson distribution with mean lambda. |
| `func PoissonProcessCount(rng *rand.Rand, rate, T float64) int` | PoissonProcessCount returns the number of events of a homogeneous Poisson process of the given rate observed on the interval [0, T]. |
| `func PoissonSample(rng *rand.Rand, lambda float64) int` | PoissonSample returns a sample from the Poisson distribution with mean lambda. |
| `func RademacherSample(rng *rand.Rand) int` | RademacherSample returns +1 or -1, each with probability 1/2. |
| `func RandomWalk2D(rng *rand.Rand, n int) (xs, ys []int)` | RandomWalk2D returns the x and y coordinates of an n-step random walk on the integer plane, where each step moves one unit in one of the four axis… |
| `func RandomWalkMean(n int, p float64) float64` | RandomWalkMean returns the analytic mean n*(2p-1) of the endpoint of an n-step simple random walk with up-probability p. |
| `func RandomWalkVariance(n int, p float64) float64` | RandomWalkVariance returns the analytic variance 4*n*p*(1-p) of the endpoint of an n-step simple random walk with up-probability p. |
| `func RayleighSample(rng *rand.Rand, sigma float64) float64` | RayleighSample returns a sample from the Rayleigh distribution with the given scale parameter sigma. |
| `func RealizedVariance(p Path) float64` | RealizedVariance returns the sum of squared log-returns of a positive-valued path, an estimator of the integrated variance. |
| `func RealizedVolatility(p Path) float64` | RealizedVolatility returns the square root of the realized variance, annualized by dividing the total horizon out. |
| `func ReflectingRandomWalk(rng *rand.Rand, n int, p float64, lo, hi int) []int` | ReflectingRandomWalk returns an n-step simple random walk confined to the closed interval [lo, hi] by reflecting boundaries: a step that would leave… |
| `func ReflectionPrincipleMaxProb(a, t, sigma float64) float64` | ReflectionPrincipleMaxProb returns P(max_{0<=s<=t} W(s) >= a) for a Brownian motion with volatility sigma and no drift, using the reflection… |
| `func SampleCovariance(x, y []float64) float64` | SampleCovariance returns the unbiased sample covariance of x and y, which must have the same length. |
| `func SampleMean(x []float64) float64` | SampleMean returns the arithmetic mean of x, or NaN for an empty slice. |
| `func SampleStdDev(x []float64) float64` | SampleStdDev returns the unbiased sample standard deviation of x. |
| `func SampleVariance(x []float64) float64` | SampleVariance returns the unbiased (Bessel-corrected) sample variance of x. |
| `func ShuffleFloat(rng *rand.Rand, x []float64)` | ShuffleFloat randomly permutes the elements of x in place using the Fisher-Yates shuffle. |
| `func SimpleRandomWalk(rng *rand.Rand, n int, p float64) []int` | SimpleRandomWalk returns the integer positions of an n-step simple random walk that starts at 0 and takes a +1 step with probability p and a -1 step… |
| `func StandardError(x []float64) float64` | StandardError returns the standard error of the mean of x. |
| `func StandardNormalSample(rng *rand.Rand) float64` | StandardNormalSample returns a sample from the standard normal distribution N(0, 1). |
| `func StoppingTime(p Path, predicate func(t, x float64) bool) (index int, time, value float64, ok bool)` | StoppingTime returns the first index, time and value at which the predicate holds along the path. |
| `func StudentTSample(rng *rand.Rand, nu float64) float64` | StudentTSample returns a sample from Student's t distribution with nu degrees of freedom. |
| `func SuperposedPoisson(a, b []float64) []float64` | SuperposedPoisson merges two sorted arrival sequences into a single sorted sequence. |
| `func SymmetricRandomWalk(rng *rand.Rand, n int) []int` | SymmetricRandomWalk returns a simple symmetric random walk (p = 1/2). |
| `func ThinnedPoisson(rng *rand.Rand, arrivals []float64, p float64) []float64` | ThinnedPoisson returns the sub-sequence of arrivals kept when each arrival is retained independently with probability p (Poisson thinning). |
| `func UniformIntSample(rng *rand.Rand, a, b int) int` | UniformIntSample returns a uniformly distributed integer in the closed interval [a, b]. |
| `func UniformSample(rng *rand.Rand, a, b float64) float64` | UniformSample returns a sample from the continuous uniform distribution on the half-open interval [a, b). |
| `func UniformVector(rng *rand.Rand, n int, a, b float64) []float64` | UniformVector returns a slice of n independent uniform variates on [a, b). |
| `func WaldExpectedStopping(expectedDisplacement, meanStep float64) float64` | WaldExpectedStopping returns E`T` = E`S_T`/E`step` from Wald's identity, given the expected total displacement at stopping and the mean single-step… |
| `func WeibullSample(rng *rand.Rand, k, lambda float64) float64` | WeibullSample returns a sample from the Weibull distribution with shape k and scale lambda. |

### Types

| Type | What it is |
| --- | --- |
| `DiffusionFunc` | DiffusionFunc is the diffusion coefficient b(t, x) of a scalar Ito SDE. |
| `DriftFunc` | DriftFunc is the drift coefficient a(t, x) of a scalar Ito SDE dX = a(t, X) dt + b(t, X) dW. |
| `GillespieResult` | GillespieResult holds the output of a stochastic simulation: the event times and the species-count vector after each event (States[0] is the initial… |
| `Path` | Path represents a sample path of a continuous-time process, stored as equal length slices of Times and Values sampled on a (usually regular) grid. |
| `Reaction` | Reaction describes a single reaction channel of a stochastic reaction network: a stoichiometry vector giving the change applied to the species count… |
| `ReactionNetwork` | ReactionNetwork is a collection of reaction channels acting on NumSpecies molecular species. |
| `SDE` | SDE bundles the drift and diffusion coefficients of a scalar Ito stochastic differential equation dX = Drift(t, X) dt + Diffusion(t, X) dW. |

<details>
<summary><code>GillespieResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Gillespie(rng *rand.Rand, net ReactionNetwork, initial []int, tMax float64) GillespieResult` | Gillespie runs the exact stochastic simulation algorithm (Gillespie's direct method) for the network starting from the given initial state until time… |
| `func GillespieMaxSteps(rng *rand.Rand, net ReactionNetwork, initial []int, tMax float64, maxSteps int) GillespieResult` | GillespieMaxSteps is like Gillespie but also stops after maxSteps reaction events, guarding against unbounded simulations. |
| `func TauLeaping(rng *rand.Rand, net ReactionNetwork, initial []int, tMax, tau float64) GillespieResult` | TauLeaping runs an explicit tau-leaping approximation of the network, firing each reaction a Poisson(propensity*tau) number of times per leap of size… |
| `func (r GillespieResult) FinalState() []int` | FinalState returns the last recorded species-count vector of the result. |
| `func (r GillespieResult) FinalTime() float64` | FinalTime returns the last recorded time of the result. |
| `func (r GillespieResult) NumEvents() int` | NumEvents returns the number of reaction events recorded (excluding the initial state). |
| `func (r GillespieResult) SpeciesTrajectory(i int) Path` | SpeciesTrajectory returns the times and counts of species index i across the simulation as a Path (counts cast to float64). |

</details>

<details>
<summary><code>Path</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BrownianBridge(rng *rand.Rand, a, b, T float64, nSteps int) Path` | BrownianBridge returns a Brownian bridge on [0, T] pinned to a at time 0 and b at time T, sampled on nSteps+1 grid points. |
| `func BrownianMotion(rng *rand.Rand, T float64, nSteps int) Path` | BrownianMotion returns a standard Brownian motion (Wiener process) on [0, T] sampled on nSteps+1 equally spaced grid points, starting at 0. |
| `func BrownianMotionScaled(rng *rand.Rand, sigma, T float64, nSteps int) Path` | BrownianMotionScaled returns a Brownian motion with volatility sigma (that is, variance sigma^2*t) on [0, T], starting at 0. |
| `func BrownianMotionWithDrift(rng *rand.Rand, mu, sigma, T float64, nSteps int) Path` | BrownianMotionWithDrift returns a Brownian motion with drift mu and volatility sigma, X(t) = mu*t + sigma*W(t), on [0, T] sampled on nSteps+1 grid… |
| `func CIRProcess(rng *rand.Rand, x0, a, b, sigma, T float64, nSteps int) Path` | CIRProcess simulates the Cox-Ingersoll-Ross process dX = a*(b - X) dt + sigma*sqrt(X) dW on [0, T] starting at x0 using a full-truncation Euler… |
| `func CompoundPoissonPath(rng *rand.Rand, rate, T float64, nSteps int, jump func(*rand.Rand) float64) Path` | CompoundPoissonPath returns the piecewise-constant sample path of a compound Poisson process on [0, T] sampled on nSteps+1 grid points. |
| `func ContinuousTimeRandomWalk(rng *rand.Rand, rate, T float64, jump func(*rand.Rand) float64) Path` | ContinuousTimeRandomWalk returns the sample path of a continuous-time random walk that waits an Exponential(rate) time between jumps and adds a jump… |
| `func CorrelatedBrownianMotion(rng *rand.Rand, corr [][]float64, T float64, nSteps int) ([]Path, error)` | CorrelatedBrownianMotion returns len(corr) Brownian motions on [0, T] with the given instantaneous correlation matrix, sampled on nSteps+1 grid… |
| `func CountingPath(arrivals []float64, T float64, nSteps int) Path` | CountingPath returns the counting process N(t) of the given sorted arrival times sampled on a regular grid of nSteps+1 points over [0, T]. |
| `func EulerMaruyama(rng *rand.Rand, sde SDE, x0, T float64, nSteps int) Path` | EulerMaruyama integrates the SDE on [0, T] starting at x0 with the Euler-Maruyama scheme over nSteps steps and returns the sample path. |
| `func EulerMaruyamaSystem(rng *rand.Rand, drift, diffusion []func(t float64, x []float64) float64, x0 []float64, T float64, nSteps int) []Path` | EulerMaruyamaSystem integrates a d-dimensional Ito SDE system with diagonal noise, dX_k = drift_k(t, X) dt + diffusion_k(t, X) dW_k, using… |
| `func FractionalBrownianMotion(rng *rand.Rand, hurst, T float64, nSteps int) (Path, error)` | FractionalBrownianMotion returns a fractional Brownian motion with Hurst exponent H (0 < H < 1) on [0, T] sampled on nSteps+1 grid points, using the… |
| `func GaussianRandomWalk(rng *rand.Rand, n int, x0, mu, sigma float64) Path` | GaussianRandomWalk returns an n-step random walk with independent normal increments of mean mu and standard deviation sigma, as a Path with unit time… |
| `func GeometricBrownianMotion(rng *rand.Rand, s0, mu, sigma, T float64, nSteps int) Path` | GeometricBrownianMotion returns an exact simulation of geometric Brownian motion dS = mu*S dt + sigma*S dW on [0, T] starting at S0, sampled on… |
| `func LevyFlight(rng *rand.Rand, n int, alpha float64) Path` | LevyFlight returns an n-step random walk in one dimension whose step sizes are symmetric with heavy Pareto-distributed magnitudes of tail index alpha… |
| `func Milstein(rng *rand.Rand, sde SDE, dDiffusion DiffusionFunc, x0, T float64, nSteps int) Path` | Milstein integrates the SDE with the Milstein scheme, which adds the first-order Ito correction 0.5*b*b'*(dW^2 - dt). |
| `func MilsteinNumeric(rng *rand.Rand, sde SDE, x0, T float64, nSteps int, h float64) Path` | MilsteinNumeric is like Milstein but approximates the derivative of the diffusion coefficient by a central finite difference with step h. |
| `func MultiBrownianMotion(rng *rand.Rand, dim int, T float64, nSteps int) []Path` | MultiBrownianMotion returns dim independent standard Brownian motions on [0, T] sampled on the same grid of nSteps+1 points. |
| `func NewPath(times, values []float64) Path` | NewPath constructs a Path from the given times and values. |
| `func OrnsteinUhlenbeck(rng *rand.Rand, x0, theta, mu, sigma, T float64, nSteps int) Path` | OrnsteinUhlenbeck returns an exact simulation of the Ornstein-Uhlenbeck process dX = theta*(mu - X) dt + sigma dW on [0, T] starting at x0, sampled… |
| `func PoissonCountingProcess(rng *rand.Rand, rate, T float64, nSteps int) Path` | PoissonCountingProcess simulates a homogeneous Poisson process of the given rate on [0, T] and returns its counting path sampled on nSteps+1 grid… |
| `func SimpleRandomWalkPath(rng *rand.Rand, n int, p float64) Path` | SimpleRandomWalkPath returns a simple random walk as a Path with unit time steps. |
| `func SimulateEnsemble(seed int64, sde SDE, x0, T float64, nSteps, paths int) []Path` | SimulateEnsemble runs EulerMaruyama for the given SDE paths times, each with an independent stream derived from the base seed, and returns the paths. |
| `func StandardBrownianBridge(rng *rand.Rand, nSteps int) Path` | StandardBrownianBridge returns a standard Brownian bridge pinned to 0 at both endpoints of [0, 1]. |
| `func StratonovichHeun(rng *rand.Rand, sde SDE, x0, T float64, nSteps int) Path` | StratonovichHeun integrates the Stratonovich SDE dX = a dt + b o dW with the stochastic Heun (predictor-corrector) scheme, which converges to the… |
| `func TamedEulerMaruyama(rng *rand.Rand, sde SDE, x0, T float64, nSteps int) Path` | TamedEulerMaruyama integrates the SDE with the tamed Euler scheme, which bounds the drift term to remain stable for super-linearly growing drift. |
| `func VasicekProcess(rng *rand.Rand, r0, a, b, sigma, T float64, nSteps int) Path` | VasicekProcess is an alias for OrnsteinUhlenbeck using interest-rate terminology: dr = a*(b - r) dt + sigma dW. |
| `func (p Path) Antithetic() Path` | Antithetic returns the path reflected about its starting value, that is the path with value v(0) - (v(t) - v(0)) = 2*v(0) - v(t). |
| `func (p Path) ArgMax() int` | ArgMax returns the index of the first maximal value along the path. |
| `func (p Path) ArgMin() int` | ArgMin returns the index of the first minimal value along the path. |
| `func (p Path) At(t float64) float64` | At returns the value of the path at time t using linear interpolation between grid points. |
| `func (p Path) Copy() Path` | Copy returns a deep copy of the path. |
| `func (p Path) CrossesLevel(level float64) bool` | CrossesLevel reports whether the path ever reaches or crosses the given level relative to its starting value. |
| `func (p Path) Empty() bool` | Empty reports whether the path has no points. |
| `func (p Path) EndTime() float64` | EndTime returns the last time of the path, or 0 for an empty path. |
| `func (p Path) Final() float64` | Final returns the last value of the path, or 0 for an empty path. |
| `func (p Path) Increments() []float64` | Increments returns the successive differences of the path values. |
| `func (p Path) IndexToHit(level float64) (int, bool)` | IndexToHit returns the first index at which the path reaches or crosses the given level, measured relative to the starting value. |
| `func (p Path) Initial() float64` | Initial returns the first value of the path, or 0 for an empty path. |
| `func (p Path) Len() int` | Len returns the number of points in the path. |
| `func (p Path) Max() float64` | Max returns the maximum value along the path. |
| `func (p Path) MaxDrawdown() float64` | MaxDrawdown returns the largest peak-to-trough drop of the path, defined as the maximum over t of running-max minus value. |
| `func (p Path) Mean() float64` | Mean returns the arithmetic mean of the path values. |
| `func (p Path) Min() float64` | Min returns the minimum value along the path. |
| `func (p Path) OccupationFraction(lo, hi float64) float64` | OccupationFraction returns the fraction of grid points whose value lies in the closed interval [lo, hi]. |
| `func (p Path) QuadraticVariation() float64` | QuadraticVariation returns the sum of squared increments, an estimate of the accumulated quadratic variation of the path. |
| `func (p Path) Range() float64` | Range returns the difference between the maximum and minimum values. |
| `func (p Path) RunningMax() []float64` | RunningMax returns the running maximum of the path values. |
| `func (p Path) RunningMean() []float64` | RunningMean returns the running arithmetic mean of the path values. |
| `func (p Path) RunningMin() []float64` | RunningMin returns the running minimum of the path values. |
| `func (p Path) Scale(c float64) Path` | Scale returns a new path whose values are multiplied by c. |
| `func (p Path) Shift(c float64) Path` | Shift returns a new path whose values are shifted by the additive constant c. |
| `func (p Path) StartTime() float64` | StartTime returns the first time of the path, or 0 for an empty path. |
| `func (p Path) StdDev() float64` | StdDev returns the sample standard deviation of the path values. |
| `func (p Path) TimeToHit(level float64) (float64, bool)` | TimeToHit returns the time of the first index at which the path reaches or crosses the given level. |
| `func (p Path) TotalVariation() float64` | TotalVariation returns the sum of the absolute increments of the path. |
| `func (p Path) Variance() float64` | Variance returns the sample variance (with Bessel's correction) of the path values. |

</details>

<details>
<summary><code>Reaction</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func MassActionReaction(stoich, reactants []int, k float64) Reaction` | MassActionReaction builds a mass-action Reaction with the given stoichiometry, reactant multiplicities and rate constant. |
| `func NewReaction(stoich []int, propensity func(state []int) float64) Reaction` | NewReaction builds a Reaction from an explicit stoichiometry vector and propensity function. |

</details>

<details>
<summary><code>ReactionNetwork</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BirthDeathNetwork(birth, death float64) ReactionNetwork` | BirthDeathNetwork returns a one-species linear birth-death network with per-capita birth rate birth and per-capita death rate death. |
| `func ImmigrationDeathNetwork(immigration, death float64) ReactionNetwork` | ImmigrationDeathNetwork returns a one-species network with constant immigration rate immigration and per-capita death rate death. |
| `func LotkaVolterraNetwork(alpha, beta, gamma float64) ReactionNetwork` | LotkaVolterraNetwork returns the stochastic Lotka-Volterra predator-prey network on species (prey, predator) with prey birth rate alpha (per prey),… |
| `func SIRNetwork(beta, gamma float64) ReactionNetwork` | SIRNetwork returns the stochastic SIR epidemic network on species (S, I, R) with infection rate constant beta and recovery rate gamma. |
| `func (net ReactionNetwork) Propensities(state []int) []float64` | Propensities returns the propensity of each reaction at the given state. |
| `func (net ReactionNetwork) TotalPropensity(state []int) float64` | TotalPropensity returns the sum of the propensities of all reactions in the network at the given state. |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/stochastic).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
