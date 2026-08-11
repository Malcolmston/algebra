# reliability — provides reliability engineering and survival-analysis primitives built entirely

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/reliability.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/reliability)

Package reliability provides reliability engineering and survival-analysis
primitives built entirely on the Go standard library.

The package covers four broad areas.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/reliability"
```

## Usage

This is the package's own `ExampleAvailability_repairableItem`, so it compiles and its output is
asserted on every `go test ./reliability/`.

```go
a := reliability.Availability(1000, 10)
	fmt.Printf("availability   = %.9f\n", a)
	fmt.Printf("unavailability = %.9f\n", reliability.Unavailability(1000, 10))
	fmt.Printf("downtime/year  = %.6f h\n", reliability.ExpectedDowntime(a, 8760))
```

```
availability   = 0.990099010
unavailability = 0.009900990
downtime/year  = 86.732673 h
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AchievedAvailability(mtbm, mmt float64) float64` | AchievedAvailability returns MTBM/(MTBM+MMT), the availability based on the mean time between maintenance actions (both corrective and preventive)… |
| `func AllocateFailureRateSeriesEqual(lambdaTotal float64, n int) float64` | AllocateFailureRateSeriesEqual apportions an overall series-system failure rate equally across n subsystems, returning the per-subsystem rate… |
| `func AllocateReliabilitySeriesEqual(target float64, n int) float64` | AllocateReliabilitySeriesEqual applies the equal-apportionment technique to a series system: to achieve an overall target reliability across n… |
| `func Availability(mtbf, mttr float64) float64` | Availability returns the steady-state availability MTBF/(MTBF+MTTR) of a repairable item from its mean time between failures and mean time to repair. |
| `func AvailabilityFromRates(lambda, mu float64) float64` | AvailabilityFromRates returns the steady-state availability from a failure rate lambda and repair rate mu, an alias of SteadyStateAvailability. |
| `func AvailabilityKofN(k, n int, a float64) float64` | AvailabilityKofN returns the availability of a k-out-of-n system of identical units each with availability a, using the same binomial sum as… |
| `func AvailabilityParallelSystem(avails ...float64) float64` | AvailabilityParallelSystem returns the availability of a parallel system, 1-∏(1-aᵢ). |
| `func AvailabilitySeriesSystem(avails ...float64) float64` | AvailabilitySeriesSystem returns the availability of a series system, the product of its component availabilities. |
| `func BathtubCDF(t, a, b, c, e float64) float64` | BathtubCDF returns the failure probability F(t)=1-R(t) of the additive bathtub model. |
| `func BathtubCumulativeHazard(t, a, b, c, e float64) float64` | BathtubCumulativeHazard returns the cumulative hazard of the additive bathtub model, |
| `func BathtubHazard(t, a, b, c, e float64) float64` | BathtubHazard evaluates an additive bathtub-shaped hazard rate at time t>=0: |
| `func BathtubPDF(t, a, b, c, e float64) float64` | BathtubPDF returns the probability density h(t)R(t) of the additive bathtub model. |
| `func BathtubReliability(t, a, b, c, e float64) float64` | BathtubReliability returns the survival probability R(t)=exp(-H(t)) of the additive bathtub model. |
| `func BridgeReliability(r1, r2, r3, r4, r5 float64) float64` | BridgeReliability returns the exact reliability of the classic five-component bridge network, evaluated by conditioning (pivotal decomposition) on… |
| `func ChiSquarePValue1DF(x float64) float64` | ChiSquarePValue1DF returns the upper-tail probability P(X>=x) for a chi-square random variable with one degree of freedom, useful for interpreting a… |
| `func CriticalityImportance(birnbaum, rComp, rSys float64) float64` | CriticalityImportance returns the criticality (Lambert) importance of a component: its Birnbaum importance weighted by its own unreliability and… |
| `func DesignLifeForReliability(rTarget, lambda float64) float64` | DesignLifeForReliability returns the maximum operating time at which a constant-failure-rate item still meets a target reliability rTarget:… |
| `func EmpiricalCDF(times []float64, t float64) float64` | EmpiricalCDF returns the fraction of a sample of failure times that are less than or equal to t, an estimate of F(t) for complete data. |
| `func EmpiricalHazardRate(times []float64, t, dt float64) float64` | EmpiricalHazardRate estimates the average hazard rate over the interval (t, t+dt] from complete data as (survivors past t that fail within dt)… |
| `func EmpiricalReliability(times []float64, t float64) float64` | EmpiricalReliability returns the fraction of a sample of failure times that exceed t, an unbiased estimate of R(t) for complete data. |
| `func ExpectedDowntime(availability, mission float64) float64` | ExpectedDowntime returns the expected downtime over a mission of the given duration for an item with the supplied availability:… |
| `func ExpectedUptime(availability, mission float64) float64` | ExpectedUptime returns the expected uptime over a mission of the given duration for an item with the supplied availability: availability·mission. |
| `func ExponentialCDF(t, lambda float64) float64` | ExponentialCDF returns the failure probability F(t)=1-e^{-λt} for the exponential distribution with failure rate lambda>0 at time t>=0. |
| `func ExponentialConditionalReliability(t, x, lambda float64) float64` | ExponentialConditionalReliability returns the probability that a component surviving to age t survives an additional time x. |
| `func ExponentialCumulativeHazard(t, lambda float64) float64` | ExponentialCumulativeHazard returns the cumulative hazard H(t)=λt of the exponential distribution. |
| `func ExponentialHazard(t, lambda float64) float64` | ExponentialHazard returns the (constant) hazard rate h(t)=λ of the exponential distribution. |
| `func ExponentialMTTF(lambda float64) float64` | ExponentialMTTF returns the mean time to failure 1/λ, identical to the mean. |
| `func ExponentialMean(lambda float64) float64` | ExponentialMean returns the mean lifetime 1/λ. |
| `func ExponentialMeanResidualLife(t, lambda float64) float64` | ExponentialMeanResidualLife returns the mean remaining life given survival to age t, which is the constant 1/λ for the exponential distribution. |
| `func ExponentialMedian(lambda float64) float64` | ExponentialMedian returns the median lifetime ln(2)/λ. |
| `func ExponentialPDF(t, lambda float64) float64` | ExponentialPDF returns the probability density f(t)=λe^{-λt} of the exponential lifetime distribution with failure rate lambda>0 at time t>=0. |
| `func ExponentialQuantile(p, lambda float64) float64` | ExponentialQuantile returns the time t at which F(t)=p, i.e. |
| `func ExponentialRateFromMTTF(mttf float64) float64` | ExponentialRateFromMTTF recovers the failure rate λ=1/MTTF. |
| `func ExponentialRateFromReliability(t, r float64) float64` | ExponentialRateFromReliability recovers the failure rate λ from a reliability observation R at time t>0: λ=-ln(R)/t. |
| `func ExponentialReliability(t, lambda float64) float64` | ExponentialReliability returns the survival probability R(t)=e^{-λt} for the exponential distribution with failure rate lambda>0 at time t>=0. |
| `func ExponentialStdDev(lambda float64) float64` | ExponentialStdDev returns the standard deviation 1/λ of the exponential distribution. |
| `func ExponentialVariance(lambda float64) float64` | ExponentialVariance returns the variance 1/λ² of the exponential distribution. |
| `func FIT(lambdaPerHour float64) float64` | FIT returns the failure-in-time rate: the number of failures expected per one billion (10^9) device-hours for a constant hourly failure rate. |
| `func FITToRate(fit float64) float64` | FITToRate converts a failure-in-time value (failures per 10^9 hours) back to an hourly failure rate. |
| `func FailureRate(f, r float64) float64` | FailureRate returns the instantaneous hazard rate f/R from a probability density f and a reliability (survival) value r. |
| `func FailureRateFromMTBF(mtbf float64) float64` | FailureRateFromMTBF returns the constant failure rate 1/MTBF. |
| `func FailuresPerMillionHours(lambdaPerHour float64) float64` | FailuresPerMillionHours converts a failure rate expressed per unit time into the FIT-like measure of expected failures per one million hours. |
| `func GammaCDF(t, k, theta float64) float64` | GammaCDF returns the failure probability F(t)=P(k,t/θ), the regularized lower incomplete gamma function. |
| `func GammaCumulativeHazard(t, k, theta float64) float64` | GammaCumulativeHazard returns the cumulative hazard H(t)=-ln R(t). |
| `func GammaHazard(t, k, theta float64) float64` | GammaHazard returns the hazard rate f(t)/R(t) of the gamma distribution. |
| `func GammaMTTF(k, theta float64) float64` | GammaMTTF returns the mean time to failure, identical to the mean. |
| `func GammaMean(k, theta float64) float64` | GammaMean returns the mean lifetime kθ. |
| `func GammaMeanResidualLife(t, k, theta float64) float64` | GammaMeanResidualLife returns the mean remaining life given survival to age t, computed by numerically integrating the reliability from t to infinity. |
| `func GammaMode(k, theta float64) float64` | GammaMode returns the mode of the gamma distribution. |
| `func GammaPDF(t, k, theta float64) float64` | GammaPDF returns the probability density of the gamma lifetime distribution with shape k>0 and scale theta>0 at time t>=0:… |
| `func GammaQuantile(p, k, theta float64) float64` | GammaQuantile returns the time at which F(t)=p, computed by inverting the regularized incomplete gamma function. |
| `func GammaReliability(t, k, theta float64) float64` | GammaReliability returns the survival probability R(t)=Q(k,t/θ), the regularized upper incomplete gamma function. |
| `func GammaStdDev(k, theta float64) float64` | GammaStdDev returns the standard deviation of the gamma distribution. |
| `func GammaVariance(k, theta float64) float64` | GammaVariance returns the variance kθ². |
| `func GompertzCDF(t, a, b float64) float64` | GompertzCDF returns the failure probability F(t)=1-R(t). |
| `func GompertzCumulativeHazard(t, a, b float64) float64` | GompertzCumulativeHazard returns the cumulative hazard H(t)=(a/b)(e^{bt}-1) of the Gompertz distribution. |
| `func GompertzHazard(t, a, b float64) float64` | GompertzHazard returns the exponentially increasing hazard rate h(t)=a·e^{bt} of the Gompertz lifetime distribution, with a>0 and b>0. |
| `func GompertzMakehamCDF(t, lambda, a, b float64) float64` | GompertzMakehamCDF returns the failure probability F(t)=1-R(t). |
| `func GompertzMakehamCumulativeHazard(t, lambda, a, b float64) float64` | GompertzMakehamCumulativeHazard returns the cumulative hazard H(t)=λt+(a/b)(e^{bt}-1). |
| `func GompertzMakehamHazard(t, lambda, a, b float64) float64` | GompertzMakehamHazard returns the hazard rate h(t)=λ+a·e^{bt} of the Gompertz–Makeham law, which adds an age-independent Makeham term lambda>=0… |
| `func GompertzMakehamPDF(t, lambda, a, b float64) float64` | GompertzMakehamPDF returns the probability density h(t)R(t). |
| `func GompertzMakehamReliability(t, lambda, a, b float64) float64` | GompertzMakehamReliability returns the survival probability R(t)=exp(-H(t)) of the Gompertz–Makeham law. |
| `func GompertzMean(a, b float64) float64` | GompertzMean returns the mean lifetime, computed by numerically integrating the reliability from 0 to infinity. |
| `func GompertzMedian(a, b float64) float64` | GompertzMedian returns the median lifetime, the 0.5 quantile. |
| `func GompertzPDF(t, a, b float64) float64` | GompertzPDF returns the probability density h(t)R(t) of the Gompertz distribution. |
| `func GompertzQuantile(p, a, b float64) float64` | GompertzQuantile returns the time at which F(t)=p: t=(1/b)ln(1-(b/a)ln(1-p)). |
| `func GompertzReliability(t, a, b float64) float64` | GompertzReliability returns the survival probability R(t)=exp(-(a/b)(e^{bt}-1)). |
| `func HazardFromReliability(r float64) float64` | HazardFromReliability converts a reliability value r into its cumulative hazard H=-ln(r). |
| `func HazardRatio(h1, h2 float64) float64` | HazardRatio returns the ratio of two hazard rates h1/h2, the constant proportional-hazards multiplier of group 1 relative to group 2. |
| `func HjorthCDF(t, delta, beta, theta float64) float64` | HjorthCDF returns the failure probability F(t)=1-R(t) of the Hjorth model. |
| `func HjorthCumulativeHazard(t, delta, beta, theta float64) float64` | HjorthCumulativeHazard returns the cumulative hazard of the Hjorth model, |
| `func HjorthHazard(t, delta, beta, theta float64) float64` | HjorthHazard evaluates the Hjorth-distribution hazard rate at time t>=0: |
| `func HjorthPDF(t, delta, beta, theta float64) float64` | HjorthPDF returns the probability density h(t)R(t) of the Hjorth model. |
| `func HjorthReliability(t, delta, beta, theta float64) float64` | HjorthReliability returns the survival probability R(t)=exp(-δt²/2)·(1+θt)^{-β/θ}. |
| `func InherentAvailability(mttf, mttr float64) float64` | InherentAvailability returns MTTF/(MTTF+MTTR), the availability that accounts only for corrective maintenance under idealized support conditions. |
| `func InstantaneousAvailability(t, lambda, mu float64) float64` | InstantaneousAvailability returns the availability at time t of a two-state Markov model that starts operational:… |
| `func KofNReliability(k, n int, r float64) float64` | KofNReliability returns the reliability of a k-out-of-n system of identical components, each with reliability r, that works when at least k of its n… |
| `func KofNReliabilityGeneral(k int, rs []float64) (float64, error)` | KofNReliabilityGeneral returns the exact reliability of a k-out-of-n system whose n components may have different reliabilities rs. |
| `func LognormalCDF(t, mu, sigma float64) float64` | LognormalCDF returns the failure probability F(t)=Φ((ln t-μ)/σ). |
| `func LognormalCumulativeHazard(t, mu, sigma float64) float64` | LognormalCumulativeHazard returns the cumulative hazard H(t)=-ln R(t). |
| `func LognormalHazard(t, mu, sigma float64) float64` | LognormalHazard returns the hazard rate f(t)/R(t) of the lognormal distribution, which rises from 0, peaks, and then decreases. |
| `func LognormalMTTF(mu, sigma float64) float64` | LognormalMTTF returns the mean time to failure, identical to the mean. |
| `func LognormalMean(mu, sigma float64) float64` | LognormalMean returns the mean lifetime exp(μ+σ²/2). |
| `func LognormalMeanResidualLife(t, mu, sigma float64) float64` | LognormalMeanResidualLife returns the mean remaining life given survival to age t, computed by numerically integrating the reliability from t to… |
| `func LognormalMedian(mu, sigma float64) float64` | LognormalMedian returns the median lifetime e^{μ}. |
| `func LognormalMode(mu, sigma float64) float64` | LognormalMode returns the mode of the lognormal distribution, e^{μ-σ²}. |
| `func LognormalPDF(t, mu, sigma float64) float64` | LognormalPDF returns the probability density of the lognormal lifetime distribution whose logarithm is normal with mean mu and standard deviation… |
| `func LognormalQuantile(p, mu, sigma float64) float64` | LognormalQuantile returns the time at which F(t)=p: t=exp(μ+σΦ^{-1}(p)). |
| `func LognormalReliability(t, mu, sigma float64) float64` | LognormalReliability returns the survival probability R(t)=1-Φ((ln t-μ)/σ). |
| `func LognormalStdDev(mu, sigma float64) float64` | LognormalStdDev returns the standard deviation of the lognormal distribution. |
| `func LognormalVariance(mu, sigma float64) float64` | LognormalVariance returns the variance (e^{σ²}-1)e^{2μ+σ²}. |
| `func MTBF(mttf, mttr float64) float64` | MTBF returns the mean time between failures MTTF+MTTR of a repairable item. |
| `func MTBFFromData(totalOperatingTime float64, failures int) float64` | MTBFFromData estimates the mean time between failures as the total operating time divided by the number of failures observed. |
| `func MTBFFromFailureRate(lambda float64) float64` | MTBFFromFailureRate returns the mean time between failures 1/λ for a constant failure rate lambda>0. |
| `func MTTFFromData(times []float64) (float64, error)` | MTTFFromData estimates the mean time to failure as the arithmetic mean of a set of (complete, uncensored) failure times. |
| `func MTTFFromReliability(r func(float64) float64, scale float64) float64` | MTTFFromReliability computes the mean time to failure by numerically integrating a reliability function R over [0,∞). |
| `func MTTR(mu float64) float64` | MTTR returns the mean time to repair 1/μ from a constant repair rate mu>0. |
| `func MaxSeriesForReliability(target, r float64) int` | MaxSeriesForReliability returns the largest number of identical series components, each with reliability r in (0,1], whose product reliability still… |
| `func MeanRank(i, n int) float64` | MeanRank returns the mean plotting position i/(n+1) of the i-th of n ordered failures. |
| `func MedianRank(i, n int) float64` | MedianRank returns Bernard's approximation to the median rank of the i-th of n ordered failures, (i-0.3)/(n+0.4), commonly used to plot empirical… |
| `func MinParallelForReliability(target, r float64) int` | MinParallelForReliability returns the smallest number of identical parallel components, each with reliability r in (0,1), needed to reach a system… |
| `func NormalCDF(t, mu, sigma float64) float64` | NormalCDF returns the failure probability F(t)=Φ((t-μ)/σ). |
| `func NormalCumulativeHazard(t, mu, sigma float64) float64` | NormalCumulativeHazard returns the cumulative hazard H(t)=-ln R(t). |
| `func NormalHazard(t, mu, sigma float64) float64` | NormalHazard returns the hazard rate f(t)/R(t) of the normal distribution, which increases monotonically with t. |
| `func NormalMTTF(mu, sigma float64) float64` | NormalMTTF returns the mean time to failure, equal to mu for the normal model. |
| `func NormalPDF(t, mu, sigma float64) float64` | NormalPDF returns the probability density of the normal lifetime distribution with mean mu and standard deviation sigma>0 at time t. |
| `func NormalQuantile(p, mu, sigma float64) float64` | NormalQuantile returns the time at which F(t)=p: t=μ+σΦ^{-1}(p). |
| `func NormalReliability(t, mu, sigma float64) float64` | NormalReliability returns the survival probability R(t)=1-Φ((t-μ)/σ). |
| `func OperationalAvailability(uptime, downtime float64) float64` | OperationalAvailability returns uptime/(uptime+downtime), the availability observed in operation including all sources of downtime. |
| `func ParallelFailureProbability(rs ...float64) float64` | ParallelFailureProbability returns the probability that a parallel system fails, the product ∏(1-rᵢ). |
| `func ParallelReliability(rs ...float64) float64` | ParallelReliability returns the reliability of a parallel (active-redundant) system, 1-∏(1-rᵢ). |
| `func ParallelReliabilityIdentical(r float64, n int) float64` | ParallelReliabilityIdentical returns the reliability 1-(1-r)ⁿ of a parallel system of n identical components each with reliability r. |
| `func ParallelSeriesReliability(strings [][]float64) float64` | ParallelSeriesReliability returns the reliability of a system built from series strings connected in parallel: each inner slice is a series string… |
| `func ParallelSystemMTTF(lambda float64, n int) float64` | ParallelSystemMTTF returns the mean time to failure of a parallel system of n identical exponential units with failure rate lambda, which is the… |
| `func ProportionalHazardsReliability(r0, hr float64) float64` | ProportionalHazardsReliability returns the reliability of an item whose hazard is a constant multiple hr of a baseline hazard: R(t)=R0(t)^{hr}, where… |
| `func RayleighCDF(t, sigma float64) float64` | RayleighCDF returns the failure probability F(t)=1-e^{-t²/(2σ²)}. |
| `func RayleighCumulativeHazard(t, sigma float64) float64` | RayleighCumulativeHazard returns the cumulative hazard H(t)=t²/(2σ²). |
| `func RayleighHazard(t, sigma float64) float64` | RayleighHazard returns the linearly increasing hazard rate h(t)=t/σ². |
| `func RayleighMTTF(sigma float64) float64` | RayleighMTTF returns the mean time to failure, identical to the mean. |
| `func RayleighMean(sigma float64) float64` | RayleighMean returns the mean lifetime σ√(π/2). |
| `func RayleighMedian(sigma float64) float64` | RayleighMedian returns the median lifetime σ√(2ln2). |
| `func RayleighMode(sigma float64) float64` | RayleighMode returns the mode of the Rayleigh distribution, equal to σ. |
| `func RayleighPDF(t, sigma float64) float64` | RayleighPDF returns the probability density of the Rayleigh lifetime distribution with scale sigma>0 at time t>=0: f(t)=(t/σ²)e^{-t²/(2σ²)}. |
| `func RayleighQuantile(p, sigma float64) float64` | RayleighQuantile returns the time at which F(t)=p: t=σ√(-2ln(1-p)). |
| `func RayleighReliability(t, sigma float64) float64` | RayleighReliability returns the survival probability R(t)=e^{-t²/(2σ²)}. |
| `func RayleighStdDev(sigma float64) float64` | RayleighStdDev returns the standard deviation of the Rayleigh distribution. |
| `func RayleighVariance(sigma float64) float64` | RayleighVariance returns the variance (2-π/2)σ². |
| `func RedundancyGain(r float64) float64` | RedundancyGain returns the improvement in reliability obtained by placing a single spare component of reliability r in active parallel with an… |
| `func ReliabilityAtTime(t, lambda float64) float64` | ReliabilityAtTime returns the reliability of a constant-failure-rate item at time t: R(t)=e^{-λt}. |
| `func ReliabilityFromHazard(h float64) float64` | ReliabilityFromHazard converts a cumulative hazard H>=0 into a reliability value R=e^{-H}. |
| `func ReliabilityImportanceBirnbaum(rWorking, rFailed float64) float64` | ReliabilityImportanceBirnbaum returns the Birnbaum reliability importance of a component in a system, the partial derivative of system reliability… |
| `func RepairRateFromMTTR(mttr float64) float64` | RepairRateFromMTTR returns the repair rate μ=1/MTTR. |
| `func SeriesFailureProbability(rs ...float64) float64` | SeriesFailureProbability returns the probability that a series system fails, 1 minus its reliability. |
| `func SeriesParallelReliability(groups [][]float64) float64` | SeriesParallelReliability returns the reliability of a system built from parallel groups connected in series: each inner slice is a parallel… |
| `func SeriesReliability(rs ...float64) float64` | SeriesReliability returns the reliability of a series system, the product of the component reliabilities. |
| `func SeriesReliabilityIdentical(r float64, n int) float64` | SeriesReliabilityIdentical returns the reliability rⁿ of a series system of n identical components each with reliability r. |
| `func SeriesSystemFailureRate(lambdas ...float64) float64` | SeriesSystemFailureRate returns the total failure rate Σλᵢ of a series system of independent constant-rate components. |
| `func SeriesSystemMTTF(lambdas ...float64) float64` | SeriesSystemMTTF returns the mean time to failure 1/Σλᵢ of a series system of independent constant-rate components. |
| `func SortedFailureTimes(times []float64) []float64` | SortedFailureTimes returns a sorted copy of a slice of failure times. |
| `func StandardNormalCDF(x float64) float64` | StandardNormalCDF returns Φ(x), the standard-normal cumulative distribution function. |
| `func StandardNormalQuantile(p float64) float64` | StandardNormalQuantile returns Φ^{-1}(p), the standard-normal quantile (probit) function. |
| `func StandbyRedundancyMTTF(lambda float64, n int) float64` | StandbyRedundancyMTTF returns the mean time to failure n/λ of a cold-standby system of n identical exponential units with per-unit failure rate… |
| `func StandbyRedundancyReliability(t, lambda float64, n int) float64` | StandbyRedundancyReliability returns the reliability at time t of a system of n identical exponential units (one active, n-1 in cold standby with… |
| `func SteadyStateAvailability(lambda, mu float64) float64` | SteadyStateAvailability returns the limiting availability μ/(λ+μ) of a two-state (up/down) Markov model with constant failure rate lambda and… |
| `func Unavailability(mtbf, mttr float64) float64` | Unavailability returns 1 minus the availability computed from MTBF and MTTR. |
| `func WeibullB10Life(k, lambda float64) float64` | WeibullB10Life returns the B10 life, the time by which 10% of the population has failed (the 0.10 quantile). |
| `func WeibullCDF(t, k, lambda float64) float64` | WeibullCDF returns the failure probability F(t)=1-e^{-(t/λ)^k}. |
| `func WeibullCharacteristicLife(lambda float64) float64` | WeibullCharacteristicLife returns the scale parameter λ, the time by which 63.2% of a Weibull-distributed population has failed. |
| `func WeibullConditionalReliability(t, x, k, lambda float64) float64` | WeibullConditionalReliability returns the probability that a component of age t survives an additional time x: R(t+x)/R(t). |
| `func WeibullCumulativeHazard(t, k, lambda float64) float64` | WeibullCumulativeHazard returns the cumulative hazard H(t)=(t/λ)^k. |
| `func WeibullHazard(t, k, lambda float64) float64` | WeibullHazard returns the hazard rate h(t)=(k/λ)(t/λ)^{k-1}. |
| `func WeibullMTTF(k, lambda float64) float64` | WeibullMTTF returns the mean time to failure, identical to the mean. |
| `func WeibullMean(k, lambda float64) float64` | WeibullMean returns the mean lifetime λΓ(1+1/k). |
| `func WeibullMeanResidualLife(t, k, lambda float64) float64` | WeibullMeanResidualLife returns the mean remaining life given survival to age t, computed by numerically integrating the reliability from t to… |
| `func WeibullMedian(k, lambda float64) float64` | WeibullMedian returns the median lifetime λ(ln2)^{1/k}. |
| `func WeibullMode(k, lambda float64) float64` | WeibullMode returns the mode of the Weibull distribution. |
| `func WeibullPDF(t, k, lambda float64) float64` | WeibullPDF returns the probability density of the two-parameter Weibull distribution with shape k>0 and scale lambda>0 at time t>=0:… |
| `func WeibullQuantile(p, k, lambda float64) float64` | WeibullQuantile returns the time at which F(t)=p: t=λ(-ln(1-p))^{1/k}. |
| `func WeibullReliability(t, k, lambda float64) float64` | WeibullReliability returns the survival probability R(t)=e^{-(t/λ)^k}. |
| `func WeibullStdDev(k, lambda float64) float64` | WeibullStdDev returns the standard deviation of the Weibull distribution. |
| `func WeibullVariance(k, lambda float64) float64` | WeibullVariance returns the variance λ²[Γ(1+2/k)-Γ(1+1/k)²]. |

### Types

| Type | What it is |
| --- | --- |
| `KMEstimate` | KMEstimate holds the result of a Kaplan–Meier product-limit survival estimate. |
| `LogRankResult` | LogRankResult holds the outcome of a two-sample log-rank test comparing the survival of group 1 against group 2. |
| `NAEstimate` | NAEstimate holds the result of a Nelson–Aalen cumulative-hazard estimate. |

<details>
<summary><code>KMEstimate</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func KaplanMeier(times []float64, events []bool) (*KMEstimate, error)` | KaplanMeier computes the Kaplan–Meier product-limit estimator of the survival function from right-censored data. |
| `func (e *KMEstimate) CumulativeHazardAt(t float64) float64` | CumulativeHazardAt returns the Kaplan–Meier-derived cumulative hazard -ln S(t). |
| `func (e *KMEstimate) MedianSurvival() float64` | MedianSurvival returns the smallest failure time at which the estimated survival drops to 0.5 or below. |
| `func (e *KMEstimate) Quantile(p float64) float64` | Quantile returns the smallest failure time at which the estimated survival drops to 1-p or below, i.e. |
| `func (e *KMEstimate) RestrictedMeanSurvival(tau float64) float64` | RestrictedMeanSurvival returns the restricted mean survival time over [0,tau], the area under the Kaplan–Meier curve up to the horizon tau. |
| `func (e *KMEstimate) SurvivalAt(t float64) float64` | SurvivalAt returns the Kaplan–Meier survival estimate S(t). |

</details>

<details>
<summary><code>LogRankResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func LogRankTest(times1 []float64, events1 []bool, times2 []float64, events2 []bool) (*LogRankResult, error)` | LogRankTest performs the two-sample log-rank test on right-censored survival data. |

</details>

<details>
<summary><code>NAEstimate</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NelsonAalen(times []float64, events []bool) (*NAEstimate, error)` | NelsonAalen computes the Nelson–Aalen estimator of the cumulative hazard function from right-censored data, together with its variance. |
| `func (e *NAEstimate) CumulativeHazardAt(t float64) float64` | CumulativeHazardAt returns the Nelson–Aalen cumulative hazard estimate H(t), a right-continuous step function equal to 0 before the first failure. |
| `func (e *NAEstimate) SurvivalAt(t float64) float64` | SurvivalAt returns the survival estimate exp(-H(t)) implied by the Nelson–Aalen cumulative hazard. |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/reliability).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
