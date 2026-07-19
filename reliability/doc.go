// Package reliability provides reliability engineering and survival-analysis
// primitives built entirely on the Go standard library.
//
// The package covers four broad areas.
//
// Lifetime distributions: closed-form probability density, cumulative
// distribution, reliability (survival), hazard-rate, cumulative-hazard,
// quantile and summary-statistic functions for the exponential, Weibull,
// lognormal, gamma, Rayleigh, normal, Gompertz and Gompertz–Makeham models,
// together with mean-time-to-failure (MTTF), conditional reliability and
// mean-residual-life helpers. A flexible additive [BathtubHazard] model and a
// two-term [HjorthHazard] model reproduce the classic decreasing/constant/
// increasing "bathtub" failure-rate curve.
//
// System reliability: reliability-block-diagram combinators for series,
// parallel (active redundancy), k-out-of-n and bridge structures, exact
// evaluation for non-identical components, standby-redundancy formulas and
// redundancy-allocation helpers such as [MinParallelForReliability].
//
// Availability: inherent, operational and steady-state availability from
// MTBF/MTTF/MTTR data or from failure/repair rates, plus availability of
// series and parallel repairable systems.
//
// Survival analysis: the [KaplanMeier] product-limit survival estimator with
// Greenwood variance, the [NelsonAalen] cumulative-hazard estimator, median
// and restricted-mean survival, log-rank testing and hazard-ratio estimation
// from right-censored data.
//
// Unless stated otherwise, times and rates must be non-negative and shape and
// scale parameters must be strictly positive; functions return NaN for
// arguments outside their natural domain rather than panicking. All routines
// are deterministic and depend only on the standard library.
package reliability
