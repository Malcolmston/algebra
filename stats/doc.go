// Package stats is a small, dependency-free statistics and probability
// toolkit written entirely with the Go standard library (only math and sort).
// It operates on plain []float64 data and is completely standalone: it does
// not import the parent algebra package.
//
// The package is organised into four groups.
//
// # Descriptive statistics
//
// Summaries of a data set: [Mean], [Median], [Mode], [Variance] and
// [PopVariance], [StdDev] and [PopStdDev], [Min], [Max], [Range], [Sum],
// [Product], [Quantile], [Percentile], [IQR], [Skewness], [Kurtosis],
// [GeometricMean], [HarmonicMean], [Covariance], [Correlation], [ZScore] and
// [WeightedMean].
//
// Functions that need a defined result on empty or too-small inputs return
// NaN rather than panicking, so callers can test with [math.IsNaN].
//
// # Combinatorics
//
// [Factorial], [Choose] (nCr) and [Perm] (nPr) are computed in the log-gamma
// domain so they stay finite and accurate for large arguments.
//
// # Probability distributions
//
// Each distribution is a small value struct carrying its parameters, with
// methods for the density/mass function, the cumulative distribution
// function, the mean and the variance, and (where feasible) a quantile:
//
//   - [Normal]      — PDF, CDF, Quantile, Mean, Variance.
//   - [Binomial]    — PMF, CDF, Mean, Variance.
//   - [Poisson]     — PMF, CDF, Mean, Variance.
//   - [Uniform]     — PDF, CDF, Quantile, Mean, Variance.
//   - [Exponential] — PDF, CDF, Quantile, Mean, Variance.
//   - [StudentT]    — PDF, CDF, Mean, Variance.
//   - [ChiSquared]  — PDF, CDF, Mean, Variance.
//   - [Gamma]       — PDF, CDF, Mean, Variance.
//
// [NormalPDF] and [NormalCDF] are convenience wrappers around the [Normal]
// methods. No random sampling is provided.
//
// # Regression
//
// [LinearRegression] fits an ordinary-least-squares line and reports its
// slope, intercept and coefficient of determination (R²).
//
// # Conventions
//
// The sample variants ([Variance], [StdDev], [Covariance]) use the unbiased
// n-1 (Bessel) denominator; the population variants use n. [Skewness] and
// [Kurtosis] use the population moment estimators, and [Kurtosis] returns
// excess kurtosis (0 for a normal distribution).
package stats
