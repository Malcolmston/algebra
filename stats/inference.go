package stats

import (
	"math"
	"math/rand"
)

// Quantile returns the p-quantile (inverse CDF) of the Student's t
// distribution: the value x such that CDF(x) = p, for p in (0, 1). It is
// computed by monotone bisection on the (strictly increasing) CDF, starting
// from a small symmetric bracket that is doubled outward until it straddles p.
// It returns -Inf for p <= 0, +Inf for p >= 1, and NaN when p is NaN. This
// single definition underlies the small-sample confidence intervals in this
// file.
func (t StudentT) Quantile(p float64) float64 {
	if math.IsNaN(p) {
		return math.NaN()
	}
	if p <= 0 {
		return math.Inf(-1)
	}
	if p >= 1 {
		return math.Inf(1)
	}
	// Widen a symmetric bracket outward until CDF(lo) <= p <= CDF(hi).
	lo, hi := -1.0, 1.0
	for i := 0; i < 100 && t.CDF(lo) > p; i++ {
		lo *= 2
	}
	for i := 0; i < 100 && t.CDF(hi) < p; i++ {
		hi *= 2
	}
	// Bisection: a fixed iteration count keeps the result deterministic and
	// converges the bracket to well below floating-point resolution.
	for i := 0; i < 100; i++ {
		mid := (lo + hi) / 2
		if t.CDF(mid) < p {
			lo = mid
		} else {
			hi = mid
		}
	}
	return (lo + hi) / 2
}

// statsClamp01 clamps x to the closed unit interval [0, 1].
func statsClamp01(x float64) float64 {
	if x < 0 {
		return 0
	}
	if x > 1 {
		return 1
	}
	return x
}

// MeanConfidenceInterval returns the two-sided confidence interval [lo, hi] for
// the population mean of xs at the given confidence level (for example 0.95 for
// a 95% interval). The interval is centred on the sample mean with half-width
// StudentT{n-1}.Quantile((1+level)/2) * StdDev(xs)/sqrt(n), using the Student's
// t distribution with n-1 degrees of freedom. It returns (NaN, NaN) when xs has
// fewer than two elements or level is not in the open interval (0, 1).
func MeanConfidenceInterval(xs []float64, level float64) (lo, hi float64) {
	n := len(xs)
	if n < 2 || level <= 0 || level >= 1 {
		return math.NaN(), math.NaN()
	}
	m := Mean(xs)
	se := StdDev(xs) / math.Sqrt(float64(n))
	crit := StudentT{Nu: float64(n - 1)}.Quantile((1 + level) / 2)
	h := crit * se
	return m - h, m + h
}

// MeanConfidenceIntervalZ returns the two-sided confidence interval [lo, hi] for
// the population mean of xs at the given confidence level using the normal (z)
// approximation with a known population standard deviation sigma. The bounds are
// the (1-level)/2 and (1+level)/2 quantiles of the sampling distribution
// Normal{Mean(xs), sigma/sqrt(n)}, evaluated with normQuantile via
// Normal.Quantile. It returns (NaN, NaN) when xs is empty, sigma is not
// positive, or level is not in the open interval (0, 1).
func MeanConfidenceIntervalZ(xs []float64, sigma, level float64) (lo, hi float64) {
	n := len(xs)
	if n < 1 || sigma <= 0 || level <= 0 || level >= 1 {
		return math.NaN(), math.NaN()
	}
	samp := Normal{Mu: Mean(xs), Sigma: sigma / math.Sqrt(float64(n))}
	alpha := (1 - level) / 2
	return samp.Quantile(alpha), samp.Quantile(1 - alpha)
}

// ProportionConfidenceInterval returns the two-sided Wilson score confidence
// interval [lo, hi] for a binomial proportion given successes out of n trials at
// the given confidence level. The Wilson interval keeps its coverage near the
// boundaries far better than the plain Wald interval; the bounds are clamped to
// [0, 1]. It returns (NaN, NaN) when n <= 0, successes is outside [0, n], or
// level is not in the open interval (0, 1).
func ProportionConfidenceInterval(successes, n int, level float64) (lo, hi float64) {
	if n <= 0 || successes < 0 || successes > n || level <= 0 || level >= 1 {
		return math.NaN(), math.NaN()
	}
	z := normQuantile((1 + level) / 2)
	z2 := z * z
	nf := float64(n)
	phat := float64(successes) / nf
	denom := 1 + z2/nf
	center := (phat + z2/(2*nf)) / denom
	half := z * math.Sqrt(phat*(1-phat)/nf+z2/(4*nf*nf)) / denom
	return statsClamp01(center - half), statsClamp01(center + half)
}

// BootstrapResult holds the outcome of a bootstrap resampling analysis:
// Estimate is the statistic evaluated on the original sample, StdErr is the
// standard error estimated as the standard deviation of the resampled
// statistics, and Lo and Hi are the endpoints of the percentile confidence
// interval.
type BootstrapResult struct {
	Estimate float64 // Estimate is the statistic on the original sample.
	StdErr   float64 // StdErr is the bootstrap standard error.
	Lo       float64 // Lo is the lower confidence-interval bound.
	Hi       float64 // Hi is the upper confidence-interval bound.
}

// Bootstrap performs a nonparametric bootstrap of the statistic stat over xs. It
// draws iters resamples with replacement via Resample, evaluates stat on each,
// and reports Estimate = stat(xs), StdErr = StdDev of the replicates, and a
// percentile confidence interval whose Lo and Hi are the (1-level)/2 and
// (1+level)/2 quantiles of the replicates (via Quantile). Randomness is drawn
// exclusively from r, so a seeded *rand.Rand yields fully deterministic results.
// The replicate slice is allocated once up front and filled by index rather than
// grown, so the hot loop performs no per-iteration reslicing. It returns an
// all-NaN result when xs is empty, stat is nil, iters <= 0, or level is not in
// the open interval (0, 1).
func Bootstrap(xs []float64, stat func([]float64) float64, iters int, level float64, r *rand.Rand) BootstrapResult {
	if len(xs) == 0 || stat == nil || iters <= 0 || level <= 0 || level >= 1 {
		nan := math.NaN()
		return BootstrapResult{Estimate: nan, StdErr: nan, Lo: nan, Hi: nan}
	}
	reps := make([]float64, iters)
	for i := range reps {
		reps[i] = stat(Resample(xs, r))
	}
	alpha := (1 - level) / 2
	return BootstrapResult{
		Estimate: stat(xs),
		StdErr:   StdDev(reps),
		Lo:       Quantile(reps, alpha),
		Hi:       Quantile(reps, 1-alpha),
	}
}

// BootstrapMean is a convenience wrapper around Bootstrap that resamples the
// sample mean of xs, returning its bootstrap estimate, standard error, and
// percentile confidence interval at the given level. Results are deterministic
// for a seeded r.
func BootstrapMean(xs []float64, iters int, level float64, r *rand.Rand) BootstrapResult {
	return Bootstrap(xs, Mean, iters, level, r)
}
