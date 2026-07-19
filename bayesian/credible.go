package bayesian

import (
	"math"
)

// Interval is a closed real interval [Lower, Upper].
type Interval struct {
	Lower float64
	Upper float64
}

// Width returns Upper − Lower.
func (iv Interval) Width() float64 { return iv.Upper - iv.Lower }

// Contains reports whether x lies within the closed interval.
func (iv Interval) Contains(x float64) bool { return x >= iv.Lower && x <= iv.Upper }

// quantiler is any distribution exposing an inverse CDF.
type quantiler interface {
	Quantile(p float64) float64
}

// EqualTailedInterval returns the equal-tailed credible interval of mass level
// (for example 0.95) for any distribution with a Quantile method: it cuts
// (1−level)/2 from each tail.
func EqualTailedInterval(d quantiler, level float64) Interval {
	a := (1 - level) / 2
	return Interval{Lower: d.Quantile(a), Upper: d.Quantile(1 - a)}
}

// BetaCredibleInterval returns the equal-tailed credible interval of the given
// mass for a Beta distribution.
func BetaCredibleInterval(d Beta, level float64) Interval {
	return EqualTailedInterval(d, level)
}

// GammaCredibleInterval returns the equal-tailed credible interval of the given
// mass for a Gamma distribution.
func GammaCredibleInterval(d Gamma, level float64) Interval {
	return EqualTailedInterval(d, level)
}

// NormalCredibleInterval returns the equal-tailed credible interval of the
// given mass for a Normal distribution.
func NormalCredibleInterval(d Normal, level float64) Interval {
	return EqualTailedInterval(d, level)
}

// InverseGammaCredibleInterval returns the equal-tailed credible interval of
// the given mass for an Inverse-Gamma distribution.
func InverseGammaCredibleInterval(d InverseGamma, level float64) Interval {
	return EqualTailedInterval(d, level)
}

// StudentTCredibleInterval returns the equal-tailed credible interval of the
// given mass for a Student-t distribution.
func StudentTCredibleInterval(d StudentT, level float64) Interval {
	return EqualTailedInterval(d, level)
}

// densityQuantiler couples a density with a quantile function for HDI search.
type densityQuantiler interface {
	PDF(x float64) float64
	Quantile(p float64) float64
}

// HighestDensityInterval returns the highest-density credible interval of the
// given mass for a unimodal continuous distribution exposing PDF and Quantile
// methods. It searches over the lower-tail probability that minimizes the
// interval width. For symmetric distributions the result matches the
// equal-tailed interval.
func HighestDensityInterval(d densityQuantiler, level float64) Interval {
	if level <= 0 || level >= 1 {
		return Interval{Lower: d.Quantile(0), Upper: d.Quantile(1)}
	}
	best := Interval{Lower: d.Quantile((1 - level) / 2), Upper: d.Quantile(1 - (1-level)/2)}
	bestW := best.Width()
	// Golden-section-like scan of the lower tail mass in [0, 1-level].
	const steps = 400
	maxLow := 1 - level
	for i := 0; i <= steps; i++ {
		lowP := maxLow * float64(i) / float64(steps)
		lo := d.Quantile(lowP)
		hi := d.Quantile(lowP + level)
		w := hi - lo
		if w < bestW {
			bestW = w
			best = Interval{Lower: lo, Upper: hi}
		}
	}
	return best
}

// BetaHDI returns the highest-density credible interval of the given mass for a
// Beta distribution.
func BetaHDI(d Beta, level float64) Interval {
	return HighestDensityInterval(d, level)
}

// GammaHDI returns the highest-density credible interval of the given mass for
// a Gamma distribution.
func GammaHDI(d Gamma, level float64) Interval {
	return HighestDensityInterval(d, level)
}

// ProbabilityGreater returns P(X > threshold) for a distribution exposing a
// CDF method.
func ProbabilityGreater(d interface{ CDF(float64) float64 }, threshold float64) float64 {
	return 1 - d.CDF(threshold)
}

// ProbabilityBetween returns P(lo < X ≤ hi) for a distribution exposing a CDF.
func ProbabilityBetween(d interface{ CDF(float64) float64 }, lo, hi float64) float64 {
	return math.Max(0, d.CDF(hi)-d.CDF(lo))
}

// ProbabilityBetaGreater returns the posterior probability that a Beta random
// variable exceeds threshold, P(θ > t).
func ProbabilityBetaGreater(d Beta, threshold float64) float64 {
	return 1 - d.CDF(threshold)
}

// ProbabilityBetaExceedsBeta returns the probability P(X > Y) that a draw from
// Beta a exceeds an independent draw from Beta b, computed by numerical
// integration over the density of X times the CDF-complement of Y. This is the
// core quantity in Bayesian A/B testing.
func ProbabilityBetaExceedsBeta(a, b Beta) float64 {
	// P(X>Y) = ∫₀¹ f_a(x) F_b(x) dx via Simpson's rule.
	const n = 2000
	h := 1.0 / n
	var sum float64
	for i := 0; i <= n; i++ {
		x := float64(i) * h
		w := 4.0
		if i == 0 || i == n {
			w = 1
		} else if i%2 == 0 {
			w = 2
		}
		sum += w * a.PDF(x) * b.CDF(x)
	}
	return sum * h / 3
}
