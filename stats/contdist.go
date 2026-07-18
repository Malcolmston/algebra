package stats

import "math"

// This file adds four continuous probability distributions to the package:
// Beta, FDist, LogNormal and Weibull. They follow the same value-struct style
// as distributions.go, exposing PDF, CDF, Quantile, Mean and Variance methods.
//
// Performance/numerical technique: every density is evaluated in log space and
// exponentiated with a single math.Exp call. Accumulating the normalizing
// constant and the variable part as a sum of logarithms (via gammaLn) avoids
// the overflow and underflow that a direct product of large gamma factors and
// small powers would suffer, and replaces several transcendental evaluations
// with one. The Beta and FDist quantiles, which have no elementary closed
// form, are found by monotone bisection on the struct's own CDF; because each
// CDF is strictly increasing on its support, bisection converges to full double
// precision in a fixed number of iterations and is fully deterministic.

// statsBisectCDFBeta returns the value x in [0, 1] at which Beta b's CDF equals
// p, for p in (0, 1). It exploits the strict monotonicity of the CDF and uses a
// fixed iteration count so the result is deterministic.
func statsBisectCDFBeta(b Beta, p float64) float64 {
	lo, hi := 0.0, 1.0
	for i := 0; i < 200; i++ {
		mid := (lo + hi) / 2
		if b.CDF(mid) < p {
			lo = mid
		} else {
			hi = mid
		}
	}
	return (lo + hi) / 2
}

// statsBisectCDFFDist returns the value x >= 0 at which FDist f's CDF equals p,
// for p in (0, 1). The upper bracket is doubled until it straddles the target,
// after which a fixed number of bisection steps refines the root; the fixed
// counts make the result deterministic.
func statsBisectCDFFDist(f FDist, p float64) float64 {
	lo, hi := 0.0, 1.0
	for f.CDF(hi) < p && hi < 1e300 {
		hi *= 2
	}
	for i := 0; i < 200; i++ {
		mid := (lo + hi) / 2
		if f.CDF(mid) < p {
			lo = mid
		} else {
			hi = mid
		}
	}
	return (lo + hi) / 2
}

// Beta is a beta distribution on the interval [0, 1] with shape parameters
// Alpha (> 0) and Beta (> 0). It is the conjugate prior of the Bernoulli and
// binomial success probability and, more generally, models a bounded quantity
// whose two shape parameters trade probability mass between the endpoints.
type Beta struct {
	Alpha float64
	Beta  float64
}

// PDF returns the probability density at x. It is 0 outside [0, 1] and NaN when
// either shape parameter is not positive. The interior is evaluated from the
// logarithm lnBeta = gammaLn(Alpha)+gammaLn(Beta)-gammaLn(Alpha+Beta) plus the
// variable terms (Alpha-1)·ln x and (Beta-1)·ln(1-x), exponentiated once.
func (b Beta) PDF(x float64) float64 {
	a, bb := b.Alpha, b.Beta
	if a <= 0 || bb <= 0 {
		return math.NaN()
	}
	if x < 0 || x > 1 {
		return 0
	}
	// Endpoint handling avoids the indeterminate 0·(-Inf) that the log-space
	// form would otherwise produce at x = 0 or x = 1.
	if x == 0 {
		switch {
		case a < 1:
			return math.Inf(1)
		case a == 1:
			return math.Exp(gammaLn(a+bb) - gammaLn(a) - gammaLn(bb))
		default:
			return 0
		}
	}
	if x == 1 {
		switch {
		case bb < 1:
			return math.Inf(1)
		case bb == 1:
			return math.Exp(gammaLn(a+bb) - gammaLn(a) - gammaLn(bb))
		default:
			return 0
		}
	}
	lnBeta := gammaLn(a) + gammaLn(bb) - gammaLn(a+bb)
	lnP := (a-1)*math.Log(x) + (bb-1)*math.Log(1-x) - lnBeta
	return math.Exp(lnP)
}

// CDF returns the cumulative probability P(X <= x), the regularized incomplete
// beta function I_x(Alpha, Beta). It is 0 for x <= 0, 1 for x >= 1, and NaN when
// either shape parameter is not positive.
func (b Beta) CDF(x float64) float64 {
	if b.Alpha <= 0 || b.Beta <= 0 {
		return math.NaN()
	}
	return regularizedIncompleteBeta(b.Alpha, b.Beta, x)
}

// Quantile returns the inverse CDF: the value x such that CDF(x) = p, for p in
// [0, 1]. It is found by monotone bisection on the CDF. Quantile returns NaN if
// p is outside [0, 1] or if either shape parameter is not positive.
func (b Beta) Quantile(p float64) float64 {
	if b.Alpha <= 0 || b.Beta <= 0 || p < 0 || p > 1 || math.IsNaN(p) {
		return math.NaN()
	}
	if p == 0 {
		return 0
	}
	if p == 1 {
		return 1
	}
	return statsBisectCDFBeta(b, p)
}

// Mean returns the mean of the distribution, Alpha/(Alpha+Beta).
func (b Beta) Mean() float64 {
	if b.Alpha <= 0 || b.Beta <= 0 {
		return math.NaN()
	}
	return b.Alpha / (b.Alpha + b.Beta)
}

// Variance returns the variance of the distribution,
// Alpha·Beta / ((Alpha+Beta)²·(Alpha+Beta+1)).
func (b Beta) Variance() float64 {
	if b.Alpha <= 0 || b.Beta <= 0 {
		return math.NaN()
	}
	s := b.Alpha + b.Beta
	return b.Alpha * b.Beta / (s * s * (s + 1))
}

// FDist is Fisher's F distribution with degrees of freedom D1 and D2 (both > 0),
// the ratio of two independent chi-squared variates each divided by its degrees
// of freedom. It is the null distribution of the variance ratio in an analysis
// of variance.
type FDist struct {
	D1 float64
	D2 float64
}

// PDF returns the probability density at x. It is 0 for x < 0 and NaN when
// either degrees-of-freedom parameter is not positive. The interior is built
// from a sum of logarithms and exponentiated once.
func (f FDist) PDF(x float64) float64 {
	d1, d2 := f.D1, f.D2
	if d1 <= 0 || d2 <= 0 {
		return math.NaN()
	}
	if x < 0 {
		return 0
	}
	if x == 0 {
		switch {
		case d1 < 2:
			return math.Inf(1)
		case d1 == 2:
			return 1
		default:
			return 0
		}
	}
	lnBeta := gammaLn(d1/2) + gammaLn(d2/2) - gammaLn((d1+d2)/2)
	lnNum := (d1/2)*math.Log(d1) + (d2/2)*math.Log(d2) + (d1/2-1)*math.Log(x)
	lnDen := ((d1 + d2) / 2) * math.Log(d1*x+d2)
	return math.Exp(lnNum - lnDen - lnBeta)
}

// CDF returns the cumulative probability P(X <= x), the regularized incomplete
// beta function I_{D1·x/(D1·x+D2)}(D1/2, D2/2). It is 0 for x < 0 and NaN when
// either degrees-of-freedom parameter is not positive.
func (f FDist) CDF(x float64) float64 {
	d1, d2 := f.D1, f.D2
	if d1 <= 0 || d2 <= 0 {
		return math.NaN()
	}
	if x < 0 {
		return 0
	}
	return regularizedIncompleteBeta(d1/2, d2/2, d1*x/(d1*x+d2))
}

// Quantile returns the inverse CDF: the value x such that CDF(x) = p, for p in
// [0, 1]. It is found by monotone bisection on the CDF after bracketing the
// root. Quantile returns NaN if p is outside [0, 1] or if either
// degrees-of-freedom parameter is not positive.
func (f FDist) Quantile(p float64) float64 {
	if f.D1 <= 0 || f.D2 <= 0 || p < 0 || p > 1 || math.IsNaN(p) {
		return math.NaN()
	}
	if p == 0 {
		return 0
	}
	if p == 1 {
		return math.Inf(1)
	}
	return statsBisectCDFFDist(f, p)
}

// Mean returns the mean of the distribution, D2/(D2-2) for D2 > 2, and NaN
// otherwise.
func (f FDist) Mean() float64 {
	if f.D1 <= 0 || f.D2 <= 0 {
		return math.NaN()
	}
	if f.D2 > 2 {
		return f.D2 / (f.D2 - 2)
	}
	return math.NaN()
}

// Variance returns the variance of the distribution,
// 2·D2²·(D1+D2-2) / (D1·(D2-2)²·(D2-4)) for D2 > 4. It is +Inf for 2 < D2 <= 4
// and NaN otherwise (the variance does not exist).
func (f FDist) Variance() float64 {
	if f.D1 <= 0 || f.D2 <= 0 {
		return math.NaN()
	}
	d1, d2 := f.D1, f.D2
	switch {
	case d2 > 4:
		return 2 * d2 * d2 * (d1 + d2 - 2) / (d1 * (d2 - 2) * (d2 - 2) * (d2 - 4))
	case d2 > 2:
		return math.Inf(1)
	default:
		return math.NaN()
	}
}

// LogNormal is a log-normal distribution: a positive variable whose natural
// logarithm is normally distributed with mean Mu and standard deviation Sigma
// (Sigma > 0). Mu and Sigma are the parameters of that underlying normal, not
// the mean and standard deviation of the log-normal variable itself.
type LogNormal struct {
	Mu    float64
	Sigma float64
}

// PDF returns the probability density at x, computed as the underlying normal
// density at ln x divided by x. It is 0 for x <= 0 and NaN when Sigma is not
// positive.
func (l LogNormal) PDF(x float64) float64 {
	if l.Sigma <= 0 {
		return math.NaN()
	}
	if x <= 0 {
		return 0
	}
	return Normal{Mu: l.Mu, Sigma: l.Sigma}.PDF(math.Log(x)) / x
}

// CDF returns the cumulative probability P(X <= x), the underlying normal CDF
// evaluated at ln x. It is 0 for x <= 0 and NaN when Sigma is not positive.
func (l LogNormal) CDF(x float64) float64 {
	if l.Sigma <= 0 {
		return math.NaN()
	}
	if x <= 0 {
		return 0
	}
	return Normal{Mu: l.Mu, Sigma: l.Sigma}.CDF(math.Log(x))
}

// Quantile returns the inverse CDF for p in [0, 1], exp(Mu + Sigma·Φ⁻¹(p)) where
// Φ⁻¹ is the standard normal quantile. It returns NaN if p is outside [0, 1] or
// if Sigma is not positive.
func (l LogNormal) Quantile(p float64) float64 {
	if l.Sigma <= 0 || p < 0 || p > 1 || math.IsNaN(p) {
		return math.NaN()
	}
	return math.Exp(l.Mu + l.Sigma*normQuantile(p))
}

// Mean returns the mean of the distribution, exp(Mu + Sigma²/2).
func (l LogNormal) Mean() float64 {
	if l.Sigma <= 0 {
		return math.NaN()
	}
	return math.Exp(l.Mu + l.Sigma*l.Sigma/2)
}

// Variance returns the variance of the distribution,
// (exp(Sigma²)-1)·exp(2·Mu + Sigma²).
func (l LogNormal) Variance() float64 {
	if l.Sigma <= 0 {
		return math.NaN()
	}
	s2 := l.Sigma * l.Sigma
	return (math.Exp(s2) - 1) * math.Exp(2*l.Mu+s2)
}

// Weibull is a Weibull distribution with shape parameter Shape (k, > 0) and
// scale parameter Scale (λ, > 0). With Shape = 1 it reduces to an exponential
// distribution with rate 1/Scale; it is widely used in reliability and
// survival analysis.
type Weibull struct {
	Shape float64
	Scale float64
}

// PDF returns the probability density at x,
// (k/λ)·(x/λ)^{k-1}·exp(-(x/λ)^k). It is 0 for x < 0 and NaN when either
// parameter is not positive.
func (w Weibull) PDF(x float64) float64 {
	k, lam := w.Shape, w.Scale
	if k <= 0 || lam <= 0 {
		return math.NaN()
	}
	if x < 0 {
		return 0
	}
	if x == 0 {
		switch {
		case k < 1:
			return math.Inf(1)
		case k == 1:
			return 1 / lam
		default:
			return 0
		}
	}
	z := x / lam
	return (k / lam) * math.Pow(z, k-1) * math.Exp(-math.Pow(z, k))
}

// CDF returns the cumulative probability P(X <= x), 1 - exp(-(x/λ)^k). It is 0
// for x < 0 and NaN when either parameter is not positive.
func (w Weibull) CDF(x float64) float64 {
	if w.Shape <= 0 || w.Scale <= 0 {
		return math.NaN()
	}
	if x < 0 {
		return 0
	}
	return 1 - math.Exp(-math.Pow(x/w.Scale, w.Shape))
}

// Quantile returns the inverse CDF for p in [0, 1], λ·(-ln(1-p))^{1/k}. It
// returns NaN if p is outside [0, 1] or if either parameter is not positive.
func (w Weibull) Quantile(p float64) float64 {
	if w.Shape <= 0 || w.Scale <= 0 || p < 0 || p > 1 || math.IsNaN(p) {
		return math.NaN()
	}
	if p == 1 {
		return math.Inf(1)
	}
	return w.Scale * math.Pow(-math.Log(1-p), 1/w.Shape)
}

// Mean returns the mean of the distribution, λ·Γ(1 + 1/k).
func (w Weibull) Mean() float64 {
	if w.Shape <= 0 || w.Scale <= 0 {
		return math.NaN()
	}
	return w.Scale * math.Gamma(1+1/w.Shape)
}

// Variance returns the variance of the distribution,
// λ²·(Γ(1 + 2/k) - Γ(1 + 1/k)²).
func (w Weibull) Variance() float64 {
	if w.Shape <= 0 || w.Scale <= 0 {
		return math.NaN()
	}
	g1 := math.Gamma(1 + 1/w.Shape)
	g2 := math.Gamma(1 + 2/w.Shape)
	return w.Scale * w.Scale * (g2 - g1*g1)
}
