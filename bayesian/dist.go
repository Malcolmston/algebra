package bayesian

import (
	"errors"
	"math"
)

// ErrParam is returned by constructors and updates when a distribution
// parameter is outside its valid domain.
var ErrParam = errors.New("bayesian: invalid distribution parameter")

// ------------------------------------------------------------------
// Beta distribution
// ------------------------------------------------------------------

// Beta is a Beta distribution with shape parameters Alpha, Beta > 0. It is the
// conjugate prior for the success probability of Bernoulli and Binomial
// sampling and for the rate of a geometric/negative-binomial model.
type Beta struct {
	Alpha float64
	Beta  float64
}

// NewBeta constructs a Beta distribution, returning ErrParam if either shape is
// non-positive.
func NewBeta(alpha, beta float64) (Beta, error) {
	if alpha <= 0 || beta <= 0 {
		return Beta{}, ErrParam
	}
	return Beta{Alpha: alpha, Beta: beta}, nil
}

// Mean returns the expected value α/(α+β).
func (d Beta) Mean() float64 { return d.Alpha / (d.Alpha + d.Beta) }

// Variance returns αβ/((α+β)²(α+β+1)).
func (d Beta) Variance() float64 {
	s := d.Alpha + d.Beta
	return d.Alpha * d.Beta / (s * s * (s + 1))
}

// StdDev returns the standard deviation, the square root of the variance.
func (d Beta) StdDev() float64 { return math.Sqrt(d.Variance()) }

// Mode returns the mode (α−1)/(α+β−2). It is defined only when α,β > 1; for
// other shapes it returns NaN because the mode lies at a boundary or is not
// unique.
func (d Beta) Mode() float64 {
	if d.Alpha > 1 && d.Beta > 1 {
		return (d.Alpha - 1) / (d.Alpha + d.Beta - 2)
	}
	return math.NaN()
}

// Skewness returns the distribution's skewness.
func (d Beta) Skewness() float64 {
	a, b := d.Alpha, d.Beta
	return 2 * (b - a) * math.Sqrt(a+b+1) / ((a + b + 2) * math.Sqrt(a*b))
}

// PDF returns the probability density at x in [0,1].
func (d Beta) PDF(x float64) float64 {
	if x < 0 || x > 1 {
		return 0
	}
	return math.Exp(d.LogPDF(x))
}

// LogPDF returns the natural logarithm of the density at x in (0,1). It returns
// math.Inf(-1) outside the open unit interval.
func (d Beta) LogPDF(x float64) float64 {
	if x <= 0 || x >= 1 {
		return math.Inf(-1)
	}
	return (d.Alpha-1)*math.Log(x) + (d.Beta-1)*math.Log(1-x) - LogBeta(d.Alpha, d.Beta)
}

// CDF returns the cumulative distribution function P(X ≤ x).
func (d Beta) CDF(x float64) float64 {
	return RegularizedIncompleteBeta(x, d.Alpha, d.Beta)
}

// Quantile returns the p-quantile (inverse CDF) for p in [0,1].
func (d Beta) Quantile(p float64) float64 {
	return InverseRegularizedIncompleteBeta(p, d.Alpha, d.Beta)
}

// MeanLog returns E[ln X] = ψ(α) − ψ(α+β).
func (d Beta) MeanLog() float64 {
	return Digamma(d.Alpha) - Digamma(d.Alpha+d.Beta)
}

// Entropy returns the differential entropy of the distribution in nats.
func (d Beta) Entropy() float64 {
	a, b := d.Alpha, d.Beta
	return LogBeta(a, b) - (a-1)*Digamma(a) - (b-1)*Digamma(b) + (a+b-2)*Digamma(a+b)
}

// ------------------------------------------------------------------
// Gamma distribution (shape / rate parameterization)
// ------------------------------------------------------------------

// Gamma is a Gamma distribution with Shape > 0 and Rate > 0 (rate = 1/scale).
// It is the conjugate prior for the rate of a Poisson or exponential model and
// for the precision of a normal likelihood.
type Gamma struct {
	Shape float64
	Rate  float64
}

// NewGamma constructs a Gamma distribution in the shape/rate parameterization,
// returning ErrParam if either parameter is non-positive.
func NewGamma(shape, rate float64) (Gamma, error) {
	if shape <= 0 || rate <= 0 {
		return Gamma{}, ErrParam
	}
	return Gamma{Shape: shape, Rate: rate}, nil
}

// Scale returns the scale parameter 1/Rate.
func (d Gamma) Scale() float64 { return 1 / d.Rate }

// Mean returns Shape/Rate.
func (d Gamma) Mean() float64 { return d.Shape / d.Rate }

// Variance returns Shape/Rate².
func (d Gamma) Variance() float64 { return d.Shape / (d.Rate * d.Rate) }

// StdDev returns the standard deviation.
func (d Gamma) StdDev() float64 { return math.Sqrt(d.Variance()) }

// Mode returns (Shape−1)/Rate for Shape ≥ 1, and 0 otherwise.
func (d Gamma) Mode() float64 {
	if d.Shape >= 1 {
		return (d.Shape - 1) / d.Rate
	}
	return 0
}

// Skewness returns 2/√Shape.
func (d Gamma) Skewness() float64 { return 2 / math.Sqrt(d.Shape) }

// PDF returns the density at x ≥ 0.
func (d Gamma) PDF(x float64) float64 {
	if x < 0 {
		return 0
	}
	return math.Exp(d.LogPDF(x))
}

// LogPDF returns the natural logarithm of the density at x > 0.
func (d Gamma) LogPDF(x float64) float64 {
	if x <= 0 {
		return math.Inf(-1)
	}
	return d.Shape*math.Log(d.Rate) + (d.Shape-1)*math.Log(x) - d.Rate*x - LogGamma(d.Shape)
}

// CDF returns the cumulative distribution function P(X ≤ x).
func (d Gamma) CDF(x float64) float64 {
	if x <= 0 {
		return 0
	}
	return RegularizedGammaP(d.Shape, d.Rate*x)
}

// Quantile returns the p-quantile for p in [0,1].
func (d Gamma) Quantile(p float64) float64 {
	return InverseRegularizedGammaP(p, d.Shape) / d.Rate
}

// MeanLog returns E[ln X] = ψ(Shape) − ln(Rate).
func (d Gamma) MeanLog() float64 {
	return Digamma(d.Shape) - math.Log(d.Rate)
}

// Entropy returns the differential entropy in nats.
func (d Gamma) Entropy() float64 {
	k := d.Shape
	return k - math.Log(d.Rate) + LogGamma(k) + (1-k)*Digamma(k)
}

// ------------------------------------------------------------------
// Inverse-Gamma distribution
// ------------------------------------------------------------------

// InverseGamma is an Inverse-Gamma distribution with Shape > 0 and Scale > 0.
// If Y ~ Gamma(Shape, rate=Scale) then 1/Y ~ InverseGamma(Shape, Scale). It is
// the conjugate prior for the variance of a normal likelihood.
type InverseGamma struct {
	Shape float64
	Scale float64
}

// NewInverseGamma constructs an Inverse-Gamma distribution, returning ErrParam
// for non-positive parameters.
func NewInverseGamma(shape, scale float64) (InverseGamma, error) {
	if shape <= 0 || scale <= 0 {
		return InverseGamma{}, ErrParam
	}
	return InverseGamma{Shape: shape, Scale: scale}, nil
}

// Mean returns Scale/(Shape−1) for Shape > 1, else +Inf.
func (d InverseGamma) Mean() float64 {
	if d.Shape > 1 {
		return d.Scale / (d.Shape - 1)
	}
	return math.Inf(1)
}

// Variance returns the variance for Shape > 2, else +Inf.
func (d InverseGamma) Variance() float64 {
	a := d.Shape
	if a > 2 {
		return d.Scale * d.Scale / ((a - 1) * (a - 1) * (a - 2))
	}
	return math.Inf(1)
}

// Mode returns Scale/(Shape+1).
func (d InverseGamma) Mode() float64 { return d.Scale / (d.Shape + 1) }

// PDF returns the density at x > 0.
func (d InverseGamma) PDF(x float64) float64 {
	if x <= 0 {
		return 0
	}
	return math.Exp(d.LogPDF(x))
}

// LogPDF returns the natural logarithm of the density at x > 0.
func (d InverseGamma) LogPDF(x float64) float64 {
	if x <= 0 {
		return math.Inf(-1)
	}
	a, b := d.Shape, d.Scale
	return a*math.Log(b) - LogGamma(a) - (a+1)*math.Log(x) - b/x
}

// CDF returns the cumulative distribution function P(X ≤ x). For an
// Inverse-Gamma this equals Q(Shape, Scale/x), the upper regularized gamma.
func (d InverseGamma) CDF(x float64) float64 {
	if x <= 0 {
		return 0
	}
	return RegularizedGammaQ(d.Shape, d.Scale/x)
}

// Quantile returns the p-quantile for p in [0,1].
func (d InverseGamma) Quantile(p float64) float64 {
	if p <= 0 {
		return 0
	}
	if p >= 1 {
		return math.Inf(1)
	}
	// x = Scale / InvGammaP(1-p, Shape)
	return d.Scale / InverseRegularizedGammaP(1-p, d.Shape)
}

// ------------------------------------------------------------------
// Normal distribution
// ------------------------------------------------------------------

// Normal is a univariate normal (Gaussian) distribution with mean Mu and
// standard deviation Sigma > 0.
type Normal struct {
	Mu    float64
	Sigma float64
}

// NewNormal constructs a Normal distribution, returning ErrParam if Sigma ≤ 0.
func NewNormal(mu, sigma float64) (Normal, error) {
	if sigma <= 0 {
		return Normal{}, ErrParam
	}
	return Normal{Mu: mu, Sigma: sigma}, nil
}

// Mean returns Mu.
func (d Normal) Mean() float64 { return d.Mu }

// Variance returns Sigma².
func (d Normal) Variance() float64 { return d.Sigma * d.Sigma }

// StdDev returns Sigma.
func (d Normal) StdDev() float64 { return d.Sigma }

// PDF returns the density at x.
func (d Normal) PDF(x float64) float64 {
	return StdNormalPDF((x-d.Mu)/d.Sigma) / d.Sigma
}

// LogPDF returns the natural logarithm of the density at x.
func (d Normal) LogPDF(x float64) float64 {
	z := (x - d.Mu) / d.Sigma
	return -0.5*z*z - math.Log(d.Sigma) - 0.5*math.Log(2*math.Pi)
}

// CDF returns the cumulative distribution function P(X ≤ x).
func (d Normal) CDF(x float64) float64 {
	return StdNormalCDF((x - d.Mu) / d.Sigma)
}

// Quantile returns the p-quantile for p in [0,1].
func (d Normal) Quantile(p float64) float64 {
	return d.Mu + d.Sigma*StdNormalQuantile(p)
}

// Entropy returns the differential entropy in nats.
func (d Normal) Entropy() float64 {
	return 0.5 * math.Log(2*math.Pi*math.E*d.Sigma*d.Sigma)
}

// ------------------------------------------------------------------
// Student-t distribution (location-scale)
// ------------------------------------------------------------------

// StudentT is a location-scale Student's t distribution with Nu > 0 degrees of
// freedom, location Loc and scale Scale > 0. It arises as the posterior
// predictive and the marginal for the mean in a normal model with unknown
// variance.
type StudentT struct {
	Nu    float64
	Loc   float64
	Scale float64
}

// NewStudentT constructs a location-scale Student-t distribution, returning
// ErrParam for non-positive Nu or Scale.
func NewStudentT(nu, loc, scale float64) (StudentT, error) {
	if nu <= 0 || scale <= 0 {
		return StudentT{}, ErrParam
	}
	return StudentT{Nu: nu, Loc: loc, Scale: scale}, nil
}

// Mean returns Loc for Nu > 1, else NaN.
func (d StudentT) Mean() float64 {
	if d.Nu > 1 {
		return d.Loc
	}
	return math.NaN()
}

// Variance returns the variance for Nu > 2, +Inf for 1 < Nu ≤ 2, else NaN.
func (d StudentT) Variance() float64 {
	if d.Nu > 2 {
		return d.Scale * d.Scale * d.Nu / (d.Nu - 2)
	}
	if d.Nu > 1 {
		return math.Inf(1)
	}
	return math.NaN()
}

// StdDev returns the standard deviation.
func (d StudentT) StdDev() float64 { return math.Sqrt(d.Variance()) }

// PDF returns the density at x.
func (d StudentT) PDF(x float64) float64 {
	return math.Exp(d.LogPDF(x))
}

// LogPDF returns the natural logarithm of the density at x.
func (d StudentT) LogPDF(x float64) float64 {
	nu := d.Nu
	t := (x - d.Loc) / d.Scale
	lead := LogGamma((nu+1)/2) - LogGamma(nu/2) - 0.5*math.Log(nu*math.Pi) - math.Log(d.Scale)
	return lead - (nu+1)/2*math.Log(1+t*t/nu)
}

// CDF returns the cumulative distribution function P(X ≤ x).
func (d StudentT) CDF(x float64) float64 {
	t := (x - d.Loc) / d.Scale
	xx := d.Nu / (d.Nu + t*t)
	ib := 0.5 * RegularizedIncompleteBeta(xx, d.Nu/2, 0.5)
	if t > 0 {
		return 1 - ib
	}
	return ib
}

// Quantile returns the p-quantile for p in [0,1].
func (d StudentT) Quantile(p float64) float64 {
	if p <= 0 {
		return math.Inf(-1)
	}
	if p >= 1 {
		return math.Inf(1)
	}
	// Standard t quantile via inverse incomplete beta.
	var t float64
	if p < 0.5 {
		xx := InverseRegularizedIncompleteBeta(2*p, d.Nu/2, 0.5)
		t = -math.Sqrt(d.Nu * (1 - xx) / xx)
	} else if p > 0.5 {
		xx := InverseRegularizedIncompleteBeta(2*(1-p), d.Nu/2, 0.5)
		t = math.Sqrt(d.Nu * (1 - xx) / xx)
	}
	return d.Loc + d.Scale*t
}
