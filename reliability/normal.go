package reliability

import "math"

// NormalPDF returns the probability density of the normal lifetime
// distribution with mean mu and standard deviation sigma>0 at time t. Although
// the normal distribution has support over the whole real line, it is a common
// wear-out lifetime model when mu is large relative to sigma.
func NormalPDF(t, mu, sigma float64) float64 {
	if sigma <= 0 {
		return math.NaN()
	}
	return normalPDF((t-mu)/sigma) / sigma
}

// NormalCDF returns the failure probability F(t)=Φ((t-μ)/σ).
func NormalCDF(t, mu, sigma float64) float64 {
	if sigma <= 0 {
		return math.NaN()
	}
	return normalCDF((t - mu) / sigma)
}

// NormalReliability returns the survival probability R(t)=1-Φ((t-μ)/σ).
func NormalReliability(t, mu, sigma float64) float64 {
	if sigma <= 0 {
		return math.NaN()
	}
	return 1 - normalCDF((t-mu)/sigma)
}

// NormalHazard returns the hazard rate f(t)/R(t) of the normal distribution,
// which increases monotonically with t.
func NormalHazard(t, mu, sigma float64) float64 {
	if sigma <= 0 {
		return math.NaN()
	}
	r := NormalReliability(t, mu, sigma)
	if r <= 0 {
		return math.NaN()
	}
	return NormalPDF(t, mu, sigma) / r
}

// NormalCumulativeHazard returns the cumulative hazard H(t)=-ln R(t).
func NormalCumulativeHazard(t, mu, sigma float64) float64 {
	r := NormalReliability(t, mu, sigma)
	if math.IsNaN(r) || r <= 0 {
		return math.NaN()
	}
	return -math.Log(r)
}

// NormalQuantile returns the time at which F(t)=p:
// t=μ+σΦ^{-1}(p).
func NormalQuantile(p, mu, sigma float64) float64 {
	if sigma <= 0 || p < 0 || p > 1 {
		return math.NaN()
	}
	return mu + sigma*normalQuantile(p)
}

// NormalMTTF returns the mean time to failure, equal to mu for the normal
// model.
func NormalMTTF(mu, sigma float64) float64 {
	if sigma <= 0 {
		return math.NaN()
	}
	return mu
}

// StandardNormalCDF returns Φ(x), the standard-normal cumulative distribution
// function.
func StandardNormalCDF(x float64) float64 {
	return normalCDF(x)
}

// StandardNormalQuantile returns Φ^{-1}(p), the standard-normal quantile
// (probit) function.
func StandardNormalQuantile(p float64) float64 {
	return normalQuantile(p)
}
