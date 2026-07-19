package reliability

import "math"

// LognormalPDF returns the probability density of the lognormal lifetime
// distribution whose logarithm is normal with mean mu and standard deviation
// sigma>0, at time t>0.
func LognormalPDF(t, mu, sigma float64) float64 {
	if sigma <= 0 || t < 0 {
		return math.NaN()
	}
	if t == 0 {
		return 0
	}
	z := (math.Log(t) - mu) / sigma
	return normalPDF(z) / (t * sigma)
}

// LognormalCDF returns the failure probability F(t)=Φ((ln t-μ)/σ).
func LognormalCDF(t, mu, sigma float64) float64 {
	if sigma <= 0 || t < 0 {
		return math.NaN()
	}
	if t == 0 {
		return 0
	}
	return normalCDF((math.Log(t) - mu) / sigma)
}

// LognormalReliability returns the survival probability R(t)=1-Φ((ln t-μ)/σ).
func LognormalReliability(t, mu, sigma float64) float64 {
	if sigma <= 0 || t < 0 {
		return math.NaN()
	}
	if t == 0 {
		return 1
	}
	return 1 - normalCDF((math.Log(t)-mu)/sigma)
}

// LognormalHazard returns the hazard rate f(t)/R(t) of the lognormal
// distribution, which rises from 0, peaks, and then decreases.
func LognormalHazard(t, mu, sigma float64) float64 {
	if sigma <= 0 || t < 0 {
		return math.NaN()
	}
	r := LognormalReliability(t, mu, sigma)
	if r <= 0 {
		return math.NaN()
	}
	return LognormalPDF(t, mu, sigma) / r
}

// LognormalCumulativeHazard returns the cumulative hazard H(t)=-ln R(t).
func LognormalCumulativeHazard(t, mu, sigma float64) float64 {
	r := LognormalReliability(t, mu, sigma)
	if math.IsNaN(r) || r <= 0 {
		return math.NaN()
	}
	return -math.Log(r)
}

// LognormalQuantile returns the time at which F(t)=p:
// t=exp(μ+σΦ^{-1}(p)).
func LognormalQuantile(p, mu, sigma float64) float64 {
	if sigma <= 0 || p < 0 || p > 1 {
		return math.NaN()
	}
	return math.Exp(mu + sigma*normalQuantile(p))
}

// LognormalMean returns the mean lifetime exp(μ+σ²/2).
func LognormalMean(mu, sigma float64) float64 {
	if sigma <= 0 {
		return math.NaN()
	}
	return math.Exp(mu + 0.5*sigma*sigma)
}

// LognormalMTTF returns the mean time to failure, identical to the mean.
func LognormalMTTF(mu, sigma float64) float64 {
	return LognormalMean(mu, sigma)
}

// LognormalVariance returns the variance (e^{σ²}-1)e^{2μ+σ²}.
func LognormalVariance(mu, sigma float64) float64 {
	if sigma <= 0 {
		return math.NaN()
	}
	s2 := sigma * sigma
	return (math.Exp(s2) - 1) * math.Exp(2*mu+s2)
}

// LognormalStdDev returns the standard deviation of the lognormal distribution.
func LognormalStdDev(mu, sigma float64) float64 {
	v := LognormalVariance(mu, sigma)
	if math.IsNaN(v) {
		return math.NaN()
	}
	return math.Sqrt(v)
}

// LognormalMedian returns the median lifetime e^{μ}.
func LognormalMedian(mu, sigma float64) float64 {
	if sigma <= 0 {
		return math.NaN()
	}
	return math.Exp(mu)
}

// LognormalMode returns the mode of the lognormal distribution, e^{μ-σ²}.
func LognormalMode(mu, sigma float64) float64 {
	if sigma <= 0 {
		return math.NaN()
	}
	return math.Exp(mu - sigma*sigma)
}

// LognormalMeanResidualLife returns the mean remaining life given survival to
// age t, computed by numerically integrating the reliability from t to
// infinity.
func LognormalMeanResidualLife(t, mu, sigma float64) float64 {
	if sigma <= 0 || t < 0 {
		return math.NaN()
	}
	rt := LognormalReliability(t, mu, sigma)
	if rt <= 0 {
		return math.NaN()
	}
	step := math.Exp(mu) * math.Max(1, sigma) / 4
	area := integrateToInfinity(func(u float64) float64 {
		return LognormalReliability(u, mu, sigma)
	}, t, step)
	return area / rt
}
