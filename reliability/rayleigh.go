package reliability

import "math"

// RayleighPDF returns the probability density of the Rayleigh lifetime
// distribution with scale sigma>0 at time t>=0:
// f(t)=(t/σ²)e^{-t²/(2σ²)}. The Rayleigh model is the Weibull model with
// shape k=2 and a linearly increasing hazard rate.
func RayleighPDF(t, sigma float64) float64 {
	if sigma <= 0 || t < 0 {
		return math.NaN()
	}
	s2 := sigma * sigma
	return (t / s2) * math.Exp(-t*t/(2*s2))
}

// RayleighCDF returns the failure probability F(t)=1-e^{-t²/(2σ²)}.
func RayleighCDF(t, sigma float64) float64 {
	if sigma <= 0 || t < 0 {
		return math.NaN()
	}
	return -math.Expm1(-t * t / (2 * sigma * sigma))
}

// RayleighReliability returns the survival probability R(t)=e^{-t²/(2σ²)}.
func RayleighReliability(t, sigma float64) float64 {
	if sigma <= 0 || t < 0 {
		return math.NaN()
	}
	return math.Exp(-t * t / (2 * sigma * sigma))
}

// RayleighHazard returns the linearly increasing hazard rate h(t)=t/σ².
func RayleighHazard(t, sigma float64) float64 {
	if sigma <= 0 || t < 0 {
		return math.NaN()
	}
	return t / (sigma * sigma)
}

// RayleighCumulativeHazard returns the cumulative hazard H(t)=t²/(2σ²).
func RayleighCumulativeHazard(t, sigma float64) float64 {
	if sigma <= 0 || t < 0 {
		return math.NaN()
	}
	return t * t / (2 * sigma * sigma)
}

// RayleighQuantile returns the time at which F(t)=p:
// t=σ√(-2ln(1-p)).
func RayleighQuantile(p, sigma float64) float64 {
	if sigma <= 0 || p < 0 || p >= 1 {
		return math.NaN()
	}
	return sigma * math.Sqrt(-2*math.Log1p(-p))
}

// RayleighMean returns the mean lifetime σ√(π/2).
func RayleighMean(sigma float64) float64 {
	if sigma <= 0 {
		return math.NaN()
	}
	return sigma * math.Sqrt(math.Pi/2)
}

// RayleighMTTF returns the mean time to failure, identical to the mean.
func RayleighMTTF(sigma float64) float64 {
	return RayleighMean(sigma)
}

// RayleighVariance returns the variance (2-π/2)σ².
func RayleighVariance(sigma float64) float64 {
	if sigma <= 0 {
		return math.NaN()
	}
	return (2 - math.Pi/2) * sigma * sigma
}

// RayleighStdDev returns the standard deviation of the Rayleigh distribution.
func RayleighStdDev(sigma float64) float64 {
	v := RayleighVariance(sigma)
	if math.IsNaN(v) {
		return math.NaN()
	}
	return math.Sqrt(v)
}

// RayleighMedian returns the median lifetime σ√(2ln2).
func RayleighMedian(sigma float64) float64 {
	if sigma <= 0 {
		return math.NaN()
	}
	return sigma * math.Sqrt(2*math.Ln2)
}

// RayleighMode returns the mode of the Rayleigh distribution, equal to σ.
func RayleighMode(sigma float64) float64 {
	if sigma <= 0 {
		return math.NaN()
	}
	return sigma
}
