package reliability

import "math"

// GammaPDF returns the probability density of the gamma lifetime distribution
// with shape k>0 and scale theta>0 at time t>=0:
// f(t)=t^{k-1}e^{-t/θ}/(θ^k Γ(k)).
func GammaPDF(t, k, theta float64) float64 {
	if k <= 0 || theta <= 0 || t < 0 {
		return math.NaN()
	}
	if t == 0 {
		if k < 1 {
			return math.Inf(1)
		}
		if k == 1 {
			return 1 / theta
		}
		return 0
	}
	return math.Exp((k-1)*math.Log(t) - t/theta - k*math.Log(theta) - lgamma(k))
}

// GammaCDF returns the failure probability F(t)=P(k,t/θ), the regularized
// lower incomplete gamma function.
func GammaCDF(t, k, theta float64) float64 {
	if k <= 0 || theta <= 0 || t < 0 {
		return math.NaN()
	}
	return regGammaP(k, t/theta)
}

// GammaReliability returns the survival probability R(t)=Q(k,t/θ), the
// regularized upper incomplete gamma function.
func GammaReliability(t, k, theta float64) float64 {
	if k <= 0 || theta <= 0 || t < 0 {
		return math.NaN()
	}
	return regGammaQ(k, t/theta)
}

// GammaHazard returns the hazard rate f(t)/R(t) of the gamma distribution.
func GammaHazard(t, k, theta float64) float64 {
	if k <= 0 || theta <= 0 || t < 0 {
		return math.NaN()
	}
	r := GammaReliability(t, k, theta)
	if r <= 0 {
		return math.NaN()
	}
	return GammaPDF(t, k, theta) / r
}

// GammaCumulativeHazard returns the cumulative hazard H(t)=-ln R(t).
func GammaCumulativeHazard(t, k, theta float64) float64 {
	r := GammaReliability(t, k, theta)
	if math.IsNaN(r) || r <= 0 {
		return math.NaN()
	}
	return -math.Log(r)
}

// GammaQuantile returns the time at which F(t)=p, computed by inverting the
// regularized incomplete gamma function.
func GammaQuantile(p, k, theta float64) float64 {
	if k <= 0 || theta <= 0 || p < 0 || p > 1 {
		return math.NaN()
	}
	return theta * invRegGammaP(k, p)
}

// GammaMean returns the mean lifetime kθ.
func GammaMean(k, theta float64) float64 {
	if k <= 0 || theta <= 0 {
		return math.NaN()
	}
	return k * theta
}

// GammaMTTF returns the mean time to failure, identical to the mean.
func GammaMTTF(k, theta float64) float64 {
	return GammaMean(k, theta)
}

// GammaVariance returns the variance kθ².
func GammaVariance(k, theta float64) float64 {
	if k <= 0 || theta <= 0 {
		return math.NaN()
	}
	return k * theta * theta
}

// GammaStdDev returns the standard deviation of the gamma distribution.
func GammaStdDev(k, theta float64) float64 {
	v := GammaVariance(k, theta)
	if math.IsNaN(v) {
		return math.NaN()
	}
	return math.Sqrt(v)
}

// GammaMode returns the mode of the gamma distribution. For k>=1 it is
// (k-1)θ; for k<1 the density is unbounded at 0 and the mode is 0.
func GammaMode(k, theta float64) float64 {
	if k <= 0 || theta <= 0 {
		return math.NaN()
	}
	if k < 1 {
		return 0
	}
	return (k - 1) * theta
}

// GammaMeanResidualLife returns the mean remaining life given survival to age
// t, computed by numerically integrating the reliability from t to infinity.
func GammaMeanResidualLife(t, k, theta float64) float64 {
	if k <= 0 || theta <= 0 || t < 0 {
		return math.NaN()
	}
	rt := GammaReliability(t, k, theta)
	if rt <= 0 {
		return math.NaN()
	}
	step := k * theta / 4
	if step <= 0 {
		step = theta
	}
	area := integrateToInfinity(func(u float64) float64 {
		return GammaReliability(u, k, theta)
	}, t, step)
	return area / rt
}
