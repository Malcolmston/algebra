package reliability

import "math"

// ExponentialPDF returns the probability density f(t)=λe^{-λt} of the
// exponential lifetime distribution with failure rate lambda>0 at time t>=0.
func ExponentialPDF(t, lambda float64) float64 {
	if lambda <= 0 || t < 0 {
		return math.NaN()
	}
	return lambda * math.Exp(-lambda*t)
}

// ExponentialCDF returns the failure probability F(t)=1-e^{-λt} for the
// exponential distribution with failure rate lambda>0 at time t>=0.
func ExponentialCDF(t, lambda float64) float64 {
	if lambda <= 0 || t < 0 {
		return math.NaN()
	}
	return -math.Expm1(-lambda * t)
}

// ExponentialReliability returns the survival probability R(t)=e^{-λt} for the
// exponential distribution with failure rate lambda>0 at time t>=0.
func ExponentialReliability(t, lambda float64) float64 {
	if lambda <= 0 || t < 0 {
		return math.NaN()
	}
	return math.Exp(-lambda * t)
}

// ExponentialHazard returns the (constant) hazard rate h(t)=λ of the
// exponential distribution.
func ExponentialHazard(t, lambda float64) float64 {
	if lambda <= 0 || t < 0 {
		return math.NaN()
	}
	return lambda
}

// ExponentialCumulativeHazard returns the cumulative hazard H(t)=λt of the
// exponential distribution.
func ExponentialCumulativeHazard(t, lambda float64) float64 {
	if lambda <= 0 || t < 0 {
		return math.NaN()
	}
	return lambda * t
}

// ExponentialQuantile returns the time t at which F(t)=p, i.e.
// t=-ln(1-p)/λ, for p in [0,1).
func ExponentialQuantile(p, lambda float64) float64 {
	if lambda <= 0 || p < 0 || p >= 1 {
		return math.NaN()
	}
	return -math.Log1p(-p) / lambda
}

// ExponentialMean returns the mean lifetime 1/λ.
func ExponentialMean(lambda float64) float64 {
	if lambda <= 0 {
		return math.NaN()
	}
	return 1 / lambda
}

// ExponentialMTTF returns the mean time to failure 1/λ, identical to the mean.
func ExponentialMTTF(lambda float64) float64 {
	return ExponentialMean(lambda)
}

// ExponentialVariance returns the variance 1/λ² of the exponential
// distribution.
func ExponentialVariance(lambda float64) float64 {
	if lambda <= 0 {
		return math.NaN()
	}
	return 1 / (lambda * lambda)
}

// ExponentialStdDev returns the standard deviation 1/λ of the exponential
// distribution.
func ExponentialStdDev(lambda float64) float64 {
	if lambda <= 0 {
		return math.NaN()
	}
	return 1 / lambda
}

// ExponentialMedian returns the median lifetime ln(2)/λ.
func ExponentialMedian(lambda float64) float64 {
	if lambda <= 0 {
		return math.NaN()
	}
	return math.Ln2 / lambda
}

// ExponentialConditionalReliability returns the probability that a component
// surviving to age t survives an additional time x. By the memoryless
// property this equals e^{-λx} independent of t.
func ExponentialConditionalReliability(t, x, lambda float64) float64 {
	if lambda <= 0 || t < 0 || x < 0 {
		return math.NaN()
	}
	return math.Exp(-lambda * x)
}

// ExponentialMeanResidualLife returns the mean remaining life given survival to
// age t, which is the constant 1/λ for the exponential distribution.
func ExponentialMeanResidualLife(t, lambda float64) float64 {
	if lambda <= 0 || t < 0 {
		return math.NaN()
	}
	return 1 / lambda
}

// ExponentialRateFromReliability recovers the failure rate λ from a reliability
// observation R at time t>0: λ=-ln(R)/t.
func ExponentialRateFromReliability(t, r float64) float64 {
	if t <= 0 || r <= 0 || r > 1 {
		return math.NaN()
	}
	return -math.Log(r) / t
}

// ExponentialRateFromMTTF recovers the failure rate λ=1/MTTF.
func ExponentialRateFromMTTF(mttf float64) float64 {
	if mttf <= 0 {
		return math.NaN()
	}
	return 1 / mttf
}
