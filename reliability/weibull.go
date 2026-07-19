package reliability

import "math"

// WeibullPDF returns the probability density of the two-parameter Weibull
// distribution with shape k>0 and scale lambda>0 at time t>=0:
// f(t)=(k/λ)(t/λ)^{k-1} e^{-(t/λ)^k}.
func WeibullPDF(t, k, lambda float64) float64 {
	if k <= 0 || lambda <= 0 || t < 0 {
		return math.NaN()
	}
	if t == 0 {
		if k < 1 {
			return math.Inf(1)
		}
		if k == 1 {
			return 1 / lambda
		}
		return 0
	}
	z := t / lambda
	return (k / lambda) * math.Pow(z, k-1) * math.Exp(-math.Pow(z, k))
}

// WeibullCDF returns the failure probability F(t)=1-e^{-(t/λ)^k}.
func WeibullCDF(t, k, lambda float64) float64 {
	if k <= 0 || lambda <= 0 || t < 0 {
		return math.NaN()
	}
	return -math.Expm1(-math.Pow(t/lambda, k))
}

// WeibullReliability returns the survival probability R(t)=e^{-(t/λ)^k}.
func WeibullReliability(t, k, lambda float64) float64 {
	if k <= 0 || lambda <= 0 || t < 0 {
		return math.NaN()
	}
	return math.Exp(-math.Pow(t/lambda, k))
}

// WeibullHazard returns the hazard rate h(t)=(k/λ)(t/λ)^{k-1}. For k<1 the
// hazard decreases (infant mortality), k=1 gives a constant rate and k>1 gives
// an increasing (wear-out) rate.
func WeibullHazard(t, k, lambda float64) float64 {
	if k <= 0 || lambda <= 0 || t < 0 {
		return math.NaN()
	}
	if t == 0 {
		if k < 1 {
			return math.Inf(1)
		}
		if k == 1 {
			return 1 / lambda
		}
		return 0
	}
	return (k / lambda) * math.Pow(t/lambda, k-1)
}

// WeibullCumulativeHazard returns the cumulative hazard H(t)=(t/λ)^k.
func WeibullCumulativeHazard(t, k, lambda float64) float64 {
	if k <= 0 || lambda <= 0 || t < 0 {
		return math.NaN()
	}
	return math.Pow(t/lambda, k)
}

// WeibullQuantile returns the time at which F(t)=p:
// t=λ(-ln(1-p))^{1/k}.
func WeibullQuantile(p, k, lambda float64) float64 {
	if k <= 0 || lambda <= 0 || p < 0 || p >= 1 {
		return math.NaN()
	}
	return lambda * math.Pow(-math.Log1p(-p), 1/k)
}

// WeibullMean returns the mean lifetime λΓ(1+1/k).
func WeibullMean(k, lambda float64) float64 {
	if k <= 0 || lambda <= 0 {
		return math.NaN()
	}
	return lambda * math.Gamma(1+1/k)
}

// WeibullMTTF returns the mean time to failure, identical to the mean.
func WeibullMTTF(k, lambda float64) float64 {
	return WeibullMean(k, lambda)
}

// WeibullVariance returns the variance λ²[Γ(1+2/k)-Γ(1+1/k)²].
func WeibullVariance(k, lambda float64) float64 {
	if k <= 0 || lambda <= 0 {
		return math.NaN()
	}
	g1 := math.Gamma(1 + 1/k)
	g2 := math.Gamma(1 + 2/k)
	return lambda * lambda * (g2 - g1*g1)
}

// WeibullStdDev returns the standard deviation of the Weibull distribution.
func WeibullStdDev(k, lambda float64) float64 {
	v := WeibullVariance(k, lambda)
	if math.IsNaN(v) {
		return math.NaN()
	}
	return math.Sqrt(v)
}

// WeibullMedian returns the median lifetime λ(ln2)^{1/k}.
func WeibullMedian(k, lambda float64) float64 {
	if k <= 0 || lambda <= 0 {
		return math.NaN()
	}
	return lambda * math.Pow(math.Ln2, 1/k)
}

// WeibullMode returns the mode of the Weibull distribution. For k<=1 the mode
// is 0; for k>1 it is λ((k-1)/k)^{1/k}.
func WeibullMode(k, lambda float64) float64 {
	if k <= 0 || lambda <= 0 {
		return math.NaN()
	}
	if k <= 1 {
		return 0
	}
	return lambda * math.Pow((k-1)/k, 1/k)
}

// WeibullConditionalReliability returns the probability that a component of age
// t survives an additional time x: R(t+x)/R(t).
func WeibullConditionalReliability(t, x, k, lambda float64) float64 {
	if k <= 0 || lambda <= 0 || t < 0 || x < 0 {
		return math.NaN()
	}
	return math.Exp(math.Pow(t/lambda, k) - math.Pow((t+x)/lambda, k))
}

// WeibullMeanResidualLife returns the mean remaining life given survival to age
// t, computed by numerically integrating the reliability from t to infinity.
func WeibullMeanResidualLife(t, k, lambda float64) float64 {
	if k <= 0 || lambda <= 0 || t < 0 {
		return math.NaN()
	}
	rt := WeibullReliability(t, k, lambda)
	if rt <= 0 {
		return math.NaN()
	}
	step := lambda / 8
	area := integrateToInfinity(func(u float64) float64 {
		return WeibullReliability(u, k, lambda)
	}, t, step)
	return area / rt
}

// WeibullCharacteristicLife returns the scale parameter λ, the time by which
// 63.2% of a Weibull-distributed population has failed.
func WeibullCharacteristicLife(lambda float64) float64 {
	if lambda <= 0 {
		return math.NaN()
	}
	return lambda
}

// WeibullB10Life returns the B10 life, the time by which 10% of the population
// has failed (the 0.10 quantile).
func WeibullB10Life(k, lambda float64) float64 {
	return WeibullQuantile(0.10, k, lambda)
}
