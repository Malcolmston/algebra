package reliability

import "math"

// GompertzHazard returns the exponentially increasing hazard rate
// h(t)=a·e^{bt} of the Gompertz lifetime distribution, with a>0 and b>0.
func GompertzHazard(t, a, b float64) float64 {
	if a <= 0 || b <= 0 || t < 0 {
		return math.NaN()
	}
	return a * math.Exp(b*t)
}

// GompertzCumulativeHazard returns the cumulative hazard
// H(t)=(a/b)(e^{bt}-1) of the Gompertz distribution.
func GompertzCumulativeHazard(t, a, b float64) float64 {
	if a <= 0 || b <= 0 || t < 0 {
		return math.NaN()
	}
	return (a / b) * math.Expm1(b*t)
}

// GompertzReliability returns the survival probability
// R(t)=exp(-(a/b)(e^{bt}-1)).
func GompertzReliability(t, a, b float64) float64 {
	h := GompertzCumulativeHazard(t, a, b)
	if math.IsNaN(h) {
		return math.NaN()
	}
	return math.Exp(-h)
}

// GompertzCDF returns the failure probability F(t)=1-R(t).
func GompertzCDF(t, a, b float64) float64 {
	r := GompertzReliability(t, a, b)
	if math.IsNaN(r) {
		return math.NaN()
	}
	return 1 - r
}

// GompertzPDF returns the probability density h(t)R(t) of the Gompertz
// distribution.
func GompertzPDF(t, a, b float64) float64 {
	h := GompertzHazard(t, a, b)
	r := GompertzReliability(t, a, b)
	if math.IsNaN(h) || math.IsNaN(r) {
		return math.NaN()
	}
	return h * r
}

// GompertzQuantile returns the time at which F(t)=p:
// t=(1/b)ln(1-(b/a)ln(1-p)).
func GompertzQuantile(p, a, b float64) float64 {
	if a <= 0 || b <= 0 || p < 0 || p >= 1 {
		return math.NaN()
	}
	inner := 1 - (b/a)*math.Log1p(-p)
	if inner <= 0 {
		return math.Inf(1)
	}
	return math.Log(inner) / b
}

// GompertzMedian returns the median lifetime, the 0.5 quantile.
func GompertzMedian(a, b float64) float64 {
	return GompertzQuantile(0.5, a, b)
}

// GompertzMean returns the mean lifetime, computed by numerically integrating
// the reliability from 0 to infinity.
func GompertzMean(a, b float64) float64 {
	if a <= 0 || b <= 0 {
		return math.NaN()
	}
	step := 1.0 / b
	return integrateToInfinity(func(u float64) float64 {
		return GompertzReliability(u, a, b)
	}, 0, step)
}

// GompertzMakehamHazard returns the hazard rate h(t)=λ+a·e^{bt} of the
// Gompertz–Makeham law, which adds an age-independent Makeham term lambda>=0
// to the Gompertz hazard.
func GompertzMakehamHazard(t, lambda, a, b float64) float64 {
	if lambda < 0 || a <= 0 || b <= 0 || t < 0 {
		return math.NaN()
	}
	return lambda + a*math.Exp(b*t)
}

// GompertzMakehamCumulativeHazard returns the cumulative hazard
// H(t)=λt+(a/b)(e^{bt}-1).
func GompertzMakehamCumulativeHazard(t, lambda, a, b float64) float64 {
	if lambda < 0 || a <= 0 || b <= 0 || t < 0 {
		return math.NaN()
	}
	return lambda*t + (a/b)*math.Expm1(b*t)
}

// GompertzMakehamReliability returns the survival probability R(t)=exp(-H(t))
// of the Gompertz–Makeham law.
func GompertzMakehamReliability(t, lambda, a, b float64) float64 {
	h := GompertzMakehamCumulativeHazard(t, lambda, a, b)
	if math.IsNaN(h) {
		return math.NaN()
	}
	return math.Exp(-h)
}

// GompertzMakehamCDF returns the failure probability F(t)=1-R(t).
func GompertzMakehamCDF(t, lambda, a, b float64) float64 {
	r := GompertzMakehamReliability(t, lambda, a, b)
	if math.IsNaN(r) {
		return math.NaN()
	}
	return 1 - r
}

// GompertzMakehamPDF returns the probability density h(t)R(t).
func GompertzMakehamPDF(t, lambda, a, b float64) float64 {
	h := GompertzMakehamHazard(t, lambda, a, b)
	r := GompertzMakehamReliability(t, lambda, a, b)
	if math.IsNaN(h) || math.IsNaN(r) {
		return math.NaN()
	}
	return h * r
}
