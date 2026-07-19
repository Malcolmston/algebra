package stochastic

import "math"

// normalCDF returns the standard normal cumulative distribution function
// Phi(x) = P(Z <= x).
func normalCDF(x float64) float64 {
	return 0.5 * math.Erfc(-x/math.Sqrt2)
}

// normalPDF returns the standard normal probability density function.
func normalPDF(x float64) float64 {
	return math.Exp(-0.5*x*x) / math.Sqrt(2*math.Pi)
}

// NormalCDF returns the cumulative distribution function of N(mu, sigma^2) at x.
func NormalCDF(x, mu, sigma float64) float64 {
	return normalCDF((x - mu) / sigma)
}

// NormalPDF returns the probability density function of N(mu, sigma^2) at x.
func NormalPDF(x, mu, sigma float64) float64 {
	return normalPDF((x-mu)/sigma) / sigma
}

// FirstPassageTime returns the first time index and time at which a path reaches
// or crosses the given level, relative to its starting value. It returns
// ok=false if the level is never reached. It is a package-level convenience
// wrapper around Path.TimeToHit.
func FirstPassageTime(p Path, level float64) (float64, bool) {
	return p.TimeToHit(level)
}

// StoppingTime returns the first index, time and value at which the predicate
// holds along the path. It returns ok=false if the predicate never holds.
func StoppingTime(p Path, predicate func(t, x float64) bool) (index int, time, value float64, ok bool) {
	for i := range p.Values {
		if predicate(p.Times[i], p.Values[i]) {
			return i, p.Times[i], p.Values[i], true
		}
	}
	return -1, math.NaN(), math.NaN(), false
}

// LastPassageTime returns the last time at which the path value equals or
// crosses level (detected as a sign change of value-level between consecutive
// points, or an exact hit). It returns ok=false if there is no crossing.
func LastPassageTime(p Path, level float64) (float64, bool) {
	found := false
	var t float64
	for i := range p.Values {
		if p.Values[i] == level {
			t = p.Times[i]
			found = true
		}
		if i > 0 {
			if (p.Values[i-1]-level)*(p.Values[i]-level) < 0 {
				t = p.Times[i]
				found = true
			}
		}
	}
	return t, found
}

// OccupationTime returns the total time the linearly interpolated path spends in
// the closed interval [lo, hi], computed by summing the durations of grid
// sub-intervals whose midpoint lies in the band.
func OccupationTime(p Path, lo, hi float64) float64 {
	total := 0.0
	for i := 1; i < len(p.Values); i++ {
		mid := 0.5 * (p.Values[i-1] + p.Values[i])
		if mid >= lo && mid <= hi {
			total += p.Times[i] - p.Times[i-1]
		}
	}
	return total
}

// ReflectionPrincipleMaxProb returns P(max_{0<=s<=t} W(s) >= a) for a Brownian
// motion with volatility sigma and no drift, using the reflection principle:
// 2*(1 - Phi(a/(sigma*sqrt(t)))) for a >= 0.
func ReflectionPrincipleMaxProb(a, t, sigma float64) float64 {
	if a <= 0 {
		return 1
	}
	if t <= 0 {
		return 0
	}
	return 2 * (1 - normalCDF(a/(sigma*math.Sqrt(t))))
}

// BrownianMaxCDF returns the CDF P(max_{0<=s<=t} W(s) <= a) of the running
// maximum of a driftless Brownian motion with volatility sigma.
func BrownianMaxCDF(a, t, sigma float64) float64 {
	return 1 - ReflectionPrincipleMaxProb(a, t, sigma)
}

// BrownianHittingProbabilityBeforeT returns the probability that a driftless
// Brownian motion with volatility sigma reaches level a (a != 0) at some time in
// [0, t]. It equals the reflection-principle maximum (or minimum) probability.
func BrownianHittingProbabilityBeforeT(a, t, sigma float64) float64 {
	return ReflectionPrincipleMaxProb(math.Abs(a), t, sigma)
}

// InverseGaussianPDF returns the density of the inverse Gaussian distribution
// with mean mu and shape lambda at x > 0. This is the first-passage-time law of
// a Brownian motion with positive drift.
func InverseGaussianPDF(x, mu, lambda float64) float64 {
	if x <= 0 {
		return 0
	}
	c := math.Sqrt(lambda / (2 * math.Pi * x * x * x))
	e := -lambda * (x - mu) * (x - mu) / (2 * mu * mu * x)
	return c * math.Exp(e)
}

// InverseGaussianCDF returns the cumulative distribution function of the inverse
// Gaussian distribution with mean mu and shape lambda at x > 0.
func InverseGaussianCDF(x, mu, lambda float64) float64 {
	if x <= 0 {
		return 0
	}
	s := math.Sqrt(lambda / x)
	a := s * (x/mu - 1)
	b := -s * (x/mu + 1)
	return normalCDF(a) + math.Exp(2*lambda/mu)*normalCDF(b)
}

// InverseGaussianMean returns the mean mu of the inverse Gaussian distribution.
func InverseGaussianMean(mu, lambda float64) float64 { return mu }

// InverseGaussianVariance returns the variance mu^3/lambda of the inverse
// Gaussian distribution.
func InverseGaussianVariance(mu, lambda float64) float64 {
	return mu * mu * mu / lambda
}

// FirstPassageDriftMean returns the mean first-passage time a/drift of a
// Brownian motion X(t) = drift*t + sigma*W(t) to a level a > 0, valid for
// drift > 0.
func FirstPassageDriftMean(a, drift float64) float64 {
	if drift <= 0 {
		return math.Inf(1)
	}
	return a / drift
}

// FirstPassageDriftCDF returns P(T_a <= t), the probability that a Brownian
// motion with drift and volatility sigma reaches level a > 0 by time t, using
// the inverse Gaussian law with mean a/drift and shape a^2/sigma^2.
func FirstPassageDriftCDF(a, drift, sigma, t float64) float64 {
	if drift <= 0 || a <= 0 {
		return math.NaN()
	}
	mu := a / drift
	lambda := a * a / (sigma * sigma)
	return InverseGaussianCDF(t, mu, lambda)
}

// WaldExpectedStopping returns E[T] = E[S_T]/E[step] from Wald's identity, given
// the expected total displacement at stopping and the mean single-step
// increment (which must be non-zero).
func WaldExpectedStopping(expectedDisplacement, meanStep float64) float64 {
	if meanStep == 0 {
		return math.Inf(1)
	}
	return expectedDisplacement / meanStep
}
