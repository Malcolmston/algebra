package reliability

import "math"

// BathtubHazard evaluates an additive bathtub-shaped hazard rate at time t>=0:
//
//	h(t) = a·e^{-b t} + c + e·t
//
// The first term (a>=0, b>0) models decreasing infant-mortality failures, the
// constant c>=0 the useful-life period, and the linear term e>=0 the
// increasing wear-out phase. Suitable parameter choices reproduce the classic
// bathtub curve.
func BathtubHazard(t, a, b, c, e float64) float64 {
	if a < 0 || b <= 0 || c < 0 || e < 0 || t < 0 {
		return math.NaN()
	}
	return a*math.Exp(-b*t) + c + e*t
}

// BathtubCumulativeHazard returns the cumulative hazard of the additive
// bathtub model,
//
//	H(t) = (a/b)(1-e^{-b t}) + c·t + (e/2)·t².
func BathtubCumulativeHazard(t, a, b, c, e float64) float64 {
	if a < 0 || b <= 0 || c < 0 || e < 0 || t < 0 {
		return math.NaN()
	}
	return (a/b)*(-math.Expm1(-b*t)) + c*t + 0.5*e*t*t
}

// BathtubReliability returns the survival probability R(t)=exp(-H(t)) of the
// additive bathtub model.
func BathtubReliability(t, a, b, c, e float64) float64 {
	h := BathtubCumulativeHazard(t, a, b, c, e)
	if math.IsNaN(h) {
		return math.NaN()
	}
	return math.Exp(-h)
}

// BathtubCDF returns the failure probability F(t)=1-R(t) of the additive
// bathtub model.
func BathtubCDF(t, a, b, c, e float64) float64 {
	r := BathtubReliability(t, a, b, c, e)
	if math.IsNaN(r) {
		return math.NaN()
	}
	return 1 - r
}

// BathtubPDF returns the probability density h(t)R(t) of the additive bathtub
// model.
func BathtubPDF(t, a, b, c, e float64) float64 {
	h := BathtubHazard(t, a, b, c, e)
	r := BathtubReliability(t, a, b, c, e)
	if math.IsNaN(h) || math.IsNaN(r) {
		return math.NaN()
	}
	return h * r
}

// HjorthHazard evaluates the Hjorth-distribution hazard rate at time t>=0:
//
//	h(t) = δ·t + β/(1+θ·t)
//
// with δ>=0, β>=0 and θ>0. Depending on the parameters the hazard can be
// increasing, decreasing, or bathtub-shaped, making the Hjorth model a compact
// alternative to the additive bathtub curve.
func HjorthHazard(t, delta, beta, theta float64) float64 {
	if delta < 0 || beta < 0 || theta <= 0 || t < 0 {
		return math.NaN()
	}
	return delta*t + beta/(1+theta*t)
}

// HjorthCumulativeHazard returns the cumulative hazard of the Hjorth model,
//
//	H(t) = δ·t²/2 + (β/θ)·ln(1+θ·t).
func HjorthCumulativeHazard(t, delta, beta, theta float64) float64 {
	if delta < 0 || beta < 0 || theta <= 0 || t < 0 {
		return math.NaN()
	}
	return 0.5*delta*t*t + (beta/theta)*math.Log1p(theta*t)
}

// HjorthReliability returns the survival probability
// R(t)=exp(-δt²/2)·(1+θt)^{-β/θ}.
func HjorthReliability(t, delta, beta, theta float64) float64 {
	h := HjorthCumulativeHazard(t, delta, beta, theta)
	if math.IsNaN(h) {
		return math.NaN()
	}
	return math.Exp(-h)
}

// HjorthCDF returns the failure probability F(t)=1-R(t) of the Hjorth model.
func HjorthCDF(t, delta, beta, theta float64) float64 {
	r := HjorthReliability(t, delta, beta, theta)
	if math.IsNaN(r) {
		return math.NaN()
	}
	return 1 - r
}

// HjorthPDF returns the probability density h(t)R(t) of the Hjorth model.
func HjorthPDF(t, delta, beta, theta float64) float64 {
	h := HjorthHazard(t, delta, beta, theta)
	r := HjorthReliability(t, delta, beta, theta)
	if math.IsNaN(h) || math.IsNaN(r) {
		return math.NaN()
	}
	return h * r
}
