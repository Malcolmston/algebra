package quasirandom

import "math"

// Integrand is a scalar function of a point in R^d used by the integration
// helpers.
type Integrand func(x []float64) float64

// IntegrateUnitCube estimates the integral of f over the unit cube [0,1]^d as
// the equal-weight average of f over the supplied quasi-random points. It
// returns an error when no points are given.
func IntegrateUnitCube(f Integrand, points [][]float64) (float64, error) {
	if len(points) == 0 {
		return 0, ErrEmptyPointSet
	}
	var sum float64
	for _, p := range points {
		sum += f(p)
	}
	return sum / float64(len(points)), nil
}

// IntegrateBox estimates the integral of f over the axis-aligned box
// [lower,upper] by mapping the unit-cube points into the box and multiplying the
// average by the box volume. lower and upper must have the same length as each
// point. It returns an error when the shapes disagree or no points are given.
func IntegrateBox(f Integrand, lower, upper []float64, points [][]float64) (float64, error) {
	if len(points) == 0 {
		return 0, ErrEmptyPointSet
	}
	d := len(points[0])
	if len(lower) != d || len(upper) != d {
		return 0, ErrDimension
	}
	vol := 1.0
	for k := 0; k < d; k++ {
		vol *= upper[k] - lower[k]
	}
	mapped := make([]float64, d)
	var sum float64
	for _, p := range points {
		if len(p) != d {
			return 0, ErrRaggedPointSet
		}
		for k := 0; k < d; k++ {
			mapped[k] = lower[k] + p[k]*(upper[k]-lower[k])
		}
		sum += f(mapped)
	}
	return vol * sum / float64(len(points)), nil
}

// QMCIntegrateHalton estimates the integral of f over [0,1]^dim using the first
// n points of the Halton sequence. It returns an error when dim < 1 or n < 1.
func QMCIntegrateHalton(f Integrand, dim, n int) (float64, error) {
	if n < 1 {
		return 0, ErrNonPositive
	}
	pts, err := HaltonSequence(dim, n)
	if err != nil {
		return 0, err
	}
	return IntegrateUnitCube(f, pts)
}

// QMCIntegrateSobol estimates the integral of f over [0,1]^dim using the first
// n Sobol points after skipping the origin. It returns an error when dim is out
// of range or n < 1.
func QMCIntegrateSobol(f Integrand, dim, n int) (float64, error) {
	if n < 1 {
		return 0, ErrNonPositive
	}
	pts, err := SobolSequenceSkip(dim, 1, n)
	if err != nil {
		return 0, err
	}
	return IntegrateUnitCube(f, pts)
}

// QMCIntegrateFaure estimates the integral of f over [0,1]^dim using the first
// n points of the Faure sequence. It returns an error when dim < 1 or n < 1.
func QMCIntegrateFaure(f Integrand, dim, n int) (float64, error) {
	if n < 1 {
		return 0, ErrNonPositive
	}
	pts, err := FaureSequence(dim, n)
	if err != nil {
		return 0, err
	}
	return IntegrateUnitCube(f, pts)
}

// QMCIntegrateHammersley estimates the integral of f over [0,1]^dim using the
// n-point Hammersley set. It returns an error when dim < 1 or n < 1.
func QMCIntegrateHammersley(f Integrand, dim, n int) (float64, error) {
	if n < 1 {
		return 0, ErrNonPositive
	}
	pts, err := HammersleySet(dim, n)
	if err != nil {
		return 0, err
	}
	return IntegrateUnitCube(f, pts)
}

// SampleMean returns the arithmetic mean of the values, or zero for an empty
// slice.
func SampleMean(values []float64) float64 {
	if len(values) == 0 {
		return 0
	}
	var s float64
	for _, v := range values {
		s += v
	}
	return s / float64(len(values))
}

// SampleVariance returns the unbiased (Bessel-corrected) sample variance of the
// values, or zero when fewer than two values are given.
func SampleVariance(values []float64) float64 {
	n := len(values)
	if n < 2 {
		return 0
	}
	m := SampleMean(values)
	var s float64
	for _, v := range values {
		d := v - m
		s += d * d
	}
	return s / float64(n-1)
}

// MonteCarloStandardError returns the estimated standard error of the sample
// mean, sqrt(SampleVariance/n). It is the classical error estimate for plain
// Monte Carlo and an (over-)conservative one for quasi-Monte-Carlo.
func MonteCarloStandardError(values []float64) float64 {
	n := len(values)
	if n < 2 {
		return 0
	}
	return math.Sqrt(SampleVariance(values) / float64(n))
}

// EvaluateSamples returns the slice of function values f(p) for each point p,
// the raw samples underlying the QMC average and its error estimates.
func EvaluateSamples(f Integrand, points [][]float64) []float64 {
	out := make([]float64, len(points))
	for i, p := range points {
		out[i] = f(p)
	}
	return out
}
