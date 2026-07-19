package fuzzy

import (
	"errors"
	"math"
)

// ErrNoArea is returned by defuzzifiers when the aggregated fuzzy set has zero
// total membership so no representative crisp value exists.
var ErrNoArea = errors.New("fuzzy: cannot defuzzify a set with zero area")

// Centroid returns the center of gravity of the discrete fuzzy set,
// sum(X[i]*Mu[i]) / sum(Mu[i]). It returns ErrNoArea when every grade is 0.
func (s Set) Centroid() (float64, error) {
	num := 0.0
	den := 0.0
	for i := range s.X {
		num += s.X[i] * s.Mu[i]
		den += s.Mu[i]
	}
	if den == 0 {
		return 0, ErrNoArea
	}
	return num / den, nil
}

// CentroidTrapz returns the continuous centroid of the fuzzy set treating the
// membership as piecewise linear between universe points and integrating with
// the trapezoidal rule, integral(x*mu dx) / integral(mu dx). It returns
// ErrNoArea when the enclosed area is 0.
func (s Set) CentroidTrapz() (float64, error) {
	num := 0.0
	den := 0.0
	for i := 0; i+1 < len(s.X); i++ {
		x0, x1 := s.X[i], s.X[i+1]
		m0, m1 := s.Mu[i], s.Mu[i+1]
		h := x1 - x0
		// area of trapezoid
		area := (m0 + m1) / 2 * h
		// moment of trapezoid about origin (exact for linear segment)
		moment := h * (x0*(2*m0+m1) + x1*(m0+2*m1)) / 6
		num += moment
		den += area
	}
	if den == 0 {
		return 0, ErrNoArea
	}
	return num / den, nil
}

// Bisector returns the bisector of area of the discrete fuzzy set, the universe
// point that splits the total membership sum into two (as equal as possible)
// halves. It returns ErrNoArea when every grade is 0.
func (s Set) Bisector() (float64, error) {
	total := s.Cardinality()
	if total == 0 {
		return 0, ErrNoArea
	}
	half := total / 2
	acc := 0.0
	for i := range s.X {
		acc += s.Mu[i]
		if acc >= half {
			return s.X[i], nil
		}
	}
	return s.X[len(s.X)-1], nil
}

// maxIndices returns the indices at which the membership attains its maximum
// (within a small tolerance) together with that maximum value.
func (s Set) maxIndices() ([]int, float64) {
	h := s.Height()
	if h == 0 {
		return nil, 0
	}
	const eps = 1e-12
	var idx []int
	for i, m := range s.Mu {
		if m >= h-eps {
			idx = append(idx, i)
		}
	}
	return idx, h
}

// MeanOfMaxima returns the mean of maxima (MOM), the average of the universe
// points at which the membership attains its maximum. It returns ErrNoArea when
// every grade is 0.
func (s Set) MeanOfMaxima() (float64, error) {
	idx, h := s.maxIndices()
	if h == 0 {
		return 0, ErrNoArea
	}
	sum := 0.0
	for _, i := range idx {
		sum += s.X[i]
	}
	return sum / float64(len(idx)), nil
}

// SmallestOfMaxima returns the smallest of maxima (SOM), the least universe
// point at which the membership attains its maximum. It returns ErrNoArea when
// every grade is 0.
func (s Set) SmallestOfMaxima() (float64, error) {
	idx, h := s.maxIndices()
	if h == 0 {
		return 0, ErrNoArea
	}
	best := math.Inf(1)
	for _, i := range idx {
		if s.X[i] < best {
			best = s.X[i]
		}
	}
	return best, nil
}

// LargestOfMaxima returns the largest of maxima (LOM), the greatest universe
// point at which the membership attains its maximum. It returns ErrNoArea when
// every grade is 0.
func (s Set) LargestOfMaxima() (float64, error) {
	idx, h := s.maxIndices()
	if h == 0 {
		return 0, ErrNoArea
	}
	best := math.Inf(-1)
	for _, i := range idx {
		if s.X[i] > best {
			best = s.X[i]
		}
	}
	return best, nil
}

// WeightedAverage returns the weighted average defuzzification
// sum(X[i]*Mu[i]) / sum(Mu[i]). For a discrete set it coincides with Centroid;
// it is provided as the standard name used by Sugeno style aggregation. It
// returns ErrNoArea when every grade is 0.
func (s Set) WeightedAverage() (float64, error) { return s.Centroid() }

// Defuzz names a defuzzification strategy.
type Defuzz int

const (
	// DefuzzCentroid selects the center of gravity.
	DefuzzCentroid Defuzz = iota
	// DefuzzBisector selects the bisector of area.
	DefuzzBisector
	// DefuzzMOM selects the mean of maxima.
	DefuzzMOM
	// DefuzzSOM selects the smallest of maxima.
	DefuzzSOM
	// DefuzzLOM selects the largest of maxima.
	DefuzzLOM
)

// Defuzzify applies the named defuzzification strategy to the set and returns
// the crisp result. It returns ErrNoArea when the set has zero area.
func (s Set) Defuzzify(method Defuzz) (float64, error) {
	switch method {
	case DefuzzCentroid:
		return s.Centroid()
	case DefuzzBisector:
		return s.Bisector()
	case DefuzzMOM:
		return s.MeanOfMaxima()
	case DefuzzSOM:
		return s.SmallestOfMaxima()
	case DefuzzLOM:
		return s.LargestOfMaxima()
	default:
		return s.Centroid()
	}
}

// String returns the canonical name of the defuzzification strategy.
func (d Defuzz) String() string {
	switch d {
	case DefuzzCentroid:
		return "centroid"
	case DefuzzBisector:
		return "bisector"
	case DefuzzMOM:
		return "mom"
	case DefuzzSOM:
		return "som"
	case DefuzzLOM:
		return "lom"
	default:
		return "unknown"
	}
}

// CentroidMF returns the centroid of the membership function mf sampled at n
// points over the closed interval [a, b] using the trapezoidal rule. It returns
// ErrNoArea when the enclosed area is 0.
func CentroidMF(mf MF, a, b float64, n int) (float64, error) {
	return FromMF(mf, Linspace(a, b, n)).CentroidTrapz()
}
