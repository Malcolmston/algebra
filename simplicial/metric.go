package simplicial

import "math"

// Metric is a distance function between two points of equal length. The Rips and
// Čech constructors accept a Metric so that complexes can be built from any
// dissimilarity, not only the Euclidean one.
type Metric func(a, b []float64) float64

// EuclideanDistance returns the L2 distance between two equal-length points. It
// returns NaN if the lengths differ.
func EuclideanDistance(a, b []float64) float64 {
	if len(a) != len(b) {
		return math.NaN()
	}
	var s float64
	for i := range a {
		d := a[i] - b[i]
		s += d * d
	}
	return math.Sqrt(s)
}

// SquaredEuclideanDistance returns the squared L2 distance, avoiding the square
// root. It returns NaN if the lengths differ.
func SquaredEuclideanDistance(a, b []float64) float64 {
	if len(a) != len(b) {
		return math.NaN()
	}
	var s float64
	for i := range a {
		d := a[i] - b[i]
		s += d * d
	}
	return s
}

// ManhattanDistance returns the L1 (taxicab) distance. It returns NaN if the
// lengths differ.
func ManhattanDistance(a, b []float64) float64 {
	if len(a) != len(b) {
		return math.NaN()
	}
	var s float64
	for i := range a {
		s += math.Abs(a[i] - b[i])
	}
	return s
}

// ChebyshevDistance returns the L∞ (maximum coordinate) distance. It returns NaN
// if the lengths differ.
func ChebyshevDistance(a, b []float64) float64 {
	if len(a) != len(b) {
		return math.NaN()
	}
	var m float64
	for i := range a {
		if d := math.Abs(a[i] - b[i]); d > m {
			m = d
		}
	}
	return m
}

// MinkowskiDistance returns the L_p distance for p ≥ 1. For p = 1 it is the
// Manhattan distance, for p = 2 the Euclidean distance. It returns NaN if the
// lengths differ.
func MinkowskiDistance(a, b []float64, p float64) float64 {
	if len(a) != len(b) {
		return math.NaN()
	}
	var s float64
	for i := range a {
		s += math.Pow(math.Abs(a[i]-b[i]), p)
	}
	return math.Pow(s, 1/p)
}
