package dynamical

import "math"

// DetectPeriod attempts to determine the period of the orbit of f starting at
// x0. It first iterates transient times to settle onto the attractor, records
// the resulting reference point, then searches for the smallest period p in
// 1..maxPeriod such that the p-th further iterate returns to within tol of the
// reference point. It returns that period, or 0 if none is found within
// maxPeriod (for example on a chaotic or very long orbit).
func DetectPeriod(f Map1D, x0 float64, transient, maxPeriod int, tol float64) int {
	x := x0
	for i := 0; i < transient; i++ {
		x = f(x)
	}
	ref := x
	y := x
	for p := 1; p <= maxPeriod; p++ {
		y = f(y)
		if math.Abs(y-ref) <= tol {
			return p
		}
	}
	return 0
}

// IsPeriodicPoint reports whether x0 is a periodic point of f with period
// exactly period, that is whether f^period(x0) returns to within tol of x0
// while no smaller positive multiple of the iterate does. A period <= 0
// returns false.
func IsPeriodicPoint(f Map1D, x0 float64, period int, tol float64) bool {
	if period <= 0 {
		return false
	}
	x := x0
	for p := 1; p <= period; p++ {
		x = f(x)
		back := math.Abs(x-x0) <= tol
		if p < period && back {
			return false
		}
		if p == period {
			return back
		}
	}
	return false
}

// DetectPeriod2D is the two-dimensional analogue of [DetectPeriod] for a map of
// the plane, measuring return in Euclidean distance.
func DetectPeriod2D(f Map2D, p0 Point2D, transient, maxPeriod int, tol float64) int {
	p := p0
	for i := 0; i < transient; i++ {
		p = f(p)
	}
	ref := p
	q := p
	for k := 1; k <= maxPeriod; k++ {
		q = f(q)
		if q.Sub(ref).Norm() <= tol {
			return k
		}
	}
	return 0
}
