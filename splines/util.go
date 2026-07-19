package splines

import "math"

// LinearInterpolate returns the piecewise-linear interpolant of the samples
// (x[i], y[i]) evaluated at q, clamping to the end values outside the data
// range. The abscissae x must be strictly increasing with len(x) >= 2 and
// len(x) == len(y).
func LinearInterpolate(x, y []float64, q float64) (float64, error) {
	if len(x) != len(y) {
		return 0, ErrLenMismatch
	}
	if len(x) < 2 {
		return 0, ErrTooFewPoints
	}
	if !strictlyIncreasing(x) {
		return 0, ErrNotIncreasing
	}
	if q <= x[0] {
		return y[0], nil
	}
	if q >= x[len(x)-1] {
		return y[len(y)-1], nil
	}
	i := searchInterval(x, q)
	t := (q - x[i]) / (x[i+1] - x[i])
	return y[i] + t*(y[i+1]-y[i]), nil
}

// NearestInterpolate returns the sample ordinate whose abscissa is closest to q
// (nearest-neighbour interpolation). Ties choose the left sample.
func NearestInterpolate(x, y []float64, q float64) (float64, error) {
	if len(x) != len(y) {
		return 0, ErrLenMismatch
	}
	if len(x) < 1 {
		return 0, ErrTooFewPoints
	}
	if len(x) == 1 {
		return y[0], nil
	}
	if !strictlyIncreasing(x) {
		return 0, ErrNotIncreasing
	}
	i := searchInterval(x, q)
	if math.Abs(q-x[i]) <= math.Abs(x[i+1]-q) {
		return y[i], nil
	}
	return y[i+1], nil
}

// StepInterpolate returns a left-continuous piecewise-constant interpolant: for
// q in [x[i], x[i+1]) it returns y[i], clamping outside the range.
func StepInterpolate(x, y []float64, q float64) (float64, error) {
	if len(x) != len(y) {
		return 0, ErrLenMismatch
	}
	if len(x) < 1 {
		return 0, ErrTooFewPoints
	}
	if len(x) == 1 {
		return y[0], nil
	}
	if !strictlyIncreasing(x) {
		return 0, ErrNotIncreasing
	}
	if q < x[0] {
		return y[0], nil
	}
	if q >= x[len(x)-1] {
		return y[len(y)-1], nil
	}
	i := searchInterval(x, q)
	return y[i], nil
}

// SampleCurve returns n+1 points sampled at uniformly spaced parameter values
// across the curve's domain (n >= 1).
func SampleCurve(c Curve, n int) []Vec {
	if n < 1 {
		n = 1
	}
	lo, hi := c.Domain()
	pts := make([]Vec, n+1)
	for i := 0; i <= n; i++ {
		u := lo + (hi-lo)*float64(i)/float64(n)
		pts[i] = c.Eval(u)
	}
	return pts
}

// PolylineLength returns the total length of the polyline through pts, the sum
// of the Euclidean distances between consecutive points.
func PolylineLength(pts []Vec) float64 {
	var s float64
	for i := 1; i < len(pts); i++ {
		s += pts[i].Dist(pts[i-1])
	}
	return s
}

// ChordLengths returns the cumulative chord-length parameter values of pts,
// starting at zero; the final entry equals the total polyline length.
func ChordLengths(pts []Vec) []float64 {
	t := make([]float64, len(pts))
	for i := 1; i < len(pts); i++ {
		t[i] = t[i-1] + pts[i].Dist(pts[i-1])
	}
	return t
}

// CentripetalParams returns parameter values for pts using the exponent alpha
// on successive chord lengths: t[0]=0 and t[i]=t[i-1]+|P[i]-P[i-1]|^alpha. With
// alpha=0.5 this yields the centripetal parameterisation commonly used to
// choose interpolation knots.
func CentripetalParams(pts []Vec, alpha float64) []float64 {
	t := make([]float64, len(pts))
	for i := 1; i < len(pts); i++ {
		t[i] = t[i-1] + math.Pow(pts[i].Dist(pts[i-1]), alpha)
	}
	return t
}

// NormalizedChordParams returns chord-length parameters rescaled to the unit
// interval [0,1]. If all points coincide it returns evenly spaced values.
func NormalizedChordParams(pts []Vec) []float64 {
	t := ChordLengths(pts)
	total := t[len(t)-1]
	if total == 0 {
		for i := range t {
			t[i] = float64(i) / float64(len(t)-1)
		}
		return t
	}
	for i := range t {
		t[i] /= total
	}
	return t
}
