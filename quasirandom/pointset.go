package quasirandom

import "math"

// BoundingBox returns the per-coordinate minimum and maximum of a point set.
// It returns an error for an empty or ragged set.
func BoundingBox(points [][]float64) (lo, hi []float64, err error) {
	if len(points) == 0 {
		return nil, nil, ErrEmptyPointSet
	}
	d := len(points[0])
	lo = make([]float64, d)
	hi = make([]float64, d)
	copy(lo, points[0])
	copy(hi, points[0])
	for _, p := range points {
		if len(p) != d {
			return nil, nil, ErrRaggedPointSet
		}
		for k := 0; k < d; k++ {
			if p[k] < lo[k] {
				lo[k] = p[k]
			}
			if p[k] > hi[k] {
				hi[k] = p[k]
			}
		}
	}
	return lo, hi, nil
}

// Centroid returns the coordinate-wise mean of a point set. It returns an error
// for an empty or ragged set.
func Centroid(points [][]float64) ([]float64, error) {
	if len(points) == 0 {
		return nil, ErrEmptyPointSet
	}
	d := len(points[0])
	c := make([]float64, d)
	for _, p := range points {
		if len(p) != d {
			return nil, ErrRaggedPointSet
		}
		for k := 0; k < d; k++ {
			c[k] += p[k]
		}
	}
	for k := 0; k < d; k++ {
		c[k] /= float64(len(points))
	}
	return c, nil
}

// CoordinateMean returns the mean of coordinate k over the point set. It returns
// an error for an empty set or an out-of-range coordinate index.
func CoordinateMean(points [][]float64, k int) (float64, error) {
	if len(points) == 0 {
		return 0, ErrEmptyPointSet
	}
	if k < 0 || k >= len(points[0]) {
		return 0, ErrDimension
	}
	var s float64
	for _, p := range points {
		s += p[k]
	}
	return s / float64(len(points)), nil
}

// CoordinateVariance returns the population variance of coordinate k over the
// point set. It returns an error for an empty set or an out-of-range index.
func CoordinateVariance(points [][]float64, k int) (float64, error) {
	m, err := CoordinateMean(points, k)
	if err != nil {
		return 0, err
	}
	var s float64
	for _, p := range points {
		d := p[k] - m
		s += d * d
	}
	return s / float64(len(points)), nil
}

// CoordinateMeans returns the vector of per-coordinate means of a point set. It
// returns an error for an empty or ragged set.
func CoordinateMeans(points [][]float64) ([]float64, error) {
	return Centroid(points)
}

// CoordinateVariances returns the vector of per-coordinate population variances
// of a point set. It returns an error for an empty or ragged set.
func CoordinateVariances(points [][]float64) ([]float64, error) {
	if len(points) == 0 {
		return nil, ErrEmptyPointSet
	}
	d := len(points[0])
	out := make([]float64, d)
	for k := 0; k < d; k++ {
		v, err := CoordinateVariance(points, k)
		if err != nil {
			return nil, err
		}
		out[k] = v
	}
	return out, nil
}

// EuclideanDistance returns the Euclidean distance between two equal-length
// points. It returns an error when the lengths disagree.
func EuclideanDistance(a, b []float64) (float64, error) {
	if len(a) != len(b) {
		return 0, ErrDimension
	}
	var s float64
	for k := range a {
		d := a[k] - b[k]
		s += d * d
	}
	return math.Sqrt(s), nil
}

// ToroidalDistance returns the Euclidean distance between two points on the unit
// torus, where each coordinate difference is taken to be at most one half. It
// returns an error when the lengths disagree.
func ToroidalDistance(a, b []float64) (float64, error) {
	if len(a) != len(b) {
		return 0, ErrDimension
	}
	var s float64
	for k := range a {
		d := math.Abs(a[k] - b[k])
		if d > 0.5 {
			d = 1 - d
		}
		s += d * d
	}
	return math.Sqrt(s), nil
}

// MinimumDistance returns the smallest pairwise Euclidean distance in a point
// set (twice the separation distance). It returns an error for a set with fewer
// than two points.
func MinimumDistance(points [][]float64) (float64, error) {
	n := len(points)
	if n < 2 {
		return 0, ErrEmptyPointSet
	}
	best := math.Inf(1)
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			d, err := EuclideanDistance(points[i], points[j])
			if err != nil {
				return 0, err
			}
			if d < best {
				best = d
			}
		}
	}
	return best, nil
}

// SeparationDistance returns half the minimum pairwise distance, the largest
// radius of non-overlapping balls centered at the points. It returns an error
// for a set with fewer than two points.
func SeparationDistance(points [][]float64) (float64, error) {
	d, err := MinimumDistance(points)
	if err != nil {
		return 0, err
	}
	return d / 2, nil
}

// CountInBox returns the number of points lying in the half-open anchored box
// [0,t) whose upper corner is t. It returns an error when a point and t have
// different lengths.
func CountInBox(points [][]float64, t []float64) (int, error) {
	cnt := 0
	for _, p := range points {
		if len(p) != len(t) {
			return 0, ErrDimension
		}
		inside := true
		for k := range t {
			if p[k] >= t[k] {
				inside = false
				break
			}
		}
		if inside {
			cnt++
		}
	}
	return cnt, nil
}

// LocalDiscrepancy returns the signed local discrepancy of a point set at the
// anchored box [0,t): the fraction of points inside minus the box volume. It
// returns an error for an empty set or a corner of the wrong length.
func LocalDiscrepancy(points [][]float64, t []float64) (float64, error) {
	if len(points) == 0 {
		return 0, ErrEmptyPointSet
	}
	cnt, err := CountInBox(points, t)
	if err != nil {
		return 0, err
	}
	vol := 1.0
	for _, tk := range t {
		vol *= tk
	}
	return float64(cnt)/float64(len(points)) - vol, nil
}

// Project returns the projection of a point set onto the coordinate axes listed
// in dims, producing a lower-dimensional point set. It returns an error when a
// requested index is out of range.
func Project(points [][]float64, dims []int) ([][]float64, error) {
	if len(points) == 0 {
		return nil, ErrEmptyPointSet
	}
	d := len(points[0])
	for _, k := range dims {
		if k < 0 || k >= d {
			return nil, ErrDimension
		}
	}
	out := make([][]float64, len(points))
	for i, p := range points {
		q := make([]float64, len(dims))
		for j, k := range dims {
			q[j] = p[k]
		}
		out[i] = q
	}
	return out, nil
}

// ScaleToBox maps a point set from the unit cube into the axis-aligned box
// [lower,upper]. It returns an error when the shapes disagree.
func ScaleToBox(points [][]float64, lower, upper []float64) ([][]float64, error) {
	if len(points) == 0 {
		return nil, ErrEmptyPointSet
	}
	d := len(points[0])
	if len(lower) != d || len(upper) != d {
		return nil, ErrDimension
	}
	out := make([][]float64, len(points))
	for i, p := range points {
		q := make([]float64, d)
		for k := 0; k < d; k++ {
			q[k] = lower[k] + p[k]*(upper[k]-lower[k])
		}
		out[i] = q
	}
	return out, nil
}

// ReflectPoint returns the reflection of a single point about the center of the
// unit cube, mapping each coordinate x to 1-x.
func ReflectPoint(p []float64) []float64 {
	out := make([]float64, len(p))
	for k, x := range p {
		out[k] = 1 - x
	}
	return out
}

// ReflectPoints returns the reflection of every point of a set about the center
// of the unit cube.
func ReflectPoints(points [][]float64) [][]float64 {
	out := make([][]float64, len(points))
	for i, p := range points {
		out[i] = ReflectPoint(p)
	}
	return out
}

// BoxToUnit maps a point set from the axis-aligned box [lower,upper] back into
// the unit cube, the inverse of ScaleToBox. It returns an error when the shapes
// disagree or a side has zero length.
func BoxToUnit(points [][]float64, lower, upper []float64) ([][]float64, error) {
	if len(points) == 0 {
		return nil, ErrEmptyPointSet
	}
	d := len(points[0])
	if len(lower) != d || len(upper) != d {
		return nil, ErrDimension
	}
	out := make([][]float64, len(points))
	for i, p := range points {
		q := make([]float64, d)
		for k := 0; k < d; k++ {
			w := upper[k] - lower[k]
			if w == 0 {
				return nil, ErrDimension
			}
			q[k] = (p[k] - lower[k]) / w
		}
		out[i] = q
	}
	return out, nil
}

// FillDistance returns the fill distance (covering radius) of a point set with
// respect to a regular test grid of res divisions per dimension: the largest
// distance from any grid node to its nearest sample point. It returns an error
// for an empty set or res < 1.
func FillDistance(points [][]float64, res int) (float64, error) {
	n, d, err := validatePoints(points)
	if err != nil {
		return 0, err
	}
	if res < 1 {
		return 0, ErrNonPositive
	}
	_ = n
	idx := make([]int, d)
	t := make([]float64, d)
	worst := 0.0
	for {
		for k := 0; k < d; k++ {
			t[k] = (float64(idx[k]) + 0.5) / float64(res)
		}
		nearest := math.Inf(1)
		for _, p := range points {
			dist, _ := EuclideanDistance(p, t)
			if dist < nearest {
				nearest = dist
			}
		}
		if nearest > worst {
			worst = nearest
		}
		carry := true
		for k := 0; k < d && carry; k++ {
			idx[k]++
			if idx[k] < res {
				carry = false
			} else {
				idx[k] = 0
			}
		}
		if carry {
			break
		}
	}
	return worst, nil
}
