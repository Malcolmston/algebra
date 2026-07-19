package simplicial

import "math"

// PointCloud is a finite set of points in Euclidean space, all of the same
// ambient dimension. It is the geometric input to the Vietoris–Rips and Čech
// constructions.
type PointCloud struct {
	points [][]float64
}

// NewPointCloud returns a point cloud from the given points. Each point is
// copied. The points should all have the same length; mixed lengths are allowed
// but the metric routines will report NaN distances between mismatched points.
func NewPointCloud(points ...[]float64) *PointCloud {
	cp := make([][]float64, len(points))
	for i, p := range points {
		cp[i] = append([]float64(nil), p...)
	}
	return &PointCloud{points: cp}
}

// Len returns the number of points in the cloud.
func (pc *PointCloud) Len() int { return len(pc.points) }

// Dim returns the ambient dimension (the length of the first point), or 0 for an
// empty cloud.
func (pc *PointCloud) Dim() int {
	if len(pc.points) == 0 {
		return 0
	}
	return len(pc.points[0])
}

// Point returns a copy of the i-th point.
func (pc *PointCloud) Point(i int) []float64 {
	return append([]float64(nil), pc.points[i]...)
}

// Add appends a copy of point p to the cloud and returns its index.
func (pc *PointCloud) Add(p []float64) int {
	pc.points = append(pc.points, append([]float64(nil), p...))
	return len(pc.points) - 1
}

// Distance returns the Euclidean distance between points i and j.
func (pc *PointCloud) Distance(i, j int) float64 {
	return EuclideanDistance(pc.points[i], pc.points[j])
}

// DistanceWith returns the distance between points i and j under the given
// metric.
func (pc *PointCloud) DistanceWith(i, j int, metric Metric) float64 {
	return metric(pc.points[i], pc.points[j])
}

// DistanceMatrix returns the symmetric matrix of pairwise Euclidean distances.
func (pc *PointCloud) DistanceMatrix() [][]float64 {
	return pc.DistanceMatrixWith(EuclideanDistance)
}

// DistanceMatrixWith returns the symmetric matrix of pairwise distances under
// the given metric.
func (pc *PointCloud) DistanceMatrixWith(metric Metric) [][]float64 {
	n := len(pc.points)
	d := make([][]float64, n)
	for i := range d {
		d[i] = make([]float64, n)
	}
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			dij := metric(pc.points[i], pc.points[j])
			d[i][j] = dij
			d[j][i] = dij
		}
	}
	return d
}

// Diameter returns the largest pairwise Euclidean distance in the cloud, or 0
// for a cloud with fewer than two points.
func (pc *PointCloud) Diameter() float64 {
	var m float64
	for i := 0; i < len(pc.points); i++ {
		for j := i + 1; j < len(pc.points); j++ {
			if d := pc.Distance(i, j); d > m {
				m = d
			}
		}
	}
	return m
}

// Centroid returns the coordinatewise mean of the points, or nil for an empty
// cloud.
func (pc *PointCloud) Centroid() []float64 {
	if len(pc.points) == 0 {
		return nil
	}
	d := pc.Dim()
	c := make([]float64, d)
	for _, p := range pc.points {
		for i := 0; i < d && i < len(p); i++ {
			c[i] += p[i]
		}
	}
	for i := range c {
		c[i] /= float64(len(pc.points))
	}
	return c
}

// GridPoints returns a point cloud of the integer lattice points of an
// n-dimensional box: for dims = [a,b,c] it produces a·b·c points at all integer
// coordinates 0..a−1 × 0..b−1 × 0..c−1.
func GridPoints(dims ...int) *PointCloud {
	pc := &PointCloud{}
	if len(dims) == 0 {
		return pc
	}
	total := 1
	for _, d := range dims {
		if d <= 0 {
			return pc
		}
		total *= d
	}
	for idx := 0; idx < total; idx++ {
		p := make([]float64, len(dims))
		rem := idx
		for k := len(dims) - 1; k >= 0; k-- {
			p[k] = float64(rem % dims[k])
			rem /= dims[k]
		}
		pc.points = append(pc.points, p)
	}
	return pc
}

// CirclePoints returns n points equally spaced on the circle of the given radius
// centred at the origin in the plane.
func CirclePoints(n int, radius float64) *PointCloud {
	pc := &PointCloud{}
	for i := 0; i < n; i++ {
		theta := 2 * math.Pi * float64(i) / float64(n)
		pc.points = append(pc.points, []float64{radius * math.Cos(theta), radius * math.Sin(theta)})
	}
	return pc
}
