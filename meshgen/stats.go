package meshgen

import (
	"math"
	"sort"
)

// MeshStats summarises area and quality statistics over an entire mesh.
type MeshStats struct {
	NumVertices    int
	NumTriangles   int
	NumEdges       int
	NumBoundary    int     // number of boundary edges
	TotalArea      float64 // sum of triangle areas
	MinArea        float64 // smallest triangle area
	MaxArea        float64 // largest triangle area
	MeanArea       float64 // mean triangle area
	MinAngleDeg    float64 // smallest interior angle over all triangles, degrees
	MaxAngleDeg    float64 // largest interior angle over all triangles, degrees
	MeanMinAngle   float64 // mean of per-triangle minimum angles, degrees
	MinQuality     float64 // smallest per-triangle radius ratio
	MeanQuality    float64 // mean per-triangle radius ratio
	MinEdgeLength  float64 // shortest edge length
	MaxEdgeLength  float64 // longest edge length
	MeanEdgeLength float64
}

// Stats computes area and quality statistics for the mesh.
func (m *Mesh) Stats() MeshStats {
	s := MeshStats{
		NumVertices:  len(m.Vertices),
		NumTriangles: len(m.Triangles),
	}
	edges := m.Edges()
	s.NumEdges = len(edges)
	s.NumBoundary = len(m.BoundaryEdges())
	if len(m.Triangles) > 0 {
		minArea := math.Inf(1)
		maxArea := math.Inf(-1)
		minAng := math.Inf(1)
		maxAng := math.Inf(-1)
		minQ := math.Inf(1)
		var sumArea, sumMinAng, sumQ float64
		for i := range m.Triangles {
			a, b, c := m.TriangleVertices(i)
			q := TriangleQualityOf(a, b, c)
			sumArea += q.Area
			if q.Area < minArea {
				minArea = q.Area
			}
			if q.Area > maxArea {
				maxArea = q.Area
			}
			if q.MinAngle < minAng {
				minAng = q.MinAngle
			}
			if q.MaxAngle > maxAng {
				maxAng = q.MaxAngle
			}
			if q.RadiusRatio < minQ {
				minQ = q.RadiusRatio
			}
			sumMinAng += q.MinAngle
			sumQ += q.RadiusRatio
		}
		n := float64(len(m.Triangles))
		s.TotalArea = sumArea
		s.MinArea = minArea
		s.MaxArea = maxArea
		s.MeanArea = sumArea / n
		s.MinAngleDeg = minAng * 180 / math.Pi
		s.MaxAngleDeg = maxAng * 180 / math.Pi
		s.MeanMinAngle = (sumMinAng / n) * 180 / math.Pi
		s.MinQuality = minQ
		s.MeanQuality = sumQ / n
	}
	if len(edges) > 0 {
		lens := m.EdgeLengths()
		s.MinEdgeLength = MinFloat(lens)
		s.MaxEdgeLength = MaxFloat(lens)
		s.MeanEdgeLength = MeanFloat(lens)
	}
	return s
}

// TriangleAreas returns the area of every triangle in triangle order.
func (m *Mesh) TriangleAreas() []float64 {
	out := make([]float64, len(m.Triangles))
	for i := range m.Triangles {
		out[i] = m.TriangleArea(i)
	}
	return out
}

// TriangleMinAngles returns the minimum interior angle (radians) of every
// triangle in triangle order.
func (m *Mesh) TriangleMinAngles() []float64 {
	out := make([]float64, len(m.Triangles))
	for i := range m.Triangles {
		a, b, c := m.TriangleVertices(i)
		out[i] = TriMinAngle(a, b, c)
	}
	return out
}

// TriangleQualities returns the radius-ratio quality of every triangle in
// triangle order.
func (m *Mesh) TriangleQualities() []float64 {
	out := make([]float64, len(m.Triangles))
	for i := range m.Triangles {
		a, b, c := m.TriangleVertices(i)
		out[i] = TriRadiusRatio(a, b, c)
	}
	return out
}

// AngleHistogram bins the interior angles of all triangles (in degrees) into
// the given number of equal-width bins spanning [0, 180]. It returns the bin
// counts. The number of bins must be positive.
func (m *Mesh) AngleHistogram(bins int) []int {
	if bins <= 0 {
		return nil
	}
	h := make([]int, bins)
	width := 180.0 / float64(bins)
	for i := range m.Triangles {
		a, b, c := m.TriangleVertices(i)
		for _, ang := range TriAngles(a, b, c) {
			deg := ang * 180 / math.Pi
			idx := int(deg / width)
			if idx >= bins {
				idx = bins - 1
			}
			if idx < 0 {
				idx = 0
			}
			h[idx]++
		}
	}
	return h
}

// StdDev returns the population standard deviation of the slice, or 0 for a
// slice with fewer than two elements.
func StdDev(xs []float64) float64 {
	n := len(xs)
	if n < 2 {
		return 0
	}
	mean := MeanFloat(xs)
	var s float64
	for _, x := range xs {
		d := x - mean
		s += d * d
	}
	return math.Sqrt(s / float64(n))
}

// Median returns the median of the slice, or 0 for an empty slice. The input
// is not modified.
func Median(xs []float64) float64 {
	n := len(xs)
	if n == 0 {
		return 0
	}
	cp := make([]float64, n)
	copy(cp, xs)
	sort.Float64s(cp)
	if n%2 == 1 {
		return cp[n/2]
	}
	return (cp[n/2-1] + cp[n/2]) / 2
}
