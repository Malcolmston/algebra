package meshgen

import "math"

// TriangleQuality bundles several shape-quality measures of a single triangle.
type TriangleQuality struct {
	Area         float64    // triangle area
	Perimeter    float64    // sum of edge lengths
	MinAngle     float64    // smallest interior angle, radians
	MaxAngle     float64    // largest interior angle, radians
	Angles       [3]float64 // interior angles at vertices a, b, c
	MinEdge      float64    // shortest edge length
	MaxEdge      float64    // longest edge length
	Inradius     float64    // radius of the inscribed circle
	Circumradius float64    // radius of the circumscribed circle
	AspectRatio  float64    // longest edge / (2*inradius)
	RadiusRatio  float64    // 2*inradius/circumradius, in [0,1]; 1 is equilateral
	EdgeRatio    float64    // longest edge / shortest edge, >= 1
	Shape        float64    // 4*sqrt(3)*area / sum(edge^2), in [0,1]; 1 is equilateral
}

// TriangleQualityOf computes the quality measures of the triangle (a, b, c).
func TriangleQualityOf(a, b, c Vec2) TriangleQuality {
	la := b.Distance(c)
	lb := a.Distance(c)
	lc := a.Distance(b)
	area := TriangleArea(a, b, c)
	per := la + lb + lc
	angs := TriAngles(a, b, c)
	minA := math.Min(angs[0], math.Min(angs[1], angs[2]))
	maxA := math.Max(angs[0], math.Max(angs[1], angs[2]))
	minE := math.Min(la, math.Min(lb, lc))
	maxE := math.Max(la, math.Max(lb, lc))
	inr := 0.0
	if per > 0 {
		inr = area / (per / 2)
	}
	var circ float64
	if r, err := Circumradius(a, b, c); err == nil {
		circ = r
	}
	aspect := math.Inf(1)
	if inr > 0 {
		aspect = maxE / (2 * inr)
	}
	radRatio := 0.0
	if circ > 0 {
		radRatio = 2 * inr / circ
	}
	edgeRatio := math.Inf(1)
	if minE > 0 {
		edgeRatio = maxE / minE
	}
	shape := 0.0
	sumSq := la*la + lb*lb + lc*lc
	if sumSq > 0 {
		shape = 4 * math.Sqrt(3) * area / sumSq
	}
	return TriangleQuality{
		Area:         area,
		Perimeter:    per,
		MinAngle:     minA,
		MaxAngle:     maxA,
		Angles:       angs,
		MinEdge:      minE,
		MaxEdge:      maxE,
		Inradius:     inr,
		Circumradius: circ,
		AspectRatio:  aspect,
		RadiusRatio:  radRatio,
		EdgeRatio:    edgeRatio,
		Shape:        shape,
	}
}

// TriAspectRatio returns the aspect ratio (longest edge over twice the
// inradius) of the triangle; 1 for an equilateral triangle up to scaling.
func TriAspectRatio(a, b, c Vec2) float64 { return TriangleQualityOf(a, b, c).AspectRatio }

// TriRadiusRatio returns the radius ratio 2*inradius/circumradius, which lies
// in [0,1] and equals 1 for an equilateral triangle.
func TriRadiusRatio(a, b, c Vec2) float64 { return TriangleQualityOf(a, b, c).RadiusRatio }

// TriShapeRegularity returns 4*sqrt(3)*area / (sum of squared edge lengths),
// a scale-invariant quality in [0,1] that is 1 for an equilateral triangle.
func TriShapeRegularity(a, b, c Vec2) float64 { return TriangleQualityOf(a, b, c).Shape }

// TriEdgeRatio returns the ratio of the longest to the shortest edge length.
func TriEdgeRatio(a, b, c Vec2) float64 { return TriangleQualityOf(a, b, c).EdgeRatio }

// IsWellShaped reports whether the triangle's minimum angle is at least
// minAngleDeg degrees.
func IsWellShaped(a, b, c Vec2, minAngleDeg float64) bool {
	return TriMinAngle(a, b, c) >= minAngleDeg*math.Pi/180
}

// IsSliver reports whether the triangle is a sliver: it has positive area but a
// minimum angle below minAngleDeg degrees.
func IsSliver(a, b, c Vec2, minAngleDeg float64) bool {
	if TriangleArea(a, b, c) <= 0 {
		return false
	}
	return TriMinAngle(a, b, c) < minAngleDeg*math.Pi/180
}

// TriangleQualityAt returns the quality of triangle i of the mesh.
func (m *Mesh) TriangleQualityAt(i int) TriangleQuality {
	a, b, c := m.TriangleVertices(i)
	return TriangleQualityOf(a, b, c)
}

// MinAngleDeg returns the smallest interior angle over all triangles of the
// mesh, in degrees. An empty mesh returns 0.
func (m *Mesh) MinAngleDeg() float64 {
	if len(m.Triangles) == 0 {
		return 0
	}
	minA := math.Inf(1)
	for i := range m.Triangles {
		a, b, c := m.TriangleVertices(i)
		if v := TriMinAngle(a, b, c); v < minA {
			minA = v
		}
	}
	return minA * 180 / math.Pi
}

// MaxAngleDeg returns the largest interior angle over all triangles of the
// mesh, in degrees. An empty mesh returns 0.
func (m *Mesh) MaxAngleDeg() float64 {
	if len(m.Triangles) == 0 {
		return 0
	}
	maxA := math.Inf(-1)
	for i := range m.Triangles {
		a, b, c := m.TriangleVertices(i)
		if v := TriMaxAngle(a, b, c); v > maxA {
			maxA = v
		}
	}
	return maxA * 180 / math.Pi
}

// CountSlivers returns the number of triangles whose minimum angle is below
// minAngleDeg degrees.
func (m *Mesh) CountSlivers(minAngleDeg float64) int {
	n := 0
	for i := range m.Triangles {
		a, b, c := m.TriangleVertices(i)
		if IsSliver(a, b, c, minAngleDeg) {
			n++
		}
	}
	return n
}
