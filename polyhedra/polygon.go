package polyhedra

import (
	"errors"
	"math"
)

// ErrDegenerate is returned by routines that cannot produce a meaningful result
// because their input is degenerate (collinear, coincident, or empty).
var ErrDegenerate = errors.New("polyhedra: degenerate input")

// TriangleNormal returns the (unnormalised) normal (b-a) × (c-a) of the triangle
// (a, b, c). Its length is twice the triangle area and its direction follows the
// right-hand rule around the vertex order.
func TriangleNormal(a, b, c Vec3) Vec3 {
	return b.Sub(a).Cross(c.Sub(a))
}

// TriangleUnitNormal returns the unit normal of triangle (a, b, c) and true, or
// the zero vector and false when the triangle is degenerate.
func TriangleUnitNormal(a, b, c Vec3) (Vec3, bool) {
	return TriangleNormal(a, b, c).Normalize()
}

// TriangleArea returns the area of the triangle with vertices a, b and c.
func TriangleArea(a, b, c Vec3) float64 {
	return 0.5 * TriangleNormal(a, b, c).Len()
}

// TriangleCentroid returns the centroid (arithmetic mean of vertices) of
// triangle (a, b, c).
func TriangleCentroid(a, b, c Vec3) Vec3 {
	return a.Add(b).Add(c).Div(3)
}

// SignedTetraVolume returns the signed volume of the tetrahedron with vertices
// a, b, c and d, equal to one sixth of (b-a) · ((c-a) × (d-a)). The sign is
// positive when (b-a, c-a, d-a) form a right-handed system.
func SignedTetraVolume(a, b, c, d Vec3) float64 {
	return b.Sub(a).Triple(c.Sub(a), d.Sub(a)) / 6
}

// TetraVolume returns the (non-negative) volume of the tetrahedron a, b, c, d.
func TetraVolume(a, b, c, d Vec3) float64 {
	return math.Abs(SignedTetraVolume(a, b, c, d))
}

// NewellNormal returns the (unnormalised) normal of the planar or near-planar
// polygon given by the ordered loop of vertices, computed with Newell's method.
// The result is robust for slightly non-planar loops and its length equals twice
// the polygon area. It returns the zero vector for fewer than three vertices.
func NewellNormal(loop []Vec3) Vec3 {
	var n Vec3
	m := len(loop)
	if m < 3 {
		return Vec3{}
	}
	for i := 0; i < m; i++ {
		cur := loop[i]
		nxt := loop[(i+1)%m]
		n.X += (cur.Y - nxt.Y) * (cur.Z + nxt.Z)
		n.Y += (cur.Z - nxt.Z) * (cur.X + nxt.X)
		n.Z += (cur.X - nxt.X) * (cur.Y + nxt.Y)
	}
	return n
}

// PolygonArea returns the area of the planar polygon given by the ordered loop
// of vertices, using Newell's normal so the result is exact for planar loops and
// well behaved for slightly non-planar ones.
func PolygonArea(loop []Vec3) float64 {
	return 0.5 * NewellNormal(loop).Len()
}

// PolygonUnitNormal returns the unit normal of the polygon loop (following the
// right-hand rule around the vertex order) and true, or the zero vector and
// false when the polygon is degenerate.
func PolygonUnitNormal(loop []Vec3) (Vec3, bool) {
	return NewellNormal(loop).Normalize()
}

// PolygonCentroid returns the area-weighted centroid of the planar polygon given
// by the ordered loop, and true. It falls back to the vertex average (and true)
// for a degenerate loop with nonzero vertex count, and returns false only for an
// empty loop. The computation triangulates the polygon as a fan from its first
// vertex.
func PolygonCentroid(loop []Vec3) (Vec3, bool) {
	m := len(loop)
	if m == 0 {
		return Vec3{}, false
	}
	if m < 3 {
		c, _ := Centroid(loop)
		return c, true
	}
	var (
		area float64
		acc  Vec3
	)
	a := loop[0]
	for i := 1; i+1 < m; i++ {
		b := loop[i]
		c := loop[i+1]
		w := TriangleArea(a, b, c)
		area += w
		acc = acc.Add(TriangleCentroid(a, b, c).Scale(w))
	}
	if area < Eps*Eps {
		c, _ := Centroid(loop)
		return c, true
	}
	return acc.Div(area), true
}

// PolygonPerimeter returns the total edge length of the closed polygon loop.
func PolygonPerimeter(loop []Vec3) float64 {
	m := len(loop)
	if m < 2 {
		return 0
	}
	var p float64
	for i := 0; i < m; i++ {
		p += loop[i].Distance(loop[(i+1)%m])
	}
	return p
}

// IsPlanar reports whether every vertex of the loop lies within tol of the plane
// through the first vertex with the polygon's Newell normal. Loops of three or
// fewer vertices are trivially planar.
func IsPlanar(loop []Vec3, tol float64) bool {
	if len(loop) <= 3 {
		return true
	}
	n, ok := PolygonUnitNormal(loop)
	if !ok {
		return false
	}
	base := loop[0]
	for _, p := range loop[1:] {
		if math.Abs(p.Sub(base).Dot(n)) > tol {
			return false
		}
	}
	return true
}

// RegularPolygonArea returns the area of a regular n-gon with side length a.
// It returns 0 for n < 3.
func RegularPolygonArea(n int, a float64) float64 {
	if n < 3 {
		return 0
	}
	return float64(n) * a * a / (4 * math.Tan(math.Pi/float64(n)))
}

// RegularPolygonCircumradius returns the circumradius of a regular n-gon with
// side length a (the distance from center to a vertex). It returns 0 for n < 3.
func RegularPolygonCircumradius(n int, a float64) float64 {
	if n < 3 {
		return 0
	}
	return a / (2 * math.Sin(math.Pi/float64(n)))
}

// RegularPolygonInradius returns the inradius (apothem) of a regular n-gon with
// side length a (the distance from center to an edge midpoint). It returns 0 for
// n < 3.
func RegularPolygonInradius(n int, a float64) float64 {
	if n < 3 {
		return 0
	}
	return a / (2 * math.Tan(math.Pi/float64(n)))
}

// RegularPolygonInteriorAngle returns the interior angle in radians of a regular
// n-gon. It returns 0 for n < 3.
func RegularPolygonInteriorAngle(n int) float64 {
	if n < 3 {
		return 0
	}
	return float64(n-2) * math.Pi / float64(n)
}
