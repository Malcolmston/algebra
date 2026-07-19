package voronoi

import (
	"errors"
	"math"
)

// ErrCollinear is returned by circumcircle constructions when the three input
// points are collinear and no unique circle exists.
var ErrCollinear = errors.New("voronoi: points are collinear")

// Orient2D returns twice the signed area of the triangle (a, b, c). The result
// is positive when a, b, c are in counterclockwise order, negative when
// clockwise, and zero when collinear. This is the standard planar orientation
// determinant.
func Orient2D(a, b, c Point) float64 {
	return (b.X-a.X)*(c.Y-a.Y) - (b.Y-a.Y)*(c.X-a.X)
}

// Orientation classifies the turn a -> b -> c as +1 (counterclockwise),
// -1 (clockwise) or 0 (collinear) using the tolerance eps on the signed area.
func Orientation(a, b, c Point, eps float64) int {
	d := Orient2D(a, b, c)
	switch {
	case d > eps:
		return 1
	case d < -eps:
		return -1
	default:
		return 0
	}
}

// IsCCW reports whether the triangle (a, b, c) is oriented counterclockwise.
func IsCCW(a, b, c Point) bool { return Orient2D(a, b, c) > 0 }

// IsCW reports whether the triangle (a, b, c) is oriented clockwise.
func IsCW(a, b, c Point) bool { return Orient2D(a, b, c) < 0 }

// Collinear reports whether a, b and c are collinear within tolerance eps
// (applied to twice the signed area).
func Collinear(a, b, c Point, eps float64) bool {
	return math.Abs(Orient2D(a, b, c)) <= eps
}

// InCircle returns a positive value when d lies inside the circle through
// a, b, c (assumed counterclockwise), negative when outside, and zero when the
// four points are concyclic. It evaluates the classic 4x4 in-circle
// determinant reduced to a 3x3 form.
func InCircle(a, b, c, d Point) float64 {
	ax := a.X - d.X
	ay := a.Y - d.Y
	bx := b.X - d.X
	by := b.Y - d.Y
	cx := c.X - d.X
	cy := c.Y - d.Y

	az := ax*ax + ay*ay
	bz := bx*bx + by*by
	cz := cx*cx + cy*cy

	return az*(bx*cy-cx*by) -
		bz*(ax*cy-cx*ay) +
		cz*(ax*by-bx*ay)
}

// InCircleTest reports whether d lies strictly inside the circumcircle of the
// triangle (a, b, c). The orientation of a, b, c is handled automatically.
func InCircleTest(a, b, c, d Point) bool {
	v := InCircle(a, b, c, d)
	if Orient2D(a, b, c) < 0 {
		v = -v
	}
	return v > 0
}

// Circumcenter returns the centre of the unique circle passing through a, b
// and c. It returns ErrCollinear if the points are collinear.
func Circumcenter(a, b, c Point) (Point, error) {
	d := 2 * (a.X*(b.Y-c.Y) + b.X*(c.Y-a.Y) + c.X*(a.Y-b.Y))
	if d == 0 {
		return Point{}, ErrCollinear
	}
	a2 := a.NormSq()
	b2 := b.NormSq()
	c2 := c.NormSq()
	ux := (a2*(b.Y-c.Y) + b2*(c.Y-a.Y) + c2*(a.Y-b.Y)) / d
	uy := (a2*(c.X-b.X) + b2*(a.X-c.X) + c2*(b.X-a.X)) / d
	return Point{ux, uy}, nil
}

// Circumradius returns the radius of the circumcircle of the triangle
// (a, b, c). It returns ErrCollinear if the points are collinear.
func Circumradius(a, b, c Point) (float64, error) {
	center, err := Circumcenter(a, b, c)
	if err != nil {
		return 0, err
	}
	return center.Distance(a), nil
}

// Circumcircle returns the circumscribed circle of the triangle (a, b, c).
// It returns ErrCollinear if the points are collinear.
func Circumcircle(a, b, c Point) (Circle, error) {
	center, err := Circumcenter(a, b, c)
	if err != nil {
		return Circle{}, err
	}
	return Circle{Center: center, Radius: center.Distance(a)}, nil
}

// TriangleSignedArea returns the signed area of the triangle (a, b, c),
// positive for counterclockwise orientation.
func TriangleSignedArea(a, b, c Point) float64 {
	return Orient2D(a, b, c) / 2
}

// TriangleArea returns the (non-negative) area of the triangle (a, b, c).
func TriangleArea(a, b, c Point) float64 {
	return math.Abs(Orient2D(a, b, c)) / 2
}

// TrianglePerimeter returns the perimeter of the triangle (a, b, c).
func TrianglePerimeter(a, b, c Point) float64 {
	return a.Distance(b) + b.Distance(c) + c.Distance(a)
}

// TriangleCentroid returns the centroid (barycentre) of the triangle (a, b, c).
func TriangleCentroid(a, b, c Point) Point {
	return Point{(a.X + b.X + c.X) / 3, (a.Y + b.Y + c.Y) / 3}
}

// Incenter returns the incentre of the triangle (a, b, c): the centre of its
// inscribed circle, weighted by the opposite side lengths.
func Incenter(a, b, c Point) Point {
	la := b.Distance(c)
	lb := a.Distance(c)
	lc := a.Distance(b)
	s := la + lb + lc
	if s == 0 {
		return a
	}
	return Point{
		(la*a.X + lb*b.X + lc*c.X) / s,
		(la*a.Y + lb*b.Y + lc*c.Y) / s,
	}
}

// Inradius returns the radius of the inscribed circle of the triangle
// (a, b, c).
func Inradius(a, b, c Point) float64 {
	s := TrianglePerimeter(a, b, c) / 2
	if s == 0 {
		return 0
	}
	return TriangleArea(a, b, c) / s
}

// PointInTriangle reports whether p lies inside or on the boundary of the
// triangle (a, b, c), using tolerance eps on the orientation sign.
func PointInTriangle(p, a, b, c Point, eps float64) bool {
	d1 := Orient2D(a, b, p)
	d2 := Orient2D(b, c, p)
	d3 := Orient2D(c, a, p)
	hasNeg := d1 < -eps || d2 < -eps || d3 < -eps
	hasPos := d1 > eps || d2 > eps || d3 > eps
	return !(hasNeg && hasPos)
}

// Barycentric returns the barycentric coordinates (u, v, w) of p with respect
// to the triangle (a, b, c), so that p = u*a + v*b + w*c and u+v+w = 1. The
// boolean is false when the triangle is degenerate.
func Barycentric(p, a, b, c Point) (u, v, w float64, ok bool) {
	den := Orient2D(a, b, c)
	if den == 0 {
		return 0, 0, 0, false
	}
	u = Orient2D(p, b, c) / den
	v = Orient2D(a, p, c) / den
	w = Orient2D(a, b, p) / den
	return u, v, w, true
}

// SegmentsIntersect reports whether the closed segments p1p2 and q1q2 share at
// least one point. Collinear overlapping segments count as intersecting.
func SegmentsIntersect(p1, p2, q1, q2 Point) bool {
	d1 := Orient2D(q1, q2, p1)
	d2 := Orient2D(q1, q2, p2)
	d3 := Orient2D(p1, p2, q1)
	d4 := Orient2D(p1, p2, q2)
	if ((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) &&
		((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0)) {
		return true
	}
	if d1 == 0 && onSegment(q1, q2, p1) {
		return true
	}
	if d2 == 0 && onSegment(q1, q2, p2) {
		return true
	}
	if d3 == 0 && onSegment(p1, p2, q1) {
		return true
	}
	if d4 == 0 && onSegment(p1, p2, q2) {
		return true
	}
	return false
}

// onSegment reports whether r lies on the segment pq given that p, q, r are
// known to be collinear.
func onSegment(p, q, r Point) bool {
	return math.Min(p.X, q.X) <= r.X && r.X <= math.Max(p.X, q.X) &&
		math.Min(p.Y, q.Y) <= r.Y && r.Y <= math.Max(p.Y, q.Y)
}
