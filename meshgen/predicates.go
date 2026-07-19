package meshgen

import (
	"errors"
	"math"
)

// ErrCollinear is returned by circumcircle constructions when the three input
// points are collinear and no unique circle exists.
var ErrCollinear = errors.New("meshgen: points are collinear")

// ErrDegenerate is returned when an input configuration is degenerate for the
// requested operation.
var ErrDegenerate = errors.New("meshgen: degenerate input")

// Orient2D returns twice the signed area of the triangle (a, b, c). It is
// positive when a, b, c are counterclockwise, negative when clockwise, and
// zero when the three points are collinear.
func Orient2D(a, b, c Vec2) float64 {
	return (b.X-a.X)*(c.Y-a.Y) - (b.Y-a.Y)*(c.X-a.X)
}

// Orientation classifies the turn a -> b -> c as +1 (counterclockwise),
// -1 (clockwise) or 0 (collinear) using the tolerance eps on the signed area.
func Orientation(a, b, c Vec2, eps float64) int {
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
func IsCCW(a, b, c Vec2) bool { return Orient2D(a, b, c) > 0 }

// IsCW reports whether the triangle (a, b, c) is oriented clockwise.
func IsCW(a, b, c Vec2) bool { return Orient2D(a, b, c) < 0 }

// Collinear reports whether a, b and c are collinear within tolerance eps
// applied to twice the signed area.
func Collinear(a, b, c Vec2, eps float64) bool {
	return math.Abs(Orient2D(a, b, c)) <= eps
}

// SignedArea2 returns the signed area of the triangle (a, b, c); positive for
// counterclockwise orientation.
func SignedArea2(a, b, c Vec2) float64 { return Orient2D(a, b, c) / 2 }

// TriangleArea returns the non-negative area of the triangle (a, b, c).
func TriangleArea(a, b, c Vec2) float64 { return math.Abs(Orient2D(a, b, c)) / 2 }

// InCircle returns a positive value when d lies inside the circle through
// a, b, c (assumed counterclockwise), negative when outside, and zero when the
// four points are concyclic.
func InCircle(a, b, c, d Vec2) float64 {
	ax := a.X - d.X
	ay := a.Y - d.Y
	bx := b.X - d.X
	by := b.Y - d.Y
	cx := c.X - d.X
	cy := c.Y - d.Y

	az := ax*ax + ay*ay
	bz := bx*bx + by*by
	cz := cx*cx + cy*cy

	return ax*(by*cz-bz*cy) - ay*(bx*cz-bz*cx) + az*(bx*cy-by*cx)
}

// InCircleCCW returns the in-circle determinant with the first three points
// re-ordered so they are counterclockwise, so that a positive result always
// means d is strictly inside the circle regardless of the input winding.
func InCircleCCW(a, b, c, d Vec2) float64 {
	if Orient2D(a, b, c) < 0 {
		a, b = b, a
	}
	return InCircle(a, b, c, d)
}

// Circumcenter returns the center of the circle passing through a, b and c.
// It returns ErrCollinear when the points are collinear.
func Circumcenter(a, b, c Vec2) (Vec2, error) {
	d := 2 * (a.X*(b.Y-c.Y) + b.X*(c.Y-a.Y) + c.X*(a.Y-b.Y))
	if math.Abs(d) < 1e-300 {
		return Vec2{}, ErrCollinear
	}
	a2 := a.X*a.X + a.Y*a.Y
	b2 := b.X*b.X + b.Y*b.Y
	c2 := c.X*c.X + c.Y*c.Y
	ux := (a2*(b.Y-c.Y) + b2*(c.Y-a.Y) + c2*(a.Y-b.Y)) / d
	uy := (a2*(c.X-b.X) + b2*(a.X-c.X) + c2*(b.X-a.X)) / d
	return Vec2{ux, uy}, nil
}

// Circumradius returns the radius of the circle passing through a, b and c.
// It returns ErrCollinear when the points are collinear.
func Circumradius(a, b, c Vec2) (float64, error) {
	ce, err := Circumcenter(a, b, c)
	if err != nil {
		return 0, err
	}
	return ce.Distance(a), nil
}

// Circumcircle returns the center and radius of the circle through a, b and c.
func Circumcircle(a, b, c Vec2) (center Vec2, radius float64, err error) {
	center, err = Circumcenter(a, b, c)
	if err != nil {
		return Vec2{}, 0, err
	}
	return center, center.Distance(a), nil
}

// Incenter returns the incenter of the triangle (a, b, c), the center of its
// inscribed circle.
func Incenter(a, b, c Vec2) Vec2 {
	la := b.Distance(c)
	lb := a.Distance(c)
	lc := a.Distance(b)
	s := la + lb + lc
	if s == 0 {
		return a
	}
	return a.Scale(la).Add(b.Scale(lb)).Add(c.Scale(lc)).Div(s)
}

// Inradius returns the radius of the inscribed circle of the triangle (a,b,c).
func Inradius(a, b, c Vec2) float64 {
	area := TriangleArea(a, b, c)
	s := (a.Distance(b) + b.Distance(c) + c.Distance(a)) / 2
	if s == 0 {
		return 0
	}
	return area / s
}

// Barycentric returns the barycentric coordinates (u, v, w) of p with respect
// to the triangle (a, b, c), so that p = u*a + v*b + w*c and u+v+w = 1. The
// boolean is false when the triangle is degenerate.
func Barycentric(p, a, b, c Vec2) (u, v, w float64, ok bool) {
	den := Orient2D(a, b, c)
	if den == 0 {
		return 0, 0, 0, false
	}
	u = Orient2D(p, b, c) / den
	v = Orient2D(a, p, c) / den
	w = Orient2D(a, b, p) / den
	return u, v, w, true
}

// PointInTriangle reports whether p lies inside or on the boundary of the
// triangle (a, b, c), using tolerance eps.
func PointInTriangle(p, a, b, c Vec2, eps float64) bool {
	u, v, w, ok := Barycentric(p, a, b, c)
	if !ok {
		return false
	}
	return u >= -eps && v >= -eps && w >= -eps
}

// ClosestPointOnSegment returns the point on the segment a-b closest to p.
func ClosestPointOnSegment(p, a, b Vec2) Vec2 {
	ab := b.Sub(a)
	den := ab.NormSq()
	if den == 0 {
		return a
	}
	t := p.Sub(a).Dot(ab) / den
	if t < 0 {
		t = 0
	} else if t > 1 {
		t = 1
	}
	return a.Add(ab.Scale(t))
}

// DistancePointSegment returns the distance from p to the segment a-b.
func DistancePointSegment(p, a, b Vec2) float64 {
	return p.Distance(ClosestPointOnSegment(p, a, b))
}

// SegmentsProperlyIntersect reports whether the open segments a-b and c-d cross
// at a single interior point. Shared endpoints and collinear overlaps return
// false.
func SegmentsProperlyIntersect(a, b, c, d Vec2) bool {
	d1 := Orient2D(c, d, a)
	d2 := Orient2D(c, d, b)
	d3 := Orient2D(a, b, c)
	d4 := Orient2D(a, b, d)
	return ((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) &&
		((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0))
}

// SegmentIntersection returns the intersection point of the lines through a-b
// and c-d. The boolean is false when the lines are parallel.
func SegmentIntersection(a, b, c, d Vec2) (Vec2, bool) {
	r := b.Sub(a)
	s := d.Sub(c)
	den := r.Cross(s)
	if den == 0 {
		return Vec2{}, false
	}
	t := c.Sub(a).Cross(s) / den
	return a.Add(r.Scale(t)), true
}

// PolygonSignedArea returns the signed area of the polygon given by its ordered
// vertices; positive for counterclockwise orientation.
func PolygonSignedArea(poly []Vec2) float64 {
	n := len(poly)
	if n < 3 {
		return 0
	}
	var s float64
	for i := 0; i < n; i++ {
		j := (i + 1) % n
		s += poly[i].Cross(poly[j])
	}
	return s / 2
}

// PolygonArea returns the non-negative area of the polygon.
func PolygonArea(poly []Vec2) float64 { return math.Abs(PolygonSignedArea(poly)) }

// PolygonIsCCW reports whether the polygon is oriented counterclockwise.
func PolygonIsCCW(poly []Vec2) bool { return PolygonSignedArea(poly) > 0 }

// PolygonPerimeter returns the total edge length of the closed polygon.
func PolygonPerimeter(poly []Vec2) float64 {
	n := len(poly)
	if n < 2 {
		return 0
	}
	var s float64
	for i := 0; i < n; i++ {
		s += poly[i].Distance(poly[(i+1)%n])
	}
	return s
}

// PolygonCentroid returns the area centroid of a simple polygon. It falls back
// to the vertex average for a degenerate (zero-area) polygon.
func PolygonCentroid(poly []Vec2) Vec2 {
	n := len(poly)
	if n == 0 {
		return Vec2{}
	}
	a := PolygonSignedArea(poly)
	if a == 0 {
		return CentroidVec2(poly)
	}
	var cx, cy float64
	for i := 0; i < n; i++ {
		j := (i + 1) % n
		cr := poly[i].Cross(poly[j])
		cx += (poly[i].X + poly[j].X) * cr
		cy += (poly[i].Y + poly[j].Y) * cr
	}
	return Vec2{cx / (6 * a), cy / (6 * a)}
}

// PointInPolygon reports whether p lies inside the simple polygon using the
// even-odd (ray casting) rule. Points exactly on an edge may return either
// result.
func PointInPolygon(p Vec2, poly []Vec2) bool {
	n := len(poly)
	if n < 3 {
		return false
	}
	inside := false
	j := n - 1
	for i := 0; i < n; i++ {
		pi, pj := poly[i], poly[j]
		if (pi.Y > p.Y) != (pj.Y > p.Y) {
			x := (pj.X-pi.X)*(p.Y-pi.Y)/(pj.Y-pi.Y) + pi.X
			if p.X < x {
				inside = !inside
			}
		}
		j = i
	}
	return inside
}

// Orient3D returns six times the signed volume of the tetrahedron (a,b,c,d).
// It is positive when d is below the plane of a, b, c as seen with a, b, c
// counterclockwise.
func Orient3D(a, b, c, d Vec3) float64 {
	return a.Sub(d).Dot(b.Sub(d).Cross(c.Sub(d)))
}

// TetraVolume returns the non-negative volume of the tetrahedron (a,b,c,d).
func TetraVolume(a, b, c, d Vec3) float64 { return math.Abs(Orient3D(a, b, c, d)) / 6 }

// TriangleArea3 returns the area of the triangle (a, b, c) in space.
func TriangleArea3(a, b, c Vec3) float64 {
	return b.Sub(a).Cross(c.Sub(a)).Norm() / 2
}

// TriangleNormal3 returns the unit normal of the triangle (a, b, c), following
// the right-hand rule for the given vertex order.
func TriangleNormal3(a, b, c Vec3) Vec3 {
	return b.Sub(a).Cross(c.Sub(a)).Normalize()
}
