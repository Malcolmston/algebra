package voronoi

import "math"

// PolygonSignedArea returns the signed area of the simple polygon given by its
// vertices in order (via the shoelace formula). It is positive for
// counterclockwise winding and negative for clockwise.
func PolygonSignedArea(poly []Point) float64 {
	n := len(poly)
	if n < 3 {
		return 0
	}
	var s float64
	for i := 0; i < n; i++ {
		j := (i + 1) % n
		s += poly[i].X*poly[j].Y - poly[j].X*poly[i].Y
	}
	return s / 2
}

// PolygonArea returns the (non-negative) area of the simple polygon.
func PolygonArea(poly []Point) float64 {
	return math.Abs(PolygonSignedArea(poly))
}

// PolygonPerimeter returns the perimeter of the polygon given by its vertices
// in order (treated as closed).
func PolygonPerimeter(poly []Point) float64 {
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

// PolygonCentroid returns the area centroid of the simple polygon. For a
// degenerate (zero-area) polygon it falls back to the vertex average.
func PolygonCentroid(poly []Point) Point {
	n := len(poly)
	if n == 0 {
		return Point{}
	}
	a := PolygonSignedArea(poly)
	if a == 0 {
		return Centroid(poly)
	}
	var cx, cy float64
	for i := 0; i < n; i++ {
		j := (i + 1) % n
		cross := poly[i].X*poly[j].Y - poly[j].X*poly[i].Y
		cx += (poly[i].X + poly[j].X) * cross
		cy += (poly[i].Y + poly[j].Y) * cross
	}
	f := 1 / (6 * a)
	return Point{cx * f, cy * f}
}

// PolygonIsCCW reports whether the polygon's vertices wind counterclockwise.
func PolygonIsCCW(poly []Point) bool {
	return PolygonSignedArea(poly) > 0
}

// PointInPolygon reports whether p lies inside the simple polygon, using the
// even-odd (ray casting) rule. Points on the boundary may be classified either
// way depending on rounding.
func PointInPolygon(p Point, poly []Point) bool {
	n := len(poly)
	if n < 3 {
		return false
	}
	inside := false
	j := n - 1
	for i := 0; i < n; i++ {
		pi := poly[i]
		pj := poly[j]
		if (pi.Y > p.Y) != (pj.Y > p.Y) {
			xint := (pj.X-pi.X)*(p.Y-pi.Y)/(pj.Y-pi.Y) + pi.X
			if p.X < xint {
				inside = !inside
			}
		}
		j = i
	}
	return inside
}

// PointInConvexPolygon reports whether p lies inside or on the convex polygon
// whose vertices are given in counterclockwise order, within tolerance eps.
func PointInConvexPolygon(p Point, poly []Point, eps float64) bool {
	n := len(poly)
	if n < 3 {
		return false
	}
	for i := 0; i < n; i++ {
		if Orient2D(poly[i], poly[(i+1)%n], p) < -eps {
			return false
		}
	}
	return true
}

// PolygonWindingNumber returns the winding number of the polygon around p. A
// nonzero value means p is enclosed.
func PolygonWindingNumber(p Point, poly []Point) int {
	n := len(poly)
	wn := 0
	for i := 0; i < n; i++ {
		a := poly[i]
		b := poly[(i+1)%n]
		if a.Y <= p.Y {
			if b.Y > p.Y && Orient2D(a, b, p) > 0 {
				wn++
			}
		} else {
			if b.Y <= p.Y && Orient2D(a, b, p) < 0 {
				wn--
			}
		}
	}
	return wn
}

// PolygonBoundingBox returns the axis-aligned bounding box of the polygon
// vertices. It returns ErrEmpty for an empty polygon.
func PolygonBoundingBox(poly []Point) (Rect, error) {
	return BoundingBox(poly)
}
