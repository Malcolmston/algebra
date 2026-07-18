package geom2d

import (
	"math"
	"sort"
)

// Polygon is a planar polygon described by an ordered list of vertices. The
// polygon is treated as closed: the edge from the last vertex back to the
// first is implicit and must not be repeated in the vertex slice. Vertices may
// be ordered clockwise or counter-clockwise; routines that care about the
// winding direction document their expectations individually.
type Polygon []Point2

// Circle is the set of points at exactly Radius distance from Center. A
// non-negative Radius is expected; a zero radius degenerates to a single point.
type Circle struct {
	Center Point2
	Radius float64
}

// Triangle is the triangle with vertices A, B and C.
type Triangle struct {
	A, B, C Point2
}

// BoundingBox is an axis-aligned rectangle spanning from Min (lower-left) to
// Max (upper-right). A valid box satisfies Min.X <= Max.X and Min.Y <= Max.Y.
type BoundingBox struct {
	Min, Max Point2
}

// geom2dcross returns the z-component of the cross product of vectors (a-o) and
// (b-o). It is positive when o->a->b makes a counter-clockwise (left) turn,
// negative for a clockwise (right) turn and zero when the three points are
// collinear.
func geom2dcross(o, a, b Point2) float64 {
	return (a.X-o.X)*(b.Y-o.Y) - (a.Y-o.Y)*(b.X-o.X)
}

// geom2dcross2 returns the cross product of the edge vectors (a1-a0) and
// (b1-b0). It is used by the rotating-calipers routine to compare edge
// directions.
func geom2dcross2(a0, a1, b0, b1 Point2) float64 {
	return (a1.X-a0.X)*(b1.Y-b0.Y) - (a1.Y-a0.Y)*(b1.X-b0.X)
}

// geom2dpointLess orders points lexicographically by X then Y. It provides the
// total order required by the monotone-chain convex-hull construction.
func geom2dpointLess(p, q Point2) bool {
	if p.X != q.X {
		return p.X < q.X
	}
	return p.Y < q.Y
}

// NewCircle returns the circle centred at c with the given radius.
func NewCircle(c Point2, radius float64) Circle {
	return Circle{Center: c, Radius: radius}
}

// NewTriangle returns the triangle with the given vertices.
func NewTriangle(a, b, c Point2) Triangle {
	return Triangle{A: a, B: b, C: c}
}

// NewBoundingBox returns the axis-aligned box spanning the two opposite corners
// p and q. The corners may be supplied in any order; the result is normalised
// so that Min holds the componentwise minimum and Max the componentwise
// maximum.
func NewBoundingBox(p, q Point2) BoundingBox {
	return BoundingBox{
		Min: Point2{X: math.Min(p.X, q.X), Y: math.Min(p.Y, q.Y)},
		Max: Point2{X: math.Max(p.X, q.X), Y: math.Max(p.Y, q.Y)},
	}
}

// ---------------------------------------------------------------------------
// Polygon: area, perimeter, centroid
// ---------------------------------------------------------------------------

// PolygonSignedArea returns the signed area of the closed polygon defined by
// pts using the shoelace (Gauss) formula. The result is positive when the
// vertices are ordered counter-clockwise and negative when they are ordered
// clockwise. Fewer than three vertices yield zero.
func PolygonSignedArea(pts []Point2) float64 {
	n := len(pts)
	if n < 3 {
		return 0
	}
	var sum float64
	for i := 0; i < n; i++ {
		j := i + 1
		if j == n {
			j = 0
		}
		sum += pts[i].X*pts[j].Y - pts[j].X*pts[i].Y
	}
	return sum / 2
}

// PolygonArea returns the non-negative area enclosed by the polygon pts,
// independent of vertex winding order.
func PolygonArea(pts []Point2) float64 {
	return math.Abs(PolygonSignedArea(pts))
}

// PolygonPerimeter returns the total edge length of the closed polygon pts,
// including the closing edge from the last vertex back to the first. Fewer than
// two vertices yield zero.
func PolygonPerimeter(pts []Point2) float64 {
	n := len(pts)
	if n < 2 {
		return 0
	}
	var per float64
	for i := 0; i < n; i++ {
		j := i + 1
		if j == n {
			j = 0
		}
		per += geom2dhypot(pts[j].X-pts[i].X, pts[j].Y-pts[i].Y)
	}
	return per
}

// PolygonCentroid returns the area centroid (centre of mass of a uniform
// lamina) of the closed polygon pts. For a degenerate polygon whose signed area
// is (near) zero it falls back to the arithmetic mean of the vertices. The
// second result reports whether the true area-weighted centroid was used.
func PolygonCentroid(pts []Point2) (Point2, bool) {
	n := len(pts)
	if n == 0 {
		return Point2{}, false
	}
	a := PolygonSignedArea(pts)
	if math.Abs(a) <= geom2dEps {
		var sx, sy float64
		for _, p := range pts {
			sx += p.X
			sy += p.Y
		}
		return Point2{X: sx / float64(n), Y: sy / float64(n)}, false
	}
	var cx, cy float64
	for i := 0; i < n; i++ {
		j := i + 1
		if j == n {
			j = 0
		}
		cross := pts[i].X*pts[j].Y - pts[j].X*pts[i].Y
		cx += (pts[i].X + pts[j].X) * cross
		cy += (pts[i].Y + pts[j].Y) * cross
	}
	f := 1.0 / (6.0 * a)
	return Point2{X: cx * f, Y: cy * f}, true
}

// SignedArea returns the signed shoelace area of the polygon; positive for
// counter-clockwise winding, negative for clockwise.
func (poly Polygon) SignedArea() float64 { return PolygonSignedArea(poly) }

// Area returns the non-negative area enclosed by the polygon.
func (poly Polygon) Area() float64 { return PolygonArea(poly) }

// Perimeter returns the total boundary length of the polygon.
func (poly Polygon) Perimeter() float64 { return PolygonPerimeter(poly) }

// Centroid returns the area centroid of the polygon. The boolean reports
// whether an area-weighted centroid was computed (false when the polygon is
// degenerate and the vertex average was returned instead).
func (poly Polygon) Centroid() (Point2, bool) { return PolygonCentroid(poly) }

// NumVertices returns the number of vertices in the polygon.
func (poly Polygon) NumVertices() int { return len(poly) }

// IsCounterClockwise reports whether the polygon vertices are ordered
// counter-clockwise (positive signed area).
func (poly Polygon) IsCounterClockwise() bool { return poly.SignedArea() > 0 }

// IsClockwise reports whether the polygon vertices are ordered clockwise
// (negative signed area).
func (poly Polygon) IsClockwise() bool { return poly.SignedArea() < 0 }

// Reverse returns a new polygon with the vertex order reversed, flipping the
// winding direction. The receiver is not modified.
func (poly Polygon) Reverse() Polygon {
	n := len(poly)
	out := make(Polygon, n)
	for i := 0; i < n; i++ {
		out[i] = poly[n-1-i]
	}
	return out
}

// EnsureCCW returns a polygon with the same shape guaranteed to be wound
// counter-clockwise, reversing the receiver only when necessary. The result may
// alias the receiver when it is already counter-clockwise.
func (poly Polygon) EnsureCCW() Polygon {
	if poly.SignedArea() < 0 {
		return poly.Reverse()
	}
	return poly
}

// EnsureCW returns a polygon with the same shape guaranteed to be wound
// clockwise, reversing the receiver only when necessary.
func (poly Polygon) EnsureCW() Polygon {
	if poly.SignedArea() > 0 {
		return poly.Reverse()
	}
	return poly
}

// Edges returns the polygon's boundary edges as segments, including the closing
// edge from the last vertex back to the first. A polygon with fewer than two
// vertices yields no edges.
func (poly Polygon) Edges() []Segment {
	n := len(poly)
	if n < 2 {
		return nil
	}
	edges := make([]Segment, n)
	for i := 0; i < n; i++ {
		j := i + 1
		if j == n {
			j = 0
		}
		edges[i] = Segment{A: poly[i], B: poly[j]}
	}
	return edges
}

// Translate returns a copy of the polygon with every vertex displaced by v.
func (poly Polygon) Translate(v Vec2) Polygon {
	out := make(Polygon, len(poly))
	for i, p := range poly {
		out[i] = Point2{X: p.X + v.X, Y: p.Y + v.Y}
	}
	return out
}

// ScaleAbout returns a copy of the polygon scaled by factor about the fixed
// point center. A factor greater than one enlarges the polygon; a negative
// factor also reflects it through center.
func (poly Polygon) ScaleAbout(factor float64, center Point2) Polygon {
	out := make(Polygon, len(poly))
	for i, p := range poly {
		out[i] = Point2{
			X: center.X + (p.X-center.X)*factor,
			Y: center.Y + (p.Y-center.Y)*factor,
		}
	}
	return out
}

// RotateAbout returns a copy of the polygon rotated by theta radians (positive
// counter-clockwise) about the fixed point center.
func (poly Polygon) RotateAbout(theta float64, center Point2) Polygon {
	s, c := math.Sincos(theta)
	out := make(Polygon, len(poly))
	for i, p := range poly {
		dx := p.X - center.X
		dy := p.Y - center.Y
		out[i] = Point2{
			X: center.X + dx*c - dy*s,
			Y: center.Y + dx*s + dy*c,
		}
	}
	return out
}

// ---------------------------------------------------------------------------
// Polygon: point containment and convexity
// ---------------------------------------------------------------------------

// PointInPolygon reports whether point p lies inside the closed polygon pts
// using the ray-casting (even-odd / crossing-number) rule. Points exactly on
// the boundary are reported inconsistently by this rule and should be tested
// separately with OnPolygonBoundary if that matters. The polygon may be
// non-convex but must be simple (non-self-intersecting) for a meaningful
// result.
func PointInPolygon(p Point2, pts []Point2) bool {
	n := len(pts)
	if n < 3 {
		return false
	}
	inside := false
	j := n - 1
	for i := 0; i < n; i++ {
		yi, yj := pts[i].Y, pts[j].Y
		if (yi > p.Y) != (yj > p.Y) {
			xint := (pts[j].X-pts[i].X)*(p.Y-yi)/(yj-yi) + pts[i].X
			if p.X < xint {
				inside = !inside
			}
		}
		j = i
	}
	return inside
}

// WindingNumber returns the winding number of the closed polygon pts around
// point p: the signed number of times the polygon boundary travels
// counter-clockwise around p. A non-zero winding number means p is enclosed.
// This is robust for self-intersecting polygons where the even-odd rule is
// ambiguous.
func WindingNumber(p Point2, pts []Point2) int {
	n := len(pts)
	if n < 3 {
		return 0
	}
	wn := 0
	for i := 0; i < n; i++ {
		j := i + 1
		if j == n {
			j = 0
		}
		if pts[i].Y <= p.Y {
			if pts[j].Y > p.Y && geom2dcross(pts[i], pts[j], p) > 0 {
				wn++
			}
		} else {
			if pts[j].Y <= p.Y && geom2dcross(pts[i], pts[j], p) < 0 {
				wn--
			}
		}
	}
	return wn
}

// PointInPolygonWinding reports whether p lies inside the closed polygon pts
// using the non-zero winding rule (the winding number is non-zero).
func PointInPolygonWinding(p Point2, pts []Point2) bool {
	return WindingNumber(p, pts) != 0
}

// OnPolygonBoundary reports whether p lies within eps of any edge of the closed
// polygon pts.
func OnPolygonBoundary(p Point2, pts []Point2, eps float64) bool {
	n := len(pts)
	if n < 2 {
		return false
	}
	for i := 0; i < n; i++ {
		j := i + 1
		if j == n {
			j = 0
		}
		if PointSegmentDistance(p, Segment{A: pts[i], B: pts[j]}) <= eps {
			return true
		}
	}
	return false
}

// IsConvexPolygon reports whether the closed polygon pts is convex. All turns
// along the boundary must have the same orientation (collinear vertices are
// permitted). Fewer than three vertices are considered non-convex.
func IsConvexPolygon(pts []Point2) bool {
	n := len(pts)
	if n < 3 {
		return false
	}
	var sign int
	for i := 0; i < n; i++ {
		a := pts[i]
		b := pts[(i+1)%n]
		c := pts[(i+2)%n]
		cr := geom2dcross(a, b, c)
		s := geom2dsign(cr, geom2dEps)
		if s == 0 {
			continue
		}
		if sign == 0 {
			sign = s
		} else if s != sign {
			return false
		}
	}
	return true
}

// ContainsPoint reports whether p lies inside the polygon using the ray-casting
// rule.
func (poly Polygon) ContainsPoint(p Point2) bool { return PointInPolygon(p, poly) }

// ContainsPointWinding reports whether p lies inside the polygon using the
// non-zero winding rule.
func (poly Polygon) ContainsPointWinding(p Point2) bool { return PointInPolygonWinding(p, poly) }

// WindingNumber returns the winding number of the polygon around p.
func (poly Polygon) WindingNumber(p Point2) int { return WindingNumber(p, poly) }

// OnBoundary reports whether p lies within eps of the polygon boundary.
func (poly Polygon) OnBoundary(p Point2, eps float64) bool { return OnPolygonBoundary(p, poly, eps) }

// IsConvex reports whether the polygon is convex.
func (poly Polygon) IsConvex() bool { return IsConvexPolygon(poly) }

// ---------------------------------------------------------------------------
// Bounding box
// ---------------------------------------------------------------------------

// BoundingBoxOf returns the smallest axis-aligned box containing all the given
// points. The second result is false when no points are supplied.
func BoundingBoxOf(pts ...Point2) (BoundingBox, bool) {
	if len(pts) == 0 {
		return BoundingBox{}, false
	}
	minX, minY := pts[0].X, pts[0].Y
	maxX, maxY := pts[0].X, pts[0].Y
	for _, p := range pts[1:] {
		minX = math.Min(minX, p.X)
		minY = math.Min(minY, p.Y)
		maxX = math.Max(maxX, p.X)
		maxY = math.Max(maxY, p.Y)
	}
	return BoundingBox{Min: Point2{X: minX, Y: minY}, Max: Point2{X: maxX, Y: maxY}}, true
}

// PolygonBoundingBox returns the axis-aligned bounding box of the polygon pts.
// The second result is false when pts is empty.
func PolygonBoundingBox(pts []Point2) (BoundingBox, bool) {
	return BoundingBoxOf(pts...)
}

// BoundingBox returns the axis-aligned bounding box of the polygon. The second
// result is false when the polygon has no vertices.
func (poly Polygon) BoundingBox() (BoundingBox, bool) { return BoundingBoxOf(poly...) }

// Width returns the horizontal extent of the box.
func (b BoundingBox) Width() float64 { return b.Max.X - b.Min.X }

// Height returns the vertical extent of the box.
func (b BoundingBox) Height() float64 { return b.Max.Y - b.Min.Y }

// Area returns the area of the box.
func (b BoundingBox) Area() float64 { return b.Width() * b.Height() }

// Perimeter returns the perimeter of the box.
func (b BoundingBox) Perimeter() float64 { return 2 * (b.Width() + b.Height()) }

// Center returns the centre point of the box.
func (b BoundingBox) Center() Point2 {
	return Point2{X: (b.Min.X + b.Max.X) / 2, Y: (b.Min.Y + b.Max.Y) / 2}
}

// ContainsPoint reports whether p lies within the closed box (boundary
// included).
func (b BoundingBox) ContainsPoint(p Point2) bool {
	return p.X >= b.Min.X && p.X <= b.Max.X && p.Y >= b.Min.Y && p.Y <= b.Max.Y
}

// Intersects reports whether the two boxes overlap (sharing an edge counts as
// intersecting).
func (b BoundingBox) Intersects(o BoundingBox) bool {
	return b.Min.X <= o.Max.X && b.Max.X >= o.Min.X &&
		b.Min.Y <= o.Max.Y && b.Max.Y >= o.Min.Y
}

// Union returns the smallest box containing both b and o.
func (b BoundingBox) Union(o BoundingBox) BoundingBox {
	return BoundingBox{
		Min: Point2{X: math.Min(b.Min.X, o.Min.X), Y: math.Min(b.Min.Y, o.Min.Y)},
		Max: Point2{X: math.Max(b.Max.X, o.Max.X), Y: math.Max(b.Max.Y, o.Max.Y)},
	}
}

// Intersection returns the overlapping box of b and o. The second result is
// false when the two boxes do not intersect.
func (b BoundingBox) Intersection(o BoundingBox) (BoundingBox, bool) {
	if !b.Intersects(o) {
		return BoundingBox{}, false
	}
	return BoundingBox{
		Min: Point2{X: math.Max(b.Min.X, o.Min.X), Y: math.Max(b.Min.Y, o.Min.Y)},
		Max: Point2{X: math.Min(b.Max.X, o.Max.X), Y: math.Min(b.Max.Y, o.Max.Y)},
	}, true
}

// Corners returns the four corners of the box in counter-clockwise order
// starting from Min (lower-left).
func (b BoundingBox) Corners() [4]Point2 {
	return [4]Point2{
		{X: b.Min.X, Y: b.Min.Y},
		{X: b.Max.X, Y: b.Min.Y},
		{X: b.Max.X, Y: b.Max.Y},
		{X: b.Min.X, Y: b.Max.Y},
	}
}

// ---------------------------------------------------------------------------
// Convex hull (Andrew's monotone chain)
// ---------------------------------------------------------------------------

// ConvexHull returns the convex hull of the input points as a polygon wound
// counter-clockwise, using Andrew's monotone-chain algorithm in O(n log n)
// time. Interior and collinear boundary points are omitted; the returned hull
// contains only true vertices. Duplicate input points are tolerated. Fewer than
// three unique points return the (de-duplicated, sorted) input unchanged.
func ConvexHull(pts []Point2) Polygon {
	n := len(pts)
	if n < 3 {
		out := make(Polygon, n)
		copy(out, pts)
		sort.Slice(out, func(i, j int) bool { return geom2dpointLess(out[i], out[j]) })
		return out
	}
	sorted := make([]Point2, n)
	copy(sorted, pts)
	sort.Slice(sorted, func(i, j int) bool { return geom2dpointLess(sorted[i], sorted[j]) })

	hull := make([]Point2, 0, 2*n)
	// Lower hull.
	for _, p := range sorted {
		for len(hull) >= 2 && geom2dcross(hull[len(hull)-2], hull[len(hull)-1], p) <= 0 {
			hull = hull[:len(hull)-1]
		}
		hull = append(hull, p)
	}
	// Upper hull.
	lower := len(hull) + 1
	for i := n - 2; i >= 0; i-- {
		p := sorted[i]
		for len(hull) >= lower && geom2dcross(hull[len(hull)-2], hull[len(hull)-1], p) <= 0 {
			hull = hull[:len(hull)-1]
		}
		hull = append(hull, p)
	}
	// Drop the last point because it repeats the first.
	return Polygon(hull[:len(hull)-1])
}

// ---------------------------------------------------------------------------
// Sutherland-Hodgman polygon clipping
// ---------------------------------------------------------------------------

// geom2dinsideHalfPlane reports whether p lies on the inside (left) side of the
// directed line a->b, boundary inclusive.
func geom2dinsideHalfPlane(p, a, b Point2) bool {
	return geom2dcross(a, b, p) >= -geom2dEps
}

// geom2dlineIntersectParam returns the intersection point of the segment p->q
// with the infinite line through a and b. The caller guarantees the segment
// crosses the line.
func geom2dlineIntersectParam(p, q, a, b Point2) Point2 {
	// Solve for t where p + t(q-p) lies on line a->b.
	rx, ry := q.X-p.X, q.Y-p.Y
	// Normal of line a->b is (-(b.Y-a.Y), b.X-a.X).
	nx, ny := -(b.Y - a.Y), b.X-a.X
	denom := rx*nx + ry*ny
	if denom == 0 {
		return p
	}
	t := ((a.X-p.X)*nx + (a.Y-p.Y)*ny) / denom
	return Point2{X: p.X + t*rx, Y: p.Y + t*ry}
}

// ClipPolygonHalfPlane clips the subject polygon against the half-plane lying
// to the left of the directed line through a and b (the inside, for a
// counter-clockwise convex boundary). It returns the portion of the subject on
// the inside as a new polygon, which may be empty.
func ClipPolygonHalfPlane(subject []Point2, a, b Point2) Polygon {
	n := len(subject)
	if n == 0 {
		return Polygon{}
	}
	out := make(Polygon, 0, n+1)
	prev := subject[n-1]
	prevIn := geom2dinsideHalfPlane(prev, a, b)
	for _, cur := range subject {
		curIn := geom2dinsideHalfPlane(cur, a, b)
		if curIn {
			if !prevIn {
				out = append(out, geom2dlineIntersectParam(prev, cur, a, b))
			}
			out = append(out, cur)
		} else if prevIn {
			out = append(out, geom2dlineIntersectParam(prev, cur, a, b))
		}
		prev = cur
		prevIn = curIn
	}
	return out
}

// ClipPolygon clips the subject polygon against a convex clip polygon using the
// Sutherland-Hodgman algorithm and returns the intersection as a new polygon.
// The clip polygon must be convex and wound counter-clockwise; the subject may
// be any simple polygon. The result may be empty when the polygons do not
// overlap.
func ClipPolygon(subject, clip []Point2) Polygon {
	if len(clip) < 3 {
		out := make(Polygon, len(subject))
		copy(out, subject)
		return out
	}
	result := make(Polygon, len(subject))
	copy(result, subject)
	m := len(clip)
	for i := 0; i < m; i++ {
		a := clip[i]
		b := clip[(i+1)%m]
		result = ClipPolygonHalfPlane(result, a, b)
		if len(result) == 0 {
			break
		}
	}
	return result
}

// ---------------------------------------------------------------------------
// Rotating calipers: convex diameter
// ---------------------------------------------------------------------------

// ConvexDiameterPoints returns the two points of the input set that are
// farthest apart (the diameter of the point set) together with that distance,
// computed via a convex hull followed by the rotating-calipers method in
// O(n log n) time. When fewer than two points are supplied the two results are
// equal and the distance is zero.
func ConvexDiameterPoints(pts []Point2) (Point2, Point2, float64) {
	h := ConvexHull(pts)
	n := len(h)
	switch n {
	case 0:
		return Point2{}, Point2{}, 0
	case 1:
		return h[0], h[0], 0
	case 2:
		return h[0], h[1], Distance(h[0], h[1])
	}
	best := -1.0
	var bi, bj int
	j := 1
	for i := 0; i < n; i++ {
		ni := (i + 1) % n
		for {
			nj := (j + 1) % n
			if geom2dcross2(h[i], h[ni], h[j], h[nj]) > geom2dEps {
				j = nj
			} else {
				break
			}
		}
		if d := DistanceSq(h[i], h[j]); d > best {
			best, bi, bj = d, i, j
		}
		if d := DistanceSq(h[ni], h[j]); d > best {
			best, bi, bj = d, ni, j
		}
	}
	return h[bi], h[bj], math.Sqrt(best)
}

// ConvexDiameter returns the diameter (maximum pairwise distance) of the input
// point set, computed with rotating calipers over the convex hull.
func ConvexDiameter(pts []Point2) float64 {
	_, _, d := ConvexDiameterPoints(pts)
	return d
}

// ---------------------------------------------------------------------------
// Triangle
// ---------------------------------------------------------------------------

// SignedArea returns the signed area of the triangle; positive when A, B, C are
// in counter-clockwise order, negative when clockwise.
func (t Triangle) SignedArea() float64 {
	return ((t.B.X-t.A.X)*(t.C.Y-t.A.Y) - (t.C.X-t.A.X)*(t.B.Y-t.A.Y)) / 2
}

// Area returns the non-negative area of the triangle.
func (t Triangle) Area() float64 { return math.Abs(t.SignedArea()) }

// Perimeter returns the sum of the triangle's three side lengths.
func (t Triangle) Perimeter() float64 {
	return Distance(t.A, t.B) + Distance(t.B, t.C) + Distance(t.C, t.A)
}

// Centroid returns the triangle's centroid (the average of its vertices, which
// coincides with the centre of mass).
func (t Triangle) Centroid() Point2 {
	return Point2{X: (t.A.X + t.B.X + t.C.X) / 3, Y: (t.A.Y + t.B.Y + t.C.Y) / 3}
}

// BoundingBox returns the axis-aligned bounding box of the triangle.
func (t Triangle) BoundingBox() BoundingBox {
	b, _ := BoundingBoxOf(t.A, t.B, t.C)
	return b
}

// Barycentric returns the barycentric coordinates (u, v, w) of point p with
// respect to the triangle, satisfying p = u*A + v*B + w*C and u+v+w = 1. The
// final boolean is false for a degenerate (zero-area) triangle, for which the
// coordinates are undefined and returned as zero.
func (t Triangle) Barycentric(p Point2) (u, v, w float64, ok bool) {
	v0 := t.B.Sub(t.A)
	v1 := t.C.Sub(t.A)
	v2 := p.Sub(t.A)
	den := v0.X*v1.Y - v1.X*v0.Y
	if math.Abs(den) <= geom2dEps {
		return 0, 0, 0, false
	}
	v = (v2.X*v1.Y - v1.X*v2.Y) / den
	w = (v0.X*v2.Y - v2.X*v0.Y) / den
	u = 1 - v - w
	return u, v, w, true
}

// PointFromBarycentric returns the Cartesian point with barycentric coordinates
// (u, v, w) relative to the triangle. The coordinates are used as supplied and
// are not required to sum to one.
func (t Triangle) PointFromBarycentric(u, v, w float64) Point2 {
	return Point2{
		X: u*t.A.X + v*t.B.X + w*t.C.X,
		Y: u*t.A.Y + v*t.B.Y + w*t.C.Y,
	}
}

// ContainsPoint reports whether p lies inside or on the boundary of the
// triangle, using barycentric coordinates. A degenerate triangle contains no
// points.
func (t Triangle) ContainsPoint(p Point2) bool {
	u, v, w, ok := t.Barycentric(p)
	if !ok {
		return false
	}
	return u >= -geom2dEps && v >= -geom2dEps && w >= -geom2dEps
}

// Circumcenter returns the centre of the triangle's circumscribed circle (the
// point equidistant from all three vertices). The second result is false for a
// degenerate triangle whose vertices are collinear.
func (t Triangle) Circumcenter() (Point2, bool) {
	ax, ay := t.A.X, t.A.Y
	bx, by := t.B.X, t.B.Y
	cx, cy := t.C.X, t.C.Y
	d := 2 * (ax*(by-cy) + bx*(cy-ay) + cx*(ay-by))
	if math.Abs(d) <= geom2dEps {
		return Point2{}, false
	}
	a2 := ax*ax + ay*ay
	b2 := bx*bx + by*by
	c2 := cx*cx + cy*cy
	ux := (a2*(by-cy) + b2*(cy-ay) + c2*(ay-by)) / d
	uy := (a2*(cx-bx) + b2*(ax-cx) + c2*(bx-ax)) / d
	return Point2{X: ux, Y: uy}, true
}

// Circumradius returns the radius of the triangle's circumscribed circle. The
// second result is false for a degenerate triangle.
func (t Triangle) Circumradius() (float64, bool) {
	a := Distance(t.B, t.C)
	b := Distance(t.A, t.C)
	c := Distance(t.A, t.B)
	area := t.Area()
	if area <= geom2dEps {
		return 0, false
	}
	return a * b * c / (4 * area), true
}

// Circumcircle returns the circle passing through all three vertices of the
// triangle. The second result is false for a degenerate triangle.
func (t Triangle) Circumcircle() (Circle, bool) {
	center, ok := t.Circumcenter()
	if !ok {
		return Circle{}, false
	}
	return Circle{Center: center, Radius: Distance(center, t.A)}, true
}

// Incenter returns the centre of the triangle's inscribed circle, the
// intersection of the angle bisectors. The second result is false for a
// degenerate triangle whose perimeter is (near) zero.
func (t Triangle) Incenter() (Point2, bool) {
	a := Distance(t.B, t.C) // side opposite A
	b := Distance(t.A, t.C) // side opposite B
	c := Distance(t.A, t.B) // side opposite C
	per := a + b + c
	if per <= geom2dEps {
		return Point2{}, false
	}
	return Point2{
		X: (a*t.A.X + b*t.B.X + c*t.C.X) / per,
		Y: (a*t.A.Y + b*t.B.Y + c*t.C.Y) / per,
	}, true
}

// Inradius returns the radius of the triangle's inscribed circle, equal to the
// area divided by the semiperimeter. The second result is false for a
// degenerate triangle.
func (t Triangle) Inradius() (float64, bool) {
	per := t.Perimeter()
	if per <= geom2dEps {
		return 0, false
	}
	return 2 * t.Area() / per, true
}

// Incircle returns the largest circle contained in the triangle, tangent to all
// three sides. The second result is false for a degenerate triangle.
func (t Triangle) Incircle() (Circle, bool) {
	center, ok := t.Incenter()
	if !ok {
		return Circle{}, false
	}
	r, ok := t.Inradius()
	if !ok {
		return Circle{}, false
	}
	return Circle{Center: center, Radius: r}, true
}

// Orthocenter returns the point where the triangle's three altitudes meet. The
// second result is false for a degenerate triangle.
func (t Triangle) Orthocenter() (Point2, bool) {
	// Orthocenter H = A + B + C - 2*O, where O is the circumcenter.
	o, ok := t.Circumcenter()
	if !ok {
		return Point2{}, false
	}
	return Point2{
		X: t.A.X + t.B.X + t.C.X - 2*o.X,
		Y: t.A.Y + t.B.Y + t.C.Y - 2*o.Y,
	}, true
}

// TriangleCircumcircle returns the circle through points a, b and c. The second
// result is false when the points are collinear.
func TriangleCircumcircle(a, b, c Point2) (Circle, bool) {
	return Triangle{A: a, B: b, C: c}.Circumcircle()
}

// TriangleIncircle returns the inscribed circle of the triangle with vertices
// a, b and c. The second result is false for a degenerate triangle.
func TriangleIncircle(a, b, c Point2) (Circle, bool) {
	return Triangle{A: a, B: b, C: c}.Incircle()
}

// Barycentric returns the barycentric coordinates of point p with respect to
// the triangle a, b, c. See Triangle.Barycentric for details.
func Barycentric(p, a, b, c Point2) (u, v, w float64, ok bool) {
	return Triangle{A: a, B: b, C: c}.Barycentric(p)
}

// ---------------------------------------------------------------------------
// Circle
// ---------------------------------------------------------------------------

// Area returns the area enclosed by the circle.
func (c Circle) Area() float64 { return math.Pi * c.Radius * c.Radius }

// Circumference returns the perimeter (circumference) of the circle.
func (c Circle) Circumference() float64 { return 2 * math.Pi * c.Radius }

// Diameter returns the diameter of the circle.
func (c Circle) Diameter() float64 { return 2 * c.Radius }

// ContainsPoint reports whether p lies inside or on the circle.
func (c Circle) ContainsPoint(p Point2) bool {
	return DistanceSq(p, c.Center) <= c.Radius*c.Radius+geom2dEps
}

// OnBoundary reports whether p lies within eps of the circle boundary.
func (c Circle) OnBoundary(p Point2, eps float64) bool {
	return math.Abs(Distance(p, c.Center)-c.Radius) <= eps
}

// BoundingBox returns the axis-aligned bounding box that tightly encloses the
// circle.
func (c Circle) BoundingBox() BoundingBox {
	return BoundingBox{
		Min: Point2{X: c.Center.X - c.Radius, Y: c.Center.Y - c.Radius},
		Max: Point2{X: c.Center.X + c.Radius, Y: c.Center.Y + c.Radius},
	}
}

// PointAtAngle returns the point on the circle at the given angle in radians,
// measured counter-clockwise from the positive x-axis.
func (c Circle) PointAtAngle(theta float64) Point2 {
	s, co := math.Sincos(theta)
	return Point2{X: c.Center.X + c.Radius*co, Y: c.Center.Y + c.Radius*s}
}

// CircleFrom3Points returns the unique circle passing through three points. The
// second result is false when the points are collinear (no finite circle
// exists).
func CircleFrom3Points(a, b, c Point2) (Circle, bool) {
	return TriangleCircumcircle(a, b, c)
}

// CircleLineIntersection returns the intersection points of circle c with the
// infinite line l. It returns zero, one (tangent) or two points, ordered along
// the direction of the line.
func CircleLineIntersection(c Circle, l Line) []Point2 {
	d := l.B.Sub(l.A)
	dd := d.NormSq()
	if dd == 0 {
		return nil
	}
	// Parameter of the foot of perpendicular from the centre.
	f := c.Center.Sub(l.A)
	t := (f.X*d.X + f.Y*d.Y) / dd
	foot := Point2{X: l.A.X + t*d.X, Y: l.A.Y + t*d.Y}
	dist2 := DistanceSq(foot, c.Center)
	r2 := c.Radius * c.Radius
	if dist2 > r2+geom2dEps {
		return nil
	}
	if dist2 >= r2-geom2dEps {
		return []Point2{foot}
	}
	h := math.Sqrt((r2 - dist2) / dd)
	p1 := Point2{X: foot.X - h*d.X, Y: foot.Y - h*d.Y}
	p2 := Point2{X: foot.X + h*d.X, Y: foot.Y + h*d.Y}
	return []Point2{p1, p2}
}

// CircleSegmentIntersection returns the intersection points of circle c with
// the finite segment s, restricted to points that actually lie on the segment.
func CircleSegmentIntersection(c Circle, s Segment) []Point2 {
	pts := CircleLineIntersection(c, Line{A: s.A, B: s.B})
	if len(pts) == 0 {
		return nil
	}
	out := make([]Point2, 0, len(pts))
	for _, p := range pts {
		if OnSegment(p, s, 1e-7) {
			out = append(out, p)
		}
	}
	return out
}

// CircleCircleIntersection returns the intersection points of two circles. It
// returns zero points when the circles are separate or one strictly contains
// the other, one point when they are tangent, and two points when they cross.
// Coincident circles (infinitely many intersections) return zero points with a
// false second result; all other cases return true.
func CircleCircleIntersection(c1, c2 Circle) ([]Point2, bool) {
	dx := c2.Center.X - c1.Center.X
	dy := c2.Center.Y - c1.Center.Y
	d := geom2dhypot(dx, dy)
	if d <= geom2dEps && math.Abs(c1.Radius-c2.Radius) <= geom2dEps {
		// Coincident circles.
		return nil, false
	}
	if d > c1.Radius+c2.Radius+geom2dEps {
		return nil, true // separate
	}
	if d < math.Abs(c1.Radius-c2.Radius)-geom2dEps {
		return nil, true // one inside the other
	}
	if d <= geom2dEps {
		return nil, true
	}
	a := (c1.Radius*c1.Radius - c2.Radius*c2.Radius + d*d) / (2 * d)
	h2 := c1.Radius*c1.Radius - a*a
	if h2 < 0 {
		h2 = 0
	}
	h := math.Sqrt(h2)
	mid := Point2{X: c1.Center.X + a*dx/d, Y: c1.Center.Y + a*dy/d}
	if h <= geom2dEps {
		return []Point2{mid}, true
	}
	ox := -h * dy / d
	oy := h * dx / d
	return []Point2{
		{X: mid.X + ox, Y: mid.Y + oy},
		{X: mid.X - ox, Y: mid.Y - oy},
	}, true
}
