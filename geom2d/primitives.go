// Package geom2d provides primitives for 2D computational geometry.
//
// It defines the value types [Point2] (a location in the plane) and [Vec2]
// (a displacement), together with the supporting [Line] and [Segment] types.
// On top of these it offers the standard toolbox of planar geometry: vector
// arithmetic (add, subtract, scale, dot, cross), norms and angles, rotation
// and reflection, projection, linear interpolation, point/point,
// point/line and point/segment distances, closest-point queries, orientation
// and collinearity predicates, and line and segment intersection.
//
// All computation is done in float64 using only the Go standard library.
// Routines are deterministic: given identical inputs they return identical
// results. Predicates that must tolerate floating-point noise accept or use an
// explicit epsilon so callers control the trade-off between robustness and
// strictness.
package geom2d

import "math"

// geom2dEps is the default absolute tolerance used by the approximate
// predicates in this package. It is comfortably above double-precision
// rounding noise for coordinates of moderate magnitude.
const geom2dEps = 1e-9

// Point2 is a point (location) in the Euclidean plane.
type Point2 struct {
	X, Y float64
}

// Vec2 is a vector (displacement or direction) in the Euclidean plane.
type Vec2 struct {
	X, Y float64
}

// Line is the infinite straight line passing through the two distinct points
// A and B. The direction of the line is B-A.
type Line struct {
	A, B Point2
}

// Segment is the finite straight segment with endpoints A and B.
type Segment struct {
	A, B Point2
}

// geom2dclamp returns x limited to the closed interval [lo, hi].
func geom2dclamp(x, lo, hi float64) float64 {
	if x < lo {
		return lo
	}
	if x > hi {
		return hi
	}
	return x
}

// geom2dsign returns -1, 0 or +1 according to the sign of x, treating any
// value with magnitude at or below eps as zero.
func geom2dsign(x, eps float64) int {
	if x > eps {
		return 1
	}
	if x < -eps {
		return -1
	}
	return 0
}

// geom2dhypot returns the Euclidean length sqrt(x*x+y*y) without spurious
// overflow or underflow.
func geom2dhypot(x, y float64) float64 {
	return math.Hypot(x, y)
}

// NewPoint2 returns the point (x, y).
func NewPoint2(x, y float64) Point2 { return Point2{X: x, Y: y} }

// NewVec2 returns the vector (x, y).
func NewVec2(x, y float64) Vec2 { return Vec2{X: x, Y: y} }

// PointFromPolar returns the point at distance r and angle theta (radians,
// measured counter-clockwise from the positive X axis) from the origin.
func PointFromPolar(r, theta float64) Point2 {
	return Point2{X: r * math.Cos(theta), Y: r * math.Sin(theta)}
}

// VecFromPolar returns the vector of magnitude r pointing at angle theta
// (radians, counter-clockwise from the positive X axis).
func VecFromPolar(r, theta float64) Vec2 {
	return Vec2{X: r * math.Cos(theta), Y: r * math.Sin(theta)}
}

// Vec returns the position vector of p (the displacement from the origin to p).
func (p Point2) Vec() Vec2 { return Vec2{X: p.X, Y: p.Y} }

// AsPoint returns v interpreted as the point at position v from the origin.
func (v Vec2) AsPoint() Point2 { return Point2{X: v.X, Y: v.Y} }

// Add returns the point p displaced by the vector v.
func (p Point2) Add(v Vec2) Point2 { return Point2{X: p.X + v.X, Y: p.Y + v.Y} }

// Sub returns the vector from q to p (that is, p-q).
func (p Point2) Sub(q Point2) Vec2 { return Vec2{X: p.X - q.X, Y: p.Y - q.Y} }

// Equal reports whether p and q are exactly equal in both coordinates.
func (p Point2) Equal(q Point2) bool { return p.X == q.X && p.Y == q.Y }

// ApproxEqual reports whether p and q are within eps (absolute) in each
// coordinate.
func (p Point2) ApproxEqual(q Point2, eps float64) bool {
	return math.Abs(p.X-q.X) <= eps && math.Abs(p.Y-q.Y) <= eps
}

// Add returns the vector sum u+v.
func (u Vec2) Add(v Vec2) Vec2 { return Vec2{X: u.X + v.X, Y: u.Y + v.Y} }

// Sub returns the vector difference u-v.
func (u Vec2) Sub(v Vec2) Vec2 { return Vec2{X: u.X - v.X, Y: u.Y - v.Y} }

// Scale returns the vector u multiplied by the scalar s.
func (u Vec2) Scale(s float64) Vec2 { return Vec2{X: u.X * s, Y: u.Y * s} }

// Neg returns the vector -u.
func (u Vec2) Neg() Vec2 { return Vec2{X: -u.X, Y: -u.Y} }

// Dot returns the dot (scalar) product u·v = u.X*v.X + u.Y*v.Y.
func (u Vec2) Dot(v Vec2) float64 { return u.X*v.X + u.Y*v.Y }

// Cross returns the 2D cross product u×v = u.X*v.Y - u.Y*v.X, the signed area
// of the parallelogram spanned by u and v. It is positive when v lies
// counter-clockwise from u.
func (u Vec2) Cross(v Vec2) float64 { return u.X*v.Y - u.Y*v.X }

// Hadamard returns the component-wise product (u.X*v.X, u.Y*v.Y).
func (u Vec2) Hadamard(v Vec2) Vec2 { return Vec2{X: u.X * v.X, Y: u.Y * v.Y} }

// Norm returns the Euclidean length |u| of the vector.
func (u Vec2) Norm() float64 { return geom2dhypot(u.X, u.Y) }

// NormSq returns the squared Euclidean length |u|², avoiding the square root.
func (u Vec2) NormSq() float64 { return u.X*u.X + u.Y*u.Y }

// IsZero reports whether u is the zero vector within eps: its length is at
// most eps.
func (u Vec2) IsZero(eps float64) bool { return u.NormSq() <= eps*eps }

// Normalize returns the unit vector in the direction of u and its original
// length. If u is the zero vector it returns the zero vector and length 0.
func (u Vec2) Normalize() (Vec2, float64) {
	n := u.Norm()
	if n == 0 {
		return Vec2{}, 0
	}
	return Vec2{X: u.X / n, Y: u.Y / n}, n
}

// Unit returns the unit vector in the direction of u, or the zero vector if u
// has zero length.
func (u Vec2) Unit() Vec2 {
	v, _ := u.Normalize()
	return v
}

// Perp returns u rotated 90° counter-clockwise: (-u.Y, u.X).
func (u Vec2) Perp() Vec2 { return Vec2{X: -u.Y, Y: u.X} }

// PerpCW returns u rotated 90° clockwise: (u.Y, -u.X).
func (u Vec2) PerpCW() Vec2 { return Vec2{X: u.Y, Y: -u.X} }

// Angle returns the direction of u in radians in (-π, π], measured
// counter-clockwise from the positive X axis. The zero vector yields 0.
func (u Vec2) Angle() float64 { return math.Atan2(u.Y, u.X) }

// Rotate returns u rotated by theta radians counter-clockwise about the
// origin.
func (u Vec2) Rotate(theta float64) Vec2 {
	s, c := math.Sincos(theta)
	return Vec2{X: u.X*c - u.Y*s, Y: u.X*s + u.Y*c}
}

// AngleBetween returns the unsigned angle in [0, π] between u and v. If either
// vector has zero length it returns 0.
func (u Vec2) AngleBetween(v Vec2) float64 {
	du, dv := u.Norm(), v.Norm()
	if du == 0 || dv == 0 {
		return 0
	}
	c := geom2dclamp(u.Dot(v)/(du*dv), -1, 1)
	return math.Acos(c)
}

// SignedAngle returns the signed angle in (-π, π] from u to v, positive when v
// is counter-clockwise from u.
func (u Vec2) SignedAngle(v Vec2) float64 {
	return math.Atan2(u.Cross(v), u.Dot(v))
}

// Project returns the vector projection of u onto v (the component of u
// parallel to v). If v is the zero vector it returns the zero vector.
func (u Vec2) Project(v Vec2) Vec2 {
	d := v.NormSq()
	if d == 0 {
		return Vec2{}
	}
	return v.Scale(u.Dot(v) / d)
}

// Reject returns the component of u orthogonal to v, that is u minus the
// projection of u onto v.
func (u Vec2) Reject(v Vec2) Vec2 { return u.Sub(u.Project(v)) }

// ScalarProject returns the signed length of the projection of u onto v.
// If v is the zero vector it returns 0.
func (u Vec2) ScalarProject(v Vec2) float64 {
	n := v.Norm()
	if n == 0 {
		return 0
	}
	return u.Dot(v) / n
}

// Reflect returns u reflected across the line through the origin with
// direction axis. If axis is the zero vector u is returned unchanged.
func (u Vec2) Reflect(axis Vec2) Vec2 {
	d := axis.NormSq()
	if d == 0 {
		return u
	}
	return u.Project(axis).Scale(2).Sub(u)
}

// ReflectAcrossNormal returns u reflected as if bouncing off a surface whose
// unit normal is n: u - 2(u·n)n. The normal should be a unit vector for a
// length-preserving reflection.
func (u Vec2) ReflectAcrossNormal(n Vec2) Vec2 {
	return u.Sub(n.Scale(2 * u.Dot(n)))
}

// Lerp returns the linear interpolation (1-t)*u + t*v.
func (u Vec2) Lerp(v Vec2, t float64) Vec2 {
	return Vec2{X: u.X + (v.X-u.X)*t, Y: u.Y + (v.Y-u.Y)*t}
}

// Equal reports whether u and v are exactly equal in both components.
func (u Vec2) Equal(v Vec2) bool { return u.X == v.X && u.Y == v.Y }

// ApproxEqual reports whether u and v are within eps (absolute) in each
// component.
func (u Vec2) ApproxEqual(v Vec2, eps float64) bool {
	return math.Abs(u.X-v.X) <= eps && math.Abs(u.Y-v.Y) <= eps
}

// Add returns the vector sum u+v.
func Add(u, v Vec2) Vec2 { return u.Add(v) }

// Sub returns the vector difference u-v.
func Sub(u, v Vec2) Vec2 { return u.Sub(v) }

// Scale returns the vector u multiplied by the scalar s.
func Scale(u Vec2, s float64) Vec2 { return u.Scale(s) }

// Dot returns the dot product u·v.
func Dot(u, v Vec2) float64 { return u.Dot(v) }

// Cross returns the 2D cross product u×v.
func Cross(u, v Vec2) float64 { return u.Cross(v) }

// Distance returns the Euclidean distance between points p and q.
func Distance(p, q Point2) float64 { return geom2dhypot(p.X-q.X, p.Y-q.Y) }

// DistanceSq returns the squared Euclidean distance between points p and q,
// avoiding the square root.
func DistanceSq(p, q Point2) float64 {
	dx, dy := p.X-q.X, p.Y-q.Y
	return dx*dx + dy*dy
}

// DistanceTo returns the Euclidean distance from p to q.
func (p Point2) DistanceTo(q Point2) float64 { return Distance(p, q) }

// DistanceSqTo returns the squared Euclidean distance from p to q.
func (p Point2) DistanceSqTo(q Point2) float64 { return DistanceSq(p, q) }

// Midpoint returns the point halfway between p and q.
func Midpoint(p, q Point2) Point2 {
	return Point2{X: (p.X + q.X) / 2, Y: (p.Y + q.Y) / 2}
}

// Lerp returns the point interpolated between p and q by parameter t:
// (1-t)*p + t*q. t=0 gives p and t=1 gives q.
func Lerp(p, q Point2, t float64) Point2 {
	return Point2{X: p.X + (q.X-p.X)*t, Y: p.Y + (q.Y-p.Y)*t}
}

// Centroid returns the arithmetic mean (centroid) of the given points. It
// returns the origin when no points are supplied.
func Centroid(pts ...Point2) Point2 {
	if len(pts) == 0 {
		return Point2{}
	}
	var sx, sy float64
	for _, p := range pts {
		sx += p.X
		sy += p.Y
	}
	n := float64(len(pts))
	return Point2{X: sx / n, Y: sy / n}
}

// RotateAround returns p rotated by theta radians counter-clockwise about the
// pivot point c.
func RotateAround(p, c Point2, theta float64) Point2 {
	return c.Add(p.Sub(c).Rotate(theta))
}

// SignedArea2 returns twice the signed area of triangle abc, equal to the
// cross product (b-a)×(c-a). It is positive when a, b, c are in
// counter-clockwise order and negative when clockwise.
func SignedArea2(a, b, c Point2) float64 {
	return b.Sub(a).Cross(c.Sub(a))
}

// TriangleArea returns the (non-negative) area of triangle abc.
func TriangleArea(a, b, c Point2) float64 {
	return math.Abs(SignedArea2(a, b, c)) / 2
}

// Orientation returns +1 if the points a, b, c make a counter-clockwise turn,
// -1 if they make a clockwise turn, and 0 if they are collinear, judged within
// eps applied to the signed area.
func Orientation(a, b, c Point2, eps float64) int {
	return geom2dsign(SignedArea2(a, b, c), eps)
}

// CCW reports whether the points a, b, c are in strict counter-clockwise
// order, using eps to reject near-collinear triples.
func CCW(a, b, c Point2, eps float64) bool {
	return SignedArea2(a, b, c) > eps
}

// Collinear reports whether the points a, b, c lie on a common line, judged by
// twice the triangle area being at most eps in magnitude.
func Collinear(a, b, c Point2, eps float64) bool {
	return math.Abs(SignedArea2(a, b, c)) <= eps
}

// AngleAt returns the interior angle in [0, π] at vertex b of the path a-b-c,
// that is the angle between the vectors b→a and b→c.
func AngleAt(a, b, c Point2) float64 {
	return a.Sub(b).AngleBetween(c.Sub(b))
}

// Slope returns the slope (dy/dx) of the line through p and q and whether it is
// finite. For a vertical line it returns (0, false).
func Slope(p, q Point2) (float64, bool) {
	dx := q.X - p.X
	if dx == 0 {
		return 0, false
	}
	return (q.Y - p.Y) / dx, true
}

// Vec returns the displacement vector from A to B along the line.
func (l Line) Vec() Vec2 { return l.B.Sub(l.A) }

// Direction returns the unit direction vector of the line (from A toward B).
func (l Line) Direction() Vec2 { return l.Vec().Unit() }

// PointAt returns the point A + t*(B-A) on the line, so t=0 gives A and t=1
// gives B.
func (l Line) PointAt(t float64) Point2 { return l.A.Add(l.Vec().Scale(t)) }

// Vec returns the displacement vector from A to B along the segment.
func (s Segment) Vec() Vec2 { return s.B.Sub(s.A) }

// Direction returns the unit direction vector of the segment (from A toward B).
func (s Segment) Direction() Vec2 { return s.Vec().Unit() }

// Length returns the Euclidean length of the segment.
func (s Segment) Length() float64 { return Distance(s.A, s.B) }

// LengthSq returns the squared Euclidean length of the segment.
func (s Segment) LengthSq() float64 { return DistanceSq(s.A, s.B) }

// Midpoint returns the midpoint of the segment.
func (s Segment) Midpoint() Point2 { return Midpoint(s.A, s.B) }

// PointAt returns the point A + t*(B-A); t is not clamped, so values outside
// [0, 1] lie on the segment's supporting line beyond its endpoints.
func (s Segment) PointAt(t float64) Point2 { return s.A.Add(s.Vec().Scale(t)) }

// AsLine returns the infinite line supporting the segment.
func (s Segment) AsLine() Line { return Line{A: s.A, B: s.B} }

// geom2dfootParam returns the parameter t such that A + t*(B-A) is the foot of
// the perpendicular from p to the line through a and b. For a segment, clamping
// t to [0, 1] yields the closest point on the segment. A degenerate
// (zero-length) direction yields t=0.
func geom2dfootParam(a, b, p Point2) float64 {
	d := b.Sub(a)
	dd := d.NormSq()
	if dd == 0 {
		return 0
	}
	return p.Sub(a).Dot(d) / dd
}

// ClosestPointOnLine returns the point on the infinite line l nearest to p
// (the foot of the perpendicular from p). If l is degenerate it returns l.A.
func ClosestPointOnLine(p Point2, l Line) Point2 {
	t := geom2dfootParam(l.A, l.B, p)
	return l.PointAt(t)
}

// ClosestPointOnSegment returns the point on segment s nearest to p, obtained
// by projecting p onto the segment's line and clamping to the endpoints.
func ClosestPointOnSegment(p Point2, s Segment) Point2 {
	t := geom2dclamp(geom2dfootParam(s.A, s.B, p), 0, 1)
	return s.PointAt(t)
}

// ClosestPoint returns the point on the line nearest to p.
func (l Line) ClosestPoint(p Point2) Point2 { return ClosestPointOnLine(p, l) }

// ClosestPoint returns the point on the segment nearest to p.
func (s Segment) ClosestPoint(p Point2) Point2 { return ClosestPointOnSegment(p, s) }

// PointLineDistance returns the perpendicular distance from p to the infinite
// line l. If l is degenerate it returns the distance from p to l.A.
func PointLineDistance(p Point2, l Line) float64 {
	d := l.Vec()
	dd := d.Norm()
	if dd == 0 {
		return Distance(p, l.A)
	}
	return math.Abs(d.Cross(p.Sub(l.A))) / dd
}

// SignedPointLineDistance returns the signed perpendicular distance from p to
// line l. The result is positive when p lies to the left of the direction A→B
// and negative to the right. A degenerate line yields the unsigned distance to
// l.A.
func SignedPointLineDistance(p Point2, l Line) float64 {
	d := l.Vec()
	dd := d.Norm()
	if dd == 0 {
		return Distance(p, l.A)
	}
	return d.Cross(p.Sub(l.A)) / dd
}

// PointSegmentDistance returns the distance from p to the closest point on
// segment s.
func PointSegmentDistance(p Point2, s Segment) float64 {
	return Distance(p, ClosestPointOnSegment(p, s))
}

// Distance returns the perpendicular distance from p to the line.
func (l Line) Distance(p Point2) float64 { return PointLineDistance(p, l) }

// Distance returns the distance from p to the closest point on the segment.
func (s Segment) Distance(p Point2) float64 { return PointSegmentDistance(p, s) }

// Contains reports whether p lies on the infinite line l within tolerance eps
// (measured as perpendicular distance).
func (l Line) Contains(p Point2, eps float64) bool {
	return PointLineDistance(p, l) <= eps
}

// Contains reports whether p lies on segment s within tolerance eps (measured
// as distance to the closest point on the segment).
func (s Segment) Contains(p Point2, eps float64) bool {
	return PointSegmentDistance(p, s) <= eps
}

// ProjectPointOntoLine returns the foot of the perpendicular from p to line l.
// It is an alias for [ClosestPointOnLine].
func ProjectPointOntoLine(p Point2, l Line) Point2 {
	return ClosestPointOnLine(p, l)
}

// ReflectPointAcrossLine returns p reflected across the infinite line l. The
// reflection maps p to the point on the far side of l at equal perpendicular
// distance.
func ReflectPointAcrossLine(p Point2, l Line) Point2 {
	foot := ClosestPointOnLine(p, l)
	return foot.Add(foot.Sub(p))
}

// PerpendicularBisector returns the line that passes through the midpoint of
// p and q at right angles to the segment pq. The returned line's A is the
// midpoint and its direction is perpendicular to q-p.
func PerpendicularBisector(p, q Point2) Line {
	m := Midpoint(p, q)
	return Line{A: m, B: m.Add(q.Sub(p).Perp())}
}

// OnSegment reports whether point p, already assumed to be collinear with the
// segment's endpoints within eps, lies within the segment's extent. Use
// [Segment.Contains] for a full membership test.
func OnSegment(p Point2, s Segment, eps float64) bool {
	if !Collinear(s.A, s.B, p, eps) {
		return false
	}
	minX, maxX := math.Min(s.A.X, s.B.X), math.Max(s.A.X, s.B.X)
	minY, maxY := math.Min(s.A.Y, s.B.Y), math.Max(s.A.Y, s.B.Y)
	return p.X >= minX-eps && p.X <= maxX+eps && p.Y >= minY-eps && p.Y <= maxY+eps
}

// LineIntersection returns the intersection point of the infinite lines l1 and
// l2 and reports whether it exists. Parallel (or coincident) lines return
// (zero, false); eps bounds the cross product below which the lines are treated
// as parallel.
func LineIntersection(l1, l2 Line, eps float64) (Point2, bool) {
	r := l1.Vec()
	s := l2.Vec()
	denom := r.Cross(s)
	if math.Abs(denom) <= eps {
		return Point2{}, false
	}
	t := l2.A.Sub(l1.A).Cross(s) / denom
	return l1.PointAt(t), true
}

// SegmentsIntersect reports whether segments s1 and s2 share at least one
// point. It handles the collinear-overlap and shared-endpoint cases, using eps
// for the orientation and on-segment tests.
func SegmentsIntersect(s1, s2 Segment, eps float64) bool {
	p1, q1 := s1.A, s1.B
	p2, q2 := s2.A, s2.B
	o1 := Orientation(p1, q1, p2, eps)
	o2 := Orientation(p1, q1, q2, eps)
	o3 := Orientation(p2, q2, p1, eps)
	o4 := Orientation(p2, q2, q1, eps)
	if o1 != o2 && o3 != o4 {
		return true
	}
	if o1 == 0 && OnSegment(p2, s1, eps) {
		return true
	}
	if o2 == 0 && OnSegment(q2, s1, eps) {
		return true
	}
	if o3 == 0 && OnSegment(p1, s2, eps) {
		return true
	}
	if o4 == 0 && OnSegment(q1, s2, eps) {
		return true
	}
	return false
}

// SegmentIntersection returns the unique intersection point of segments s1 and
// s2 and reports whether such a single point exists. Non-intersecting,
// parallel, or collinearly overlapping segments return (zero, false); use
// [SegmentsIntersect] to detect overlap without a unique point.
func SegmentIntersection(s1, s2 Segment, eps float64) (Point2, bool) {
	p, r := s1.A, s1.Vec()
	q, s := s2.A, s2.Vec()
	denom := r.Cross(s)
	if math.Abs(denom) <= eps {
		return Point2{}, false
	}
	qp := q.Sub(p)
	t := qp.Cross(s) / denom
	u := qp.Cross(r) / denom
	if t < -eps || t > 1+eps || u < -eps || u > 1+eps {
		return Point2{}, false
	}
	return s1.PointAt(geom2dclamp(t, 0, 1)), true
}

// LineSegmentIntersection returns the intersection point of infinite line l
// with segment s and reports whether it exists. A parallel line, or a crossing
// that falls outside the segment, returns (zero, false).
func LineSegmentIntersection(l Line, s Segment, eps float64) (Point2, bool) {
	r := l.Vec()
	sd := s.Vec()
	denom := r.Cross(sd)
	if math.Abs(denom) <= eps {
		return Point2{}, false
	}
	u := s.A.Sub(l.A).Cross(r) / denom
	if u < -eps || u > 1+eps {
		return Point2{}, false
	}
	return s.PointAt(geom2dclamp(u, 0, 1)), true
}

// SegmentSegmentDistance returns the minimum distance between segments s1 and
// s2, which is zero when they intersect.
func SegmentSegmentDistance(s1, s2 Segment, eps float64) float64 {
	if SegmentsIntersect(s1, s2, eps) {
		return 0
	}
	d := PointSegmentDistance(s1.A, s2)
	d = math.Min(d, PointSegmentDistance(s1.B, s2))
	d = math.Min(d, PointSegmentDistance(s2.A, s1))
	d = math.Min(d, PointSegmentDistance(s2.B, s1))
	return d
}
