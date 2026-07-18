package geom3d

import "math"

// Plane is an oriented plane in space, represented by the equation
// Normal·x = D. Points on the positive side of the plane (in the direction of
// Normal) have Normal·x > D. Normal is not required to be unit length, but the
// distance queries assume it is; call [Plane.Normalize] to enforce this.
type Plane struct {
	Normal Vec3
	D      float64
}

// PlaneFromPointNormal returns the plane passing through point p with the given
// normal. The normal is normalized so that distance queries are metric.
func PlaneFromPointNormal(p, normal Vec3) Plane {
	n := normal.Unit()
	return Plane{Normal: n, D: n.Dot(p)}
}

// PlaneFromPoints returns the plane through the three points a, b and c and
// true, or the zero plane and false if the points are collinear. The normal
// follows the right-hand rule from (b-a)×(c-a) and is unit length.
func PlaneFromPoints(a, b, c Vec3) (Plane, bool) {
	n, ln := b.Sub(a).Cross(c.Sub(a)).Normalize()
	if ln == 0 {
		return Plane{}, false
	}
	return Plane{Normal: n, D: n.Dot(a)}, true
}

// Normalize returns the plane rescaled so that Normal has unit length while
// describing the same set of points. If the normal has zero length the plane
// is returned unchanged.
func (pl Plane) Normalize() Plane {
	n, ln := pl.Normal.Normalize()
	if ln == 0 {
		return pl
	}
	return Plane{Normal: n, D: pl.D / ln}
}

// SignedDistance returns the signed distance from point p to the plane. It is
// positive on the side the normal points toward and negative on the other. The
// plane's normal is assumed to be unit length.
func (pl Plane) SignedDistance(p Vec3) float64 {
	return pl.Normal.Dot(p) - pl.D
}

// Distance returns the unsigned (absolute) distance from point p to the plane.
// The plane's normal is assumed to be unit length.
func (pl Plane) Distance(p Vec3) float64 {
	return math.Abs(pl.SignedDistance(p))
}

// ClosestPoint returns the orthogonal projection of point p onto the plane, the
// point on the plane nearest to p. The plane's normal is assumed to be unit
// length.
func (pl Plane) ClosestPoint(p Vec3) Vec3 {
	return p.Sub(pl.Normal.Scale(pl.SignedDistance(p)))
}

// Side reports which side of the plane point p lies on: +1 for the positive
// (normal) side, -1 for the negative side, and 0 if p lies on the plane within
// the absolute tolerance eps. The plane's normal is assumed to be unit length.
func (pl Plane) Side(p Vec3, eps float64) int {
	d := pl.SignedDistance(p)
	switch {
	case d > eps:
		return 1
	case d < -eps:
		return -1
	default:
		return 0
	}
}

// Contains reports whether point p lies on the plane within the absolute
// tolerance eps. The plane's normal is assumed to be unit length.
func (pl Plane) Contains(p Vec3, eps float64) bool {
	return math.Abs(pl.SignedDistance(p)) <= eps
}

// Line3 is the infinite straight line passing through Point in the direction
// Dir. Dir need not be unit length but must be non-zero for the geometric
// queries to be meaningful.
type Line3 struct {
	Point Vec3
	Dir   Vec3
}

// LineThrough returns the line passing through the two distinct points a and b,
// with direction b-a.
func LineThrough(a, b Vec3) Line3 {
	return Line3{Point: a, Dir: b.Sub(a)}
}

// At returns the point Point + t*Dir on the line.
func (l Line3) At(t float64) Vec3 {
	return l.Point.Add(l.Dir.Scale(t))
}

// Param returns the parameter t of the point on the line closest to p, so that
// [Line3.At](t) is the closest point. It returns 0 if the line's direction has
// zero length.
func (l Line3) Param(p Vec3) float64 {
	dd := l.Dir.LengthSq()
	if dd == 0 {
		return 0
	}
	return p.Sub(l.Point).Dot(l.Dir) / dd
}

// ClosestPoint returns the point on the line nearest to p.
func (l Line3) ClosestPoint(p Vec3) Vec3 {
	return l.At(l.Param(p))
}

// Distance returns the perpendicular distance from point p to the line. It
// returns the distance from p to Point if the line's direction has zero length.
func (l Line3) Distance(p Vec3) float64 {
	dd := l.Dir.LengthSq()
	if dd == 0 {
		return p.Distance(l.Point)
	}
	// |(p-Point) x Dir| / |Dir|
	return p.Sub(l.Point).Cross(l.Dir).Length() / math.Sqrt(dd)
}

// ClosestPointsBetweenLines returns the pair of points, one on line l1 and one
// on line l2, that are closest to each other, together with true. If the lines
// are parallel it returns an arbitrary closest pair (the foot of l1.Point on
// l2, paired with l1.Point) and false.
func ClosestPointsBetweenLines(l1, l2 Line3) (p1, p2 Vec3, ok bool) {
	d1 := l1.Dir
	d2 := l2.Dir
	r := l1.Point.Sub(l2.Point)
	a := d1.Dot(d1)
	b := d1.Dot(d2)
	c := d2.Dot(d2)
	e := d1.Dot(r)
	f := d2.Dot(r)
	denom := a*c - b*b
	if math.Abs(denom) <= geom3dEps {
		// Parallel lines: fix s=0, solve for t.
		s := 0.0
		var t float64
		if c != 0 {
			t = f / c
		}
		return l1.At(s), l2.At(t), false
	}
	s := (b*f - c*e) / denom
	t := (a*f - b*e) / denom
	return l1.At(s), l2.At(t), true
}

// Ray is a half-line starting at Origin and extending in the direction Dir.
// Dir need not be unit length; parameters returned by the intersection tests
// are expressed in units of Dir's length.
type Ray struct {
	Origin Vec3
	Dir    Vec3
}

// NewRay returns the ray from origin with the given direction. The direction is
// normalized so that intersection parameters are distances along the ray.
func NewRay(origin, dir Vec3) Ray {
	return Ray{Origin: origin, Dir: dir.Unit()}
}

// At returns the point Origin + t*Dir on the ray.
func (r Ray) At(t float64) Vec3 {
	return r.Origin.Add(r.Dir.Scale(t))
}

// IntersectPlane returns the parameter t at which the ray meets the plane and
// true, or 0 and false if the ray is parallel to the plane (or points away and
// never reaches it). A returned t is negative if the intersection lies behind
// the ray's origin. The plane's normal is assumed to be unit length.
func (r Ray) IntersectPlane(pl Plane) (t float64, ok bool) {
	denom := pl.Normal.Dot(r.Dir)
	if math.Abs(denom) <= geom3dEps {
		return 0, false
	}
	t = (pl.D - pl.Normal.Dot(r.Origin)) / denom
	return t, true
}
