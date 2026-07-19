package hull3d

import (
	"errors"
	"fmt"
	"math"
)

// Plane represents an oriented plane in R^3 as the set of points p satisfying
// Normal·p = Offset. The Normal need not be unit length; the positive
// half-space is {p : Normal·p > Offset}.
type Plane struct {
	Normal Vec3
	Offset float64
}

// NewPlane returns the plane with the given normal and offset.
func NewPlane(normal Vec3, offset float64) Plane { return Plane{normal, offset} }

// PlaneFromPointNormal returns the plane through point p with the given normal.
func PlaneFromPointNormal(p, normal Vec3) Plane {
	return Plane{normal, normal.Dot(p)}
}

// PlaneFromPoints returns the plane through the three points a, b, c, oriented
// so that its normal is (b-a)×(c-a). It returns an error if the points are
// collinear (degenerate).
func PlaneFromPoints(a, b, c Vec3) (Plane, error) {
	n := b.Sub(a).Cross(c.Sub(a))
	if n.LengthSq() < 1e-300 {
		return Plane{}, errors.New("hull3d: collinear points do not define a plane")
	}
	return Plane{n, n.Dot(a)}, nil
}

// Normalized returns an equivalent plane whose normal has unit length, so that
// SignedDistance reports true Euclidean distance. It returns an error if the
// normal is the zero vector.
func (p Plane) Normalized() (Plane, error) {
	n := p.Normal.Length()
	if n < 1e-300 {
		return Plane{}, errors.New("hull3d: degenerate plane with zero normal")
	}
	return Plane{p.Normal.Scale(1 / n), p.Offset / n}, nil
}

// Eval returns Normal·q - Offset, which is positive for q in the positive
// half-space, zero on the plane, and negative below.
func (p Plane) Eval(q Vec3) float64 { return p.Normal.Dot(q) - p.Offset }

// SignedDistance returns the signed distance from q to the plane, positive on
// the side the normal points to. The magnitude is exact only when the normal is
// unit length (see [Plane.Normalized]).
func (p Plane) SignedDistance(q Vec3) float64 {
	n := p.Normal.Length()
	if n < 1e-300 {
		return 0
	}
	return p.Eval(q) / n
}

// Distance returns the unsigned Euclidean distance from q to the plane.
func (p Plane) Distance(q Vec3) float64 { return math.Abs(p.SignedDistance(q)) }

// Side returns +1 if q lies in the positive half-space beyond tolerance eps, -1
// if it lies in the negative half-space, and 0 if it is within eps of the plane.
func (p Plane) Side(q Vec3, eps float64) int {
	d := p.SignedDistance(q)
	if d > eps {
		return 1
	}
	if d < -eps {
		return -1
	}
	return 0
}

// Contains reports whether q lies on the plane to within tolerance eps.
func (p Plane) Contains(q Vec3, eps float64) bool { return p.Side(q, eps) == 0 }

// Project returns the orthogonal projection of q onto the plane. It returns an
// error if the plane is degenerate.
func (p Plane) Project(q Vec3) (Vec3, error) {
	pn, err := p.Normalized()
	if err != nil {
		return Vec3{}, err
	}
	return q.Sub(pn.Normal.Scale(pn.Eval(q))), nil
}

// Flip returns the plane describing the same point set with the opposite
// orientation (negated normal and offset).
func (p Plane) Flip() Plane { return Plane{p.Normal.Neg(), -p.Offset} }

// PointOnPlane returns some point lying on the plane (the foot of the
// perpendicular from the origin). It returns an error if the plane is
// degenerate.
func (p Plane) PointOnPlane() (Vec3, error) {
	d := p.Normal.LengthSq()
	if d < 1e-300 {
		return Vec3{}, errors.New("hull3d: degenerate plane")
	}
	return p.Normal.Scale(p.Offset / d), nil
}

// String renders the plane's equation.
func (p Plane) String() string {
	return fmt.Sprintf("{n=%s, d=%g}", p.Normal, p.Offset)
}

// Orient3D returns a value whose sign tells the orientation of point d relative
// to the plane through a, b, c: positive when d lies on the positive side of the
// oriented triangle (a,b,c) (that is, above it with respect to the right-hand
// rule normal (b-a)×(c-a)), negative below, and zero when the four points are
// coplanar. The magnitude is six times the signed volume of the tetrahedron
// (a,b,c,d).
func Orient3D(a, b, c, d Vec3) float64 {
	return Triple(a.Sub(d), b.Sub(d), c.Sub(d))
}

// OrientSign returns the sign of [Orient3D] with a coplanarity tolerance eps:
// +1, -1, or 0 when |Orient3D| <= eps.
func OrientSign(a, b, c, d Vec3, eps float64) int {
	v := Orient3D(a, b, c, d)
	if v > eps {
		return 1
	}
	if v < -eps {
		return -1
	}
	return 0
}

// Coplanar reports whether the four points a, b, c, d are coplanar to within
// tolerance eps on the orientation determinant.
func Coplanar(a, b, c, d Vec3, eps float64) bool {
	return math.Abs(Orient3D(a, b, c, d)) <= eps
}

// Collinear reports whether the three points a, b, c are collinear to within
// tolerance eps on the squared area of the triangle they span.
func Collinear(a, b, c Vec3, eps float64) bool {
	return b.Sub(a).Cross(c.Sub(a)).LengthSq() <= eps*eps
}

// InSphere returns a value whose sign tells whether point e lies inside the
// sphere through a, b, c, d. When (a,b,c,d) is positively oriented
// ([Orient3D] > 0), a positive result means e is strictly inside the
// circumsphere, negative means outside, and zero means cocircular. Callers that
// do not control the orientation of the base tetrahedron should combine this
// with the sign of [Orient3D].
func InSphere(a, b, c, d, e Vec3) float64 {
	ax, ay, az := a.X-e.X, a.Y-e.Y, a.Z-e.Z
	bx, by, bz := b.X-e.X, b.Y-e.Y, b.Z-e.Z
	cx, cy, cz := c.X-e.X, c.Y-e.Y, c.Z-e.Z
	dx, dy, dz := d.X-e.X, d.Y-e.Y, d.Z-e.Z

	ab := ax*by - bx*ay
	bc := bx*cy - cx*by
	cd := cx*dy - dx*cy
	da := dx*ay - ax*dy
	ac := ax*cy - cx*ay
	bd := bx*dy - dx*by

	alift := ax*ax + ay*ay + az*az
	blift := bx*bx + by*by + bz*bz
	clift := cx*cx + cy*cy + cz*cz
	dlift := dx*dx + dy*dy + dz*dz

	abc := az*bc - bz*ac + cz*ab
	bcd := bz*cd - cz*bd + dz*bc
	cda := cz*da + dz*ac + az*cd
	dab := dz*ab + az*bd + bz*da

	return (dlift*abc - clift*dab) + (blift*cda - alift*bcd)
}

// SupportingPlane returns the plane through the triangle (a,b,c) oriented so
// that its positive half-space does not contain the interior reference point
// interior; that is, the returned normal points away from interior. It is used
// to give hull faces outward-facing normals. It returns an error if the three
// points are collinear.
func SupportingPlane(a, b, c, interior Vec3) (Plane, error) {
	pl, err := PlaneFromPoints(a, b, c)
	if err != nil {
		return Plane{}, err
	}
	if pl.Eval(interior) > 0 {
		pl = pl.Flip()
	}
	return pl, nil
}
