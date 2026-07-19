package polyhedra

import (
	"fmt"
	"math"
)

// Eps is the default absolute tolerance used by the approximate predicates in
// this package. It sits comfortably above double-precision rounding noise for
// coordinates of moderate magnitude.
const Eps = 1e-9

// Vec3 is a vector in three-dimensional Euclidean space. It is used both as a
// position (a point) and as a displacement or direction, depending on context.
type Vec3 struct {
	X, Y, Z float64
}

// NewVec3 returns the vector with the given components.
func NewVec3(x, y, z float64) Vec3 {
	return Vec3{X: x, Y: y, Z: z}
}

// ZeroVec3 returns the zero vector.
func ZeroVec3() Vec3 { return Vec3{} }

// XAxis returns the unit vector (1, 0, 0).
func XAxis() Vec3 { return Vec3{1, 0, 0} }

// YAxis returns the unit vector (0, 1, 0).
func YAxis() Vec3 { return Vec3{0, 1, 0} }

// ZAxis returns the unit vector (0, 0, 1).
func ZAxis() Vec3 { return Vec3{0, 0, 1} }

// Add returns the component-wise sum v + w.
func (v Vec3) Add(w Vec3) Vec3 { return Vec3{v.X + w.X, v.Y + w.Y, v.Z + w.Z} }

// Sub returns the component-wise difference v - w.
func (v Vec3) Sub(w Vec3) Vec3 { return Vec3{v.X - w.X, v.Y - w.Y, v.Z - w.Z} }

// Scale returns the vector v scaled by the scalar s.
func (v Vec3) Scale(s float64) Vec3 { return Vec3{v.X * s, v.Y * s, v.Z * s} }

// Div returns v with each component divided by the scalar s.
func (v Vec3) Div(s float64) Vec3 { return Vec3{v.X / s, v.Y / s, v.Z / s} }

// Neg returns the additive inverse -v.
func (v Vec3) Neg() Vec3 { return Vec3{-v.X, -v.Y, -v.Z} }

// Dot returns the Euclidean dot (inner) product v · w.
func (v Vec3) Dot(w Vec3) float64 { return v.X*w.X + v.Y*w.Y + v.Z*w.Z }

// Cross returns the cross product v × w, a vector orthogonal to both operands
// whose length equals the area of the parallelogram they span.
func (v Vec3) Cross(w Vec3) Vec3 {
	return Vec3{
		v.Y*w.Z - v.Z*w.Y,
		v.Z*w.X - v.X*w.Z,
		v.X*w.Y - v.Y*w.X,
	}
}

// LenSq returns the squared Euclidean length of v.
func (v Vec3) LenSq() float64 { return v.Dot(v) }

// Len returns the Euclidean length (2-norm) of v.
func (v Vec3) Len() float64 { return math.Sqrt(v.Dot(v)) }

// Normalize returns the unit vector in the direction of v and true, or the zero
// vector and false when v has length below Eps.
func (v Vec3) Normalize() (Vec3, bool) {
	n := v.Len()
	if n < Eps {
		return Vec3{}, false
	}
	return v.Div(n), true
}

// Unit returns the unit vector in the direction of v, or the zero vector if v is
// (numerically) zero. It is the panic-free convenience form of Normalize.
func (v Vec3) Unit() Vec3 {
	u, ok := v.Normalize()
	if !ok {
		return Vec3{}
	}
	return u
}

// DistanceSq returns the squared Euclidean distance between points v and w.
func (v Vec3) DistanceSq(w Vec3) float64 { return v.Sub(w).LenSq() }

// Distance returns the Euclidean distance between points v and w.
func (v Vec3) Distance(w Vec3) float64 { return v.Sub(w).Len() }

// Equal reports whether v and w are equal within absolute tolerance tol in each
// component.
func (v Vec3) Equal(w Vec3, tol float64) bool {
	return math.Abs(v.X-w.X) <= tol && math.Abs(v.Y-w.Y) <= tol && math.Abs(v.Z-w.Z) <= tol
}

// Lerp returns the linear interpolation (1-t)*v + t*w.
func (v Vec3) Lerp(w Vec3, t float64) Vec3 { return v.Add(w.Sub(v).Scale(t)) }

// Midpoint returns the point halfway between v and w.
func (v Vec3) Midpoint(w Vec3) Vec3 { return v.Add(w).Scale(0.5) }

// Angle returns the unsigned angle in radians between v and w, in [0, π]. It
// returns 0 when either operand is (numerically) zero.
func (v Vec3) Angle(w Vec3) float64 {
	uv, ok1 := v.Normalize()
	uw, ok2 := w.Normalize()
	if !ok1 || !ok2 {
		return 0
	}
	d := uv.Dot(uw)
	if d > 1 {
		d = 1
	} else if d < -1 {
		d = -1
	}
	return math.Acos(d)
}

// Project returns the vector projection of v onto w. It returns the zero vector
// when w is (numerically) zero.
func (v Vec3) Project(w Vec3) Vec3 {
	d := w.LenSq()
	if d < Eps*Eps {
		return Vec3{}
	}
	return w.Scale(v.Dot(w) / d)
}

// Reject returns the component of v orthogonal to w, that is v minus its
// projection onto w.
func (v Vec3) Reject(w Vec3) Vec3 { return v.Sub(v.Project(w)) }

// Reflect returns v reflected across the plane through the origin with unit
// normal n, i.e. v - 2 (v·n) n. The normal is assumed to be of unit length.
func (v Vec3) Reflect(n Vec3) Vec3 { return v.Sub(n.Scale(2 * v.Dot(n))) }

// IsZero reports whether every component of v is within Eps of zero.
func (v Vec3) IsZero() bool {
	return math.Abs(v.X) <= Eps && math.Abs(v.Y) <= Eps && math.Abs(v.Z) <= Eps
}

// Hadamard returns the component-wise (Hadamard) product of v and w.
func (v Vec3) Hadamard(w Vec3) Vec3 { return Vec3{v.X * w.X, v.Y * w.Y, v.Z * w.Z} }

// MinComp returns the component-wise minimum of v and w.
func (v Vec3) MinComp(w Vec3) Vec3 {
	return Vec3{math.Min(v.X, w.X), math.Min(v.Y, w.Y), math.Min(v.Z, w.Z)}
}

// MaxComp returns the component-wise maximum of v and w.
func (v Vec3) MaxComp(w Vec3) Vec3 {
	return Vec3{math.Max(v.X, w.X), math.Max(v.Y, w.Y), math.Max(v.Z, w.Z)}
}

// Abs returns the component-wise absolute value of v.
func (v Vec3) Abs() Vec3 { return Vec3{math.Abs(v.X), math.Abs(v.Y), math.Abs(v.Z)} }

// MaxAbsComp returns the largest absolute value among the components of v.
func (v Vec3) MaxAbsComp() float64 {
	return math.Max(math.Abs(v.X), math.Max(math.Abs(v.Y), math.Abs(v.Z)))
}

// Triple returns the scalar triple product v · (w × u), equal to the signed
// volume of the parallelepiped spanned by the three vectors.
func (v Vec3) Triple(w, u Vec3) float64 { return v.Dot(w.Cross(u)) }

// String renders v in a compact, deterministic textual form.
func (v Vec3) String() string {
	return fmt.Sprintf("(%g, %g, %g)", v.X, v.Y, v.Z)
}

// Sum returns the vector sum of all points in ps (the zero vector when empty).
func Sum(ps []Vec3) Vec3 {
	var s Vec3
	for _, p := range ps {
		s = s.Add(p)
	}
	return s
}

// Centroid returns the arithmetic mean of the points in ps and true, or the
// zero vector and false when ps is empty.
func Centroid(ps []Vec3) (Vec3, bool) {
	if len(ps) == 0 {
		return Vec3{}, false
	}
	return Sum(ps).Div(float64(len(ps))), true
}
