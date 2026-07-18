// Package geom3d provides primitives for 3D computational geometry.
//
// The core value type is [Vec3], a three-component vector used both as a
// location (point) and as a displacement (direction). On top of it the
// package offers the standard toolbox of spatial geometry: vector arithmetic
// (add, subtract, scale, dot, cross, triple product), norms, distances,
// angles, normalization, projection/rejection and reflection; 3x3 matrix
// algebra ([Mat3]) with determinant, inverse, transpose and multiplication;
// rotation matrices about the coordinate axes and about an arbitrary axis
// (Rodrigues' rotation formula); conversion between Euler angles and rotation
// matrices; planes ([Plane]) and lines ([Line3]) with closest-point and
// distance queries; ray casting ([Ray]) against planes, triangles
// (Moeller-Trumbore), axis-aligned boxes ([AABB]) and spheres ([Sphere]);
// barycentric coordinates; and intersection tests between boxes and spheres.
//
// All computation is done in float64 using only the Go standard library.
// Routines are deterministic: given identical inputs they return identical
// results. Predicates that must tolerate floating-point noise accept or use an
// explicit epsilon so callers control the trade-off between robustness and
// strictness.
package geom3d

import "math"

// geom3dEps is the default absolute tolerance used by the approximate
// predicates in this package. It is comfortably above double-precision
// rounding noise for coordinates of moderate magnitude.
const geom3dEps = 1e-9

// Vec3 is a vector in three-dimensional Euclidean space. It is used both as a
// position (a point) and as a displacement or direction, depending on context.
type Vec3 struct {
	X, Y, Z float64
}

// NewVec3 returns the vector with the given components.
func NewVec3(x, y, z float64) Vec3 {
	return Vec3{X: x, Y: y, Z: z}
}

// Add returns the component-wise sum v+w.
func (v Vec3) Add(w Vec3) Vec3 {
	return Vec3{v.X + w.X, v.Y + w.Y, v.Z + w.Z}
}

// Sub returns the component-wise difference v-w.
func (v Vec3) Sub(w Vec3) Vec3 {
	return Vec3{v.X - w.X, v.Y - w.Y, v.Z - w.Z}
}

// Scale returns v with every component multiplied by the scalar s.
func (v Vec3) Scale(s float64) Vec3 {
	return Vec3{v.X * s, v.Y * s, v.Z * s}
}

// Neg returns the additive inverse -v.
func (v Vec3) Neg() Vec3 {
	return Vec3{-v.X, -v.Y, -v.Z}
}

// Div returns v with every component divided by the scalar s. It panics if s
// is zero.
func (v Vec3) Div(s float64) Vec3 {
	if s == 0 {
		panic("geom3d: division of Vec3 by zero")
	}
	return Vec3{v.X / s, v.Y / s, v.Z / s}
}

// Hadamard returns the component-wise (Hadamard) product of v and w.
func (v Vec3) Hadamard(w Vec3) Vec3 {
	return Vec3{v.X * w.X, v.Y * w.Y, v.Z * w.Z}
}

// Dot returns the dot (inner) product v·w.
func (v Vec3) Dot(w Vec3) float64 {
	return v.X*w.X + v.Y*w.Y + v.Z*w.Z
}

// Cross returns the cross (vector) product v×w. The result is orthogonal to
// both v and w and follows the right-hand rule.
func (v Vec3) Cross(w Vec3) Vec3 {
	return Vec3{
		v.Y*w.Z - v.Z*w.Y,
		v.Z*w.X - v.X*w.Z,
		v.X*w.Y - v.Y*w.X,
	}
}

// LengthSq returns the squared Euclidean length v·v. It avoids the square root
// of [Vec3.Length] and is preferred when only relative magnitudes matter.
func (v Vec3) LengthSq() float64 {
	return v.X*v.X + v.Y*v.Y + v.Z*v.Z
}

// Length returns the Euclidean length (magnitude) of v.
func (v Vec3) Length() float64 {
	return math.Sqrt(v.LengthSq())
}

// Normalize returns the unit vector pointing in the same direction as v,
// together with the original length. If v has zero length the zero vector and
// a length of 0 are returned.
func (v Vec3) Normalize() (Vec3, float64) {
	n := v.Length()
	if n == 0 {
		return Vec3{}, 0
	}
	return Vec3{v.X / n, v.Y / n, v.Z / n}, n
}

// Unit returns the unit vector in the direction of v, or the zero vector if v
// has zero length. It is a convenience wrapper over [Vec3.Normalize] that
// discards the length.
func (v Vec3) Unit() Vec3 {
	u, _ := v.Normalize()
	return u
}

// DistanceSq returns the squared Euclidean distance between the points v and w.
func (v Vec3) DistanceSq(w Vec3) float64 {
	return v.Sub(w).LengthSq()
}

// Distance returns the Euclidean distance between the points v and w.
func (v Vec3) Distance(w Vec3) float64 {
	return v.Sub(w).Length()
}

// Lerp returns the linear interpolation (1-t)*v + t*w. With t=0 it returns v
// and with t=1 it returns w; values outside [0,1] extrapolate.
func (v Vec3) Lerp(w Vec3, t float64) Vec3 {
	return Vec3{
		v.X + (w.X-v.X)*t,
		v.Y + (w.Y-v.Y)*t,
		v.Z + (w.Z-v.Z)*t,
	}
}

// Angle returns the unsigned angle in radians between v and w, in the range
// [0, pi]. It returns 0 if either vector has zero length.
func (v Vec3) Angle(w Vec3) float64 {
	d := math.Sqrt(v.LengthSq() * w.LengthSq())
	if d == 0 {
		return 0
	}
	c := v.Dot(w) / d
	c = geom3dclamp(c, -1, 1)
	return math.Acos(c)
}

// Project returns the vector projection of v onto w, that is the component of
// v parallel to w. It returns the zero vector if w has zero length.
func (v Vec3) Project(w Vec3) Vec3 {
	d := w.LengthSq()
	if d == 0 {
		return Vec3{}
	}
	return w.Scale(v.Dot(w) / d)
}

// Reject returns the vector rejection of v from w, that is the component of v
// orthogonal to w: v minus its projection onto w.
func (v Vec3) Reject(w Vec3) Vec3 {
	return v.Sub(v.Project(w))
}

// Reflect returns v reflected across the plane whose unit normal is n. The
// caller is responsible for supplying a unit-length normal; the reflection is
// v - 2*(v·n)*n.
func (v Vec3) Reflect(n Vec3) Vec3 {
	return v.Sub(n.Scale(2 * v.Dot(n)))
}

// Min returns the component-wise minimum of v and w.
func (v Vec3) Min(w Vec3) Vec3 {
	return Vec3{math.Min(v.X, w.X), math.Min(v.Y, w.Y), math.Min(v.Z, w.Z)}
}

// Max returns the component-wise maximum of v and w.
func (v Vec3) Max(w Vec3) Vec3 {
	return Vec3{math.Max(v.X, w.X), math.Max(v.Y, w.Y), math.Max(v.Z, w.Z)}
}

// Abs returns the component-wise absolute value of v.
func (v Vec3) Abs() Vec3 {
	return Vec3{math.Abs(v.X), math.Abs(v.Y), math.Abs(v.Z)}
}

// IsZero reports whether every component of v is zero within the package
// default tolerance.
func (v Vec3) IsZero() bool {
	return math.Abs(v.X) <= geom3dEps &&
		math.Abs(v.Y) <= geom3dEps &&
		math.Abs(v.Z) <= geom3dEps
}

// Equal reports whether v and w are equal within the absolute tolerance eps,
// compared component-wise.
func (v Vec3) Equal(w Vec3, eps float64) bool {
	return math.Abs(v.X-w.X) <= eps &&
		math.Abs(v.Y-w.Y) <= eps &&
		math.Abs(v.Z-w.Z) <= eps
}

// ScalarTriple returns the scalar triple product a·(b×c), which equals the
// signed volume of the parallelepiped spanned by a, b and c.
func ScalarTriple(a, b, c Vec3) float64 {
	return a.Dot(b.Cross(c))
}

// geom3dclamp returns x limited to the closed interval [lo, hi].
func geom3dclamp(x, lo, hi float64) float64 {
	if x < lo {
		return lo
	}
	if x > hi {
		return hi
	}
	return x
}
