// Package hypercomplex implements arithmetic over several hypercomplex number
// systems that extend the real and complex numbers.
//
// The package provides:
//
//   - Quaternion: Hamilton's quaternions, the four-dimensional normed division
//     algebra. Includes multiplication, conjugation, norm, inverse, the
//     exponential and logarithm maps, power, spherical linear interpolation
//     (Slerp), rotation of 3-vectors, and conversions to and from axis-angle,
//     Euler angles, and rotation matrices.
//   - DualQuaternion: the tensor of the dual numbers with the quaternions,
//     used to represent rigid-body transforms (screw motions) in a single
//     algebraic object.
//   - Octonion: the eight-dimensional non-associative normed division algebra,
//     built from quaternions by the Cayley-Dickson construction.
//   - SplitComplex: the split-complex (hyperbolic) numbers x + y*j with
//     j*j = +1, the algebra underlying Lorentz boosts in 1+1 dimensions.
//
// Every routine uses only the Go standard library and is deterministic. Angles
// are expressed in radians throughout, and 3-vectors are represented by the
// [Vec3] type defined in this package.
package hypercomplex

import "math"

// Vec3 is a three-dimensional real vector. It is used for the vector (imaginary)
// part of quaternions and as the argument and result of rotations.
type Vec3 struct {
	X, Y, Z float64
}

// V3 constructs a [Vec3] from its three components.
func V3(x, y, z float64) Vec3 {
	return Vec3{X: x, Y: y, Z: z}
}

// Add returns the component-wise sum v + w.
func (v Vec3) Add(w Vec3) Vec3 {
	return Vec3{v.X + w.X, v.Y + w.Y, v.Z + w.Z}
}

// Sub returns the component-wise difference v - w.
func (v Vec3) Sub(w Vec3) Vec3 {
	return Vec3{v.X - w.X, v.Y - w.Y, v.Z - w.Z}
}

// Scale returns the vector v scaled by the real factor s.
func (v Vec3) Scale(s float64) Vec3 {
	return Vec3{v.X * s, v.Y * s, v.Z * s}
}

// Neg returns the additive inverse -v.
func (v Vec3) Neg() Vec3 {
	return Vec3{-v.X, -v.Y, -v.Z}
}

// Dot returns the Euclidean dot product v·w.
func (v Vec3) Dot(w Vec3) float64 {
	return v.X*w.X + v.Y*w.Y + v.Z*w.Z
}

// Cross returns the cross product v×w.
func (v Vec3) Cross(w Vec3) Vec3 {
	return Vec3{
		v.Y*w.Z - v.Z*w.Y,
		v.Z*w.X - v.X*w.Z,
		v.X*w.Y - v.Y*w.X,
	}
}

// Norm returns the Euclidean length |v|.
func (v Vec3) Norm() float64 {
	return math.Sqrt(v.Dot(v))
}

// Normalize returns v scaled to unit length. If v is the zero vector it is
// returned unchanged.
func (v Vec3) Normalize() Vec3 {
	n := v.Norm()
	if n == 0 {
		return v
	}
	return v.Scale(1 / n)
}

// Equal reports whether v and w agree component-wise to within the absolute
// tolerance tol.
func (v Vec3) Equal(w Vec3, tol float64) bool {
	return math.Abs(v.X-w.X) <= tol &&
		math.Abs(v.Y-w.Y) <= tol &&
		math.Abs(v.Z-w.Z) <= tol
}
