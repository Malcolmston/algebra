package physics

import "math"

// Vec3 is a three-dimensional Cartesian vector with float64 components. It is a
// small value type: all methods take and return Vec3 by value, so ordinary
// vector arithmetic performs no heap allocation and the compiler can keep
// intermediate vectors in registers. Vec3 is used throughout the package for
// positions, velocities, forces and field samples.
type Vec3 struct {
	// X is the component along the first Cartesian axis.
	X float64
	// Y is the component along the second Cartesian axis.
	Y float64
	// Z is the component along the third Cartesian axis.
	Z float64
}

// NewVec3 returns the vector with components (x, y, z).
func NewVec3(x, y, z float64) Vec3 { return Vec3{x, y, z} }

// Add returns the component-wise sum a + b.
func (a Vec3) Add(b Vec3) Vec3 { return Vec3{a.X + b.X, a.Y + b.Y, a.Z + b.Z} }

// Sub returns the component-wise difference a - b.
func (a Vec3) Sub(b Vec3) Vec3 { return Vec3{a.X - b.X, a.Y - b.Y, a.Z - b.Z} }

// Scale returns the vector a with every component multiplied by the scalar s.
func (a Vec3) Scale(s float64) Vec3 { return Vec3{a.X * s, a.Y * s, a.Z * s} }

// Neg returns the vector pointing opposite to a, that is a scaled by -1.
func (a Vec3) Neg() Vec3 { return Vec3{-a.X, -a.Y, -a.Z} }

// Dot returns the scalar (inner) product a·b = aₓbₓ + a_yb_y + a_zb_z.
func (a Vec3) Dot(b Vec3) float64 { return a.X*b.X + a.Y*b.Y + a.Z*b.Z }

// Cross returns the vector (cross) product a×b, a vector orthogonal to both a
// and b whose magnitude equals the area of the parallelogram they span and
// whose direction follows the right-hand rule.
func (a Vec3) Cross(b Vec3) Vec3 {
	return Vec3{
		a.Y*b.Z - a.Z*b.Y,
		a.Z*b.X - a.X*b.Z,
		a.X*b.Y - a.Y*b.X,
	}
}

// Norm returns the Euclidean length (magnitude) |a| = √(a·a) of the vector.
func (a Vec3) Norm() float64 { return math.Sqrt(a.X*a.X + a.Y*a.Y + a.Z*a.Z) }

// Norm2 returns the squared Euclidean length a·a. It avoids the square root of
// [Vec3.Norm] and is preferred when only comparing magnitudes.
func (a Vec3) Norm2() float64 { return a.X*a.X + a.Y*a.Y + a.Z*a.Z }

// Normalize returns the unit vector pointing in the direction of a. If a is the
// zero vector it is returned unchanged, since its direction is undefined.
func (a Vec3) Normalize() Vec3 {
	n := a.Norm()
	if n == 0 {
		return a
	}
	inv := 1 / n
	return Vec3{a.X * inv, a.Y * inv, a.Z * inv}
}

// Distance returns the Euclidean distance |a - b| between the points a and b.
func (a Vec3) Distance(b Vec3) float64 {
	dx, dy, dz := a.X-b.X, a.Y-b.Y, a.Z-b.Z
	return math.Sqrt(dx*dx + dy*dy + dz*dz)
}

// Lerp returns the linear interpolation (1-t)·a + t·b between a and b. t = 0
// yields a and t = 1 yields b; values outside [0, 1] extrapolate.
func (a Vec3) Lerp(b Vec3, t float64) Vec3 {
	return Vec3{
		a.X + (b.X-a.X)*t,
		a.Y + (b.Y-a.Y)*t,
		a.Z + (b.Z-a.Z)*t,
	}
}

// Angle returns the unsigned angle in radians between a and b, in the range
// [0, π]. If either vector is the zero vector the result is 0.
func (a Vec3) Angle(b Vec3) float64 {
	na, nb := a.Norm(), b.Norm()
	if na == 0 || nb == 0 {
		return 0
	}
	c := physicsClampUnit(a.Dot(b) / (na * nb))
	return math.Acos(c)
}

// physicsClampUnit clamps x to the closed interval [-1, 1], guarding acos and
// asin against rounding that pushes a cosine slightly outside its valid range.
func physicsClampUnit(x float64) float64 {
	if x > 1 {
		return 1
	}
	if x < -1 {
		return -1
	}
	return x
}
