// Package diffgeo implements classical differential geometry of curves and
// surfaces in three-dimensional Euclidean space using only the Go standard
// library.
//
// The package is organized around two functional descriptions of geometric
// objects:
//
//   - A [Curve] is a function t ↦ (x, y, z) tracing a parametric space curve.
//     From it the package computes the Frenet apparatus: unit tangent, principal
//     normal and binormal ([FrenetFrame]), signed [Curvature] and [Torsion], the
//     curvature vector, radius of curvature and [ArcLength].
//
//   - A [Surface] is a function (u, v) ↦ (x, y, z) describing a parametric
//     patch. From it the package computes the [FirstForm] and [SecondForm]
//     (the first and second fundamental forms), the unit [SurfaceNormal],
//     [GaussianCurvature], [MeanCurvature], [PrincipalCurvatures],
//     [NormalCurvature], the [AreaElement] and total [SurfaceArea], the
//     [Christoffel] symbols of the second kind, and numerically integrated
//     geodesics ([GeodesicPath]) and [ParallelTransport].
//
// Derivatives are approximated by central finite differences with carefully
// chosen step sizes, so every routine works for an arbitrary caller-supplied
// parametrization without symbolic differentiation. Ready-made
// parametrizations are provided as constructors ([Line], [Circle], [Ellipse],
// [Helix], [Sphere], [Cylinder], [Torus], [Graph], [SurfacePlane]) so that
// results can be checked against the well-known closed-form curvatures of these
// standard shapes.
//
// All computation is done in float64. Routines are deterministic: identical
// inputs yield identical outputs. Approximate predicates take an explicit
// epsilon so the caller controls the robustness/strictness trade-off.
package diffgeo

import "math"

// Vec3 is a vector in three-dimensional Euclidean space. It doubles as a point
// (a position) and as a displacement or direction, depending on context.
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

// Dot returns the Euclidean inner product v·w.
func (v Vec3) Dot(w Vec3) float64 {
	return v.X*w.X + v.Y*w.Y + v.Z*w.Z
}

// Cross returns the cross product v×w, a vector orthogonal to both operands.
func (v Vec3) Cross(w Vec3) Vec3 {
	return Vec3{
		v.Y*w.Z - v.Z*w.Y,
		v.Z*w.X - v.X*w.Z,
		v.X*w.Y - v.Y*w.X,
	}
}

// Norm returns the Euclidean length |v|.
func (v Vec3) Norm() float64 {
	return math.Sqrt(v.X*v.X + v.Y*v.Y + v.Z*v.Z)
}

// Norm2 returns the squared Euclidean length |v|², avoiding the square root.
func (v Vec3) Norm2() float64 {
	return v.X*v.X + v.Y*v.Y + v.Z*v.Z
}

// Normalize returns the unit vector pointing along v. If v is the zero vector
// (within [Eps]) it returns the zero vector unchanged.
func (v Vec3) Normalize() Vec3 {
	n := v.Norm()
	if n < Eps {
		return Vec3{}
	}
	return v.Scale(1 / n)
}

// Distance returns the Euclidean distance between the points v and w.
func (v Vec3) Distance(w Vec3) float64 {
	return v.Sub(w).Norm()
}

// Lerp returns the linear interpolation (1-t)·v + t·w.
func (v Vec3) Lerp(w Vec3, t float64) Vec3 {
	return v.Scale(1 - t).Add(w.Scale(t))
}

// Angle returns the unsigned angle in radians between v and w, in [0, π].
// If either operand is the zero vector it returns 0.
func (v Vec3) Angle(w Vec3) float64 {
	d := v.Norm() * w.Norm()
	if d < Eps {
		return 0
	}
	c := v.Dot(w) / d
	if c > 1 {
		c = 1
	} else if c < -1 {
		c = -1
	}
	return math.Acos(c)
}

// IsZero reports whether every component of v is within eps of zero.
func (v Vec3) IsZero(eps float64) bool {
	return math.Abs(v.X) <= eps && math.Abs(v.Y) <= eps && math.Abs(v.Z) <= eps
}

// Equal reports whether v and w agree component-wise within the tolerance eps.
func (v Vec3) Equal(w Vec3, eps float64) bool {
	return math.Abs(v.X-w.X) <= eps &&
		math.Abs(v.Y-w.Y) <= eps &&
		math.Abs(v.Z-w.Z) <= eps
}

// Eps is the default absolute tolerance used by this package's approximate
// predicates and degeneracy checks. It sits comfortably above double-precision
// rounding noise for quantities of moderate magnitude.
const Eps = 1e-9
