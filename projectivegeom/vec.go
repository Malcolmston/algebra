package projectivegeom

import (
	"fmt"
	"math"
)

// Eps is the default absolute tolerance used by the approximate predicates in
// this package. It is comfortably above double-precision rounding noise for
// coordinates of moderate magnitude.
const Eps = 1e-9

// Vec3 is a 3-component real vector. It is the common carrier for homogeneous
// points and lines of RP^2 and for ordinary vectors of R^3.
type Vec3 struct {
	X, Y, Z float64
}

// NewVec3 returns the vector [x y z].
func NewVec3(x, y, z float64) Vec3 { return Vec3{x, y, z} }

// Add returns the component-wise sum v+w.
func (v Vec3) Add(w Vec3) Vec3 { return Vec3{v.X + w.X, v.Y + w.Y, v.Z + w.Z} }

// Sub returns the component-wise difference v-w.
func (v Vec3) Sub(w Vec3) Vec3 { return Vec3{v.X - w.X, v.Y - w.Y, v.Z - w.Z} }

// Scale returns v scaled by the scalar s.
func (v Vec3) Scale(s float64) Vec3 { return Vec3{v.X * s, v.Y * s, v.Z * s} }

// Neg returns the negation -v.
func (v Vec3) Neg() Vec3 { return Vec3{-v.X, -v.Y, -v.Z} }

// Dot returns the Euclidean dot product v·w.
func (v Vec3) Dot(w Vec3) float64 { return v.X*w.X + v.Y*w.Y + v.Z*w.Z }

// Cross returns the vector cross product v×w. In RP^2 this is simultaneously
// the line joining two points and the point meeting two lines.
func (v Vec3) Cross(w Vec3) Vec3 {
	return Vec3{
		v.Y*w.Z - v.Z*w.Y,
		v.Z*w.X - v.X*w.Z,
		v.X*w.Y - v.Y*w.X,
	}
}

// Norm2 returns the squared Euclidean norm v·v.
func (v Vec3) Norm2() float64 { return v.Dot(v) }

// Norm returns the Euclidean norm |v|.
func (v Vec3) Norm() float64 { return math.Sqrt(v.Norm2()) }

// Normalized returns v scaled to unit Euclidean length and true, or the zero
// vector and false when v is (numerically) zero.
func (v Vec3) Normalized() (Vec3, bool) {
	n := v.Norm()
	if n < Eps {
		return Vec3{}, false
	}
	return v.Scale(1 / n), true
}

// IsZero reports whether every component of v is within tol of zero.
func (v Vec3) IsZero(tol float64) bool {
	return math.Abs(v.X) <= tol && math.Abs(v.Y) <= tol && math.Abs(v.Z) <= tol
}

// ApproxEqual reports whether v and w agree component-wise within tol.
func (v Vec3) ApproxEqual(w Vec3, tol float64) bool {
	return math.Abs(v.X-w.X) <= tol && math.Abs(v.Y-w.Y) <= tol && math.Abs(v.Z-w.Z) <= tol
}

// Parallel reports whether v and w are scalar multiples of one another within
// tol, i.e. their cross product vanishes.
func (v Vec3) Parallel(w Vec3, tol float64) bool {
	return v.Cross(w).IsZero(tol)
}

// Max returns the largest component of v.
func (v Vec3) Max() float64 { return math.Max(v.X, math.Max(v.Y, v.Z)) }

// MaxAbs returns the largest absolute value among the components of v.
func (v Vec3) MaxAbs() float64 {
	return math.Max(math.Abs(v.X), math.Max(math.Abs(v.Y), math.Abs(v.Z)))
}

// Get returns component i (0=X, 1=Y, 2=Z). It panics on an out-of-range index.
func (v Vec3) Get(i int) float64 {
	switch i {
	case 0:
		return v.X
	case 1:
		return v.Y
	case 2:
		return v.Z
	}
	panic("projectivegeom: Vec3 index out of range")
}

// AngleBetween returns the unsigned angle in radians between v and w treated as
// direction vectors in R^3, in the range [0, pi].
func (v Vec3) AngleBetween(w Vec3) float64 {
	nv, nw := v.Norm(), w.Norm()
	if nv < Eps || nw < Eps {
		return 0
	}
	c := v.Dot(w) / (nv * nw)
	c = math.Max(-1, math.Min(1, c))
	return math.Acos(c)
}

// String renders v in bracketed homogeneous form.
func (v Vec3) String() string { return fmt.Sprintf("[%g %g %g]", v.X, v.Y, v.Z) }

// Vec4 is a 4-component real vector, the carrier for homogeneous points and
// planes of RP^3.
type Vec4 struct {
	X, Y, Z, W float64
}

// NewVec4 returns the vector [x y z w].
func NewVec4(x, y, z, w float64) Vec4 { return Vec4{x, y, z, w} }

// Add returns the component-wise sum v+w.
func (v Vec4) Add(w Vec4) Vec4 { return Vec4{v.X + w.X, v.Y + w.Y, v.Z + w.Z, v.W + w.W} }

// Sub returns the component-wise difference v-w.
func (v Vec4) Sub(w Vec4) Vec4 { return Vec4{v.X - w.X, v.Y - w.Y, v.Z - w.Z, v.W - w.W} }

// Scale returns v scaled by the scalar s.
func (v Vec4) Scale(s float64) Vec4 { return Vec4{v.X * s, v.Y * s, v.Z * s, v.W * s} }

// Neg returns the negation -v.
func (v Vec4) Neg() Vec4 { return Vec4{-v.X, -v.Y, -v.Z, -v.W} }

// Dot returns the Euclidean dot product v·w.
func (v Vec4) Dot(w Vec4) float64 { return v.X*w.X + v.Y*w.Y + v.Z*w.Z + v.W*w.W }

// Norm2 returns the squared Euclidean norm v·v.
func (v Vec4) Norm2() float64 { return v.Dot(v) }

// Norm returns the Euclidean norm |v|.
func (v Vec4) Norm() float64 { return math.Sqrt(v.Norm2()) }

// Normalized returns v scaled to unit Euclidean length and true, or the zero
// vector and false when v is (numerically) zero.
func (v Vec4) Normalized() (Vec4, bool) {
	n := v.Norm()
	if n < Eps {
		return Vec4{}, false
	}
	return v.Scale(1 / n), true
}

// IsZero reports whether every component of v is within tol of zero.
func (v Vec4) IsZero(tol float64) bool {
	return math.Abs(v.X) <= tol && math.Abs(v.Y) <= tol &&
		math.Abs(v.Z) <= tol && math.Abs(v.W) <= tol
}

// ApproxEqual reports whether v and w agree component-wise within tol.
func (v Vec4) ApproxEqual(w Vec4, tol float64) bool {
	return math.Abs(v.X-w.X) <= tol && math.Abs(v.Y-w.Y) <= tol &&
		math.Abs(v.Z-w.Z) <= tol && math.Abs(v.W-w.W) <= tol
}

// MaxAbs returns the largest absolute value among the components of v.
func (v Vec4) MaxAbs() float64 {
	return math.Max(math.Max(math.Abs(v.X), math.Abs(v.Y)),
		math.Max(math.Abs(v.Z), math.Abs(v.W)))
}

// Get returns component i (0=X, 1=Y, 2=Z, 3=W). It panics on an out-of-range
// index.
func (v Vec4) Get(i int) float64 {
	switch i {
	case 0:
		return v.X
	case 1:
		return v.Y
	case 2:
		return v.Z
	case 3:
		return v.W
	}
	panic("projectivegeom: Vec4 index out of range")
}

// String renders v in bracketed homogeneous form.
func (v Vec4) String() string { return fmt.Sprintf("[%g %g %g %g]", v.X, v.Y, v.Z, v.W) }

// Cross4 returns the generalized cross product of three 4-vectors: the unique
// (up to scale) vector orthogonal to a, b and c. It is used to intersect three
// planes at a point or to span a plane through three points in RP^3.
func Cross4(a, b, c Vec4) Vec4 {
	// Component i is (-1)^i times the 3x3 minor obtained by deleting column i
	// from the 3x4 matrix whose rows are a, b, c.
	m := func(c0, c1, c2 int) float64 {
		r := [3][4]float64{
			{a.X, a.Y, a.Z, a.W},
			{b.X, b.Y, b.Z, b.W},
			{c.X, c.Y, c.Z, c.W},
		}
		return det3(
			r[0][c0], r[0][c1], r[0][c2],
			r[1][c0], r[1][c1], r[1][c2],
			r[2][c0], r[2][c1], r[2][c2],
		)
	}
	return Vec4{
		X: +m(1, 2, 3),
		Y: -m(0, 2, 3),
		Z: +m(0, 1, 3),
		W: -m(0, 1, 2),
	}
}

// det3 returns the determinant of the 3x3 matrix given in row-major order.
func det3(a, b, c, d, e, f, g, h, i float64) float64 {
	return a*(e*i-f*h) - b*(d*i-f*g) + c*(d*h-e*g)
}
