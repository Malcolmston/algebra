package meshgen

import (
	"fmt"
	"math"
)

// Vec2 is a point or vector in the Euclidean plane.
type Vec2 struct {
	X, Y float64
}

// NewVec2 returns the Vec2 with the given coordinates.
func NewVec2(x, y float64) Vec2 { return Vec2{X: x, Y: y} }

// Vec2Zero is the origin (0, 0).
func Vec2Zero() Vec2 { return Vec2{0, 0} }

// Add returns the componentwise sum p + q.
func (p Vec2) Add(q Vec2) Vec2 { return Vec2{p.X + q.X, p.Y + q.Y} }

// Sub returns the componentwise difference p - q.
func (p Vec2) Sub(q Vec2) Vec2 { return Vec2{p.X - q.X, p.Y - q.Y} }

// Scale returns p scaled by the factor s.
func (p Vec2) Scale(s float64) Vec2 { return Vec2{p.X * s, p.Y * s} }

// Div returns p divided componentwise by the nonzero scalar s.
func (p Vec2) Div(s float64) Vec2 { return Vec2{p.X / s, p.Y / s} }

// Neg returns the vector reflected through the origin, -p.
func (p Vec2) Neg() Vec2 { return Vec2{-p.X, -p.Y} }

// Dot returns the dot product p . q.
func (p Vec2) Dot(q Vec2) float64 { return p.X*q.X + p.Y*q.Y }

// Cross returns the scalar 2-D cross product p.X*q.Y - p.Y*q.X.
func (p Vec2) Cross(q Vec2) float64 { return p.X*q.Y - p.Y*q.X }

// NormSq returns the squared Euclidean norm |p|^2.
func (p Vec2) NormSq() float64 { return p.X*p.X + p.Y*p.Y }

// Norm returns the Euclidean norm |p|.
func (p Vec2) Norm() float64 { return math.Hypot(p.X, p.Y) }

// DistanceSq returns the squared Euclidean distance between p and q.
func (p Vec2) DistanceSq(q Vec2) float64 { return p.Sub(q).NormSq() }

// Distance returns the Euclidean distance between p and q.
func (p Vec2) Distance(q Vec2) float64 { return p.Sub(q).Norm() }

// Normalize returns the unit vector in the direction of p. The zero vector is
// returned unchanged.
func (p Vec2) Normalize() Vec2 {
	n := p.Norm()
	if n == 0 {
		return p
	}
	return p.Scale(1 / n)
}

// Perp returns the vector p rotated 90 degrees counterclockwise, (-p.Y, p.X).
func (p Vec2) Perp() Vec2 { return Vec2{-p.Y, p.X} }

// Rotate returns p rotated counterclockwise by theta radians about the origin.
func (p Vec2) Rotate(theta float64) Vec2 {
	s, c := math.Sincos(theta)
	return Vec2{p.X*c - p.Y*s, p.X*s + p.Y*c}
}

// RotateAbout returns p rotated counterclockwise by theta radians about c.
func (p Vec2) RotateAbout(c Vec2, theta float64) Vec2 {
	return p.Sub(c).Rotate(theta).Add(c)
}

// Lerp returns the linear interpolation (1-t)*p + t*q.
func (p Vec2) Lerp(q Vec2, t float64) Vec2 {
	return Vec2{p.X + (q.X-p.X)*t, p.Y + (q.Y-p.Y)*t}
}

// Midpoint returns the midpoint of the segment p-q.
func (p Vec2) Midpoint(q Vec2) Vec2 { return Vec2{(p.X + q.X) / 2, (p.Y + q.Y) / 2} }

// Hadamard returns the componentwise product of p and q.
func (p Vec2) Hadamard(q Vec2) Vec2 { return Vec2{p.X * q.X, p.Y * q.Y} }

// Min returns the componentwise minimum of p and q.
func (p Vec2) Min(q Vec2) Vec2 { return Vec2{math.Min(p.X, q.X), math.Min(p.Y, q.Y)} }

// Max returns the componentwise maximum of p and q.
func (p Vec2) Max(q Vec2) Vec2 { return Vec2{math.Max(p.X, q.X), math.Max(p.Y, q.Y)} }

// Abs returns the componentwise absolute value of p.
func (p Vec2) Abs() Vec2 { return Vec2{math.Abs(p.X), math.Abs(p.Y)} }

// Angle returns the angle of p measured counterclockwise from the +X axis, in
// the range (-pi, pi].
func (p Vec2) Angle() float64 { return math.Atan2(p.Y, p.X) }

// AngleTo returns the signed angle from p to q in radians, in (-pi, pi].
func (p Vec2) AngleTo(q Vec2) float64 {
	return math.Atan2(p.Cross(q), p.Dot(q))
}

// IsZero reports whether p is exactly the zero vector.
func (p Vec2) IsZero() bool { return p.X == 0 && p.Y == 0 }

// Equal reports whether p and q are exactly equal.
func (p Vec2) Equal(q Vec2) bool { return p.X == q.X && p.Y == q.Y }

// ApproxEqual reports whether p and q agree to within absolute tolerance eps.
func (p Vec2) ApproxEqual(q Vec2, eps float64) bool {
	return math.Abs(p.X-q.X) <= eps && math.Abs(p.Y-q.Y) <= eps
}

// IsFinite reports whether both coordinates of p are finite.
func (p Vec2) IsFinite() bool {
	return !math.IsInf(p.X, 0) && !math.IsInf(p.Y, 0) &&
		!math.IsNaN(p.X) && !math.IsNaN(p.Y)
}

// String returns a human-readable representation of p.
func (p Vec2) String() string { return fmt.Sprintf("(%g, %g)", p.X, p.Y) }

// CentroidVec2 returns the arithmetic mean of the given points. It returns the
// origin for an empty slice.
func CentroidVec2(pts []Vec2) Vec2 {
	if len(pts) == 0 {
		return Vec2{}
	}
	var s Vec2
	for _, p := range pts {
		s = s.Add(p)
	}
	return s.Div(float64(len(pts)))
}

// BoundingBox2 returns the axis-aligned minimum and maximum corners enclosing
// the given points. It returns two zero vectors for an empty slice.
func BoundingBox2(pts []Vec2) (min, max Vec2) {
	if len(pts) == 0 {
		return Vec2{}, Vec2{}
	}
	min, max = pts[0], pts[0]
	for _, p := range pts[1:] {
		min = min.Min(p)
		max = max.Max(p)
	}
	return min, max
}
