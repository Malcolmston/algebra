package meshgen

import (
	"fmt"
	"math"
)

// Vec3 is a point or vector in three-dimensional Euclidean space.
type Vec3 struct {
	X, Y, Z float64
}

// NewVec3 returns the Vec3 with the given coordinates.
func NewVec3(x, y, z float64) Vec3 { return Vec3{X: x, Y: y, Z: z} }

// Vec3Zero is the origin (0, 0, 0).
func Vec3Zero() Vec3 { return Vec3{} }

// Add returns the componentwise sum p + q.
func (p Vec3) Add(q Vec3) Vec3 { return Vec3{p.X + q.X, p.Y + q.Y, p.Z + q.Z} }

// Sub returns the componentwise difference p - q.
func (p Vec3) Sub(q Vec3) Vec3 { return Vec3{p.X - q.X, p.Y - q.Y, p.Z - q.Z} }

// Scale returns p scaled by the factor s.
func (p Vec3) Scale(s float64) Vec3 { return Vec3{p.X * s, p.Y * s, p.Z * s} }

// Div returns p divided componentwise by the nonzero scalar s.
func (p Vec3) Div(s float64) Vec3 { return Vec3{p.X / s, p.Y / s, p.Z / s} }

// Neg returns the vector reflected through the origin, -p.
func (p Vec3) Neg() Vec3 { return Vec3{-p.X, -p.Y, -p.Z} }

// Dot returns the dot product p . q.
func (p Vec3) Dot(q Vec3) float64 { return p.X*q.X + p.Y*q.Y + p.Z*q.Z }

// Cross returns the vector cross product p x q.
func (p Vec3) Cross(q Vec3) Vec3 {
	return Vec3{
		p.Y*q.Z - p.Z*q.Y,
		p.Z*q.X - p.X*q.Z,
		p.X*q.Y - p.Y*q.X,
	}
}

// NormSq returns the squared Euclidean norm |p|^2.
func (p Vec3) NormSq() float64 { return p.X*p.X + p.Y*p.Y + p.Z*p.Z }

// Norm returns the Euclidean norm |p|.
func (p Vec3) Norm() float64 { return math.Sqrt(p.NormSq()) }

// DistanceSq returns the squared Euclidean distance between p and q.
func (p Vec3) DistanceSq(q Vec3) float64 { return p.Sub(q).NormSq() }

// Distance returns the Euclidean distance between p and q.
func (p Vec3) Distance(q Vec3) float64 { return p.Sub(q).Norm() }

// Normalize returns the unit vector in the direction of p. The zero vector is
// returned unchanged.
func (p Vec3) Normalize() Vec3 {
	n := p.Norm()
	if n == 0 {
		return p
	}
	return p.Scale(1 / n)
}

// Lerp returns the linear interpolation (1-t)*p + t*q.
func (p Vec3) Lerp(q Vec3, t float64) Vec3 {
	return Vec3{
		p.X + (q.X-p.X)*t,
		p.Y + (q.Y-p.Y)*t,
		p.Z + (q.Z-p.Z)*t,
	}
}

// Midpoint returns the midpoint of the segment p-q.
func (p Vec3) Midpoint(q Vec3) Vec3 {
	return Vec3{(p.X + q.X) / 2, (p.Y + q.Y) / 2, (p.Z + q.Z) / 2}
}

// Min returns the componentwise minimum of p and q.
func (p Vec3) Min(q Vec3) Vec3 {
	return Vec3{math.Min(p.X, q.X), math.Min(p.Y, q.Y), math.Min(p.Z, q.Z)}
}

// Max returns the componentwise maximum of p and q.
func (p Vec3) Max(q Vec3) Vec3 {
	return Vec3{math.Max(p.X, q.X), math.Max(p.Y, q.Y), math.Max(p.Z, q.Z)}
}

// TripleProduct returns the scalar triple product p . (q x r), equal to the
// signed volume of the parallelepiped spanned by the three vectors.
func (p Vec3) TripleProduct(q, r Vec3) float64 { return p.Dot(q.Cross(r)) }

// AngleBetween returns the unsigned angle in radians between p and q, in
// the range [0, pi].
func (p Vec3) AngleBetween(q Vec3) float64 {
	d := p.Norm() * q.Norm()
	if d == 0 {
		return 0
	}
	c := p.Dot(q) / d
	if c > 1 {
		c = 1
	} else if c < -1 {
		c = -1
	}
	return math.Acos(c)
}

// Equal reports whether p and q are exactly equal.
func (p Vec3) Equal(q Vec3) bool { return p.X == q.X && p.Y == q.Y && p.Z == q.Z }

// ApproxEqual reports whether p and q agree to within absolute tolerance eps.
func (p Vec3) ApproxEqual(q Vec3, eps float64) bool {
	return math.Abs(p.X-q.X) <= eps &&
		math.Abs(p.Y-q.Y) <= eps &&
		math.Abs(p.Z-q.Z) <= eps
}

// IsFinite reports whether every coordinate of p is finite.
func (p Vec3) IsFinite() bool {
	for _, c := range [3]float64{p.X, p.Y, p.Z} {
		if math.IsInf(c, 0) || math.IsNaN(c) {
			return false
		}
	}
	return true
}

// String returns a human-readable representation of p.
func (p Vec3) String() string { return fmt.Sprintf("(%g, %g, %g)", p.X, p.Y, p.Z) }

// CentroidVec3 returns the arithmetic mean of the given points. It returns the
// origin for an empty slice.
func CentroidVec3(pts []Vec3) Vec3 {
	if len(pts) == 0 {
		return Vec3{}
	}
	var s Vec3
	for _, p := range pts {
		s = s.Add(p)
	}
	return s.Div(float64(len(pts)))
}

// BoundingBox3 returns the axis-aligned minimum and maximum corners enclosing
// the given points. It returns two zero vectors for an empty slice.
func BoundingBox3(pts []Vec3) (min, max Vec3) {
	if len(pts) == 0 {
		return Vec3{}, Vec3{}
	}
	min, max = pts[0], pts[0]
	for _, p := range pts[1:] {
		min = min.Min(p)
		max = max.Max(p)
	}
	return min, max
}
