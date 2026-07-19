package hull3d

import (
	"errors"
	"fmt"
	"math"
	"sort"
)

// Vec3 is a three-dimensional vector over float64, used interchangeably as a
// point or a displacement in R^3.
type Vec3 struct {
	X, Y, Z float64
}

// NewVec3 returns the vector with the given components.
func NewVec3(x, y, z float64) Vec3 { return Vec3{x, y, z} }

// Zero returns the zero vector (0, 0, 0).
func Zero() Vec3 { return Vec3{} }

// UnitX returns the first standard basis vector (1, 0, 0).
func UnitX() Vec3 { return Vec3{1, 0, 0} }

// UnitY returns the second standard basis vector (0, 1, 0).
func UnitY() Vec3 { return Vec3{0, 1, 0} }

// UnitZ returns the third standard basis vector (0, 0, 1).
func UnitZ() Vec3 { return Vec3{0, 0, 1} }

// Splat returns the vector with every component equal to s.
func Splat(s float64) Vec3 { return Vec3{s, s, s} }

// Add returns the component-wise sum a + b.
func (a Vec3) Add(b Vec3) Vec3 { return Vec3{a.X + b.X, a.Y + b.Y, a.Z + b.Z} }

// Sub returns the component-wise difference a - b.
func (a Vec3) Sub(b Vec3) Vec3 { return Vec3{a.X - b.X, a.Y - b.Y, a.Z - b.Z} }

// Neg returns the additive inverse -a.
func (a Vec3) Neg() Vec3 { return Vec3{-a.X, -a.Y, -a.Z} }

// Scale returns the scalar multiple s*a.
func (a Vec3) Scale(s float64) Vec3 { return Vec3{a.X * s, a.Y * s, a.Z * s} }

// Mul returns the component-wise (Hadamard) product of a and b.
func (a Vec3) Mul(b Vec3) Vec3 { return Vec3{a.X * b.X, a.Y * b.Y, a.Z * b.Z} }

// Div returns the component-wise quotient a / b. It returns an error if any
// component of b is zero.
func (a Vec3) Div(b Vec3) (Vec3, error) {
	if b.X == 0 || b.Y == 0 || b.Z == 0 {
		return Vec3{}, errors.New("hull3d: division by zero component")
	}
	return Vec3{a.X / b.X, a.Y / b.Y, a.Z / b.Z}, nil
}

// Dot returns the Euclidean inner product a·b.
func (a Vec3) Dot(b Vec3) float64 { return a.X*b.X + a.Y*b.Y + a.Z*b.Z }

// Cross returns the vector cross product a×b.
func (a Vec3) Cross(b Vec3) Vec3 {
	return Vec3{
		a.Y*b.Z - a.Z*b.Y,
		a.Z*b.X - a.X*b.Z,
		a.X*b.Y - a.Y*b.X,
	}
}

// Triple returns the scalar triple product a·(b×c), the signed volume of the
// parallelepiped spanned by a, b and c.
func Triple(a, b, c Vec3) float64 { return a.Dot(b.Cross(c)) }

// LengthSq returns the squared Euclidean length of a.
func (a Vec3) LengthSq() float64 { return a.Dot(a) }

// Length returns the Euclidean length of a.
func (a Vec3) Length() float64 { return math.Sqrt(a.LengthSq()) }

// Norm is an alias for Length.
func (a Vec3) Norm() float64 { return a.Length() }

// L1 returns the Manhattan (taxicab) norm |x|+|y|+|z|.
func (a Vec3) L1() float64 { return math.Abs(a.X) + math.Abs(a.Y) + math.Abs(a.Z) }

// LInf returns the Chebyshev (maximum) norm max(|x|,|y|,|z|).
func (a Vec3) LInf() float64 {
	return math.Max(math.Abs(a.X), math.Max(math.Abs(a.Y), math.Abs(a.Z)))
}

// DistanceSq returns the squared Euclidean distance between a and b.
func (a Vec3) DistanceSq(b Vec3) float64 { return a.Sub(b).LengthSq() }

// Distance returns the Euclidean distance between a and b.
func (a Vec3) Distance(b Vec3) float64 { return a.Sub(b).Length() }

// Normalize returns a unit vector in the direction of a. It returns an error if
// a is the zero vector (has length below the machine epsilon).
func (a Vec3) Normalize() (Vec3, error) {
	n := a.Length()
	if n < 1e-300 || math.IsNaN(n) {
		return Vec3{}, errors.New("hull3d: cannot normalize zero-length vector")
	}
	return a.Scale(1 / n), nil
}

// NormalizeOr returns a unit vector in the direction of a, or the fallback
// vector if a has effectively zero length.
func (a Vec3) NormalizeOr(fallback Vec3) Vec3 {
	if u, err := a.Normalize(); err == nil {
		return u
	}
	return fallback
}

// Lerp returns the linear interpolation (1-t)*a + t*b.
func (a Vec3) Lerp(b Vec3, t float64) Vec3 { return a.Add(b.Sub(a).Scale(t)) }

// Midpoint returns the arithmetic mean of a and b.
func (a Vec3) Midpoint(b Vec3) Vec3 { return a.Add(b).Scale(0.5) }

// Abs returns the component-wise absolute value of a.
func (a Vec3) Abs() Vec3 { return Vec3{math.Abs(a.X), math.Abs(a.Y), math.Abs(a.Z)} }

// Min returns the component-wise minimum of a and b.
func (a Vec3) Min(b Vec3) Vec3 {
	return Vec3{math.Min(a.X, b.X), math.Min(a.Y, b.Y), math.Min(a.Z, b.Z)}
}

// Max returns the component-wise maximum of a and b.
func (a Vec3) Max(b Vec3) Vec3 {
	return Vec3{math.Max(a.X, b.X), math.Max(a.Y, b.Y), math.Max(a.Z, b.Z)}
}

// MinComponent returns the smallest of the three components of a.
func (a Vec3) MinComponent() float64 { return math.Min(a.X, math.Min(a.Y, a.Z)) }

// MaxComponent returns the largest of the three components of a.
func (a Vec3) MaxComponent() float64 { return math.Max(a.X, math.Max(a.Y, a.Z)) }

// MaxAxis returns the index (0, 1 or 2) of the component of a with the largest
// absolute value.
func (a Vec3) MaxAxis() int {
	ax, ay, az := math.Abs(a.X), math.Abs(a.Y), math.Abs(a.Z)
	if ax >= ay && ax >= az {
		return 0
	}
	if ay >= az {
		return 1
	}
	return 2
}

// Get returns component i of a for i in {0,1,2}. It panics for other indices,
// matching the convention of indexing a fixed-length container.
func (a Vec3) Get(i int) float64 {
	switch i {
	case 0:
		return a.X
	case 1:
		return a.Y
	case 2:
		return a.Z
	}
	panic("hull3d: Vec3 index out of range")
}

// With returns a copy of a with component i replaced by v.
func (a Vec3) With(i int, v float64) Vec3 {
	switch i {
	case 0:
		a.X = v
	case 1:
		a.Y = v
	case 2:
		a.Z = v
	default:
		panic("hull3d: Vec3 index out of range")
	}
	return a
}

// IsZero reports whether every component of a is exactly zero.
func (a Vec3) IsZero() bool { return a.X == 0 && a.Y == 0 && a.Z == 0 }

// IsFinite reports whether all components of a are finite (neither Inf nor NaN).
func (a Vec3) IsFinite() bool {
	return !math.IsInf(a.X, 0) && !math.IsInf(a.Y, 0) && !math.IsInf(a.Z, 0) &&
		!math.IsNaN(a.X) && !math.IsNaN(a.Y) && !math.IsNaN(a.Z)
}

// ApproxEqual reports whether a and b agree to within an absolute tolerance eps
// in every component.
func (a Vec3) ApproxEqual(b Vec3, eps float64) bool {
	return math.Abs(a.X-b.X) <= eps && math.Abs(a.Y-b.Y) <= eps && math.Abs(a.Z-b.Z) <= eps
}

// Equal reports whether a and b are exactly equal component-wise.
func (a Vec3) Equal(b Vec3) bool { return a == b }

// AngleBetween returns the unsigned angle in radians between a and b, in the
// range [0, pi]. It returns an error if either vector has zero length.
func (a Vec3) AngleBetween(b Vec3) (float64, error) {
	la, lb := a.Length(), b.Length()
	if la < 1e-300 || lb < 1e-300 {
		return 0, errors.New("hull3d: angle undefined for zero-length vector")
	}
	c := a.Dot(b) / (la * lb)
	c = math.Max(-1, math.Min(1, c))
	return math.Acos(c), nil
}

// ProjectOnto returns the orthogonal projection of a onto the direction b. It
// returns an error if b is the zero vector.
func (a Vec3) ProjectOnto(b Vec3) (Vec3, error) {
	d := b.LengthSq()
	if d < 1e-300 {
		return Vec3{}, errors.New("hull3d: cannot project onto zero vector")
	}
	return b.Scale(a.Dot(b) / d), nil
}

// RejectFrom returns the component of a orthogonal to b (a minus its projection
// onto b). It returns an error if b is the zero vector.
func (a Vec3) RejectFrom(b Vec3) (Vec3, error) {
	p, err := a.ProjectOnto(b)
	if err != nil {
		return Vec3{}, err
	}
	return a.Sub(p), nil
}

// Reflect returns a reflected across the plane through the origin with unit
// normal n. If n is not unit length the result is scaled accordingly.
func (a Vec3) Reflect(n Vec3) Vec3 { return a.Sub(n.Scale(2 * a.Dot(n))) }

// AnyPerpendicular returns a non-zero vector orthogonal to a. It returns an
// error if a is the zero vector.
func (a Vec3) AnyPerpendicular() (Vec3, error) {
	if a.IsZero() {
		return Vec3{}, errors.New("hull3d: no perpendicular to zero vector")
	}
	if math.Abs(a.X) <= math.Abs(a.Y) && math.Abs(a.X) <= math.Abs(a.Z) {
		return a.Cross(UnitX()), nil
	}
	if math.Abs(a.Y) <= math.Abs(a.Z) {
		return a.Cross(UnitY()), nil
	}
	return a.Cross(UnitZ()), nil
}

// OrthonormalBasis returns two unit vectors u, v that together with a
// normalized copy of a form a right-handed orthonormal basis. It returns an
// error if a is the zero vector.
func (a Vec3) OrthonormalBasis() (u, v Vec3, err error) {
	w, err := a.Normalize()
	if err != nil {
		return Vec3{}, Vec3{}, err
	}
	p, _ := a.AnyPerpendicular()
	u, _ = p.Normalize()
	v = w.Cross(u)
	return u, v, nil
}

// Clamp returns a with each component clamped to the box [lo, hi].
func (a Vec3) Clamp(lo, hi Vec3) Vec3 {
	return Vec3{
		math.Max(lo.X, math.Min(hi.X, a.X)),
		math.Max(lo.Y, math.Min(hi.Y, a.Y)),
		math.Max(lo.Z, math.Min(hi.Z, a.Z)),
	}
}

// String renders a in a compact parenthesised form.
func (a Vec3) String() string { return fmt.Sprintf("(%g, %g, %g)", a.X, a.Y, a.Z) }

// Centroid returns the arithmetic mean of the given points. It returns the zero
// vector for an empty slice.
func Centroid(pts []Vec3) Vec3 {
	if len(pts) == 0 {
		return Vec3{}
	}
	var s Vec3
	for _, p := range pts {
		s = s.Add(p)
	}
	return s.Scale(1 / float64(len(pts)))
}

// SumVec3 returns the sum of the given points.
func SumVec3(pts []Vec3) Vec3 {
	var s Vec3
	for _, p := range pts {
		s = s.Add(p)
	}
	return s
}

// BoundingBox returns the axis-aligned minimum and maximum corners enclosing
// the given points. It returns an error for an empty slice.
func BoundingBox(pts []Vec3) (min, max Vec3, err error) {
	if len(pts) == 0 {
		return Vec3{}, Vec3{}, errors.New("hull3d: bounding box of empty set")
	}
	min, max = pts[0], pts[0]
	for _, p := range pts[1:] {
		min = min.Min(p)
		max = max.Max(p)
	}
	return min, max, nil
}

// FarthestPoint returns the index and value of the point in pts that maximises
// the dot product with direction d (the support point in direction d). It
// returns an error for an empty slice.
func FarthestPoint(pts []Vec3, d Vec3) (int, Vec3, error) {
	if len(pts) == 0 {
		return -1, Vec3{}, errors.New("hull3d: farthest point of empty set")
	}
	best := 0
	bestDot := pts[0].Dot(d)
	for i := 1; i < len(pts); i++ {
		if v := pts[i].Dot(d); v > bestDot {
			bestDot, best = v, i
		}
	}
	return best, pts[best], nil
}

// SortByAxis returns a copy of pts sorted ascending by component axis (0, 1 or
// 2). The input slice is not modified.
func SortByAxis(pts []Vec3, axis int) []Vec3 {
	out := append([]Vec3(nil), pts...)
	sort.SliceStable(out, func(i, j int) bool { return out[i].Get(axis) < out[j].Get(axis) })
	return out
}

// DedupPoints returns the points of pts with near-duplicates (within Euclidean
// distance eps) removed, preserving first-seen order. It is O(n^2) and intended
// for the modest point counts typical of hull construction.
func DedupPoints(pts []Vec3, eps float64) []Vec3 {
	var out []Vec3
	for _, p := range pts {
		dup := false
		for _, q := range out {
			if p.DistanceSq(q) <= eps*eps {
				dup = true
				break
			}
		}
		if !dup {
			out = append(out, p)
		}
	}
	return out
}
