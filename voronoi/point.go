package voronoi

import (
	"fmt"
	"math"
	"sort"
)

// Point is a location in the Euclidean plane.
type Point struct {
	X, Y float64
}

// NewPoint returns the Point with the given coordinates.
func NewPoint(x, y float64) Point { return Point{X: x, Y: y} }

// Origin is the point (0, 0).
func Origin() Point { return Point{0, 0} }

// Add returns the componentwise sum p + q.
func (p Point) Add(q Point) Point { return Point{p.X + q.X, p.Y + q.Y} }

// Sub returns the componentwise difference p - q.
func (p Point) Sub(q Point) Point { return Point{p.X - q.X, p.Y - q.Y} }

// Scale returns p scaled by the factor s.
func (p Point) Scale(s float64) Point { return Point{p.X * s, p.Y * s} }

// Neg returns the point reflected through the origin, -p.
func (p Point) Neg() Point { return Point{-p.X, -p.Y} }

// Dot returns the dot product p . q, treating both points as vectors.
func (p Point) Dot(q Point) float64 { return p.X*q.X + p.Y*q.Y }

// Cross returns the scalar (2D) cross product p x q = p.X*q.Y - p.Y*q.X.
func (p Point) Cross(q Point) float64 { return p.X*q.Y - p.Y*q.X }

// NormSq returns the squared Euclidean norm |p|^2.
func (p Point) NormSq() float64 { return p.X*p.X + p.Y*p.Y }

// Norm returns the Euclidean norm |p|.
func (p Point) Norm() float64 { return math.Hypot(p.X, p.Y) }

// DistanceSq returns the squared Euclidean distance between p and q.
func (p Point) DistanceSq(q Point) float64 {
	dx := p.X - q.X
	dy := p.Y - q.Y
	return dx*dx + dy*dy
}

// Distance returns the Euclidean distance between p and q.
func (p Point) Distance(q Point) float64 {
	return math.Hypot(p.X-q.X, p.Y-q.Y)
}

// Midpoint returns the point halfway between p and q.
func (p Point) Midpoint(q Point) Point {
	return Point{(p.X + q.X) / 2, (p.Y + q.Y) / 2}
}

// Lerp returns the linear interpolation (1-t)*p + t*q.
func (p Point) Lerp(q Point, t float64) Point {
	return Point{p.X + t*(q.X-p.X), p.Y + t*(q.Y-p.Y)}
}

// Perp returns p rotated by +90 degrees, i.e. the vector (-p.Y, p.X).
func (p Point) Perp() Point { return Point{-p.Y, p.X} }

// Angle returns the angle of the vector p measured from the positive X axis,
// in the range (-pi, pi].
func (p Point) Angle() float64 { return math.Atan2(p.Y, p.X) }

// AngleTo returns the angle of the vector q - p from the positive X axis.
func (p Point) AngleTo(q Point) float64 {
	return math.Atan2(q.Y-p.Y, q.X-p.X)
}

// Rotate returns p rotated about the origin by theta radians (counterclockwise).
func (p Point) Rotate(theta float64) Point {
	s, c := math.Sincos(theta)
	return Point{p.X*c - p.Y*s, p.X*s + p.Y*c}
}

// RotateAbout returns p rotated by theta radians about the centre c.
func (p Point) RotateAbout(c Point, theta float64) Point {
	return p.Sub(c).Rotate(theta).Add(c)
}

// Normalize returns the unit vector in the direction of p. If p is the zero
// vector it is returned unchanged.
func (p Point) Normalize() Point {
	n := p.Norm()
	if n == 0 {
		return p
	}
	return Point{p.X / n, p.Y / n}
}

// Equal reports whether p and q are exactly equal.
func (p Point) Equal(q Point) bool { return p.X == q.X && p.Y == q.Y }

// ApproxEqual reports whether p and q are within eps (per coordinate) of each
// other.
func (p Point) ApproxEqual(q Point, eps float64) bool {
	return math.Abs(p.X-q.X) <= eps && math.Abs(p.Y-q.Y) <= eps
}

// IsFinite reports whether both coordinates of p are finite.
func (p Point) IsFinite() bool {
	return !math.IsInf(p.X, 0) && !math.IsInf(p.Y, 0) &&
		!math.IsNaN(p.X) && !math.IsNaN(p.Y)
}

// String returns a human-readable representation of p.
func (p Point) String() string {
	return fmt.Sprintf("(%g, %g)", p.X, p.Y)
}

// Centroid returns the arithmetic mean of the given points. It returns the
// origin for an empty slice.
func Centroid(pts []Point) Point {
	if len(pts) == 0 {
		return Point{}
	}
	var sx, sy float64
	for _, p := range pts {
		sx += p.X
		sy += p.Y
	}
	n := float64(len(pts))
	return Point{sx / n, sy / n}
}

// DedupePoints returns the points with exact duplicates removed, preserving the
// order of first appearance.
func DedupePoints(pts []Point) []Point {
	seen := make(map[Point]struct{}, len(pts))
	out := make([]Point, 0, len(pts))
	for _, p := range pts {
		if _, ok := seen[p]; ok {
			continue
		}
		seen[p] = struct{}{}
		out = append(out, p)
	}
	return out
}

// SortPointsXY returns a copy of pts sorted lexicographically by X then Y.
func SortPointsXY(pts []Point) []Point {
	out := make([]Point, len(pts))
	copy(out, pts)
	sort.Slice(out, func(i, j int) bool {
		if out[i].X != out[j].X {
			return out[i].X < out[j].X
		}
		return out[i].Y < out[j].Y
	})
	return out
}

// ClonePoints returns a fresh copy of the given slice.
func ClonePoints(pts []Point) []Point {
	out := make([]Point, len(pts))
	copy(out, pts)
	return out
}
