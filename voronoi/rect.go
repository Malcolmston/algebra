package voronoi

import (
	"errors"
	"fmt"
	"math"
)

// ErrEmpty is returned by aggregate constructions when given no points.
var ErrEmpty = errors.New("voronoi: empty point set")

// Rect is an axis-aligned rectangle described by its minimum and maximum
// corners. It is well formed when Min.X <= Max.X and Min.Y <= Max.Y.
type Rect struct {
	Min, Max Point
}

// NewRect returns the axis-aligned rectangle spanning the two corners,
// normalizing so that Min holds the smaller coordinates.
func NewRect(a, b Point) Rect {
	return Rect{
		Min: Point{math.Min(a.X, b.X), math.Min(a.Y, b.Y)},
		Max: Point{math.Max(a.X, b.X), math.Max(a.Y, b.Y)},
	}
}

// BoundingBox returns the axis-aligned bounding box of the given points. It
// returns ErrEmpty when pts is empty.
func BoundingBox(pts []Point) (Rect, error) {
	if len(pts) == 0 {
		return Rect{}, ErrEmpty
	}
	mn := pts[0]
	mx := pts[0]
	for _, p := range pts[1:] {
		mn.X = math.Min(mn.X, p.X)
		mn.Y = math.Min(mn.Y, p.Y)
		mx.X = math.Max(mx.X, p.X)
		mx.Y = math.Max(mx.Y, p.Y)
	}
	return Rect{Min: mn, Max: mx}, nil
}

// Width returns the horizontal extent of the rectangle.
func (r Rect) Width() float64 { return r.Max.X - r.Min.X }

// Height returns the vertical extent of the rectangle.
func (r Rect) Height() float64 { return r.Max.Y - r.Min.Y }

// Area returns the area of the rectangle.
func (r Rect) Area() float64 { return r.Width() * r.Height() }

// Perimeter returns the perimeter of the rectangle.
func (r Rect) Perimeter() float64 { return 2 * (r.Width() + r.Height()) }

// Center returns the centre point of the rectangle.
func (r Rect) Center() Point {
	return Point{(r.Min.X + r.Max.X) / 2, (r.Min.Y + r.Max.Y) / 2}
}

// Contains reports whether p lies inside or on the rectangle, within tolerance
// eps.
func (r Rect) Contains(p Point, eps float64) bool {
	return p.X >= r.Min.X-eps && p.X <= r.Max.X+eps &&
		p.Y >= r.Min.Y-eps && p.Y <= r.Max.Y+eps
}

// Expand returns the rectangle grown outward by delta on every side. Negative
// delta shrinks it.
func (r Rect) Expand(delta float64) Rect {
	return Rect{
		Min: Point{r.Min.X - delta, r.Min.Y - delta},
		Max: Point{r.Max.X + delta, r.Max.Y + delta},
	}
}

// Corners returns the four corners of the rectangle in counterclockwise order
// starting from Min.
func (r Rect) Corners() [4]Point {
	return [4]Point{
		{r.Min.X, r.Min.Y},
		{r.Max.X, r.Min.Y},
		{r.Max.X, r.Max.Y},
		{r.Min.X, r.Max.Y},
	}
}

// Union returns the smallest rectangle containing both r and s.
func (r Rect) Union(s Rect) Rect {
	return Rect{
		Min: Point{math.Min(r.Min.X, s.Min.X), math.Min(r.Min.Y, s.Min.Y)},
		Max: Point{math.Max(r.Max.X, s.Max.X), math.Max(r.Max.Y, s.Max.Y)},
	}
}

// Diagonal returns the length of the rectangle's diagonal.
func (r Rect) Diagonal() float64 { return math.Hypot(r.Width(), r.Height()) }

// String returns a human-readable representation of the rectangle.
func (r Rect) String() string {
	return fmt.Sprintf("Rect{min: %s, max: %s}", r.Min, r.Max)
}
