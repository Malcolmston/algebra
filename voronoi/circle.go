package voronoi

import (
	"fmt"
	"math"
)

// Circle is a circle defined by its centre and radius.
type Circle struct {
	Center Point
	Radius float64
}

// NewCircle returns the circle with the given centre and radius.
func NewCircle(center Point, radius float64) Circle {
	return Circle{Center: center, Radius: radius}
}

// CircleThrough returns the unique circle passing through the three points.
// It returns ErrCollinear if the points are collinear.
func CircleThrough(a, b, c Point) (Circle, error) {
	return Circumcircle(a, b, c)
}

// CircleFromDiameter returns the circle whose diameter is the segment ab.
func CircleFromDiameter(a, b Point) Circle {
	return Circle{Center: a.Midpoint(b), Radius: a.Distance(b) / 2}
}

// Area returns the area enclosed by the circle.
func (c Circle) Area() float64 { return math.Pi * c.Radius * c.Radius }

// Circumference returns the circle's circumference.
func (c Circle) Circumference() float64 { return 2 * math.Pi * c.Radius }

// Diameter returns the circle's diameter.
func (c Circle) Diameter() float64 { return 2 * c.Radius }

// Contains reports whether p lies inside or on the circle, within tolerance
// eps applied to the radius.
func (c Circle) Contains(p Point, eps float64) bool {
	return c.Center.Distance(p) <= c.Radius+eps
}

// ContainsStrict reports whether p lies strictly inside the circle, within
// tolerance eps.
func (c Circle) ContainsStrict(p Point, eps float64) bool {
	return c.Center.Distance(p) < c.Radius-eps
}

// OnBoundary reports whether p lies on the circle boundary within tolerance
// eps.
func (c Circle) OnBoundary(p Point, eps float64) bool {
	return math.Abs(c.Center.Distance(p)-c.Radius) <= eps
}

// PointAt returns the point on the circle at angle theta (radians) measured
// counterclockwise from the positive X axis.
func (c Circle) PointAt(theta float64) Point {
	s, co := math.Sincos(theta)
	return Point{c.Center.X + c.Radius*co, c.Center.Y + c.Radius*s}
}

// Sample returns n points spaced evenly around the circle, starting at angle 0.
// It returns nil for n <= 0.
func (c Circle) Sample(n int) []Point {
	if n <= 0 {
		return nil
	}
	pts := make([]Point, n)
	for i := 0; i < n; i++ {
		pts[i] = c.PointAt(2 * math.Pi * float64(i) / float64(n))
	}
	return pts
}

// Intersects reports whether two circles overlap or touch.
func (c Circle) Intersects(o Circle) bool {
	d := c.Center.Distance(o.Center)
	return d <= c.Radius+o.Radius && d >= math.Abs(c.Radius-o.Radius)
}

// BoundingBox returns the axis-aligned bounding box of the circle.
func (c Circle) BoundingBox() Rect {
	return Rect{
		Min: Point{c.Center.X - c.Radius, c.Center.Y - c.Radius},
		Max: Point{c.Center.X + c.Radius, c.Center.Y + c.Radius},
	}
}

// String returns a human-readable representation of the circle.
func (c Circle) String() string {
	return fmt.Sprintf("Circle{center: %s, r: %g}", c.Center, c.Radius)
}

// MinimumEnclosingCircle returns the smallest circle that contains all of the
// given points. It uses Welzl's randomized-free incremental algorithm with a
// deterministic point order, running in expected linear time. It returns a
// zero-radius circle for a single point and the empty Circle for no points.
func MinimumEnclosingCircle(pts []Point) Circle {
	if len(pts) == 0 {
		return Circle{}
	}
	p := ClonePoints(pts)
	var c Circle
	c = Circle{Center: p[0], Radius: 0}
	for i := 1; i < len(p); i++ {
		if c.Contains(p[i], Eps) {
			continue
		}
		c = Circle{Center: p[i], Radius: 0}
		for j := 0; j < i; j++ {
			if c.Contains(p[j], Eps) {
				continue
			}
			c = CircleFromDiameter(p[i], p[j])
			for k := 0; k < j; k++ {
				if c.Contains(p[k], Eps) {
					continue
				}
				if cc, err := Circumcircle(p[i], p[j], p[k]); err == nil {
					c = cc
				}
			}
		}
	}
	return c
}

// EnclosingCircleRadius returns the radius of the minimum enclosing circle of
// the given points.
func EnclosingCircleRadius(pts []Point) float64 {
	return MinimumEnclosingCircle(pts).Radius
}
