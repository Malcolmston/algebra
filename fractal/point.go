package fractal

import "math"

// Point2D is a point in the Euclidean plane with coordinates X and Y.
type Point2D struct {
	X, Y float64
}

// Add returns the vector sum p+q.
func (p Point2D) Add(q Point2D) Point2D {
	return Point2D{p.X + q.X, p.Y + q.Y}
}

// Sub returns the vector difference p-q.
func (p Point2D) Sub(q Point2D) Point2D {
	return Point2D{p.X - q.X, p.Y - q.Y}
}

// Scale returns p with both coordinates multiplied by s.
func (p Point2D) Scale(s float64) Point2D {
	return Point2D{p.X * s, p.Y * s}
}

// Rotate returns p rotated counterclockwise about the origin by angle radians.
func (p Point2D) Rotate(angle float64) Point2D {
	sn, cs := math.Sincos(angle)
	return Point2D{p.X*cs - p.Y*sn, p.X*sn + p.Y*cs}
}

// Dist returns the Euclidean distance between p and q.
func (p Point2D) Dist(q Point2D) float64 {
	return math.Hypot(p.X-q.X, p.Y-q.Y)
}

// Norm returns the Euclidean length of p treated as a vector from the origin.
func (p Point2D) Norm() float64 {
	return math.Hypot(p.X, p.Y)
}

// Midpoint returns the point halfway between p and q.
func (p Point2D) Midpoint(q Point2D) Point2D {
	return Point2D{(p.X + q.X) / 2, (p.Y + q.Y) / 2}
}

// Lerp returns the linear interpolation between p (t=0) and q (t=1).
func (p Point2D) Lerp(q Point2D, t float64) Point2D {
	return Point2D{p.X + (q.X-p.X)*t, p.Y + (q.Y-p.Y)*t}
}

// BoundingBox returns the axis-aligned bounding box of pts as its minimum and
// maximum corners. It returns two zero points when pts is empty.
func BoundingBox(pts []Point2D) (min, max Point2D) {
	if len(pts) == 0 {
		return Point2D{}, Point2D{}
	}
	min, max = pts[0], pts[0]
	for _, p := range pts[1:] {
		if p.X < min.X {
			min.X = p.X
		}
		if p.Y < min.Y {
			min.Y = p.Y
		}
		if p.X > max.X {
			max.X = p.X
		}
		if p.Y > max.Y {
			max.Y = p.Y
		}
	}
	return min, max
}
