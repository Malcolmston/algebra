package tilings

import (
	"math"
)

// Epsilon is the default absolute tolerance used by approximate comparisons
// throughout the package.
const Epsilon = 1e-9

// Phi is the golden ratio (1+sqrt(5))/2, the inflation factor of the Penrose
// substitution tilings.
const Phi = 1.618033988749894848204586834365638

// Deg2Rad converts an angle in degrees to radians.
func Deg2Rad(deg float64) float64 { return deg * math.Pi / 180 }

// Rad2Deg converts an angle in radians to degrees.
func Rad2Deg(rad float64) float64 { return rad * 180 / math.Pi }

// Point is a point in the plane, also used as a two-dimensional vector.
type Point struct {
	X, Y float64
}

// NewPoint returns the point (x, y).
func NewPoint(x, y float64) Point { return Point{x, y} }

// Polar returns the point at radius r and angle theta (radians) from the origin.
func Polar(r, theta float64) Point {
	return Point{r * math.Cos(theta), r * math.Sin(theta)}
}

// Origin returns the point (0, 0).
func Origin() Point { return Point{0, 0} }

// Add returns the vector sum p+q.
func (p Point) Add(q Point) Point { return Point{p.X + q.X, p.Y + q.Y} }

// Sub returns the vector difference p-q.
func (p Point) Sub(q Point) Point { return Point{p.X - q.X, p.Y - q.Y} }

// Neg returns the additive inverse -p.
func (p Point) Neg() Point { return Point{-p.X, -p.Y} }

// Scale returns p scaled by the factor s.
func (p Point) Scale(s float64) Point { return Point{p.X * s, p.Y * s} }

// Dot returns the Euclidean dot product p·q.
func (p Point) Dot(q Point) float64 { return p.X*q.X + p.Y*q.Y }

// Cross returns the scalar cross product p.X*q.Y - p.Y*q.X.
func (p Point) Cross(q Point) float64 { return p.X*q.Y - p.Y*q.X }

// Norm returns the Euclidean length of p.
func (p Point) Norm() float64 { return math.Hypot(p.X, p.Y) }

// Norm2 returns the squared Euclidean length of p.
func (p Point) Norm2() float64 { return p.X*p.X + p.Y*p.Y }

// Normalize returns the unit vector in the direction of p; the zero vector is
// returned unchanged.
func (p Point) Normalize() Point {
	n := p.Norm()
	if n == 0 {
		return p
	}
	return Point{p.X / n, p.Y / n}
}

// Perp returns p rotated by +90 degrees, i.e. the left-hand perpendicular.
func (p Point) Perp() Point { return Point{-p.Y, p.X} }

// Angle returns the direction of p measured counter-clockwise from the positive
// x-axis, in radians in (-pi, pi].
func (p Point) Angle() float64 { return math.Atan2(p.Y, p.X) }

// Rotate returns p rotated about the origin by theta radians.
func (p Point) Rotate(theta float64) Point {
	c, s := math.Cos(theta), math.Sin(theta)
	return Point{p.X*c - p.Y*s, p.X*s + p.Y*c}
}

// RotateAbout returns p rotated about the point c by theta radians.
func (p Point) RotateAbout(c Point, theta float64) Point {
	return p.Sub(c).Rotate(theta).Add(c)
}

// Distance returns the Euclidean distance between p and q.
func (p Point) Distance(q Point) float64 { return p.Sub(q).Norm() }

// Distance2 returns the squared Euclidean distance between p and q.
func (p Point) Distance2(q Point) float64 { return p.Sub(q).Norm2() }

// Midpoint returns the midpoint of p and q.
func (p Point) Midpoint(q Point) Point {
	return Point{(p.X + q.X) / 2, (p.Y + q.Y) / 2}
}

// Lerp returns the linear interpolation (1-t)*p + t*q.
func (p Point) Lerp(q Point, t float64) Point {
	return Point{p.X + (q.X-p.X)*t, p.Y + (q.Y-p.Y)*t}
}

// IsZero reports whether p is within eps of the origin.
func (p Point) IsZero(eps float64) bool { return p.Norm() <= eps }

// ApproxEqual reports whether p and q are within eps of one another.
func (p Point) ApproxEqual(q Point, eps float64) bool {
	return p.Sub(q).Norm() <= eps
}

// AngleBetween returns the unsigned angle between vectors u and v in [0, pi].
func AngleBetween(u, v Point) float64 {
	nu, nv := u.Norm(), v.Norm()
	if nu == 0 || nv == 0 {
		return 0
	}
	c := u.Dot(v) / (nu * nv)
	if c > 1 {
		c = 1
	} else if c < -1 {
		c = -1
	}
	return math.Acos(c)
}

// SignedAngle returns the signed angle from u to v in (-pi, pi], positive for a
// counter-clockwise turn.
func SignedAngle(u, v Point) float64 {
	return math.Atan2(u.Cross(v), u.Dot(v))
}

// Centroid returns the arithmetic mean of the given points; the empty input
// yields the origin.
func Centroid(pts ...Point) Point {
	if len(pts) == 0 {
		return Point{}
	}
	var s Point
	for _, p := range pts {
		s = s.Add(p)
	}
	return s.Scale(1 / float64(len(pts)))
}

// PolygonArea returns the signed area of the polygon with the given vertices in
// order; a positive value indicates counter-clockwise orientation.
func PolygonArea(vertices []Point) float64 {
	n := len(vertices)
	if n < 3 {
		return 0
	}
	var a float64
	for i := 0; i < n; i++ {
		j := (i + 1) % n
		a += vertices[i].Cross(vertices[j])
	}
	return a / 2
}

// PolygonCentroid returns the area centroid of the simple polygon with the
// given vertices in order; for a degenerate (zero-area) polygon it returns the
// vertex centroid.
func PolygonCentroid(vertices []Point) Point {
	n := len(vertices)
	if n == 0 {
		return Point{}
	}
	a := PolygonArea(vertices)
	if math.Abs(a) < Epsilon {
		return Centroid(vertices...)
	}
	var cx, cy float64
	for i := 0; i < n; i++ {
		j := (i + 1) % n
		cr := vertices[i].Cross(vertices[j])
		cx += (vertices[i].X + vertices[j].X) * cr
		cy += (vertices[i].Y + vertices[j].Y) * cr
	}
	return Point{cx / (6 * a), cy / (6 * a)}
}

// PolygonContains reports whether p lies inside the simple polygon with the
// given vertices, using the even-odd ray-casting rule. Points on the boundary
// may be reported either way.
func PolygonContains(vertices []Point, p Point) bool {
	n := len(vertices)
	if n < 3 {
		return false
	}
	inside := false
	j := n - 1
	for i := 0; i < n; i++ {
		vi, vj := vertices[i], vertices[j]
		if (vi.Y > p.Y) != (vj.Y > p.Y) {
			x := (vj.X-vi.X)*(p.Y-vi.Y)/(vj.Y-vi.Y) + vi.X
			if p.X < x {
				inside = !inside
			}
		}
		j = i
	}
	return inside
}

// RegularPolygon returns the vertices of a regular n-gon centred at c with
// circumradius r; the first vertex is placed at angle theta0 and vertices run
// counter-clockwise. It returns nil for n < 3.
func RegularPolygon(n int, c Point, r, theta0 float64) []Point {
	if n < 3 {
		return nil
	}
	out := make([]Point, n)
	for k := 0; k < n; k++ {
		a := theta0 + 2*math.Pi*float64(k)/float64(n)
		out[k] = c.Add(Polar(r, a))
	}
	return out
}

// NormalizeAngle returns theta reduced to the half-open interval [0, 2*pi).
func NormalizeAngle(theta float64) float64 {
	two := 2 * math.Pi
	t := math.Mod(theta, two)
	if t < 0 {
		t += two
	}
	return t
}

// NormalizeAngleSigned returns theta reduced to the interval (-pi, pi].
func NormalizeAngleSigned(theta float64) float64 {
	t := NormalizeAngle(theta)
	if t > math.Pi {
		t -= 2 * math.Pi
	}
	return t
}

// approxEqualScalar reports whether a and b are within eps.
func approxEqualScalar(a, b, eps float64) bool { return math.Abs(a-b) <= eps }
