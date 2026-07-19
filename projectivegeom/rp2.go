package projectivegeom

import (
	"errors"
	"math"
)

// ErrDegenerate is returned by constructions whose inputs are degenerate, for
// example joining two coincident points or meeting two identical lines.
var ErrDegenerate = errors.New("projectivegeom: degenerate configuration")

// ErrAtInfinity is returned when an affine value is requested from a point that
// lies on the line at infinity.
var ErrAtInfinity = errors.New("projectivegeom: point at infinity")

// Point is a point of the real projective plane RP^2 in homogeneous
// coordinates [X Y Z], defined up to a non-zero scalar. The finite point (x, y)
// corresponds to [x y 1]; points with Z = 0 lie on the line at infinity.
type Point struct {
	V Vec3
}

// Line is a line of RP^2 in homogeneous (dual) coordinates [A B C] denoting the
// locus A*X + B*Y + C*Z = 0, defined up to a non-zero scalar.
type Line struct {
	V Vec3
}

// NewPoint returns the homogeneous point [x y z].
func NewPoint(x, y, z float64) Point { return Point{Vec3{x, y, z}} }

// PointFromAffine returns the finite point [x y 1].
func PointFromAffine(x, y float64) Point { return Point{Vec3{x, y, 1}} }

// PointAtInfinity returns the point at infinity in the direction (dx, dy),
// namely [dx dy 0].
func PointAtInfinity(dx, dy float64) Point { return Point{Vec3{dx, dy, 0}} }

// NewLine returns the homogeneous line [a b c] for a*x + b*y + c = 0.
func NewLine(a, b, c float64) Line { return Line{Vec3{a, b, c}} }

// LineAtInfinity returns the line at infinity [0 0 1].
func LineAtInfinity() Line { return Line{Vec3{0, 0, 1}} }

// IsValid reports whether the point is a non-zero homogeneous triple.
func (p Point) IsValid() bool { return !p.V.IsZero(Eps) }

// IsValid reports whether the line is a non-zero homogeneous triple.
func (l Line) IsValid() bool { return !l.V.IsZero(Eps) }

// IsAtInfinity reports whether p lies on the line at infinity (Z within tol of
// zero after scaling by the largest component).
func (p Point) IsAtInfinity(tol float64) bool {
	s := p.V.MaxAbs()
	if s < Eps {
		return false
	}
	return math.Abs(p.V.Z) <= tol*s
}

// Affine returns the Euclidean coordinates (x, y) = (X/Z, Y/Z), or an error
// when p lies on the line at infinity.
func (p Point) Affine() (x, y float64, err error) {
	if p.IsAtInfinity(Eps) {
		return 0, 0, ErrAtInfinity
	}
	return p.V.X / p.V.Z, p.V.Y / p.V.Z, nil
}

// Normalized returns p rescaled so that Z = 1 when finite, or scaled to unit
// Euclidean length when at infinity. The second result is false only for the
// zero triple.
func (p Point) Normalized() (Point, bool) {
	if !p.IsValid() {
		return Point{}, false
	}
	if math.Abs(p.V.Z) > Eps*p.V.MaxAbs() {
		return Point{p.V.Scale(1 / p.V.Z)}, true
	}
	u, _ := p.V.Normalized()
	return Point{u}, true
}

// Normalized returns l rescaled to unit Euclidean length. The second result is
// false only for the zero triple.
func (l Line) Normalized() (Line, bool) {
	u, ok := l.V.Normalized()
	return Line{u}, ok
}

// Equal reports whether p and q denote the same projective point, i.e. their
// homogeneous triples are parallel within tol.
func (p Point) Equal(q Point, tol float64) bool {
	sp, sq := p.V.MaxAbs(), q.V.MaxAbs()
	if sp < Eps || sq < Eps {
		return sp < Eps && sq < Eps
	}
	return p.V.Scale(1/sp).Parallel(q.V.Scale(1/sq), tol)
}

// Equal reports whether l and m denote the same projective line.
func (l Line) Equal(m Line, tol float64) bool {
	return Point{l.V}.Equal(Point{m.V}, tol)
}

// Join returns the line through the two points p and q. When p and q coincide
// the result is the zero triple and ErrDegenerate is returned.
func Join(p, q Point) (Line, error) {
	c := p.V.Cross(q.V)
	if c.IsZero(Eps) {
		return Line{}, ErrDegenerate
	}
	return Line{c}, nil
}

// Meet returns the intersection point of the two lines l and m. When l and m
// are the same line the result is the zero triple and ErrDegenerate is
// returned.
func Meet(l, m Line) (Point, error) {
	c := l.V.Cross(m.V)
	if c.IsZero(Eps) {
		return Point{}, ErrDegenerate
	}
	return Point{c}, nil
}

// LineThrough is an alias for Join reading naturally as "the line through p and
// q".
func LineThrough(p, q Point) (Line, error) { return Join(p, q) }

// Intersect is an alias for Meet reading naturally as "the intersection of l
// and m".
func Intersect(l, m Line) (Point, error) { return Meet(l, m) }

// Incidence returns the (scale-dependent) value p.V·l.V, which is zero exactly
// when p lies on l.
func Incidence(p Point, l Line) float64 { return p.V.Dot(l.V) }

// OnLine reports whether p lies on l within tol, using coordinates scaled by
// their largest magnitudes so the test is independent of homogeneous scale.
func OnLine(p Point, l Line, tol float64) bool {
	sp, sl := p.V.MaxAbs(), l.V.MaxAbs()
	if sp < Eps || sl < Eps {
		return false
	}
	return math.Abs(p.V.Dot(l.V)/(sp*sl)) <= tol
}

// Collinear reports whether the three points are collinear within tol, i.e. the
// determinant of their homogeneous coordinates vanishes after scale
// normalization.
func Collinear(p, q, r Point, tol float64) bool {
	sp, sq, sr := p.V.MaxAbs(), q.V.MaxAbs(), r.V.MaxAbs()
	if sp < Eps || sq < Eps || sr < Eps {
		return false
	}
	a := p.V.Scale(1 / sp)
	b := q.V.Scale(1 / sq)
	c := r.V.Scale(1 / sr)
	return math.Abs(a.Dot(b.Cross(c))) <= tol
}

// Concurrent reports whether the three lines pass through a common point within
// tol (the dual of Collinear).
func Concurrent(l, m, n Line, tol float64) bool {
	return Collinear(Point{l.V}, Point{m.V}, Point{n.V}, tol)
}

// Triangle holds three points assumed to be in general position (no two equal,
// not collinear). Its sides are the joins of opposite vertex pairs.
type Triangle struct {
	A, B, C Point
}

// Sides returns the three side lines BC, CA, AB of the triangle. It returns
// ErrDegenerate when the vertices are not in general position.
func (t Triangle) Sides() (bc, ca, ab Line, err error) {
	if bc, err = Join(t.B, t.C); err != nil {
		return
	}
	if ca, err = Join(t.C, t.A); err != nil {
		return
	}
	ab, err = Join(t.A, t.B)
	return
}

// IsDegenerate reports whether the three vertices are collinear (or coincident)
// within tol.
func (t Triangle) IsDegenerate(tol float64) bool {
	return Collinear(t.A, t.B, t.C, tol)
}

// Centroid returns the affine centroid of the triangle's finite vertices. It
// returns ErrAtInfinity if any vertex lies at infinity.
func (t Triangle) Centroid() (Point, error) {
	ax, ay, err := t.A.Affine()
	if err != nil {
		return Point{}, err
	}
	bx, by, err := t.B.Affine()
	if err != nil {
		return Point{}, err
	}
	cx, cy, err := t.C.Affine()
	if err != nil {
		return Point{}, err
	}
	return PointFromAffine((ax+bx+cx)/3, (ay+by+cy)/3), nil
}

// LineAngle returns the direction angle of the finite line [a b c] in radians,
// in the range (-pi/2, pi/2], measured as atan2(-a, b) folded into that range.
func (l Line) LineAngle() float64 {
	ang := math.Atan2(-l.V.X, l.V.Y)
	for ang > math.Pi/2 {
		ang -= math.Pi
	}
	for ang <= -math.Pi/2 {
		ang += math.Pi
	}
	return ang
}

// ParallelLines reports whether l and m are parallel in the affine sense, i.e.
// they meet on the line at infinity, within tol.
func ParallelLines(l, m Line, tol float64) bool {
	p := l.V.Cross(m.V)
	if p.IsZero(Eps) {
		return true // identical lines are trivially parallel
	}
	return math.Abs(p.Z) <= tol*Vec3{p.X, p.Y, 0}.Norm()
}

// PerpendicularLines reports whether the finite lines l and m are perpendicular
// within tol (using their affine normals (a, b)).
func PerpendicularLines(l, m Line, tol float64) bool {
	return math.Abs(l.V.X*m.V.X+l.V.Y*m.V.Y) <=
		tol*math.Hypot(l.V.X, l.V.Y)*math.Hypot(m.V.X, m.V.Y)
}
