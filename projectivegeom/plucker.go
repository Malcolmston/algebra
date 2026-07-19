package projectivegeom

import "math"

// LinePluecker holds the six Pluecker coordinates of a line in RP^3, in the
// order (P01, P02, P03, P12, P13, P23) where Pij = p_i*q_j - p_j*q_i for two
// points p, q on the line. The coordinates are defined up to a common non-zero
// scale and satisfy the Grassmann-Pluecker relation
// P01*P23 - P02*P13 + P03*P12 = 0.
type LinePluecker struct {
	P01, P02, P03, P12, P13, P23 float64
}

// LineFromPoints returns the Pluecker coordinates of the line through two
// distinct points of RP^3. It returns ErrDegenerate when the points coincide.
func LineFromPoints(p, q SPoint) (LinePluecker, error) {
	a := [4]float64{p.V.X, p.V.Y, p.V.Z, p.V.W}
	b := [4]float64{q.V.X, q.V.Y, q.V.Z, q.V.W}
	l := LinePluecker{
		P01: a[0]*b[1] - a[1]*b[0],
		P02: a[0]*b[2] - a[2]*b[0],
		P03: a[0]*b[3] - a[3]*b[0],
		P12: a[1]*b[2] - a[2]*b[1],
		P13: a[1]*b[3] - a[3]*b[1],
		P23: a[2]*b[3] - a[3]*b[2],
	}
	if l.MaxAbs() < Eps {
		return LinePluecker{}, ErrDegenerate
	}
	return l, nil
}

// LineFromPlanes returns the Pluecker coordinates of the line where two distinct
// planes meet, obtained from their dual coordinates by the Hodge duality. It
// returns ErrDegenerate when the planes coincide.
func LineFromPlanes(s, t SPlane) (LinePluecker, error) {
	a := [4]float64{s.V.X, s.V.Y, s.V.Z, s.V.W}
	b := [4]float64{t.V.X, t.V.Y, t.V.Z, t.V.W}
	d01 := a[0]*b[1] - a[1]*b[0]
	d02 := a[0]*b[2] - a[2]*b[0]
	d03 := a[0]*b[3] - a[3]*b[0]
	d12 := a[1]*b[2] - a[2]*b[1]
	d13 := a[1]*b[3] - a[3]*b[1]
	d23 := a[2]*b[3] - a[3]*b[2]
	l := LinePluecker{
		P01: d23,
		P02: -d13,
		P03: d12,
		P12: d03,
		P13: -d02,
		P23: d01,
	}
	if l.MaxAbs() < Eps {
		return LinePluecker{}, ErrDegenerate
	}
	return l, nil
}

// MaxAbs returns the largest absolute coordinate of the line.
func (l LinePluecker) MaxAbs() float64 {
	m := math.Abs(l.P01)
	for _, v := range []float64{l.P02, l.P03, l.P12, l.P13, l.P23} {
		if a := math.Abs(v); a > m {
			m = a
		}
	}
	return m
}

// GrassmannRelation returns P01*P23 - P02*P13 + P03*P12, which vanishes exactly
// for genuine lines.
func (l LinePluecker) GrassmannRelation() float64 {
	return l.P01*l.P23 - l.P02*l.P13 + l.P03*l.P12
}

// IsLine reports whether the coordinates satisfy the Grassmann-Pluecker relation
// within tol and are non-zero, hence represent an actual line.
func (l LinePluecker) IsLine(tol float64) bool {
	s := l.MaxAbs()
	if s < Eps {
		return false
	}
	return math.Abs(l.GrassmannRelation()) <= tol*s*s
}

// Direction returns a direction vector of the line, taken from its (X,W),(Y,W),
// (Z,W) coordinates. For a line at infinity the vector is zero.
func (l LinePluecker) Direction() Vec3 { return Vec3{l.P03, l.P13, l.P23} }

// Moment returns the moment vector of the line about the origin, from its
// spatial coordinates.
func (l LinePluecker) Moment() Vec3 { return Vec3{l.P12, -l.P02, l.P01} }

// ContainsPoint reports whether the point x lies on the line within tol, using
// the wedge condition x ∧ L = 0.
func (l LinePluecker) ContainsPoint(x SPoint, tol float64) bool {
	c := [4]float64{x.V.X, x.V.Y, x.V.Z, x.V.W}
	f := [4]float64{
		c[1]*l.P23 - c[2]*l.P13 + c[3]*l.P12,
		c[0]*l.P23 - c[2]*l.P03 + c[3]*l.P02,
		c[0]*l.P13 - c[1]*l.P03 + c[3]*l.P01,
		c[0]*l.P12 - c[1]*l.P02 + c[2]*l.P01,
	}
	scale := (1 + x.V.MaxAbs()) * (1 + l.MaxAbs())
	for _, v := range f {
		if math.Abs(v) > tol*scale {
			return false
		}
	}
	return true
}

// ReciprocalProduct returns the bilinear side operator of two lines,
// P01*Q23 + P23*Q01 - P02*Q13 - P13*Q02 + P03*Q12 + P12*Q03. It vanishes exactly
// when the two lines are coplanar (they meet, possibly at infinity).
func ReciprocalProduct(l, m LinePluecker) float64 {
	return l.P01*m.P23 + l.P23*m.P01 -
		l.P02*m.P13 - l.P13*m.P02 +
		l.P03*m.P12 + l.P12*m.P03
}

// LinesMeet reports whether two lines of RP^3 are coplanar (intersect) within
// tol.
func LinesMeet(l, m LinePluecker, tol float64) bool {
	scale := (1 + l.MaxAbs()) * (1 + m.MaxAbs())
	return math.Abs(ReciprocalProduct(l, m)) <= tol*scale
}

// Equal reports whether two lines have proportional Pluecker coordinates within
// tol.
func (l LinePluecker) Equal(m LinePluecker, tol float64) bool {
	a := []float64{l.P01, l.P02, l.P03, l.P12, l.P13, l.P23}
	b := []float64{m.P01, m.P02, m.P03, m.P12, m.P13, m.P23}
	na, nb := l.MaxAbs(), m.MaxAbs()
	if na < Eps || nb < Eps {
		return na < Eps && nb < Eps
	}
	for i := range a {
		for j := i + 1; j < len(a); j++ {
			if math.Abs(a[i]/na*b[j]/nb-a[j]/na*b[i]/nb) > tol {
				return false
			}
		}
	}
	return true
}
