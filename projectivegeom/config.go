package projectivegeom

// PappusConfig holds two triples of collinear points, one triple on each of two
// lines, as required by the theorem of Pappus.
type PappusConfig struct {
	A, B, C    Point // three points on the first line
	A2, B2, C2 Point // three points on the second line
}

// CrossPoints returns the three Pappus cross points
//
//	X = (A B2) ∩ (A2 B),  Y = (B C2) ∩ (B2 C),  Z = (C A2) ∩ (C2 A),
//
// which the theorem of Pappus asserts are collinear whenever A, B, C lie on one
// line and A2, B2, C2 on another. It returns an error when any join or meet is
// degenerate.
func (p PappusConfig) CrossPoints() (x, y, z Point, err error) {
	meetJoins := func(u1, v1, u2, v2 Point) (Point, error) {
		l1, e := Join(u1, v1)
		if e != nil {
			return Point{}, e
		}
		l2, e := Join(u2, v2)
		if e != nil {
			return Point{}, e
		}
		return Meet(l1, l2)
	}
	if x, err = meetJoins(p.A, p.B2, p.A2, p.B); err != nil {
		return
	}
	if y, err = meetJoins(p.B, p.C2, p.B2, p.C); err != nil {
		return
	}
	z, err = meetJoins(p.C, p.A2, p.C2, p.A)
	return
}

// PappusLine returns the line containing the three Pappus cross points, the
// Pappus line of the configuration. It returns an error when the cross points
// are not defined or happen to coincide.
func (p PappusConfig) PappusLine() (Line, error) {
	x, y, _, err := p.CrossPoints()
	if err != nil {
		return Line{}, err
	}
	return Join(x, y)
}

// Holds reports whether the theorem of Pappus holds numerically for this
// configuration, i.e. the three cross points are collinear within tol. It
// returns false on any degeneracy.
func (p PappusConfig) Holds(tol float64) bool {
	x, y, z, err := p.CrossPoints()
	if err != nil {
		return false
	}
	return Collinear(x, y, z, tol)
}

// TrianglePair holds two triangles for the theorem of Desargues.
type TrianglePair struct {
	T1, T2 Triangle
}

// PerspectiveCenter returns the center of perspectivity of the two triangles,
// the common point of the lines A A', B B', C C'. It returns an error when the
// triangles are not perspective from a point or a join is degenerate.
func (tp TrianglePair) PerspectiveCenter() (Point, error) {
	l1, err := Join(tp.T1.A, tp.T2.A)
	if err != nil {
		return Point{}, err
	}
	l2, err := Join(tp.T1.B, tp.T2.B)
	if err != nil {
		return Point{}, err
	}
	return Meet(l1, l2)
}

// IsPerspectiveFromPoint reports whether the three lines joining corresponding
// vertices are concurrent within tol.
func (tp TrianglePair) IsPerspectiveFromPoint(tol float64) bool {
	l1, err := Join(tp.T1.A, tp.T2.A)
	if err != nil {
		return false
	}
	l2, err := Join(tp.T1.B, tp.T2.B)
	if err != nil {
		return false
	}
	l3, err := Join(tp.T1.C, tp.T2.C)
	if err != nil {
		return false
	}
	return Concurrent(l1, l2, l3, tol)
}

// AxisPoints returns the three intersections of corresponding sides,
//
//	P = (AB) ∩ (A'B'),  Q = (BC) ∩ (B'C'),  R = (CA) ∩ (C'A'),
//
// which Desargues' theorem asserts are collinear exactly when the triangles are
// perspective from a point. It returns an error on any degeneracy.
func (tp TrianglePair) AxisPoints() (p, q, r Point, err error) {
	sideMeet := func(u1, v1, u2, v2 Point) (Point, error) {
		l1, e := Join(u1, v1)
		if e != nil {
			return Point{}, e
		}
		l2, e := Join(u2, v2)
		if e != nil {
			return Point{}, e
		}
		return Meet(l1, l2)
	}
	if p, err = sideMeet(tp.T1.A, tp.T1.B, tp.T2.A, tp.T2.B); err != nil {
		return
	}
	if q, err = sideMeet(tp.T1.B, tp.T1.C, tp.T2.B, tp.T2.C); err != nil {
		return
	}
	r, err = sideMeet(tp.T1.C, tp.T1.A, tp.T2.C, tp.T2.A)
	return
}

// DesarguesAxis returns the axis of perspectivity, the line through the three
// AxisPoints. It returns an error when those points are undefined or coincide.
func (tp TrianglePair) DesarguesAxis() (Line, error) {
	p, q, _, err := tp.AxisPoints()
	if err != nil {
		return Line{}, err
	}
	return Join(p, q)
}

// IsPerspectiveFromLine reports whether the three AxisPoints are collinear
// within tol.
func (tp TrianglePair) IsPerspectiveFromLine(tol float64) bool {
	p, q, r, err := tp.AxisPoints()
	if err != nil {
		return false
	}
	return Collinear(p, q, r, tol)
}
