package projectivegeom

// DualOfPoint returns the line whose homogeneous coordinates equal those of the
// point p. Under the duality of RP^2 this exchanges the roles of points and
// lines; combined constructions transform into their duals.
func DualOfPoint(p Point) Line { return Line{p.V} }

// DualOfLine returns the point whose homogeneous coordinates equal those of the
// line l, the dual of DualOfPoint.
func DualOfLine(l Line) Point { return Point{l.V} }

// DualPointOfPlane returns the point of the dual space RP^3 whose coordinates
// equal those of the plane s.
func DualPointOfPlane(s SPlane) SPoint { return SPoint{s.V} }

// DualPlaneOfPoint returns the plane whose coordinates equal those of the point
// p, the dual of DualPointOfPlane.
func DualPlaneOfPoint(p SPoint) SPlane { return SPlane{p.V} }

// PolarPoint is a synonym for Conic.Polar reading the point as the pole; it
// returns the polar line of p with respect to q.
func PolarPoint(q Conic, p Point) Line { return q.Polar(p) }

// ConjugatePoints reports whether p1 and p2 are conjugate with respect to the
// conic, i.e. each lies on the polar of the other, within tol. The relation is
// symmetric because the conic matrix is symmetric.
func ConjugatePoints(q Conic, p1, p2 Point, tol float64) bool {
	val := q.M.Quad(p1.V, p2.V)
	s := q.M.MaxAbs() * p1.V.MaxAbs() * p2.V.MaxAbs()
	if s < Eps {
		return false
	}
	return absF(val)/s <= tol
}

// ConjugateLines reports whether two lines are conjugate with respect to the
// non-degenerate conic, i.e. the pole of each lies on the other, within tol.
func ConjugateLines(q Conic, l1, l2 Line, tol float64) bool {
	adj := q.M.Adjugate()
	val := adj.Quad(l1.V, l2.V)
	s := adj.MaxAbs() * l1.V.MaxAbs() * l2.V.MaxAbs()
	if s < Eps {
		return false
	}
	return absF(val)/s <= tol
}

func absF(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}
