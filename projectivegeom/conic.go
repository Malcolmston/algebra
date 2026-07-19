package projectivegeom

import (
	"errors"
	"math"
)

// ErrNotOnConic is returned when a tangent is requested at a point that does
// not lie on the conic.
var ErrNotOnConic = errors.New("projectivegeom: point is not on the conic")

// Conic is a plane conic represented by a symmetric 3x3 matrix M, denoting the
// locus of homogeneous points x with x^T M x = 0. The point conic is
//
//	a*X^2 + b*X*Y + c*Y^2 + d*X*Z + e*Y*Z + f*Z^2 = 0
//
// with M = [[a, b/2, d/2], [b/2, c, e/2], [d/2, e/2, f]].
type Conic struct {
	M Mat3
}

// NewConic builds the conic with the given point-equation coefficients a, b, c,
// d, e, f.
func NewConic(a, b, c, d, e, f float64) Conic {
	return Conic{NewMat3(a, b/2, d/2, b/2, c, e/2, d/2, e/2, f)}
}

// ConicFromMatrix wraps a matrix as a conic, replacing it with its symmetric
// part so that only the quadratic form matters.
func ConicFromMatrix(m Mat3) Conic { return Conic{m.Symmetrize()} }

// Coeffs returns the point-equation coefficients (a, b, c, d, e, f).
func (q Conic) Coeffs() (a, b, c, d, e, f float64) {
	return q.M[0][0], 2 * q.M[0][1], q.M[1][1], 2 * q.M[0][2], 2 * q.M[1][2], q.M[2][2]
}

// Evaluate returns the value x^T M x of the conic's quadratic form at p, which
// is zero exactly when p lies on the conic.
func (q Conic) Evaluate(p Point) float64 { return q.M.Quad(p.V, p.V) }

// OnConic reports whether p lies on the conic within tol, using scale-invariant
// normalization.
func (q Conic) OnConic(p Point, tol float64) bool {
	sp := p.V.MaxAbs()
	sm := q.M.MaxAbs()
	if sp < Eps || sm < Eps {
		return false
	}
	return math.Abs(q.Evaluate(p)/(sp*sp*sm)) <= tol
}

// ConicFromFivePoints returns the unique conic through the five given points,
// computed from the 5x5 minors of their monomial matrix. It returns
// ErrConfiguration when the five points fail to determine a conic (for example
// four of them coincide) yielding a numerically zero matrix.
func ConicFromFivePoints(pts [5]Point) (Conic, error) {
	// Monomials in the order [X^2, XY, Y^2, XZ, YZ, Z^2].
	rows := make([][]float64, 5)
	for i, p := range pts {
		x, y, z := p.V.X, p.V.Y, p.V.Z
		rows[i] = []float64{x * x, x * y, y * y, x * z, y * z, z * z}
	}
	coeff := make([]float64, 6)
	for j := 0; j < 6; j++ {
		minor := make([][]float64, 5)
		for i := 0; i < 5; i++ {
			r := make([]float64, 0, 5)
			for k := 0; k < 6; k++ {
				if k == j {
					continue
				}
				r = append(r, rows[i][k])
			}
			minor[i] = r
		}
		d, _ := detN(minor)
		if j%2 == 1 {
			d = -d
		}
		coeff[j] = d
	}
	q := NewConic(coeff[0], coeff[1], coeff[2], coeff[3], coeff[4], coeff[5])
	if q.M.MaxAbs() < Eps {
		return Conic{}, ErrConfiguration
	}
	return q, nil
}

// Polar returns the polar line of the point p with respect to the conic, namely
// the line M*p. When p lies on the conic the polar is the tangent at p.
func (q Conic) Polar(p Point) Line { return Line{q.M.MulVec(p.V)} }

// Pole returns the pole of the line l with respect to the conic, the point
// whose polar is l. It uses the adjugate of M so that the pole is well defined
// (possibly at infinity) whenever the conic is non-degenerate. It returns
// ErrDegenerate when the conic matrix is singular.
func (q Conic) Pole(l Line) (Point, error) {
	adj := q.M.Adjugate()
	p := adj.MulVec(l.V)
	if p.IsZero(Eps) {
		return Point{}, ErrDegenerate
	}
	return Point{p}, nil
}

// TangentAt returns the tangent line to the conic at a point p that lies on it,
// or ErrNotOnConic when p is not on the conic within a mild tolerance.
func (q Conic) TangentAt(p Point) (Line, error) {
	if !q.OnConic(p, 1e-6) {
		return Line{}, ErrNotOnConic
	}
	pl := q.Polar(p)
	if pl.V.IsZero(Eps) {
		return Line{}, ErrDegenerate
	}
	return pl, nil
}

// twoPointsOnLine returns two distinct points that span the line l.
func twoPointsOnLine(l Line) (Point, Point) {
	cands := []Vec3{
		l.V.Cross(Vec3{1, 0, 0}),
		l.V.Cross(Vec3{0, 1, 0}),
		l.V.Cross(Vec3{0, 0, 1}),
	}
	// Choose the two candidates of largest norm (they are the best conditioned
	// and guaranteed independent when l is non-zero).
	best, second := 0, 1
	if cands[1].Norm() > cands[0].Norm() {
		best, second = 1, 0
	}
	for i := 2; i < 3; i++ {
		if cands[i].Norm() > cands[best].Norm() {
			second = best
			best = i
		} else if cands[i].Norm() > cands[second].Norm() {
			second = i
		}
	}
	return Point{cands[best]}, Point{cands[second]}
}

// IntersectLine returns the intersection points of the line l with the conic:
// zero points when the line misses the conic, one (returned twice) when it is
// tangent, or two when it is secant. Complex intersections are omitted.
func (q Conic) IntersectLine(l Line) []Point {
	P, Q := twoPointsOnLine(l)
	alpha := q.M.Quad(Q.V, Q.V)
	beta := 2 * q.M.Quad(P.V, Q.V)
	gamma := q.M.Quad(P.V, P.V)
	scale := q.M.MaxAbs() * math.Max(P.V.Norm2(), Q.V.Norm2())
	if scale < Eps {
		return nil
	}
	if math.Abs(alpha) < 1e-12*scale {
		// Linear: beta t + gamma = 0.
		if math.Abs(beta) < 1e-12*scale {
			return nil
		}
		t := -gamma / beta
		return []Point{{P.V.Add(Q.V.Scale(t))}}
	}
	disc := beta*beta - 4*alpha*gamma
	if disc < -1e-9*scale*scale {
		return nil
	}
	if disc < 0 {
		disc = 0
	}
	sq := math.Sqrt(disc)
	t1 := (-beta + sq) / (2 * alpha)
	t2 := (-beta - sq) / (2 * alpha)
	p1 := Point{P.V.Add(Q.V.Scale(t1))}
	p2 := Point{P.V.Add(Q.V.Scale(t2))}
	if sq < 1e-9*(1+math.Abs(beta)) {
		return []Point{p1, p1}
	}
	return []Point{p1, p2}
}

// TangentsFrom returns the two tangent lines to the conic through an external
// point p, found as the joins of p with the two contact points where the polar
// of p meets the conic. It returns ErrDegenerate when p lies on the conic (a
// single tangent) or when the construction fails.
func (q Conic) TangentsFrom(p Point) ([2]Line, error) {
	polar := q.Polar(p)
	contacts := q.IntersectLine(polar)
	if len(contacts) < 2 {
		return [2]Line{}, ErrDegenerate
	}
	t1, err := Join(p, contacts[0])
	if err != nil {
		return [2]Line{}, err
	}
	t2, err := Join(p, contacts[1])
	if err != nil {
		return [2]Line{}, err
	}
	return [2]Line{t1, t2}, nil
}

// Center returns the center of the conic, the pole of the line at infinity. For
// a parabola the center lies at infinity. It returns ErrDegenerate when the
// conic matrix is singular.
func (q Conic) Center() (Point, error) { return q.Pole(LineAtInfinity()) }

// Det returns the determinant of the conic matrix. Its vanishing marks a
// degenerate conic.
func (q Conic) Det() float64 { return q.M.Det() }

// Discriminant returns the determinant of the upper-left 2x2 block, a*c -
// (b/2)^2, whose sign classifies a non-degenerate conic as an ellipse
// (positive), parabola (zero) or hyperbola (negative).
func (q Conic) Discriminant() float64 {
	return q.M[0][0]*q.M[1][1] - q.M[0][1]*q.M[0][1]
}

// IsDegenerate reports whether the conic matrix is singular within tol,
// measured relative to the cube of its largest entry.
func (q Conic) IsDegenerate(tol float64) bool {
	s := q.M.MaxAbs()
	if s < Eps {
		return true
	}
	return math.Abs(q.Det()) <= tol*s*s*s
}

// ConicType classifies the affine shape of a conic.
type ConicType int

const (
	// ConicDegenerate marks a singular conic (line pair, double line or point).
	ConicDegenerate ConicType = iota
	// ConicEllipse marks a real or imaginary ellipse.
	ConicEllipse
	// ConicParabola marks a parabola.
	ConicParabola
	// ConicHyperbola marks a hyperbola.
	ConicHyperbola
)

// String returns the conic type name.
func (t ConicType) String() string {
	switch t {
	case ConicEllipse:
		return "ellipse"
	case ConicParabola:
		return "parabola"
	case ConicHyperbola:
		return "hyperbola"
	default:
		return "degenerate"
	}
}

// Classify returns the affine type of the conic. A near-singular matrix is
// reported as ConicDegenerate; otherwise the sign of the discriminant selects
// ellipse, parabola or hyperbola. The tol governs both tests relative to the
// matrix scale.
func (q Conic) Classify(tol float64) ConicType {
	if q.IsDegenerate(tol) {
		return ConicDegenerate
	}
	s := q.M.MaxAbs()
	disc := q.Discriminant()
	if math.Abs(disc) <= tol*s*s {
		return ConicParabola
	}
	if disc > 0 {
		return ConicEllipse
	}
	return ConicHyperbola
}

// IsEllipse reports whether the conic classifies as an ellipse within tol.
func (q Conic) IsEllipse(tol float64) bool { return q.Classify(tol) == ConicEllipse }

// IsParabola reports whether the conic classifies as a parabola within tol.
func (q Conic) IsParabola(tol float64) bool { return q.Classify(tol) == ConicParabola }

// IsHyperbola reports whether the conic classifies as a hyperbola within tol.
func (q Conic) IsHyperbola(tol float64) bool { return q.Classify(tol) == ConicHyperbola }

// Rank returns the numerical rank (0..3) of the conic matrix, distinguishing a
// proper conic (rank 3) from a line pair (rank 2) or double line (rank 1). The
// tol is relative to the matrix scale.
func (q Conic) Rank(tol float64) int {
	return matrixRank3(q.M, tol)
}

// Transform returns the image of the conic under the homography h. A point x
// lies on the image conic when h^{-1} x lies on the original, giving image
// matrix h^{-T} M h^{-1}. It returns ErrSingular when h is not invertible.
func (q Conic) Transform(h Homography) (Conic, error) {
	inv, ok := h.M.Inverse()
	if !ok {
		return Conic{}, ErrSingular
	}
	m := inv.Transpose().Mul(q.M).Mul(inv)
	return Conic{m.Symmetrize()}, nil
}

// Dual returns the dual (line) conic, whose points are the tangent lines of the
// original. Its matrix is the adjugate of M. It returns ErrDegenerate when the
// adjugate vanishes.
func (q Conic) Dual() (Conic, error) {
	adj := q.M.Adjugate().Symmetrize()
	if adj.MaxAbs() < Eps {
		return Conic{}, ErrDegenerate
	}
	return Conic{adj}, nil
}

// UnitCircle returns the conic X^2 + Y^2 - Z^2 = 0, the unit circle.
func UnitCircle() Conic { return NewConic(1, 0, 1, 0, 0, -1) }

// CircleConic returns the conic for the circle of radius r centered at (cx, cy).
func CircleConic(cx, cy, r float64) Conic {
	return NewConic(1, 0, 1, -2*cx, -2*cy, cx*cx+cy*cy-r*r)
}

// matrixRank3 returns the numerical rank of a 3x3 matrix via row reduction with
// a scale-relative pivot threshold.
func matrixRank3(m Mat3, tol float64) int {
	s := m.MaxAbs()
	if s < Eps {
		return 0
	}
	a := [3][3]float64{m[0], m[1], m[2]}
	rank := 0
	thr := tol * s
	for col := 0; col < 3 && rank < 3; col++ {
		piv := -1
		best := thr
		for r := rank; r < 3; r++ {
			if v := math.Abs(a[r][col]); v > best {
				best, piv = v, r
			}
		}
		if piv < 0 {
			continue
		}
		a[rank], a[piv] = a[piv], a[rank]
		inv := 1 / a[rank][col]
		for r := 0; r < 3; r++ {
			if r == rank {
				continue
			}
			f := a[r][col] * inv
			for c := 0; c < 3; c++ {
				a[r][c] -= f * a[rank][c]
			}
		}
		rank++
	}
	return rank
}
