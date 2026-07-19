package projectivegeom

import "math"

// SPoint is a point of real projective space RP^3 in homogeneous coordinates
// [X Y Z W], defined up to a non-zero scalar. The finite point (x, y, z)
// corresponds to [x y z 1].
type SPoint struct {
	V Vec4
}

// SPlane is a plane of RP^3 in homogeneous (dual) coordinates [A B C D] denoting
// the locus A*X + B*Y + C*Z + D*W = 0.
type SPlane struct {
	V Vec4
}

// NewSPoint returns the homogeneous point [x y z w].
func NewSPoint(x, y, z, w float64) SPoint { return SPoint{Vec4{x, y, z, w}} }

// SPointFromAffine returns the finite point [x y z 1].
func SPointFromAffine(x, y, z float64) SPoint { return SPoint{Vec4{x, y, z, 1}} }

// NewSPlane returns the homogeneous plane [a b c d] for a*x + b*y + c*z + d = 0.
func NewSPlane(a, b, c, d float64) SPlane { return SPlane{Vec4{a, b, c, d}} }

// PlaneAtInfinity returns the plane at infinity [0 0 0 1] of RP^3.
func PlaneAtInfinity() SPlane { return SPlane{Vec4{0, 0, 0, 1}} }

// IsValid reports whether the point is a non-zero homogeneous 4-tuple.
func (p SPoint) IsValid() bool { return !p.V.IsZero(Eps) }

// IsValid reports whether the plane is a non-zero homogeneous 4-tuple.
func (s SPlane) IsValid() bool { return !s.V.IsZero(Eps) }

// IsAtInfinity reports whether the point lies on the plane at infinity within
// tol.
func (p SPoint) IsAtInfinity(tol float64) bool {
	s := p.V.MaxAbs()
	if s < Eps {
		return false
	}
	return math.Abs(p.V.W) <= tol*s
}

// Affine returns the Euclidean coordinates (x, y, z) = (X/W, Y/W, Z/W), or
// ErrAtInfinity when the point lies on the plane at infinity.
func (p SPoint) Affine() (x, y, z float64, err error) {
	if p.IsAtInfinity(Eps) {
		return 0, 0, 0, ErrAtInfinity
	}
	return p.V.X / p.V.W, p.V.Y / p.V.W, p.V.Z / p.V.W, nil
}

// Equal reports whether the two points denote the same projective point.
func (p SPoint) Equal(q SPoint, tol float64) bool { return vec4Parallel(p.V, q.V, tol) }

// Equal reports whether the two planes denote the same projective plane.
func (s SPlane) Equal(t SPlane, tol float64) bool { return vec4Parallel(s.V, t.V, tol) }

// IncidencePointPlane returns the value p.V·s.V, which is zero exactly when the
// point lies on the plane.
func IncidencePointPlane(p SPoint, s SPlane) float64 { return p.V.Dot(s.V) }

// OnPlane reports whether the point lies on the plane within tol using
// scale-invariant normalization.
func OnPlane(p SPoint, s SPlane, tol float64) bool {
	sp, ss := p.V.MaxAbs(), s.V.MaxAbs()
	if sp < Eps || ss < Eps {
		return false
	}
	return math.Abs(p.V.Dot(s.V)/(sp*ss)) <= tol
}

// PlaneThrough3Points returns the plane through three points of RP^3. It returns
// ErrDegenerate when the points are collinear (the plane is undetermined).
func PlaneThrough3Points(p, q, r SPoint) (SPlane, error) {
	v := Cross4(p.V, q.V, r.V)
	if v.IsZero(Eps) {
		return SPlane{}, ErrDegenerate
	}
	return SPlane{v}, nil
}

// PointOf3Planes returns the common point of three planes of RP^3. It returns
// ErrDegenerate when the planes share a common line or are otherwise
// dependent.
func PointOf3Planes(a, b, c SPlane) (SPoint, error) {
	v := Cross4(a.V, b.V, c.V)
	if v.IsZero(Eps) {
		return SPoint{}, ErrDegenerate
	}
	return SPoint{v}, nil
}

// Coplanar reports whether four points of RP^3 are coplanar within tol, i.e. the
// determinant of their coordinates vanishes after scale normalization.
func Coplanar(p, q, r, s SPoint, tol float64) bool {
	scaleAndRow := func(v Vec4) []float64 {
		n := v.MaxAbs()
		if n < Eps {
			n = 1
		}
		return []float64{v.X / n, v.Y / n, v.Z / n, v.W / n}
	}
	d, _ := detN([][]float64{
		scaleAndRow(p.V), scaleAndRow(q.V), scaleAndRow(r.V), scaleAndRow(s.V),
	})
	return math.Abs(d) <= tol
}

// Concurrent4Planes reports whether four planes of RP^3 pass through a common
// point within tol (the dual of Coplanar).
func Concurrent4Planes(a, b, c, d SPlane, tol float64) bool {
	return Coplanar(SPoint{a.V}, SPoint{b.V}, SPoint{c.V}, SPoint{d.V}, tol)
}

// Collineation3 is a projective transformation of RP^3, an invertible 4x4
// matrix acting on homogeneous point coordinates.
type Collineation3 struct {
	M Mat4
}

// IdentityCollineation3 returns the identity transformation of RP^3.
func IdentityCollineation3() Collineation3 { return Collineation3{Identity4()} }

// Apply returns the image of the point p under the collineation.
func (c Collineation3) Apply(p SPoint) SPoint { return SPoint{c.M.MulVec(p.V)} }

// ApplyPlane returns the image of the plane s, transformed by the inverse
// transpose so that incidence is preserved. It returns ErrSingular when the
// matrix is not invertible.
func (c Collineation3) ApplyPlane(s SPlane) (SPlane, error) {
	inv, ok := c.M.Inverse()
	if !ok {
		return SPlane{}, ErrSingular
	}
	return SPlane{inv.Transpose().MulVec(s.V)}, nil
}

// Inverse returns the inverse collineation, or ErrSingular when the matrix is
// not invertible.
func (c Collineation3) Inverse() (Collineation3, error) {
	inv, ok := c.M.Inverse()
	if !ok {
		return Collineation3{}, ErrSingular
	}
	return Collineation3{inv}, nil
}

// Compose3 returns the collineation g followed by h.
func Compose3(h, g Collineation3) Collineation3 { return Collineation3{h.M.Mul(g.M)} }

// Translation3 returns the affine translation of RP^3 by (dx, dy, dz).
func Translation3(dx, dy, dz float64) Collineation3 {
	m := Identity4()
	m[0][3] = dx
	m[1][3] = dy
	m[2][3] = dz
	return Collineation3{m}
}

// vec4Parallel reports whether two 4-vectors are scalar multiples within tol,
// after normalizing each by its largest-magnitude component.
func vec4Parallel(a, b Vec4, tol float64) bool {
	na, nb := a.MaxAbs(), b.MaxAbs()
	if na < Eps || nb < Eps {
		return na < Eps && nb < Eps
	}
	u := a.Scale(1 / na)
	v := b.Scale(1 / nb)
	// Cross-difference test: u_i v_j - u_j v_i == 0 for all i<j.
	comps := [4]float64{u.X, u.Y, u.Z, u.W}
	compv := [4]float64{v.X, v.Y, v.Z, v.W}
	for i := 0; i < 4; i++ {
		for j := i + 1; j < 4; j++ {
			if math.Abs(comps[i]*compv[j]-comps[j]*compv[i]) > tol {
				return false
			}
		}
	}
	return true
}
