package projectivegeom

import (
	"errors"
	"math"
)

// ErrSingular is returned when a matrix that must be invertible is singular.
var ErrSingular = errors.New("projectivegeom: singular matrix")

// ErrConfiguration is returned by fitting routines whose input points are in a
// forbidden special position (for example three of four collinear).
var ErrConfiguration = errors.New("projectivegeom: invalid point configuration")

// Homography is a projective transformation (collineation) of RP^2, represented
// by an invertible 3x3 matrix acting on homogeneous point coordinates. Two
// matrices differing by a non-zero scalar denote the same homography.
type Homography struct {
	M Mat3
}

// IdentityHomography returns the identity transformation.
func IdentityHomography() Homography { return Homography{Identity3()} }

// NewHomography wraps a 3x3 matrix as a homography without checking
// invertibility.
func NewHomography(m Mat3) Homography { return Homography{m} }

// Apply returns the image of the point p under the homography.
func (h Homography) Apply(p Point) Point { return Point{h.M.MulVec(p.V)} }

// ApplyLine returns the image of the line l under the homography. Lines
// transform by the inverse transpose so that incidence is preserved. It returns
// ErrSingular when the homography matrix is not invertible.
func (h Homography) ApplyLine(l Line) (Line, error) {
	inv, ok := h.M.Inverse()
	if !ok {
		return Line{}, ErrSingular
	}
	return Line{inv.Transpose().MulVec(l.V)}, nil
}

// Compose returns the homography g followed by h, whose matrix is h.M * g.M so
// that Compose(h, g).Apply(p) equals h.Apply(g.Apply(p)).
func Compose(h, g Homography) Homography { return Homography{h.M.Mul(g.M)} }

// Then returns the homography that applies h and then g, i.e. Compose(g, h).
func (h Homography) Then(g Homography) Homography { return Compose(g, h) }

// Inverse returns the inverse homography, or ErrSingular when the matrix is not
// invertible.
func (h Homography) Inverse() (Homography, error) {
	inv, ok := h.M.Inverse()
	if !ok {
		return Homography{}, ErrSingular
	}
	return Homography{inv}, nil
}

// Determinant returns the determinant of the homography matrix. Its sign
// indicates whether the map is orientation preserving.
func (h Homography) Determinant() float64 { return h.M.Det() }

// IsInvertible reports whether the homography matrix is numerically invertible.
func (h Homography) IsInvertible() bool { return math.Abs(h.M.Det()) > Eps*Eps }

// ApproxEqual reports whether h and g denote the same homography, comparing
// their matrices after scaling each by its largest-magnitude entry.
func (h Homography) ApproxEqual(g Homography, tol float64) bool {
	sh, sg := h.M.MaxAbs(), g.M.MaxAbs()
	if sh < Eps || sg < Eps {
		return false
	}
	a := h.M.Scale(1 / sh)
	b := g.M.Scale(1 / sg)
	return a.ApproxEqual(b, tol) || a.ApproxEqual(b.Scale(-1), tol)
}

// Translation returns the affine translation by (dx, dy).
func Translation(dx, dy float64) Homography {
	return Homography{NewMat3(1, 0, dx, 0, 1, dy, 0, 0, 1)}
}

// ScalingXY returns the affine scaling by sx along x and sy along y about the
// origin.
func ScalingXY(sx, sy float64) Homography {
	return Homography{NewMat3(sx, 0, 0, 0, sy, 0, 0, 0, 1)}
}

// Scaling returns the uniform affine scaling by s about the origin.
func Scaling(s float64) Homography { return ScalingXY(s, s) }

// Rotation returns the Euclidean rotation by theta radians about the origin.
func Rotation(theta float64) Homography {
	c, s := math.Cos(theta), math.Sin(theta)
	return Homography{NewMat3(c, -s, 0, s, c, 0, 0, 0, 1)}
}

// RotationAbout returns the Euclidean rotation by theta radians about the
// finite point (cx, cy).
func RotationAbout(theta, cx, cy float64) Homography {
	return Compose(Translation(cx, cy), Compose(Rotation(theta), Translation(-cx, -cy)))
}

// Similarity returns the direct similarity that scales by s, rotates by theta
// and then translates by (dx, dy).
func Similarity(s, theta, dx, dy float64) Homography {
	return Compose(Translation(dx, dy), Compose(Rotation(theta), Scaling(s)))
}

// Affine returns the affine map [[a b tx],[c d ty],[0 0 1]].
func Affine(a, b, c, d, tx, ty float64) Homography {
	return Homography{NewMat3(a, b, tx, c, d, ty, 0, 0, 1)}
}

// ShearX returns the affine shear x' = x + k*y.
func ShearX(k float64) Homography { return Homography{NewMat3(1, k, 0, 0, 1, 0, 0, 0, 1)} }

// ShearY returns the affine shear y' = y + k*x.
func ShearY(k float64) Homography { return Homography{NewMat3(1, 0, 0, k, 1, 0, 0, 0, 1)} }

// ReflectionX returns the Euclidean reflection across the x-axis.
func ReflectionX() Homography { return ScalingXY(1, -1) }

// ReflectionY returns the Euclidean reflection across the y-axis.
func ReflectionY() Homography { return ScalingXY(-1, 1) }

// IsAffine reports whether the homography fixes the line at infinity, i.e. its
// bottom row is proportional to (0, 0, 1) within tol.
func (h Homography) IsAffine(tol float64) bool {
	s := h.M.MaxAbs()
	if s < Eps {
		return false
	}
	return math.Abs(h.M[2][0]) <= tol*s && math.Abs(h.M[2][1]) <= tol*s
}

// IsSimilarity reports whether the homography is affine and its linear 2x2 part
// is a scaled orthogonal matrix within tol (a Euclidean similarity).
func (h Homography) IsSimilarity(tol float64) bool {
	if !h.IsAffine(tol) {
		return false
	}
	a, b := h.M[0][0], h.M[0][1]
	c, d := h.M[1][0], h.M[1][1]
	// Columns of the linear part must be orthogonal and of equal length.
	col1 := a*a + c*c
	col2 := b*b + d*d
	dot := a*b + c*d
	scale := math.Max(col1, col2)
	if scale < Eps {
		return false
	}
	return math.Abs(col1-col2) <= tol*scale && math.Abs(dot) <= tol*scale
}

// basisToPoints returns the matrix taking the projective basis (e1, e2, e3, u)
// with u = [1 1 1] to the four given points, which must be in general position
// (no three collinear). It returns ErrConfiguration otherwise.
func basisToPoints(p1, p2, p3, p4 Point) (Mat3, error) {
	base := Mat3{
		{p1.V.X, p2.V.X, p3.V.X},
		{p1.V.Y, p2.V.Y, p3.V.Y},
		{p1.V.Z, p2.V.Z, p3.V.Z},
	}
	inv, ok := base.Inverse()
	if !ok {
		return Mat3{}, ErrConfiguration
	}
	lam := inv.MulVec(p4.V)
	if math.Abs(lam.X) < Eps || math.Abs(lam.Y) < Eps || math.Abs(lam.Z) < Eps {
		return Mat3{}, ErrConfiguration
	}
	return Mat3{
		{base[0][0] * lam.X, base[0][1] * lam.Y, base[0][2] * lam.Z},
		{base[1][0] * lam.X, base[1][1] * lam.Y, base[1][2] * lam.Z},
		{base[2][0] * lam.X, base[2][1] * lam.Y, base[2][2] * lam.Z},
	}, nil
}

// HomographyFromPoints returns the unique homography mapping the four source
// points to the four destination points, given no three points in either set
// are collinear. It returns ErrConfiguration when either set is degenerate.
func HomographyFromPoints(src, dst [4]Point) (Homography, error) {
	a, err := basisToPoints(src[0], src[1], src[2], src[3])
	if err != nil {
		return Homography{}, err
	}
	b, err := basisToPoints(dst[0], dst[1], dst[2], dst[3])
	if err != nil {
		return Homography{}, err
	}
	aInv, ok := a.Inverse()
	if !ok {
		return Homography{}, ErrConfiguration
	}
	return Homography{b.Mul(aInv)}, nil
}

// FixedPointsReal returns the real fixed points of the homography, that is the
// points p with h.Apply(p) proportional to p. These are the eigenvectors of the
// matrix with real eigenvalues. Up to three are returned; complex eigenvalue
// pairs are omitted.
func (h Homography) FixedPointsReal(tol float64) []Point {
	vals := realEigenvalues3(h.M)
	var out []Point
	for _, lam := range vals {
		v, ok := eigenvector3(h.M, lam, tol)
		if !ok {
			continue
		}
		dup := false
		for _, q := range out {
			if (Point{v}).Equal(Point{q.V}, tol*10) {
				dup = true
				break
			}
		}
		if !dup {
			out = append(out, Point{v})
		}
	}
	return out
}
