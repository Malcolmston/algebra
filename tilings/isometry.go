package tilings

import (
	"fmt"
	"math"
)

// Isometry is a planar isometry, the affine map x -> M x + t whose linear part
// M = [[M00, M01], [M10, M11]] is orthogonal (determinant +/-1). Orientation is
// preserved when the determinant is +1 (identity, translation, rotation) and
// reversed when it is -1 (reflection, glide reflection).
type Isometry struct {
	M00, M01, M10, M11 float64
	Tx, Ty             float64
}

// IsometryKind enumerates the five kinds of planar isometry.
type IsometryKind int

const (
	// KindIdentity is the identity map.
	KindIdentity IsometryKind = iota
	// KindTranslation is a nonzero translation.
	KindTranslation
	// KindRotation is a rotation about a fixed point by a nonzero angle.
	KindRotation
	// KindReflection is a reflection across a line.
	KindReflection
	// KindGlideReflection is a reflection across a line composed with a nonzero
	// translation parallel to that line.
	KindGlideReflection
)

// String returns the name of the isometry kind.
func (k IsometryKind) String() string {
	switch k {
	case KindIdentity:
		return "identity"
	case KindTranslation:
		return "translation"
	case KindRotation:
		return "rotation"
	case KindReflection:
		return "reflection"
	case KindGlideReflection:
		return "glide reflection"
	default:
		return "unknown"
	}
}

// Identity returns the identity isometry.
func Identity() Isometry {
	return Isometry{M00: 1, M01: 0, M10: 0, M11: 1}
}

// Translation returns the translation by the vector (dx, dy).
func Translation(dx, dy float64) Isometry {
	return Isometry{M00: 1, M01: 0, M10: 0, M11: 1, Tx: dx, Ty: dy}
}

// TranslationVec returns the translation by the vector v.
func TranslationVec(v Point) Isometry { return Translation(v.X, v.Y) }

// Rotation returns the rotation about the origin by theta radians.
func Rotation(theta float64) Isometry {
	c, s := math.Cos(theta), math.Sin(theta)
	return Isometry{M00: c, M01: -s, M10: s, M11: c}
}

// RotationAbout returns the rotation about the point center by theta radians.
func RotationAbout(center Point, theta float64) Isometry {
	r := Rotation(theta)
	// x -> R(x-c)+c = R x + (c - R c)
	rc := r.applyLinear(center)
	r.Tx = center.X - rc.X
	r.Ty = center.Y - rc.Y
	return r
}

// Reflection returns the reflection across the line through the origin whose
// direction makes angle axis (radians) with the positive x-axis.
func Reflection(axis float64) Isometry {
	c, s := math.Cos(2*axis), math.Sin(2*axis)
	return Isometry{M00: c, M01: s, M10: s, M11: -c}
}

// ReflectionLine returns the reflection across the line through point p whose
// direction makes angle axis (radians) with the positive x-axis.
func ReflectionLine(p Point, axis float64) Isometry {
	r := Reflection(axis)
	rp := r.applyLinear(p)
	r.Tx = p.X - rp.X
	r.Ty = p.Y - rp.Y
	return r
}

// GlideReflection returns the glide reflection whose mirror line passes through
// point p in the direction axis (radians) and whose glide vector has signed
// length glide along that direction.
func GlideReflection(p Point, axis, glide float64) Isometry {
	base := ReflectionLine(p, axis)
	g := Polar(glide, axis)
	base.Tx += g.X
	base.Ty += g.Y
	return base
}

// applyLinear applies only the linear part of the isometry to v.
func (a Isometry) applyLinear(v Point) Point {
	return Point{a.M00*v.X + a.M01*v.Y, a.M10*v.X + a.M11*v.Y}
}

// Apply returns the image of the point p under the isometry.
func (a Isometry) Apply(p Point) Point {
	return Point{a.M00*p.X + a.M01*p.Y + a.Tx, a.M10*p.X + a.M11*p.Y + a.Ty}
}

// ApplyVec returns the image of the vector v under the linear part only,
// ignoring the translation.
func (a Isometry) ApplyVec(v Point) Point { return a.applyLinear(v) }

// Compose returns the isometry a∘inner, i.e. the map applying inner first and
// then a.
func (a Isometry) Compose(inner Isometry) Isometry {
	return Isometry{
		M00: a.M00*inner.M00 + a.M01*inner.M10,
		M01: a.M00*inner.M01 + a.M01*inner.M11,
		M10: a.M10*inner.M00 + a.M11*inner.M10,
		M11: a.M10*inner.M01 + a.M11*inner.M11,
		Tx:  a.M00*inner.Tx + a.M01*inner.Ty + a.Tx,
		Ty:  a.M10*inner.Tx + a.M11*inner.Ty + a.Ty,
	}
}

// Then returns the isometry b∘a, i.e. the map applying a first and then b.
func (a Isometry) Then(b Isometry) Isometry { return b.Compose(a) }

// Determinant returns the determinant of the linear part, either +1 or -1 for a
// true isometry.
func (a Isometry) Determinant() float64 { return a.M00*a.M11 - a.M01*a.M10 }

// IsDirect reports whether the isometry preserves orientation (determinant +1).
func (a Isometry) IsDirect() bool { return a.Determinant() > 0 }

// Inverse returns the inverse isometry.
func (a Isometry) Inverse() Isometry {
	det := a.Determinant()
	// Inverse of orthogonal linear part.
	i00 := a.M11 / det
	i01 := -a.M01 / det
	i10 := -a.M10 / det
	i11 := a.M00 / det
	return Isometry{
		M00: i00, M01: i01, M10: i10, M11: i11,
		Tx: -(i00*a.Tx + i01*a.Ty),
		Ty: -(i10*a.Tx + i11*a.Ty),
	}
}

// Translation returns the translation part (Tx, Ty) of the isometry as a point.
func (a Isometry) Translation() Point { return Point{a.Tx, a.Ty} }

// LinearAngle returns atan2(M10, M00), the rotation angle of the linear part
// for a direct isometry and, for an opposite isometry, twice the mirror-axis
// angle.
func (a Isometry) LinearAngle() float64 { return math.Atan2(a.M10, a.M00) }

// RotationAngle returns the rotation angle in (-pi, pi] for a direct isometry.
// For an opposite isometry it returns 0.
func (a Isometry) RotationAngle() float64 {
	if !a.IsDirect() {
		return 0
	}
	return NormalizeAngleSigned(math.Atan2(a.M10, a.M00))
}

// ApproxEqual reports whether a and b agree to within eps in every entry.
func (a Isometry) ApproxEqual(b Isometry, eps float64) bool {
	return approxEqualScalar(a.M00, b.M00, eps) &&
		approxEqualScalar(a.M01, b.M01, eps) &&
		approxEqualScalar(a.M10, b.M10, eps) &&
		approxEqualScalar(a.M11, b.M11, eps) &&
		approxEqualScalar(a.Tx, b.Tx, eps) &&
		approxEqualScalar(a.Ty, b.Ty, eps)
}

// IsIdentity reports whether the isometry is the identity to within eps.
func (a Isometry) IsIdentity(eps float64) bool {
	return a.ApproxEqual(Identity(), eps)
}

// Order returns the smallest positive integer k such that a^k is the identity
// to within eps, searching up to max. It returns 0 if no such k exists in range
// (for example a nontrivial translation or glide reflection).
func (a Isometry) Order(max int, eps float64) int {
	cur := a
	for k := 1; k <= max; k++ {
		if cur.IsIdentity(eps) {
			return k
		}
		cur = a.Compose(cur)
	}
	return 0
}

// FixedPoint returns a point fixed by the isometry and reports whether one
// exists. Identities, rotations and reflections have fixed points; nonzero
// translations and glide reflections do not.
func (a Isometry) FixedPoint() (Point, bool) {
	// Solve (I-M) x = t.
	d00 := 1 - a.M00
	d01 := -a.M01
	d10 := -a.M10
	d11 := 1 - a.M11
	det := d00*d11 - d01*d10
	if math.Abs(det) < Epsilon {
		return Point{}, false
	}
	x := (d11*a.Tx - d01*a.Ty) / det
	y := (-d10*a.Tx + d00*a.Ty) / det
	return Point{x, y}, true
}

// Classification bundles the invariant data recovered from an isometry.
type Classification struct {
	// Kind is the kind of isometry.
	Kind IsometryKind
	// Angle is the rotation angle in (-pi, pi] for a rotation, and 0 otherwise.
	Angle float64
	// Center is the fixed point of a rotation.
	Center Point
	// AxisPoint is a point on the mirror line of a reflection or glide.
	AxisPoint Point
	// AxisAngle is the direction of the mirror line (radians) for a reflection
	// or glide.
	AxisAngle float64
	// Glide is the signed glide length along the axis for a glide reflection.
	Glide float64
	// TransVec is the translation vector for a translation.
	TransVec Point
}

// Classify decomposes the isometry into its kind and invariant data using the
// tolerance eps.
func (a Isometry) Classify(eps float64) Classification {
	if a.IsDirect() {
		theta := a.RotationAngle()
		if approxEqualScalar(theta, 0, eps) {
			if math.Hypot(a.Tx, a.Ty) <= eps {
				return Classification{Kind: KindIdentity}
			}
			return Classification{Kind: KindTranslation, TransVec: Point{a.Tx, a.Ty}}
		}
		c, _ := a.FixedPoint()
		return Classification{Kind: KindRotation, Angle: theta, Center: c}
	}
	// Opposite isometry: mirror axis angle is half the linear angle.
	axis := math.Atan2(a.M01, a.M00) / 2
	u := Point{math.Cos(axis), math.Sin(axis)}
	t := Point{a.Tx, a.Ty}
	glide := t.Dot(u)
	perp := t.Sub(u.Scale(glide))
	axisPoint := perp.Scale(0.5)
	if math.Abs(glide) <= eps {
		return Classification{Kind: KindReflection, AxisPoint: axisPoint, AxisAngle: NormalizeAngleSigned(axis)}
	}
	return Classification{Kind: KindGlideReflection, AxisPoint: axisPoint, AxisAngle: NormalizeAngleSigned(axis), Glide: glide}
}

// Kind returns just the kind of the isometry using tolerance eps.
func (a Isometry) Kind(eps float64) IsometryKind { return a.Classify(eps).Kind }

// String returns a short human-readable description of the isometry.
func (a Isometry) String() string {
	cl := a.Classify(Epsilon)
	switch cl.Kind {
	case KindIdentity:
		return "identity"
	case KindTranslation:
		return fmt.Sprintf("translation(%.4g,%.4g)", cl.TransVec.X, cl.TransVec.Y)
	case KindRotation:
		return fmt.Sprintf("rotation(%.4g deg about %.4g,%.4g)", Rad2Deg(cl.Angle), cl.Center.X, cl.Center.Y)
	case KindReflection:
		return fmt.Sprintf("reflection(axis %.4g deg through %.4g,%.4g)", Rad2Deg(cl.AxisAngle), cl.AxisPoint.X, cl.AxisPoint.Y)
	default:
		return fmt.Sprintf("glide(axis %.4g deg, glide %.4g)", Rad2Deg(cl.AxisAngle), cl.Glide)
	}
}

// ComposeAll returns the composition of the given isometries applied left to
// right: the first element is applied first. The empty input yields the
// identity.
func ComposeAll(seq ...Isometry) Isometry {
	result := Identity()
	for _, a := range seq {
		result = a.Compose(result)
	}
	return result
}
