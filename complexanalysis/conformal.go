package complexanalysis

import "math/cmplx"

// Mobius represents a Mobius (linear fractional) transformation
// z -> (A*z + B) / (C*z + D). Such maps are the conformal automorphisms of the
// extended complex plane (the Riemann sphere).
type Mobius struct {
	A, B, C, D complex128
}

// NewMobius constructs the Mobius transformation with the given coefficients.
func NewMobius(a, b, c, d complex128) Mobius {
	return Mobius{A: a, B: b, C: c, D: d}
}

// IdentityMobius returns the identity transformation z -> z.
func IdentityMobius() Mobius {
	return Mobius{A: 1, B: 0, C: 0, D: 1}
}

// Apply evaluates the transformation at z. When the denominator vanishes it
// returns a complex infinity, matching the point at infinity on the sphere.
func (m Mobius) Apply(z complex128) complex128 {
	denom := m.C*z + m.D
	if denom == 0 {
		return cmplx.Inf()
	}
	return (m.A*z + m.B) / denom
}

// Determinant returns A*D - B*C. A Mobius map is invertible exactly when this is
// non-zero.
func (m Mobius) Determinant() complex128 {
	return m.A*m.D - m.B*m.C
}

// Compose returns the transformation m(other(z)), the composition of m after
// other. Its coefficient matrix is the product of the two coefficient matrices.
func (m Mobius) Compose(other Mobius) Mobius {
	return Mobius{
		A: m.A*other.A + m.B*other.C,
		B: m.A*other.B + m.B*other.D,
		C: m.C*other.A + m.D*other.C,
		D: m.C*other.B + m.D*other.D,
	}
}

// Inverse returns the inverse transformation z -> (D*z - B) / (-C*z + A).
func (m Mobius) Inverse() Mobius {
	return Mobius{A: m.D, B: -m.B, C: -m.C, D: m.A}
}

// Normalize returns an equivalent transformation whose coefficient matrix has
// determinant 1, dividing all coefficients by a square root of the current
// determinant. A zero-determinant map is returned unchanged.
func (m Mobius) Normalize() Mobius {
	det := m.Determinant()
	if det == 0 {
		return m
	}
	s := cmplx.Sqrt(det)
	return Mobius{A: m.A / s, B: m.B / s, C: m.C / s, D: m.D / s}
}

// FixedPoints returns the fixed points of the transformation, the solutions of
// C*z^2 + (D-A)*z - B = 0. For a non-identity map this is one or two points;
// when C == 0 the point at infinity is a fixed point and only the finite ones
// are returned. The identity map (and any map fixing every point) returns nil.
func (m Mobius) FixedPoints() []complex128 {
	if m.C == 0 {
		if m.A == m.D {
			// z -> z + B/D: only infinity is fixed unless B == 0 (identity).
			return nil
		}
		// (A z + B)/D = z  =>  z = B/(D-A).
		return []complex128{m.B / (m.D - m.A)}
	}
	a := m.C
	b := m.D - m.A
	c := -m.B
	disc := cmplx.Sqrt(b*b - 4*a*c)
	z1 := (-b + disc) / (2 * a)
	z2 := (-b - disc) / (2 * a)
	if z1 == z2 {
		return []complex128{z1}
	}
	return []complex128{z1, z2}
}

// TranslationMap returns the Mobius map z -> z + b.
func TranslationMap(b complex128) Mobius {
	return Mobius{A: 1, B: b, C: 0, D: 1}
}

// ScalingMap returns the Mobius map z -> a*z. For |a| = 1 this is a pure
// rotation and RotationMap is a convenient alias.
func ScalingMap(a complex128) Mobius {
	return Mobius{A: a, B: 0, C: 0, D: 1}
}

// RotationMap returns the rotation z -> e^{i*theta} * z about the origin.
func RotationMap(theta float64) Mobius {
	return ScalingMap(cmplx.Rect(1, theta))
}

// InversionMap returns the Mobius map z -> 1/z.
func InversionMap() Mobius {
	return Mobius{A: 0, B: 1, C: 1, D: 0}
}

// MobiusFromPoints returns the unique Mobius transformation that sends the three
// distinct points z1, z2, z3 to w1, w2, w3 respectively. It is built by
// composing the standard map taking z1, z2, z3 to 0, 1, infinity with the
// inverse of the analogous map for the w points.
func MobiusFromPoints(z1, z2, z3, w1, w2, w3 complex128) Mobius {
	return complexanalysisToStandard(w1, w2, w3).Inverse().Compose(complexanalysisToStandard(z1, z2, z3))
}

// complexanalysisToStandard returns the Mobius map sending p1, p2, p3 to
// 0, 1, infinity, whose formula is the cross-ratio (z-p1)(p2-p3)/((z-p3)(p2-p1)).
func complexanalysisToStandard(p1, p2, p3 complex128) Mobius {
	return Mobius{
		A: p2 - p3,
		B: -p1 * (p2 - p3),
		C: p2 - p1,
		D: -p3 * (p2 - p1),
	}
}

// CrossRatio returns the cross-ratio (z-z1)(z2-z3) / ((z-z3)(z2-z1)) of the four
// points. The cross-ratio is invariant under every Mobius transformation.
func CrossRatio(z, z1, z2, z3 complex128) complex128 {
	return ((z - z1) * (z2 - z3)) / ((z - z3) * (z2 - z1))
}

// CayleyTransform returns the value of the Cayley transform (z-i)/(z+i), which
// maps the upper half-plane conformally onto the open unit disk.
func CayleyTransform(z complex128) complex128 {
	i := complex(0, 1)
	return (z - i) / (z + i)
}

// InverseCayleyTransform returns i*(1+w)/(1-w), the inverse of the Cayley
// transform, mapping the unit disk back to the upper half-plane.
func InverseCayleyTransform(w complex128) complex128 {
	i := complex(0, 1)
	return i * (1 + w) / (1 - w)
}

// JoukowskiMap returns the Joukowski map (z + 1/z)/2, which maps circles about
// the origin to confocal ellipses and is used in classical airfoil theory. It
// returns a complex infinity at z = 0.
func JoukowskiMap(z complex128) complex128 {
	if z == 0 {
		return cmplx.Inf()
	}
	return (z + 1/z) / 2
}
