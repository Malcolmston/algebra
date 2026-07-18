package hypercomplex

import "math"

// DualQuaternion is a dual quaternion Real + Dual*epsilon, where epsilon is the
// dual unit with epsilon^2 = 0 and Real and Dual are ordinary [Quaternion]
// values. A unit dual quaternion encodes a rigid-body transform: its real part
// carries the rotation and its dual part carries the translation.
type DualQuaternion struct {
	Real Quaternion
	Dual Quaternion
}

// DualQuat constructs a dual quaternion from its real and dual quaternion parts.
func DualQuat(real, dual Quaternion) DualQuaternion {
	return DualQuaternion{Real: real, Dual: dual}
}

// IdentityDualQuat returns the identity dual quaternion (identity rotation, zero
// translation).
func IdentityDualQuat() DualQuaternion {
	return DualQuaternion{Real: IdentityQuat()}
}

// DualQuatFromRotationTranslation returns the unit dual quaternion for the rigid
// transform that first rotates by the (unit) quaternion rot and then translates
// by t. The rotation is normalized internally.
func DualQuatFromRotationTranslation(rot Quaternion, t Vec3) DualQuaternion {
	r := rot.Normalize()
	// Dual part = 0.5 * (0,t) * r.
	tq := QuatFromVector(t)
	dual := tq.Mul(r).Scale(0.5)
	return DualQuaternion{Real: r, Dual: dual}
}

// DualQuatFromTranslation returns the unit dual quaternion representing a pure
// translation by t with no rotation.
func DualQuatFromTranslation(t Vec3) DualQuaternion {
	return DualQuatFromRotationTranslation(IdentityQuat(), t)
}

// Add returns the component-wise sum of the two dual quaternions.
func (d DualQuaternion) Add(e DualQuaternion) DualQuaternion {
	return DualQuaternion{Real: d.Real.Add(e.Real), Dual: d.Dual.Add(e.Dual)}
}

// Scale returns d with both parts multiplied by the real factor s.
func (d DualQuaternion) Scale(s float64) DualQuaternion {
	return DualQuaternion{Real: d.Real.Scale(s), Dual: d.Dual.Scale(s)}
}

// Mul returns the dual quaternion product d*e. Because epsilon^2 = 0, the
// product is (Rd*Re) + (Rd*De + Dd*Re)*epsilon. Composing two rigid transforms
// corresponds to multiplying their unit dual quaternions.
func (d DualQuaternion) Mul(e DualQuaternion) DualQuaternion {
	return DualQuaternion{
		Real: d.Real.Mul(e.Real),
		Dual: d.Real.Mul(e.Dual).Add(d.Dual.Mul(e.Real)),
	}
}

// Conjugate returns the quaternion conjugate applied to both parts,
// (Real*, Dual*). This is the conjugate used to transform points via the
// sandwich product.
func (d DualQuaternion) Conjugate() DualQuaternion {
	return DualQuaternion{Real: d.Real.Conj(), Dual: d.Dual.Conj()}
}

// DualNumberConjugate returns the dual-number conjugate (Real, -Dual), which
// negates the dual part while leaving the real part unchanged.
func (d DualQuaternion) DualNumberConjugate() DualQuaternion {
	return DualQuaternion{Real: d.Real, Dual: d.Dual.Neg()}
}

// Norm returns the norm of a unit dual quaternion's real part, |Real|. For a
// valid rigid transform this equals 1.
func (d DualQuaternion) Norm() float64 {
	return d.Real.Norm()
}

// Normalize returns d rescaled so that its real part has unit norm, dividing
// both parts by |Real|. If the real part is zero, d is returned unchanged.
func (d DualQuaternion) Normalize() DualQuaternion {
	n := d.Real.Norm()
	if n == 0 {
		return d
	}
	return d.Scale(1 / n)
}

// Rotation returns the rotation component of d as a unit quaternion (its
// normalized real part).
func (d DualQuaternion) Rotation() Quaternion {
	return d.Real.Normalize()
}

// Translation returns the translation vector encoded by the unit dual
// quaternion d, recovered as 2 * Dual * Real*.
func (d DualQuaternion) Translation() Vec3 {
	r := d.Real.Normalize()
	// Scale the dual part consistently with the normalized real part.
	scale := 1.0
	if n := d.Real.Norm(); n != 0 {
		scale = 1 / n
	}
	dual := d.Dual.Scale(scale)
	t := dual.Mul(r.Conj()).Scale(2)
	return t.Vector()
}

// TransformPoint applies the rigid transform encoded by the unit dual
// quaternion d to the point p, returning the rotated-then-translated point.
func (d DualQuaternion) TransformPoint(p Vec3) Vec3 {
	rot := d.Rotation()
	return rot.RotateVector(p).Add(d.Translation())
}

// Equal reports whether d and e agree in both parts to within the absolute
// tolerance tol.
func (d DualQuaternion) Equal(e DualQuaternion, tol float64) bool {
	return d.Real.Equal(e.Real, tol) && d.Dual.Equal(e.Dual, tol)
}

// ScrewAngle returns the rotation angle in radians of the rigid transform
// encoded by d, in the range [0, pi]. It is a convenience wrapper over the
// axis-angle decomposition of the rotation part.
func (d DualQuaternion) ScrewAngle() float64 {
	_, angle := d.Rotation().ToAxisAngle()
	return math.Mod(angle, 2*math.Pi)
}
