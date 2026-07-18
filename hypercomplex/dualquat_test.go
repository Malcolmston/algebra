package hypercomplex

import (
	"math"
	"testing"
)

const dTol = 1e-9

func TestDualQuatPureTranslation(t *testing.T) {
	d := DualQuatFromTranslation(V3(1, 2, 3))
	if got := d.Translation(); !got.Equal(V3(1, 2, 3), dTol) {
		t.Errorf("translation = %+v want (1,2,3)", got)
	}
	// A pure translation leaves orientation unchanged.
	if got := d.Rotation(); !got.Equal(IdentityQuat(), dTol) {
		t.Errorf("rotation = %+v want identity", got)
	}
	// Transforming a point just adds the translation.
	if got := d.TransformPoint(V3(5, 5, 5)); !got.Equal(V3(6, 7, 8), dTol) {
		t.Errorf("transform = %+v want (6,7,8)", got)
	}
}

func TestDualQuatRotationTranslation(t *testing.T) {
	rot := QuatFromAxisAngle(V3(0, 0, 1), math.Pi/2)
	tr := V3(1, 0, 0)
	d := DualQuatFromRotationTranslation(rot, tr)
	// Recover rotation and translation.
	if got := d.Rotation(); !got.Equal(rot, dTol) && !got.Equal(rot.Neg(), dTol) {
		t.Errorf("rotation = %+v want %+v", got, rot)
	}
	if got := d.Translation(); !got.Equal(tr, dTol) {
		t.Errorf("translation = %+v want %+v", got, tr)
	}
	// Transform x-axis: rotate (1,0,0)->(0,1,0), then translate +x -> (1,1,0).
	if got := d.TransformPoint(V3(1, 0, 0)); !got.Equal(V3(1, 1, 0), dTol) {
		t.Errorf("transform = %+v want (1,1,0)", got)
	}
}

func TestDualQuatCompose(t *testing.T) {
	// Composition of two transforms equals transforming sequentially.
	d1 := DualQuatFromRotationTranslation(QuatFromAxisAngle(V3(0, 0, 1), math.Pi/2), V3(1, 0, 0))
	d2 := DualQuatFromTranslation(V3(0, 0, 2))
	comp := d1.Mul(d2)
	p := V3(2, 0, 0)
	// d1 applied after d2: first d2 then d1 under this multiplication order.
	want := d1.TransformPoint(d2.TransformPoint(p))
	got := comp.TransformPoint(p)
	if !got.Equal(want, 1e-9) {
		t.Errorf("composed transform = %+v want %+v", got, want)
	}
}

func TestDualQuatIdentity(t *testing.T) {
	id := IdentityDualQuat()
	p := V3(7, -3, 2)
	if got := id.TransformPoint(p); !got.Equal(p, dTol) {
		t.Errorf("identity transform = %+v want %+v", got, p)
	}
	if math.Abs(id.Norm()-1) > dTol {
		t.Errorf("identity norm = %v want 1", id.Norm())
	}
}

func TestDualQuatConjugates(t *testing.T) {
	d := DualQuatFromRotationTranslation(QuatFromAxisAngle(V3(1, 1, 0), 0.5), V3(1, 2, 3))
	// Applying the transform then its inverse (quaternion conjugate composition)
	// returns the identity rotation; check the dual-number conjugate negates dual.
	dc := d.DualNumberConjugate()
	if !dc.Dual.Equal(d.Dual.Neg(), dTol) {
		t.Errorf("dual-number conjugate dual = %+v want negation", dc.Dual)
	}
	if !dc.Real.Equal(d.Real, dTol) {
		t.Errorf("dual-number conjugate real changed")
	}
}
