package complexanalysis

import (
	"math/cmplx"
	"testing"
)

func TestMobiusBasics(t *testing.T) {
	id := IdentityMobius()
	if !closeC(id.Apply(complex(3, 4)), complex(3, 4), testTol) {
		t.Error("identity should fix every point")
	}
	m := NewMobius(2, 1, 1, 3) // (2z+1)/(z+3)
	if !closeC(m.Apply(0), complex(1.0/3, 0), testTol) {
		t.Error("Apply wrong")
	}
	if !closeC(m.Determinant(), 5, testTol) {
		t.Error("Determinant wrong")
	}
	// m composed with its inverse is the identity map.
	comp := m.Compose(m.Inverse())
	z := complex(0.3, 1.2)
	if !closeC(comp.Apply(z), z, 1e-9) {
		t.Error("m . m^-1 should be identity map")
	}
	// Normalize keeps the action but sets determinant magnitude to 1.
	nm := m.Normalize()
	if !closeC(nm.Determinant(), 1, 1e-9) {
		t.Error("Normalize determinant should be 1")
	}
	if !closeC(nm.Apply(z), m.Apply(z), 1e-9) {
		t.Error("Normalize must preserve the map")
	}
}

func TestMobiusFromPoints(t *testing.T) {
	m := MobiusFromPoints(0, 1, -1, 1, complex(0, 1), complex(0, -1))
	if !closeC(m.Apply(0), 1, 1e-9) {
		t.Error("0 should map to 1")
	}
	if !closeC(m.Apply(1), complex(0, 1), 1e-9) {
		t.Error("1 should map to i")
	}
	if !closeC(m.Apply(-1), complex(0, -1), 1e-9) {
		t.Error("-1 should map to -i")
	}
}

func TestCrossRatioInvariance(t *testing.T) {
	z, z1, z2, z3 := complex(0.5, 0.5), complex(0, 0), complex(1, 0), complex(0, 1)
	m := NewMobius(complex(1, 1), 2, complex(0, 1), 3)
	before := CrossRatio(z, z1, z2, z3)
	after := CrossRatio(m.Apply(z), m.Apply(z1), m.Apply(z2), m.Apply(z3))
	if !closeC(before, after, 1e-9) {
		t.Errorf("cross-ratio not invariant: %v vs %v", before, after)
	}
}

func TestCayleyAndJoukowski(t *testing.T) {
	if !closeC(CayleyTransform(complex(0, 1)), 0, testTol) {
		t.Error("Cayley(i) should be 0")
	}
	if !closeC(InverseCayleyTransform(0), complex(0, 1), testTol) {
		t.Error("InverseCayley(0) should be i")
	}
	// Round trip.
	z := complex(0.4, 2.0)
	if !closeC(InverseCayleyTransform(CayleyTransform(z)), z, 1e-9) {
		t.Error("Cayley round trip failed")
	}
	if !closeC(JoukowskiMap(1), 1, testTol) {
		t.Error("Joukowski(1) should be 1")
	}
	if !closeC(JoukowskiMap(complex(0, 1)), 0, testTol) {
		t.Error("Joukowski(i) should be 0")
	}
}

func TestFixedPointsAndHelpers(t *testing.T) {
	// z -> 1/z has fixed points +/-1.
	fp := InversionMap().FixedPoints()
	if len(fp) != 2 {
		t.Fatalf("expected 2 fixed points, got %d", len(fp))
	}
	got1, got2 := fp[0], fp[1]
	if !((closeC(got1, 1, 1e-9) && closeC(got2, -1, 1e-9)) ||
		(closeC(got1, -1, 1e-9) && closeC(got2, 1, 1e-9))) {
		t.Errorf("fixed points wrong: %v", fp)
	}
	if !closeC(TranslationMap(3).Apply(1), 4, testTol) {
		t.Error("TranslationMap wrong")
	}
	if !closeC(ScalingMap(2).Apply(complex(1, 1)), complex(2, 2), testTol) {
		t.Error("ScalingMap wrong")
	}
	if !closeC(RotationMap(0).Apply(complex(1, 0)), 1, testTol) {
		t.Error("RotationMap wrong")
	}
	if !closeC(InversionMap().Apply(2), 0.5, testTol) {
		t.Error("InversionMap wrong")
	}
	_ = cmplx.Abs
}
