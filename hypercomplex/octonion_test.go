package hypercomplex

import (
	"math"
	"testing"
)

const oTol = 1e-9

func TestOctonionNormMultiplicative(t *testing.T) {
	o := Oct(1, -2, 3, 0.5, -1.5, 2, 0.25, -0.75)
	p := Oct(0.3, 1, -1, 2, 0.5, -0.5, 1.25, 3)
	lhs := o.Mul(p).Norm()
	rhs := o.Norm() * p.Norm()
	if math.Abs(lhs-rhs) > 1e-9 {
		t.Errorf("|o*p| = %v want |o||p| = %v", lhs, rhs)
	}
}

func TestOctonionConjInverse(t *testing.T) {
	o := Oct(1, 2, 3, 4, 5, 6, 7, 8)
	// o * conj(o) = |o|^2 (real).
	prod := o.Mul(o.Conj())
	want := OctFromScalar(o.NormSq())
	if !prod.Equal(want, 1e-9) {
		t.Errorf("o*conj(o) = %+v want %+v", prod, want)
	}
	// o * o^-1 = 1.
	if got := o.Mul(o.Inverse()); !got.Equal(IdentityOct(), 1e-9) {
		t.Errorf("o*o^-1 = %+v want identity", got)
	}
}

func TestOctonionUnitBasisProducts(t *testing.T) {
	// e1*e1 = -1 (imaginary units square to -1).
	e1 := Oct(0, 1, 0, 0, 0, 0, 0, 0)
	if got := e1.Mul(e1); !got.Equal(OctFromScalar(-1), oTol) {
		t.Errorf("e1*e1 = %+v want -1", got)
	}
	// Anticommutativity of distinct imaginary units: e1*e2 = -(e2*e1).
	e2 := Oct(0, 0, 1, 0, 0, 0, 0, 0)
	a := e1.Mul(e2)
	b := e2.Mul(e1)
	if !a.Equal(b.Neg(), oTol) {
		t.Errorf("e1*e2 = %+v but e2*e1 = %+v; not anticommutative", a, b)
	}
}

func TestOctonionAlternativeButNotAssociative(t *testing.T) {
	o := Oct(1, 2, 0, -1, 0.5, 0, 3, 1)
	p := Oct(0, 1, -2, 0.5, 1, 2, 0, -1)
	q := Oct(2, 0, 1, 1, -1, 0.5, 2, 0)
	// Alternative: associator vanishes when two arguments coincide.
	if got := Associator(o, o, p); !got.Equal(Octonion{}, 1e-9) {
		t.Errorf("associator(o,o,p) = %+v want 0", got)
	}
	if got := Associator(o, p, p); !got.Equal(Octonion{}, 1e-9) {
		t.Errorf("associator(o,p,p) = %+v want 0", got)
	}
	// Not associative in general: pick a triple with non-zero associator.
	assoc := Associator(o, p, q)
	if assoc.Equal(Octonion{}, 1e-9) {
		t.Errorf("expected non-zero associator for independent triple")
	}
}

func TestOctonionCommutator(t *testing.T) {
	e1 := Oct(0, 1, 0, 0, 0, 0, 0, 0)
	e2 := Oct(0, 0, 1, 0, 0, 0, 0, 0)
	// [e1,e2] = e1*e2 - e2*e1 = 2*(e1*e2) since anticommuting.
	c := Commutator(e1, e2)
	want := e1.Mul(e2).Scale(2)
	if !c.Equal(want, oTol) {
		t.Errorf("commutator = %+v want %+v", c, want)
	}
}

func BenchmarkOctonionMul(b *testing.B) {
	o := Oct(1, -2, 3, 0.5, -1.5, 2, 0.25, -0.75)
	p := Oct(0.3, 1, -1, 2, 0.5, -0.5, 1.25, 3)
	var sink Octonion
	for i := 0; i < b.N; i++ {
		sink = o.Mul(p)
	}
	_ = sink
}
