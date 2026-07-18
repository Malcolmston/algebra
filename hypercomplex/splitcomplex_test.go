package hypercomplex

import (
	"math"
	"testing"
)

const sTol = 1e-9

func TestSplitComplexMul(t *testing.T) {
	// j*j = 1.
	j := Split(0, 1)
	if got := j.Mul(j); !got.Equal(Split(1, 0), sTol) {
		t.Errorf("j*j = %+v want 1", got)
	}
	// (2+3j)(1+4j) = (2+12) + (8+3)j = 14 + 11j.
	got := Split(2, 3).Mul(Split(1, 4))
	if !got.Equal(Split(14, 11), sTol) {
		t.Errorf("(2+3j)(1+4j) = %+v want 14+11j", got)
	}
}

func TestSplitComplexModulusMultiplicative(t *testing.T) {
	z := Split(3, 1)
	w := Split(2, 5)
	// ModulusSq is multiplicative: N(zw) = N(z)N(w).
	lhs := z.Mul(w).ModulusSq()
	rhs := z.ModulusSq() * w.ModulusSq()
	if math.Abs(lhs-rhs) > sTol {
		t.Errorf("N(zw) = %v want %v", lhs, rhs)
	}
	// N(3+1j) = 9 - 1 = 8.
	if got := z.ModulusSq(); math.Abs(got-8) > sTol {
		t.Errorf("N(z) = %v want 8", got)
	}
}

func TestSplitComplexInverse(t *testing.T) {
	z := Split(3, 1)
	inv, ok := z.Inverse()
	if !ok {
		t.Fatal("expected invertible")
	}
	if got := z.Mul(inv); !got.Equal(Split(1, 0), sTol) {
		t.Errorf("z*z^-1 = %+v want 1", got)
	}
	// Null element 1+j is a zero divisor.
	null := Split(1, 1)
	if !null.IsZeroDivisor(sTol) {
		t.Errorf("1+j should be a zero divisor")
	}
	if _, ok := null.Inverse(); ok {
		t.Errorf("1+j should not be invertible")
	}
}

func TestSplitComplexExp(t *testing.T) {
	// exp(0 + phi j) = cosh phi + j sinh phi.
	phi := 0.75
	got := Split(0, phi).Exp()
	want := Split(math.Cosh(phi), math.Sinh(phi))
	if !got.Equal(want, sTol) {
		t.Errorf("exp(phi j) = %+v want %+v", got, want)
	}
	// exp preserves modulus 1 on the boost axis: N(exp(phi j)) = 1.
	if m := got.ModulusSq(); math.Abs(m-1) > sTol {
		t.Errorf("N(exp(phi j)) = %v want 1", m)
	}
}

func TestSplitComplexModulusArgumentRoundTrip(t *testing.T) {
	z := SplitFromModulusArgument(2.5, 0.6)
	if math.Abs(z.Abs()-2.5) > sTol {
		t.Errorf("Abs = %v want 2.5", z.Abs())
	}
	if math.Abs(z.Argument()-0.6) > sTol {
		t.Errorf("Argument = %v want 0.6", z.Argument())
	}
}

func TestSplitComplexLightCone(t *testing.T) {
	z := Split(3, 1)
	u, v := z.LightCone()
	if math.Abs(u-4) > sTol || math.Abs(v-2) > sTol {
		t.Errorf("lightcone = (%v,%v) want (4,2)", u, v)
	}
	// Multiplication is component-wise in the null basis.
	w := Split(2, 5)
	u2, v2 := w.LightCone()
	pu, pv := z.Mul(w).LightCone()
	if math.Abs(pu-u*u2) > sTol || math.Abs(pv-v*v2) > sTol {
		t.Errorf("null-basis product = (%v,%v) want (%v,%v)", pu, pv, u*u2, v*v2)
	}
	// Round trip.
	if got := SplitFromLightCone(u, v); !got.Equal(z, sTol) {
		t.Errorf("lightcone round trip = %+v want %+v", got, z)
	}
}
