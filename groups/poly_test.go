package groups

import (
	"math"
	"testing"
)

func TestPolyArithmetic(t *testing.T) {
	a := Poly{1, 2, 3} // 3x^2+2x+1
	b := Poly{0, 1}    // x
	if !PolyEqual(PolyAdd(a, b), Poly{1, 3, 3}) {
		t.Errorf("PolyAdd wrong: %v", PolyAdd(a, b))
	}
	if !PolyEqual(PolySub(a, a), Poly{}) {
		t.Errorf("PolySub self should be zero")
	}
	// (x+1)(x-1) = x^2 - 1
	prod := PolyMul(Poly{1, 1}, Poly{-1, 1})
	if !PolyEqual(prod, Poly{-1, 0, 1}) {
		t.Errorf("PolyMul=%v want x^2-1", prod)
	}
	if a.Degree() != 2 {
		t.Errorf("Degree want 2")
	}
}

func TestPolyDivMod(t *testing.T) {
	// (x^2 - 1) / (x - 1) = x + 1, remainder 0
	q, r := PolyDivMod(Poly{-1, 0, 1}, Poly{-1, 1})
	if !PolyEqual(q, Poly{1, 1}) || !r.IsZero() {
		t.Errorf("DivMod q=%v r=%v want x+1, 0", q, r)
	}
	// x^2 + 1 divided by x gives q=x, r=1
	q, r = PolyDivMod(Poly{1, 0, 1}, Poly{0, 1})
	if !PolyEqual(q, Poly{0, 1}) || !PolyEqual(r, Poly{1}) {
		t.Errorf("DivMod q=%v r=%v want x, 1", q, r)
	}
	// verify a = q*b + r on a random-ish case
	a := Poly{2, -3, 0, 1} // x^3 - 3x + 2
	b := Poly{-1, 1}       // x - 1
	q, r = PolyDivMod(a, b)
	recon := PolyAdd(PolyMul(q, b), r)
	if !PolyEqual(recon, a) {
		t.Errorf("a != q*b+r: %v", recon)
	}
}

func TestPolyGCD(t *testing.T) {
	// gcd(x^2-1, x^2-2x+1) = x-1 (monic)
	a := Poly{-1, 0, 1} // x^2-1
	b := Poly{1, -2, 1} // (x-1)^2
	g := PolyGCD(a, b)
	if !PolyEqual(g, Poly{-1, 1}) {
		t.Errorf("PolyGCD=%v want x-1", g)
	}
	// gcd with monic normalization: gcd(2x^2-2, x-1) = x-1
	g2 := PolyGCD(Poly{-2, 0, 2}, Poly{-1, 1})
	if !PolyEqual(g2, Poly{-1, 1}) {
		t.Errorf("PolyGCD monic=%v want x-1", g2)
	}
}

func TestPolyEvalDeriv(t *testing.T) {
	a := Poly{1, 2, 3} // 3x^2+2x+1
	if math.Abs(PolyEval(a, 2)-17) > 1e-9 {
		t.Errorf("PolyEval(2) want 17 got %v", PolyEval(a, 2))
	}
	d := PolyDerivative(a) // 6x+2
	if !PolyEqual(d, Poly{2, 6}) {
		t.Errorf("PolyDerivative=%v want 6x+2", d)
	}
	if PolyDerivative(Poly{5}) != nil && !PolyDerivative(Poly{5}).IsZero() {
		t.Errorf("derivative of constant should be zero")
	}
}

func TestPolyMonicString(t *testing.T) {
	m := PolyMonic(Poly{2, 4}) // 4x+2 -> x+0.5
	if math.Abs(m.LeadingCoeff()-1) > 1e-9 {
		t.Errorf("PolyMonic leading coeff want 1")
	}
	if (Poly{}).String() != "0" {
		t.Errorf("zero poly string")
	}
	if (Poly{-1, 0, 1}).String() != "x^2 - 1" {
		t.Errorf("string got %q", (Poly{-1, 0, 1}).String())
	}
}
