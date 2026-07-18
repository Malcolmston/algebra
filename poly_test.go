package algebra

import (
	"testing"
)

// paReconstruct multiplies out a factor list (Base^Mult) and returns the
// resulting polynomial for comparison with the original.
func paReconstruct(v string, fs []PolyFactor) *Poly {
	acc := NewPoly(v, Int(1))
	for _, f := range fs {
		for k := 0; k < f.Mult; k++ {
			acc = acc.MulP(f.Base)
		}
	}
	return acc
}

// paEqualCoeffs reports whether two polynomials have identical coefficients.
func paEqualCoeffs(a, b *Poly) bool {
	if len(a.coeffs) != len(b.coeffs) {
		return false
	}
	for i := range a.coeffs {
		if a.coeffs[i].Cmp(b.coeffs[i]) != 0 {
			return false
		}
	}
	return true
}

func TestPolyBasics(t *testing.T) {
	p := NewPoly("x", Int(1), Int(-2), Int(1)) // 1 - 2x + x^2
	if got := p.String(); got != "x^2 - 2*x + 1" {
		t.Errorf("String = %q, want %q", got, "x^2 - 2*x + 1")
	}
	if p.Degree() != 2 {
		t.Errorf("Degree = %d, want 2", p.Degree())
	}
	if p.Var() != "x" {
		t.Errorf("Var = %q, want x", p.Var())
	}
	if !p.Coeff(1).Equal(Int(-2)) {
		t.Errorf("Coeff(1) = %v, want -2", p.Coeff(1))
	}
	if !p.Coeff(5).Equal(Int(0)) {
		t.Errorf("Coeff(5) = %v, want 0", p.Coeff(5))
	}
	if !p.LeadingCoeff().Equal(Int(1)) {
		t.Errorf("LeadingCoeff = %v, want 1", p.LeadingCoeff())
	}
	if p.IsZero() {
		t.Error("p reported zero")
	}
	z := NewPoly("x", Int(0), Int(0))
	if !z.IsZero() || z.Degree() != -1 || z.String() != "0" {
		t.Errorf("zero poly wrong: zero=%v deg=%d str=%q", z.IsZero(), z.Degree(), z.String())
	}
}

func TestPolyEval(t *testing.T) {
	p := NewPoly("x", Int(1), Int(-2), Int(1)) // (x-1)^2
	cases := []struct {
		x, want Expr
	}{
		{Int(3), Int(4)},
		{Int(1), Int(0)},
		{Int(0), Int(1)},
		{Rat(1, 2), Rat(1, 4)},
	}
	for _, c := range cases {
		if got := p.Eval(c.x); !got.Equal(c.want) {
			t.Errorf("Eval(%v) = %v, want %v", c.x, got, c.want)
		}
	}
}

func TestPolyArithmetic(t *testing.T) {
	a := NewPoly("x", Int(1), Int(1))  // x + 1
	b := NewPoly("x", Int(-1), Int(1)) // x - 1

	if got := a.MulP(b); !paEqualCoeffs(got, NewPoly("x", Int(-1), Int(0), Int(1))) {
		t.Errorf("MulP = %s, want x^2 - 1", got)
	}
	if got := a.Add(b); !paEqualCoeffs(got, NewPoly("x", Int(0), Int(2))) {
		t.Errorf("Add = %s, want 2*x", got)
	}
	if got := a.Sub(b); !paEqualCoeffs(got, NewPoly("x", Int(2))) {
		t.Errorf("Sub = %s, want 2", got)
	}
	if got := a.Scale(Int(3)); !paEqualCoeffs(got, NewPoly("x", Int(3), Int(3))) {
		t.Errorf("Scale = %s, want 3*x + 3", got)
	}
}

func TestPolyDivMod(t *testing.T) {
	// (x^2 - 1) / (x - 1) = x + 1 remainder 0.
	num := NewPoly("x", Int(-1), Int(0), Int(1))
	den := NewPoly("x", Int(-1), Int(1))
	quo, rem, err := num.DivMod(den)
	if err != nil {
		t.Fatal(err)
	}
	if !paEqualCoeffs(quo, NewPoly("x", Int(1), Int(1))) || !rem.IsZero() {
		t.Errorf("DivMod = (%s, %s), want (x + 1, 0)", quo, rem)
	}
	// (x^2 + 1) / (x - 1) = x + 1 remainder 2.
	num2 := NewPoly("x", Int(1), Int(0), Int(1))
	quo2, rem2, _ := num2.DivMod(den)
	if !paEqualCoeffs(quo2, NewPoly("x", Int(1), Int(1))) || !paEqualCoeffs(rem2, NewPoly("x", Int(2))) {
		t.Errorf("DivMod = (%s, %s), want (x + 1, 2)", quo2, rem2)
	}
	// Division by zero polynomial errors.
	if _, _, err := num.DivMod(NewPoly("x", Int(0))); err == nil {
		t.Error("DivMod by zero did not error")
	}
}

func TestPolyDerivativeMonic(t *testing.T) {
	p := NewPoly("x", Int(1), Int(-2), Int(1)) // x^2 - 2x + 1
	if got := p.Derivative(); !paEqualCoeffs(got, NewPoly("x", Int(-2), Int(2))) {
		t.Errorf("Derivative = %s, want 2*x - 2", got)
	}
	q := NewPoly("x", Int(4), Int(2)) // 2x + 4
	if got := q.Monic(); !paEqualCoeffs(got, NewPoly("x", Int(2), Int(1))) {
		t.Errorf("Monic = %s, want x + 2", got)
	}
}

func TestPolyGCDLCM(t *testing.T) {
	// gcd((x-1)(x+1), (x-1)^2) = x - 1.
	a := NewPoly("x", Int(-1), Int(0), Int(1)) // x^2 - 1
	b := NewPoly("x", Int(1), Int(-2), Int(1)) // x^2 - 2x + 1
	if got := PolyGCD(a, b); !paEqualCoeffs(got, NewPoly("x", Int(-1), Int(1))) {
		t.Errorf("PolyGCD = %s, want x - 1", got)
	}
	// lcm(x-1, x+1) = x^2 - 1.
	l := PolyLCM(NewPoly("x", Int(-1), Int(1)), NewPoly("x", Int(1), Int(1)))
	if !paEqualCoeffs(l, NewPoly("x", Int(-1), Int(0), Int(1))) {
		t.Errorf("PolyLCM = %s, want x^2 - 1", l)
	}
}

func TestPolySquareFree(t *testing.T) {
	// x^3 - 3x + 2 = (x-1)^2 (x+2).
	p := NewPoly("x", Int(2), Int(-3), Int(0), Int(1))
	fs := p.SquareFree()
	if !paEqualCoeffs(paReconstruct("x", fs), p) {
		t.Errorf("SquareFree reconstruction mismatch: %v", fs)
	}
	want := map[string]int{"x - 1": 2, "x + 2": 1}
	got := map[string]int{}
	for _, f := range fs {
		got[f.Base.String()] = f.Mult
	}
	for k, v := range want {
		if got[k] != v {
			t.Errorf("SquareFree factor %q mult = %d, want %d (all=%v)", k, got[k], v, got)
		}
	}
}

func TestPolyFactor(t *testing.T) {
	// x^3 - 3x + 2 = (x-1)^2 (x+2).
	p := NewPoly("x", Int(2), Int(-3), Int(0), Int(1))
	fs := p.Factor()
	if !paEqualCoeffs(paReconstruct("x", fs), p) {
		t.Errorf("Factor reconstruction mismatch: %v", fs)
	}
	got := map[string]int{}
	for _, f := range fs {
		got[f.Base.String()] = f.Mult
	}
	if got["x - 1"] != 2 || got["x + 2"] != 1 {
		t.Errorf("Factor = %v, want (x-1)^2 (x+2)", got)
	}

	// 2x^2 - 2 = 2 (x-1)(x+1): leading coefficient returned as a constant factor.
	q := NewPoly("x", Int(-2), Int(0), Int(2))
	qf := q.Factor()
	if !paEqualCoeffs(paReconstruct("x", qf), q) {
		t.Errorf("Factor reconstruction mismatch for 2x^2-2: %v", qf)
	}

	// x^2 - 2 is irreducible over the rationals and returned as-is.
	irr := NewPoly("x", Int(-2), Int(0), Int(1))
	irf := irr.Factor()
	if len(irf) != 1 || irf[0].Base.String() != "x^2 - 2" || irf[0].Mult != 1 {
		t.Errorf("Factor(x^2-2) = %v, want single (x^2 - 2)", irf)
	}
}

func TestPolyResultant(t *testing.T) {
	// res(x-2, x-3) = -1.
	r := NewPoly("x", Int(-2), Int(1)).Resultant(NewPoly("x", Int(-3), Int(1)))
	if !r.Equal(Int(-1)) {
		t.Errorf("Resultant(x-2, x-3) = %v, want -1", r)
	}
	// Shared factor => resultant 0.
	shared := NewPoly("x", Int(-1), Int(0), Int(1)).Resultant(NewPoly("x", Int(-1), Int(1)))
	if !shared.Equal(Int(0)) {
		t.Errorf("Resultant with shared factor = %v, want 0", shared)
	}
}

func TestPolyDiscriminant(t *testing.T) {
	cases := []struct {
		p    *Poly
		want Expr
	}{
		{NewPoly("x", Int(6), Int(-5), Int(1)), Int(1)},  // x^2 - 5x + 6, disc = 1
		{NewPoly("x", Int(1), Int(0), Int(1)), Int(-4)},  // x^2 + 1, disc = -4
		{NewPoly("x", Int(-4), Int(0), Int(1)), Int(16)}, // x^2 - 4, disc = 16
	}
	for _, c := range cases {
		if got := c.p.Discriminant(); !got.Equal(c.want) {
			t.Errorf("Discriminant(%s) = %v, want %v", c.p, got, c.want)
		}
	}
}

func TestPartialFractionsDistinct(t *testing.T) {
	// 1 / ((x-1)(x-2)) = -1/(x-1) + 1/(x-2).
	num := NewPoly("x", Int(1))
	den := NewPoly("x", Int(2), Int(-3), Int(1)) // x^2 - 3x + 2
	poly, terms, err := PartialFractions(num, den)
	if err != nil {
		t.Fatal(err)
	}
	if !poly.IsZero() {
		t.Errorf("polynomial part = %s, want 0", poly)
	}
	if len(terms) != 2 {
		t.Fatalf("got %d terms, want 2: %+v", len(terms), terms)
	}
	if terms[0].Denom.String() != "x - 1" || !terms[0].Coeff.Equal(Int(-1)) || terms[0].Power != 1 {
		t.Errorf("term0 = %v/(%s)^%d, want -1/(x-1)", terms[0].Coeff, terms[0].Denom, terms[0].Power)
	}
	if terms[1].Denom.String() != "x - 2" || !terms[1].Coeff.Equal(Int(1)) || terms[1].Power != 1 {
		t.Errorf("term1 = %v/(%s)^%d, want 1/(x-2)", terms[1].Coeff, terms[1].Denom, terms[1].Power)
	}
}

func TestPartialFractionsRepeated(t *testing.T) {
	// 1 / (x-1)^2 = 1/(x-1)^2.
	num := NewPoly("x", Int(1))
	den := NewPoly("x", Int(1), Int(-2), Int(1)) // (x-1)^2
	poly, terms, err := PartialFractions(num, den)
	if err != nil {
		t.Fatal(err)
	}
	if !poly.IsZero() {
		t.Errorf("polynomial part = %s, want 0", poly)
	}
	if len(terms) != 1 {
		t.Fatalf("got %d terms, want 1: %+v", len(terms), terms)
	}
	if terms[0].Denom.String() != "x - 1" || terms[0].Power != 2 || !terms[0].Coeff.Equal(Int(1)) {
		t.Errorf("term = %v/(%s)^%d, want 1/(x-1)^2", terms[0].Coeff, terms[0].Denom, terms[0].Power)
	}
}

func TestPartialFractionsImproper(t *testing.T) {
	// (x^2) / (x^2 - 1) = 1 + 1/(x^2-1) split into partial fractions.
	num := NewPoly("x", Int(0), Int(0), Int(1))
	den := NewPoly("x", Int(-1), Int(0), Int(1))
	poly, terms, err := PartialFractions(num, den)
	if err != nil {
		t.Fatal(err)
	}
	if !paEqualCoeffs(poly, NewPoly("x", Int(1))) {
		t.Errorf("polynomial part = %s, want 1", poly)
	}
	if len(terms) != 2 {
		t.Fatalf("got %d terms, want 2", len(terms))
	}
}

func TestApartExpr(t *testing.T) {
	x := Sym("x")
	// e = 1 / ((x-1)(x-2)).
	e := Mul(
		Pow(Add(x, Int(-1)), Int(-1)),
		Pow(Add(x, Int(-2)), Int(-1)),
	)
	got, err := ApartExpr(e, x)
	if err != nil {
		t.Fatal(err)
	}
	// The decomposition must be numerically equal to e at several points.
	for _, xv := range []int64{3, 4, 5, 10} {
		want, err1 := Evalf(Subs(e, x, Int(xv)))
		have, err2 := Evalf(Subs(got, x, Int(xv)))
		if err1 != nil || err2 != nil {
			t.Fatalf("eval error: %v %v", err1, err2)
		}
		if diff := want - have; diff > 1e-9 || diff < -1e-9 {
			t.Errorf("ApartExpr at x=%d: got %v, want %v", xv, have, want)
		}
	}
}

func TestPolyFromRoundTrip(t *testing.T) {
	x := Sym("x")
	e := Add(Pow(x, Int(2)), Mul(Int(-2), x), Int(1)) // x^2 - 2x + 1
	p, err := PolyFrom(e, x)
	if err != nil {
		t.Fatal(err)
	}
	if !paEqualCoeffs(p, NewPoly("x", Int(1), Int(-2), Int(1))) {
		t.Errorf("PolyFrom = %s, want x^2 - 2x + 1", p)
	}
	// Round trip through Expr and back.
	p2, err := PolyFrom(p.Expr(), x)
	if err != nil {
		t.Fatal(err)
	}
	if !paEqualCoeffs(p, p2) {
		t.Errorf("Expr round trip mismatch: %s vs %s", p, p2)
	}
}

func TestPolyFromErrors(t *testing.T) {
	x := Sym("x")
	// Non-polynomial (negative power) must error.
	if _, err := PolyFrom(Pow(x, Int(-1)), x); err == nil {
		t.Error("PolyFrom(x^-1) did not error")
	}
	// Non-rational coefficient must error.
	if _, err := PolyFrom(Mul(Pi, x), x); err == nil {
		t.Error("PolyFrom(pi*x) did not error")
	}
	// Non-symbol target must error.
	if _, err := PolyFrom(x, Int(2)); err == nil {
		t.Error("PolyFrom with non-symbol did not error")
	}
}

func TestPolyExprBuild(t *testing.T) {
	// Verify Expr uses exact rational coefficients.
	p := NewPoly("x", Rat(1, 2), Int(0), Int(1)) // x^2 + 1/2
	want := Add(Pow(Sym("x"), Int(2)), Rat(1, 2))
	if !Simplify(p.Expr()).Equal(Simplify(want)) {
		t.Errorf("Expr = %v, want %v", p.Expr(), want)
	}
}
