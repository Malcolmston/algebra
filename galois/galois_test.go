package galois

import (
	"fmt"
	"math/big"
	"math/rand"
	"testing"
)

func bigEq(a, b *big.Int) bool { return a.Cmp(b) == 0 }

// ---- prime-field modular arithmetic ----

func TestModularArithmetic(t *testing.T) {
	p := big.NewInt(7)
	tests := []struct {
		name string
		got  *big.Int
		want int64
	}{
		{"add", AddMod(big.NewInt(5), big.NewInt(4), p), 2},
		{"sub", SubMod(big.NewInt(2), big.NewInt(5), p), 4},
		{"mul", MulMod(big.NewInt(3), big.NewInt(6), p), 4},
		{"neg", NegMod(big.NewInt(3), p), 4},
	}
	for _, tc := range tests {
		if !bigEq(tc.got, big.NewInt(tc.want)) {
			t.Errorf("%s: got %v want %d", tc.name, tc.got, tc.want)
		}
	}
}

func TestInvAndPowMod(t *testing.T) {
	p := big.NewInt(13)
	for a := int64(1); a < 13; a++ {
		inv, err := InvMod(big.NewInt(a), p)
		if err != nil {
			t.Fatalf("InvMod(%d): %v", a, err)
		}
		if !bigEq(MulMod(big.NewInt(a), inv, p), big1) {
			t.Errorf("a*inv != 1 for a=%d", a)
		}
	}
	got, err := PowMod(big.NewInt(2), big.NewInt(10), big.NewInt(13))
	if err != nil {
		t.Fatal(err)
	}
	if !bigEq(got, big.NewInt(10)) { // 2^10 = 1024 = 78*13+10
		t.Errorf("2^10 mod 13 = %v want 10", got)
	}
	// negative exponent
	neg, err := PowMod(big.NewInt(2), big.NewInt(-1), big.NewInt(13))
	if err != nil {
		t.Fatal(err)
	}
	if !bigEq(neg, big.NewInt(7)) { // inverse of 2 mod 13 is 7
		t.Errorf("2^-1 mod 13 = %v want 7", neg)
	}
}

func TestLegendreAndSqrt(t *testing.T) {
	p := big.NewInt(7)
	cases := []struct {
		a    int64
		leg  int
		isQR bool
	}{
		{1, 1, true}, {2, 1, true}, {3, -1, false}, {4, 1, true}, {0, 0, true},
	}
	for _, c := range cases {
		if got := Legendre(big.NewInt(c.a), p); got != c.leg {
			t.Errorf("Legendre(%d,7)=%d want %d", c.a, got, c.leg)
		}
		if got := IsQuadraticResidue(big.NewInt(c.a), p); got != c.isQR {
			t.Errorf("IsQuadraticResidue(%d,7)=%v want %v", c.a, got, c.isQR)
		}
	}
	// square roots for several primes
	for _, p := range []int64{7, 11, 13, 17, 101} {
		P := big.NewInt(p)
		for a := int64(0); a < p; a++ {
			if !IsQuadraticResidue(big.NewInt(a), P) {
				continue
			}
			r, err := SqrtMod(big.NewInt(a), P)
			if err != nil {
				t.Fatalf("SqrtMod(%d,%d): %v", a, p, err)
			}
			if !bigEq(MulMod(r, r, P), big.NewInt(a)) {
				t.Errorf("SqrtMod(%d,%d)=%v, square=%v", a, p, r, MulMod(r, r, P))
			}
		}
	}
}

func TestOrderAndPrimitiveRoot(t *testing.T) {
	ord, err := MultiplicativeOrder(big.NewInt(2), big.NewInt(7))
	if err != nil {
		t.Fatal(err)
	}
	if !bigEq(ord, big.NewInt(3)) {
		t.Errorf("ord_7(2)=%v want 3", ord)
	}
	g, err := PrimitiveRoot(big.NewInt(7))
	if err != nil {
		t.Fatal(err)
	}
	if !bigEq(g, big.NewInt(3)) {
		t.Errorf("PrimitiveRoot(7)=%v want 3", g)
	}
	if !IsPrimitiveRoot(big.NewInt(3), big.NewInt(7)) {
		t.Error("3 should be a primitive root mod 7")
	}
	if IsPrimitiveRoot(big.NewInt(2), big.NewInt(7)) {
		t.Error("2 should not be a primitive root mod 7")
	}
}

// ---- integer utilities ----

func TestIntUtilities(t *testing.T) {
	if got := EulerPhi(big.NewInt(12)); !bigEq(got, big.NewInt(4)) {
		t.Errorf("phi(12)=%v want 4", got)
	}
	if MobiusMu(big.NewInt(30)) != -1 {
		t.Error("mu(30) want -1")
	}
	if MobiusMu(big.NewInt(12)) != 0 {
		t.Error("mu(12) want 0")
	}
	f := FactorInt(big.NewInt(360)) // 2^3 * 3^2 * 5
	want := []PrimePower{{big.NewInt(2), 3}, {big.NewInt(3), 2}, {big.NewInt(5), 1}}
	if len(f) != 3 {
		t.Fatalf("FactorInt(360) len=%d", len(f))
	}
	for i, pp := range f {
		if !bigEq(pp.Prime, want[i].Prime) || pp.Exp != want[i].Exp {
			t.Errorf("factor %d = %v^%d want %v^%d", i, pp.Prime, pp.Exp, want[i].Prime, want[i].Exp)
		}
	}
	if NumDivisors(big.NewInt(360)) != 24 {
		t.Errorf("NumDivisors(360)=%d want 24", NumDivisors(big.NewInt(360)))
	}
	g, x, y := ExtendedGcd(big.NewInt(240), big.NewInt(46))
	if !bigEq(g, big.NewInt(2)) {
		t.Errorf("gcd=%v want 2", g)
	}
	chk := new(big.Int).Add(new(big.Int).Mul(x, big.NewInt(240)), new(big.Int).Mul(y, big.NewInt(46)))
	if !bigEq(chk, g) {
		t.Errorf("bezout mismatch: %v", chk)
	}
}

// ---- polynomial arithmetic ----

func TestPolyArithmetic(t *testing.T) {
	p := big.NewInt(7)
	a := NewPoly(p, 1, 2, 3) // 3x^2+2x+1
	b := NewPoly(p, 4, 5)    // 5x+4
	sum := a.Add(b)
	if !sum.Equal(NewPoly(p, 5, 0, 3)) {
		t.Errorf("sum=%v", sum)
	}
	prod := a.Mul(b)
	// (3x^2+2x+1)(5x+4) = 15x^3+12x^2+10x^2+8x+5x+4 = 15x^3+22x^2+13x+4
	// mod 7: 15->1, 22->1, 13->6, 4 => x^3 + x^2 + 6x + 4
	if !prod.Equal(NewPoly(p, 4, 6, 1, 1)) {
		t.Errorf("prod=%v", prod)
	}
	q, r, err := a.DivMod(b)
	if err != nil {
		t.Fatal(err)
	}
	// verify a = q*b + r
	if !q.Mul(b).Add(r).Equal(a) {
		t.Errorf("division identity failed: q=%v r=%v", q, r)
	}
	if r.Degree() >= b.Degree() {
		t.Errorf("remainder degree too large: %v", r)
	}
}

func TestPolyGcdAndDerivative(t *testing.T) {
	p := big.NewInt(5)
	// (x+1)^2(x+2) and (x+1)(x+3) share gcd (x+1)
	f := NewPoly(p, 1, 1).Pow(2).Mul(NewPoly(p, 2, 1))
	g := NewPoly(p, 1, 1).Mul(NewPoly(p, 3, 1))
	gcd := f.Gcd(g)
	if !gcd.Equal(NewPoly(p, 1, 1)) {
		t.Errorf("gcd=%v want x+1", gcd)
	}
	// extended gcd cofactor identity
	G, s, tt := f.ExtendedGcd(g)
	lhs := f.Mul(s).Add(g.Mul(tt))
	if !lhs.Equal(G) {
		t.Errorf("extended gcd identity failed: %v != %v", lhs, G)
	}
	// derivative of x^3+x+1 over GF(5) is 3x^2+1
	d := NewPoly(p, 1, 1, 0, 1).Derivative()
	if !d.Equal(NewPoly(p, 1, 0, 3)) {
		t.Errorf("derivative=%v", d)
	}
}

func TestPolyPowMod(t *testing.T) {
	p := big.NewInt(2)
	f := NewPoly(p, 1, 1, 0, 1) // x^3+x+1
	x := PolyX(p)
	// x^8 mod (x^3+x+1) should equal x (since field GF(8), x^(2^3)=x)
	got, err := x.PowMod(big.NewInt(8), f)
	if err != nil {
		t.Fatal(err)
	}
	if !got.Equal(x) {
		t.Errorf("x^8 mod f = %v want x", got)
	}
}

// ---- irreducibility & primitivity ----

func TestIrreducibility(t *testing.T) {
	p := big.NewInt(2)
	irr := NewPoly(p, 1, 1, 0, 1) // x^3+x+1 irreducible
	red := NewPoly(p, 0, 1, 0, 1) // x^3+x = x(x^2+1) reducible
	if !IsIrreducible(irr) {
		t.Error("x^3+x+1 should be irreducible over GF(2)")
	}
	if IsIrreducible(red) {
		t.Error("x^3+x should be reducible over GF(2)")
	}
	// count matches enumeration
	for _, n := range []int{2, 3, 4} {
		want := NumberOfIrreducibles(p, n)
		all := AllIrreduciblePolys(p, n)
		if int64(len(all)) != want.Int64() {
			t.Errorf("GF(2) deg %d: enumerated %d, formula %v", n, len(all), want)
		}
	}
}

func TestFindPolys(t *testing.T) {
	p := big.NewInt(3)
	f, err := FindIrreducible(p, 4)
	if err != nil {
		t.Fatal(err)
	}
	if !IsIrreducible(f) || f.Degree() != 4 {
		t.Errorf("FindIrreducible gave non-irreducible or wrong degree: %v", f)
	}
	prim, err := FindPrimitivePoly(big.NewInt(2), 4)
	if err != nil {
		t.Fatal(err)
	}
	if !IsPrimitivePoly(prim) {
		t.Errorf("FindPrimitivePoly not primitive: %v", prim)
	}
	conway, err := ConwayStylePoly(big.NewInt(2), 3)
	if err != nil {
		t.Fatal(err)
	}
	if !IsPrimitivePoly(conway) {
		t.Errorf("ConwayStylePoly not primitive: %v", conway)
	}
}

// ---- factorization ----

func TestFactorization(t *testing.T) {
	// (x+1)^2 over GF(2) from x^2+1
	f := NewPoly(big.NewInt(2), 1, 0, 1)
	lead, factors := Factor(f)
	if !bigEq(lead, big1) {
		t.Errorf("lead=%v", lead)
	}
	if len(factors) != 1 || factors[0].Exp != 2 || !factors[0].Factor.Equal(NewPoly(big.NewInt(2), 1, 1)) {
		t.Errorf("Factor(x^2+1 over GF2)=%v", factors)
	}
	// x^2+1 over GF(5) = (x+2)(x+3)
	g := NewPoly(big.NewInt(5), 1, 0, 1)
	_, gf := Factor(g)
	if DistinctFactorCount(g) != 2 || FactorCount(g) != 2 {
		t.Errorf("Factor(x^2+1 over GF5) distinct=%d total=%d", DistinctFactorCount(g), FactorCount(g))
	}
	// reconstruct product
	prod := PolyOne(big.NewInt(5))
	for _, pf := range gf {
		prod = prod.Mul(pf.Factor.Pow(pf.Exp))
	}
	if !prod.Equal(g) {
		t.Errorf("reconstruction failed: %v != %v", prod, g)
	}
	// a larger random-ish product over GF(3)
	p := big.NewInt(3)
	base := NewPoly(p, 2, 1).Mul(NewPoly(p, 1, 0, 1)).Mul(NewPoly(p, 1, 1)) // (x+2)(x^2+1)(x+1)
	_, bf := Factor(base)
	rec := PolyOne(p)
	for _, pf := range bf {
		rec = rec.Mul(pf.Factor.Pow(pf.Exp))
	}
	if !rec.Equal(base.Monic()) {
		t.Errorf("GF(3) reconstruction failed: %v != %v", rec, base.Monic())
	}
	if !IsSquareFree(NewPoly(p, 1, 1).Mul(NewPoly(p, 2, 1))) {
		t.Error("(x+1)(x+2) should be square free")
	}
	if IsSquareFree(NewPoly(p, 1, 1).Pow(2)) {
		t.Error("(x+1)^2 should not be square free")
	}
}

// ---- field arithmetic in GF(p^n) ----

func gf8(t *testing.T) *Field {
	t.Helper()
	f, err := NewField(big.NewInt(2), NewPoly(big.NewInt(2), 1, 1, 0, 1))
	if err != nil {
		t.Fatal(err)
	}
	return f
}

func TestFieldBasics(t *testing.T) {
	f := gf8(t)
	if !bigEq(f.Order(), big.NewInt(8)) || f.Degree() != 3 {
		t.Errorf("order=%v degree=%d", f.Order(), f.Degree())
	}
	x := f.Element(0, 1) // element x
	cube, err := x.Pow(big.NewInt(3))
	if err != nil {
		t.Fatal(err)
	}
	if !cube.Equal(f.Element(1, 1)) { // x^3 = x+1
		t.Errorf("x^3=%v want x+1", cube)
	}
	// x has order 7 (primitive)
	if !bigEq(cube.Field.Order(), big.NewInt(8)) {
		t.Fatal("field mismatch")
	}
	seven, _ := x.Pow(big.NewInt(7))
	if !seven.IsOne() {
		t.Errorf("x^7=%v want 1", seven)
	}
	if !x.IsPrimitive() {
		t.Error("x should be primitive in GF(8) with x^3+x+1")
	}
}

func TestFieldInverseAllNonzero(t *testing.T) {
	for _, spec := range []struct {
		p    int64
		mod  []int64
		size int64
	}{
		{2, []int64{1, 1, 0, 1}, 8},     // GF(8)
		{3, []int64{1, 0, 1}, 9},        // GF(9), x^2+1
		{2, []int64{1, 1, 0, 0, 1}, 16}, // GF(16), x^4+x+1
	} {
		f, err := NewField(big.NewInt(spec.p), NewPoly(big.NewInt(spec.p), spec.mod...))
		if err != nil {
			t.Fatalf("field %v: %v", spec, err)
		}
		for _, a := range f.NonzeroElements() {
			inv, err := a.Inv()
			if err != nil {
				t.Fatalf("Inv(%v): %v", a, err)
			}
			if !a.Mul(inv).IsOne() {
				t.Errorf("a*inv != 1 for a=%v in %v", a, f)
			}
		}
	}
}

func TestFrobeniusTraceNorm(t *testing.T) {
	// GF(4): x^2+x+1
	f, err := NewField(big.NewInt(2), NewPoly(big.NewInt(2), 1, 1, 1))
	if err != nil {
		t.Fatal(err)
	}
	x := f.Element(0, 1)
	// Frobenius additive
	a := f.Element(1, 1)
	b := f.Element(0, 1)
	lhs := a.Add(b).Frobenius()
	rhs := a.Frobenius().Add(b.Frobenius())
	if !lhs.Equal(rhs) {
		t.Errorf("Frobenius not additive: %v vs %v", lhs, rhs)
	}
	if !bigEq(x.Trace(), big1) {
		t.Errorf("Tr(x)=%v want 1", x.Trace())
	}
	if !bigEq(x.Norm(), big1) {
		t.Errorf("N(x)=%v want 1", x.Norm())
	}
	// minimal poly of x is x^2+x+1 (the modulus)
	mp := x.MinimalPoly()
	if !mp.Equal(NewPoly(big.NewInt(2), 1, 1, 1)) {
		t.Errorf("MinimalPoly(x)=%v want x^2+x+1", mp)
	}
	// x is a root of its minimal polynomial in the field
	acc := f.Zero()
	for i, c := range mp.Coeff {
		xi, _ := x.Pow(big.NewInt(int64(i)))
		acc = acc.Add(xi.Mul(f.Element(c.Int64())))
	}
	if !acc.IsZero() {
		t.Errorf("x not a root of its minimal polynomial: %v", acc)
	}
}

func TestConjugates(t *testing.T) {
	f := gf8(t)
	x := f.Element(0, 1)
	conj := x.Conjugates()
	if len(conj) != 3 {
		t.Fatalf("expected 3 conjugates, got %d", len(conj))
	}
	// conjugates are x, x^2, x^4
	x2, _ := x.Pow(big.NewInt(2))
	x4, _ := x.Pow(big.NewInt(4))
	if !conj[0].Equal(x) || !conj[1].Equal(x2) || !conj[2].Equal(x4) {
		t.Errorf("unexpected conjugates: %v", conj)
	}
}

func TestPrimitiveElementAndOrders(t *testing.T) {
	f := gf8(t)
	g, err := f.PrimitiveElement()
	if err != nil {
		t.Fatal(err)
	}
	ord, _ := g.Order()
	if !bigEq(ord, big.NewInt(7)) {
		t.Errorf("primitive element order=%v want 7", ord)
	}
	// all nonzero orders divide 7
	for _, a := range f.NonzeroElements() {
		o, _ := a.Order()
		if new(big.Int).Mod(big.NewInt(7), o).Sign() != 0 {
			t.Errorf("order %v does not divide 7 for %v", o, a)
		}
	}
}

func TestDiscreteLog(t *testing.T) {
	// prime field GF(7), base 3 (primitive root)
	f, _ := NewPrimeField(big.NewInt(7))
	base := f.Element(3)
	target := f.Element(5)
	x, err := DiscreteLog(base, target)
	if err != nil {
		t.Fatal(err)
	}
	chk, _ := base.Pow(x)
	if !chk.Equal(target) {
		t.Errorf("3^%v = %v want 5", x, chk)
	}
	// extension field GF(8)
	ef := gf8(t)
	g, _ := ef.PrimitiveElement()
	for k := int64(0); k < 7; k++ {
		h, _ := g.Pow(big.NewInt(k))
		lg, err := DiscreteLog(g, h)
		if err != nil {
			t.Fatalf("DiscreteLog: %v", err)
		}
		back, _ := g.Pow(lg)
		if !back.Equal(h) {
			t.Errorf("discrete log roundtrip failed at k=%d", k)
		}
	}
}

func TestSubfields(t *testing.T) {
	f, err := NewFieldFromDegree(big.NewInt(2), 6)
	if err != nil {
		t.Fatal(err)
	}
	subs := f.Subfields()
	want := []int{1, 2, 3, 6}
	if len(subs) != len(want) {
		t.Fatalf("subfields=%v want %v", subs, want)
	}
	for i := range want {
		if subs[i] != want[i] {
			t.Errorf("subfields=%v want %v", subs, want)
		}
	}
	if IsSubfieldDegree(6, 4) {
		t.Error("GF(2^4) is not a subfield of GF(2^6)")
	}
}

func TestCyclotomicCosets(t *testing.T) {
	// binary cyclotomic cosets mod 15 (q=2, n=4)
	cosets := CyclotomicCosets(big.NewInt(2), 4)
	total := 0
	for _, c := range cosets {
		total += len(c)
	}
	if total != 15 {
		t.Errorf("cosets cover %d elements want 15", total)
	}
	c1 := CyclotomicCoset(big.NewInt(2), 1, 4)
	// {1,2,4,8}
	want := []int64{1, 2, 4, 8}
	if len(c1) != len(want) {
		t.Fatalf("coset of 1 = %v want %v", c1, want)
	}
	for i := range want {
		if c1[i] != want[i] {
			t.Errorf("coset of 1 = %v want %v", c1, want)
		}
	}
}

func TestNewFieldRejectsReducible(t *testing.T) {
	// x^2+1 is reducible over GF(2) since (x+1)^2 = x^2+1
	_, err := NewField(big.NewInt(2), NewPoly(big.NewInt(2), 1, 0, 1))
	if err == nil {
		t.Error("expected error for reducible modulus")
	}
}

func TestMatrixNullSpace(t *testing.T) {
	p := big.NewInt(5)
	m := NewMatModP(p, 2, 3)
	m.Set(0, 0, big.NewInt(1))
	m.Set(0, 1, big.NewInt(2))
	m.Set(0, 2, big.NewInt(3))
	m.Set(1, 0, big.NewInt(2))
	m.Set(1, 1, big.NewInt(4))
	m.Set(1, 2, big.NewInt(6))
	if m.Rank() != 1 {
		t.Errorf("rank=%d want 1", m.Rank())
	}
	ns := m.NullSpace()
	for _, v := range ns {
		// M*v == 0
		for i := 0; i < m.Rows; i++ {
			acc := big.NewInt(0)
			for j := 0; j < m.Cols; j++ {
				acc.Add(acc, new(big.Int).Mul(m.Data[i][j], v[j]))
			}
			acc.Mod(acc, p)
			if acc.Sign() != 0 {
				t.Errorf("null vector fails: %v", v)
			}
		}
	}
	if len(ns) != 2 { // cols - rank = 3 - 1
		t.Errorf("nullspace dim=%d want 2", len(ns))
	}
}

func TestRandomElementDeterministic(t *testing.T) {
	f := gf8(t)
	r := rand.New(rand.NewSource(42))
	seen := map[string]bool{}
	for i := 0; i < 50; i++ {
		e := f.RandomElement(r)
		if e.ToInt().Cmp(f.Order()) >= 0 {
			t.Errorf("random element out of range: %v", e)
		}
		seen[e.String()] = true
	}
	if len(seen) < 2 {
		t.Error("random element generator appears constant")
	}
}

// ---- runnable example ----

func ExampleField() {
	// Build GF(2^3) = GF(2)[x]/(x^3 + x + 1).
	f, _ := NewField(big.NewInt(2), NewPoly(big.NewInt(2), 1, 1, 0, 1))
	fmt.Println(f)

	x := f.Element(0, 1) // the residue class of x
	cube, _ := x.Pow(big.NewInt(3))
	fmt.Println("x^3 =", cube) // x^3 = x + 1

	ord, _ := x.Order()
	fmt.Println("order of x =", ord) // 7, so x is primitive

	fmt.Println("Tr(x) =", x.Trace())
	fmt.Println("N(x)  =", x.Norm())
	// Output:
	// GF(2^3)
	// x^3 = x + 1
	// order of x = 7
	// Tr(x) = 0
	// N(x)  = 1
}

func ExampleDiscreteLog() {
	f, _ := NewPrimeField(big.NewInt(7))
	base := f.Element(3) // 3 is a primitive root mod 7
	x, _ := DiscreteLog(base, f.Element(5))
	fmt.Printf("3^%v = 5 (mod 7)\n", x)
	// Output:
	// 3^5 = 5 (mod 7)
}
