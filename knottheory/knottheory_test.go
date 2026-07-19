package knottheory

import (
	"fmt"
	"math"
	"testing"
)

// --- Laurent polynomial arithmetic ---

func TestLaurentArithmetic(t *testing.T) {
	a := NewLaurent(-1, []int{1, 0, 2}) // X^-1 + 2X
	b := NewLaurent(0, []int{3, -1})    // 3 - X
	if got := a.Add(b); !got.Equal(NewLaurent(-1, []int{1, 3, 1})) {
		t.Errorf("Add: got %s", got)
	}
	if got := a.Sub(a); !got.IsZero() {
		t.Errorf("Sub self: got %s want 0", got)
	}
	// (X^-1 + 2X)(3 - X) = 3X^-1 - 1 + 6X - 2X^2
	want := NewLaurent(-1, []int{3, -1, 6, -2})
	if got := a.Mul(b); !got.Equal(want) {
		t.Errorf("Mul: got %s want %s", got, want)
	}
	if got := a.MinDegree(); got != -1 {
		t.Errorf("MinDegree: got %d", got)
	}
	if got := a.MaxDegree(); got != 1 {
		t.Errorf("MaxDegree: got %d", got)
	}
	if got := a.SpanWidth(); got != 2 {
		t.Errorf("SpanWidth: got %d", got)
	}
	if got := a.NumTerms(); got != 2 {
		t.Errorf("NumTerms: got %d", got)
	}
	m := Monomial(5, -3)
	if got := m.Mul(Monomial(2, 3)); !got.Equal(LaurentConst(10)) {
		t.Errorf("monomial mul: got %s", got)
	}
	if !OneLaurent().IsOne() || ZeroLaurent().IsOne() {
		t.Errorf("IsOne classification wrong")
	}
}

func TestLaurentPowEvalReverse(t *testing.T) {
	p := NewLaurent(0, []int{1, 1}) // 1 + X
	if got := p.Pow(3); !got.Equal(NewLaurent(0, []int{1, 3, 3, 1})) {
		t.Errorf("Pow: got %s", got)
	}
	if got := p.Eval(2); math.Abs(got-3) > 1e-12 {
		t.Errorf("Eval: got %v want 3", got)
	}
	q := NewLaurent(-2, []int{1, 0, 5}) // X^-2 + 5
	if got := q.Reverse(); !got.Equal(NewLaurent(0, []int{5, 0, 1})) {
		t.Errorf("Reverse: got %s", got)
	}
	sym := NewLaurent(-1, []int{1, -1, 1})
	if !sym.IsPalindromic() {
		t.Errorf("expected palindromic")
	}
	if NewLaurent(0, []int{1, 2}).IsPalindromic() {
		t.Errorf("1+2X should not be palindromic")
	}
	if got := q.SubstitutePow(2); !got.Equal(NewLaurent(-4, []int{1, 0, 0, 0, 5})) {
		t.Errorf("SubstitutePow: got %s", got)
	}
	if got := sym.EvalUnit(1); got != 1 {
		t.Errorf("EvalUnit(1): got %d", got)
	}
	if got := sym.EvalUnit(-1); got != -3 {
		t.Errorf("EvalUnit(-1): got %d", got)
	}
}

func TestLaurentDivExact(t *testing.T) {
	// (t^6 - 1) / (t^3 - 1) = t^3 + 1
	num := Monomial(1, 6).Sub(OneLaurent())
	den := Monomial(1, 3).Sub(OneLaurent())
	q, ok := num.DivExact(den)
	if !ok || !q.Equal(NewLaurent(0, []int{1, 0, 0, 1})) {
		t.Errorf("DivExact: ok=%v got %s", ok, q)
	}
	if _, ok := Monomial(1, 2).DivExact(NewLaurent(0, []int{1, 1, 1})); ok {
		t.Errorf("expected inexact division to report ok=false")
	}
}

// --- Permutations ---

func TestPermutation(t *testing.T) {
	p, err := NewPermutation([]int{1, 2, 0})
	if err != nil {
		t.Fatal(err)
	}
	if p.Order() != 3 {
		t.Errorf("Order: got %d want 3", p.Order())
	}
	if p.NumCycles() != 1 {
		t.Errorf("NumCycles: got %d want 1", p.NumCycles())
	}
	if p.Sign() != 1 {
		t.Errorf("Sign of 3-cycle should be +1, got %d", p.Sign())
	}
	if !p.Compose(p.Inverse()).IsIdentity() {
		t.Errorf("p * p^-1 should be identity")
	}
	if _, err := NewPermutation([]int{0, 0}); err == nil {
		t.Errorf("expected error for non-bijection")
	}
	// Transposition decomposition recomposes to p.
	rec := IdentityPermutation(p.Size())
	for _, tr := range p.Transpositions() {
		swap, _ := TranspositionPermutation(p.Size(), tr[0], tr[1])
		rec = rec.Compose(swap)
	}
	if !rec.Equal(p) {
		t.Errorf("transposition recomposition: got %v want %v", rec, p)
	}
	if p.TranspositionCount() != 2 {
		t.Errorf("TranspositionCount: got %d want 2", p.TranspositionCount())
	}
}

// --- Braids ---

func TestBraidBasics(t *testing.T) {
	b := MustBraid(3, 1, 2, -1)
	if b.Strands() != 3 || b.Length() != 3 {
		t.Errorf("strands/length wrong")
	}
	if b.ExponentSum() != 1 {
		t.Errorf("ExponentSum: got %d want 1", b.ExponentSum())
	}
	if bb, _ := b.Concat(b.Inverse()); !bb.FreeReduce().IsTrivial() {
		t.Errorf("b * b^-1 should freely reduce to identity")
	}
	if _, err := NewBraid(3, 3); err == nil {
		t.Errorf("generator 3 is out of range for B_3")
	}
	if BraidGroupRank(5) != 4 {
		t.Errorf("BraidGroupRank(5) should be 4")
	}
	// Artin far commutation preserves the permutation.
	x := MustBraid(4, 1, 3)
	y := MustBraid(4, 3, 1)
	if !x.Permutation().Equal(y.Permutation()) {
		t.Errorf("sigma_1 sigma_3 and sigma_3 sigma_1 should share a permutation")
	}
}

func TestBraidClosureComponents(t *testing.T) {
	tests := []struct {
		name     string
		braid    Braid
		strands  int
		comps    int
		expSum   int
		positive bool
	}{
		{"trefoil", MustBraid(2, 1, 1, 1), 2, 1, 3, true},
		{"hopf", MustBraid(2, 1, 1), 2, 2, 2, true},
		{"unknot2", MustBraid(2, 1), 2, 1, 1, true},
		{"fig8", MustBraid(3, 1, -2, 1, -2), 3, 1, 0, false},
		{"pure", MustBraid(3, 1, 1), 3, 3, 2, true},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			if tc.braid.NumComponents() != tc.comps {
				t.Errorf("components: got %d want %d", tc.braid.NumComponents(), tc.comps)
			}
			if tc.braid.ExponentSum() != tc.expSum {
				t.Errorf("exponent sum: got %d want %d", tc.braid.ExponentSum(), tc.expSum)
			}
			if tc.braid.IsPositive() != tc.positive {
				t.Errorf("positivity: got %v want %v", tc.braid.IsPositive(), tc.positive)
			}
		})
	}
}

func TestTorusBraidAndTwists(t *testing.T) {
	// T(2,3) braid is sigma_1^3.
	b, err := TorusBraid(2, 3)
	if err != nil {
		t.Fatal(err)
	}
	if got := b.Word(); len(got) != 3 || got[0] != 1 {
		t.Errorf("TorusBraid(2,3) word: got %v", got)
	}
	// Full twist on 3 strands has exponent sum 3*(3-1)=6 and is a pure braid.
	ft := FullTwist(3)
	if ft.ExponentSum() != 6 {
		t.Errorf("full twist exponent sum: got %d want 6", ft.ExponentSum())
	}
	if !ft.IsPureBraid() {
		t.Errorf("full twist should be a pure braid")
	}
	// Garside half twist on 4 strands has length 4*3/2 = 6.
	if GarsideHalfTwist(4).Length() != 6 {
		t.Errorf("garside half twist length: got %d want 6", GarsideHalfTwist(4).Length())
	}
	// Its square is the full twist permutation.
	half := GarsideHalfTwist(4)
	sq, _ := half.Concat(half)
	if !sq.Permutation().IsIdentity() {
		t.Errorf("square of half twist should be a pure braid")
	}
}

// --- Burau-based Alexander polynomial of a braid closure ---

func TestBraidAlexander(t *testing.T) {
	tests := []struct {
		name  string
		braid Braid
		want  Laurent
	}{
		{"trefoil", MustBraid(2, 1, 1, 1), NewLaurent(-1, []int{1, -1, 1})},
		{"cinquefoil", MustBraid(2, 1, 1, 1, 1, 1), NewLaurent(-2, []int{1, -1, 1, -1, 1})},
		{"figure8", MustBraid(3, 1, -2, 1, -2), NewLaurent(-1, []int{-1, 3, -1})},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			got, err := tc.braid.AlexanderPolynomial()
			if err != nil {
				t.Fatal(err)
			}
			if !got.Equal(tc.want) {
				t.Errorf("Alexander: got %s want %s", got.StringVar("t"), tc.want.StringVar("t"))
			}
			if !got.IsPalindromic() {
				t.Errorf("Alexander polynomial should be palindromic: %s", got.StringVar("t"))
			}
		})
	}
	// Burau of the identity braid is the identity matrix.
	if got := IdentityBraid(3).ReducedBurau(); got.Rows() != 2 || !got.At(0, 0).IsOne() {
		t.Errorf("reduced Burau of identity should be identity")
	}
}

// --- PD codes, Kauffman bracket and Jones polynomial ---

func TestKauffmanAndJones(t *testing.T) {
	tests := []struct {
		name      string
		pd        PDCode
		writhe    int
		bracket   Laurent
		jonesSqrt Laurent
	}{
		{
			name:      "unknot",
			pd:        UnknotPD(),
			writhe:    0,
			bracket:   OneLaurent(),
			jonesSqrt: OneLaurent(),
		},
		{
			name:      "trefoil",
			pd:        TrefoilPD(),
			writhe:    3,
			bracket:   NewLaurent(-7, []int{1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, -1}), // A^-7 - A^-3 - A^5
			jonesSqrt: NewLaurent(2, []int{1, 0, 0, 0, 1, 0, -1}),                     // u^2 + u^6 - u^8
		},
		{
			name:      "figure8",
			pd:        FigureEightPD(),
			writhe:    0,
			jonesSqrt: NewLaurent(-4, []int{1, 0, -1, 0, 1, 0, -1, 0, 1}),
		},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			if tc.pd.Writhe() != tc.writhe {
				t.Errorf("writhe: got %d want %d", tc.pd.Writhe(), tc.writhe)
			}
			if !tc.bracket.IsZero() {
				if got := tc.pd.KauffmanBracket(); !got.Equal(tc.bracket) {
					t.Errorf("bracket: got %s want %s", got.StringVar("A"), tc.bracket.StringVar("A"))
				}
			}
			if got := tc.pd.JonesPolynomialSqrt(); !got.Equal(tc.jonesSqrt) {
				t.Errorf("jones: got %s want %s", got.StringVar("u"), tc.jonesSqrt.StringVar("u"))
			}
		})
	}
}

func TestJonesInT(t *testing.T) {
	// Figure-eight Jones in t: t^-2 - t^-1 + 1 - t + t^2.
	j, err := FigureEightPD().JonesPolynomial()
	if err != nil {
		t.Fatal(err)
	}
	want := NewLaurent(-2, []int{1, -1, 1, -1, 1})
	if !j.Equal(want) {
		t.Errorf("figure8 Jones(t): got %s want %s", j.StringVar("t"), want.StringVar("t"))
	}
	// Trefoil Jones in t: -t^4 + t^3 + t.
	jt, err := TrefoilPD().JonesPolynomial()
	if err != nil {
		t.Fatal(err)
	}
	if !jt.Equal(NewLaurent(1, []int{1, 0, 1, -1})) {
		t.Errorf("trefoil Jones(t): got %s", jt.StringVar("t"))
	}
	// A 2-component link has half-integer powers, reported as an error.
	if _, err := HopfLinkPD(true).JonesPolynomial(); err == nil {
		t.Errorf("Hopf link Jones(t) should report half-integer powers")
	}
	// Unknot Jones is 1.
	if u, _ := UnknotPD().JonesPolynomial(); !u.IsOne() {
		t.Errorf("unknot Jones should be 1")
	}
	// Mirror sends Jones(t) to Jones(1/t).
	if !jt.Reverse().Equal(mustJones(t, LeftTrefoilPD())) {
		t.Errorf("mirror Jones mismatch")
	}
}

func mustJones(t *testing.T, pd PDCode) Laurent {
	t.Helper()
	j, err := pd.JonesPolynomial()
	if err != nil {
		t.Fatal(err)
	}
	return j
}

// --- Alexander polynomial and determinant from a diagram ---

func TestPDAlexanderAndDeterminant(t *testing.T) {
	tests := []struct {
		name string
		pd   PDCode
		alex Laurent
		det  int
		arf  int
	}{
		{"trefoil", TrefoilPD(), NewLaurent(-1, []int{1, -1, 1}), 3, 1},
		{"figure8", FigureEightPD(), NewLaurent(-1, []int{-1, 3, -1}), 5, 1},
		{"cinquefoil", CinquefoilPD(), NewLaurent(-2, []int{1, -1, 1, -1, 1}), 5, 1},
		{"unknot", UnknotPD(), OneLaurent(), 1, 0},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			got, err := tc.pd.AlexanderPolynomial()
			if err != nil {
				t.Fatal(err)
			}
			if !got.Equal(tc.alex) {
				t.Errorf("Alexander: got %s want %s", got.StringVar("t"), tc.alex.StringVar("t"))
			}
			if tc.pd.KnotDeterminant() != tc.det {
				t.Errorf("determinant: got %d want %d", tc.pd.KnotDeterminant(), tc.det)
			}
			if tc.pd.ArfInvariant() != tc.arf {
				t.Errorf("Arf: got %d want %d", tc.pd.ArfInvariant(), tc.arf)
			}
		})
	}
	// The diagram Alexander polynomial agrees with the braid Burau one.
	braidAlex, _ := MustBraid(2, 1, 1, 1).AlexanderPolynomial()
	pdAlex, _ := TrefoilPD().AlexanderPolynomial()
	if !braidAlex.Equal(pdAlex) {
		t.Errorf("braid and PD Alexander disagree: %s vs %s", braidAlex.StringVar("t"), pdAlex.StringVar("t"))
	}
	// 3-colourability of the trefoil (determinant divisible by 3).
	if !TrefoilPD().IsThreeColorable() {
		t.Errorf("trefoil should be 3-colourable")
	}
	if FigureEightPD().IsThreeColorable() {
		t.Errorf("figure-eight should not be 3-colourable")
	}
	if !FigureEightPD().IsPColorable(5) {
		t.Errorf("figure-eight should be 5-colourable")
	}
}

// --- Gauss codes, writhe and linking numbers ---

func TestGaussCodeAndLinking(t *testing.T) {
	// Trefoil signed Gauss code: O1+ U2+ O3+ U1+ O2+ U3+.
	entries := []GaussEntry{
		{Crossing: 1, Over: true, Sign: 1},
		{Crossing: 2, Over: false, Sign: 1},
		{Crossing: 3, Over: true, Sign: 1},
		{Crossing: 1, Over: false, Sign: 1},
		{Crossing: 2, Over: true, Sign: 1},
		{Crossing: 3, Over: false, Sign: 1},
	}
	gc, err := NewGaussCode(entries)
	if err != nil {
		t.Fatal(err)
	}
	if gc.CrossingNumber() != 3 {
		t.Errorf("crossing number: got %d want 3", gc.CrossingNumber())
	}
	if gc.Writhe() != 3 {
		t.Errorf("writhe: got %d want 3", gc.Writhe())
	}
	if gc.Mirror().Writhe() != -3 {
		t.Errorf("mirror writhe: got %d want -3", gc.Mirror().Writhe())
	}
	// Hopf link: two components sharing two crossings of the same sign.
	compA := GaussCode{Entries: []GaussEntry{
		{Crossing: 1, Over: true, Sign: 1},
		{Crossing: 2, Over: true, Sign: 1},
	}}
	compB := GaussCode{Entries: []GaussEntry{
		{Crossing: 1, Over: false, Sign: 1},
		{Crossing: 2, Over: false, Sign: 1},
	}}
	d, err := NewDiagram(compA, compB)
	if err != nil {
		t.Fatal(err)
	}
	if d.NumComponents() != 2 {
		t.Errorf("components: got %d want 2", d.NumComponents())
	}
	lk, err := d.LinkingNumber(0, 1)
	if err != nil {
		t.Fatal(err)
	}
	if lk != 1 {
		t.Errorf("linking number: got %d want 1", lk)
	}
	if d.Writhe() != 2 {
		t.Errorf("diagram writhe: got %d want 2", d.Writhe())
	}
	// Reversing every sign negates the linking number.
	nA := compA.Mirror()
	nB := compB.Mirror()
	dn, _ := NewDiagram(nA, nB)
	if lk2, _ := dn.LinkingNumber(0, 1); lk2 != -1 {
		t.Errorf("mirror linking number: got %d want -1", lk2)
	}
	// Bad code: crossing appearing twice as over.
	if _, err := NewGaussCode([]GaussEntry{
		{Crossing: 1, Over: true, Sign: 1},
		{Crossing: 1, Over: true, Sign: 1},
	}); err == nil {
		t.Errorf("expected validation error for two over-passes")
	}
}

// --- Reidemeister recognisers ---

func TestReidemeister(t *testing.T) {
	// A single kink: crossing 1 appears at adjacent positions.
	kink := GaussCode{Entries: []GaussEntry{
		{Crossing: 1, Over: true, Sign: 1},
		{Crossing: 1, Over: false, Sign: 1},
		{Crossing: 2, Over: true, Sign: -1},
		{Crossing: 3, Over: false, Sign: 1},
		{Crossing: 2, Over: false, Sign: -1},
		{Crossing: 3, Over: true, Sign: 1},
	}}
	kinks := kink.ReidemeisterIKinks()
	if len(kinks) != 1 || kinks[0] != 1 {
		t.Errorf("ReidemeisterIKinks: got %v want [1]", kinks)
	}
	if !kink.IsReidemeisterIReducible() {
		t.Errorf("expected type I reducible")
	}
	reduced, ok := kink.RemoveReidemeisterIKink(1)
	if !ok || reduced.CrossingNumber() != 2 {
		t.Errorf("RemoveReidemeisterIKink: ok=%v crossings=%d", ok, reduced.CrossingNumber())
	}
	// A reducible type II bigon: crossings 1,2 form an over-arc and an under-arc.
	bigon := GaussCode{Entries: []GaussEntry{
		{Crossing: 1, Over: true, Sign: 1},
		{Crossing: 2, Over: true, Sign: -1},
		{Crossing: 3, Over: true, Sign: 1},
		{Crossing: 2, Over: false, Sign: -1},
		{Crossing: 1, Over: false, Sign: 1},
		{Crossing: 3, Over: false, Sign: 1},
	}}
	pairs := bigon.ReidemeisterIIReducible()
	found := false
	for _, p := range pairs {
		if p == [2]int{1, 2} {
			found = true
		}
	}
	if !found {
		t.Errorf("ReidemeisterIIReducible: got %v want to contain [1 2]", pairs)
	}
	if red2, ok := bigon.RemoveReidemeisterII(1, 2); !ok || red2.CrossingNumber() != 1 {
		t.Errorf("RemoveReidemeisterII: ok=%v crossings=%d", ok, red2.CrossingNumber())
	}
	// The trefoil has no type I or type II reductions.
	tref := GaussCode{Entries: []GaussEntry{
		{Crossing: 1, Over: true, Sign: 1},
		{Crossing: 2, Over: false, Sign: 1},
		{Crossing: 3, Over: true, Sign: 1},
		{Crossing: 1, Over: false, Sign: 1},
		{Crossing: 2, Over: true, Sign: 1},
		{Crossing: 3, Over: false, Sign: 1},
	}}
	if !tref.IsReduced() {
		t.Errorf("trefoil diagram should be reduced")
	}
	if len(tref.ReidemeisterIIITriangles()) == 0 {
		t.Errorf("trefoil should present a type III triangle of chords")
	}
}

// --- Torus knot invariants ---

func TestTorusInvariants(t *testing.T) {
	tests := []struct {
		p, q      int
		crossing  int
		genus     int
		bridge    int
		det       int
		component int
	}{
		{2, 3, 3, 1, 2, 3, 1},
		{2, 5, 5, 2, 2, 5, 1},
		{2, 7, 7, 3, 2, 7, 1},
		{3, 4, 8, 3, 3, 3, 1},
		{3, 5, 10, 4, 3, 1, 1},
		{2, 4, 0, 0, 0, 0, 2}, // a link, not a knot
	}
	for _, tc := range tests {
		name := fmt.Sprintf("T(%d,%d)", tc.p, tc.q)
		t.Run(name, func(t *testing.T) {
			if TorusLinkComponents(tc.p, tc.q) != tc.component {
				t.Errorf("components: got %d want %d", TorusLinkComponents(tc.p, tc.q), tc.component)
			}
			if tc.component != 1 {
				return // remaining checks apply to knots
			}
			if TorusKnotCrossingNumber(tc.p, tc.q) != tc.crossing {
				t.Errorf("crossing: got %d want %d", TorusKnotCrossingNumber(tc.p, tc.q), tc.crossing)
			}
			if TorusKnotGenus(tc.p, tc.q) != tc.genus {
				t.Errorf("genus: got %d want %d", TorusKnotGenus(tc.p, tc.q), tc.genus)
			}
			if TorusKnotBridgeNumber(tc.p, tc.q) != tc.bridge {
				t.Errorf("bridge: got %d want %d", TorusKnotBridgeNumber(tc.p, tc.q), tc.bridge)
			}
			if TorusKnotDeterminant(tc.p, tc.q) != tc.det {
				t.Errorf("determinant: got %d want %d", TorusKnotDeterminant(tc.p, tc.q), tc.det)
			}
			if TorusKnotUnknottingNumber(tc.p, tc.q) != tc.genus {
				t.Errorf("unknotting number should equal genus")
			}
		})
	}
}

func TestTorusPolynomialsMatchDiagrams(t *testing.T) {
	// Torus formula for T(2,3) matches the trefoil diagram.
	tj, err := TorusKnotJones(2, 3)
	if err != nil {
		t.Fatal(err)
	}
	dj, err := TrefoilPD().JonesPolynomial()
	if err != nil {
		t.Fatal(err)
	}
	if !tj.Equal(dj) {
		t.Errorf("T(2,3) Jones: formula %s vs diagram %s", tj.StringVar("t"), dj.StringVar("t"))
	}
	ta, _ := TorusKnotAlexander(2, 3)
	da, _ := TrefoilPD().AlexanderPolynomial()
	if !ta.Equal(da) {
		t.Errorf("T(2,3) Alexander: formula %s vs diagram %s", ta.StringVar("t"), da.StringVar("t"))
	}
	// T(2,5) matches the cinquefoil.
	tj5, _ := TorusKnotJones(2, 5)
	dj5, _ := CinquefoilPD().JonesPolynomial()
	if !tj5.Equal(dj5) {
		t.Errorf("T(2,5) Jones: formula %s vs diagram %s", tj5.StringVar("t"), dj5.StringVar("t"))
	}
	// Torus Alexander is always palindromic.
	for _, pq := range [][2]int{{3, 4}, {3, 5}, {2, 9}} {
		a, err := TorusKnotAlexander(pq[0], pq[1])
		if err != nil {
			t.Fatal(err)
		}
		if !a.IsPalindromic() {
			t.Errorf("T(%d,%d) Alexander not palindromic: %s", pq[0], pq[1], a.StringVar("t"))
		}
	}
	// Coprimality is required.
	if _, err := TorusKnotJones(2, 4); err == nil {
		t.Errorf("TorusKnotJones(2,4) should error (not coprime)")
	}
}

// --- Matrix determinant over the Laurent ring ---

func TestLaurentMatrixDeterminant(t *testing.T) {
	m := NewLaurentMatrix(2, 2)
	m.Set(0, 0, Monomial(1, 1)) // t
	m.Set(0, 1, OneLaurent())
	m.Set(1, 0, OneLaurent())
	m.Set(1, 1, Monomial(1, -1)) // t^-1
	// det = t*t^-1 - 1 = 0
	if got := m.Determinant(); !got.IsZero() {
		t.Errorf("determinant: got %s want 0", got)
	}
	id := IdentityLaurentMatrix(3)
	if got := id.Determinant(); !got.IsOne() {
		t.Errorf("identity determinant should be 1, got %s", got)
	}
}

// Example_trefoilInvariants computes several invariants of the right-handed
// trefoil and prints them.
func Example_trefoilInvariants() {
	pd := TrefoilPD()
	jones, _ := pd.JonesPolynomial()
	alex, _ := pd.AlexanderPolynomial()
	fmt.Println("writhe:", pd.Writhe())
	fmt.Println("jones:", jones.StringVar("t"))
	fmt.Println("alexander:", alex.StringVar("t"))
	fmt.Println("determinant:", pd.KnotDeterminant())
	fmt.Println("3-colorable:", pd.IsThreeColorable())
	// Output:
	// writhe: 3
	// jones: -t^4 + t^3 + t
	// alexander: t - 1 + t^-1
	// determinant: 3
	// 3-colorable: true
}

// ExampleBraid_AlexanderPolynomial recovers the Alexander polynomial of the
// figure-eight knot from a braid whose closure is that knot.
func ExampleBraid_AlexanderPolynomial() {
	b := MustBraid(3, 1, -2, 1, -2)
	alex, _ := b.AlexanderPolynomial()
	fmt.Println("components:", b.NumComponents())
	fmt.Println("alexander:", alex.StringVar("t"))
	// Output:
	// components: 1
	// alexander: -t + 3 - t^-1
}
