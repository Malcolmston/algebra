package cellular

import (
	"fmt"
	"math"
	"testing"
)

const tol = 1e-9

func approx(a, b float64) bool { return math.Abs(a-b) <= tol }

func TestElementaryApply(t *testing.T) {
	// Rule 110 truth table: neighbourhood p -> output bit of 110 = 0b01101110.
	r := Rule110()
	tests := []struct {
		l, c, rt, want int
	}{
		{1, 1, 1, 0},
		{1, 1, 0, 1},
		{1, 0, 1, 1},
		{1, 0, 0, 0},
		{0, 1, 1, 1},
		{0, 1, 0, 1},
		{0, 0, 1, 1},
		{0, 0, 0, 0},
	}
	for _, tc := range tests {
		if got := r.ApplyLCR(tc.l, tc.c, tc.rt); got != tc.want {
			t.Errorf("rule110(%d%d%d)=%d, want %d", tc.l, tc.c, tc.rt, got, tc.want)
		}
		if got := r.Apply([]int{tc.l, tc.c, tc.rt}); got != tc.want {
			t.Errorf("rule110.Apply(%d%d%d)=%d, want %d", tc.l, tc.c, tc.rt, got, tc.want)
		}
	}
}

func TestElementaryConstruction(t *testing.T) {
	if _, err := NewElementaryRule(256); err == nil {
		t.Error("expected error for rule 256")
	}
	if _, err := NewElementaryRule(-1); err == nil {
		t.Error("expected error for rule -1")
	}
	r := MustElementaryRule(110)
	if r.Number() != 110 || r.States() != 2 || r.Radius() != 1 {
		t.Errorf("rule110 metadata wrong: %d %d %d", r.Number(), r.States(), r.Radius())
	}
	bitsTable := RuleToBits(r)
	back, err := RuleFromBits(bitsTable)
	if err != nil {
		t.Fatal(err)
	}
	if back != r {
		t.Errorf("RuleFromBits round trip: got %d, want %d", back, r)
	}
}

func TestElementaryConjugates(t *testing.T) {
	tests := []struct {
		rule, mirror, comp, mc, classSize int
	}{
		{110, 124, 137, 193, 4},
		{90, 90, 165, 165, 2},
		{30, 86, 135, 149, 4},
		{60, 102, 195, 153, 4},
		{150, 150, 150, 150, 1},
		{184, 226, 226, 184, 2},
	}
	for _, tc := range tests {
		r := ElementaryRule(tc.rule)
		if int(r.MirrorRule()) != tc.mirror {
			t.Errorf("rule %d mirror = %d, want %d", tc.rule, r.MirrorRule(), tc.mirror)
		}
		if int(r.ComplementRule()) != tc.comp {
			t.Errorf("rule %d complement = %d, want %d", tc.rule, r.ComplementRule(), tc.comp)
		}
		if int(r.MirrorComplementRule()) != tc.mc {
			t.Errorf("rule %d mirror-complement = %d, want %d", tc.rule, r.MirrorComplementRule(), tc.mc)
		}
		if r.EquivalenceClassSize() != tc.classSize {
			t.Errorf("rule %d class size = %d, want %d", tc.rule, r.EquivalenceClassSize(), tc.classSize)
		}
	}
	// Double mirror and double complement are the identity.
	for n := 0; n < 256; n++ {
		r := ElementaryRule(n)
		if r.MirrorRule().MirrorRule() != r {
			t.Fatalf("mirror not involutive at %d", n)
		}
		if r.ComplementRule().ComplementRule() != r {
			t.Fatalf("complement not involutive at %d", n)
		}
	}
}

func TestInequivalentCount(t *testing.T) {
	if got := len(InequivalentElementaryRules()); got != 88 {
		t.Errorf("inequivalent rules = %d, want 88", got)
	}
	if got := len(AllElementaryRules()); got != 256 {
		t.Errorf("all rules = %d, want 256", got)
	}
	// Every rule's canonical representative must itself be canonical.
	for n := 0; n < 256; n++ {
		m := ElementaryRule(n).MinimalEquivalent()
		if !m.IsCanonical() {
			t.Errorf("canonical of %d is %d which is not canonical", n, m)
		}
	}
}

func TestAdditiveRules(t *testing.T) {
	want := map[int]bool{0: true, 60: true, 90: true, 102: true, 150: true, 170: true, 204: true, 240: true}
	got := AdditiveElementaryRules()
	if len(got) != len(want) {
		t.Fatalf("additive rules = %v, want the 8 linear rules", got)
	}
	for _, r := range got {
		if !want[int(r)] {
			t.Errorf("rule %d reported additive but should not be", int(r))
		}
	}
	// Rule 90 and 150 are additive; 110 and 30 are not.
	if !Rule90().IsAdditive() || !Rule150().IsAdditive() {
		t.Error("rule 90/150 should be additive")
	}
	if Rule110().IsAdditive() || Rule30().IsAdditive() {
		t.Error("rule 110/30 should not be additive")
	}
	// Rule 255 is affine (constant 1) but not strictly additive.
	if ElementaryRule(255).IsAdditive() {
		t.Error("rule 255 is affine, not additive")
	}
}

func TestRuleProperties(t *testing.T) {
	if !Rule150().IsTotalistic() {
		t.Error("rule 150 should be totalistic")
	}
	if Rule110().IsTotalistic() {
		t.Error("rule 110 should not be totalistic")
	}
	if !Rule90().IsOuterTotalistic() {
		t.Error("rule 90 should be outer-totalistic")
	}
	if !Rule90().IsSymmetric() {
		t.Error("rule 90 should be symmetric")
	}
	if Rule110().IsSymmetric() {
		t.Error("rule 110 should not be symmetric")
	}
	if !ElementaryRule(0).IsQuiescent() || ElementaryRule(1).IsQuiescent() {
		t.Error("quiescence check wrong")
	}
	if !approx(Rule110().LambdaParameter(), 0.625) {
		t.Errorf("rule110 lambda = %v, want 0.625", Rule110().LambdaParameter())
	}
	if Rule110().ActiveCount() != 5 {
		t.Errorf("rule110 active count = %d, want 5", Rule110().ActiveCount())
	}
}

func TestSierpinski(t *testing.T) {
	// Rule 90 from a single seed produces Pascal's triangle mod 2 (Sierpinski).
	rows := ElementaryEvolveSeed(90, 9, 4)
	want := []string{
		"....#....",
		"...#.#...",
		"..#...#..",
		".#.#.#.#.",
		"#.......#",
	}
	got := DiagramString(rows)
	if got != joinLines(want) {
		t.Errorf("rule90 Sierpinski mismatch:\n%s\nwant:\n%s", got, joinLines(want))
	}
	// Light cone half-width grows by exactly one per step.
	widths := LightConeWidth(rows)
	for i, w := range widths {
		if w != i {
			t.Errorf("light cone width at step %d = %d, want %d", i, w, i)
		}
	}
}

func joinLines(lines []string) string {
	out := ""
	for i, l := range lines {
		if i > 0 {
			out += "\n"
		}
		out += l
	}
	return out
}

func TestBoundaryConditions(t *testing.T) {
	state := []int{1, 0, 0}
	if got := cellAt(state, -1, Periodic); got != 0 {
		t.Errorf("periodic -1 = %d, want 0 (last cell)", got)
	}
	if got := cellAt(state, 3, Periodic); got != 1 {
		t.Errorf("periodic 3 = %d, want 1 (first cell)", got)
	}
	if got := cellAt(state, -1, FixedOne); got != 1 {
		t.Errorf("fixed-one -1 = %d, want 1", got)
	}
	if got := cellAt(state, -1, FixedZero); got != 0 {
		t.Errorf("fixed-zero -1 = %d, want 0", got)
	}
	if got := cellAt(state, -1, Reflect); got != 1 {
		t.Errorf("reflect -1 = %d, want 1 (edge duplicated)", got)
	}
	// Rule 90 XOR on a ring conserves parity check: periodic evolution stays
	// within the ring width.
	init := SingleSeedState(7)
	rows := Evolve1D(Rule90(), init, 3, Periodic)
	for _, r := range rows {
		if len(r) != 7 {
			t.Fatalf("periodic evolution changed width to %d", len(r))
		}
	}
}

func TestClassHeuristic(t *testing.T) {
	tests := []struct {
		rule, class int
	}{
		{0, 1},
		{255, 1},
		{4, 2},
		{90, 3},
		{30, 3},
		{150, 3},
	}
	for _, tc := range tests {
		if got := ElementaryRule(tc.rule).Class(); got != tc.class {
			t.Errorf("Class(%d) = %d, want %d", tc.rule, got, tc.class)
		}
	}
	if ClassName(1) != "uniform" || ClassName(4) != "complex" {
		t.Error("ClassName wrong")
	}
}

func TestTotalisticRule(t *testing.T) {
	// A binary radius-1 totalistic rule with code 6 = digits [0,1,1,0]:
	// sum 0->0, 1->1, 2->1, 3->0. This is the "majority-off / xor-ish" rule.
	tr, err := NewTotalisticRule(2, 1, 6)
	if err != nil {
		t.Fatal(err)
	}
	if tr.MaxSum() != 3 {
		t.Errorf("MaxSum = %d, want 3", tr.MaxSum())
	}
	tab := tr.Table()
	wantTab := []int{0, 1, 1, 0}
	for i := range wantTab {
		if tab[i] != wantTab[i] {
			t.Errorf("table[%d] = %d, want %d", i, tab[i], wantTab[i])
		}
	}
	cases := []struct {
		nb   []int
		want int
	}{
		{[]int{0, 0, 0}, 0},
		{[]int{1, 0, 0}, 1},
		{[]int{1, 1, 0}, 1},
		{[]int{1, 1, 1}, 0},
	}
	for _, c := range cases {
		if got := tr.Apply(c.nb); got != c.want {
			t.Errorf("Apply(%v) = %d, want %d", c.nb, got, c.want)
		}
	}
	if _, err := NewTotalisticRule(1, 1, 0); err == nil {
		t.Error("expected error for k<2")
	}
}

func TestOuterTotalistic(t *testing.T) {
	// Reconstruct Conway-like behaviour is 2-D; here just check a 1-D
	// outer-totalistic rule round-trips through its table.
	o, err := NewOuterTotalisticRule(2, 1, 6)
	if err != nil {
		t.Fatal(err)
	}
	if o.OuterMax() != 2 {
		t.Errorf("OuterMax = %d, want 2", o.OuterMax())
	}
	// index = centre*(outerMax+1)+outerSum, table from code 6 = [0,1,1,0,0,0]
	// with 6 entries (k=2, outerMax=2 -> 2*3=6 entries).
	if len(o.Table()) != 6 {
		t.Errorf("table len = %d, want 6", len(o.Table()))
	}
}

func TestGeneralRule(t *testing.T) {
	// A general k=2 r=1 rule built from rule 110's table must evolve identically.
	g := ElementaryAsGeneral(Rule110())
	init := RandomBinaryState(30, 7)
	a := IterateState(Rule110(), init, 10, Periodic)
	b := IterateState(g, init, 10, Periodic)
	if !EqualState(a, b) {
		t.Error("GeneralRule does not match elementary rule 110")
	}
	if !approx(LangtonLambda(g), 0.625) {
		t.Errorf("LangtonLambda = %v, want 0.625", LangtonLambda(g))
	}
	if _, err := NewGeneralRule(2, 1, []int{0, 1}); err == nil {
		t.Error("expected error for wrong table length")
	}
}

func TestEncoding(t *testing.T) {
	code, err := EncodeBaseK([]int{1, 0, 1}, 2)
	if err != nil || code != 5 {
		t.Errorf("EncodeBaseK([1 0 1],2) = %d, %v; want 5", code, err)
	}
	digits, err := DecodeBaseK(5, 3, 2)
	if err != nil {
		t.Fatal(err)
	}
	if digits[0] != 1 || digits[1] != 0 || digits[2] != 1 {
		t.Errorf("DecodeBaseK = %v, want [1 0 1]", digits)
	}
	if _, err := DecodeBaseK(8, 3, 2); err == nil {
		t.Error("expected overflow error for 8 in 3 bits")
	}
	if NumRules(2, 1).Int64() != 256 {
		t.Errorf("NumRules(2,1) = %v, want 256", NumRules(2, 1))
	}
	if NumTotalisticRules(2, 1).Int64() != 16 {
		t.Errorf("NumTotalisticRules(2,1) = %v, want 16", NumTotalisticRules(2, 1))
	}
}

func TestEntropy(t *testing.T) {
	if !approx(BinaryEntropy(0.5), 1.0) {
		t.Errorf("BinaryEntropy(0.5) = %v, want 1", BinaryEntropy(0.5))
	}
	if !approx(BinaryEntropy(0), 0) || !approx(BinaryEntropy(1), 0) {
		t.Error("BinaryEntropy at 0/1 should be 0")
	}
	if !approx(SpatialEntropy([]int{0, 1, 0, 1}), 1.0) {
		t.Errorf("balanced entropy = %v, want 1", SpatialEntropy([]int{0, 1, 0, 1}))
	}
	if !approx(SpatialEntropy([]int{1, 1, 1, 1}), 0) {
		t.Error("uniform entropy should be 0")
	}
	// Uniform distribution over 4 symbols has entropy 2 bits.
	if !approx(ShannonEntropyBits([]float64{0.25, 0.25, 0.25, 0.25}), 2.0) {
		t.Error("4-symbol uniform entropy should be 2")
	}
}

func TestDamageSpread(t *testing.T) {
	// Rule 90 (linear XOR) spreads damage; identity-like rule 204 does not.
	init := RandomBinaryState(101, 42)
	dmg := DamageSpread(Rule90(), init, 50, 20, Periodic)
	if dmg[0] != 1 {
		t.Errorf("initial damage = %d, want 1", dmg[0])
	}
	if dmg[20] <= dmg[0] {
		t.Errorf("rule 90 damage did not spread: %v", dmg)
	}
	// Rule 204 is the identity; damage stays exactly 1.
	stay := DamageSpread(ElementaryRule(204), init, 50, 20, Periodic)
	for i, d := range stay {
		if d != 1 {
			t.Errorf("identity rule damage at %d = %d, want 1", i, d)
		}
	}
}

func TestSecondOrderReversible(t *testing.T) {
	prev := RandomBinaryState(40, 11)
	cur := RandomBinaryState(40, 22)
	for _, ruleNum := range []int{90, 30, 110, 150} {
		so, err := NewSecondOrderCA(ElementaryRule(ruleNum), prev, cur, Periodic)
		if err != nil {
			t.Fatal(err)
		}
		snapPrev := CloneState(so.Prev)
		snapCur := CloneState(so.Cur)
		for i := 0; i < 25; i++ {
			so.Step()
		}
		for i := 0; i < 25; i++ {
			so.StepBack()
		}
		if !EqualState(so.Cur, snapCur) || !EqualState(so.Prev, snapPrev) {
			t.Errorf("rule %d second-order CA not reversible", ruleNum)
		}
	}
	// Length mismatch is rejected.
	if _, err := NewSecondOrderCA(Rule90(), []int{0}, []int{0, 1}, Periodic); err == nil {
		t.Error("expected length-mismatch error")
	}
}

func TestMargolus(t *testing.T) {
	for name, m := range map[string]MargolusRule{
		"critters": CrittersRule(),
		"rot180":   Rotate180Rule(),
	} {
		if !m.IsPermutation() {
			t.Errorf("%s is not a permutation", name)
		}
		inv, err := m.Inverse()
		if err != nil {
			t.Fatal(err)
		}
		// Both example rules are involutions.
		if inv != m {
			t.Errorf("%s inverse differs from itself", name)
		}
	}
	// A non-permutation is detected.
	var bad MargolusRule
	if bad.IsPermutation() {
		t.Error("all-zero Margolus rule should not be a permutation")
	}

	// Margolus evolution followed by inverse evolution restores the grid.
	g, _ := NewGridFrom([][]int{
		{1, 0, 1, 0},
		{0, 1, 0, 0},
		{1, 1, 0, 1},
		{0, 0, 1, 0},
	})
	crit := CrittersRule()
	step, err := MargolusStep(crit, g, true)
	if err != nil {
		t.Fatal(err)
	}
	back, err := MargolusStep(crit, step, true)
	if err != nil {
		t.Fatal(err)
	}
	if !back.Equal(g) {
		t.Error("involutive Margolus step should restore grid")
	}
	// Odd dimensions are rejected.
	odd := NewGrid(3, 4)
	if _, err := MargolusStep(crit, odd, true); err == nil {
		t.Error("expected error for odd dimensions")
	}
}

func TestConwayLife(t *testing.T) {
	l := Conway()
	if l.String() != "B3/S23" {
		t.Errorf("Conway string = %q, want B3/S23", l.String())
	}
	// Blinker oscillates with period 2.
	blinker := CenteredPattern(Blinker(), 5, 5)
	period, pre := FindPeriod(l, blinker, Moore, FixedZero, 20)
	if period != 2 || pre != 0 {
		t.Errorf("blinker period = %d pre = %d, want 2, 0", period, pre)
	}
	// Block is a still life.
	block := CenteredPattern(Block(), 4, 4)
	bp, _ := FindPeriod(l, block, Moore, FixedZero, 5)
	if bp != 1 {
		t.Errorf("block period = %d, want 1", bp)
	}
	// A glider conserves population (5 cells) each generation on a torus.
	g := GliderOn(12, 12)
	frames := LifeEvolve(l, g, 16, Moore, Periodic)
	for i, f := range frames {
		if f.Population() != 5 {
			t.Errorf("glider population at step %d = %d, want 5", i, f.Population())
		}
	}
	// A glider returns to its own shape after 4 steps, translated by (1,1).
	after4 := frames[4]
	orig := LiveCells(g)
	shifted := LiveCells(after4)
	if len(orig) != len(shifted) {
		t.Fatal("glider changed cell count")
	}
	for i := range orig {
		if shifted[i][0]-orig[i][0] != 1 || shifted[i][1]-orig[i][1] != 1 {
			t.Errorf("glider not translated by (1,1) at cell %d: %v -> %v", i, orig[i], shifted[i])
			break
		}
	}
}

func TestRuleStringParsing(t *testing.T) {
	tests := []struct {
		in, canonical string
	}{
		{"B3/S23", "B3/S23"},
		{"b36/s23", "B36/S23"},
		{"23/3", "B3/S23"},
		{"B2/S", "B2/S"},
	}
	for _, tc := range tests {
		l, err := ParseRuleString(tc.in)
		if err != nil {
			t.Fatalf("ParseRuleString(%q): %v", tc.in, err)
		}
		if l.String() != tc.canonical {
			t.Errorf("ParseRuleString(%q).String() = %q, want %q", tc.in, l.String(), tc.canonical)
		}
	}
	for _, bad := range []string{"", "B3", "B9/S2", "X3/S2", "3/2/1"} {
		if _, err := ParseRuleString(bad); err == nil {
			t.Errorf("ParseRuleString(%q) should error", bad)
		}
	}
}

func TestNeighbourCounts(t *testing.T) {
	g, _ := NewGridFrom([][]int{
		{1, 1, 1},
		{1, 0, 1},
		{1, 1, 1},
	})
	if got := CountNeighbours(g, 1, 1, Moore, FixedZero); got != 8 {
		t.Errorf("Moore neighbours of centre = %d, want 8", got)
	}
	if got := CountNeighbours(g, 1, 1, VonNeumann, FixedZero); got != 4 {
		t.Errorf("von Neumann neighbours of centre = %d, want 4", got)
	}
	// Toroidal corner sees wrapped neighbours.
	corner, _ := NewGridFrom([][]int{
		{1, 0},
		{0, 1},
	})
	if got := CountNeighbours(corner, 0, 0, Moore, Periodic); got == 0 {
		t.Error("toroidal Moore count should wrap around")
	}
}

func TestGridBasics(t *testing.T) {
	g := NewGrid(3, 4)
	g.Set(1, 2, 1)
	if g.At(1, 2) != 1 || g.Population() != 1 {
		t.Error("Set/At/Population inconsistent")
	}
	if g.At(-1, 0) != 0 {
		t.Error("out-of-range At should be 0")
	}
	c := g.Clone()
	c.Set(0, 0, 1)
	if g.Population() != 1 {
		t.Error("Clone should be independent")
	}
	minR, minC, maxR, maxC, ok := BoundingBox(g)
	if !ok || minR != 1 || minC != 2 || maxR != 1 || maxC != 2 {
		t.Errorf("BoundingBox = %d,%d,%d,%d ok=%v", minR, minC, maxR, maxC, ok)
	}
}

func TestStateHelpers(t *testing.T) {
	s := SingleSeedState(5)
	if Population(s) != 1 || s[2] != 1 {
		t.Error("SingleSeedState wrong")
	}
	if !approx(Density(s), 0.2) {
		t.Errorf("Density = %v, want 0.2", Density(s))
	}
	r1 := RandomBinaryState(50, 99)
	r2 := RandomBinaryState(50, 99)
	if !EqualState(r1, r2) {
		t.Error("RandomBinaryState must be reproducible for equal seeds")
	}
	if EqualState(r1, RandomBinaryState(50, 100)) {
		t.Error("different seeds should generally differ")
	}
	parsed := StateFromString("#.#.", '#')
	if !EqualState(parsed, []int{1, 0, 1, 0}) {
		t.Errorf("StateFromString = %v", parsed)
	}
	if StateToString([]int{1, 0, 1}, '.', '#') != "#.#" {
		t.Error("StateToString wrong")
	}
	hd, err := HammingDistance([]int{0, 1, 0}, []int{1, 1, 0})
	if err != nil || hd != 1 {
		t.Errorf("HammingDistance = %d, %v; want 1", hd, err)
	}
	if _, err := HammingDistance([]int{0}, []int{0, 1}); err == nil {
		t.Error("HammingDistance length mismatch should error")
	}
}

func TestSpacetimeUtilities(t *testing.T) {
	rows := ElementaryEvolveSeed(90, 7, 3)
	h, w := DiagramDimensions(rows)
	if h != 4 || w != 7 {
		t.Errorf("dimensions = %d,%d, want 4,7", h, w)
	}
	tr := TransposeRows(rows)
	if len(tr) != 7 || len(tr[0]) != 4 {
		t.Errorf("transpose dimensions wrong: %d x %d", len(tr), len(tr[0]))
	}
	round := ParseDiagram(DiagramString(rows), '#')
	for i := range rows {
		if !EqualState(rows[i], round[i]) {
			t.Errorf("ParseDiagram round trip failed at row %d", i)
		}
	}
	if CountOnes(rows) != Population(rows[0])+Population(rows[1])+Population(rows[2])+Population(rows[3]) {
		t.Error("CountOnes mismatch")
	}
}

func ExampleElementaryRule_evolve() {
	// Rule 90 grows the Sierpinski triangle from a single central cell.
	rows := ElementaryEvolveSeed(Rule90(), 9, 4)
	fmt.Println(DiagramString(rows))
	// Output:
	// ....#....
	// ...#.#...
	// ..#...#..
	// .#.#.#.#.
	// #.......#
}

func ExampleConway() {
	// One generation turns a horizontal blinker into a vertical one.
	blinker, _ := NewGridFrom([][]int{
		{0, 0, 0},
		{1, 1, 1},
		{0, 0, 0},
	})
	next := LifeStep(Conway(), blinker, Moore, FixedZero)
	fmt.Println(next)
	// Output:
	// .#.
	// .#.
	// .#.
}
