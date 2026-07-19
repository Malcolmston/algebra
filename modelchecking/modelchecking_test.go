package modelchecking

import (
	"fmt"
	"math/rand"
	"reflect"
	"testing"
)

// ---------- helpers ----------

func eqInts(a, b []int) bool {
	if len(a) == 0 && len(b) == 0 {
		return true
	}
	return reflect.DeepEqual(a, b)
}

// sampleModel builds a fixed 4-state Kripke structure used across tests.
//
//	s0 {}    -> s1, s2
//	s1 {p}   -> s3
//	s2 {}    -> s2
//	s3 {p,q} -> s3
func sampleModel() *Kripke {
	k := NewKripke(4)
	k.SetInitial(0)
	k.AddTransition(0, 1)
	k.AddTransition(0, 2)
	k.AddTransition(1, 3)
	k.AddTransition(2, 2)
	k.AddTransition(3, 3)
	k.AddLabel(1, "p")
	k.AddLabel(3, "p")
	k.AddLabel(3, "q")
	return k
}

// ---------- StateSet ----------

func TestStateSet(t *testing.T) {
	s := NewStateSet(10)
	for _, e := range []int{1, 3, 5, 9} {
		s.Add(e)
	}
	if s.Len() != 4 {
		t.Fatalf("len=%d", s.Len())
	}
	if !eqInts(s.Elements(), []int{1, 3, 5, 9}) {
		t.Fatalf("elements=%v", s.Elements())
	}
	if !s.Contains(5) || s.Contains(4) {
		t.Fatal("contains")
	}
	comp := s.Complement()
	if !eqInts(comp.Elements(), []int{0, 2, 4, 6, 7, 8}) {
		t.Fatalf("complement=%v", comp.Elements())
	}
	t2 := StateSetFromSlice(10, []int{3, 9, 7})
	if !eqInts(s.Union(t2).Elements(), []int{1, 3, 5, 7, 9}) {
		t.Fatalf("union=%v", s.Union(t2).Elements())
	}
	if !eqInts(s.Intersect(t2).Elements(), []int{3, 9}) {
		t.Fatalf("intersect=%v", s.Intersect(t2).Elements())
	}
	if !eqInts(s.Difference(t2).Elements(), []int{1, 5}) {
		t.Fatalf("difference=%v", s.Difference(t2).Elements())
	}
	if !StateSetFromSlice(10, []int{3}).IsSubset(s) {
		t.Fatal("subset")
	}
	if s.IsSubset(t2) {
		t.Fatal("not subset")
	}
	if !s.Intersects(t2) {
		t.Fatal("intersects")
	}
	if s.String() != "{1, 3, 5, 9}" {
		t.Fatalf("string=%q", s.String())
	}
	full := FullStateSet(3)
	if !eqInts(full.Elements(), []int{0, 1, 2}) {
		t.Fatal("full")
	}
}

// ---------- Kripke ----------

func TestKripke(t *testing.T) {
	k := sampleModel()
	if k.NumStates() != 4 {
		t.Fatal("num states")
	}
	if !eqInts(k.Successors(0), []int{1, 2}) {
		t.Fatalf("succ0=%v", k.Successors(0))
	}
	if !eqInts(k.Predecessors(3), []int{1, 3}) {
		t.Fatalf("pred3=%v", k.Predecessors(3))
	}
	if !k.Holds(3, "q") || k.Holds(1, "q") {
		t.Fatal("holds")
	}
	if !reflect.DeepEqual(k.Label(3), []string{"p", "q"}) {
		t.Fatalf("label3=%v", k.Label(3))
	}
	if !reflect.DeepEqual(k.Propositions(), []string{"p", "q"}) {
		t.Fatalf("props=%v", k.Propositions())
	}
	if !k.IsTotal() {
		t.Fatal("total")
	}
	if !eqInts(k.Reachable().Elements(), []int{0, 1, 2, 3}) {
		t.Fatalf("reach=%v", k.Reachable().Elements())
	}
	if k.NumTransitions() != 5 {
		t.Fatalf("edges=%d", k.NumTransitions())
	}
	if err := k.Validate(); err != nil {
		t.Fatal(err)
	}
	// deadlock and MakeTotal
	d := NewKripke(2)
	d.AddTransition(0, 1)
	if d.IsTotal() {
		t.Fatal("should have deadlock")
	}
	if n := d.MakeTotal(); n != 1 {
		t.Fatalf("added=%d", n)
	}
	if !d.IsTotal() {
		t.Fatal("now total")
	}
}

// ---------- LTL syntax ----------

func TestLTLParseAndNNF(t *testing.T) {
	tests := []struct {
		in  string
		out string
	}{
		{"p", "p"},
		{"!p", "!p"},
		{"p & q | r", "(p & q) | r"},
		{"G p", "G p"},
		{"p U q", "p U q"},
		{"X X p", "X X p"},
		{"F G p", "F G p"},
	}
	for _, tc := range tests {
		f, err := ParseLTL(tc.in)
		if err != nil {
			t.Fatalf("parse %q: %v", tc.in, err)
		}
		if got := normalizeSpaces(f.String()); got != tc.out {
			t.Errorf("parse %q printed %q want %q", tc.in, got, tc.out)
		}
	}
	// NNF pushes negation to literals only.
	f := MustParseLTL("!(p U q)").NNF()
	if got := normalizeSpaces(f.String()); got != "!p R !q" {
		t.Errorf("NNF !(pUq) = %q", got)
	}
	g := MustParseLTL("!G p").NNF()
	if got := normalizeSpaces(g.String()); got != "true U !p" {
		t.Errorf("NNF !Gp = %q", got)
	}
	if _, err := ParseLTL("p &"); err == nil {
		t.Error("expected parse error")
	}
	if !MustParseLTL("F p").Equal(LTLEventually(LTLVar("p"))) {
		t.Error("equal")
	}
	if MustParseLTL("F p").Atoms()[0] != "p" {
		t.Error("atoms")
	}
}

func TestLTLSimplify(t *testing.T) {
	cases := []struct{ in, out string }{
		{"p & true", "p"},
		{"p | false", "p"},
		{"p & false", "false"},
		{"!!p", "p"},
		{"F F p", "F p"},
		{"true U p", "true U p"},
	}
	for _, c := range cases {
		got := normalizeSpaces(MustParseLTL(c.in).Simplify().String())
		if got != c.out {
			t.Errorf("simplify %q = %q want %q", c.in, got, c.out)
		}
	}
}

// ---------- CTL model checking ----------

func TestCTLCheck(t *testing.T) {
	k := sampleModel()
	cases := []struct {
		f    string
		want []int
	}{
		{"true", []int{0, 1, 2, 3}},
		{"false", nil},
		{"p", []int{1, 3}},
		{"q", []int{3}},
		{"!p", []int{0, 2}},
		{"p & q", []int{3}},
		{"EX p", []int{0, 1, 3}},
		{"AX p", []int{1, 3}},
		{"EF q", []int{0, 1, 3}},
		{"EG p", []int{1, 3}},
		{"AF q", []int{1, 3}},
		{"AG (q -> p)", []int{0, 1, 2, 3}},
		{"E[true U q]", []int{0, 1, 3}},
		{"A[p U q]", []int{1, 3}},
		{"EG true", []int{0, 1, 2, 3}},
	}
	for _, c := range cases {
		set, err := CTLCheck(k, MustParseCTL(c.f))
		if err != nil {
			t.Fatalf("%s: %v", c.f, err)
		}
		if !eqInts(set.Elements(), c.want) {
			t.Errorf("CTLCheck %s = %v want %v", c.f, set.Elements(), c.want)
		}
	}
	// initial-state model checking
	ok, _ := CTLModelCheck(k, MustParseCTL("EF q"))
	if !ok {
		t.Error("EF q should hold at s0")
	}
	ok, _ = CTLModelCheck(k, MustParseCTL("AG p"))
	if ok {
		t.Error("AG p should fail at s0")
	}
}

func TestSatOperatorsDirect(t *testing.T) {
	k := sampleModel()
	p := k.LabelStateSet("p")
	q := k.LabelStateSet("q")
	if !eqInts(SatEX(k, p).Elements(), []int{0, 1, 3}) {
		t.Errorf("SatEX p=%v", SatEX(k, p).Elements())
	}
	if !eqInts(SatEU(k, p, q).Elements(), []int{1, 3}) {
		t.Errorf("SatEU=%v", SatEU(k, p, q).Elements())
	}
	if !eqInts(SatEG(k, p).Elements(), []int{1, 3}) {
		t.Errorf("SatEG=%v", SatEG(k, p).Elements())
	}
	if !eqInts(SatEF(k, q).Elements(), []int{0, 1, 3}) {
		t.Errorf("SatEF=%v", SatEF(k, q).Elements())
	}
	// AR / ER greatest fixpoints sanity: E[false R q] = EG q
	er := SatER(k, NewStateSet(k.n), q)
	eg := SatEG(k, q)
	if !er.Equal(eg) {
		t.Errorf("E[false R q]=%v EG q=%v", er.Elements(), eg.Elements())
	}
}

// ---------- Büchi / emptiness ----------

func TestBuchiEmptiness(t *testing.T) {
	// Automaton with one accepting state on a self loop: nonempty.
	b := NewBuchi(2)
	b.SetInitial(0)
	b.SetAccepting(1)
	b.AddEdge(0, 1, TrueGuard())
	b.AddEdge(1, 1, TrueGuard())
	if b.IsEmpty() {
		t.Error("should be nonempty")
	}
	if b.IsEmptySCC() {
		t.Error("SCC: should be nonempty")
	}
	l, ok := b.AcceptingLasso()
	if !ok || l.LoopStart() != 1 {
		t.Errorf("lasso=%v ok=%v", l, ok)
	}
	// Accepting state not on a cycle: empty.
	b2 := NewBuchi(2)
	b2.SetInitial(0)
	b2.SetAccepting(1)
	b2.AddEdge(0, 1, TrueGuard())
	if !b2.IsEmpty() {
		t.Error("should be empty (no cycle through accepting)")
	}
	if !b2.IsEmptySCC() {
		t.Error("SCC empty")
	}
}

func TestDegeneralize(t *testing.T) {
	// GBA with two acceptance sets {0} and {1} over a 2-cycle 0<->1.
	g := NewGenBuchi(2)
	g.SetInitial(0)
	g.AddEdge(0, 1, TrueGuard())
	g.AddEdge(1, 0, TrueGuard())
	g.AddAcceptanceSet(StateSetFromSlice(2, []int{0}))
	g.AddAcceptanceSet(StateSetFromSlice(2, []int{1}))
	b := Degeneralize(g)
	if b.IsEmpty() {
		t.Error("2-cycle visiting both sets should be nonempty")
	}
}

// ---------- LTL decision problems ----------

func TestLTLDecision(t *testing.T) {
	valid := []string{"p | !p", "G p -> p", "G p -> F p", "G p -> X p"}
	for _, s := range valid {
		if !LTLValid(MustParseLTL(s)) {
			t.Errorf("%s should be valid", s)
		}
	}
	invalid := []string{"F p -> G p", "p", "X p -> p"}
	for _, s := range invalid {
		if LTLValid(MustParseLTL(s)) {
			t.Errorf("%s should not be valid", s)
		}
	}
	unsat := []string{"p & !p", "G p & F !p", "false"}
	for _, s := range unsat {
		if ok, _ := LTLSatisfiable(MustParseLTL(s)); ok {
			t.Errorf("%s should be unsatisfiable", s)
		}
	}
	sat := []string{"F p & G q", "p U q", "G F p"}
	for _, s := range sat {
		if ok, _ := LTLSatisfiable(MustParseLTL(s)); !ok {
			t.Errorf("%s should be satisfiable", s)
		}
	}
	equiv := [][2]string{
		{"F p", "true U p"},
		{"G p", "!F !p"},
		{"!X p", "X !p"},
		{"F F p", "F p"},
		{"p W q", "(p U q) | G p"},
	}
	for _, e := range equiv {
		if !LTLEquivalent(MustParseLTL(e[0]), MustParseLTL(e[1])) {
			t.Errorf("%s != %s", e[0], e[1])
		}
	}
	if LTLEquivalent(MustParseLTL("F p"), MustParseLTL("G p")) {
		t.Error("Fp != Gp")
	}
}

// ---------- LTL model checking + counterexample ----------

func TestLTLModelCheck(t *testing.T) {
	k := sampleModel()
	// G(q -> p) holds on every path (q implies p everywhere).
	holds, cex, err := LTLModelCheck(k, MustParseLTL("G (q -> p)"))
	if err != nil {
		t.Fatal(err)
	}
	if !holds || cex != nil {
		t.Errorf("G(q->p) holds=%v cex=%v", holds, cex)
	}
	// G p fails from s0 (s0 has no p).
	holds, cex, _ = LTLModelCheck(k, MustParseLTL("G p"))
	if holds || cex == nil {
		t.Fatalf("G p should fail with counterexample")
	}
	path := cex.StatePath()
	if len(path) == 0 || path[0] != 0 {
		t.Errorf("counterexample should start at s0: %v", path)
	}
	// F q fails: the path 0->2->2... never reaches q.
	holds, cex, _ = LTLModelCheck(k, MustParseLTL("F q"))
	if holds || cex == nil {
		t.Fatal("F q should fail")
	}
	// consistency with CTL universal duals
	dualCases := [][2]string{{"G p", "AG p"}, {"F q", "AF q"}, {"X p", "AX p"}, {"p U q", "A[p U q]"}}
	for _, dc := range dualCases {
		lh, _, _ := LTLModelCheck(k, MustParseLTL(dc[0]))
		ch, _ := CTLModelCheck(k, MustParseCTL(dc[1]))
		if lh != ch {
			t.Errorf("%s(%v) vs %s(%v)", dc[0], lh, dc[1], ch)
		}
	}
}

func TestLTLvsCTLRandom(t *testing.T) {
	rng := rand.New(rand.NewSource(7))
	props := []string{"p", "q"}
	pairs := [][2]string{{"G p", "AG p"}, {"F q", "AF q"}, {"X q", "AX q"}, {"p U q", "A[p U q]"}}
	for iter := 0; iter < 200; iter++ {
		n := 3 + rng.Intn(3)
		k := NewKripke(n)
		k.SetInitial(0)
		for s := 0; s < n; s++ {
			for d := 0; d < 1+rng.Intn(2); d++ {
				k.AddTransition(s, rng.Intn(n))
			}
			for _, p := range props {
				if rng.Intn(2) == 0 {
					k.AddLabel(s, p)
				}
			}
		}
		k.MakeTotal()
		for _, pr := range pairs {
			lh, _, err := LTLModelCheck(k, MustParseLTL(pr[0]))
			if err != nil {
				t.Fatal(err)
			}
			ch, _ := CTLModelCheck(k, MustParseCTL(pr[1]))
			if lh != ch {
				t.Fatalf("iter %d %s=%v %s=%v\n%s", iter, pr[0], lh, pr[1], ch, k)
			}
		}
	}
}

// ---------- fairness ----------

func TestFairness(t *testing.T) {
	// s0 -> s0, s0 -> s1, s1 -> s1 ; p only at s0.
	k := NewKripke(2)
	k.SetInitial(0)
	k.AddTransition(0, 0)
	k.AddTransition(0, 1)
	k.AddTransition(1, 1)
	k.AddLabel(0, "p")
	fc := NewFairness(StateSetFromSlice(2, []int{1})) // s1 infinitely often
	// Without fairness EG p = {0}; with fairness no fair path stays in p.
	if !eqInts(SatEG(k, k.LabelStateSet("p")).Elements(), []int{0}) {
		t.Errorf("SatEG p=%v", SatEG(k, k.LabelStateSet("p")).Elements())
	}
	fe := FairEG(k, k.LabelStateSet("p"), fc)
	if !fe.IsEmpty() {
		t.Errorf("FairEG p should be empty, got %v", fe.Elements())
	}
	fs := FairStates(k, fc)
	if !eqInts(fs.Elements(), []int{0, 1}) {
		t.Errorf("FairStates=%v", fs.Elements())
	}
	// FairCTLCheck of EG p under fairness is empty; ordinary is {0}.
	fset, err := FairCTLCheck(k, MustParseCTL("EG p"), fc)
	if err != nil {
		t.Fatal(err)
	}
	if !fset.IsEmpty() {
		t.Errorf("fair EG p = %v", fset.Elements())
	}
	// With no fairness constraints, fair check equals ordinary check.
	none := NewFairness()
	for _, s := range []string{"EG p", "EF q", "EX p"} {
		_ = s
	}
	ord, _ := CTLCheck(k, MustParseCTL("EG p"))
	fairNone, _ := FairCTLCheck(k, MustParseCTL("EG p"), none)
	if !ord.Equal(fairNone) {
		t.Errorf("no-fairness mismatch: %v vs %v", ord.Elements(), fairNone.Elements())
	}
}

// ---------- bounded model checking ----------

func TestBMC(t *testing.T) {
	// chain 0->1->2->3->3
	k := NewKripke(4)
	k.SetInitial(0)
	k.AddTransition(0, 1)
	k.AddTransition(1, 2)
	k.AddTransition(2, 3)
	k.AddTransition(3, 3)
	for i := 0; i < 3; i++ {
		k.AddLabel(i, "safe")
	}
	// reachability
	if p, ok := BoundedReachable(k, StateSetFromSlice(4, []int{3}), 3); !ok || !eqInts(p, []int{0, 1, 2, 3}) {
		t.Errorf("reach path=%v ok=%v", p, ok)
	}
	if _, ok := BoundedReachable(k, StateSetFromSlice(4, []int{3}), 2); ok {
		t.Error("should not reach s3 within 2 steps")
	}
	// invariant "safe" violated at s3 within bound 3
	cex, found := BMCInvariant(k, k.LabelStateSet("safe"), 3)
	if !found || cex.Lasso.Prefix[len(cex.Lasso.Prefix)-1] != 3 {
		t.Errorf("expected invariant violation ending at s3: %v", cex)
	}
	// within bound 2 no violation
	if !BMCReachInvariant(k, k.LabelStateSet("safe"), 2) {
		t.Error("no violation within 2 steps")
	}
	// unrolling
	u := Unroll(k, 2)
	if !eqInts(u.ReachableWithin().Elements(), []int{0, 1, 2}) {
		t.Errorf("unroll=%v", u.ReachableWithin().Elements())
	}
	if u.Depth() != 2 {
		t.Errorf("depth=%d", u.Depth())
	}
	// lasso exists (self loop at s3)
	if _, ok := BMCFindLasso(k, 10); !ok {
		t.Error("lasso should exist")
	}
	// EG safe bounded: safe holds at 0,1,2 but not 3, so no infinite safe lasso
	if _, ok := BMCExistsGlobally(k, k.LabelStateSet("safe"), 10); ok {
		t.Error("no infinite path staying in safe")
	}
}

// ---------- counterexample extraction ----------

func TestCounterexampleExtraction(t *testing.T) {
	k := sampleModel()
	// AG p fails at s0; extract path to a !p state.
	cex, ok := AGCounterexample(k, k.LabelStateSet("p"))
	if !ok {
		t.Fatal("AG p should fail")
	}
	last := cex.Lasso.Prefix[len(cex.Lasso.Prefix)-1]
	if k.Holds(last, "p") {
		t.Errorf("counterexample should end at a !p state, ended at s%d", last)
	}
	// AF q fails at s0 via the s2 self-loop.
	af, ok := AFCounterexample(k, k.LabelStateSet("q"))
	if !ok {
		t.Fatal("AF q should fail")
	}
	for _, s := range af.Lasso.States() {
		if k.Holds(s, "q") {
			t.Errorf("AF counterexample visits q at s%d", s)
		}
	}
	// EF q witness from s0
	if p, ok := EFWitness(k, 0, k.LabelStateSet("q")); !ok || p[len(p)-1] != 3 {
		t.Errorf("EF witness=%v ok=%v", p, ok)
	}
	// EG p witness from s1 (loops in p)
	l, ok := EGWitness(k, 1, k.LabelStateSet("p"))
	if !ok {
		t.Fatal("EG p witness from s1")
	}
	for _, s := range l.States() {
		if !k.Holds(s, "p") {
			t.Errorf("EG witness leaves p at s%d", s)
		}
	}
	// EU witness
	if p, ok := EUWitness(k, 1, k.LabelStateSet("p"), k.LabelStateSet("q")); !ok || p[len(p)-1] != 3 {
		t.Errorf("EU witness=%v ok=%v", p, ok)
	}
}

// ---------- guards ----------

func TestGuard(t *testing.T) {
	g := NewGuard([]string{"a"}, []string{"b"})
	if !g.Satisfies(map[string]bool{"a": true}) {
		t.Error("a & !b satisfied by {a}")
	}
	if g.Satisfies(map[string]bool{"a": true, "b": true}) {
		t.Error("a & !b not satisfied by {a,b}")
	}
	if TrueGuard().Contradictory() {
		t.Error("true guard not contradictory")
	}
	if !NewGuard([]string{"a"}, []string{"a"}).Contradictory() {
		t.Error("a & !a contradictory")
	}
	if g.String() != "a & !b" {
		t.Errorf("guard string=%q", g.String())
	}
}

// ---------- examples ----------

func Example() {
	k := NewKripke(3)
	k.SetInitial(0)
	k.AddTransition(0, 1)
	k.AddTransition(1, 2)
	k.AddTransition(2, 2)
	k.AddLabel(1, "p")
	k.AddLabel(2, "q")

	set, _ := CTLCheck(k, MustParseCTL("EF q"))
	fmt.Println("EF q:", set)

	ok, _ := CTLModelCheck(k, MustParseCTL("AF q"))
	fmt.Println("AF q holds:", ok)

	fmt.Println("G p -> p valid:", LTLValid(MustParseLTL("G p -> p")))
	// Output:
	// EF q: {0, 1, 2}
	// AF q holds: true
	// G p -> p valid: true
}

func ExampleLTLModelCheck() {
	k := NewKripke(2)
	k.SetInitial(0)
	k.AddTransition(0, 1)
	k.AddTransition(1, 1)
	k.AddLabel(0, "p")
	k.AddLabel(1, "p")

	holds, _, _ := LTLModelCheck(k, MustParseLTL("G p"))
	fmt.Println("G p holds:", holds)
	// Output:
	// G p holds: true
}
