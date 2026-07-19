package proofsystems

import (
	"fmt"
	"sort"
	"testing"
)

// ---------- terms & unification ----------

func TestUnifyTerms(t *testing.T) {
	x, y := NewVar("x"), NewVar("y")
	a, b := NewConst("a"), NewConst("b")
	cases := []struct {
		name    string
		s, t    Term
		unifies bool
		xImage  string // expected image of x, "" to skip
	}{
		{"vars", x, y, true, "y"},
		{"var-const", x, a, true, "a"},
		{"func-match", NewFunc("f", x, a), NewFunc("f", b, y), true, "b"},
		{"func-clash", NewFunc("f", a), NewFunc("g", a), false, ""},
		{"arity-clash", NewFunc("f", a), NewFunc("f", a, b), false, ""},
		{"occurs-check", x, NewFunc("f", x), false, ""},
		{"const-clash", a, b, false, ""},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			s, err := UnifyTerms(c.s, c.t)
			if c.unifies != (err == nil) {
				t.Fatalf("unifies=%v, err=%v", c.unifies, err)
			}
			if err != nil {
				return
			}
			if !s.ApplyTerm(c.s).Equal(s.ApplyTerm(c.t)) {
				t.Fatalf("unifier does not equate terms: %s", s)
			}
			if c.xImage != "" {
				img, _ := s.Get("x")
				if img.String() != c.xImage {
					t.Fatalf("x -> %s, want %s", img.String(), c.xImage)
				}
			}
		})
	}
}

func TestMatchTerm(t *testing.T) {
	x := NewVar("x")
	a, b := NewConst("a"), NewConst("b")
	if !Matches(NewFunc("f", x, x), NewFunc("f", a, a)) {
		t.Errorf("f(x,x) should match f(a,a)")
	}
	if Matches(NewFunc("f", x, x), NewFunc("f", a, b)) {
		t.Errorf("f(x,x) should not match f(a,b)")
	}
	// matching is one-way: the subject variables are not bound
	if Matches(NewConst("a"), NewVar("y")) {
		t.Errorf("constant should not match a free variable in one-way match")
	}
}

func TestSubstitutionCompose(t *testing.T) {
	// s = {x -> f(y)}, o = {y -> a}; (o∘s)(x) = f(a)
	s := SingletonSubstitution("x", NewFunc("f", NewVar("y")))
	o := SingletonSubstitution("y", NewConst("a"))
	comp := s.Compose(o)
	got := comp.ApplyTerm(NewVar("x"))
	if got.String() != "f(a)" {
		t.Fatalf("composition gave %s, want f(a)", got)
	}
}

// ---------- parsing, evaluation, semantics ----------

func TestParseAndEval(t *testing.T) {
	f := MustParseFormula("(A -> B) & A")
	a := Assignment{"A": true, "B": true}
	if v := MustEval(f, a); !v {
		t.Errorf("expected true")
	}
	a2 := Assignment{"A": true, "B": false}
	if v := MustEval(f, a2); v {
		t.Errorf("expected false")
	}
}

func TestSemantics(t *testing.T) {
	cases := []struct {
		src                       string
		taut, contra, sat, contin bool
	}{
		{"A | !A", true, false, true, false},
		{"A & !A", false, true, false, false},
		{"A -> A", true, false, true, false},
		{"A -> B", false, false, true, true},
		{"(A -> B) <-> (!B -> !A)", true, false, true, false},
		{"((A -> B) -> A) -> A", true, false, true, false},
	}
	for _, c := range cases {
		t.Run(c.src, func(t *testing.T) {
			f := MustParseFormula(c.src)
			if IsTautology(f) != c.taut {
				t.Errorf("IsTautology=%v want %v", IsTautology(f), c.taut)
			}
			if IsContradiction(f) != c.contra {
				t.Errorf("IsContradiction=%v want %v", IsContradiction(f), c.contra)
			}
			if IsSatisfiable(f) != c.sat {
				t.Errorf("IsSatisfiable=%v want %v", IsSatisfiable(f), c.sat)
			}
			if IsContingent(f) != c.contin {
				t.Errorf("IsContingent=%v want %v", IsContingent(f), c.contin)
			}
		})
	}
}

func TestEntailmentAndEquivalence(t *testing.T) {
	prem := []Formula{MustParseFormula("A"), MustParseFormula("A -> B")}
	if !Entails(prem, MustParseFormula("B")) {
		t.Errorf("modus ponens should entail B")
	}
	if Entails(prem, MustParseFormula("C")) {
		t.Errorf("should not entail unrelated C")
	}
	if !Equivalent(MustParseFormula("A -> B"), MustParseFormula("!A | B")) {
		t.Errorf("A->B equivalent to !A|B")
	}
	if !Equivalent(MustParseFormula("!(A & B)"), MustParseFormula("!A | !B")) {
		t.Errorf("De Morgan")
	}
	if CountModels(MustParseFormula("A | B")) != 3 {
		t.Errorf("A|B has 3 models, got %d", CountModels(MustParseFormula("A | B")))
	}
}

func TestTruthTable(t *testing.T) {
	tt, err := NewTruthTable(MustParseFormula("A & B"))
	if err != nil {
		t.Fatal(err)
	}
	if len(tt.Rows) != 4 {
		t.Fatalf("expected 4 rows, got %d", len(tt.Rows))
	}
	if tt.CountTrue() != 1 {
		t.Errorf("A&B true in 1 row, got %d", tt.CountTrue())
	}
}

// ---------- normal forms ----------

func TestNNFandCNF(t *testing.T) {
	f := MustParseFormula("!(A -> B)")
	nnf := ToNNF(f)
	// !(A->B) == A & !B
	if !Equivalent(nnf, MustParseFormula("A & !B")) {
		t.Errorf("NNF of !(A->B) wrong: %s", nnf)
	}
	cnf := ToCNF(MustParseFormula("A | (B & C)"))
	if len(cnf.Clauses) != 2 {
		t.Errorf("expected 2 clauses, got %d: %s", len(cnf.Clauses), cnf)
	}
	// clausal form must be equivalent to the original
	for _, a := range AllAssignments([]string{"A", "B", "C"}) {
		orig, _ := Eval(MustParseFormula("A | (B & C)"), a)
		if cnf.Eval(a) != orig {
			t.Errorf("CNF disagrees at %s", a)
		}
	}
}

func TestHornDetection(t *testing.T) {
	horn := NewPClause(NegPLit("A"), NegPLit("B"), PosPLit("C"))
	notHorn := NewPClause(PosPLit("A"), PosPLit("B"))
	if !IsHornClause(horn) {
		t.Errorf("expected Horn clause")
	}
	if IsHornClause(notHorn) {
		t.Errorf("two positive literals is not Horn")
	}
}

// ---------- Tseitin + DPLL ----------

func TestTseitinEquisatisfiable(t *testing.T) {
	srcs := []string{
		"A & B",
		"(A | B) & (!A | !B)",
		"A & !A",
		"(A -> B) & (B -> C) & A & !C",
		"((A | B) & (C | D)) -> (A | C)",
	}
	for _, src := range srcs {
		f := MustParseFormula(src)
		want := IsSatisfiable(f)
		got := DPLLSatisfiable(TseitinCNF(f, "_t"))
		if want != got {
			t.Errorf("%s: DPLL(Tseitin)=%v, truth-table sat=%v", src, got, want)
		}
	}
}

func TestDPLL(t *testing.T) {
	sat := ToCNF(MustParseFormula("(A | B) & (!A | B)"))
	r := DPLL(sat)
	if !r.Sat {
		t.Fatalf("expected satisfiable")
	}
	if !sat.Eval(r.Model) {
		t.Errorf("returned model does not satisfy CNF: %s", r.Model)
	}
	unsat := ToCNF(MustParseFormula("A & !A"))
	if DPLLSatisfiable(unsat) {
		t.Errorf("A & !A must be unsatisfiable")
	}
}

func TestUnitPropagation(t *testing.T) {
	n := NewPCNF(
		NewPClause(PosPLit("A")),
		NewPClause(NegPLit("A"), PosPLit("B")),
	)
	_, assign, ok := UnitPropagate(n)
	if !ok {
		t.Fatalf("no conflict expected")
	}
	if !assign["A"] || !assign["B"] {
		t.Errorf("expected A and B forced true, got %s", assign)
	}
}

// ---------- resolution & tableau ----------

func TestResolution(t *testing.T) {
	if !ResolutionValid(MustParseFormula("A | !A")) {
		t.Errorf("excluded middle should be valid by resolution")
	}
	if ResolutionValid(MustParseFormula("A")) {
		t.Errorf("bare A is not valid")
	}
	prem := []Formula{MustParseFormula("A -> B"), MustParseFormula("B -> C"), MustParseFormula("A")}
	if !ResolutionEntails(prem, MustParseFormula("C")) {
		t.Errorf("hypothetical syllogism should hold")
	}
}

func TestTableau(t *testing.T) {
	cases := []struct {
		src   string
		valid bool
	}{
		{"A | !A", true},
		{"(A -> B) -> (!B -> !A)", true},
		{"((A -> B) -> A) -> A", true},
		{"A -> B", false},
		{"A & !A -> B", true},
	}
	for _, c := range cases {
		if got := TableauValid(MustParseFormula(c.src)); got != c.valid {
			t.Errorf("%s: TableauValid=%v want %v", c.src, got, c.valid)
		}
	}
}

// ---------- first-order clausification & resolution ----------

func TestFOResolutionSyllogism(t *testing.T) {
	// forall x. Man(x) -> Mortal(x);  Man(socrates)  |=  Mortal(socrates)
	allMortal := Forall("x", Imp(Atom("Man", NewVar("x")), Atom("Mortal", NewVar("x"))))
	manSocrates := Atom("Man", NewConst("socrates"))
	mortalSocrates := Atom("Mortal", NewConst("socrates"))
	if !FOEntails([]Formula{allMortal, manSocrates}, mortalSocrates, 100) {
		t.Errorf("Socrates syllogism should be entailed")
	}
	// negative control: does not entail Mortal(plato) without Man(plato)
	if FOEntails([]Formula{allMortal, manSocrates}, Atom("Mortal", NewConst("plato")), 60) {
		t.Errorf("should not entail Mortal(plato)")
	}
}

func TestClausify(t *testing.T) {
	// exists x. forall y. P(x,y) skolemizes x to a constant
	f := Exists("x", Forall("y", Atom("P", NewVar("x"), NewVar("y"))))
	clauses := Clausify(f)
	if len(clauses) != 1 || len(clauses[0].Lits) != 1 {
		t.Fatalf("expected single unit clause, got %v", clauses)
	}
	lit := clauses[0].Lits[0]
	if lit.Pred != "P" || !lit.Args[0].IsConst() {
		t.Errorf("first argument should be a Skolem constant, got %s", lit)
	}
}

// ---------- SLD / Prolog ----------

func TestSLDResolution(t *testing.T) {
	tom, bob, ann := NewConst("tom"), NewConst("bob"), NewConst("ann")
	X, Y, Z := NewVar("x"), NewVar("y"), NewVar("z")
	p := NewProgram(
		Fact(NewFOLiteral("parent", tom, bob)),
		Fact(NewFOLiteral("parent", bob, ann)),
		Rule(NewFOLiteral("ancestor", X, Y), NewFOLiteral("parent", X, Y)),
		Rule(NewFOLiteral("ancestor", X, Y),
			NewFOLiteral("parent", X, Z), NewFOLiteral("ancestor", Z, Y)),
	)
	if !p.Query(NewFOLiteral("ancestor", tom, ann)) {
		t.Errorf("tom should be an ancestor of ann")
	}
	if p.Query(NewFOLiteral("ancestor", ann, tom)) {
		t.Errorf("ann should not be an ancestor of tom")
	}
	// enumerate descendants of tom
	ans := p.Solve([]FOLiteral{NewFOLiteral("ancestor", tom, NewVar("w"))}, DefaultSolveOptions())
	got := map[string]bool{}
	for _, s := range ans {
		img, _ := s.Get("w")
		got[img.String()] = true
	}
	if !got["bob"] || !got["ann"] || len(got) != 2 {
		t.Errorf("expected descendants {bob, ann}, got %v", keysOf(got))
	}
}

// ---------- natural deduction ----------

func TestNaturalDeduction(t *testing.T) {
	A := Prop("A")
	B := Prop("B")

	// theorem: A -> A
	idProof := ImpIntro(A, Assume(A))
	if !NDProves(idProof) {
		t.Errorf("A -> A should be provable with no open assumptions")
	}

	// modus ponens depends on assumptions A and A->B
	mp := ImpElim(Assume(A), Assume(Imp(A, B)))
	if !NDValid(mp, []Formula{A, Imp(A, B)}) {
		t.Errorf("modus ponens should be valid from its premises")
	}
	if NDProves(mp) {
		t.Errorf("modus ponens is not a closed theorem")
	}

	// theorem: (A & B) -> (B & A)
	ab := And(A, B)
	body := AndIntro(AndElim2(Assume(ab)), AndElim1(Assume(ab)))
	commProof := ImpIntro(ab, body)
	if !commProof.Concl.Equal(Imp(ab, And(B, A))) {
		t.Fatalf("unexpected conclusion %s", commProof.Concl)
	}
	if !NDProves(commProof) {
		t.Errorf("conjunction commutativity should be a theorem")
	}

	// an unsound step must be rejected
	bad := Derivation{Rule: NDAndE1, Concl: A, Premises: []Derivation{Assume(A)}}
	if _, err := CheckND(bad); err == nil {
		t.Errorf("∧E1 applied to a non-conjunction must fail")
	}
}

// ---------- sequent calculus ----------

func TestSequentProver(t *testing.T) {
	peirce := MustParseFormula("((A -> B) -> A) -> A")
	d, ok := ProveFormula(peirce)
	if !ok {
		t.Fatalf("Peirce's law should be provable in LK")
	}
	if err := CheckSequent(d); err != nil {
		t.Errorf("generated proof failed the checker: %v", err)
	}
	if !ValidSequent(Sequent{Right: []Formula{peirce}}) {
		t.Errorf("Peirce sequent should be valid")
	}
	// a non-theorem is not provable
	if _, ok := ProveFormula(Prop("A")); ok {
		t.Errorf("bare A is not a theorem")
	}
	// classic entailment sequent A, A->B |- B
	seq := Sequent{
		Left:  []Formula{Prop("A"), MustParseFormula("A -> B")},
		Right: []Formula{Prop("B")},
	}
	if !SequentProvable(seq) {
		t.Errorf("modus ponens sequent should be provable")
	}
}

func TestCheckSequentRejectsBadProof(t *testing.T) {
	// claim A |- B by a bogus axiom
	bad := SequentDerivation{
		Rule:  SeqAxiom,
		Concl: Sequent{Left: []Formula{Prop("A")}, Right: []Formula{Prop("B")}},
	}
	if err := CheckSequent(bad); err == nil {
		t.Errorf("A |- B is not an axiom and must be rejected")
	}
}

// ---------- helpers & example ----------

func keysOf(m map[string]bool) []string {
	out := make([]string, 0, len(m))
	for k := range m {
		out = append(out, k)
	}
	sort.Strings(out)
	return out
}

func ExampleEntails() {
	premises := []Formula{
		MustParseFormula("A"),
		MustParseFormula("A -> B"),
	}
	conclusion := MustParseFormula("B")
	fmt.Println(Entails(premises, conclusion))
	// Output: true
}

func ExampleUnifyTerms() {
	// unify f(x, a) with f(b, y)
	t1 := NewFunc("f", NewVar("x"), NewConst("a"))
	t2 := NewFunc("f", NewConst("b"), NewVar("y"))
	s, _ := UnifyTerms(t1, t2)
	fmt.Println(s)
	// Output: {x -> b, y -> a}
}

func ExampleProgram_Query() {
	tom, bob := NewConst("tom"), NewConst("bob")
	x, y := NewVar("x"), NewVar("y")
	p := NewProgram(
		Fact(NewFOLiteral("parent", tom, bob)),
		Rule(NewFOLiteral("ancestor", x, y), NewFOLiteral("parent", x, y)),
	)
	fmt.Println(p.Query(NewFOLiteral("ancestor", tom, bob)))
	// Output: true
}
