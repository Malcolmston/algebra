package satsolver

import (
	"fmt"
	"reflect"
	"sort"
	"testing"
)

func TestLitBasics(t *testing.T) {
	tests := []struct {
		l       Lit
		wantVar int
		wantNeg bool
		wantStr string
	}{
		{PosLit(3), 3, false, "x3"},
		{NegLit(3), 3, true, "~x3"},
		{MakeLit(5, true), 5, true, "~x5"},
		{MakeLit(5, false), 5, false, "x5"},
	}
	for _, tc := range tests {
		if got := tc.l.Var(); got != tc.wantVar {
			t.Errorf("%v.Var()=%d want %d", tc.l, got, tc.wantVar)
		}
		if got := tc.l.IsNeg(); got != tc.wantNeg {
			t.Errorf("%v.IsNeg()=%v want %v", tc.l, got, tc.wantNeg)
		}
		if got := tc.l.String(); got != tc.wantStr {
			t.Errorf("%v.String()=%q want %q", tc.l, got, tc.wantStr)
		}
		if got := tc.l.Negate().Negate(); got != tc.l {
			t.Errorf("double negate changed literal: %v", got)
		}
	}
}

func TestClauseEval(t *testing.T) {
	c := NewClause(PosLit(1), NegLit(2))
	tests := []struct {
		assign map[int]bool
		wantOr bool
		wantAn bool
	}{
		{map[int]bool{1: true, 2: true}, true, false},
		{map[int]bool{1: false, 2: true}, false, false},
		{map[int]bool{1: false, 2: false}, true, false},
		{map[int]bool{1: true, 2: false}, true, true},
	}
	for _, tc := range tests {
		if got := c.EvalOr(tc.assign); got != tc.wantOr {
			t.Errorf("EvalOr(%v)=%v want %v", tc.assign, got, tc.wantOr)
		}
		if got := c.EvalAnd(tc.assign); got != tc.wantAn {
			t.Errorf("EvalAnd(%v)=%v want %v", tc.assign, got, tc.wantAn)
		}
	}
	if !NewClause(PosLit(1), NegLit(1)).IsTautology() {
		t.Errorf("clause with x and ~x should be tautology")
	}
	if NewClause(PosLit(1), PosLit(2)).IsTautology() {
		t.Errorf("clause without complementary pair is not a tautology")
	}
}

func TestCNFEvalAndAssign(t *testing.T) {
	// (x1 | x2) & (~x1 | x3)
	f := NewCNF(
		NewClause(PosLit(1), PosLit(2)),
		NewClause(NegLit(1), PosLit(3)),
	)
	if !f.Eval(map[int]bool{1: true, 2: false, 3: true}) {
		t.Errorf("expected satisfying assignment to evaluate true")
	}
	if f.Eval(map[int]bool{1: true, 2: false, 3: false}) {
		t.Errorf("expected falsifying assignment to evaluate false")
	}
	// Conditioning on x1=true removes clause 1 and shortens clause 2.
	g := f.Assign(PosLit(1))
	if g.NumClauses() != 1 {
		t.Fatalf("after Assign(x1) expected 1 clause, got %d", g.NumClauses())
	}
	if !reflect.DeepEqual(g.Clauses[0], NewClause(PosLit(3))) {
		t.Errorf("after Assign(x1) expected clause (x3), got %v", g.Clauses[0])
	}
}

func TestDIMACSRoundTrip(t *testing.T) {
	f := NewCNF(
		NewClause(PosLit(1), NegLit(2), PosLit(3)),
		NewClause(NegLit(1), PosLit(2)),
	)
	s := f.DIMACS()
	g, err := ParseDIMACS(s)
	if err != nil {
		t.Fatalf("ParseDIMACS error: %v", err)
	}
	if !EquivalentCNF(f, g) {
		t.Errorf("DIMACS round trip changed the formula:\n%s", s)
	}
	if len(g.Clauses) != 2 {
		t.Errorf("expected 2 clauses after parse, got %d", len(g.Clauses))
	}
}

func TestParseEval(t *testing.T) {
	tests := []struct {
		src  string
		env  map[string]bool
		want bool
	}{
		{"a & b", map[string]bool{"a": true, "b": true}, true},
		{"a & b", map[string]bool{"a": true, "b": false}, false},
		{"a | b", map[string]bool{"a": false, "b": true}, true},
		{"a -> b", map[string]bool{"a": true, "b": false}, false},
		{"a -> b", map[string]bool{"a": false, "b": false}, true},
		{"a <-> b", map[string]bool{"a": true, "b": true}, true},
		{"a ^ b", map[string]bool{"a": true, "b": true}, false},
		{"~a", map[string]bool{"a": false}, true},
		{"a & b | c", map[string]bool{"a": false, "b": false, "c": true}, true},
		{"a & (b | c)", map[string]bool{"a": true, "b": false, "c": false}, false},
	}
	for _, tc := range tests {
		got, err := EvalString(tc.src, tc.env)
		if err != nil {
			t.Errorf("EvalString(%q) error: %v", tc.src, err)
			continue
		}
		if got != tc.want {
			t.Errorf("EvalString(%q, %v)=%v want %v", tc.src, tc.env, got, tc.want)
		}
	}
}

func TestParsePrecedence(t *testing.T) {
	// & binds tighter than |
	e := MustParse("a | b & c")
	want := Or{X: Variable("a"), Y: And{X: Variable("b"), Y: Variable("c")}}
	if !ExprEqual(e, want) {
		t.Errorf("precedence parse = %s want %s", e, want)
	}
}

func TestSimplify(t *testing.T) {
	tests := []struct {
		src  string
		want Expr
	}{
		{"a & a", Variable("a")},
		{"a | a", Variable("a")},
		{"a & ~a", False},
		{"a | ~a", True},
		{"a & T", Variable("a")},
		{"a & F", False},
		{"a | T", True},
		{"a | F", Variable("a")},
		{"~~a", Variable("a")},
	}
	for _, tc := range tests {
		got := Simplify(MustParse(tc.src))
		if !ExprEqual(got, tc.want) {
			t.Errorf("Simplify(%q)=%s want %s", tc.src, got, tc.want)
		}
	}
}

func TestSimplifyPreservesSemantics(t *testing.T) {
	srcs := []string{
		"(a -> b) & (b -> c)",
		"a ^ b ^ c",
		"(a <-> b) | (c & ~a)",
		"~(a & b) | (a & ~c)",
	}
	for _, s := range srcs {
		e := MustParse(s)
		if !Equivalent(e, Simplify(e)) {
			t.Errorf("Simplify changed semantics of %q", s)
		}
	}
}

func TestNNFAndCNFEquivalence(t *testing.T) {
	srcs := []string{
		"a -> (b -> c)",
		"(a <-> b)",
		"~(a | (b & c))",
		"a ^ b",
		"(a & b) | (c & ~d)",
	}
	for _, s := range srcs {
		e := MustParse(s)
		nnf := ToNNF(e)
		if !IsNNF(nnf) {
			t.Errorf("ToNNF(%q) not in NNF: %s", s, nnf)
		}
		if !Equivalent(e, nnf) {
			t.Errorf("ToNNF changed semantics of %q", s)
		}
		cnf := ToCNFExpr(e)
		if !IsCNFExpr(cnf) {
			t.Errorf("ToCNFExpr(%q) not in CNF: %s", s, cnf)
		}
		if !Equivalent(e, cnf) {
			t.Errorf("ToCNFExpr changed semantics of %q", s)
		}
		dnf := ToDNFExpr(e)
		if !IsDNFExpr(dnf) {
			t.Errorf("ToDNFExpr(%q) not in DNF: %s", s, dnf)
		}
		if !Equivalent(e, dnf) {
			t.Errorf("ToDNFExpr changed semantics of %q", s)
		}
	}
}

func TestToCNFFormulaEquivalence(t *testing.T) {
	e := MustParse("(a <-> b) & (c | ~a)")
	f, vm := ToCNFFormula(e)
	// Reconstruct and compare truth tables over the original variables.
	rebuilt := CNFToExpr(f, vm)
	if !Equivalent(e, rebuilt) {
		t.Errorf("ToCNFFormula/CNFToExpr changed semantics: %s vs %s", e, rebuilt)
	}
}

func TestTseitinEquisatisfiable(t *testing.T) {
	tests := []struct {
		src     string
		wantSat bool
	}{
		{"(a <-> b) & (b <-> c)", true},
		{"a & ~a", false},
		{"(a | b) & (~a | b) & (a | ~b) & (~a | ~b)", false},
		{"a ^ b ^ c", true},
		{"(a -> b) & (b -> c) & a & ~c", false},
	}
	for _, tc := range tests {
		e := MustParse(tc.src)
		f := TseitinCNF(e)
		if got := IsSatisfiableCNF(f); got != tc.wantSat {
			t.Errorf("Tseitin sat of %q = %v want %v", tc.src, got, tc.wantSat)
		}
		// Cross-check against direct semantic satisfiability.
		if IsSatisfiable(e) != tc.wantSat {
			t.Errorf("IsSatisfiable(%q)=%v want %v", tc.src, IsSatisfiable(e), tc.wantSat)
		}
	}
}

func TestDPLL(t *testing.T) {
	tests := []struct {
		name    string
		f       CNF
		wantSat bool
	}{
		{
			"sat",
			NewCNF(NewClause(PosLit(1), PosLit(2)), NewClause(NegLit(1), PosLit(3))),
			true,
		},
		{
			"unsat-unit-conflict",
			NewCNF(
				NewClause(PosLit(1), PosLit(2)),
				NewClause(NegLit(1), PosLit(2)),
				NewClause(NegLit(2)),
			),
			false,
		},
		{
			"empty-formula-sat",
			NewCNF(),
			true,
		},
		{
			"empty-clause-unsat",
			NewCNF(NewClause()),
			false,
		},
	}
	for _, tc := range tests {
		model, ok := SolveCNF(tc.f)
		if ok != tc.wantSat {
			t.Errorf("%s: SolveCNF sat=%v want %v", tc.name, ok, tc.wantSat)
			continue
		}
		if ok && !tc.f.Eval(model) {
			t.Errorf("%s: returned model %v does not satisfy formula", tc.name, model)
		}
	}
}

func TestUnitPropagate(t *testing.T) {
	f := NewCNF(NewClause(PosLit(1)), NewClause(NegLit(1), PosLit(2)))
	g, assign, ok := UnitPropagate(f, map[int]bool{})
	if !ok {
		t.Fatalf("expected no conflict")
	}
	if assign[1] != true || assign[2] != true {
		t.Errorf("unit propagation assignment = %v, want x1=x2=true", assign)
	}
	if !g.IsEmpty() {
		t.Errorf("expected all clauses satisfied, got %s", g)
	}
}

func TestPureLiteral(t *testing.T) {
	// x1 is pure positive, x2 pure negative.
	f := NewCNF(NewClause(PosLit(1), PosLit(2).Negate()), NewClause(PosLit(1), NegLit(3)), NewClause(PosLit(3), NegLit(2)))
	pures := f.PureLiterals()
	got := map[Lit]bool{}
	for _, l := range pures {
		got[l] = true
	}
	if !got[PosLit(1)] {
		t.Errorf("expected x1 pure positive in %v", pures)
	}
	if !got[NegLit(2)] {
		t.Errorf("expected ~x2 pure negative in %v", pures)
	}
}

func TestResolve(t *testing.T) {
	c1 := NewClause(PosLit(1), PosLit(2))
	c2 := NewClause(NegLit(1), PosLit(3))
	r, ok := Resolve(c1, c2, 1)
	if !ok {
		t.Fatalf("expected resolvable clauses")
	}
	want := NewClause(PosLit(2), PosLit(3)).Sorted()
	if !reflect.DeepEqual(r, want) {
		t.Errorf("Resolve = %v want %v", r, want)
	}
	if _, ok := Resolve(c1, NewClause(PosLit(3)), 1); ok {
		t.Errorf("expected non-resolvable on variable 1")
	}
}

func TestSemantics(t *testing.T) {
	tests := []struct {
		src     string
		taut    bool
		contra  bool
		sat     bool
		count   int
		numVars int
	}{
		{"a | ~a", true, false, true, 2, 1},
		{"a & ~a", false, true, false, 0, 1},
		{"a & b", false, false, true, 1, 2},
		{"a -> b", false, false, true, 3, 2},
		{"a ^ b", false, false, true, 2, 2},
	}
	for _, tc := range tests {
		e := MustParse(tc.src)
		if got := IsTautology(e); got != tc.taut {
			t.Errorf("IsTautology(%q)=%v want %v", tc.src, got, tc.taut)
		}
		if got := IsContradiction(e); got != tc.contra {
			t.Errorf("IsContradiction(%q)=%v want %v", tc.src, got, tc.contra)
		}
		if got := IsSatisfiable(e); got != tc.sat {
			t.Errorf("IsSatisfiable(%q)=%v want %v", tc.src, got, tc.sat)
		}
		if got := CountModels(e); got != tc.count {
			t.Errorf("CountModels(%q)=%d want %d", tc.src, got, tc.count)
		}
	}
}

func TestEquivalentAndEntails(t *testing.T) {
	if !Equivalent(MustParse("a -> b"), MustParse("~a | b")) {
		t.Errorf("a->b should be equivalent to ~a|b")
	}
	if !Equivalent(MustParse("~(a & b)"), MustParse("~a | ~b")) {
		t.Errorf("De Morgan equivalence failed")
	}
	if !Entails(MustParse("a & b"), MustParse("a")) {
		t.Errorf("a&b should entail a")
	}
	if Entails(MustParse("a"), MustParse("a & b")) {
		t.Errorf("a should not entail a&b")
	}
}

func TestTruthTable(t *testing.T) {
	e := MustParse("a & b")
	tt := NewTruthTable(e)
	if tt.NumRows() != 4 {
		t.Fatalf("expected 4 rows, got %d", tt.NumRows())
	}
	if got := tt.Minterms(); !reflect.DeepEqual(got, []int{3}) {
		t.Errorf("minterms of a&b = %v want [3]", got)
	}
	if got := Maxterms(e); !reflect.DeepEqual(got, []int{0, 1, 2}) {
		t.Errorf("maxterms of a&b = %v want [0 1 2]", got)
	}
	// Index/env round trip.
	vars := []string{"a", "b"}
	for i := 0; i < 4; i++ {
		env := IndexToEnv(vars, i)
		if EnvToIndex(vars, env) != i {
			t.Errorf("index round trip failed at %d", i)
		}
	}
}

func TestBDDSatCountMatchesModels(t *testing.T) {
	srcs := []string{
		"a & b | ~a & c",
		"(a <-> b) & (c | d)",
		"a ^ b ^ c",
		"a -> (b & c)",
	}
	for _, s := range srcs {
		e := MustParse(s)
		b, root := NewBDDFromExpr(e)
		if got, want := b.SatCount(root), CountModels(e); got != want {
			t.Errorf("BDD SatCount(%q)=%d want %d", s, got, want)
		}
		// Evaluate BDD against expression on all assignments.
		vars := Vars(e)
		n := len(vars)
		for mask := 0; mask < (1 << n); mask++ {
			env := indexEnv(vars, mask)
			if b.Eval(root, env) != e.Eval(env) {
				t.Errorf("BDD eval mismatch for %q at %v", s, env)
			}
		}
	}
}

func TestBDDTautologyContradiction(t *testing.T) {
	b := NewBDD([]string{"a", "b"})
	taut := b.Or(b.Var("a"), b.Not(b.Var("a")))
	if !b.IsTautology(taut) {
		t.Errorf("a | ~a should be BDD true terminal")
	}
	contra := b.And(b.Var("a"), b.Not(b.Var("a")))
	if !b.IsContradiction(contra) {
		t.Errorf("a & ~a should be BDD false terminal")
	}
	// Canonicity: two syntactically different but equal functions share a node.
	f1 := b.And(b.Var("a"), b.Var("b"))
	f2 := b.Not(b.Or(b.Not(b.Var("a")), b.Not(b.Var("b"))))
	if !b.Equal(f1, f2) {
		t.Errorf("De Morgan-equal functions should be the same BDD node")
	}
}

func TestBDDAnySat(t *testing.T) {
	e := MustParse("a & ~b & c")
	b, root := NewBDDFromExpr(e)
	env, ok := b.AnySat(root)
	if !ok {
		t.Fatalf("expected satisfiable")
	}
	if !e.Eval(env) {
		t.Errorf("AnySat returned non-satisfying env %v", env)
	}
}

func TestQuineMcCluskeyNotC(t *testing.T) {
	// f(a,b,c) true on even minterms => f = ~c.
	minterms := []int{0, 2, 4, 6}
	cover := MinimizeSOP(minterms, nil, 3)
	if len(cover) != 1 {
		t.Fatalf("expected single implicant cover, got %d: %v", len(cover), cover)
	}
	if got := cover[0].String(); got != "--0" {
		t.Errorf("cover implicant = %q want \"--0\"", got)
	}
	names := []string{"a", "b", "c"}
	if got := SOPString(cover, names); got != "~c" {
		t.Errorf("SOPString = %q want \"~c\"", got)
	}
	// Reconstructed expression must equal ~c.
	expr := MinimizeSOPExpr(minterms, nil, names)
	if !Equivalent(expr, Not{X: Variable("c")}) {
		t.Errorf("minimised expr %s not equivalent to ~c", expr)
	}
}

func TestQuineMcCluskeyXor(t *testing.T) {
	// f(a,b) = a xor b : minterms 1 (01) and 2 (10).
	minterms := []int{1, 2}
	cover := MinimizeSOP(minterms, nil, 2)
	if len(cover) != 2 {
		t.Fatalf("expected 2 implicants for xor, got %d", len(cover))
	}
	expr := MinimizeSOPExpr(minterms, []int{}, []string{"a", "b"})
	if !Equivalent(expr, Xor{X: Variable("a"), Y: Variable("b")}) {
		t.Errorf("minimised %s not equivalent to a xor b", expr)
	}
}

func TestQuineMcCluskeyDontCares(t *testing.T) {
	// Classic: minterms {1,3,7,11,15}, don't-cares {0,2,5} over 4 vars.
	minterms := []int{1, 3, 7, 11, 15}
	dc := []int{0, 2, 5}
	cover := MinimizeSOP(minterms, dc, 4)
	// Every required minterm must be covered and no off-set minterm covered.
	offset := map[int]bool{}
	for m := 0; m < 16; m++ {
		offset[m] = true
	}
	for _, m := range append(append([]int{}, minterms...), dc...) {
		offset[m] = false
	}
	for _, m := range minterms {
		found := false
		for _, im := range cover {
			if im.Covers(m) {
				found = true
				break
			}
		}
		if !found {
			t.Errorf("minterm %d not covered", m)
		}
	}
	for m, off := range offset {
		if !off {
			continue
		}
		for _, im := range cover {
			if im.Covers(m) {
				t.Errorf("off-set minterm %d wrongly covered by %s", m, im)
			}
		}
	}
}

func TestKarnaughMap(t *testing.T) {
	k := NewKarnaughMap([]int{0, 2, 4, 6}, nil, 3)
	if k.RowVars != 1 || k.ColVars != 2 {
		t.Errorf("unexpected K-map shape: rows=%d cols=%d", k.RowVars, k.ColVars)
	}
	// Count the ones; must equal number of on-set minterms.
	ones := 0
	for _, row := range k.Cells {
		for _, v := range row {
			if v == 1 {
				ones++
			}
		}
	}
	if ones != 4 {
		t.Errorf("expected 4 ones in K-map, got %d", ones)
	}
}

func TestGrayCode(t *testing.T) {
	g := GrayCode(3)
	want := []int{0, 1, 3, 2, 6, 7, 5, 4}
	if !reflect.DeepEqual(g, want) {
		t.Errorf("GrayCode(3)=%v want %v", g, want)
	}
	for i := 0; i < 16; i++ {
		if GrayDecode(GrayEncode(i)) != i {
			t.Errorf("gray round trip failed at %d", i)
		}
	}
}

func TestPopCountBitString(t *testing.T) {
	if PopCount(0b1011) != 3 {
		t.Errorf("PopCount(0b1011)=%d want 3", PopCount(0b1011))
	}
	if BitString(5, 4) != "0101" {
		t.Errorf("BitString(5,4)=%q want 0101", BitString(5, 4))
	}
}

func TestVarsSorted(t *testing.T) {
	e := MustParse("c & (a | b)")
	got := Vars(e)
	want := []string{"a", "b", "c"}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("Vars=%v want %v", got, want)
	}
	if Size(e) != 5 {
		t.Errorf("Size=%d want 5", Size(e))
	}
}

func TestSubstitute(t *testing.T) {
	e := MustParse("a & b")
	sub := Substitute(e, "a", MustParse("c | d"))
	if !Equivalent(sub, MustParse("(c | d) & b")) {
		t.Errorf("Substitute gave %s", sub)
	}
}

func TestAllSolutions(t *testing.T) {
	f := NewCNF(NewClause(PosLit(1), PosLit(2)))
	sols := AllSolutions(f)
	if CountSolutions(f) != 3 {
		t.Errorf("expected 3 solutions, got %d", CountSolutions(f))
	}
	for _, s := range sols {
		if !f.Eval(s) {
			t.Errorf("enumerated non-solution %v", s)
		}
	}
}

func TestDNFNegate(t *testing.T) {
	d := NewDNF(NewClause(PosLit(1), PosLit(2)), NewClause(NegLit(3)))
	cnf := d.Negate()
	vars := []int{1, 2, 3}
	for mask := 0; mask < 8; mask++ {
		a := maskToAssign(vars, mask)
		if d.Eval(a) == cnf.Eval(a) {
			t.Errorf("DNF and its negation agree at %v", a)
		}
	}
}

func TestVarMap(t *testing.T) {
	vm := NewVarMap([]string{"x", "y"})
	if i, _ := vm.Index("x"); i != 1 {
		t.Errorf("x index = %d want 1", i)
	}
	if vm.Add("z") != 3 {
		t.Errorf("adding z should give index 3")
	}
	if vm.Name(2) != "y" {
		t.Errorf("Name(2)=%q want y", vm.Name(2))
	}
}

// helper to keep sort imported meaningfully in assertions
func sortedInts(xs []int) []int {
	out := append([]int(nil), xs...)
	sort.Ints(out)
	return out
}

func TestPrimeImplicantsSorted(t *testing.T) {
	primes := PrimeImplicants([]int{1, 2}, nil, 2)
	if len(primes) != 2 {
		t.Fatalf("expected 2 primes for xor, got %d", len(primes))
	}
	m := sortedInts(primes[0].Minterms())
	if len(m) != 1 {
		t.Errorf("each xor prime covers 1 minterm, got %v", m)
	}
}

func ExampleIsTautology() {
	e := MustParse("a | ~a")
	fmt.Println(IsTautology(e))
	// Output: true
}

func ExampleSolveCNF() {
	// (x1 OR x2) AND (NOT x1 OR x3)
	f := NewCNF(
		NewClause(PosLit(1), PosLit(2)),
		NewClause(NegLit(1), PosLit(3)),
	)
	_, sat := SolveCNF(f)
	fmt.Println(sat)
	// Output: true
}

func ExampleMinimizeSOPExpr() {
	// Even minterms of three variables minimise to ~c.
	expr := MinimizeSOPExpr([]int{0, 2, 4, 6}, nil, []string{"a", "b", "c"})
	fmt.Println(expr)
	// Output: ~c
}

func ExampleParse() {
	e := MustParse("a -> b")
	fmt.Println(Equivalent(e, MustParse("~a | b")))
	// Output: true
}
