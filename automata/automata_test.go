package automata

import (
	"fmt"
	"testing"
)

// evenA builds a 2-state DFA over {a,b} accepting strings with an even number of
// 'a's.
func evenA() *DFA {
	d := NewDFA(2, []rune{'a', 'b'}, 0)
	d.SetTransition(0, 'a', 1)
	d.SetTransition(1, 'a', 0)
	d.SetTransition(0, 'b', 0)
	d.SetTransition(1, 'b', 1)
	d.AddAccept(0)
	return d
}

// evenA4 builds a redundant 4-state version of evenA for minimisation tests.
func evenA4() *DFA {
	d := NewDFA(4, []rune{'a', 'b'}, 0)
	// 0,2 even (accept); 1,3 odd.
	d.SetTransition(0, 'a', 1)
	d.SetTransition(0, 'b', 2)
	d.SetTransition(1, 'a', 0)
	d.SetTransition(1, 'b', 3)
	d.SetTransition(2, 'a', 3)
	d.SetTransition(2, 'b', 0)
	d.SetTransition(3, 'a', 2)
	d.SetTransition(3, 'b', 1)
	d.AddAccept(0, 2)
	return d
}

// onlyAB builds a DFA accepting exactly the string "ab".
func onlyAB() *DFA {
	d := NewDFA(4, []rune{'a', 'b'}, 0)
	d.SetTransition(0, 'a', 1)
	d.SetTransition(1, 'b', 2)
	// state 3 is dead.
	d.AddAccept(2)
	return d
}

func TestDFAAccepts(t *testing.T) {
	d := evenA()
	cases := []struct {
		in   string
		want bool
	}{
		{"", true},
		{"a", false},
		{"aa", true},
		{"ab", false},
		{"aba", true},
		{"abab", true},
		{"bbbb", true},
		{"aabbaa", true},
	}
	for _, c := range cases {
		if got := d.Accepts(c.in); got != c.want {
			t.Errorf("Accepts(%q) = %v, want %v", c.in, got, c.want)
		}
	}
}

func TestDFAValidate(t *testing.T) {
	if err := evenA().Validate(); err != nil {
		t.Fatalf("unexpected validation error: %v", err)
	}
}

func TestDFACompleteAndComplement(t *testing.T) {
	d := onlyAB()
	if d.IsComplete() {
		t.Fatal("onlyAB should be partial")
	}
	comp := d.Complement()
	// Complement of {ab} accepts everything except "ab".
	for _, s := range []string{"", "a", "b", "abb", "ba"} {
		if !comp.Accepts(s) {
			t.Errorf("complement should accept %q", s)
		}
	}
	if comp.Accepts("ab") {
		t.Error("complement should reject ab")
	}
}

func TestShortestAndCount(t *testing.T) {
	d := evenA()
	s, ok := d.ShortestAccepted()
	if !ok || s != "" {
		t.Fatalf("ShortestAccepted = %q,%v want '',true", s, ok)
	}
	if got := d.CountAcceptedWords(0); got.Int64() != 1 {
		t.Errorf("count len 0 = %v want 1", got)
	}
	if got := d.CountAcceptedWords(2); got.Int64() != 2 {
		t.Errorf("count len 2 = %v want 2", got)
	}
	if got := d.CountAcceptedWords(3); got.Int64() != 4 {
		// length-3 words over {a,b} with even a's: choose positions of a's in
		// even count: 0 a's (bbb) =1, 2 a's C(3,2)=3 => 4.
		t.Errorf("count len 3 = %v want 4", got)
	}
}

func TestFiniteAndUniversal(t *testing.T) {
	if evenA().IsFiniteLanguage() {
		t.Error("evenA language is infinite")
	}
	if !onlyAB().IsFiniteLanguage() {
		t.Error("{ab} is finite")
	}
	if evenA().IsEmptyLanguage() {
		t.Error("evenA is not empty")
	}
	// Universal DFA over {a}.
	u := NewDFA(1, []rune{'a'}, 0)
	u.SetTransition(0, 'a', 0)
	u.AddAccept(0)
	if !u.IsUniversal() {
		t.Error("u should be universal")
	}
}

func TestEpsilonClosure(t *testing.T) {
	n := NewNFA(3, []rune{'a'}, 0)
	n.AddEpsilon(0, 1)
	n.AddEpsilon(1, 2)
	n.AddTransition(0, 'a', 0)
	cl := n.EpsilonClosureState(0)
	for _, q := range []int{0, 1, 2} {
		if !cl[q] {
			t.Errorf("closure missing %d", q)
		}
	}
	if len(cl) != 3 {
		t.Errorf("closure size = %d want 3", len(cl))
	}
}

func TestSubsetConstruction(t *testing.T) {
	// NFA for strings over {a,b} that end in "ab".
	n := NewNFA(3, []rune{'a', 'b'}, 0)
	n.AddTransition(0, 'a', 0)
	n.AddTransition(0, 'b', 0)
	n.AddTransition(0, 'a', 1)
	n.AddTransition(1, 'b', 2)
	n.AddAccept(2)
	d := n.ToDFA()
	cases := map[string]bool{
		"ab": true, "aab": true, "bab": true, "abab": true,
		"a": false, "b": false, "ba": false, "abb": false, "": false,
	}
	for s, want := range cases {
		if got := d.Accepts(s); got != want {
			t.Errorf("DFA.Accepts(%q) = %v want %v", s, got, want)
		}
		if got := n.Accepts(s); got != want {
			t.Errorf("NFA.Accepts(%q) = %v want %v", s, got, want)
		}
	}
}

func TestMinimize(t *testing.T) {
	d := evenA4()
	m := d.Minimize()
	if m.NumStates != 2 {
		t.Errorf("Hopcroft states = %d want 2", m.NumStates)
	}
	mm := d.MinimizeMoore()
	if mm.NumStates != 2 {
		t.Errorf("Moore states = %d want 2", mm.NumStates)
	}
	if !Equivalent(d, m) {
		t.Error("minimised DFA not equivalent to original")
	}
	if !Equivalent(m, mm) {
		t.Error("Hopcroft and Moore disagree")
	}
	classes := d.EquivalenceClasses()
	if len(classes) != 2 {
		t.Errorf("equivalence classes = %d want 2", len(classes))
	}
}

func TestProductOps(t *testing.T) {
	even := evenA()
	// endsA: strings over {a,b} ending in 'a'.
	endsA := NewDFA(2, []rune{'a', 'b'}, 0)
	endsA.SetTransition(0, 'a', 1)
	endsA.SetTransition(0, 'b', 0)
	endsA.SetTransition(1, 'a', 1)
	endsA.SetTransition(1, 'b', 0)
	endsA.AddAccept(1)

	un := Union(even, endsA)
	in := Intersection(even, endsA)
	di := Difference(even, endsA)

	check := func(name string, d *DFA, s string, want bool) {
		if got := d.Accepts(s); got != want {
			t.Errorf("%s.Accepts(%q) = %v want %v", name, s, got, want)
		}
	}
	// "aa": even a's (2) true, ends in a true.
	check("union", un, "aa", true)
	check("inter", in, "aa", true)
	check("diff", di, "aa", false)
	// "a": even false, ends a true.
	check("union", un, "a", true)
	check("inter", in, "a", false)
	check("diff", di, "a", false)
	// "bb": even true, ends a false.
	check("union", un, "bb", true)
	check("inter", in, "bb", false)
	check("diff", di, "bb", true)
}

func TestEquivalenceAndWitness(t *testing.T) {
	even := evenA()
	if !even.IsEquivalentTo(even.Minimize()) {
		t.Error("DFA not equivalent to its minimisation")
	}
	// Universal over {a,b}.
	uni := NewDFA(1, []rune{'a', 'b'}, 0)
	uni.SetTransition(0, 'a', 0)
	uni.SetTransition(0, 'b', 0)
	uni.AddAccept(0)
	w, ok := Witness(even, uni)
	if !ok {
		t.Fatal("expected a witness distinguishing evenA from universal")
	}
	// Shortest distinguishing string is "a" (odd a's).
	if w != "a" {
		t.Errorf("witness = %q want %q", w, "a")
	}
	if Equivalent(even, uni) {
		t.Error("evenA and universal must not be equivalent")
	}
	if !Subset(even, uni) {
		t.Error("evenA should be a subset of universal")
	}
}

func TestPumping(t *testing.T) {
	d := evenA()
	p := PumpingLength(d)
	if p != 2 {
		t.Errorf("pumping length = %d want 2", p)
	}
	decomp := PumpingDecomposition(d)
	x, y, z, ok := decomp("aabb")
	if !ok {
		t.Fatal("expected a decomposition of aabb")
	}
	if x+y+z != "aabb" {
		t.Errorf("xyz = %q want aabb", x+y+z)
	}
	if len(y) == 0 {
		t.Error("|y| must be >= 1")
	}
	// Pump: x y^i z must remain accepted for several i.
	for _, i := range []int{0, 1, 2, 3} {
		w := x
		for k := 0; k < i; k++ {
			w += y
		}
		w += z
		if !d.Accepts(w) {
			t.Errorf("pumped word %q (i=%d) should be accepted", w, i)
		}
	}
}

func TestRegexParseAndMatch(t *testing.T) {
	re := MustCompile("a(b|c)*")
	cases := map[string]bool{
		"a": true, "ab": true, "ac": true, "abcbcb": true,
		"": false, "b": false, "ba": false, "bc": false,
	}
	for s, want := range cases {
		if got := re.Matches(s); got != want {
			t.Errorf("Matches(%q) = %v want %v", s, got, want)
		}
	}
}

func TestRegexAlternationConcatStar(t *testing.T) {
	// (ab|ba)* over {a,b}.
	re := MustCompile("(ab|ba)*")
	cases := map[string]bool{
		"": true, "ab": true, "ba": true, "abba": true, "abab": true, "baab": true,
		"a": false, "b": false, "aba": false, "abb": false,
	}
	for s, want := range cases {
		if got := re.Matches(s); got != want {
			t.Errorf("Matches(%q) = %v want %v", s, got, want)
		}
	}
}

func TestRegexPlusOptionalDot(t *testing.T) {
	re := MustCompile("a+b?")
	cases := map[string]bool{
		"a": true, "aa": true, "aaab": true, "ab": true,
		"": false, "b": false, "abb": false, "ba": false,
	}
	for s, want := range cases {
		if got := re.Matches(s); got != want {
			t.Errorf("Matches(%q) = %v want %v", s, got, want)
		}
	}
	// Wildcard over the alphabet {a,b}.
	dot, err := CompileWithAlphabet("a.b", []rune{'a', 'b'})
	if err != nil {
		t.Fatal(err)
	}
	for _, s := range []string{"aab", "abb"} {
		if !dot.Matches(s) {
			t.Errorf("dot should match %q", s)
		}
	}
	if dot.Matches("ab") {
		t.Error("dot pattern a.b should not match ab")
	}
}

func TestRegexError(t *testing.T) {
	if _, err := ParseRegex("a("); err == nil {
		t.Error("expected error for unbalanced parenthesis")
	}
	if _, err := ParseRegex("*a"); err == nil {
		t.Error("expected error for leading quantifier")
	}
}

func TestLangOps(t *testing.T) {
	a := LiteralNFA("ab")
	b := LiteralNFA("cd")
	con := Concatenate(a, b)
	if !con.Accepts("abcd") || con.Accepts("ab") || con.Accepts("abc") {
		t.Error("concatenation failed")
	}
	un := UnionNFA(a, b)
	if !un.Accepts("ab") || !un.Accepts("cd") || un.Accepts("abcd") {
		t.Error("union failed")
	}
	st := StarNFA(LiteralNFA("ab"))
	for _, s := range []string{"", "ab", "abab", "ababab"} {
		if !st.Accepts(s) {
			t.Errorf("star should accept %q", s)
		}
	}
	if st.Accepts("aba") {
		t.Error("star should reject aba")
	}
	pl := PlusNFA(LiteralNFA("x"))
	if pl.Accepts("") || !pl.Accepts("x") || !pl.Accepts("xxx") {
		t.Error("plus failed")
	}
	op := OptionalNFA(LiteralNFA("z"))
	if !op.Accepts("") || !op.Accepts("z") || op.Accepts("zz") {
		t.Error("optional failed")
	}
	pw := PowerNFA(LiteralNFA("ab"), 3)
	if !pw.Accepts("ababab") || pw.Accepts("abab") {
		t.Error("power failed")
	}
}

func TestReversal(t *testing.T) {
	// Language {ab}; reversal is {ba}.
	d := onlyAB()
	r := d.Reverse()
	if !r.Accepts("ba") || r.Accepts("ab") {
		t.Error("reversal of {ab} should be {ba}")
	}
}

func TestDFAToRegexRoundTrip(t *testing.T) {
	d := evenA()
	node := DFAToRegex(d)
	// Recompile the regex over the same alphabet and check equivalence on a
	// battery of strings.
	re, err := CompileWithAlphabet(node.String(), []rune{'a', 'b'})
	if err != nil {
		t.Fatalf("recompile %q: %v", node.String(), err)
	}
	for _, s := range []string{"", "a", "b", "aa", "ab", "ba", "bb", "aab", "abab", "aaa", "bbba"} {
		if re.Matches(s) != d.Accepts(s) {
			t.Errorf("round-trip mismatch on %q: regex=%v dfa=%v (regex=%s)",
				s, re.Matches(s), d.Accepts(s), node.String())
		}
	}
}

func TestNFAEpsilonElimination(t *testing.T) {
	// a?b built with epsilon edges, then eliminate epsilons.
	n := RegexToNFAMust("a?b")
	e := n.RemoveEpsilon()
	if e.HasEpsilon() {
		t.Error("epsilon should be eliminated")
	}
	for _, s := range []string{"b", "ab"} {
		if !e.Accepts(s) {
			t.Errorf("should accept %q", s)
		}
	}
	if e.Accepts("aab") || e.Accepts("") {
		t.Error("should reject aab and empty")
	}
}

// RegexToNFAMust is a small test helper.
func RegexToNFAMust(p string) *NFA {
	n, err := RegexToNFA(p)
	if err != nil {
		panic(err)
	}
	return n
}

func TestTuringMachineAnBn(t *testing.T) {
	// TM deciding {0^n 1^n : n >= 0}.
	m := NewTM("q0", "acc", "rej", '_')
	// q0: at left; cross a 0 with X.
	m.SetRule("q0", '0', "q1", 'X', Right)
	m.SetRule("q0", 'Y', "q3", 'Y', Right)
	m.SetRule("q0", '_', "acc", '_', Stay)
	// q1: scan right for a 1.
	m.SetRule("q1", '0', "q1", '0', Right)
	m.SetRule("q1", 'Y', "q1", 'Y', Right)
	m.SetRule("q1", '1', "q2", 'Y', Left)
	// q2: scan left back to X.
	m.SetRule("q2", '0', "q2", '0', Left)
	m.SetRule("q2", 'Y', "q2", 'Y', Left)
	m.SetRule("q2", 'X', "q0", 'X', Right)
	// q3: only Y's then blank.
	m.SetRule("q3", 'Y', "q3", 'Y', Right)
	m.SetRule("q3", '_', "acc", '_', Stay)

	cases := map[string]bool{
		"":       true,
		"01":     true,
		"0011":   true,
		"000111": true,
		"1":      false,
		"0":      false,
		"10":     false,
		"0111":   false,
		"0101":   false,
		"001":    false,
	}
	for in, want := range cases {
		got, err := m.Accepts(in, 10000)
		if err != nil {
			t.Fatalf("TM did not halt on %q: %v", in, err)
		}
		if got != want {
			t.Errorf("TM.Accepts(%q) = %v want %v", in, got, want)
		}
	}
}

func TestTuringMachineTransducer(t *testing.T) {
	// TM that flips a single leading bit: 0 -> 1 then accept.
	m := NewTM("s", "acc", "rej", '_')
	m.SetRule("s", '0', "acc", '1', Stay)
	m.SetRule("s", '1', "acc", '0', Stay)
	out, err := m.Compute("0", 100)
	if err != nil {
		t.Fatal(err)
	}
	if out != "1" {
		t.Errorf("Compute(0) = %q want 1", out)
	}
}

func TestPDAAnBn(t *testing.T) {
	// PDA for {a^n b^n : n >= 0}, accept by final state.
	p := NewPDA("q0", 'Z', AcceptByFinalState)
	p.AddTransition("q0", 'a', 'Z', "q0", 'A', 'Z')
	p.AddTransition("q0", 'a', 'A', "q0", 'A', 'A')
	p.AddEpsilonTransition("q0", 'Z', "q2")
	p.AddTransition("q0", 'b', 'A', "q1")
	p.AddTransition("q1", 'b', 'A', "q1")
	p.AddEpsilonTransition("q1", 'Z', "q2")
	p.AddAccept("q2")

	cases := map[string]bool{
		"":       true,
		"ab":     true,
		"aabb":   true,
		"aaabbb": true,
		"a":      false,
		"b":      false,
		"aab":    false,
		"abb":    false,
		"ba":     false,
		"aabbb":  false,
	}
	for in, want := range cases {
		if got := p.Accepts(in); got != want {
			t.Errorf("PDA.Accepts(%q) = %v want %v", in, got, want)
		}
	}
}

func TestPDAEmptyStack(t *testing.T) {
	// Balanced parentheses over {(,)} accepting by empty stack.
	p := NewPDA("q", 'Z', AcceptByEmptyStack)
	p.AddTransition("q", '(', 'Z', "q", 'O', 'Z')
	p.AddTransition("q", '(', 'O', "q", 'O', 'O')
	p.AddTransition("q", ')', 'O', "q")
	p.AddEpsilonTransition("q", 'Z', "q")
	cases := map[string]bool{
		"":     true,
		"()":   true,
		"(())": true,
		"()()": true,
		"(":    false,
		")":    false,
		"(()":  false,
		"())":  false,
		")(":   false,
	}
	for in, want := range cases {
		if got := p.Accepts(in); got != want {
			t.Errorf("PDA.Accepts(%q) = %v want %v", in, got, want)
		}
	}
}

func TestStateSet(t *testing.T) {
	s := NewStateSet(3, 1, 2)
	if s.Key() != "1,2,3" {
		t.Errorf("Key = %q want 1,2,3", s.Key())
	}
	if !s.Equal(NewStateSet(1, 2, 3)) {
		t.Error("Equal failed")
	}
	u := s.Union(NewStateSet(4))
	if u.Len() != 4 || !u.Contains(4) {
		t.Error("Union failed")
	}
	in := s.Intersect(NewStateSet(2, 5))
	if in.Len() != 1 || !in.Contains(2) {
		t.Error("Intersect failed")
	}
}

func ExampleCompile() {
	re := MustCompile("a(b|c)*")
	fmt.Println(re.Matches("abcbc"))
	fmt.Println(re.Matches("a"))
	fmt.Println(re.Matches("b"))
	// Output:
	// true
	// true
	// false
}

func ExampleNFA_ToDFA() {
	// NFA for strings over {a,b} ending in "ab".
	n := NewNFA(3, []rune{'a', 'b'}, 0)
	n.AddTransition(0, 'a', 0)
	n.AddTransition(0, 'b', 0)
	n.AddTransition(0, 'a', 1)
	n.AddTransition(1, 'b', 2)
	n.AddAccept(2)
	d := n.ToDFA().Minimize()
	fmt.Println(d.Accepts("aab"), d.Accepts("abb"))
	// Output:
	// true false
}

func ExampleDFA_CountAcceptedWords() {
	fmt.Println(evenA().CountAcceptedWords(3))
	// Output:
	// 4
}
