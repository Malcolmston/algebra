package logic

import "testing"

func TestCanonicalStrings(t *testing.T) {
	if got := DNFString(MustParse("A ^ B")); got != "!A&B | A&!B" {
		t.Errorf("DNFString(A^B) = %q", got)
	}
	if got := CNFString(MustParse("A ^ B")); got != "(A|B) & (!A|!B)" {
		t.Errorf("CNFString(A^B) = %q", got)
	}
	if got := DNFString(MustParse("A & !A")); got != "F" {
		t.Errorf("DNFString(contradiction) = %q, want F", got)
	}
	if got := CNFString(MustParse("A | !A")); got != "T" {
		t.Errorf("CNFString(tautology) = %q, want T", got)
	}
}

func TestNormalFormsEquivalent(t *testing.T) {
	srcs := []string{
		"A ^ B",
		"(A -> B) & (B -> C)",
		"A <-> (B | C)",
		"!(A & (B | !C))",
		"A nand (B nor C)",
	}
	for _, s := range srcs {
		e := MustParse(s)
		cnf := ToCNF(e)
		dnf := ToDNF(e)
		nnf := ToNNF(e)
		if !Equivalent(e, cnf) {
			t.Errorf("ToCNF(%q) not equivalent", s)
		}
		if !Equivalent(e, dnf) {
			t.Errorf("ToDNF(%q) not equivalent", s)
		}
		if !Equivalent(e, nnf) {
			t.Errorf("ToNNF(%q) not equivalent", s)
		}
		if !IsCNF(cnf) {
			t.Errorf("ToCNF(%q) not in CNF: %s", s, cnf.String())
		}
		if !IsDNF(dnf) {
			t.Errorf("ToDNF(%q) not in DNF: %s", s, dnf.String())
		}
		if !logicIsNNFForm(nnf) {
			t.Errorf("ToNNF(%q) not in NNF: %s", s, nnf.String())
		}
	}
}

// logicIsNNFForm reports whether e uses only And, Or and negations of variables
// or constants -- the shape guaranteed by ToNNF.
func logicIsNNFForm(e Expr) bool {
	switch t := e.(type) {
	case Var, Const:
		return true
	case *UnaryExpr:
		if t.Op != NotOp {
			return false
		}
		switch t.X.(type) {
		case Var, Const:
			return true
		default:
			return false
		}
	case *BinaryExpr:
		if t.Op != AndOp && t.Op != OrOp {
			return false
		}
		return logicIsNNFForm(t.L) && logicIsNNFForm(t.R)
	default:
		return false
	}
}

func TestIsCNFDNF(t *testing.T) {
	if !IsDNF(MustParse("A & B | C")) {
		t.Errorf("A&B|C should be DNF")
	}
	if IsCNF(MustParse("A & B | C")) {
		t.Errorf("A&B|C should not be CNF")
	}
	if !IsCNF(MustParse("(A | B) & C")) {
		t.Errorf("(A|B)&C should be CNF")
	}
	if IsDNF(MustParse("(A | B) & C")) {
		t.Errorf("(A|B)&C should not be DNF")
	}
	if !IsCNF(MustParse("A")) || !IsDNF(MustParse("A")) {
		t.Errorf("a single literal is both CNF and DNF")
	}
}

func TestSimplify(t *testing.T) {
	cases := []struct {
		src  string
		want string // expected String() of the simplified expression
	}{
		{"A & A", "A"},
		{"A | A", "A"},
		{"A | !A", "T"},
		{"A & !A", "F"},
		{"A -> A", "T"},
		{"A & T", "A"},
		{"A | F", "A"},
		{"A & F", "F"},
		{"A | T", "T"},
		{"!!A", "A"},
		{"A ^ A", "F"},
		{"A <-> A", "T"},
	}
	for _, c := range cases {
		got := Simplify(MustParse(c.src)).String()
		if got != c.want {
			t.Errorf("Simplify(%q) = %q, want %q", c.src, got, c.want)
		}
	}
}

func TestSimplifyPreservesMeaning(t *testing.T) {
	srcs := []string{
		"(A & B) | (A & B)",
		"A & (B | !B)",
		"(A -> B) & T",
		"!(A | F)",
		"A ^ (B ^ B)",
		"(A <-> A) & (C | !C)",
	}
	for _, s := range srcs {
		e := MustParse(s)
		if !Equivalent(e, Simplify(e)) {
			t.Errorf("Simplify(%q) changed meaning", s)
		}
	}
}
