package logic

import (
	"reflect"
	"testing"
)

func TestParseEval(t *testing.T) {
	cases := []struct {
		src  string
		env  map[string]bool
		want bool
	}{
		{"A & B", map[string]bool{"A": true, "B": true}, true},
		{"A & B", map[string]bool{"A": true, "B": false}, false},
		{"A | B", map[string]bool{"A": false, "B": false}, false},
		{"!A", map[string]bool{"A": false}, true},
		{"A -> B", map[string]bool{"A": true, "B": false}, false},
		{"A -> B", map[string]bool{"A": false, "B": false}, true},
		{"A <-> B", map[string]bool{"A": true, "B": true}, true},
		{"A ^ B", map[string]bool{"A": true, "B": true}, false},
		{"A nand B", map[string]bool{"A": true, "B": true}, false},
		{"A nor B", map[string]bool{"A": false, "B": false}, true},
		{"A xnor B", map[string]bool{"A": true, "B": true}, true},
		{"T", nil, true},
		{"F", nil, false},
		{"1 & 0", nil, false},
		{"A and (B or not C)", map[string]bool{"A": true, "B": false, "C": false}, true},
	}
	for _, c := range cases {
		got, err := EvalString(c.src, c.env)
		if err != nil {
			t.Errorf("EvalString(%q) error: %v", c.src, err)
			continue
		}
		if got != c.want {
			t.Errorf("EvalString(%q) = %v, want %v", c.src, got, c.want)
		}
	}
}

func TestParsePrecedence(t *testing.T) {
	// & binds tighter than |, which binds tighter than ->, which binds
	// tighter than <->. Verify by comparing truth tables to a fully
	// parenthesised reference.
	pairs := [][2]string{
		{"A | B & C", "A | (B & C)"},
		{"A -> B | C", "A -> (B | C)"},
		{"A <-> B -> C", "A <-> (B -> C)"},
		{"!A & B", "(!A) & B"},
		{"A ^ B & C", "A ^ (B & C)"},
	}
	for _, p := range pairs {
		e1 := MustParse(p[0])
		e2 := MustParse(p[1])
		if !Equivalent(e1, e2) {
			t.Errorf("%q not equivalent to %q", p[0], p[1])
		}
	}
}

func TestParseRoundTrip(t *testing.T) {
	// The String rendering must re-parse to an equivalent expression.
	srcs := []string{"A & B | C", "!(A -> B)", "A <-> (B ^ C)", "A nand B nor C"}
	for _, s := range srcs {
		e := MustParse(s)
		re := MustParse(e.String())
		if !Equivalent(e, re) {
			t.Errorf("round trip failed for %q -> %q", s, e.String())
		}
	}
}

func TestParseErrors(t *testing.T) {
	bad := []string{"A &", "(A", "A B", "@", ""}
	for _, s := range bad {
		if _, err := Parse(s); err == nil {
			t.Errorf("Parse(%q) expected error", s)
		}
	}
}

func TestUnboundVariable(t *testing.T) {
	e := MustParse("A & B")
	if _, err := e.Eval(map[string]bool{"A": true}); err == nil {
		t.Errorf("expected unbound-variable error")
	}
}

func TestVars(t *testing.T) {
	e := MustParse("C | A & B | A")
	got := Vars(e)
	want := []string{"A", "B", "C"}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("Vars = %v, want %v", got, want)
	}
}

func TestConstructors(t *testing.T) {
	e := NewIff(NewAnd(NewVar("A"), NewVar("B")), NewNot(NewOr(NewNot(NewVar("A")), NewNot(NewVar("B")))))
	// (A & B) <-> !(!A | !B) is a tautology (De Morgan).
	if !IsTautology(e) {
		t.Errorf("De Morgan tautology failed: %s", e.String())
	}
	if !Equivalent(NewConst(true), True) {
		t.Errorf("True constant mismatch")
	}
}
