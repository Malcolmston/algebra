package logic

import (
	"reflect"
	"testing"
)

func TestTautologyContradiction(t *testing.T) {
	cases := []struct {
		src                        string
		taut, contra, sat, conting bool
	}{
		{"A | !A", true, false, true, false},
		{"A & !A", false, true, false, false},
		{"A -> A", true, false, true, false},
		{"A -> B", false, false, true, true},
		{"(A -> B) & (B -> A) <-> (A <-> B)", true, false, true, false},
		{"A & B", false, false, true, true},
	}
	for _, c := range cases {
		e := MustParse(c.src)
		if got := IsTautology(e); got != c.taut {
			t.Errorf("IsTautology(%q) = %v, want %v", c.src, got, c.taut)
		}
		if got := IsContradiction(e); got != c.contra {
			t.Errorf("IsContradiction(%q) = %v, want %v", c.src, got, c.contra)
		}
		if got := IsSatisfiable(e); got != c.sat {
			t.Errorf("IsSatisfiable(%q) = %v, want %v", c.src, got, c.sat)
		}
		if got := IsContingency(e); got != c.conting {
			t.Errorf("IsContingency(%q) = %v, want %v", c.src, got, c.conting)
		}
	}
}

func TestModels(t *testing.T) {
	e := MustParse("A & (B | C)")
	// Satisfying assignments over sorted vars A,B,C:
	// A=1,B=0,C=1 (5); A=1,B=1,C=0 (6); A=1,B=1,C=1 (7).
	if n := CountModels(e); n != 3 {
		t.Errorf("CountModels = %d, want 3", n)
	}
	m, ok := FindModel(e)
	if !ok {
		t.Fatalf("FindModel found nothing")
	}
	if r, _ := e.Eval(m); !r {
		t.Errorf("FindModel returned non-model %v", m)
	}
	// Deterministic first model is index 5: A=1,B=0,C=1.
	wantFirst := map[string]bool{"A": true, "B": false, "C": true}
	if !reflect.DeepEqual(m, wantFirst) {
		t.Errorf("FindModel = %v, want %v", m, wantFirst)
	}
	if got := len(AllModels(e)); got != 3 {
		t.Errorf("AllModels len = %d, want 3", got)
	}
}

func TestUnsatModels(t *testing.T) {
	e := MustParse("A & !A")
	if _, ok := FindModel(e); ok {
		t.Errorf("unsat expr should have no model")
	}
	if CountModels(e) != 0 {
		t.Errorf("unsat expr should have 0 models")
	}
}

func TestEquivalentDeMorgan(t *testing.T) {
	a := MustParse("!(A & B)")
	b := MustParse("!A | !B")
	if !Equivalent(a, b) {
		t.Errorf("De Morgan equivalence failed")
	}
	if Equivalent(MustParse("A | B"), MustParse("A & B")) {
		t.Errorf("A|B should not equal A&B")
	}
}

func TestEntails(t *testing.T) {
	if !Entails(MustParse("A & B"), MustParse("A")) {
		t.Errorf("A & B should entail A")
	}
	if Entails(MustParse("A"), MustParse("A & B")) {
		t.Errorf("A should not entail A & B")
	}
	// Modus ponens: A and (A->B) entail B.
	if !Entails(MustParse("A & (A -> B)"), MustParse("B")) {
		t.Errorf("modus ponens entailment failed")
	}
}
