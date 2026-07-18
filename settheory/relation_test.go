package settheory

import (
	"reflect"
	"testing"
)

func TestRelationProperties(t *testing.T) {
	set := NewIntSet(1, 2, 3)
	// Equality on {1,2,3}: reflexive, symmetric, transitive.
	eq := RelationFromPairs([][2]int{{1, 1}, {2, 2}, {3, 3}})
	if !eq.IsReflexiveOn(set) || !eq.IsSymmetric() || !eq.IsTransitive() {
		t.Errorf("identity relation properties wrong")
	}
	if !eq.IsEquivalenceOn(set) {
		t.Errorf("identity should be an equivalence relation")
	}
	// A strict-less relation is antisymmetric, transitive, irreflexive.
	lt := RelationFromPairs([][2]int{{1, 2}, {2, 3}, {1, 3}})
	if !lt.IsAntisymmetric() || !lt.IsTransitive() || !lt.IsIrreflexiveOn(set) {
		t.Errorf("strict order properties wrong")
	}
	if lt.IsSymmetric() {
		t.Errorf("strict order must not be symmetric")
	}
}

func TestTransitiveClosureKnown(t *testing.T) {
	// Path 1->2->3->4; closure is the full reachability order.
	r := RelationFromPairs([][2]int{{1, 2}, {2, 3}, {3, 4}})
	got := r.TransitiveClosure().Pairs()
	want := []Pair{{1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4}}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("TransitiveClosure = %v, want %v", got, want)
	}
}

func TestSymmetricAndReflexiveClosure(t *testing.T) {
	r := RelationFromPairs([][2]int{{1, 2}})
	sc := r.SymmetricClosure()
	if !sc.Related(2, 1) || !sc.Related(1, 2) {
		t.Errorf("symmetric closure missing a pair")
	}
	rc := r.ReflexiveClosureOn(NewIntSet(1, 2, 3))
	for _, x := range []int{1, 2, 3} {
		if !rc.Related(x, x) {
			t.Errorf("reflexive closure missing (%d,%d)", x, x)
		}
	}
}

func TestEquivalenceClosure(t *testing.T) {
	set := NewIntSet(1, 2, 3, 4)
	// 1~2 and 3~4 induce two classes {1,2} and {3,4}.
	r := RelationFromPairs([][2]int{{1, 2}, {3, 4}})
	ec := r.EquivalenceClosureOn(set)
	if !ec.IsEquivalenceOn(set) {
		t.Fatalf("closure is not an equivalence relation")
	}
	if !ec.Related(2, 1) || !ec.Related(1, 1) || ec.Related(1, 3) {
		t.Errorf("equivalence closure classes wrong")
	}
	classes := EquivalenceClasses(ec, set)
	if classes.Len() != 2 {
		t.Errorf("expected 2 classes, got %d", classes.Len())
	}
}

func TestInverseAndCompose(t *testing.T) {
	r := RelationFromPairs([][2]int{{1, 2}, {2, 3}})
	if got := r.Inverse().Pairs(); !reflect.DeepEqual(got, []Pair{{2, 1}, {3, 2}}) {
		t.Errorf("Inverse = %v", got)
	}
	// r∘r: apply r then r again. 1->2->3 gives (1,3); 2->3->? none.
	comp := r.Compose(r)
	if !comp.Related(1, 3) || comp.Len() != 1 {
		t.Errorf("Compose = %v", comp.Pairs())
	}
}

func TestDomainRangeField(t *testing.T) {
	r := RelationFromPairs([][2]int{{1, 2}, {3, 4}})
	if got := r.Domain().Elements(); !reflect.DeepEqual(got, []int{1, 3}) {
		t.Errorf("Domain = %v", got)
	}
	if got := r.Range().Elements(); !reflect.DeepEqual(got, []int{2, 4}) {
		t.Errorf("Range = %v", got)
	}
	if got := r.Field().Elements(); !reflect.DeepEqual(got, []int{1, 2, 3, 4}) {
		t.Errorf("Field = %v", got)
	}
}
