package settheory

import (
	"reflect"
	"testing"
)

func TestIntSetOps(t *testing.T) {
	a := NewIntSet(1, 2, 3, 4)
	b := NewIntSet(3, 4, 5, 6)

	if got := a.Union(b).Elements(); !reflect.DeepEqual(got, []int{1, 2, 3, 4, 5, 6}) {
		t.Errorf("Union = %v", got)
	}
	if got := a.Intersection(b).Elements(); !reflect.DeepEqual(got, []int{3, 4}) {
		t.Errorf("Intersection = %v", got)
	}
	if got := a.Difference(b).Elements(); !reflect.DeepEqual(got, []int{1, 2}) {
		t.Errorf("Difference = %v", got)
	}
	if got := a.SymmetricDifference(b).Elements(); !reflect.DeepEqual(got, []int{1, 2, 5, 6}) {
		t.Errorf("SymmetricDifference = %v", got)
	}
	// Inclusion-exclusion: |A ∪ B| = |A| + |B| - |A ∩ B|.
	if a.Union(b).Len() != a.Len()+b.Len()-a.Intersection(b).Len() {
		t.Errorf("inclusion-exclusion violated")
	}
}

func TestIntSetPredicates(t *testing.T) {
	full := NewIntSet(1, 2, 3)
	sub := NewIntSet(1, 2)
	if !sub.IsSubsetOf(full) || !sub.IsProperSubsetOf(full) {
		t.Errorf("subset predicates failed")
	}
	if !full.IsSupersetOf(sub) {
		t.Errorf("superset failed")
	}
	if full.IsProperSubsetOf(full) {
		t.Errorf("a set is not a proper subset of itself")
	}
	if !NewIntSet(1, 2).IsDisjoint(NewIntSet(3, 4)) {
		t.Errorf("disjoint failed")
	}
	if NewIntSet(1, 2).IsDisjoint(NewIntSet(2, 3)) {
		t.Errorf("non-disjoint reported disjoint")
	}
	empty := NewIntSet()
	if !empty.IsSubsetOf(full) {
		t.Errorf("empty set must be a subset of every set")
	}
}

func TestIntSetMinMaxSum(t *testing.T) {
	s := NewIntSet(4, -2, 7, 0)
	if mn, ok := s.Min(); !ok || mn != -2 {
		t.Errorf("Min = %d,%v", mn, ok)
	}
	if mx, ok := s.Max(); !ok || mx != 7 {
		t.Errorf("Max = %d,%v", mx, ok)
	}
	if s.Sum() != 9 {
		t.Errorf("Sum = %d", s.Sum())
	}
	if _, ok := NewIntSet().Min(); ok {
		t.Errorf("Min of empty set should report !ok")
	}
}

func TestPowerSetCardinality(t *testing.T) {
	for n := 0; n <= 8; n++ {
		elems := make([]int, n)
		for i := range elems {
			elems[i] = i
		}
		ps := NewIntSet(elems...).PowerSet()
		if len(ps) != 1<<uint(n) {
			t.Errorf("|P(S)| for n=%d = %d, want %d", n, len(ps), 1<<uint(n))
		}
	}
	// The full set is the last subset produced.
	ps := NewIntSet(1, 2, 3).PowerSet()
	if !ps[len(ps)-1].Equal(NewIntSet(1, 2, 3)) {
		t.Errorf("last power-set entry should be the full set")
	}
	if !ps[0].IsEmpty() {
		t.Errorf("first power-set entry should be empty")
	}
}

func TestCartesianProduct(t *testing.T) {
	a := NewIntSet(1, 2)
	b := NewIntSet(3, 4)
	got := a.CartesianProduct(b)
	want := []Pair{{1, 3}, {1, 4}, {2, 3}, {2, 4}}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("CartesianProduct = %v, want %v", got, want)
	}
}

func TestIntSetString(t *testing.T) {
	if got := NewIntSet(3, 1, 2).String(); got != "{1, 2, 3}" {
		t.Errorf("String = %q", got)
	}
	if got := NewIntSet().String(); got != "{}" {
		t.Errorf("empty String = %q", got)
	}
}

func TestStringSetOps(t *testing.T) {
	a := NewStringSet("a", "b", "c")
	b := NewStringSet("b", "c", "d")
	if got := a.Union(b).Elements(); !reflect.DeepEqual(got, []string{"a", "b", "c", "d"}) {
		t.Errorf("Union = %v", got)
	}
	if got := a.Intersection(b).Elements(); !reflect.DeepEqual(got, []string{"b", "c"}) {
		t.Errorf("Intersection = %v", got)
	}
	if got := a.SymmetricDifference(b).Elements(); !reflect.DeepEqual(got, []string{"a", "d"}) {
		t.Errorf("SymmetricDifference = %v", got)
	}
	if len(a.PowerSet()) != 8 {
		t.Errorf("power set size = %d", len(a.PowerSet()))
	}
	if got := a.String(); got != `{"a", "b", "c"}` {
		t.Errorf("String = %q", got)
	}
}
