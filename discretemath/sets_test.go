package discretemath

import (
	"sort"
	"testing"
)

func sortedInts(s Set[int]) []int {
	out := s.Elements()
	sort.Ints(out)
	return out
}

func TestSetBasics(t *testing.T) {
	s := NewSet(1, 2, 2, 3)
	if s.Len() != 3 {
		t.Errorf("Len = %d, want 3", s.Len())
	}
	if !s.Contains(2) || s.Contains(9) {
		t.Error("Contains failed")
	}
	s.Add(4, 4)
	if s.Len() != 4 {
		t.Errorf("after Add Len = %d, want 4", s.Len())
	}
	s.Remove(1)
	if s.Contains(1) {
		t.Error("Remove failed")
	}
	c := s.Clone()
	c.Add(99)
	if s.Contains(99) {
		t.Error("Clone should be independent")
	}
}

func TestSetAlgebra(t *testing.T) {
	a := NewSet(1, 2, 3, 4)
	b := NewSet(3, 4, 5, 6)

	if got, want := sortedInts(Union(a, b)), []int{1, 2, 3, 4, 5, 6}; !equalInts(got, want) {
		t.Errorf("Union = %v, want %v", got, want)
	}
	if got, want := sortedInts(Intersection(a, b)), []int{3, 4}; !equalInts(got, want) {
		t.Errorf("Intersection = %v, want %v", got, want)
	}
	if got, want := sortedInts(Difference(a, b)), []int{1, 2}; !equalInts(got, want) {
		t.Errorf("Difference = %v, want %v", got, want)
	}
	if got, want := sortedInts(SymmetricDifference(a, b)), []int{1, 2, 5, 6}; !equalInts(got, want) {
		t.Errorf("SymmetricDifference = %v, want %v", got, want)
	}
}

func TestSetRelations(t *testing.T) {
	a := NewSet(1, 2)
	b := NewSet(1, 2, 3)
	if !a.IsSubset(b) || b.IsSubset(a) {
		t.Error("IsSubset failed")
	}
	if !b.IsSuperset(a) || a.IsSuperset(b) {
		t.Error("IsSuperset failed")
	}
	if !a.Equal(NewSet(2, 1)) || a.Equal(b) {
		t.Error("Equal failed")
	}
	if !NewSet(1, 2).IsDisjoint(NewSet(3, 4)) {
		t.Error("expected disjoint")
	}
	if NewSet(1, 2).IsDisjoint(NewSet(2, 3)) {
		t.Error("expected not disjoint")
	}
}

func TestPowerSet(t *testing.T) {
	ps := PowerSet([]int{1, 2, 3})
	if len(ps) != 8 { // 2^3
		t.Fatalf("PowerSet size = %d, want 8", len(ps))
	}
	// Check each subset length count matches binomial coefficients 1,3,3,1.
	counts := map[int]int{}
	for _, sub := range ps {
		counts[len(sub)]++
	}
	if counts[0] != 1 || counts[1] != 3 || counts[2] != 3 || counts[3] != 1 {
		t.Errorf("subset size distribution = %v", counts)
	}
	if len(PowerSet([]int{})) != 1 {
		t.Error("power set of empty set should have one element (the empty set)")
	}
}

func TestCartesianProduct(t *testing.T) {
	p := CartesianProduct([]int{1, 2}, []string{"a", "b"})
	if len(p) != 4 {
		t.Fatalf("size = %d, want 4", len(p))
	}
	want := []Pair[int, string]{
		{1, "a"}, {1, "b"}, {2, "a"}, {2, "b"},
	}
	for i, w := range want {
		if p[i] != w {
			t.Errorf("pair %d = %v, want %v", i, p[i], w)
		}
	}
}

func TestCartesianProductN(t *testing.T) {
	got := CartesianProductN([]int{0, 1}, []int{0, 1}, []int{0, 1})
	if len(got) != 8 {
		t.Fatalf("size = %d, want 8", len(got))
	}
	// Last input varies fastest: first tuple all zero, second is {0,0,1}.
	if !equalInts(got[0], []int{0, 0, 0}) || !equalInts(got[1], []int{0, 0, 1}) {
		t.Errorf("ordering wrong: %v, %v", got[0], got[1])
	}
	// Empty operand -> empty product.
	if len(CartesianProductN([]int{1, 2}, []int{})) != 0 {
		t.Error("expected empty product with empty operand")
	}
	// No operands -> single empty tuple.
	if len(CartesianProductN[int]()) != 1 {
		t.Error("expected single empty tuple with no operands")
	}
}

func equalInts(a, b []int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}
