package settheory

import (
	"reflect"
	"testing"
)

func TestDivisibilityLattice(t *testing.T) {
	p, err := NewDivisibilityPoset(12)
	if err != nil {
		t.Fatal(err)
	}
	if got := p.Elements(); !reflect.DeepEqual(got, []int{1, 2, 3, 4, 6, 12}) {
		t.Fatalf("divisors = %v", got)
	}
	// In a divisibility lattice meet = gcd and join = lcm.
	if m, ok := p.Meet(4, 6); !ok || m != 2 {
		t.Errorf("Meet(4,6) = %d,%v want 2", m, ok)
	}
	if j, ok := p.Join(4, 6); !ok || j != 12 {
		t.Errorf("Join(4,6) = %d,%v want 12", j, ok)
	}
	if m, ok := p.Meet(2, 3); !ok || m != 1 {
		t.Errorf("Meet(2,3) = %d,%v want 1", m, ok)
	}
	if !p.IsLattice() {
		t.Errorf("divisibility poset must be a lattice")
	}
	if l, ok := p.LeastElement(); !ok || l != 1 {
		t.Errorf("LeastElement = %d,%v want 1", l, ok)
	}
	if g, ok := p.GreatestElement(); !ok || g != 12 {
		t.Errorf("GreatestElement = %d,%v want 12", g, ok)
	}
	if p.Height() != 4 {
		t.Errorf("Height = %d, want 4 (chain 1|2|4|12)", p.Height())
	}
}

func TestHasseEdgesDivisibility12(t *testing.T) {
	p, _ := NewDivisibilityPoset(12)
	got := p.HasseEdges()
	want := []Pair{{1, 2}, {1, 3}, {2, 4}, {2, 6}, {3, 6}, {4, 12}, {6, 12}}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("HasseEdges = %v, want %v", got, want)
	}
}

func TestTopologicalOrderConsistency(t *testing.T) {
	p, _ := NewDivisibilityPoset(36)
	order := p.TopologicalOrder()
	idx := make(map[int]int, len(order))
	for i, x := range order {
		idx[x] = i
	}
	// Every strict order relation must be respected by the linear extension.
	for _, a := range p.Elements() {
		for _, b := range p.Elements() {
			if p.Less(a, b) && idx[a] >= idx[b] {
				t.Errorf("topological order violates %d < %d", a, b)
			}
		}
	}
	if len(order) != p.Size() {
		t.Errorf("order length %d != size %d", len(order), p.Size())
	}
}

func TestChainPoset(t *testing.T) {
	p := NewChainPoset(5, 1, 3)
	if !p.IsChain() {
		t.Errorf("chain poset should be a chain")
	}
	if !p.IsLattice() {
		t.Errorf("every finite chain is a lattice")
	}
	if p.Height() != 3 {
		t.Errorf("Height = %d, want 3", p.Height())
	}
	if !reflect.DeepEqual(p.MinimalElements(), []int{1}) {
		t.Errorf("minimal = %v", p.MinimalElements())
	}
	if !reflect.DeepEqual(p.MaximalElements(), []int{5}) {
		t.Errorf("maximal = %v", p.MaximalElements())
	}
}

func TestAntichainAndNonLattice(t *testing.T) {
	// Poset with two incomparable maximal elements (a "diamond" without a top):
	// 1 < 2, 1 < 3, so {2,3} is an antichain and has no join.
	set := NewIntSet(1, 2, 3)
	leq := RelationFromPairs([][2]int{{1, 1}, {2, 2}, {3, 3}, {1, 2}, {1, 3}})
	p, err := NewPoset(set, leq)
	if err != nil {
		t.Fatal(err)
	}
	if !p.IsAntichain([]int{2, 3}) {
		t.Errorf("{2,3} should be an antichain")
	}
	if _, ok := p.Join(2, 3); ok {
		t.Errorf("Join(2,3) should not exist")
	}
	if p.IsLattice() {
		t.Errorf("poset without a top element is not a lattice")
	}
	if !reflect.DeepEqual(p.MaximalElements(), []int{2, 3}) {
		t.Errorf("maximal = %v", p.MaximalElements())
	}
}

func TestNewPosetValidation(t *testing.T) {
	set := NewIntSet(1, 2)
	// Not transitive: 1<=2, 2<=3 missing 3 entirely and no reflexive.
	bad := RelationFromPairs([][2]int{{1, 2}})
	if _, err := NewPoset(set, bad); err == nil {
		t.Errorf("expected error for non-reflexive relation")
	}
	// References element outside the set.
	set2 := NewIntSet(1)
	outside := RelationFromPairs([][2]int{{1, 1}, {1, 9}})
	if _, err := NewPoset(set2, outside); err == nil {
		t.Errorf("expected error for out-of-set reference")
	}
	if _, err := NewDivisibilityPoset(0); err == nil {
		t.Errorf("expected error for n < 1")
	}
}

func BenchmarkTransitiveClosure(b *testing.B) {
	// Build a dense-ish random-looking but deterministic relation on 60 nodes:
	// a chain plus extra forward edges, the heaviest routine (Warshall O(n^3)).
	const n = 60
	r := make(Relation)
	for i := 0; i < n-1; i++ {
		r.Add(i, i+1)
	}
	for i := 0; i < n; i++ {
		r.Add(i, (i*7+3)%n)
	}
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = r.TransitiveClosure()
	}
}
