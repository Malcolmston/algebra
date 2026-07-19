package matroids

import (
	"fmt"
	"reflect"
	"testing"
)

// ---------------------------------------------------------------------------
// Basic rank / independence properties across representations.
// ---------------------------------------------------------------------------

func TestUniformMatroidBasics(t *testing.T) {
	m := NewUniformMatroid(2, 4)
	if got := m.Size(); got != 4 {
		t.Fatalf("Size = %d, want 4", got)
	}
	if got := FullRank(m); got != 2 {
		t.Fatalf("FullRank = %d, want 2", got)
	}
	tests := []struct {
		set  []int
		rank int
	}{
		{nil, 0},
		{[]int{0}, 1},
		{[]int{0, 1}, 2},
		{[]int{0, 1, 2}, 2},
		{[]int{0, 1, 2, 3}, 2},
		{[]int{2, 2, 2}, 1},
	}
	for _, tc := range tests {
		if got := m.Rank(tc.set); got != tc.rank {
			t.Errorf("Rank(%v) = %d, want %d", tc.set, got, tc.rank)
		}
	}
	if !Independent(m, []int{0, 1}) {
		t.Error("{0,1} should be independent")
	}
	if Independent(m, []int{0, 1, 2}) {
		t.Error("{0,1,2} should be dependent")
	}
	if n := BasisCount(m); n != 6 {
		t.Errorf("BasisCount = %d, want 6", n)
	}
	if cs := Circuits(m); len(cs) != 4 {
		t.Errorf("number of circuits = %d, want 4", len(cs))
	}
	if g := Girth(m); g != 3 {
		t.Errorf("Girth = %d, want 3", g)
	}
	if !IsUniform(m) {
		t.Error("U(2,4) should be recognised as uniform")
	}
}

func TestFreeAndTrivial(t *testing.T) {
	free := NewFreeMatroid(3)
	if FullRank(free) != 3 {
		t.Errorf("free FullRank = %d, want 3", FullRank(free))
	}
	if got := Coloops(free); !reflect.DeepEqual(got, []int{0, 1, 2}) {
		t.Errorf("free Coloops = %v, want [0 1 2]", got)
	}
	if len(Circuits(free)) != 0 {
		t.Error("free matroid has no circuits")
	}
	if Girth(free) != -1 {
		t.Errorf("free Girth = %d, want -1 (infinite)", Girth(free))
	}
	triv := NewTrivialMatroid(3)
	if FullRank(triv) != 0 {
		t.Errorf("trivial FullRank = %d, want 0", FullRank(triv))
	}
	if got := Loops(triv); !reflect.DeepEqual(got, []int{0, 1, 2}) {
		t.Errorf("trivial Loops = %v, want [0 1 2]", got)
	}
}

func TestPartitionMatroid(t *testing.T) {
	// two blocks: {0,1,2} cap 1, {3,4} cap 2
	m := NewPartitionMatroidFromBlocks(5, [][]int{{0, 1, 2}, {3, 4}}, []int{1, 2})
	tests := []struct {
		set  []int
		rank int
	}{
		{[]int{0, 1, 2}, 1},
		{[]int{0, 3, 4}, 3},
		{[]int{0, 1, 3, 4}, 3},
		{[]int{3, 4}, 2},
	}
	for _, tc := range tests {
		if got := m.Rank(tc.set); got != tc.rank {
			t.Errorf("Rank(%v) = %d, want %d", tc.set, got, tc.rank)
		}
	}
	if FullRank(m) != 3 {
		t.Errorf("FullRank = %d, want 3", FullRank(m))
	}
}

func TestGraphicMatroidK4(t *testing.T) {
	m := CycleMatroidOfCompleteGraph(4)
	if m.Size() != 6 {
		t.Fatalf("K4 edges = %d, want 6", m.Size())
	}
	if FullRank(m) != 3 {
		t.Errorf("FullRank = %d, want 3", FullRank(m))
	}
	// Number of spanning trees of K4 is 4^2 = 16 (Cayley).
	if bc := BasisCount(m); bc != 16 {
		t.Errorf("spanning trees = %d, want 16", bc)
	}
	if !IsConnected(m) {
		t.Error("K4 cycle matroid should be connected")
	}
}

func TestGraphicMatroidTriangle(t *testing.T) {
	m := CycleMatroidOfCycleGraph(3)
	if FullRank(m) != 2 {
		t.Errorf("FullRank = %d, want 2", FullRank(m))
	}
	cs := Circuits(m)
	if len(cs) != 1 || !reflect.DeepEqual(cs[0], []int{0, 1, 2}) {
		t.Errorf("circuits = %v, want [[0 1 2]]", cs)
	}
	// self loop is a loop element
	g := NewGraphicMatroid(2)
	g.AddEdge(0, 0)
	g.AddEdge(0, 1)
	if !IsLoop(g, 0) {
		t.Error("self-loop edge should be a matroid loop")
	}
	if !IsColoop(g, 1) {
		t.Error("bridge edge should be a coloop")
	}
}

func TestLinearMatroid(t *testing.T) {
	// Columns represent U(2,3): any two columns independent, all three dependent.
	m := NewLinearMatroid([][]int64{
		{1, 0, 1},
		{0, 1, 1},
	})
	if FullRank(m) != 2 {
		t.Errorf("FullRank = %d, want 2", FullRank(m))
	}
	if m.Rank([]int{0, 1, 2}) != 2 {
		t.Errorf("Rank(all) = %d, want 2", m.Rank([]int{0, 1, 2}))
	}
	if !Independent(m, []int{0, 2}) {
		t.Error("{0,2} should be independent")
	}
	if !IsCircuit(m, []int{0, 1, 2}) {
		t.Error("{0,1,2} should be a circuit")
	}
	// a zero column is a loop
	m2 := NewLinearMatroid([][]int64{
		{1, 0, 0},
		{0, 0, 1},
	})
	if !IsLoop(m2, 1) {
		t.Error("zero column should be a loop")
	}
}

func TestBinaryMatroidFano(t *testing.T) {
	// Fano plane F7: 7 columns, rank 3 over GF(2), every 2 columns independent,
	// exactly the 7 lines are the 3-element circuits.
	cols := [][]int{
		{1, 0, 0, 1, 1, 0, 1},
		{0, 1, 0, 1, 0, 1, 1},
		{0, 0, 1, 0, 1, 1, 1},
	}
	m := NewBinaryMatroid(cols)
	if m.Size() != 7 {
		t.Fatalf("Size = %d, want 7", m.Size())
	}
	if FullRank(m) != 3 {
		t.Errorf("FullRank = %d, want 3", FullRank(m))
	}
	if !IsSimple(m) {
		t.Error("Fano matroid is simple")
	}
	// The Fano matroid has exactly 7 three-element circuits (its lines) plus
	// larger circuits; check girth is 3.
	if Girth(m) != 3 {
		t.Errorf("Girth = %d, want 3", Girth(m))
	}
	if err := CheckRankAxioms(m); err != nil {
		t.Errorf("Fano rank axioms: %v", err)
	}
}

func TestTransversalMatroid(t *testing.T) {
	// element 0 -> {a,b}; 1 -> {a}; 2 -> {a}
	const a, b = 10, 11
	m := NewTransversalMatroid(3, [][]int{{a, b}, {a}, {a}})
	if FullRank(m) != 2 {
		t.Errorf("FullRank = %d, want 2", FullRank(m))
	}
	if Independent(m, []int{1, 2}) {
		t.Error("{1,2} share the only target, should be dependent")
	}
	if !Independent(m, []int{0, 1}) {
		t.Error("{0,1} should be independent")
	}
	if _, ok := m.SystemOfDistinctRepresentatives([]int{0, 1}); !ok {
		t.Error("{0,1} should have an SDR")
	}
	if _, ok := m.SystemOfDistinctRepresentatives([]int{1, 2}); ok {
		t.Error("{1,2} should not have an SDR")
	}
}

func TestExplicitMatroid(t *testing.T) {
	m := NewExplicitMatroidFromBases(3, [][]int{{0, 1}, {0, 2}, {1, 2}})
	if FullRank(m) != 2 {
		t.Errorf("FullRank = %d, want 2", FullRank(m))
	}
	if m.Rank([]int{0, 1, 2}) != 2 {
		t.Errorf("Rank(all) = %d, want 2", m.Rank([]int{0, 1, 2}))
	}
	if err := m.Validate(); err != nil {
		t.Errorf("Validate: %v", err)
	}
	// From circuits: single circuit {0,1,2} on 3 elements gives U(2,3).
	mc := NewExplicitMatroidFromCircuits(3, [][]int{{0, 1, 2}})
	if FullRank(mc) != 2 {
		t.Errorf("circuit-built FullRank = %d, want 2", FullRank(mc))
	}
	if !IsUniform(mc) {
		t.Error("should be uniform U(2,3)")
	}
	// A non-matroid family should fail the axiom check.
	if IsMatroidIndependenceSystem(2, [][]int{{}, {0}, {0, 1}}) {
		// {0,1} independent but {1} not listed -> not downward closed
		t.Error("non-downward-closed family should not be a matroid")
	}
}

// ---------------------------------------------------------------------------
// Duality, minors, direct sums.
// ---------------------------------------------------------------------------

func TestDualityUniformSelfDual(t *testing.T) {
	m := NewUniformMatroid(2, 4)
	d := Dual(m)
	if FullRank(d) != 2 {
		t.Errorf("dual FullRank = %d, want 2 (self-dual)", FullRank(d))
	}
	if !IsUniform(d) {
		t.Error("dual of U(2,4) should be uniform U(2,4)")
	}
	// bases of dual are complements of bases
	if len(Cobases(m)) != BasisCount(m) {
		t.Error("cobasis count mismatch")
	}
	if err := CheckRankAxioms(d); err != nil {
		t.Errorf("dual rank axioms: %v", err)
	}
}

func TestMinors(t *testing.T) {
	m := NewUniformMatroid(2, 4)
	// deletion of one element from U(2,4) gives U(2,3)
	del := Deletion(m, []int{3})
	if del.Size() != 3 || FullRank(del) != 2 {
		t.Errorf("deletion: size=%d rank=%d, want 3,2", del.Size(), FullRank(del))
	}
	// contraction of one non-loop element from U(2,4) gives U(1,3)
	con := Contraction(m, []int{3})
	if con.Size() != 3 || FullRank(con) != 1 {
		t.Errorf("contraction: size=%d rank=%d, want 3,1", con.Size(), FullRank(con))
	}
	if got := con.Elements(); !reflect.DeepEqual(got, []int{0, 1, 2}) {
		t.Errorf("minor Elements = %v, want [0 1 2]", got)
	}
	if err := CheckRankAxioms(con); err != nil {
		t.Errorf("contraction rank axioms: %v", err)
	}
}

func TestDirectSum(t *testing.T) {
	a := NewUniformMatroid(1, 2)
	b := NewUniformMatroid(1, 2)
	ds := DirectSum(a, b)
	if ds.Size() != 4 {
		t.Fatalf("Size = %d, want 4", ds.Size())
	}
	if FullRank(ds) != 2 {
		t.Errorf("FullRank = %d, want 2", FullRank(ds))
	}
	// {0,1} is a circuit (first summand), {0,2} independent (different parts)
	if !IsCircuit(ds, []int{0, 1}) {
		t.Error("{0,1} should be a circuit in the first summand")
	}
	if !Independent(ds, []int{0, 2}) {
		t.Error("{0,2} in different summands should be independent")
	}
	if NumComponents(ds) != 2 {
		t.Errorf("components = %d, want 2", NumComponents(ds))
	}
	if IsConnected(ds) {
		t.Error("a direct sum of two nontrivial matroids is disconnected")
	}
}

// ---------------------------------------------------------------------------
// Greedy, intersection, union.
// ---------------------------------------------------------------------------

func TestGreedy(t *testing.T) {
	m := NewUniformMatroid(2, 4)
	w := []float64{1, 4, 2, 3}
	res := Greedy(m, w)
	if !reflect.DeepEqual(res.Set, []int{1, 3}) {
		t.Errorf("Greedy set = %v, want [1 3]", res.Set)
	}
	if res.Weight != 7 {
		t.Errorf("Greedy weight = %v, want 7", res.Weight)
	}
}

func TestGreedyMSTIsMinWeightBasis(t *testing.T) {
	// triangle with edge weights 1,2,3; MST picks the two lightest edges.
	g := NewGraphicMatroidFromEdges(3, [][2]int{{0, 1}, {1, 2}, {0, 2}})
	w := []float64{1, 2, 3}
	res := GreedyMinWeightBasis(g, w)
	if res.Weight != 3 { // edges 0 and 1
		t.Errorf("MST weight = %v, want 3", res.Weight)
	}
	if !reflect.DeepEqual(res.Set, []int{0, 1}) {
		t.Errorf("MST edges = %v, want [0 1]", res.Set)
	}
}

// bipartiteEdgeMatroids builds the two partition matroids whose common
// independent sets are the matchings of a bipartite graph on the given edges.
func bipartiteEdgeMatroids(nl, nr int, edges [][2]int) (*PartitionMatroid, *PartitionMatroid) {
	n := len(edges)
	leftBlock := make([]int, n)
	rightBlock := make([]int, n)
	for i, e := range edges {
		leftBlock[i] = e[0]
		rightBlock[i] = e[1]
	}
	lc := make([]int, nl)
	for i := range lc {
		lc[i] = 1
	}
	rc := make([]int, nr)
	for i := range rc {
		rc[i] = 1
	}
	return NewPartitionMatroid(n, leftBlock, lc), NewPartitionMatroid(n, rightBlock, rc)
}

func TestMatroidIntersectionMatching(t *testing.T) {
	// K_{2,2}: edges e0=(0,0) e1=(0,1) e2=(1,0) e3=(1,1)
	edges := [][2]int{{0, 0}, {0, 1}, {1, 0}, {1, 1}}
	m1, m2 := bipartiteEdgeMatroids(2, 2, edges)
	if got := IntersectionSize(m1, m2); got != 2 {
		t.Errorf("max matching = %d, want 2", got)
	}
	common := Intersection(m1, m2)
	if len(common) != 2 {
		t.Errorf("common independent set size = %d, want 2", len(common))
	}
	// it must be independent in both
	if !Independent(m1, common) || !Independent(m2, common) {
		t.Errorf("returned set %v not common independent", common)
	}
}

func TestWeightedMatroidIntersection(t *testing.T) {
	edges := [][2]int{{0, 0}, {0, 1}, {1, 0}, {1, 1}}
	m1, m2 := bipartiteEdgeMatroids(2, 2, edges)
	// weights favour the diagonal e1,e2
	w := []float64{1, 5, 5, 1}
	if got := WeightedIntersectionWeight(m1, m2, w); got != 10 {
		t.Errorf("max weight matching = %v, want 10", got)
	}
	set := WeightedIntersection(m1, m2, w)
	if !reflect.DeepEqual(set, []int{1, 2}) {
		t.Errorf("weighted intersection set = %v, want [1 2]", set)
	}
}

func TestMatroidUnion(t *testing.T) {
	// union of two U(1,3) is U(2,3)
	a := NewUniformMatroid(1, 3)
	b := NewUniformMatroid(1, 3)
	u := Union(a, b)
	if u.Size() != 3 {
		t.Fatalf("Size = %d, want 3", u.Size())
	}
	if got := u.Rank([]int{0, 1, 2}); got != 2 {
		t.Errorf("union rank = %d, want 2", got)
	}
	if mis := u.MaxIndependentSet(); len(mis) != 2 {
		t.Errorf("union max independent size = %d, want 2", len(mis))
	}
	if !IsUniform(u) {
		t.Error("union of two U(1,3) should be U(2,3)")
	}
	// partition assigns disjoint singletons to each part
	part := u.Partition()
	if len(part) != 2 {
		t.Fatalf("partition parts = %d, want 2", len(part))
	}
	total := len(part[0]) + len(part[1])
	if total != 2 {
		t.Errorf("partition total = %d, want 2", total)
	}
}

func TestUnionPartitionable(t *testing.T) {
	// two spanning-tree matroids of a triangle: 3 edges partition into 2
	// forests? A triangle has 3 edges, arboricity 2, so it is NOT partitionable
	// into 1 forest but IS the union of two forests covering all edges only if
	// union rank == 3. Union rank of two graphic matroids of C3 = 3.
	g1 := NewGraphicMatroidFromEdges(3, [][2]int{{0, 1}, {1, 2}, {0, 2}})
	g2 := NewGraphicMatroidFromEdges(3, [][2]int{{0, 1}, {1, 2}, {0, 2}})
	if !IsPartitionable([]Matroid{g1, g2}) {
		t.Error("triangle should be partitionable into two forests")
	}
	if got := UnionRank([]Matroid{g1, g2}); got != 3 {
		t.Errorf("union rank = %d, want 3", got)
	}
}

// ---------------------------------------------------------------------------
// Closure, flats, connectivity.
// ---------------------------------------------------------------------------

func TestClosureAndFlats(t *testing.T) {
	m := NewUniformMatroid(2, 4)
	// closure of a rank-2 independent set is the whole ground set
	if got := Closure(m, []int{0, 1}); !reflect.DeepEqual(got, []int{0, 1, 2, 3}) {
		t.Errorf("Closure({0,1}) = %v, want [0 1 2 3]", got)
	}
	// closure of a single element is itself (rank 1, no parallels)
	if got := Closure(m, []int{0}); !reflect.DeepEqual(got, []int{0}) {
		t.Errorf("Closure({0}) = %v, want [0]", got)
	}
	// flats of rank 1 are the four singletons; rank-0 flat is empty set
	if got := FlatsOfRank(m, 1); len(got) != 4 {
		t.Errorf("rank-1 flats = %d, want 4", len(got))
	}
	if got := Hyperplanes(m); len(got) != 4 {
		t.Errorf("hyperplanes = %d, want 4", len(got))
	}
}

func TestConnectivityFunction(t *testing.T) {
	m := NewUniformMatroid(2, 4)
	// lambda of {0,1} in U(2,4): r({0,1}) + r({2,3}) - r(E) = 2+2-2 = 2
	if got := ConnectivityFunction(m, []int{0, 1}); got != 2 {
		t.Errorf("lambda({0,1}) = %d, want 2", got)
	}
	if IsSeparator(m, []int{0, 1}) {
		t.Error("{0,1} is not a separator of connected U(2,4)")
	}
}

// ---------------------------------------------------------------------------
// Axiom checks across all representations.
// ---------------------------------------------------------------------------

func TestRankAxiomsAcrossRepresentations(t *testing.T) {
	ms := map[string]Matroid{
		"uniform":   NewUniformMatroid(2, 5),
		"free":      NewFreeMatroid(4),
		"trivial":   NewTrivialMatroid(4),
		"partition": NewPartitionMatroidFromBlocks(5, [][]int{{0, 1, 2}, {3, 4}}, []int{1, 2}),
		"graphic":   CycleMatroidOfCompleteGraph(4),
		"linear":    NewLinearMatroid([][]int64{{1, 0, 1, 2}, {0, 1, 1, 3}}),
		"binary":    NewBinaryMatroid([][]int{{1, 0, 1}, {0, 1, 1}}),
		"transvers": NewTransversalMatroid(4, [][]int{{0, 1}, {0}, {1, 2}, {2}}),
		"explicit":  NewExplicitMatroidFromBases(4, [][]int{{0, 1}, {0, 2}, {1, 2}, {0, 3}, {1, 3}, {2, 3}}),
	}
	for name, m := range ms {
		if err := CheckRankAxioms(m); err != nil {
			t.Errorf("%s: rank axioms fail: %v", name, err)
		}
		if err := CheckRankAxioms(Dual(m)); err != nil {
			t.Errorf("%s dual: rank axioms fail: %v", name, err)
		}
	}
}

func TestSetHelpers(t *testing.T) {
	if !SetEqual(SetUnion([]int{1, 2}, []int{2, 3}), []int{1, 2, 3}) {
		t.Error("SetUnion wrong")
	}
	if !SetEqual(SetIntersection([]int{1, 2, 3}, []int{2, 3, 4}), []int{2, 3}) {
		t.Error("SetIntersection wrong")
	}
	if !SetEqual(SetDifference([]int{1, 2, 3}, []int{2}), []int{1, 3}) {
		t.Error("SetDifference wrong")
	}
	if !SetEqual(Complement(4, []int{1, 3}), []int{0, 2}) {
		t.Error("Complement wrong")
	}
	if len(Subsets([]int{0, 1, 2})) != 8 {
		t.Error("Subsets count wrong")
	}
	if len(SubsetsOfSize([]int{0, 1, 2, 3}, 2)) != 6 {
		t.Error("SubsetsOfSize count wrong")
	}
}

// ---------------------------------------------------------------------------
// Runnable examples.
// ---------------------------------------------------------------------------

func ExampleGreedy() {
	m := NewUniformMatroid(2, 4)
	weights := []float64{1, 4, 2, 3}
	res := Greedy(m, weights)
	fmt.Println(res.Set, res.Weight)
	// Output: [1 3] 7
}

func ExampleIntersection() {
	// Maximum bipartite matching in K_{2,2} via matroid intersection.
	edges := [][2]int{{0, 0}, {0, 1}, {1, 0}, {1, 1}}
	n := len(edges)
	left := make([]int, n)
	right := make([]int, n)
	for i, e := range edges {
		left[i], right[i] = e[0], e[1]
	}
	m1 := NewPartitionMatroid(n, left, []int{1, 1})
	m2 := NewPartitionMatroid(n, right, []int{1, 1})
	fmt.Println(IntersectionSize(m1, m2))
	// Output: 2
}

func ExampleDual() {
	m := NewUniformMatroid(2, 4)
	d := Dual(m)
	fmt.Println(FullRank(d))
	// Output: 2
}
