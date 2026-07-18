package settheory

import "sort"

// Poset is a finite partially ordered set: a set of integer elements together
// with a reflexive, antisymmetric and transitive "less than or equal" relation.
// The zero value is not usable; build posets with NewPoset, NewChainPoset or
// NewDivisibilityPoset.
type Poset struct {
	elements IntSet
	leq      Relation
}

// NewPoset returns a Poset over set whose order is leq. It returns an error if
// leq relates any value outside set, or if leq is not a partial order on set
// (reflexive on set, antisymmetric and transitive).
func NewPoset(set IntSet, leq Relation) (Poset, error) {
	for p := range leq {
		if _, ok := set[p.From]; !ok {
			return Poset{}, settheoryErrorf("NewPoset: relation references %d which is not in the element set", p.From)
		}
		if _, ok := set[p.To]; !ok {
			return Poset{}, settheoryErrorf("NewPoset: relation references %d which is not in the element set", p.To)
		}
	}
	if !leq.IsReflexiveOn(set) {
		return Poset{}, settheoryErrorf("NewPoset: relation is not reflexive on the element set")
	}
	if !leq.IsAntisymmetric() {
		return Poset{}, settheoryErrorf("NewPoset: relation is not antisymmetric")
	}
	if !leq.IsTransitive() {
		return Poset{}, settheoryErrorf("NewPoset: relation is not transitive")
	}
	return Poset{elements: set.Clone(), leq: leq.Clone()}, nil
}

// NewChainPoset returns the totally ordered Poset on the given elements ordered
// by ascending integer value. Duplicate elements are collapsed.
func NewChainPoset(elems ...int) Poset {
	set := NewIntSet(elems...)
	leq := make(Relation)
	sorted := set.Elements()
	for i := range sorted {
		for j := i; j < len(sorted); j++ {
			leq[Pair{From: sorted[i], To: sorted[j]}] = struct{}{}
		}
	}
	return Poset{elements: set, leq: leq}
}

// NewDivisibilityPoset returns the divisibility Poset on the positive divisors
// of n, where a ≤ b means a divides b. It returns an error for n < 1. This poset
// is always a lattice, with meet given by greatest common divisor and join by
// least common multiple.
func NewDivisibilityPoset(n int) (Poset, error) {
	if n < 1 {
		return Poset{}, settheoryErrorf("NewDivisibilityPoset: n must be positive, got %d", n)
	}
	set := make(IntSet)
	for d := 1; d <= n; d++ {
		if n%d == 0 {
			set[d] = struct{}{}
		}
	}
	leq := make(Relation)
	divs := set.Elements()
	for _, a := range divs {
		for _, b := range divs {
			if b%a == 0 {
				leq[Pair{From: a, To: b}] = struct{}{}
			}
		}
	}
	return Poset{elements: set, leq: leq}, nil
}

// Elements returns the elements of the poset in ascending order.
func (p Poset) Elements() []int {
	return p.elements.Elements()
}

// Size returns the number of elements in the poset.
func (p Poset) Size() int {
	return len(p.elements)
}

// Leq reports whether a ≤ b in the poset order.
func (p Poset) Leq(a, b int) bool {
	return p.leq.Related(a, b)
}

// Less reports whether a < b, i.e. a ≤ b and a ≠ b.
func (p Poset) Less(a, b int) bool {
	return a != b && p.leq.Related(a, b)
}

// AreComparable reports whether a ≤ b or b ≤ a.
func (p Poset) AreComparable(a, b int) bool {
	return p.leq.Related(a, b) || p.leq.Related(b, a)
}

// Covers reports whether upper covers lower: lower < upper and there is no
// element c with lower < c < upper.
func (p Poset) Covers(lower, upper int) bool {
	if !p.Less(lower, upper) {
		return false
	}
	for c := range p.elements {
		if c == lower || c == upper {
			continue
		}
		if p.Less(lower, c) && p.Less(c, upper) {
			return false
		}
	}
	return true
}

// HasseEdges returns the covering pairs of the poset as Pairs (From = lower,
// To = upper) with To covering From. Edges are sorted lexicographically. These
// are exactly the edges of the poset's Hasse diagram.
func (p Poset) HasseEdges() []Pair {
	var edges []Pair
	elems := p.elements.Elements()
	for _, lo := range elems {
		for _, hi := range elems {
			if p.Covers(lo, hi) {
				edges = append(edges, Pair{From: lo, To: hi})
			}
		}
	}
	sort.Slice(edges, func(i, j int) bool {
		if edges[i].From != edges[j].From {
			return edges[i].From < edges[j].From
		}
		return edges[i].To < edges[j].To
	})
	return edges
}

// TopologicalOrder returns a linear extension of the poset: a permutation of the
// elements in which a appears before b whenever a < b. The extension is made
// deterministic by always selecting the smallest currently-minimal element.
func (p Poset) TopologicalOrder() []int {
	remaining := p.elements.Clone()
	var order []int
	for len(remaining) > 0 {
		// Find the smallest element that is minimal among those remaining.
		cands := remaining.Elements() // ascending
		pick := cands[0]
		found := false
		for _, x := range cands {
			isMinimal := true
			for y := range remaining {
				if y != x && p.Less(y, x) {
					isMinimal = false
					break
				}
			}
			if isMinimal {
				pick = x
				found = true
				break
			}
		}
		if !found {
			// Should be impossible for a valid partial order, but guard anyway.
			pick = cands[0]
		}
		order = append(order, pick)
		delete(remaining, pick)
	}
	return order
}

// MinimalElements returns the elements that have no strictly smaller element,
// in ascending order.
func (p Poset) MinimalElements() []int {
	var out []int
	for _, x := range p.elements.Elements() {
		minimal := true
		for y := range p.elements {
			if y != x && p.Less(y, x) {
				minimal = false
				break
			}
		}
		if minimal {
			out = append(out, x)
		}
	}
	return out
}

// MaximalElements returns the elements that have no strictly larger element, in
// ascending order.
func (p Poset) MaximalElements() []int {
	var out []int
	for _, x := range p.elements.Elements() {
		maximal := true
		for y := range p.elements {
			if y != x && p.Less(x, y) {
				maximal = false
				break
			}
		}
		if maximal {
			out = append(out, x)
		}
	}
	return out
}

// LeastElement returns the unique element that is ≤ every element of the poset,
// if one exists. The boolean result reports existence.
func (p Poset) LeastElement() (int, bool) {
	for _, x := range p.elements.Elements() {
		isLeast := true
		for y := range p.elements {
			if !p.leq.Related(x, y) {
				isLeast = false
				break
			}
		}
		if isLeast {
			return x, true
		}
	}
	return 0, false
}

// GreatestElement returns the unique element that is ≥ every element of the
// poset, if one exists. The boolean result reports existence.
func (p Poset) GreatestElement() (int, bool) {
	for _, x := range p.elements.Elements() {
		isGreatest := true
		for y := range p.elements {
			if !p.leq.Related(y, x) {
				isGreatest = false
				break
			}
		}
		if isGreatest {
			return x, true
		}
	}
	return 0, false
}

// UpperBounds returns the elements that are ≥ every element of subset, in
// ascending order.
func (p Poset) UpperBounds(subset []int) []int {
	var out []int
	for _, x := range p.elements.Elements() {
		ub := true
		for _, s := range subset {
			if !p.leq.Related(s, x) {
				ub = false
				break
			}
		}
		if ub {
			out = append(out, x)
		}
	}
	return out
}

// LowerBounds returns the elements that are ≤ every element of subset, in
// ascending order.
func (p Poset) LowerBounds(subset []int) []int {
	var out []int
	for _, x := range p.elements.Elements() {
		lb := true
		for _, s := range subset {
			if !p.leq.Related(x, s) {
				lb = false
				break
			}
		}
		if lb {
			out = append(out, x)
		}
	}
	return out
}

// Join returns the least upper bound (supremum) of a and b, if it exists. The
// boolean result reports existence; the join exists when the set of common
// upper bounds has a unique least element.
func (p Poset) Join(a, b int) (int, bool) {
	ubs := p.UpperBounds([]int{a, b})
	return p.settheoryLeastOf(ubs)
}

// Meet returns the greatest lower bound (infimum) of a and b, if it exists. The
// boolean result reports existence; the meet exists when the set of common
// lower bounds has a unique greatest element.
func (p Poset) Meet(a, b int) (int, bool) {
	lbs := p.LowerBounds([]int{a, b})
	return p.settheoryGreatestOf(lbs)
}

// IsLattice reports whether every pair of elements has both a meet and a join,
// i.e. the poset is a lattice. The empty poset and any single-element poset are
// trivially lattices.
func (p Poset) IsLattice() bool {
	elems := p.elements.Elements()
	for i := range elems {
		for j := i; j < len(elems); j++ {
			if _, ok := p.Join(elems[i], elems[j]); !ok {
				return false
			}
			if _, ok := p.Meet(elems[i], elems[j]); !ok {
				return false
			}
		}
	}
	return true
}

// IsChain reports whether the poset is totally ordered: every two elements are
// comparable.
func (p Poset) IsChain() bool {
	elems := p.elements.Elements()
	for i := range elems {
		for j := i + 1; j < len(elems); j++ {
			if !p.AreComparable(elems[i], elems[j]) {
				return false
			}
		}
	}
	return true
}

// IsAntichain reports whether subset is an antichain: no two distinct elements
// of subset are comparable.
func (p Poset) IsAntichain(subset []int) bool {
	for i := range subset {
		for j := i + 1; j < len(subset); j++ {
			if p.AreComparable(subset[i], subset[j]) {
				return false
			}
		}
	}
	return true
}

// Height returns the number of elements in a longest chain of the poset. The
// empty poset has height zero. The height of a poset equals the length of its
// longest path in the covering (Hasse) digraph plus one.
func (p Poset) Height() int {
	order := p.TopologicalOrder()
	longest := make(map[int]int, len(order))
	best := 0
	for _, x := range order {
		l := 1
		for _, y := range order {
			if y == x {
				break
			}
			if p.Less(y, x) && longest[y]+1 > l {
				l = longest[y] + 1
			}
		}
		longest[x] = l
		if l > best {
			best = l
		}
	}
	return best
}

// settheoryLeastOf returns the unique least element of cands under the poset
// order, if one exists.
func (p Poset) settheoryLeastOf(cands []int) (int, bool) {
	for _, x := range cands {
		isLeast := true
		for _, y := range cands {
			if !p.leq.Related(x, y) {
				isLeast = false
				break
			}
		}
		if isLeast {
			return x, true
		}
	}
	return 0, false
}

// settheoryGreatestOf returns the unique greatest element of cands under the
// poset order, if one exists.
func (p Poset) settheoryGreatestOf(cands []int) (int, bool) {
	for _, x := range cands {
		isGreatest := true
		for _, y := range cands {
			if !p.leq.Related(y, x) {
				isGreatest = false
				break
			}
		}
		if isGreatest {
			return x, true
		}
	}
	return 0, false
}
