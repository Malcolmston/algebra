package settheory

import (
	"fmt"
	"sort"
)

// Pair is an ordered pair of integers (From, To). It is used both as an element
// of a binary Relation and as an edge in the Hasse diagram of a Poset.
type Pair struct {
	From int
	To   int
}

// String renders the pair as "(from, to)".
func (p Pair) String() string {
	return fmt.Sprintf("(%d, %d)", p.From, p.To)
}

// Relation is a binary relation on the integers, represented as a set of
// ordered pairs. The zero value is not usable; create relations with
// NewRelation.
type Relation map[Pair]struct{}

// NewRelation returns a new Relation containing the given pairs. Duplicate pairs
// are collapsed.
func NewRelation(pairs ...Pair) Relation {
	r := make(Relation, len(pairs))
	for _, p := range pairs {
		r[p] = struct{}{}
	}
	return r
}

// RelationFromPairs builds a Relation from ordered pairs supplied as two-element
// arrays [from, to]. It is a convenience for table-driven construction.
func RelationFromPairs(pairs [][2]int) Relation {
	r := make(Relation, len(pairs))
	for _, p := range pairs {
		r[Pair{From: p[0], To: p[1]}] = struct{}{}
	}
	return r
}

// Add inserts the ordered pair (a, b) into the relation.
func (r Relation) Add(a, b int) {
	r[Pair{From: a, To: b}] = struct{}{}
}

// Remove deletes the ordered pair (a, b) from the relation.
func (r Relation) Remove(a, b int) {
	delete(r, Pair{From: a, To: b})
}

// Related reports whether the ordered pair (a, b) is in the relation.
func (r Relation) Related(a, b int) bool {
	_, ok := r[Pair{From: a, To: b}]
	return ok
}

// Len returns the number of ordered pairs in the relation.
func (r Relation) Len() int {
	return len(r)
}

// Clone returns an independent copy of the relation.
func (r Relation) Clone() Relation {
	c := make(Relation, len(r))
	for p := range r {
		c[p] = struct{}{}
	}
	return c
}

// Pairs returns the ordered pairs of the relation sorted lexicographically by
// (From, To). The result is freshly allocated.
func (r Relation) Pairs() []Pair {
	out := make([]Pair, 0, len(r))
	for p := range r {
		out = append(out, p)
	}
	sort.Slice(out, func(i, j int) bool {
		if out[i].From != out[j].From {
			return out[i].From < out[j].From
		}
		return out[i].To < out[j].To
	})
	return out
}

// Equal reports whether r and other contain exactly the same pairs.
func (r Relation) Equal(other Relation) bool {
	if len(r) != len(other) {
		return false
	}
	for p := range r {
		if _, ok := other[p]; !ok {
			return false
		}
	}
	return true
}

// Domain returns the set of first coordinates appearing in the relation.
func (r Relation) Domain() IntSet {
	d := make(IntSet)
	for p := range r {
		d[p.From] = struct{}{}
	}
	return d
}

// Range returns the set of second coordinates appearing in the relation.
func (r Relation) Range() IntSet {
	rng := make(IntSet)
	for p := range r {
		rng[p.To] = struct{}{}
	}
	return rng
}

// Field returns the union of the domain and range, i.e. every value that
// participates in the relation.
func (r Relation) Field() IntSet {
	f := make(IntSet)
	for p := range r {
		f[p.From] = struct{}{}
		f[p.To] = struct{}{}
	}
	return f
}

// Inverse returns the converse relation, containing (b, a) for every pair
// (a, b) in r.
func (r Relation) Inverse() Relation {
	inv := make(Relation, len(r))
	for p := range r {
		inv[Pair{From: p.To, To: p.From}] = struct{}{}
	}
	return inv
}

// Compose returns the composition r∘other defined so that (a, c) is in the
// result whenever (a, b) is in other and (b, c) is in r for some b. This is the
// standard "apply other, then r" ordering.
func (r Relation) Compose(other Relation) Relation {
	// Index r by its From coordinate for efficient lookup.
	byFrom := make(map[int][]int)
	for p := range r {
		byFrom[p.From] = append(byFrom[p.From], p.To)
	}
	out := make(Relation)
	for p := range other {
		for _, c := range byFrom[p.To] {
			out[Pair{From: p.From, To: c}] = struct{}{}
		}
	}
	return out
}

// IsReflexiveOn reports whether (x, x) is in the relation for every x in set.
func (r Relation) IsReflexiveOn(set IntSet) bool {
	for x := range set {
		if _, ok := r[Pair{From: x, To: x}]; !ok {
			return false
		}
	}
	return true
}

// IsIrreflexiveOn reports whether (x, x) is absent for every x in set.
func (r Relation) IsIrreflexiveOn(set IntSet) bool {
	for x := range set {
		if _, ok := r[Pair{From: x, To: x}]; ok {
			return false
		}
	}
	return true
}

// IsSymmetric reports whether (b, a) is in the relation whenever (a, b) is.
func (r Relation) IsSymmetric() bool {
	for p := range r {
		if _, ok := r[Pair{From: p.To, To: p.From}]; !ok {
			return false
		}
	}
	return true
}

// IsAntisymmetric reports whether the presence of both (a, b) and (b, a)
// implies a equals b.
func (r Relation) IsAntisymmetric() bool {
	for p := range r {
		if p.From == p.To {
			continue
		}
		if _, ok := r[Pair{From: p.To, To: p.From}]; ok {
			return false
		}
	}
	return true
}

// IsTransitive reports whether (a, c) is in the relation whenever both (a, b)
// and (b, c) are.
func (r Relation) IsTransitive() bool {
	byFrom := make(map[int][]int)
	for p := range r {
		byFrom[p.From] = append(byFrom[p.From], p.To)
	}
	for p := range r {
		for _, c := range byFrom[p.To] {
			if _, ok := r[Pair{From: p.From, To: c}]; !ok {
				return false
			}
		}
	}
	return true
}

// IsEquivalenceOn reports whether the relation is reflexive on set, symmetric
// and transitive, i.e. an equivalence relation over set.
func (r Relation) IsEquivalenceOn(set IntSet) bool {
	return r.IsReflexiveOn(set) && r.IsSymmetric() && r.IsTransitive()
}

// IsPartialOrderOn reports whether the relation is reflexive on set,
// antisymmetric and transitive, i.e. a (non-strict) partial order over set.
func (r Relation) IsPartialOrderOn(set IntSet) bool {
	return r.IsReflexiveOn(set) && r.IsAntisymmetric() && r.IsTransitive()
}

// ReflexiveClosureOn returns the smallest relation containing r that is
// reflexive on set. It adds (x, x) for every x in set.
func (r Relation) ReflexiveClosureOn(set IntSet) Relation {
	out := r.Clone()
	for x := range set {
		out[Pair{From: x, To: x}] = struct{}{}
	}
	return out
}

// SymmetricClosure returns the smallest symmetric relation containing r. It adds
// (b, a) for every pair (a, b).
func (r Relation) SymmetricClosure() Relation {
	out := r.Clone()
	for p := range r {
		out[Pair{From: p.To, To: p.From}] = struct{}{}
	}
	return out
}

// TransitiveClosure returns the smallest transitive relation containing r,
// computed with Warshall's algorithm over the field of the relation. The result
// is the reachability relation of r viewed as a directed graph.
func (r Relation) TransitiveClosure() Relation {
	verts := r.Field().Elements()
	n := len(verts)
	index := make(map[int]int, n)
	for i, v := range verts {
		index[v] = i
	}
	reach := make([][]bool, n)
	for i := range reach {
		reach[i] = make([]bool, n)
	}
	for p := range r {
		reach[index[p.From]][index[p.To]] = true
	}
	for k := 0; k < n; k++ {
		for i := 0; i < n; i++ {
			if !reach[i][k] {
				continue
			}
			for j := 0; j < n; j++ {
				if reach[k][j] {
					reach[i][j] = true
				}
			}
		}
	}
	out := make(Relation)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if reach[i][j] {
				out[Pair{From: verts[i], To: verts[j]}] = struct{}{}
			}
		}
	}
	return out
}

// EquivalenceClosureOn returns the smallest equivalence relation over set that
// contains r. It is the reflexive (on set), symmetric and transitive closure.
func (r Relation) EquivalenceClosureOn(set IntSet) Relation {
	return r.SymmetricClosure().ReflexiveClosureOn(set).TransitiveClosure()
}
