package settheory

import (
	"sort"
	"strconv"
	"strings"
)

// IntSet is a finite set of int values backed by a hash map. The zero value is
// not usable; create sets with NewIntSet. A nil IntSet behaves as an empty set
// for read-only operations but must not be mutated.
type IntSet map[int]struct{}

// NewIntSet returns a new IntSet containing the given elements. Duplicate
// elements are collapsed, so NewIntSet(1, 1, 2) has cardinality two.
func NewIntSet(elems ...int) IntSet {
	s := make(IntSet, len(elems))
	for _, e := range elems {
		s[e] = struct{}{}
	}
	return s
}

// Add inserts x into the set. Adding an element already present is a no-op.
func (s IntSet) Add(x int) {
	s[x] = struct{}{}
}

// Remove deletes x from the set. Removing an absent element is a no-op.
func (s IntSet) Remove(x int) {
	delete(s, x)
}

// Contains reports whether x is a member of the set.
func (s IntSet) Contains(x int) bool {
	_, ok := s[x]
	return ok
}

// Len returns the cardinality of the set.
func (s IntSet) Len() int {
	return len(s)
}

// IsEmpty reports whether the set has no elements.
func (s IntSet) IsEmpty() bool {
	return len(s) == 0
}

// Clone returns an independent copy of the set that shares no storage with the
// receiver.
func (s IntSet) Clone() IntSet {
	c := make(IntSet, len(s))
	for e := range s {
		c[e] = struct{}{}
	}
	return c
}

// Elements returns the members of the set in ascending order. The result is a
// freshly allocated slice; mutating it does not affect the set.
func (s IntSet) Elements() []int {
	out := make([]int, 0, len(s))
	for e := range s {
		out = append(out, e)
	}
	sort.Ints(out)
	return out
}

// Equal reports whether s and other contain exactly the same elements.
func (s IntSet) Equal(other IntSet) bool {
	if len(s) != len(other) {
		return false
	}
	for e := range s {
		if _, ok := other[e]; !ok {
			return false
		}
	}
	return true
}

// Union returns a new set containing every element that is in s or in other.
func (s IntSet) Union(other IntSet) IntSet {
	u := make(IntSet, len(s)+len(other))
	for e := range s {
		u[e] = struct{}{}
	}
	for e := range other {
		u[e] = struct{}{}
	}
	return u
}

// Intersection returns a new set containing every element that is in both s and
// other.
func (s IntSet) Intersection(other IntSet) IntSet {
	small, large := s, other
	if len(large) < len(small) {
		small, large = large, small
	}
	i := make(IntSet)
	for e := range small {
		if _, ok := large[e]; ok {
			i[e] = struct{}{}
		}
	}
	return i
}

// Difference returns a new set containing every element that is in s but not in
// other.
func (s IntSet) Difference(other IntSet) IntSet {
	d := make(IntSet)
	for e := range s {
		if _, ok := other[e]; !ok {
			d[e] = struct{}{}
		}
	}
	return d
}

// SymmetricDifference returns a new set containing every element that is in
// exactly one of s and other.
func (s IntSet) SymmetricDifference(other IntSet) IntSet {
	d := make(IntSet)
	for e := range s {
		if _, ok := other[e]; !ok {
			d[e] = struct{}{}
		}
	}
	for e := range other {
		if _, ok := s[e]; !ok {
			d[e] = struct{}{}
		}
	}
	return d
}

// IsSubsetOf reports whether every element of s is also an element of other.
// The empty set is a subset of every set.
func (s IntSet) IsSubsetOf(other IntSet) bool {
	if len(s) > len(other) {
		return false
	}
	for e := range s {
		if _, ok := other[e]; !ok {
			return false
		}
	}
	return true
}

// IsSupersetOf reports whether every element of other is also an element of s.
func (s IntSet) IsSupersetOf(other IntSet) bool {
	return other.IsSubsetOf(s)
}

// IsProperSubsetOf reports whether s is a subset of other and the two sets are
// not equal.
func (s IntSet) IsProperSubsetOf(other IntSet) bool {
	return len(s) < len(other) && s.IsSubsetOf(other)
}

// IsDisjoint reports whether s and other share no elements.
func (s IntSet) IsDisjoint(other IntSet) bool {
	small, large := s, other
	if len(large) < len(small) {
		small, large = large, small
	}
	for e := range small {
		if _, ok := large[e]; ok {
			return false
		}
	}
	return true
}

// Min returns the smallest element of the set. The boolean result is false when
// the set is empty.
func (s IntSet) Min() (int, bool) {
	first := true
	var m int
	for e := range s {
		if first || e < m {
			m, first = e, false
		}
	}
	return m, !first
}

// Max returns the largest element of the set. The boolean result is false when
// the set is empty.
func (s IntSet) Max() (int, bool) {
	first := true
	var m int
	for e := range s {
		if first || e > m {
			m, first = e, false
		}
	}
	return m, !first
}

// Sum returns the sum of all elements of the set; the empty set sums to zero.
func (s IntSet) Sum() int {
	total := 0
	for e := range s {
		total += e
	}
	return total
}

// CartesianProduct returns every ordered pair (a, b) with a in s and b in
// other, sorted lexicographically. The result has len(s)*len(other) entries.
func (s IntSet) CartesianProduct(other IntSet) []Pair {
	as := s.Elements()
	bs := other.Elements()
	out := make([]Pair, 0, len(as)*len(bs))
	for _, a := range as {
		for _, b := range bs {
			out = append(out, Pair{From: a, To: b})
		}
	}
	return out
}

// PowerSet returns all 2^n subsets of the set, where n is its cardinality. The
// subsets are ordered by the binary counter over the ascending element list, so
// the empty set is first and the full set is last. PowerSet panics if the set
// has more than 20 elements, since the result would exceed one million subsets.
func (s IntSet) PowerSet() []IntSet {
	elems := s.Elements()
	n := len(elems)
	if n > 20 {
		panic("settheory: PowerSet called on a set with more than 20 elements")
	}
	total := 1 << uint(n)
	out := make([]IntSet, 0, total)
	for mask := 0; mask < total; mask++ {
		sub := make(IntSet)
		for i := 0; i < n; i++ {
			if mask&(1<<uint(i)) != 0 {
				sub[elems[i]] = struct{}{}
			}
		}
		out = append(out, sub)
	}
	return out
}

// String returns a canonical textual representation of the set such as
// "{1, 2, 3}", with elements in ascending order. The empty set renders as
// "{}".
func (s IntSet) String() string {
	elems := s.Elements()
	parts := make([]string, len(elems))
	for i, e := range elems {
		parts[i] = strconv.Itoa(e)
	}
	return "{" + strings.Join(parts, ", ") + "}"
}
