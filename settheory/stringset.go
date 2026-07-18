package settheory

import (
	"sort"
	"strconv"
	"strings"
)

// StringSet is a finite set of string values backed by a hash map. The zero
// value is not usable; create sets with NewStringSet. A nil StringSet behaves
// as an empty set for read-only operations but must not be mutated.
type StringSet map[string]struct{}

// NewStringSet returns a new StringSet containing the given elements. Duplicate
// elements are collapsed.
func NewStringSet(elems ...string) StringSet {
	s := make(StringSet, len(elems))
	for _, e := range elems {
		s[e] = struct{}{}
	}
	return s
}

// Add inserts x into the set. Adding an element already present is a no-op.
func (s StringSet) Add(x string) {
	s[x] = struct{}{}
}

// Remove deletes x from the set. Removing an absent element is a no-op.
func (s StringSet) Remove(x string) {
	delete(s, x)
}

// Contains reports whether x is a member of the set.
func (s StringSet) Contains(x string) bool {
	_, ok := s[x]
	return ok
}

// Len returns the cardinality of the set.
func (s StringSet) Len() int {
	return len(s)
}

// IsEmpty reports whether the set has no elements.
func (s StringSet) IsEmpty() bool {
	return len(s) == 0
}

// Clone returns an independent copy of the set that shares no storage with the
// receiver.
func (s StringSet) Clone() StringSet {
	c := make(StringSet, len(s))
	for e := range s {
		c[e] = struct{}{}
	}
	return c
}

// Elements returns the members of the set in ascending lexicographic order. The
// result is a freshly allocated slice; mutating it does not affect the set.
func (s StringSet) Elements() []string {
	out := make([]string, 0, len(s))
	for e := range s {
		out = append(out, e)
	}
	sort.Strings(out)
	return out
}

// Equal reports whether s and other contain exactly the same elements.
func (s StringSet) Equal(other StringSet) bool {
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
func (s StringSet) Union(other StringSet) StringSet {
	u := make(StringSet, len(s)+len(other))
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
func (s StringSet) Intersection(other StringSet) StringSet {
	small, large := s, other
	if len(large) < len(small) {
		small, large = large, small
	}
	i := make(StringSet)
	for e := range small {
		if _, ok := large[e]; ok {
			i[e] = struct{}{}
		}
	}
	return i
}

// Difference returns a new set containing every element that is in s but not in
// other.
func (s StringSet) Difference(other StringSet) StringSet {
	d := make(StringSet)
	for e := range s {
		if _, ok := other[e]; !ok {
			d[e] = struct{}{}
		}
	}
	return d
}

// SymmetricDifference returns a new set containing every element that is in
// exactly one of s and other.
func (s StringSet) SymmetricDifference(other StringSet) StringSet {
	d := make(StringSet)
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
func (s StringSet) IsSubsetOf(other StringSet) bool {
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
func (s StringSet) IsSupersetOf(other StringSet) bool {
	return other.IsSubsetOf(s)
}

// IsProperSubsetOf reports whether s is a subset of other and the two sets are
// not equal.
func (s StringSet) IsProperSubsetOf(other StringSet) bool {
	return len(s) < len(other) && s.IsSubsetOf(other)
}

// IsDisjoint reports whether s and other share no elements.
func (s StringSet) IsDisjoint(other StringSet) bool {
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

// PowerSet returns all 2^n subsets of the set, where n is its cardinality. The
// subsets are ordered by the binary counter over the ascending element list.
// PowerSet panics if the set has more than 20 elements.
func (s StringSet) PowerSet() []StringSet {
	elems := s.Elements()
	n := len(elems)
	if n > 20 {
		panic("settheory: PowerSet called on a set with more than 20 elements")
	}
	total := 1 << uint(n)
	out := make([]StringSet, 0, total)
	for mask := 0; mask < total; mask++ {
		sub := make(StringSet)
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
// "{a, b, c}", with elements quoted and in ascending order. The empty set
// renders as "{}".
func (s StringSet) String() string {
	elems := s.Elements()
	parts := make([]string, len(elems))
	for i, e := range elems {
		parts[i] = strconv.Quote(e)
	}
	return "{" + strings.Join(parts, ", ") + "}"
}
