package automata

import (
	"sort"
	"strconv"
	"strings"
)

// StateSet is a set of automaton states represented as a map from state index
// to presence. It is used throughout the package for subset construction,
// epsilon-closure computation and nondeterministic simulation.
type StateSet map[int]bool

// NewStateSet returns a StateSet containing exactly the given states.
func NewStateSet(states ...int) StateSet {
	s := make(StateSet, len(states))
	for _, q := range states {
		s[q] = true
	}
	return s
}

// Add inserts state q into the set and returns the set for chaining.
func (s StateSet) Add(q int) StateSet {
	s[q] = true
	return s
}

// Contains reports whether state q is a member of the set.
func (s StateSet) Contains(q int) bool {
	return s[q]
}

// Len returns the number of states in the set.
func (s StateSet) Len() int {
	return len(s)
}

// Clone returns an independent copy of the set.
func (s StateSet) Clone() StateSet {
	c := make(StateSet, len(s))
	for q := range s {
		c[q] = true
	}
	return c
}

// Sorted returns the members of the set as an ascending slice of ints.
func (s StateSet) Sorted() []int {
	out := make([]int, 0, len(s))
	for q := range s {
		out = append(out, q)
	}
	sort.Ints(out)
	return out
}

// Key returns a canonical string key for the set, suitable for use as a map
// index. Two StateSets have the same Key if and only if they contain the same
// states.
func (s StateSet) Key() string {
	sorted := s.Sorted()
	parts := make([]string, len(sorted))
	for i, q := range sorted {
		parts[i] = strconv.Itoa(q)
	}
	return strings.Join(parts, ",")
}

// Equal reports whether two StateSets contain exactly the same states.
func (s StateSet) Equal(other StateSet) bool {
	if len(s) != len(other) {
		return false
	}
	for q := range s {
		if !other[q] {
			return false
		}
	}
	return true
}

// Union returns a new StateSet containing every state in either operand.
func (s StateSet) Union(other StateSet) StateSet {
	out := s.Clone()
	for q := range other {
		out[q] = true
	}
	return out
}

// Intersect returns a new StateSet containing states present in both operands.
func (s StateSet) Intersect(other StateSet) StateSet {
	out := make(StateSet)
	for q := range s {
		if other[q] {
			out[q] = true
		}
	}
	return out
}

// IsEmpty reports whether the set contains no states.
func (s StateSet) IsEmpty() bool {
	return len(s) == 0
}

// sortedRunes returns a sorted, de-duplicated copy of rs.
func sortedRunes(rs []rune) []rune {
	seen := make(map[rune]bool, len(rs))
	out := make([]rune, 0, len(rs))
	for _, r := range rs {
		if !seen[r] {
			seen[r] = true
			out = append(out, r)
		}
	}
	sort.Slice(out, func(i, j int) bool { return out[i] < out[j] })
	return out
}

// mergeRunes returns the sorted union of two rune slices, de-duplicated.
func mergeRunes(a, b []rune) []rune {
	return sortedRunes(append(append([]rune{}, a...), b...))
}
