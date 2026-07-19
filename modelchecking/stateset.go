package modelchecking

import (
	"fmt"
	"sort"
	"strings"
)

// StateSet is a fixed-universe set of states drawn from {0, 1, ..., n-1}.
// It is represented as a packed bitset so that the Boolean operations used by
// the CTL fixpoint algorithms (union, intersection, complement) run in time
// proportional to n/64. The universe size n is fixed at construction and every
// binary operation requires both operands to share the same universe.
type StateSet struct {
	n     int
	words []uint64
}

// NewStateSet returns an empty set over the universe {0, ..., n-1}. It panics
// only if n is negative, which cannot arise from a valid Kripke structure.
func NewStateSet(n int) StateSet {
	if n < 0 {
		n = 0
	}
	return StateSet{n: n, words: make([]uint64, (n+63)/64)}
}

// FullStateSet returns the set containing every state of the universe
// {0, ..., n-1}.
func FullStateSet(n int) StateSet {
	s := NewStateSet(n)
	for i := 0; i < n; i++ {
		s.Add(i)
	}
	return s
}

// StateSetFromSlice builds a set over the universe {0, ..., n-1} containing
// exactly the listed elements. Out-of-range elements are ignored.
func StateSetFromSlice(n int, elems []int) StateSet {
	s := NewStateSet(n)
	for _, e := range elems {
		s.Add(e)
	}
	return s
}

// Universe returns the size n of the universe over which the set is defined.
func (s StateSet) Universe() int { return s.n }

// Add inserts state i into the set. Out-of-range indices are ignored so that
// callers may add speculatively without bounds checking.
func (s StateSet) Add(i int) {
	if i < 0 || i >= s.n {
		return
	}
	s.words[i>>6] |= 1 << uint(i&63)
}

// Remove deletes state i from the set. Out-of-range indices are ignored.
func (s StateSet) Remove(i int) {
	if i < 0 || i >= s.n {
		return
	}
	s.words[i>>6] &^= 1 << uint(i&63)
}

// Contains reports whether state i is a member of the set.
func (s StateSet) Contains(i int) bool {
	if i < 0 || i >= s.n {
		return false
	}
	return s.words[i>>6]&(1<<uint(i&63)) != 0
}

// Len returns the number of states in the set.
func (s StateSet) Len() int {
	c := 0
	for _, w := range s.words {
		c += popcount(w)
	}
	return c
}

// IsEmpty reports whether the set has no members.
func (s StateSet) IsEmpty() bool {
	for _, w := range s.words {
		if w != 0 {
			return false
		}
	}
	return true
}

// Elements returns the members of the set in ascending order.
func (s StateSet) Elements() []int {
	out := make([]int, 0, s.Len())
	for wi, w := range s.words {
		base := wi << 6
		for w != 0 {
			b := trailingZeros(w)
			out = append(out, base+b)
			w &= w - 1
		}
	}
	return out
}

// Clone returns an independent copy of the set.
func (s StateSet) Clone() StateSet {
	w := make([]uint64, len(s.words))
	copy(w, s.words)
	return StateSet{n: s.n, words: w}
}

// Equal reports whether s and t contain the same elements over the same
// universe.
func (s StateSet) Equal(t StateSet) bool {
	if s.n != t.n {
		return false
	}
	for i := range s.words {
		if s.words[i] != t.words[i] {
			return false
		}
	}
	return true
}

// Union returns the set of states in s or t. The two sets must share a
// universe; otherwise the larger universe is used and the smaller is treated as
// padded with empty words.
func (s StateSet) Union(t StateSet) StateSet {
	r := NewStateSet(maxInt(s.n, t.n))
	for i := range s.words {
		r.words[i] |= s.words[i]
	}
	for i := range t.words {
		r.words[i] |= t.words[i]
	}
	return r
}

// Intersect returns the set of states in both s and t.
func (s StateSet) Intersect(t StateSet) StateSet {
	r := NewStateSet(maxInt(s.n, t.n))
	m := minInt(len(s.words), len(t.words))
	for i := 0; i < m; i++ {
		r.words[i] = s.words[i] & t.words[i]
	}
	return r
}

// Difference returns the set of states in s but not in t.
func (s StateSet) Difference(t StateSet) StateSet {
	r := s.Clone()
	m := minInt(len(s.words), len(t.words))
	for i := 0; i < m; i++ {
		r.words[i] &^= t.words[i]
	}
	return r
}

// Complement returns the set of states of the universe not in s.
func (s StateSet) Complement() StateSet {
	r := FullStateSet(s.n)
	for i := range s.words {
		r.words[i] &^= s.words[i]
	}
	return r
}

// IsSubset reports whether every element of s is also in t.
func (s StateSet) IsSubset(t StateSet) bool {
	for wi, w := range s.words {
		var tw uint64
		if wi < len(t.words) {
			tw = t.words[wi]
		}
		if w&^tw != 0 {
			return false
		}
	}
	return true
}

// Intersects reports whether s and t share at least one element.
func (s StateSet) Intersects(t StateSet) bool {
	m := minInt(len(s.words), len(t.words))
	for i := 0; i < m; i++ {
		if s.words[i]&t.words[i] != 0 {
			return true
		}
	}
	return false
}

// Any returns some element of the set and true, or (0, false) if the set is
// empty. The element returned is the smallest for determinism.
func (s StateSet) Any() (int, bool) {
	for wi, w := range s.words {
		if w != 0 {
			return (wi << 6) + trailingZeros(w), true
		}
	}
	return 0, false
}

// String renders the set as {a, b, c} with elements in ascending order.
func (s StateSet) String() string {
	els := s.Elements()
	parts := make([]string, len(els))
	for i, e := range els {
		parts[i] = fmt.Sprintf("%d", e)
	}
	return "{" + strings.Join(parts, ", ") + "}"
}

func popcount(w uint64) int {
	c := 0
	for w != 0 {
		w &= w - 1
		c++
	}
	return c
}

func trailingZeros(w uint64) int {
	if w == 0 {
		return 64
	}
	n := 0
	for w&1 == 0 {
		w >>= 1
		n++
	}
	return n
}

func maxInt(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func minInt(a, b int) int {
	if a < b {
		return a
	}
	return b
}

// sortedUnique returns a sorted, de-duplicated copy of the input string slice.
func sortedUnique(in []string) []string {
	if len(in) == 0 {
		return nil
	}
	cp := append([]string(nil), in...)
	sort.Strings(cp)
	out := cp[:0]
	var prev string
	for i, v := range cp {
		if i == 0 || v != prev {
			out = append(out, v)
			prev = v
		}
	}
	return out
}
