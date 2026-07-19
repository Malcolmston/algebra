package matroids

import "sort"

// IntSet is a mutable set of integers implemented as a map. It is used
// internally and exported as a convenience for callers who prefer set-valued
// intermediate results. The zero value is not usable; construct sets with
// [NewIntSet].
type IntSet map[int]bool

// NewIntSet returns a new [IntSet] containing the given elements. Duplicate
// arguments are collapsed.
func NewIntSet(elems ...int) IntSet {
	s := make(IntSet, len(elems))
	for _, e := range elems {
		s[e] = true
	}
	return s
}

// Add inserts e into the set.
func (s IntSet) Add(e int) { s[e] = true }

// AddAll inserts every element of elems into the set.
func (s IntSet) AddAll(elems []int) {
	for _, e := range elems {
		s[e] = true
	}
}

// Remove deletes e from the set if present.
func (s IntSet) Remove(e int) { delete(s, e) }

// Contains reports whether e is a member of the set.
func (s IntSet) Contains(e int) bool { return s[e] }

// Len returns the number of elements in the set.
func (s IntSet) Len() int { return len(s) }

// Slice returns the elements of the set in ascending order.
func (s IntSet) Slice() []int {
	out := make([]int, 0, len(s))
	for e := range s {
		out = append(out, e)
	}
	sort.Ints(out)
	return out
}

// Clone returns an independent copy of the set.
func (s IntSet) Clone() IntSet {
	c := make(IntSet, len(s))
	for e := range s {
		c[e] = true
	}
	return c
}

// Union returns a new set containing every element of s or t.
func (s IntSet) Union(t IntSet) IntSet {
	r := s.Clone()
	for e := range t {
		r[e] = true
	}
	return r
}

// Intersect returns a new set containing the elements common to s and t.
func (s IntSet) Intersect(t IntSet) IntSet {
	r := make(IntSet)
	for e := range s {
		if t[e] {
			r[e] = true
		}
	}
	return r
}

// Difference returns a new set containing the elements of s that are not in t.
func (s IntSet) Difference(t IntSet) IntSet {
	r := make(IntSet)
	for e := range s {
		if !t[e] {
			r[e] = true
		}
	}
	return r
}

// IsSubsetOf reports whether every element of s is also in t.
func (s IntSet) IsSubsetOf(t IntSet) bool {
	for e := range s {
		if !t[e] {
			return false
		}
	}
	return true
}

// Equal reports whether s and t contain exactly the same elements.
func (s IntSet) Equal(t IntSet) bool {
	if len(s) != len(t) {
		return false
	}
	return s.IsSubsetOf(t)
}

// DistinctSorted returns the distinct elements of a in ascending order.
func DistinctSorted(a []int) []int {
	if len(a) == 0 {
		return nil
	}
	seen := make(map[int]bool, len(a))
	out := make([]int, 0, len(a))
	for _, e := range a {
		if !seen[e] {
			seen[e] = true
			out = append(out, e)
		}
	}
	sort.Ints(out)
	return out
}

// Cardinality returns the number of distinct elements in a.
func Cardinality(a []int) int {
	if len(a) == 0 {
		return 0
	}
	seen := make(map[int]bool, len(a))
	for _, e := range a {
		seen[e] = true
	}
	return len(seen)
}

// SetContains reports whether a contains e.
func SetContains(a []int, e int) bool {
	for _, x := range a {
		if x == e {
			return true
		}
	}
	return false
}

// CopySet returns a copy of a.
func CopySet(a []int) []int {
	out := make([]int, len(a))
	copy(out, a)
	return out
}

// SetUnion returns the sorted union of a and b as distinct elements.
func SetUnion(a, b []int) []int {
	s := NewIntSet(a...)
	s.AddAll(b)
	return s.Slice()
}

// SetIntersection returns the sorted intersection of a and b.
func SetIntersection(a, b []int) []int {
	sb := NewIntSet(b...)
	r := make(IntSet)
	for _, e := range a {
		if sb.Contains(e) {
			r.Add(e)
		}
	}
	return r.Slice()
}

// SetDifference returns the sorted set of elements in a but not in b.
func SetDifference(a, b []int) []int {
	sb := NewIntSet(b...)
	r := make(IntSet)
	for _, e := range a {
		if !sb.Contains(e) {
			r.Add(e)
		}
	}
	return r.Slice()
}

// SetEqual reports whether a and b contain the same distinct elements.
func SetEqual(a, b []int) bool {
	return NewIntSet(a...).Equal(NewIntSet(b...))
}

// IsSubset reports whether every distinct element of a is also in b.
func IsSubset(a, b []int) bool {
	sb := NewIntSet(b...)
	for _, e := range a {
		if !sb.Contains(e) {
			return false
		}
	}
	return true
}

// SetInsert returns the sorted set a with e added.
func SetInsert(a []int, e int) []int {
	s := NewIntSet(a...)
	s.Add(e)
	return s.Slice()
}

// SetRemove returns the sorted set a with e removed.
func SetRemove(a []int, e int) []int {
	s := NewIntSet(a...)
	s.Remove(e)
	return s.Slice()
}

// Complement returns the sorted elements of {0,...,n-1} that are not in a.
func Complement(n int, a []int) []int {
	sa := NewIntSet(a...)
	out := make([]int, 0, n)
	for e := 0; e < n; e++ {
		if !sa.Contains(e) {
			out = append(out, e)
		}
	}
	return out
}

// RangeSet returns the slice {0, 1, ..., n-1}.
func RangeSet(n int) []int {
	out := make([]int, n)
	for i := range out {
		out[i] = i
	}
	return out
}

// Subsets returns every subset of the distinct elements of a. The result has
// 2^k entries, where k is the number of distinct elements, and includes the
// empty set. It panics if a has more than 63 distinct elements.
func Subsets(a []int) [][]int {
	elems := DistinctSorted(a)
	k := len(elems)
	if k > 63 {
		panic("matroids: Subsets on more than 63 elements")
	}
	out := make([][]int, 0, 1<<uint(k))
	for mask := 0; mask < (1 << uint(k)); mask++ {
		var sub []int
		for i := 0; i < k; i++ {
			if mask&(1<<uint(i)) != 0 {
				sub = append(sub, elems[i])
			}
		}
		out = append(out, sub)
	}
	return out
}

// SubsetsOfSize returns every size-k subset of the distinct elements of a.
func SubsetsOfSize(a []int, k int) [][]int {
	elems := DistinctSorted(a)
	n := len(elems)
	var out [][]int
	if k < 0 || k > n {
		return out
	}
	idx := make([]int, k)
	for i := range idx {
		idx[i] = i
	}
	for {
		sub := make([]int, k)
		for i, j := range idx {
			sub[i] = elems[j]
		}
		out = append(out, sub)
		// advance
		i := k - 1
		for i >= 0 && idx[i] == n-k+i {
			i--
		}
		if i < 0 {
			break
		}
		idx[i]++
		for j := i + 1; j < k; j++ {
			idx[j] = idx[j-1] + 1
		}
	}
	return out
}
