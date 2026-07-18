package discretemath

// Set is a mathematical set of distinct comparable elements, represented as a
// map whose keys are the members. The zero value is not usable; construct sets
// with NewSet.
type Set[T comparable] map[T]struct{}

// Pair is an ordered two-element tuple, used as the element type of a Cartesian
// product of two (possibly heterogeneous) slices.
type Pair[T, U any] struct {
	// First is the component drawn from the first operand.
	First T
	// Second is the component drawn from the second operand.
	Second U
}

// NewSet returns a new Set containing the given items, with duplicates removed.
func NewSet[T comparable](items ...T) Set[T] {
	s := make(Set[T], len(items))
	for _, it := range items {
		s[it] = struct{}{}
	}
	return s
}

// Add inserts the given items into the set. Items already present are ignored.
func (s Set[T]) Add(items ...T) {
	for _, it := range items {
		s[it] = struct{}{}
	}
}

// Remove deletes the given items from the set. Items not present are ignored.
func (s Set[T]) Remove(items ...T) {
	for _, it := range items {
		delete(s, it)
	}
}

// Contains reports whether item is a member of the set.
func (s Set[T]) Contains(item T) bool {
	_, ok := s[item]
	return ok
}

// Len returns the number of elements in the set (its cardinality).
func (s Set[T]) Len() int {
	return len(s)
}

// Elements returns the members of the set as a slice in unspecified order.
func (s Set[T]) Elements() []T {
	out := make([]T, 0, len(s))
	for k := range s {
		out = append(out, k)
	}
	return out
}

// Clone returns an independent copy of the set.
func (s Set[T]) Clone() Set[T] {
	out := make(Set[T], len(s))
	for k := range s {
		out[k] = struct{}{}
	}
	return out
}

// Equal reports whether the set has exactly the same members as other.
func (s Set[T]) Equal(other Set[T]) bool {
	if len(s) != len(other) {
		return false
	}
	for k := range s {
		if _, ok := other[k]; !ok {
			return false
		}
	}
	return true
}

// IsSubset reports whether every element of the receiver is also in other.
func (s Set[T]) IsSubset(other Set[T]) bool {
	if len(s) > len(other) {
		return false
	}
	for k := range s {
		if _, ok := other[k]; !ok {
			return false
		}
	}
	return true
}

// IsSuperset reports whether the receiver contains every element of other.
func (s Set[T]) IsSuperset(other Set[T]) bool {
	return other.IsSubset(s)
}

// IsDisjoint reports whether the receiver and other share no common element.
func (s Set[T]) IsDisjoint(other Set[T]) bool {
	small, large := s, other
	if len(large) < len(small) {
		small, large = large, small
	}
	for k := range small {
		if _, ok := large[k]; ok {
			return false
		}
	}
	return true
}

// Union returns a new set containing every element that appears in any of the
// given sets.
func Union[T comparable](sets ...Set[T]) Set[T] {
	out := make(Set[T])
	for _, s := range sets {
		for k := range s {
			out[k] = struct{}{}
		}
	}
	return out
}

// Intersection returns a new set containing exactly the elements common to all
// of the given sets. With no arguments it returns an empty set.
func Intersection[T comparable](sets ...Set[T]) Set[T] {
	out := make(Set[T])
	if len(sets) == 0 {
		return out
	}
	// Iterate the smallest set for efficiency.
	smallest := 0
	for i := 1; i < len(sets); i++ {
		if len(sets[i]) < len(sets[smallest]) {
			smallest = i
		}
	}
	for k := range sets[smallest] {
		inAll := true
		for i, s := range sets {
			if i == smallest {
				continue
			}
			if _, ok := s[k]; !ok {
				inAll = false
				break
			}
		}
		if inAll {
			out[k] = struct{}{}
		}
	}
	return out
}

// Difference returns a new set containing the elements of a that are not in b
// (the relative complement a \ b).
func Difference[T comparable](a, b Set[T]) Set[T] {
	out := make(Set[T])
	for k := range a {
		if _, ok := b[k]; !ok {
			out[k] = struct{}{}
		}
	}
	return out
}

// SymmetricDifference returns a new set containing the elements that are in
// exactly one of a and b.
func SymmetricDifference[T comparable](a, b Set[T]) Set[T] {
	out := make(Set[T])
	for k := range a {
		if _, ok := b[k]; !ok {
			out[k] = struct{}{}
		}
	}
	for k := range b {
		if _, ok := a[k]; !ok {
			out[k] = struct{}{}
		}
	}
	return out
}

// PowerSet returns every subset of the given elements as a slice of slices. The
// input is treated as an ordered collection; each output subset preserves the
// relative order of the input, and the empty subset is included. The result has
// 2**len(elements) entries, so the input length must be modest.
func PowerSet[T any](elements []T) [][]T {
	n := len(elements)
	total := 1 << uint(n)
	out := make([][]T, 0, total)
	for mask := 0; mask < total; mask++ {
		subset := make([]T, 0, PopCount(uint64(mask)))
		for i := 0; i < n; i++ {
			if mask&(1<<uint(i)) != 0 {
				subset = append(subset, elements[i])
			}
		}
		out = append(out, subset)
	}
	return out
}

// CartesianProduct returns the Cartesian product of a and b: every ordered pair
// (x, y) with x drawn from a and y from b. The pairs are ordered with a's index
// varying slowest.
func CartesianProduct[T, U any](a []T, b []U) []Pair[T, U] {
	out := make([]Pair[T, U], 0, len(a)*len(b))
	for _, x := range a {
		for _, y := range b {
			out = append(out, Pair[T, U]{First: x, Second: y})
		}
	}
	return out
}

// CartesianProductN returns the Cartesian product of any number of homogeneous
// slices as a slice of tuples, each tuple having one element from each input
// slice in order. The last input varies fastest. If any input slice is empty the
// product is empty; with no inputs the result contains a single empty tuple.
func CartesianProductN[T any](sets ...[]T) [][]T {
	out := [][]T{{}}
	for _, set := range sets {
		next := make([][]T, 0, len(out)*len(set))
		for _, prefix := range out {
			for _, v := range set {
				tuple := make([]T, len(prefix)+1)
				copy(tuple, prefix)
				tuple[len(prefix)] = v
				next = append(next, tuple)
			}
		}
		out = next
	}
	return out
}
