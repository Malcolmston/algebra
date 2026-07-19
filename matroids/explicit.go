package matroids

import (
	"errors"
	"sort"
)

// ExplicitMatroid is a matroid given directly by the collection of its
// independent sets. Membership is tested against a canonical bit-mask
// representation, so ranks and all derived notions are exact. Ground-set size
// is limited to 62 elements by the bit-mask encoding.
type ExplicitMatroid struct {
	n     int
	indep map[uint64]bool
}

// maxExplicitElements is the largest ground-set size supported by the bit-mask
// encoding used by ExplicitMatroid.
const maxExplicitElements = 62

func maskOf(set []int, n int) uint64 {
	var mask uint64
	for _, e := range set {
		if e >= 0 && e < n {
			mask |= 1 << uint(e)
		}
	}
	return mask
}

// NewExplicitMatroidFromIndependentSets builds a matroid on n elements whose
// independent sets are exactly the given sets, downward-closed automatically:
// every subset of a listed set is also made independent, and the empty set is
// always independent. It does not verify the exchange axiom; use
// [ExplicitMatroid.Validate] or [CheckRankAxioms] for that. It panics if
// n > 62.
func NewExplicitMatroidFromIndependentSets(n int, sets [][]int) *ExplicitMatroid {
	if n > maxExplicitElements {
		panic("matroids: explicit matroid supports at most 62 elements")
	}
	indep := map[uint64]bool{0: true}
	for _, s := range sets {
		mask := maskOf(s, n)
		// add mask and all its subsets
		addSubsets(indep, mask)
	}
	return &ExplicitMatroid{n: n, indep: indep}
}

// addSubsets records mask and every subset of mask as independent.
func addSubsets(indep map[uint64]bool, mask uint64) {
	// enumerate submasks of mask
	sub := mask
	for {
		indep[sub] = true
		if sub == 0 {
			break
		}
		sub = (sub - 1) & mask
	}
}

// NewExplicitMatroidFromBases builds a matroid on n elements whose bases are
// exactly the given sets. A set is independent when it is a subset of some
// basis. It panics if n > 62 or if the bases are not all the same size.
func NewExplicitMatroidFromBases(n int, bases [][]int) *ExplicitMatroid {
	if n > maxExplicitElements {
		panic("matroids: explicit matroid supports at most 62 elements")
	}
	size := -1
	indep := map[uint64]bool{0: true}
	for _, b := range bases {
		d := DistinctSorted(b)
		if size == -1 {
			size = len(d)
		} else if len(d) != size {
			panic("matroids: bases of unequal size")
		}
		addSubsets(indep, maskOf(d, n))
	}
	return &ExplicitMatroid{n: n, indep: indep}
}

// NewExplicitMatroidFromCircuits builds a matroid on n elements from its
// circuits: a set is independent exactly when it contains no circuit. It panics
// if n > 62.
func NewExplicitMatroidFromCircuits(n int, circuits [][]int) *ExplicitMatroid {
	if n > maxExplicitElements {
		panic("matroids: explicit matroid supports at most 62 elements")
	}
	circMasks := make([]uint64, 0, len(circuits))
	for _, c := range circuits {
		circMasks = append(circMasks, maskOf(c, n))
	}
	indep := make(map[uint64]bool)
	full := (uint64(1) << uint(n)) - 1
	// enumerate every subset; mark independent if it contains no circuit
	for mask := uint64(0); mask <= full; mask++ {
		ok := true
		for _, cm := range circMasks {
			if cm != 0 && mask&cm == cm {
				ok = false
				break
			}
		}
		if ok {
			indep[mask] = true
		}
		if mask == full {
			break
		}
	}
	return &ExplicitMatroid{n: n, indep: indep}
}

// Size returns the number of ground-set elements.
func (m *ExplicitMatroid) Size() int { return m.n }

// IsIndependentSet reports whether the exact set (as a mask) is recorded as
// independent.
func (m *ExplicitMatroid) IsIndependentSet(set []int) bool {
	return m.indep[maskOf(set, m.n)]
}

// Rank returns the size of a largest independent subset of set, found by a
// downward search over the recorded independent sets.
func (m *ExplicitMatroid) Rank(set []int) int {
	mask := maskOf(set, m.n)
	if m.indep[mask] {
		return popcount(mask)
	}
	// The rank of set equals the size of the largest recorded-independent
	// submask of mask. Because independence is downward closed, a greedy
	// removal of single elements finds a maximal independent subset whose size
	// is the rank.
	best := 0
	// breadth over submasks by size using greedy augmentation from empty set.
	cur := uint64(0)
	elems := maskBits(mask, m.n)
	for _, e := range elems {
		trial := cur | (1 << uint(e))
		if m.indep[trial] {
			cur = trial
		}
	}
	best = popcount(cur)
	return best
}

// ListIndependentSets returns every recorded independent set of m, sorted by
// size then lexicographically.
func (m *ExplicitMatroid) ListIndependentSets() [][]int {
	out := make([][]int, 0, len(m.indep))
	for mask := range m.indep {
		out = append(out, maskBits(mask, m.n))
	}
	sort.SliceStable(out, func(i, j int) bool {
		if len(out[i]) != len(out[j]) {
			return len(out[i]) < len(out[j])
		}
		return lexLess(out[i], out[j])
	})
	return out
}

// Validate checks that the recorded independent sets form a matroid: the empty
// set is independent, independence is downward closed, and the exchange axiom
// holds. It returns nil on success or a descriptive error.
func (m *ExplicitMatroid) Validate() error {
	if !m.indep[0] {
		return errors.New("matroids: empty set not independent")
	}
	// downward closure
	for mask := range m.indep {
		bits := maskBits(mask, m.n)
		for _, e := range bits {
			if !m.indep[mask&^(1<<uint(e))] {
				return errors.New("matroids: independence not downward closed")
			}
		}
	}
	// exchange axiom: if |A| < |B| both independent, some e in B\A extends A
	sets := m.ListIndependentSets()
	for _, a := range sets {
		for _, b := range sets {
			if len(a) >= len(b) {
				continue
			}
			am := maskOf(a, m.n)
			bm := maskOf(b, m.n)
			extendable := false
			diff := bm &^ am
			for _, e := range maskBits(diff, m.n) {
				if m.indep[am|(1<<uint(e))] {
					extendable = true
					break
				}
			}
			if !extendable {
				return errors.New("matroids: exchange axiom violated")
			}
		}
	}
	return nil
}

// popcount returns the number of set bits in mask.
func popcount(mask uint64) int {
	c := 0
	for mask != 0 {
		mask &= mask - 1
		c++
	}
	return c
}

// maskBits returns the element indices set in mask, in ascending order.
func maskBits(mask uint64, n int) []int {
	var out []int
	for e := 0; e < n; e++ {
		if mask&(1<<uint(e)) != 0 {
			out = append(out, e)
		}
	}
	return out
}

// IsMatroidIndependenceSystem reports whether the given family of subsets of
// {0..n-1} satisfies the matroid axioms (non-empty, downward closed, exchange).
// It builds an [ExplicitMatroid] from the exact sets and validates it. Note the
// family must already be downward closed to pass, since it is taken verbatim.
func IsMatroidIndependenceSystem(n int, sets [][]int) bool {
	if n > maxExplicitElements {
		return false
	}
	indep := make(map[uint64]bool)
	for _, s := range sets {
		indep[maskOf(s, n)] = true
	}
	m := &ExplicitMatroid{n: n, indep: indep}
	return m.Validate() == nil
}
