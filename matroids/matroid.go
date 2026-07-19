package matroids

import (
	"errors"
	"sort"
)

// ErrInvalidElement is returned or panicked with when a ground-set element is
// outside the range [0, Size()).
var ErrInvalidElement = errors.New("matroids: element out of range")

// ErrDimensionMismatch is returned when two matroids that must share a ground
// set have different sizes.
var ErrDimensionMismatch = errors.New("matroids: ground-set size mismatch")

// Matroid is the interface satisfied by every matroid representation in this
// package. The ground set is {0, 1, ..., Size()-1}. Rank returns the rank of a
// subset of the ground set: the size of a largest independent subset of the
// given set. Implementations must satisfy the rank axioms; [CheckRankAxioms]
// verifies them by brute force.
type Matroid interface {
	// Size returns the number of elements in the ground set.
	Size() int
	// Rank returns the rank of the given subset. Repeated or out-of-range
	// elements in the argument must be tolerated (out-of-range elements are
	// ignored, repeats collapsed).
	Rank(set []int) int
}

// Ground returns the ground set {0, 1, ..., m.Size()-1}.
func Ground(m Matroid) []int { return RangeSet(m.Size()) }

// FullRank returns the rank of the whole ground set of m, i.e. the size of a
// basis of m.
func FullRank(m Matroid) int { return m.Rank(Ground(m)) }

// Nullity returns |set| - rank(set), the number of elements of set that are
// dependent on the rest.
func Nullity(m Matroid, set []int) int {
	d := DistinctSorted(set)
	return len(d) - m.Rank(d)
}

// Corank returns the rank of the ground set minus the rank of set, i.e. how
// many more independent elements are needed to span m.
func Corank(m Matroid, set []int) int {
	return FullRank(m) - m.Rank(set)
}

// Independent reports whether set is an independent set of m, i.e. its rank
// equals its cardinality.
func Independent(m Matroid, set []int) bool {
	d := DistinctSorted(set)
	return m.Rank(d) == len(d)
}

// Dependent reports whether set is dependent (not independent) in m.
func Dependent(m Matroid, set []int) bool { return !Independent(m, set) }

// IsSpanning reports whether set spans m, i.e. its rank equals the rank of the
// whole ground set.
func IsSpanning(m Matroid, set []int) bool {
	return m.Rank(set) == FullRank(m)
}

// IsBasis reports whether set is a basis of m: a maximal independent set,
// equivalently an independent spanning set.
func IsBasis(m Matroid, set []int) bool {
	d := DistinctSorted(set)
	r := m.Rank(d)
	return r == len(d) && r == FullRank(m)
}

// IsCircuit reports whether set is a circuit of m: a minimal dependent set. A
// set C is a circuit exactly when it is dependent but every proper subset
// C\{e} is independent.
func IsCircuit(m Matroid, set []int) bool {
	c := DistinctSorted(set)
	if len(c) == 0 {
		return false
	}
	if m.Rank(c) != len(c)-1 {
		return false
	}
	for i := range c {
		sub := make([]int, 0, len(c)-1)
		sub = append(sub, c[:i]...)
		sub = append(sub, c[i+1:]...)
		if m.Rank(sub) != len(sub) {
			return false
		}
	}
	return true
}

// Spans reports whether adding e to set does not increase its rank, i.e. e is
// in the closure of set.
func Spans(m Matroid, set []int, e int) bool {
	return m.Rank(SetInsert(set, e)) == m.Rank(DistinctSorted(set))
}

// Closure returns the closure (span) of set: every element e of the ground set
// such that rank(set ∪ {e}) = rank(set). The result is sorted and always
// contains the distinct elements of set that lie in range.
func Closure(m Matroid, set []int) []int {
	d := DistinctSorted(set)
	base := m.Rank(d)
	var out []int
	for e := 0; e < m.Size(); e++ {
		if m.Rank(SetInsert(d, e)) == base {
			out = append(out, e)
		}
	}
	return out
}

// IsClosed reports whether set equals its own closure in m; such sets are the
// flats of m.
func IsClosed(m Matroid, set []int) bool {
	d := DistinctSorted(set)
	return SetEqual(d, Closure(m, d))
}

// IsFlat is a synonym for [IsClosed].
func IsFlat(m Matroid, set []int) bool { return IsClosed(m, set) }

// IsLoop reports whether element e is a loop of m, i.e. rank({e}) = 0.
func IsLoop(m Matroid, e int) bool { return m.Rank([]int{e}) == 0 }

// IsColoop reports whether element e is a coloop (isthmus, bridge) of m, i.e.
// e lies in every basis, equivalently removing e drops the rank of the ground
// set by one.
func IsColoop(m Matroid, e int) bool {
	if e < 0 || e >= m.Size() {
		return false
	}
	full := FullRank(m)
	return m.Rank(SetRemove(Ground(m), e)) == full-1
}

// Loops returns the sorted list of loops of m.
func Loops(m Matroid) []int {
	var out []int
	for e := 0; e < m.Size(); e++ {
		if IsLoop(m, e) {
			out = append(out, e)
		}
	}
	return out
}

// Coloops returns the sorted list of coloops of m.
func Coloops(m Matroid) []int {
	var out []int
	for e := 0; e < m.Size(); e++ {
		if IsColoop(m, e) {
			out = append(out, e)
		}
	}
	return out
}

// Parallel reports whether distinct non-loop elements e and f are parallel in
// m, i.e. {e, f} is a circuit (rank 1 with both dependent on each other).
func Parallel(m Matroid, e, f int) bool {
	if e == f {
		return false
	}
	if IsLoop(m, e) || IsLoop(m, f) {
		return false
	}
	return m.Rank([]int{e, f}) == 1
}

// ParallelClasses returns the partition of the non-loop elements into maximal
// classes of mutually parallel elements. Each class is sorted, and classes are
// ordered by their smallest element.
func ParallelClasses(m Matroid) [][]int {
	n := m.Size()
	assigned := make([]bool, n)
	var out [][]int
	for e := 0; e < n; e++ {
		if assigned[e] || IsLoop(m, e) {
			continue
		}
		class := []int{e}
		assigned[e] = true
		for f := e + 1; f < n; f++ {
			if !assigned[f] && Parallel(m, e, f) {
				class = append(class, f)
				assigned[f] = true
			}
		}
		out = append(out, class)
	}
	return out
}

// IsSimple reports whether m is a simple matroid: it has no loops and no two
// distinct parallel elements.
func IsSimple(m Matroid) bool {
	if len(Loops(m)) > 0 {
		return false
	}
	for _, c := range ParallelClasses(m) {
		if len(c) > 1 {
			return false
		}
	}
	return true
}

// IndependentSets returns every independent set of m, sorted by size then
// lexicographically. It enumerates all subsets and is intended for small
// ground sets.
func IndependentSets(m Matroid) [][]int {
	return filterSubsets(m, func(s []int) bool { return Independent(m, s) })
}

// DependentSets returns every dependent set of m.
func DependentSets(m Matroid) [][]int {
	return filterSubsets(m, func(s []int) bool { return Dependent(m, s) })
}

// SpanningSets returns every spanning set of m.
func SpanningSets(m Matroid) [][]int {
	return filterSubsets(m, func(s []int) bool { return IsSpanning(m, s) })
}

// Bases returns every basis of m, sorted lexicographically.
func Bases(m Matroid) [][]int {
	full := FullRank(m)
	var out [][]int
	for _, s := range SubsetsOfSize(Ground(m), full) {
		if m.Rank(s) == full {
			out = append(out, s)
		}
	}
	return out
}

// BasisCount returns the number of bases of m.
func BasisCount(m Matroid) int { return len(Bases(m)) }

// Circuits returns every circuit of m, sorted by size then lexicographically.
func Circuits(m Matroid) [][]int {
	return filterSubsets(m, func(s []int) bool { return IsCircuit(m, s) })
}

// Girth returns the size of a smallest circuit of m. If m has no circuit (m is
// a free matroid) the girth is infinite and -1 is returned.
func Girth(m Matroid) int {
	best := -1
	for _, c := range Circuits(m) {
		if best == -1 || len(c) < best {
			best = len(c)
		}
	}
	return best
}

// FundamentalCircuit returns the unique circuit contained in basis ∪ {e},
// where basis is an independent set and e is an element not in basis whose
// addition creates a dependency. If basis ∪ {e} is independent (no circuit is
// created) the result is nil.
func FundamentalCircuit(m Matroid, basis []int, e int) []int {
	b := DistinctSorted(basis)
	if SetContains(b, e) {
		return nil
	}
	be := SetInsert(b, e)
	if Independent(m, be) {
		return nil
	}
	circuit := []int{e}
	for _, x := range b {
		// x is in the fundamental circuit iff removing x from be keeps e
		// dependent, i.e. be\{x} is independent.
		rest := SetRemove(be, x)
		if Independent(m, rest) {
			circuit = append(circuit, x)
		}
	}
	sort.Ints(circuit)
	return circuit
}

// Flats returns every flat (closed set) of m, sorted by rank then
// lexicographically.
func Flats(m Matroid) [][]int {
	flats := filterSubsets(m, func(s []int) bool { return IsClosed(m, s) })
	sort.SliceStable(flats, func(i, j int) bool {
		ri, rj := m.Rank(flats[i]), m.Rank(flats[j])
		if ri != rj {
			return ri < rj
		}
		return lexLess(flats[i], flats[j])
	})
	return flats
}

// FlatsOfRank returns every flat of m whose rank equals r.
func FlatsOfRank(m Matroid, r int) [][]int {
	var out [][]int
	for _, f := range Flats(m) {
		if m.Rank(f) == r {
			out = append(out, f)
		}
	}
	return out
}

// Hyperplanes returns the hyperplanes of m: the flats of rank FullRank(m)-1.
func Hyperplanes(m Matroid) [][]int {
	full := FullRank(m)
	if full == 0 {
		return nil
	}
	return FlatsOfRank(m, full-1)
}

// RankProfile returns a slice p where p[k] is the number of independent sets of
// m of size k, for k from 0 to FullRank(m).
func RankProfile(m Matroid) []int {
	p := make([]int, FullRank(m)+1)
	for _, s := range IndependentSets(m) {
		p[len(s)]++
	}
	return p
}

// CheckRankAxioms verifies, by brute force over all subsets, that m satisfies
// the matroid rank axioms: 0 ≤ rank(X) ≤ |X|; rank is monotone; and rank is
// submodular. It returns nil if all axioms hold and a descriptive error
// otherwise. It is intended for testing small matroids.
func CheckRankAxioms(m Matroid) error {
	n := m.Size()
	subs := Subsets(Ground(m))
	for _, x := range subs {
		rx := m.Rank(x)
		if rx < 0 || rx > len(x) {
			return errors.New("matroids: rank not bounded by cardinality")
		}
		for e := 0; e < n; e++ {
			if SetContains(x, e) {
				continue
			}
			rxe := m.Rank(SetInsert(x, e))
			if rxe < rx || rxe > rx+1 {
				return errors.New("matroids: rank fails unit increase / monotonicity")
			}
		}
	}
	// Submodularity: r(A)+r(B) >= r(A∪B)+r(A∩B).
	for _, a := range subs {
		for _, b := range subs {
			lhs := m.Rank(a) + m.Rank(b)
			rhs := m.Rank(SetUnion(a, b)) + m.Rank(SetIntersection(a, b))
			if lhs < rhs {
				return errors.New("matroids: rank fails submodularity")
			}
		}
	}
	return nil
}

// filterSubsets returns all subsets of the ground set of m satisfying pred,
// ordered by size then lexicographically.
func filterSubsets(m Matroid, pred func([]int) bool) [][]int {
	all := Subsets(Ground(m))
	var out [][]int
	for _, s := range all {
		if pred(s) {
			out = append(out, s)
		}
	}
	sort.SliceStable(out, func(i, j int) bool {
		if len(out[i]) != len(out[j]) {
			return len(out[i]) < len(out[j])
		}
		return lexLess(out[i], out[j])
	})
	return out
}

// lexLess reports whether a is lexicographically smaller than b.
func lexLess(a, b []int) bool {
	for i := 0; i < len(a) && i < len(b); i++ {
		if a[i] != b[i] {
			return a[i] < b[i]
		}
	}
	return len(a) < len(b)
}
