package matroids

import "sort"

// DualMatroid is the dual of a matroid m on the same ground set. Its rank
// function is r*(S) = |S| + r(E\S) - r(E), where r is the rank function of m
// and E is the ground set. The bases of the dual are the complements of the
// bases of m.
type DualMatroid struct {
	base     Matroid
	fullRank int
}

// Dual returns the dual matroid of m. Independence, bases, circuits and so on
// of the returned matroid are the "co-" versions for m: its circuits are the
// cocircuits of m and its bases are the complements of the bases of m.
func Dual(m Matroid) *DualMatroid {
	return &DualMatroid{base: m, fullRank: FullRank(m)}
}

// Size returns the ground-set size, the same as the primal matroid.
func (d *DualMatroid) Size() int { return d.base.Size() }

// Primal returns the underlying matroid whose dual this is.
func (d *DualMatroid) Primal() Matroid { return d.base }

// Rank returns the dual rank r*(S) = |S| + r(E\S) - r(E).
func (d *DualMatroid) Rank(set []int) int {
	s := distinctInRange(set, d.base.Size())
	comp := Complement(d.base.Size(), s)
	return len(s) + d.base.Rank(comp) - d.fullRank
}

// Cocircuits returns the cocircuits of m: the circuits of its dual, equivalently
// the minimal sets meeting every basis of m.
func Cocircuits(m Matroid) [][]int { return Circuits(Dual(m)) }

// Cobases returns the cobases of m: the bases of its dual, i.e. the complements
// of the bases of m.
func Cobases(m Matroid) [][]int {
	n := m.Size()
	var out [][]int
	for _, b := range Bases(m) {
		out = append(out, Complement(n, b))
	}
	sort.SliceStable(out, func(i, j int) bool { return lexLess(out[i], out[j]) })
	return out
}

// FundamentalCocircuit returns the fundamental cocircuit of element e with
// respect to basis: the unique cocircuit contained in (E \ basis) ∪ {e}, where
// e is an element of basis. It is the fundamental circuit of e in the dual.
func FundamentalCocircuit(m Matroid, basis []int, e int) []int {
	d := Dual(m)
	cobasis := Complement(m.Size(), basis)
	return FundamentalCircuit(d, cobasis, e)
}

// MinorMatroid is a minor of a base matroid, obtained by deleting a set D and
// contracting a set C of ground elements. Its ground set is E \ (D ∪ C),
// re-indexed to {0, ..., k-1}; [MinorMatroid.Elements] recovers the original
// labels. Its rank function is r'(T) = r(map(T) ∪ C) - r(C).
type MinorMatroid struct {
	base     Matroid
	keep     []int // keep[i] is the original label of new element i
	contract []int
	rankC    int
}

// Minor returns the minor of m obtained by deleting the elements in del and
// contracting the elements in con. del and con must be disjoint; any overlap
// panics. The remaining elements are re-indexed to {0, ..., k-1} in ascending
// order of their original labels.
func Minor(m Matroid, del, con []int) *MinorMatroid {
	d := NewIntSet(del...)
	c := NewIntSet(con...)
	for e := range c {
		if d.Contains(e) {
			panic("matroids: deletion and contraction sets overlap")
		}
	}
	var keep []int
	for e := 0; e < m.Size(); e++ {
		if !d.Contains(e) && !c.Contains(e) {
			keep = append(keep, e)
		}
	}
	conSlice := c.Slice()
	return &MinorMatroid{
		base:     m,
		keep:     keep,
		contract: conSlice,
		rankC:    m.Rank(conSlice),
	}
}

// Deletion returns the minor of m obtained by deleting the elements in del.
func Deletion(m Matroid, del []int) *MinorMatroid { return Minor(m, del, nil) }

// Contraction returns the minor of m obtained by contracting the elements in
// con.
func Contraction(m Matroid, con []int) *MinorMatroid { return Minor(m, nil, con) }

// Restriction returns the restriction of m to the elements in keep, i.e. the
// deletion of all other elements.
func Restriction(m Matroid, keep []int) *MinorMatroid {
	del := Complement(m.Size(), keep)
	return Minor(m, del, nil)
}

// Size returns the number of elements of the minor.
func (mm *MinorMatroid) Size() int { return len(mm.keep) }

// Elements returns the original ground-set labels of the minor's elements, in
// new-index order (Elements()[i] is the original label of new element i).
func (mm *MinorMatroid) Elements() []int {
	out := make([]int, len(mm.keep))
	copy(out, mm.keep)
	return out
}

// OriginalLabel maps a new element index of the minor back to its original
// label in the base matroid.
func (mm *MinorMatroid) OriginalLabel(i int) int { return mm.keep[i] }

// Rank returns the rank of the given set of new-indexed elements in the minor.
func (mm *MinorMatroid) Rank(set []int) int {
	mapped := make([]int, 0, len(set)+len(mm.contract))
	seen := make(map[int]bool, len(set))
	for _, e := range set {
		if e >= 0 && e < len(mm.keep) && !seen[e] {
			seen[e] = true
			mapped = append(mapped, mm.keep[e])
		}
	}
	mapped = append(mapped, mm.contract...)
	return mm.base.Rank(mapped) - mm.rankC
}

// DirectSumMatroid is the direct sum (1-sum) of two or more matroids on
// disjoint ground sets. The ground set is the concatenation of the parts: part
// p occupies a contiguous block of indices. Its rank is the sum of the ranks of
// the restrictions to each block.
type DirectSumMatroid struct {
	parts   []Matroid
	offsets []int // offsets[p] is the starting index of part p
	total   int
}

// DirectSum returns the direct sum of the given matroids. The ground set of the
// result is the disjoint union of the parts, with part 0 first. With no
// arguments it returns an empty matroid.
func DirectSum(parts ...Matroid) *DirectSumMatroid {
	offsets := make([]int, len(parts))
	total := 0
	for i, p := range parts {
		offsets[i] = total
		total += p.Size()
	}
	cp := make([]Matroid, len(parts))
	copy(cp, parts)
	return &DirectSumMatroid{parts: cp, offsets: offsets, total: total}
}

// Size returns the total number of elements across all parts.
func (ds *DirectSumMatroid) Size() int { return ds.total }

// NumParts returns the number of summands.
func (ds *DirectSumMatroid) NumParts() int { return len(ds.parts) }

// PartOf returns the part index and the local (within-part) index of global
// element e.
func (ds *DirectSumMatroid) PartOf(e int) (part, local int) {
	for p := len(ds.parts) - 1; p >= 0; p-- {
		if e >= ds.offsets[p] {
			return p, e - ds.offsets[p]
		}
	}
	return 0, e
}

// Rank returns the sum over parts of the rank of the elements of set falling in
// that part (translated to local indices).
func (ds *DirectSumMatroid) Rank(set []int) int {
	local := make([][]int, len(ds.parts))
	seen := make(map[int]bool, len(set))
	for _, e := range set {
		if e < 0 || e >= ds.total || seen[e] {
			continue
		}
		seen[e] = true
		p, l := ds.PartOf(e)
		local[p] = append(local[p], l)
	}
	total := 0
	for p, part := range ds.parts {
		total += part.Rank(local[p])
	}
	return total
}
