package exterioralgebra

import "sort"

// Popcount returns the number of set bits in mask, i.e. the grade of the basis
// blade that mask encodes.
func Popcount(mask uint) int {
	c := 0
	for mask != 0 {
		mask &= mask - 1
		c++
	}
	return c
}

// FullMask returns the bitmask 0b1…1 with the low n bits set, encoding the
// top-grade basis blade e_0∧…∧e_{n-1} of Λ(Rⁿ).
func FullMask(n int) uint {
	if n <= 0 {
		return 0
	}
	return (uint(1) << uint(n)) - 1
}

// MaskToIndices returns the sorted slice of basis indices whose bits are set in
// mask. The result is always in strictly increasing order.
func MaskToIndices(mask uint) []int {
	var idx []int
	for i := 0; mask != 0; i++ {
		if mask&1 == 1 {
			idx = append(idx, i)
		}
		mask >>= 1
	}
	return idx
}

// IndicesToMask converts an ordered list of basis indices into the bitmask of
// the corresponding blade together with the sign of the permutation that sorts
// the indices into increasing order. A repeated index yields a zero sign
// (the blade is zero) and ok == false; a negative index or one that is not less
// than n also yields ok == false.
//
// When ok is true the returned mask encodes the sorted blade and sign is ±1,
// so that e_{idx[0]}∧…∧e_{idx[k-1]} = sign · e_mask.
func IndicesToMask(n int, idx ...int) (mask uint, sign int, ok bool) {
	for _, i := range idx {
		if i < 0 || i >= n {
			return 0, 0, false
		}
	}
	// Count inversions to get the permutation parity, and detect repeats.
	sign = 1
	work := append([]int(nil), idx...)
	for i := 1; i < len(work); i++ {
		for j := i; j > 0 && work[j-1] > work[j]; j-- {
			work[j-1], work[j] = work[j], work[j-1]
			sign = -sign
		}
	}
	for i := 1; i < len(work); i++ {
		if work[i-1] == work[i] {
			return 0, 0, false
		}
	}
	for _, i := range work {
		mask |= uint(1) << uint(i)
	}
	return mask, sign, true
}

// reorderSign returns the sign such that e_a ∧ e_b = reorderSign(a,b) · e_{a|b}
// for disjoint blade masks a and b. It is +1 or -1 according to the parity of
// the number of transpositions needed to merge the two increasing index lists,
// and is meaningful only when a&b == 0.
func reorderSign(a, b uint) int {
	a >>= 1
	swaps := 0
	for a != 0 {
		swaps += Popcount(a & b)
		a >>= 1
	}
	if swaps&1 == 0 {
		return 1
	}
	return -1
}

// sortedMasks returns the masks of m sorted first by grade (population count)
// and then, within a grade, by the numeric value of the mask. It gives Forms a
// stable, human-readable term order.
func sortedMasks(m map[uint]float64) []uint {
	ks := make([]uint, 0, len(m))
	for k := range m {
		ks = append(ks, k)
	}
	sort.Slice(ks, func(i, j int) bool {
		gi, gj := Popcount(ks[i]), Popcount(ks[j])
		if gi != gj {
			return gi < gj
		}
		return ks[i] < ks[j]
	})
	return ks
}
