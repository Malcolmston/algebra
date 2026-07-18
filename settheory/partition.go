package settheory

import "sort"

// Partition is a partition of a finite set of integers into non-empty,
// pairwise-disjoint blocks whose union is the whole set. Blocks are represented
// as IntSets. Canonical partitions produced by this package list their blocks in
// ascending order of least element.
type Partition []IntSet

// Blocks returns the blocks of the partition sorted by ascending least element,
// each block itself presented via IntSet.Elements when rendered. The returned
// slice is freshly ordered but shares the underlying IntSet values.
func (p Partition) Blocks() []IntSet {
	out := make([]IntSet, len(p))
	copy(out, p)
	sort.Slice(out, func(i, j int) bool {
		mi, _ := out[i].Min()
		mj, _ := out[j].Min()
		return mi < mj
	})
	return out
}

// Len returns the number of blocks in the partition.
func (p Partition) Len() int {
	return len(p)
}

// Universe returns the union of all blocks, i.e. the set being partitioned.
func (p Partition) Universe() IntSet {
	u := make(IntSet)
	for _, b := range p {
		for e := range b {
			u[e] = struct{}{}
		}
	}
	return u
}

// IsValidPartitionOf reports whether p is a valid partition of set: every block
// is non-empty, blocks are pairwise disjoint, and their union equals set.
func (p Partition) IsValidPartitionOf(set IntSet) bool {
	seen := make(IntSet)
	total := 0
	for _, b := range p {
		if b.IsEmpty() {
			return false
		}
		for e := range b {
			if _, ok := seen[e]; ok {
				return false // overlap between blocks
			}
			if _, ok := set[e]; !ok {
				return false // element outside the universe
			}
			seen[e] = struct{}{}
			total++
		}
	}
	return total == len(set)
}

// EquivalenceClasses returns the partition of set into the equivalence classes
// of the relation rel, which must be an equivalence relation on set. The blocks
// are canonicalized in ascending order of least element. It panics if rel is not
// an equivalence relation on set.
func EquivalenceClasses(rel Relation, set IntSet) Partition {
	if !rel.IsEquivalenceOn(set) {
		panic("settheory: EquivalenceClasses called with a non-equivalence relation")
	}
	elems := set.Elements()
	assigned := make(map[int]int) // element -> block index
	var blocks []IntSet
	for _, x := range elems {
		if _, ok := assigned[x]; ok {
			continue
		}
		idx := len(blocks)
		block := make(IntSet)
		for _, y := range elems {
			if rel.Related(x, y) {
				block[y] = struct{}{}
				assigned[y] = idx
			}
		}
		blocks = append(blocks, block)
	}
	return Partition(blocks).Blocks()
}

// RelationFromPartition returns the equivalence relation induced by p: two
// elements are related exactly when they lie in the same block.
func RelationFromPartition(p Partition) Relation {
	r := make(Relation)
	for _, b := range p {
		elems := b.Elements()
		for _, a := range elems {
			for _, c := range elems {
				r[Pair{From: a, To: c}] = struct{}{}
			}
		}
	}
	return r
}

// Refines reports whether p is a refinement of q, i.e. every block of p is
// contained in some block of q. Both partitions are assumed to cover the same
// universe.
func (p Partition) Refines(q Partition) bool {
	for _, bp := range p {
		contained := false
		for _, bq := range q {
			if bp.IsSubsetOf(bq) {
				contained = true
				break
			}
		}
		if !contained {
			return false
		}
	}
	return true
}

// BellNumber returns the nth Bell number, the number of partitions of an n-set,
// computed exactly with the Bell triangle. It panics for negative n. Values are
// exact for n up to about 25 before int overflow becomes a concern on 64-bit
// platforms.
func BellNumber(n int) int {
	if n < 0 {
		panic("settheory: BellNumber called with negative n")
	}
	if n == 0 {
		return 1
	}
	prev := []int{1}
	for row := 1; row <= n; row++ {
		cur := make([]int, row+1)
		cur[0] = prev[len(prev)-1]
		for i := 1; i <= row; i++ {
			cur[i] = cur[i-1] + prev[i-1]
		}
		prev = cur
	}
	return prev[0]
}

// StirlingSecond returns the Stirling number of the second kind S(n, k), the
// number of ways to partition an n-set into exactly k non-empty blocks. It
// panics for negative arguments and returns zero when k exceeds n. The recurrence
// S(n, k) = k*S(n-1, k) + S(n-1, k-1) is used.
func StirlingSecond(n, k int) int {
	if n < 0 || k < 0 {
		panic("settheory: StirlingSecond called with negative argument")
	}
	if k == 0 {
		if n == 0 {
			return 1
		}
		return 0
	}
	if k > n {
		return 0
	}
	// dp[j] = S(i, j) for the current row i.
	dp := make([]int, k+1)
	dp[0] = 1 // S(0,0)=1; S(0,j)=0 for j>0
	for i := 1; i <= n; i++ {
		for j := min(i, k); j >= 1; j-- {
			dp[j] = j*dp[j] + dp[j-1]
		}
		dp[0] = 0 // S(i,0)=0 for i>=1
	}
	return dp[k]
}
