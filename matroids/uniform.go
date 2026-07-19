package matroids

// UniformMatroid is the uniform matroid U(r, n): its ground set has n elements
// and a subset is independent exactly when it has at most r elements. Its rank
// function is rank(S) = min(|S|, r).
type UniformMatroid struct {
	n int
	r int
}

// NewUniformMatroid returns the uniform matroid U(r, n) on n elements with rank
// r. It panics if n < 0 or r < 0. A rank exceeding n is clamped to n.
func NewUniformMatroid(r, n int) *UniformMatroid {
	if n < 0 || r < 0 {
		panic("matroids: negative parameter for uniform matroid")
	}
	if r > n {
		r = n
	}
	return &UniformMatroid{n: n, r: r}
}

// Size returns the number of ground-set elements, n.
func (m *UniformMatroid) Size() int { return m.n }

// UniformRank returns the parameter r of U(r, n), i.e. the rank of the whole
// matroid.
func (m *UniformMatroid) UniformRank() int { return m.r }

// Rank returns min(|distinct in-range elements of set|, r).
func (m *UniformMatroid) Rank(set []int) int {
	c := 0
	seen := make(map[int]bool, len(set))
	for _, e := range set {
		if e >= 0 && e < m.n && !seen[e] {
			seen[e] = true
			c++
		}
	}
	if c < m.r {
		return c
	}
	return m.r
}

// NewFreeMatroid returns the free matroid on n elements, U(n, n), in which
// every subset is independent.
func NewFreeMatroid(n int) *UniformMatroid { return NewUniformMatroid(n, n) }

// FreeMatroid is an alias constructor returning the free matroid on n elements;
// it is provided for symmetry with the other named constructions.
func FreeMatroid(n int) *UniformMatroid { return NewFreeMatroid(n) }

// NewTrivialMatroid returns the rank-zero matroid on n elements, U(0, n), in
// which every element is a loop and the only independent set is empty.
func NewTrivialMatroid(n int) *UniformMatroid { return NewUniformMatroid(0, n) }

// TrivialMatroid is an alias constructor returning the rank-zero matroid on n
// elements.
func TrivialMatroid(n int) *UniformMatroid { return NewTrivialMatroid(n) }

// IsUniform reports whether m has the rank profile of a uniform matroid, i.e.
// every subset of size at most r is independent, where r = FullRank(m). It is
// checked by brute force and is intended for small matroids.
func IsUniform(m Matroid) bool {
	r := FullRank(m)
	for _, s := range Subsets(Ground(m)) {
		want := len(s)
		if want > r {
			want = r
		}
		if m.Rank(s) != want {
			return false
		}
	}
	return true
}

// PartitionMatroid is a matroid whose ground set is partitioned into blocks,
// each with an independence capacity: a set is independent when it contains at
// most cap[b] elements from block b. Its rank function is the sum over blocks
// of min(|S ∩ block|, cap[b]).
type PartitionMatroid struct {
	n       int
	blockOf []int // blockOf[e] is the block index of element e
	caps    []int // caps[b] is the capacity of block b
}

// NewPartitionMatroid builds a partition matroid on n elements. blockOf must
// have length n and assign each element a block index in [0, len(caps)). caps
// gives the per-block capacity (a negative capacity is treated as 0). It panics
// on a length mismatch or an out-of-range block index.
func NewPartitionMatroid(n int, blockOf, caps []int) *PartitionMatroid {
	if len(blockOf) != n {
		panic("matroids: blockOf length must equal n")
	}
	c := make([]int, len(caps))
	for i, v := range caps {
		if v < 0 {
			v = 0
		}
		c[i] = v
	}
	bo := make([]int, n)
	for e, b := range blockOf {
		if b < 0 || b >= len(caps) {
			panic("matroids: block index out of range")
		}
		bo[e] = b
	}
	return &PartitionMatroid{n: n, blockOf: bo, caps: c}
}

// NewPartitionMatroidFromBlocks builds a partition matroid from an explicit
// list of blocks (each a slice of element indices in [0, n)) and matching
// capacities. Elements not appearing in any block become loops (they belong to
// an implicit capacity-zero block). It panics on out-of-range elements or a
// blocks/caps length mismatch.
func NewPartitionMatroidFromBlocks(n int, blocks [][]int, caps []int) *PartitionMatroid {
	if len(blocks) != len(caps) {
		panic("matroids: blocks and caps length mismatch")
	}
	blockOf := make([]int, n)
	for i := range blockOf {
		blockOf[i] = len(blocks) // implicit loop block
	}
	for b, blk := range blocks {
		for _, e := range blk {
			if e < 0 || e >= n {
				panic("matroids: element out of range in block")
			}
			blockOf[e] = b
		}
	}
	fullCaps := make([]int, len(caps)+1)
	copy(fullCaps, caps)
	fullCaps[len(caps)] = 0 // loop block capacity
	return NewPartitionMatroid(n, blockOf, fullCaps)
}

// Size returns the number of ground-set elements.
func (m *PartitionMatroid) Size() int { return m.n }

// NumBlocks returns the number of blocks in the partition.
func (m *PartitionMatroid) NumBlocks() int { return len(m.caps) }

// BlockOf returns the block index of element e.
func (m *PartitionMatroid) BlockOf(e int) int { return m.blockOf[e] }

// Capacity returns the independence capacity of block b.
func (m *PartitionMatroid) Capacity(b int) int { return m.caps[b] }

// Rank returns the sum over blocks of min(count in block, capacity).
func (m *PartitionMatroid) Rank(set []int) int {
	counts := make([]int, len(m.caps))
	seen := make(map[int]bool, len(set))
	for _, e := range set {
		if e >= 0 && e < m.n && !seen[e] {
			seen[e] = true
			counts[m.blockOf[e]]++
		}
	}
	total := 0
	for b, c := range counts {
		if c > m.caps[b] {
			c = m.caps[b]
		}
		total += c
	}
	return total
}
