package matroids

// UnionMatroid is the matroid union M1 ∨ M2 ∨ ... ∨ Mk of several matroids on a
// common ground set. A set S is independent in the union exactly when it can be
// partitioned into sets I1, ..., Ik with Ii independent in Mi. Its rank is
// computed by the matroid-union theorem, reduced to a single matroid
// intersection between the direct sum of the parts (on k disjoint copies of the
// ground set) and a partition matroid that allows each ground element in at
// most one copy.
type UnionMatroid struct {
	parts []Matroid
	n     int
	sum   *DirectSumMatroid
}

// Union returns the matroid union of the given matroids, which must all share a
// ground set (equal Size). It panics on a size mismatch. With no arguments it
// returns the rank-zero matroid on zero elements.
func Union(parts ...Matroid) *UnionMatroid {
	n := 0
	if len(parts) > 0 {
		n = parts[0].Size()
		for _, p := range parts {
			if p.Size() != n {
				panic(ErrDimensionMismatch)
			}
		}
	}
	cp := make([]Matroid, len(parts))
	copy(cp, parts)
	return &UnionMatroid{parts: cp, n: n, sum: DirectSum(cp...)}
}

// Size returns the size of the shared ground set.
func (u *UnionMatroid) Size() int { return u.n }

// NumParts returns the number of matroids being combined.
func (u *UnionMatroid) NumParts() int { return len(u.parts) }

// Rank returns the union rank of set: the maximum size of a subset of set that
// is partitionable into independent sets of the parts.
func (u *UnionMatroid) Rank(set []int) int {
	if len(u.parts) == 0 || u.n == 0 {
		return 0
	}
	allowed := make([]bool, u.n)
	for _, e := range distinctInRange(set, u.n) {
		allowed[e] = true
	}
	common := u.intersectReduction(allowed)
	return len(common)
}

// MaxIndependentSet returns a maximum-size subset of the whole ground set that
// is independent in the union (partitionable into independents of the parts).
func (u *UnionMatroid) MaxIndependentSet() []int {
	if len(u.parts) == 0 || u.n == 0 {
		return nil
	}
	allowed := make([]bool, u.n)
	for e := range allowed {
		allowed[e] = true
	}
	common := u.intersectReduction(allowed)
	res := NewIntSet()
	for _, idx := range common {
		res.Add(idx % u.n)
	}
	return res.Slice()
}

// Partition returns an assignment of a maximum union-independent set into the
// parts: result[i] is the (sorted) list of ground elements assigned to matroid
// i, each list independent in its matroid, all lists pairwise disjoint. Their
// union is a maximum union-independent set.
func (u *UnionMatroid) Partition() [][]int {
	out := make([][]int, len(u.parts))
	if len(u.parts) == 0 || u.n == 0 {
		return out
	}
	allowed := make([]bool, u.n)
	for e := range allowed {
		allowed[e] = true
	}
	common := u.intersectReduction(allowed)
	for _, idx := range common {
		part := idx / u.n
		elem := idx % u.n
		out[part] = append(out[part], elem)
	}
	for i := range out {
		if out[i] != nil {
			out[i] = DistinctSorted(out[i])
		}
	}
	return out
}

// intersectReduction runs the reduction to matroid intersection restricted to
// the ground elements marked in allowed. It returns the maximum common
// independent set (as global indices into the k·n disjoint copies).
func (u *UnionMatroid) intersectReduction(allowed []bool) []int {
	k := len(u.parts)
	N := k * u.n
	blockOf := make([]int, N)
	for i := 0; i < k; i++ {
		for e := 0; e < u.n; e++ {
			blockOf[i*u.n+e] = e
		}
	}
	caps := make([]int, u.n)
	for e := 0; e < u.n; e++ {
		if allowed[e] {
			caps[e] = 1
		}
	}
	b := NewPartitionMatroid(N, blockOf, caps)
	return Intersection(u.sum, b)
}

// UnionMaxIndependentSet returns a maximum-size set independent in the union of
// the given matroids.
func UnionMaxIndependentSet(parts []Matroid) []int {
	return Union(parts...).MaxIndependentSet()
}

// UnionRank returns the rank of the union of the given matroids (the maximum
// size of a set partitionable into independents of the parts).
func UnionRank(parts []Matroid) int {
	u := Union(parts...)
	return u.Rank(Ground(u))
}

// IsPartitionable reports whether the whole ground set of the given matroids
// can be partitioned into independent sets, one per matroid; equivalently the
// union rank equals the ground-set size.
func IsPartitionable(parts []Matroid) bool {
	if len(parts) == 0 {
		return true
	}
	u := Union(parts...)
	return u.Rank(Ground(u)) == u.Size()
}
