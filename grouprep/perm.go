package grouprep

import "sort"

// Perm is a permutation of {0,...,n-1} in one-line notation: p[i] is the image
// of i. A valid Perm of length n contains each value in [0,n) exactly once.
type Perm []int

// IdentityPerm returns the identity permutation of degree n.
func IdentityPerm(n int) Perm {
	p := make(Perm, n)
	for i := range p {
		p[i] = i
	}
	return p
}

// Degree returns the size n of the set {0,...,n-1} that p permutes.
func (p Perm) Degree() int { return len(p) }

// Apply returns the image p[i] of point i.
func (p Perm) Apply(i int) int { return p[i] }

// IsValid reports whether p is a genuine permutation of {0,...,len(p)-1}.
func (p Perm) IsValid() bool {
	seen := make([]bool, len(p))
	for _, v := range p {
		if v < 0 || v >= len(p) || seen[v] {
			return false
		}
		seen[v] = true
	}
	return true
}

// Equal reports whether p and q are the same permutation.
func (p Perm) Equal(q Perm) bool {
	if len(p) != len(q) {
		return false
	}
	for i := range p {
		if p[i] != q[i] {
			return false
		}
	}
	return true
}

// Compose returns the composition p∘q, applying q first: (p∘q)[i] = p[q[i]].
// It panics if the degrees differ.
func (p Perm) Compose(q Perm) Perm {
	if len(p) != len(q) {
		panic("grouprep: Compose requires equal degrees")
	}
	r := make(Perm, len(p))
	for i := range q {
		r[i] = p[q[i]]
	}
	return r
}

// Inverse returns the inverse permutation of p.
func (p Perm) Inverse() Perm {
	r := make(Perm, len(p))
	for i, v := range p {
		r[v] = i
	}
	return r
}

// Sign returns the sign of p: +1 if p is an even permutation and -1 if odd,
// computed from the parity of the number of transpositions.
func (p Perm) Sign() int {
	visited := make([]bool, len(p))
	sign := 1
	for i := range p {
		if visited[i] {
			continue
		}
		length := 0
		for j := i; !visited[j]; j = p[j] {
			visited[j] = true
			length++
		}
		if length%2 == 0 {
			sign = -sign
		}
	}
	return sign
}

// Order returns the multiplicative order of p, the least positive k with
// p^k = identity. It equals the least common multiple of the cycle lengths.
func (p Perm) Order() int {
	visited := make([]bool, len(p))
	order := 1
	for i := range p {
		if visited[i] {
			continue
		}
		length := 0
		for j := i; !visited[j]; j = p[j] {
			visited[j] = true
			length++
		}
		order = lcmInt(order, length)
	}
	return order
}

// CycleDecomposition returns the disjoint cycles of p as a slice of slices,
// including fixed points as length-one cycles. Each cycle begins with its
// smallest element and the cycles are ordered by that element.
func (p Perm) CycleDecomposition() [][]int {
	visited := make([]bool, len(p))
	var cycles [][]int
	for i := range p {
		if visited[i] {
			continue
		}
		var cyc []int
		for j := i; !visited[j]; j = p[j] {
			visited[j] = true
			cyc = append(cyc, j)
		}
		cycles = append(cycles, cyc)
	}
	return cycles
}

// CycleType returns the multiset of cycle lengths of p in non-increasing order,
// a partition of the degree. It is the conjugacy-class invariant for symmetric
// groups.
func (p Perm) CycleType() []int {
	var t []int
	for _, c := range p.CycleDecomposition() {
		t = append(t, len(c))
	}
	sort.Sort(sort.Reverse(sort.IntSlice(t)))
	return t
}

// PermMatrix returns the permutation matrix of p; see [PermutationMatrix].
func (p Perm) PermMatrix() Matrix {
	return PermutationMatrix([]int(p))
}

// AllPermutations returns every permutation of {0,...,n-1} in lexicographic
// order. It panics if n is negative.
func AllPermutations(n int) []Perm {
	if n < 0 {
		panic("grouprep: AllPermutations requires n >= 0")
	}
	base := IdentityPerm(n)
	var out []Perm
	var rec func(k int)
	cur := append(Perm(nil), base...)
	used := make([]bool, n)
	rec = func(k int) {
		if k == n {
			out = append(out, append(Perm(nil), cur...))
			return
		}
		for v := 0; v < n; v++ {
			if used[v] {
				continue
			}
			used[v] = true
			cur[k] = v
			rec(k + 1)
			used[v] = false
		}
	}
	rec(0)
	return out
}

// lcmInt returns the least common multiple of two non-negative integers.
func lcmInt(a, b int) int {
	if a == 0 || b == 0 {
		return 0
	}
	return a / gcdInt(a, b) * b
}

// gcdInt returns the greatest common divisor of two non-negative integers.
func gcdInt(a, b int) int {
	for b != 0 {
		a, b = b, a%b
	}
	if a < 0 {
		return -a
	}
	return a
}
