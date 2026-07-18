package groups

import (
	"fmt"
	"strings"
)

// Perm is a permutation of {0, 1, ..., n-1} in one-line notation: p[i] is the
// image of i under p. A valid Perm of length n contains each value in [0, n)
// exactly once.
type Perm []int

// Identity returns the identity permutation of degree n, mapping every i to
// itself. The degree n must be non-negative.
func Identity(n int) Perm {
	if n < 0 {
		panic("groups: Identity requires n >= 0")
	}
	p := make(Perm, n)
	for i := range p {
		p[i] = i
	}
	return p
}

// IsPermutation reports whether p is a valid permutation: a bijection of
// {0, ..., len(p)-1} with every value appearing exactly once.
func IsPermutation(p Perm) bool {
	n := len(p)
	seen := make([]bool, n)
	for _, v := range p {
		if v < 0 || v >= n || seen[v] {
			return false
		}
		seen[v] = true
	}
	return true
}

// Degree returns the size n of the set {0, ..., n-1} that p permutes, i.e.
// len(p).
func (p Perm) Degree() int {
	return len(p)
}

// Apply returns the image p[i] of the point i under the permutation p. It
// panics if i is out of range.
func (p Perm) Apply(i int) int {
	return p[i]
}

// Compose returns the composition p∘q, the permutation that applies q first
// and then p, so (Compose(p, q))[i] == p[q[i]]. Both permutations must have
// the same degree.
func Compose(p, q Perm) Perm {
	if len(p) != len(q) {
		panic("groups: Compose requires equal degrees")
	}
	r := make(Perm, len(p))
	for i := range q {
		r[i] = p[q[i]]
	}
	return r
}

// Inverse returns the inverse permutation p^-1, satisfying
// Compose(p, p.Inverse()) == Identity(p.Degree()).
func (p Perm) Inverse() Perm {
	inv := make(Perm, len(p))
	for i, v := range p {
		inv[v] = i
	}
	return inv
}

// Equal reports whether p and q are the same permutation (equal degree and
// equal images).
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

// IsIdentity reports whether p fixes every point, i.e. p[i] == i for all i.
func (p Perm) IsIdentity() bool {
	for i, v := range p {
		if v != i {
			return false
		}
	}
	return true
}

// Cycles returns the disjoint-cycle decomposition of p. Each cycle is a slice
// of points listed starting from its smallest element, and the cycles are
// ordered by their smallest element. Fixed points appear as length-one cycles.
func (p Perm) Cycles() [][]int {
	n := len(p)
	seen := make([]bool, n)
	var cycles [][]int
	for i := 0; i < n; i++ {
		if seen[i] {
			continue
		}
		var cyc []int
		j := i
		for !seen[j] {
			seen[j] = true
			cyc = append(cyc, j)
			j = p[j]
		}
		cycles = append(cycles, cyc)
	}
	return cycles
}

// CycleType returns the multiset of cycle lengths of p sorted in
// non-increasing order (a partition of the degree n). For example the cycle
// type of a product of a 3-cycle and a transposition in S5 is [3, 2].
func (p Perm) CycleType() []int {
	cycles := p.Cycles()
	lengths := make([]int, len(cycles))
	for i, c := range cycles {
		lengths[i] = len(c)
	}
	// Insertion sort into non-increasing order (small slices).
	for i := 1; i < len(lengths); i++ {
		v := lengths[i]
		j := i - 1
		for j >= 0 && lengths[j] < v {
			lengths[j+1] = lengths[j]
			j--
		}
		lengths[j+1] = v
	}
	return lengths
}

// Order returns the order of p in the symmetric group: the smallest k > 0 with
// p^k the identity. It equals the least common multiple of the cycle lengths.
// The identity permutation has order 1.
func (p Perm) Order() int {
	ord := 1
	for _, c := range p.Cycles() {
		ord = Lcm(ord, len(c))
	}
	return ord
}

// Sign returns the sign (parity) of p: +1 if p is an even permutation and -1
// if it is odd. The sign equals (-1)^(n - c) where n is the degree and c is
// the number of disjoint cycles including fixed points.
func (p Perm) Sign() int {
	n := len(p)
	c := len(p.Cycles())
	if (n-c)%2 == 0 {
		return 1
	}
	return -1
}

// IsEven reports whether p is an even permutation (sign +1).
func (p Perm) IsEven() bool {
	return p.Sign() == 1
}

// Pow returns p raised to the integer power k. Negative k uses the inverse of
// p, and Pow(0) is the identity of the same degree.
func (p Perm) Pow(k int) Perm {
	base := p
	if k < 0 {
		base = p.Inverse()
		k = -k
	}
	result := Identity(len(p))
	for i := 0; i < k; i++ {
		result = Compose(base, result)
	}
	return result
}

// NumInversions returns the number of inversions of p: pairs i < j with
// p[i] > p[j]. The parity of this count matches [Perm.Sign].
func (p Perm) NumInversions() int {
	count := 0
	for i := 0; i < len(p); i++ {
		for j := i + 1; j < len(p); j++ {
			if p[i] > p[j] {
				count++
			}
		}
	}
	return count
}

// String renders p in disjoint-cycle notation, e.g. "(0 2 1)(3 4)". The
// identity permutation renders as "()". Fixed points are omitted.
func (p Perm) String() string {
	var b strings.Builder
	nonTrivial := false
	for _, c := range p.Cycles() {
		if len(c) == 1 {
			continue
		}
		nonTrivial = true
		b.WriteByte('(')
		for i, v := range c {
			if i > 0 {
				b.WriteByte(' ')
			}
			fmt.Fprintf(&b, "%d", v)
		}
		b.WriteByte(')')
	}
	if !nonTrivial {
		return "()"
	}
	return b.String()
}

// PermFromCycles builds a permutation of degree n from disjoint cycles given
// as slices of points. Points not mentioned are fixed. It panics if n < 0, if
// any point is out of range, or if a point appears more than once across the
// cycles.
func PermFromCycles(n int, cycles [][]int) Perm {
	if n < 0 {
		panic("groups: PermFromCycles requires n >= 0")
	}
	p := Identity(n)
	used := make([]bool, n)
	for _, cyc := range cycles {
		for _, v := range cyc {
			if v < 0 || v >= n {
				panic("groups: PermFromCycles point out of range")
			}
			if used[v] {
				panic("groups: PermFromCycles cycles are not disjoint")
			}
			used[v] = true
		}
		m := len(cyc)
		for i := 0; i < m; i++ {
			p[cyc[i]] = cyc[(i+1)%m]
		}
	}
	return p
}

// Transposition returns the permutation of degree n that swaps i and j and
// fixes every other point. It panics if i or j is out of range.
func Transposition(n, i, j int) Perm {
	if i < 0 || i >= n || j < 0 || j >= n {
		panic("groups: Transposition index out of range")
	}
	p := Identity(n)
	p[i], p[j] = j, i
	return p
}

// Factorial returns n! as an int. n must be non-negative and small enough that
// n! fits in an int; Factorial(0) == 1.
func Factorial(n int) int {
	if n < 0 {
		panic("groups: Factorial requires n >= 0")
	}
	result := 1
	for i := 2; i <= n; i++ {
		result *= i
	}
	return result
}

// SymmetricGroupOrder returns |S_n| = n!, the order of the symmetric group on
// n points. n must be non-negative.
func SymmetricGroupOrder(n int) int {
	return Factorial(n)
}

// AlternatingGroupOrder returns |A_n|, the order of the alternating group on n
// points: n!/2 for n >= 2 and 1 for n in {0, 1}. n must be non-negative.
func AlternatingGroupOrder(n int) int {
	if n < 0 {
		panic("groups: AlternatingGroupOrder requires n >= 0")
	}
	if n < 2 {
		return 1
	}
	return Factorial(n) / 2
}

// SymmetricGroup returns all n! permutations of degree n in lexicographic
// order of their one-line notation. n must be non-negative and small.
func SymmetricGroup(n int) []Perm {
	if n < 0 {
		panic("groups: SymmetricGroup requires n >= 0")
	}
	var result []Perm
	cur := Identity(n)
	result = append(result, append(Perm(nil), cur...))
	for {
		next, ok := groupsNextPerm(cur)
		if !ok {
			break
		}
		cur = next
		result = append(result, append(Perm(nil), cur...))
	}
	return result
}

// AlternatingGroup returns all even permutations of degree n (the alternating
// group A_n) in lexicographic order. n must be non-negative and small.
func AlternatingGroup(n int) []Perm {
	var result []Perm
	for _, p := range SymmetricGroup(n) {
		if p.Sign() == 1 {
			result = append(result, p)
		}
	}
	return result
}

// groupsNextPerm advances p to the next permutation in lexicographic order,
// returning a fresh slice and true, or (nil, false) if p is already the last.
func groupsNextPerm(p Perm) (Perm, bool) {
	n := len(p)
	next := append(Perm(nil), p...)
	i := n - 2
	for i >= 0 && next[i] >= next[i+1] {
		i--
	}
	if i < 0 {
		return nil, false
	}
	j := n - 1
	for next[j] <= next[i] {
		j--
	}
	next[i], next[j] = next[j], next[i]
	for l, r := i+1, n-1; l < r; l, r = l+1, r-1 {
		next[l], next[r] = next[r], next[l]
	}
	return next, true
}
