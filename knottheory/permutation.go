package knottheory

import (
	"fmt"
	"strings"
)

// Permutation is a bijection of the set {0, 1, ..., n-1} represented by the
// slice of images: p[i] is the image of i. It is used to describe the
// underlying permutation of a braid and the number of components of its
// closure.
type Permutation []int

// NewPermutation returns the permutation with the given image slice after
// checking that it is a genuine bijection of {0, ..., n-1}. It returns an error
// if the slice is not a permutation.
func NewPermutation(images []int) (Permutation, error) {
	n := len(images)
	seen := make([]bool, n)
	for _, v := range images {
		if v < 0 || v >= n || seen[v] {
			return nil, fmt.Errorf("knottheory: %v is not a permutation of {0,...,%d}", images, n-1)
		}
		seen[v] = true
	}
	p := make(Permutation, n)
	copy(p, images)
	return p, nil
}

// IdentityPermutation returns the identity permutation on n points.
func IdentityPermutation(n int) Permutation {
	p := make(Permutation, n)
	for i := range p {
		p[i] = i
	}
	return p
}

// TranspositionPermutation returns the permutation on n points that swaps i and
// j and fixes everything else. It returns an error if i or j is out of range.
func TranspositionPermutation(n, i, j int) (Permutation, error) {
	if i < 0 || j < 0 || i >= n || j >= n {
		return nil, fmt.Errorf("knottheory: transposition indices %d,%d out of range for n=%d", i, j, n)
	}
	p := IdentityPermutation(n)
	p[i], p[j] = j, i
	return p, nil
}

// Size returns the number of points the permutation acts on.
func (p Permutation) Size() int { return len(p) }

// Apply returns the image of i under the permutation.
func (p Permutation) Apply(i int) int { return p[i] }

// Compose returns the permutation p∘q defined by (p∘q)(i) = p(q(i)). Both
// permutations must act on the same number of points.
func (p Permutation) Compose(q Permutation) Permutation {
	n := len(p)
	r := make(Permutation, n)
	for i := 0; i < n; i++ {
		r[i] = p[q[i]]
	}
	return r
}

// Inverse returns the inverse permutation.
func (p Permutation) Inverse() Permutation {
	r := make(Permutation, len(p))
	for i, v := range p {
		r[v] = i
	}
	return r
}

// Equal reports whether p and q are the same permutation.
func (p Permutation) Equal(q Permutation) bool {
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

// IsIdentity reports whether p fixes every point.
func (p Permutation) IsIdentity() bool {
	for i, v := range p {
		if v != i {
			return false
		}
	}
	return true
}

// Cycles returns the cycle decomposition of p as a slice of cycles, each cycle
// being the list of points visited starting from its smallest element. Fixed
// points appear as singleton cycles.
func (p Permutation) Cycles() [][]int {
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

// NumCycles returns the number of cycles in the decomposition of p, counting
// fixed points. For the permutation of a braid this equals the number of
// components of the braid closure.
func (p Permutation) NumCycles() int { return len(p.Cycles()) }

// Order returns the multiplicative order of p, that is the least positive k with
// p^k the identity. It equals the least common multiple of the cycle lengths.
func (p Permutation) Order() int {
	ord := 1
	for _, c := range p.Cycles() {
		ord = lcmInt(ord, len(c))
	}
	return ord
}

// Sign returns +1 if p is an even permutation and -1 if it is odd. The sign is
// (-1)^(n - number of cycles).
func (p Permutation) Sign() int {
	n := len(p)
	c := p.NumCycles()
	if (n-c)%2 == 0 {
		return 1
	}
	return -1
}

// Power returns p raised to the integer power k, which may be negative.
func (p Permutation) Power(k int) Permutation {
	base := p
	if k < 0 {
		base = p.Inverse()
		k = -k
	}
	r := IdentityPermutation(len(p))
	for i := 0; i < k; i++ {
		r = r.Compose(base)
	}
	return r
}

// String renders p in one-line notation, for example "[1 2 0]".
func (p Permutation) String() string {
	parts := make([]string, len(p))
	for i, v := range p {
		parts[i] = fmt.Sprintf("%d", v)
	}
	return "[" + strings.Join(parts, " ") + "]"
}

// lcmInt returns the least common multiple of a and b, or 0 if either is 0.
func lcmInt(a, b int) int {
	if a == 0 || b == 0 {
		return 0
	}
	g := gcdInt(a, b)
	if g == 0 {
		return 0
	}
	return a / g * b
}
