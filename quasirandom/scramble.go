package quasirandom

// IdentityPermutation returns the identity permutation (0,1,...,n-1). It returns
// an empty slice for non-positive n.
func IdentityPermutation(n int) []int {
	if n < 0 {
		n = 0
	}
	p := make([]int, n)
	for i := range p {
		p[i] = i
	}
	return p
}

// ReversePermutation returns the permutation (n-1,...,1,0) that reverses digit
// order, a simple scrambling of the digits {0,...,n-1}.
func ReversePermutation(n int) []int {
	if n < 0 {
		n = 0
	}
	p := make([]int, n)
	for i := range p {
		p[i] = n - 1 - i
	}
	return p
}

// AffinePermutation returns the permutation i -> (a*i + b) mod n over the digit
// set {0,...,n-1}. When a is coprime to n this is a bijection and hence a valid
// scrambling permutation; otherwise the returned slice is not a permutation.
func AffinePermutation(n, a, b int) []int {
	if n < 0 {
		n = 0
	}
	p := make([]int, n)
	for i := range p {
		v := (a*i + b) % n
		if v < 0 {
			v += n
		}
		p[i] = v
	}
	return p
}

// FaurePermutation returns Faure's deterministic scrambling permutation of the
// digit set {0,...,base-1}. The permutations are defined recursively — even
// bases interleave two copies of the half-base permutation, odd bases splice
// the fixed middle digit into the shifted previous permutation — and are the
// standard choice for improving the uniformity of scrambled Halton and
// van der Corput sequences. For base < 2 the identity of that size is returned.
func FaurePermutation(base int) []int {
	if base < 2 {
		return IdentityPermutation(base)
	}
	return faurePerm(base)
}

func faurePerm(b int) []int {
	if b == 1 {
		return []int{0}
	}
	if b == 2 {
		return []int{0, 1}
	}
	if b%2 == 0 {
		c := b / 2
		s := faurePerm(c)
		out := make([]int, b)
		for j := 0; j < c; j++ {
			out[j] = 2 * s[j]
			out[c+j] = 2*s[j] + 1
		}
		return out
	}
	c := (b - 1) / 2
	t := faurePerm(b - 1)
	for i := range t {
		if t[i] >= c {
			t[i]++
		}
	}
	out := make([]int, 0, b)
	out = append(out, t[:c]...)
	out = append(out, c)
	out = append(out, t[c:]...)
	return out
}

// IsPermutation reports whether p is a permutation of {0,...,len(p)-1}: each
// value appears exactly once and lies in range.
func IsPermutation(p []int) bool {
	seen := make([]bool, len(p))
	for _, v := range p {
		if v < 0 || v >= len(p) || seen[v] {
			return false
		}
		seen[v] = true
	}
	return true
}

// InvertPermutation returns the inverse of the permutation p, so that
// InvertPermutation(p)[p[i]] == i. It returns nil when p is not a permutation.
func InvertPermutation(p []int) []int {
	if !IsPermutation(p) {
		return nil
	}
	inv := make([]int, len(p))
	for i, v := range p {
		inv[v] = i
	}
	return inv
}

// ComposePermutations returns the permutation r with r[i] = p[q[i]], the effect
// of applying q and then p. It returns nil when the lengths differ or either
// argument is not a permutation.
func ComposePermutations(p, q []int) []int {
	if len(p) != len(q) || !IsPermutation(p) || !IsPermutation(q) {
		return nil
	}
	r := make([]int, len(p))
	for i := range q {
		r[i] = p[q[i]]
	}
	return r
}

// ApplyPermutation returns a new slice with entry i taken from xs[p[i]],
// reordering xs by the permutation p. It returns nil when the lengths differ.
func ApplyPermutation(xs []float64, p []int) []float64 {
	if len(xs) != len(p) {
		return nil
	}
	out := make([]float64, len(xs))
	for i, pi := range p {
		out[i] = xs[pi]
	}
	return out
}

// DigitScramble returns the value obtained by scrambling the base-b digits of
// the radical inverse of n with the single permutation perm applied to every
// digit position, equivalent to ScrambledRadicalInverse.
func DigitScramble(base int, n uint64, perm []int) (float64, error) {
	return ScrambledRadicalInverse(base, n, perm)
}
