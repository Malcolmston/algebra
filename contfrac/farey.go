package contfrac

// FareySequence returns the Farey sequence F_n: all reduced fractions in the
// closed interval [0, 1] with denominator at most n, listed in ascending order
// from 0/1 to 1/1. It panics for n < 1.
func FareySequence(n int64) []Frac {
	if n < 1 {
		panic("contfrac: FareySequence requires n >= 1")
	}
	out := []Frac{{0, 1}}
	a, b, c, d := int64(0), int64(1), int64(1), n
	for c <= n {
		k := (n + b) / d
		a, b, c, d = c, d, k*c-a, k*d-b
		out = append(out, Frac{a, b})
	}
	return out
}

// FareyLength returns the number of terms in the Farey sequence F_n, which is
// 1 + sum_{k=1}^{n} phi(k). It panics for n < 1.
func FareyLength(n int64) int64 {
	if n < 1 {
		panic("contfrac: FareyLength requires n >= 1")
	}
	total := int64(1)
	for k := int64(1); k <= n; k++ {
		total += EulerPhi(k)
	}
	return total
}

// FareySuccessor returns the term immediately after a/b in the Farey sequence
// F_n. The fraction a/b must be reduced with 0 <= a/b < 1 and b <= n. It panics
// if a/b == 1/1, which has no successor in F_n.
func FareySuccessor(a, b, n int64) Frac {
	a, b = ReduceFraction(a, b)
	if a == b {
		panic("contfrac: 1/1 has no Farey successor")
	}
	// Successor c/d satisfies b*c - a*d = 1 with the largest d <= n.
	// d ≡ -a^{-1} (mod b); pick the largest such d in [1, n].
	d := largestCongruent(negMod(modInverse(a, b), b), b, n)
	c := (1 + a*d) / b
	return NewFrac(c, d)
}

// FareyPredecessor returns the term immediately before a/b in the Farey
// sequence F_n. The fraction a/b must be reduced with 0 < a/b <= 1 and b <= n.
// It panics if a/b == 0/1, which has no predecessor.
func FareyPredecessor(a, b, n int64) Frac {
	a, b = ReduceFraction(a, b)
	if a == 0 {
		panic("contfrac: 0/1 has no Farey predecessor")
	}
	// Predecessor c/d satisfies a*d - b*c = 1 with the largest d <= n.
	// d ≡ a^{-1} (mod b).
	d := largestCongruent(modInverse(a, b), b, n)
	c := (a*d - 1) / b
	return NewFrac(c, d)
}

// FareyNeighbors returns the predecessor and successor of a/b in F_n. If a/b is
// an endpoint (0/1 or 1/1) the missing neighbour is returned as the endpoint
// itself.
func FareyNeighbors(a, b, n int64) (prev, next Frac) {
	a, b = ReduceFraction(a, b)
	if a == 0 {
		prev = Frac{0, 1}
	} else {
		prev = FareyPredecessor(a, b, n)
	}
	if a == b {
		next = Frac{1, 1}
	} else {
		next = FareySuccessor(a, b, n)
	}
	return prev, next
}

// FareyNext returns the term after c/d given the two consecutive Farey terms
// a/b and c/d of F_n, using the standard mediant recurrence
// k = floor((n+b)/d), next = (k*c-a)/(k*d-b).
func FareyNext(a, b, c, d, n int64) Frac {
	k := (n + b) / d
	return Frac{k*c - a, k*d - b}
}

// AreFareyNeighbors reports whether a/b and c/d are adjacent in some Farey
// sequence, i.e. |b*c - a*d| == 1.
func AreFareyNeighbors(a, b, c, d int64) bool {
	return absInt(b*c-a*d) == 1
}

// modInverse returns the modular inverse of a modulo m (m > 0), in the range
// [0, m). For m == 1 it returns 0. It assumes gcd(a, m) == 1.
func modInverse(a, m int64) int64 {
	if m == 1 {
		return 0
	}
	a = mod(a, m)
	// Extended Euclidean algorithm.
	g, x := m, int64(0)
	newg, newx := a, int64(1)
	for newg != 0 {
		q := g / newg
		g, newg = newg, g-q*newg
		x, newx = newx, x-q*newx
	}
	return mod(x, m)
}

// negMod returns (-v) reduced modulo m into [0, m).
func negMod(v, m int64) int64 {
	return mod(-v, m)
}

// largestCongruent returns the largest value d in [1, n] with d ≡ r (mod m),
// where r is already reduced into [0, m). When r == 0 the representative m is
// used so that d stays positive.
func largestCongruent(r, m, n int64) int64 {
	if r == 0 {
		r = m
	}
	if r > n {
		return r // smallest positive representative (may exceed n for tiny n)
	}
	return r + ((n-r)/m)*m
}
