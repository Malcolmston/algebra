package combin

import (
	"math/big"
	"sort"
)

// ----------------------------------------------------------------------------
// unexported helpers (all prefixed with "combin" to avoid collisions with the
// sibling file in this package)
// ----------------------------------------------------------------------------

// combinMinInt returns the smaller of a and b.
func combinMinInt(a, b int) int {
	if a < b {
		return a
	}
	return b
}

// combinCopyInts returns a fresh copy of s.
func combinCopyInts(s []int) []int {
	out := make([]int, len(s))
	copy(out, s)
	return out
}

// combinFactorial returns n! as a big.Int (0! == 1). Negative n yields 0.
func combinFactorial(n int) *big.Int {
	if n < 0 {
		return big.NewInt(0)
	}
	r := big.NewInt(1)
	t := new(big.Int)
	for i := 2; i <= n; i++ {
		r.Mul(r, t.SetInt64(int64(i)))
	}
	return r
}

// combinBinom returns the binomial coefficient C(n,k) as a big.Int.
func combinBinom(n, k int) *big.Int {
	if k < 0 || n < 0 || k > n {
		return big.NewInt(0)
	}
	if k > n-k {
		k = n - k
	}
	r := big.NewInt(1)
	t := new(big.Int)
	for i := 0; i < k; i++ {
		r.Mul(r, t.SetInt64(int64(n-i)))
		r.Div(r, t.SetInt64(int64(i+1)))
	}
	return r
}

// combinPowInt returns base**exp as a big.Int for exp >= 0.
func combinPowInt(base, exp int) *big.Int {
	return new(big.Int).Exp(big.NewInt(int64(base)), big.NewInt(int64(exp)), nil)
}

// combinDivisors returns the positive divisors of n in ascending order.
func combinDivisors(n int) []int {
	var ds []int
	for i := 1; i*i <= n; i++ {
		if n%i == 0 {
			ds = append(ds, i)
			if i != n/i {
				ds = append(ds, n/i)
			}
		}
	}
	sort.Ints(ds)
	return ds
}

// combinTotient returns Euler's totient phi(n) for n >= 1.
func combinTotient(n int) int {
	result := n
	nn := n
	for p := 2; p*p <= nn; p++ {
		if nn%p == 0 {
			for nn%p == 0 {
				nn /= p
			}
			result -= result / p
		}
	}
	if nn > 1 {
		result -= result / nn
	}
	return result
}

// combinMobius returns the Moebius function mu(n) for n >= 1.
func combinMobius(n int) int {
	if n == 1 {
		return 1
	}
	result := 1
	nn := n
	for p := 2; p*p <= nn; p++ {
		if nn%p == 0 {
			nn /= p
			if nn%p == 0 {
				return 0
			}
			result = -result
		}
	}
	if nn > 1 {
		result = -result
	}
	return result
}

// combinDerangementBig returns the number of derangements !n as a big.Int.
func combinDerangementBig(n int) *big.Int {
	if n < 0 {
		return big.NewInt(0)
	}
	if n == 0 {
		return big.NewInt(1)
	}
	if n == 1 {
		return big.NewInt(0)
	}
	prev2 := big.NewInt(1) // !0
	prev1 := big.NewInt(0) // !1
	cur := new(big.Int)
	for i := 2; i <= n; i++ {
		// !i = (i-1) * (!(i-1) + !(i-2))
		cur = new(big.Int).Add(prev1, prev2)
		cur.Mul(cur, big.NewInt(int64(i-1)))
		prev2, prev1 = prev1, cur
	}
	return new(big.Int).Set(cur)
}

// combinZigzag returns the n-th Euler up/down (zigzag) number A_n using the
// Entringer number triangle. A_0 = A_1 = 1.
func combinZigzag(n int) *big.Int {
	if n < 0 {
		return big.NewInt(0)
	}
	if n == 0 {
		return big.NewInt(1)
	}
	prev := []*big.Int{big.NewInt(1)} // Entringer row 0: E(0,0) = 1
	for r := 1; r <= n; r++ {
		cur := make([]*big.Int, r+1)
		cur[0] = big.NewInt(0)
		for k := 1; k <= r; k++ {
			cur[k] = new(big.Int).Add(cur[k-1], prev[r-k])
		}
		prev = cur
	}
	return prev[n]
}

// ----------------------------------------------------------------------------
// Integer partitions
// ----------------------------------------------------------------------------

// PartitionNumber returns p(n), the number of integer partitions of n, computed
// by an O(n^2) dynamic program. It returns 1 for n == 0 and 0 for n < 0.
func PartitionNumber(n int) int {
	if n < 0 {
		return 0
	}
	dp := make([]int, n+1)
	dp[0] = 1
	for part := 1; part <= n; part++ {
		for j := part; j <= n; j++ {
			dp[j] += dp[j-part]
		}
	}
	return dp[n]
}

// PartitionNumberTable returns the slice [p(0), p(1), ..., p(n)].
func PartitionNumberTable(n int) []int {
	if n < 0 {
		return nil
	}
	dp := make([]int, n+1)
	dp[0] = 1
	for part := 1; part <= n; part++ {
		for j := part; j <= n; j++ {
			dp[j] += dp[j-part]
		}
	}
	return dp
}

// PartitionNumberInto returns the number of partitions of n into exactly k
// positive parts, using the recurrence P(n,k) = P(n-1,k-1) + P(n-k,k).
func PartitionNumberInto(n, k int) int {
	if k < 0 || n < 0 {
		return 0
	}
	if k == 0 {
		if n == 0 {
			return 1
		}
		return 0
	}
	if n < k {
		return 0
	}
	dp := make([][]int, n+1)
	for i := range dp {
		dp[i] = make([]int, k+1)
	}
	dp[0][0] = 1
	for i := 1; i <= n; i++ {
		for j := 1; j <= k && j <= i; j++ {
			v := dp[i-1][j-1]
			if i-j >= 0 {
				v += dp[i-j][j]
			}
			dp[i][j] = v
		}
	}
	return dp[n][k]
}

// PartitionNumberAtMost returns the number of partitions of n into at most k
// parts, equivalently the number of partitions of n whose parts are all <= k.
func PartitionNumberAtMost(n, k int) int {
	if n < 0 || k < 0 {
		return 0
	}
	dp := make([]int, n+1)
	dp[0] = 1
	for part := 1; part <= k; part++ {
		for j := part; j <= n; j++ {
			dp[j] += dp[j-part]
		}
	}
	return dp[n]
}

// PartitionNumberDistinct returns the number of partitions of n into distinct
// (unequal) positive parts.
func PartitionNumberDistinct(n int) int {
	if n < 0 {
		return 0
	}
	dp := make([]int, n+1)
	dp[0] = 1
	for part := 1; part <= n; part++ {
		for j := n; j >= part; j-- {
			dp[j] += dp[j-part]
		}
	}
	return dp[n]
}

// PartitionNumberOdd returns the number of partitions of n into odd parts.
// By Euler's theorem this equals PartitionNumberDistinct(n).
func PartitionNumberOdd(n int) int {
	if n < 0 {
		return 0
	}
	dp := make([]int, n+1)
	dp[0] = 1
	for part := 1; part <= n; part += 2 {
		for j := part; j <= n; j++ {
			dp[j] += dp[j-part]
		}
	}
	return dp[n]
}

// combinPartitionsMax returns every partition of n (non-increasing) whose
// largest part is at most maxPart.
func combinPartitionsMax(n, maxPart int) [][]int {
	if n == 0 {
		return [][]int{{}}
	}
	var res [][]int
	for first := combinMinInt(n, maxPart); first >= 1; first-- {
		for _, rest := range combinPartitionsMax(n-first, first) {
			p := append([]int{first}, rest...)
			res = append(res, p)
		}
	}
	return res
}

// Partitions returns all integer partitions of n as non-increasing slices, in
// reverse lexicographic order (starting with {n} and ending with all ones).
func Partitions(n int) [][]int {
	if n < 0 {
		return nil
	}
	return combinPartitionsMax(n, n)
}

// combinPartitionsInto returns every partition of n into exactly k parts
// (non-increasing) whose largest part is at most maxPart.
func combinPartitionsInto(n, k, maxPart int) [][]int {
	if k == 0 {
		if n == 0 {
			return [][]int{{}}
		}
		return nil
	}
	var res [][]int
	hi := combinMinInt(n-(k-1), maxPart)
	for first := hi; first >= 1; first-- {
		for _, rest := range combinPartitionsInto(n-first, k-1, first) {
			p := append([]int{first}, rest...)
			res = append(res, p)
		}
	}
	return res
}

// PartitionsInto returns all partitions of n into exactly k positive parts,
// each as a non-increasing slice.
func PartitionsInto(n, k int) [][]int {
	if n < 0 || k < 0 {
		return nil
	}
	return combinPartitionsInto(n, k, n)
}

// combinPartitionsDistinct returns every partition of n into distinct parts
// (strictly decreasing) whose largest part is at most maxPart.
func combinPartitionsDistinct(n, maxPart int) [][]int {
	if n == 0 {
		return [][]int{{}}
	}
	var res [][]int
	for first := combinMinInt(n, maxPart); first >= 1; first-- {
		for _, rest := range combinPartitionsDistinct(n-first, first-1) {
			p := append([]int{first}, rest...)
			res = append(res, p)
		}
	}
	return res
}

// PartitionsDistinct returns all partitions of n into distinct positive parts,
// each as a strictly decreasing slice.
func PartitionsDistinct(n int) [][]int {
	if n < 0 {
		return nil
	}
	return combinPartitionsDistinct(n, n)
}

// NextPartition advances p (a non-increasing partition) in place to the next
// partition of the same integer in reverse lexicographic order. It returns the
// updated slice (which may be reallocated) and true, or the unchanged input and
// false when p is already the final partition (all ones).
func NextPartition(p []int) ([]int, bool) {
	sum := 0
	for _, x := range p {
		sum += x
	}
	if len(p) == sum { // all parts equal to 1 (also handles the empty partition)
		return p, false
	}
	k := len(p) - 1
	for p[k] == 1 {
		k--
	}
	rem := (len(p) - 1 - k) + 1
	p[k]--
	val := p[k]
	p = p[:k+1]
	for rem > val {
		p = append(p, val)
		rem -= val
	}
	if rem > 0 {
		p = append(p, rem)
	}
	return p, true
}

// PartitionConjugate returns the conjugate (transpose) of the partition p, which
// must be given in non-increasing order. The conjugate is also non-increasing.
func PartitionConjugate(p []int) []int {
	if len(p) == 0 {
		return []int{}
	}
	max := p[0]
	conj := make([]int, max)
	for _, part := range p {
		for i := 0; i < part; i++ {
			conj[i]++
		}
	}
	return conj
}

// PartitionRank returns the rank of the partition p (largest part minus the
// number of parts). The empty partition has rank 0.
func PartitionRank(p []int) int {
	if len(p) == 0 {
		return 0
	}
	return p[0] - len(p)
}

// PartitionDurfeeSquare returns the side length of the Durfee square of the
// partition p (the largest d with at least d parts each of size >= d). p must be
// non-increasing.
func PartitionDurfeeSquare(p []int) int {
	d := 0
	for i, part := range p {
		if part >= i+1 {
			d = i + 1
		} else {
			break
		}
	}
	return d
}

// PentagonalNumber returns the k-th (signed) pentagonal number k(3k-1)/2. It
// accepts negative k, yielding the generalized pentagonal numbers.
func PentagonalNumber(k int) int {
	return k * (3*k - 1) / 2
}

// GeneralizedPentagonal returns the i-th generalized pentagonal number in the
// order 0, 1, 2, 5, 7, 12, 15, ... generated by k = 0, 1, -1, 2, -2, 3, ...
func GeneralizedPentagonal(i int) int {
	if i <= 0 {
		return 0
	}
	m := (i + 1) / 2
	if i%2 == 0 {
		m = -m
	}
	return PentagonalNumber(m)
}

// ----------------------------------------------------------------------------
// Compositions
// ----------------------------------------------------------------------------

// CompositionNumber returns the number of compositions (ordered partitions) of
// n into positive parts, which is 2^(n-1) for n >= 1 and 1 for n == 0.
func CompositionNumber(n int) int {
	if n < 0 {
		return 0
	}
	if n == 0 {
		return 1
	}
	return 1 << (n - 1)
}

// CompositionNumberInto returns the number of compositions of n into exactly k
// positive parts, which is C(n-1, k-1).
func CompositionNumberInto(n, k int) int {
	if k == 0 {
		if n == 0 {
			return 1
		}
		return 0
	}
	return int(combinBinom(n-1, k-1).Int64())
}

// WeakCompositionNumber returns the number of weak compositions of n into
// exactly k non-negative parts, which is C(n+k-1, k-1).
func WeakCompositionNumber(n, k int) int {
	if k <= 0 {
		if n == 0 && k == 0 {
			return 1
		}
		return 0
	}
	return int(combinBinom(n+k-1, k-1).Int64())
}

// combinCompositions returns all compositions of n into positive parts.
func combinCompositions(n int) [][]int {
	if n == 0 {
		return [][]int{{}}
	}
	var res [][]int
	for first := 1; first <= n; first++ {
		for _, rest := range combinCompositions(n - first) {
			c := append([]int{first}, rest...)
			res = append(res, c)
		}
	}
	return res
}

// CompositionList returns all compositions of n into positive parts, in
// lexicographic order.
func CompositionList(n int) [][]int {
	if n < 0 {
		return nil
	}
	return combinCompositions(n)
}

// combinCompositionsInto returns all compositions of n into exactly k positive
// parts.
func combinCompositionsInto(n, k int) [][]int {
	if k == 0 {
		if n == 0 {
			return [][]int{{}}
		}
		return nil
	}
	var res [][]int
	for first := 1; first <= n-(k-1); first++ {
		for _, rest := range combinCompositionsInto(n-first, k-1) {
			c := append([]int{first}, rest...)
			res = append(res, c)
		}
	}
	return res
}

// CompositionListInto returns all compositions of n into exactly k positive
// parts.
func CompositionListInto(n, k int) [][]int {
	if n < 0 || k < 0 {
		return nil
	}
	return combinCompositionsInto(n, k)
}

// combinWeakCompositions returns all weak compositions of n into exactly k
// non-negative parts.
func combinWeakCompositions(n, k int) [][]int {
	if k == 0 {
		if n == 0 {
			return [][]int{{}}
		}
		return nil
	}
	if k == 1 {
		return [][]int{{n}}
	}
	var res [][]int
	for first := 0; first <= n; first++ {
		for _, rest := range combinWeakCompositions(n-first, k-1) {
			c := append([]int{first}, rest...)
			res = append(res, c)
		}
	}
	return res
}

// WeakCompositionList returns all weak compositions of n into exactly k
// non-negative parts, in lexicographic order.
func WeakCompositionList(n, k int) [][]int {
	if n < 0 || k < 0 {
		return nil
	}
	return combinWeakCompositions(n, k)
}

// ----------------------------------------------------------------------------
// Derangements and involutions
// ----------------------------------------------------------------------------

// DerangementNumber returns !n, the number of permutations of n elements with no
// fixed point, using the recurrence !n = (n-1)(!(n-1) + !(n-2)).
func DerangementNumber(n int) int {
	if n < 0 {
		return 0
	}
	if n == 0 {
		return 1
	}
	if n == 1 {
		return 0
	}
	p2, p1 := 1, 0
	cur := 0
	for i := 2; i <= n; i++ {
		cur = (i - 1) * (p1 + p2)
		p2, p1 = p1, cur
	}
	return cur
}

// DerangementProbability returns !n / n!, the probability that a uniformly
// random permutation of n elements is a derangement. It converges to 1/e.
func DerangementProbability(n int) float64 {
	if n < 0 {
		return 0
	}
	r := new(big.Rat).SetFrac(combinDerangementBig(n), combinFactorial(n))
	f, _ := r.Float64()
	return f
}

// DerangementList returns every derangement of {0, 1, ..., n-1} as a slice, in
// lexicographic order.
func DerangementList(n int) [][]int {
	perms := AllPermutations(combinRange(n))
	var res [][]int
	for _, p := range perms {
		fixed := false
		for i, v := range p {
			if v == i {
				fixed = true
				break
			}
		}
		if !fixed {
			res = append(res, p)
		}
	}
	combinSortIntSlices(res)
	return res
}

// InvolutionNumber returns the number of involutions (self-inverse permutations)
// of n elements, using I(n) = I(n-1) + (n-1) I(n-2).
func InvolutionNumber(n int) int {
	if n < 0 {
		return 0
	}
	if n == 0 || n == 1 {
		return 1
	}
	p2, p1 := 1, 1
	cur := 0
	for i := 2; i <= n; i++ {
		cur = p1 + (i-1)*p2
		p2, p1 = p1, cur
	}
	return cur
}

// InvolutionList returns every involution of {0, 1, ..., n-1} as a slice, in
// lexicographic order.
func InvolutionList(n int) [][]int {
	perms := AllPermutations(combinRange(n))
	var res [][]int
	for _, p := range perms {
		isInv := true
		for i, v := range p {
			if p[v] != i {
				isInv = false
				break
			}
		}
		if isInv {
			res = append(res, p)
		}
	}
	combinSortIntSlices(res)
	return res
}

// combinRange returns the slice {0, 1, ..., n-1}.
func combinRange(n int) []int {
	s := make([]int, n)
	for i := range s {
		s[i] = i
	}
	return s
}

// combinSortIntSlices sorts a slice of int slices lexicographically in place.
func combinSortIntSlices(s [][]int) {
	sort.Slice(s, func(a, b int) bool {
		x, y := s[a], s[b]
		for i := 0; i < len(x) && i < len(y); i++ {
			if x[i] != y[i] {
				return x[i] < y[i]
			}
		}
		return len(x) < len(y)
	})
}

// ----------------------------------------------------------------------------
// Permutations
// ----------------------------------------------------------------------------

// NextPermutation rearranges a into the lexicographically next greater
// permutation of its elements and returns true. If a is already the highest
// permutation it is reversed to the lowest (sorted ascending) permutation and
// false is returned. Equal elements are handled correctly.
func NextPermutation(a []int) bool {
	n := len(a)
	if n < 2 {
		return false
	}
	i := n - 2
	for i >= 0 && a[i] >= a[i+1] {
		i--
	}
	if i < 0 {
		for l, r := 0, n-1; l < r; l, r = l+1, r-1 {
			a[l], a[r] = a[r], a[l]
		}
		return false
	}
	j := n - 1
	for a[j] <= a[i] {
		j--
	}
	a[i], a[j] = a[j], a[i]
	for l, r := i+1, n-1; l < r; l, r = l+1, r-1 {
		a[l], a[r] = a[r], a[l]
	}
	return true
}

// PrevPermutation rearranges a into the lexicographically previous permutation
// of its elements and returns true. If a is already the lowest permutation it is
// reversed to the highest permutation and false is returned.
func PrevPermutation(a []int) bool {
	n := len(a)
	if n < 2 {
		return false
	}
	i := n - 2
	for i >= 0 && a[i] <= a[i+1] {
		i--
	}
	if i < 0 {
		for l, r := 0, n-1; l < r; l, r = l+1, r-1 {
			a[l], a[r] = a[r], a[l]
		}
		return false
	}
	j := n - 1
	for a[j] >= a[i] {
		j--
	}
	a[i], a[j] = a[j], a[i]
	for l, r := i+1, n-1; l < r; l, r = l+1, r-1 {
		a[l], a[r] = a[r], a[l]
	}
	return true
}

// AllPermutations returns all len(items)! permutations of items (treating equal
// values as distinct by position), each as a fresh slice.
func AllPermutations(items []int) [][]int {
	n := len(items)
	if n == 0 {
		return [][]int{{}}
	}
	cur := combinCopyInts(items)
	var res [][]int
	var rec func(k int)
	rec = func(k int) {
		if k == n {
			res = append(res, combinCopyInts(cur))
			return
		}
		for i := k; i < n; i++ {
			cur[k], cur[i] = cur[i], cur[k]
			rec(k + 1)
			cur[k], cur[i] = cur[i], cur[k]
		}
	}
	rec(0)
	return res
}

// MultisetPermutations returns all distinct permutations of items (a multiset),
// in lexicographic order, each as a fresh slice.
func MultisetPermutations(items []int) [][]int {
	a := combinCopyInts(items)
	sort.Ints(a)
	res := [][]int{combinCopyInts(a)}
	for NextPermutation(a) {
		res = append(res, combinCopyInts(a))
	}
	return res
}

// MultisetPermutationList expands the count vector counts into a multiset (value
// i repeated counts[i] times) and returns all its distinct permutations.
func MultisetPermutationList(counts []int) [][]int {
	var items []int
	for v, c := range counts {
		for j := 0; j < c; j++ {
			items = append(items, v)
		}
	}
	return MultisetPermutations(items)
}

// MultisetPermutationNumber returns the number of distinct permutations of a
// multiset whose element multiplicities are counts, i.e. the multinomial
// coefficient (sum counts)! / prod(counts[i]!).
func MultisetPermutationNumber(counts []int) *big.Int {
	total := 0
	for _, c := range counts {
		if c < 0 {
			return big.NewInt(0)
		}
		total += c
	}
	res := combinFactorial(total)
	for _, c := range counts {
		res.Div(res, combinFactorial(c))
	}
	return res
}

// PermutationParity returns the sign of the permutation a of {0,...,n-1}:
// +1 if it has an even number of inversions, -1 if odd.
func PermutationParity(a []int) int {
	inv := 0
	for i := 0; i < len(a); i++ {
		for j := i + 1; j < len(a); j++ {
			if a[i] > a[j] {
				inv++
			}
		}
	}
	if inv%2 == 0 {
		return 1
	}
	return -1
}

// PermutationRank returns the lexicographic rank (0-based) of the permutation a,
// which must be a permutation of {0, 1, ..., len(a)-1}.
func PermutationRank(a []int) *big.Int {
	n := len(a)
	rank := big.NewInt(0)
	for i := 0; i < n; i++ {
		smaller := 0
		for j := i + 1; j < n; j++ {
			if a[j] < a[i] {
				smaller++
			}
		}
		term := new(big.Int).Mul(big.NewInt(int64(smaller)), combinFactorial(n-1-i))
		rank.Add(rank, term)
	}
	return rank
}

// PermutationUnrank returns the permutation of {0, 1, ..., n-1} whose
// lexicographic rank (0-based) is rank.
func PermutationUnrank(n int, rank *big.Int) []int {
	digits := FactorialNumberSystem(rank, n)
	avail := combinRange(n)
	res := make([]int, n)
	for i := 0; i < n; i++ {
		idx := digits[i]
		res[i] = avail[idx]
		avail = append(avail[:idx], avail[idx+1:]...)
	}
	return res
}

// FactorialNumberSystem returns the factoradic (Lehmer) digits of rank as a
// length-n slice; digit i lies in [0, n-1-i] and the final digit is always 0.
func FactorialNumberSystem(rank *big.Int, n int) []int {
	digits := make([]int, n)
	r := new(big.Int).Set(rank)
	q := new(big.Int)
	for i := 0; i < n; i++ {
		f := combinFactorial(n - 1 - i)
		q.DivMod(r, f, r)
		digits[i] = int(q.Int64())
	}
	return digits
}

// ----------------------------------------------------------------------------
// Combinations and subsets
// ----------------------------------------------------------------------------

// NextCombination advances c, a strictly increasing k-combination of indices
// drawn from {0, ..., n-1}, to the next combination in lexicographic order. It
// returns true on success, or false (leaving c unchanged) when c is the last
// combination.
func NextCombination(c []int, n int) bool {
	k := len(c)
	if k == 0 {
		return false
	}
	i := k - 1
	for i >= 0 && c[i] == n-k+i {
		i--
	}
	if i < 0 {
		return false
	}
	c[i]++
	for j := i + 1; j < k; j++ {
		c[j] = c[j-1] + 1
	}
	return true
}

// CombinationList returns all k-element index subsets of {0, 1, ..., n-1}, in
// lexicographic order, each as a strictly increasing slice.
func CombinationList(n, k int) [][]int {
	if k < 0 || k > n {
		return nil
	}
	c := make([]int, k)
	for i := range c {
		c[i] = i
	}
	res := [][]int{combinCopyInts(c)}
	for NextCombination(c, n) {
		res = append(res, combinCopyInts(c))
	}
	return res
}

// CombinationListOf returns all k-element combinations of items, in the order
// induced by index position, each as a fresh slice.
func CombinationListOf(items []int, k int) [][]int {
	n := len(items)
	idx := CombinationList(n, k)
	res := make([][]int, 0, len(idx))
	for _, c := range idx {
		s := make([]int, k)
		for i, ci := range c {
			s[i] = items[ci]
		}
		res = append(res, s)
	}
	return res
}

// CombinationRank returns the colexicographic rank of the strictly increasing
// combination c in the combinatorial number system: sum of C(c[i], i+1).
func CombinationRank(c []int) *big.Int {
	rank := big.NewInt(0)
	for i, v := range c {
		rank.Add(rank, combinBinom(v, i+1))
	}
	return rank
}

// CombinationUnrank returns the k-combination (strictly increasing) whose
// colexicographic rank is rank, inverting CombinationRank.
func CombinationUnrank(k int, rank *big.Int) []int {
	c := make([]int, k)
	r := new(big.Int).Set(rank)
	for i := k; i >= 1; i-- {
		// find the largest v with C(v, i) <= r
		v := i - 1
		for combinBinom(v+1, i).Cmp(r) <= 0 {
			v++
		}
		c[i-1] = v
		r.Sub(r, combinBinom(v, i))
	}
	return c
}

// MultisetCoefficient returns the number of multisets of size k drawn from n
// distinct kinds, i.e. C(n+k-1, k).
func MultisetCoefficient(n, k int) *big.Int {
	return combinBinom(n+k-1, k)
}

// PowerSet returns every subset of items as a fresh slice, ordered by the binary
// value of the inclusion mask (the empty set first).
func PowerSet(items []int) [][]int {
	n := len(items)
	total := 1 << n
	res := make([][]int, 0, total)
	for mask := 0; mask < total; mask++ {
		var s []int
		for i := 0; i < n; i++ {
			if mask&(1<<i) != 0 {
				s = append(s, items[i])
			}
		}
		res = append(res, s)
	}
	return res
}

// PowerSetSize returns 2^n, the number of subsets of an n-element set.
func PowerSetSize(n int) *big.Int {
	return new(big.Int).Lsh(big.NewInt(1), uint(n))
}

// SubsetList returns every subset of {0, 1, ..., n-1} as an index slice, ordered
// by the binary value of the inclusion mask.
func SubsetList(n int) [][]int {
	return PowerSet(combinRange(n))
}

// ----------------------------------------------------------------------------
// Gray codes
// ----------------------------------------------------------------------------

// GrayEncode returns the reflected binary Gray code of n (n XOR (n >> 1)).
func GrayEncode(n uint64) uint64 {
	return n ^ (n >> 1)
}

// GrayDecode returns the integer whose Gray code is g, inverting GrayEncode.
func GrayDecode(g uint64) uint64 {
	g ^= g >> 32
	g ^= g >> 16
	g ^= g >> 8
	g ^= g >> 4
	g ^= g >> 2
	g ^= g >> 1
	return g
}

// GrayCodeRank returns the position (rank) of the Gray code g within the
// standard reflected Gray sequence; it is an alias for GrayDecode.
func GrayCodeRank(g uint64) uint64 {
	return GrayDecode(g)
}

// GrayCodeAt returns the Gray code at position rank in the standard reflected
// Gray sequence; it is an alias for GrayEncode.
func GrayCodeAt(rank uint64) uint64 {
	return GrayEncode(rank)
}

// GrayCodeSequence returns the 2^bits reflected Gray codes in order. Consecutive
// entries (and the last and first) differ in exactly one bit.
func GrayCodeSequence(bits int) []uint64 {
	if bits < 0 {
		return nil
	}
	total := uint64(1) << uint(bits)
	res := make([]uint64, total)
	for i := uint64(0); i < total; i++ {
		res[i] = GrayEncode(i)
	}
	return res
}

// ----------------------------------------------------------------------------
// Triangles and sequences of classical combinatorial numbers
// ----------------------------------------------------------------------------

// NarayanaTriangleRow returns row n of the Narayana triangle, the values
// N(n,1), ..., N(n,n) where N(n,k) = C(n,k) C(n,k-1) / n. The row sums to the
// n-th Catalan number.
func NarayanaTriangleRow(n int) []*big.Int {
	if n <= 0 {
		return []*big.Int{big.NewInt(1)}
	}
	row := make([]*big.Int, n)
	for k := 1; k <= n; k++ {
		v := new(big.Int).Mul(combinBinom(n, k), combinBinom(n, k-1))
		v.Div(v, big.NewInt(int64(n)))
		row[k-1] = v
	}
	return row
}

// NarayanaTriangle returns rows 0 through n of the Narayana triangle.
func NarayanaTriangle(n int) [][]*big.Int {
	rows := make([][]*big.Int, n+1)
	for r := 0; r <= n; r++ {
		rows[r] = NarayanaTriangleRow(r)
	}
	return rows
}

// MotzkinSequence returns the Motzkin numbers M(0), M(1), ..., M(n), using
// M(i) = M(i-1) + sum_{j=0}^{i-2} M(j) M(i-2-j).
func MotzkinSequence(n int) []*big.Int {
	if n < 0 {
		return nil
	}
	m := make([]*big.Int, n+1)
	m[0] = big.NewInt(1)
	if n >= 1 {
		m[1] = big.NewInt(1)
	}
	for i := 2; i <= n; i++ {
		s := new(big.Int).Set(m[i-1])
		for j := 0; j <= i-2; j++ {
			s.Add(s, new(big.Int).Mul(m[j], m[i-2-j]))
		}
		m[i] = s
	}
	return m
}

// CatalanSequence returns the Catalan numbers C(0), C(1), ..., C(n), using
// C(i+1) = C(i) (2)(2i+1) / (i+2).
func CatalanSequence(n int) []*big.Int {
	if n < 0 {
		return nil
	}
	c := make([]*big.Int, n+1)
	c[0] = big.NewInt(1)
	for i := 0; i < n; i++ {
		v := new(big.Int).Mul(c[i], big.NewInt(int64(2*(2*i+1))))
		v.Div(v, big.NewInt(int64(i+2)))
		c[i+1] = v
	}
	return c
}

// EulerianRow returns row n of the Eulerian triangle: the counts of
// permutations of n elements with 0, 1, ..., n-1 ascents. The row sums to n!.
func EulerianRow(n int) []*big.Int {
	if n <= 0 {
		return []*big.Int{big.NewInt(1)}
	}
	prev := []*big.Int{big.NewInt(1)} // row for n = 1
	for r := 2; r <= n; r++ {
		cur := make([]*big.Int, r)
		for j := 0; j < r; j++ {
			t := big.NewInt(0)
			if j < len(prev) {
				t = new(big.Int).Mul(big.NewInt(int64(j+1)), prev[j])
			}
			if j-1 >= 0 && j-1 < len(prev) {
				t.Add(t, new(big.Int).Mul(big.NewInt(int64(r-j)), prev[j-1]))
			}
			cur[j] = t
		}
		prev = cur
	}
	return prev
}

// EulerianTriangle returns rows 0 through n of the Eulerian triangle.
func EulerianTriangle(n int) [][]*big.Int {
	rows := make([][]*big.Int, n+1)
	for r := 0; r <= n; r++ {
		rows[r] = EulerianRow(r)
	}
	return rows
}

// BellTriangleFull returns rows 0 through n of the Bell (Aitken/Peirce)
// triangle. The first entry of row i equals the i-th Bell number.
func BellTriangleFull(n int) [][]*big.Int {
	if n < 0 {
		return nil
	}
	rows := make([][]*big.Int, n+1)
	rows[0] = []*big.Int{big.NewInt(1)}
	for i := 1; i <= n; i++ {
		prev := rows[i-1]
		row := make([]*big.Int, i+1)
		row[0] = new(big.Int).Set(prev[len(prev)-1])
		for j := 1; j <= i; j++ {
			row[j] = new(big.Int).Add(row[j-1], prev[j-1])
		}
		rows[i] = row
	}
	return rows
}

// ----------------------------------------------------------------------------
// Bernoulli, Euler and zigzag numbers
// ----------------------------------------------------------------------------

// BernoulliSequence returns the Bernoulli numbers B(0), B(1), ..., B(n) as exact
// rationals, using the convention B(1) = -1/2.
func BernoulliSequence(n int) []*big.Rat {
	if n < 0 {
		return nil
	}
	b := make([]*big.Rat, n+1)
	for m := 0; m <= n; m++ {
		if m == 0 {
			b[m] = big.NewRat(1, 1)
			continue
		}
		s := new(big.Rat)
		for k := 0; k < m; k++ {
			c := new(big.Rat).SetInt(combinBinom(m+1, k))
			s.Add(s, c.Mul(c, b[k]))
		}
		s.Quo(s, new(big.Rat).SetInt64(int64(m+1)))
		s.Neg(s)
		b[m] = s
	}
	return b
}

// BernoulliNumber returns the n-th Bernoulli number B(n) as an exact rational,
// using the convention B(1) = -1/2.
func BernoulliNumber(n int) *big.Rat {
	seq := BernoulliSequence(n)
	if seq == nil {
		return new(big.Rat)
	}
	return seq[n]
}

// BernoulliNumberFloat returns the n-th Bernoulli number as a float64.
func BernoulliNumberFloat(n int) float64 {
	f, _ := BernoulliNumber(n).Float64()
	return f
}

// ZigzagNumber returns the n-th Euler up/down (zigzag) number A(n): the number
// of alternating permutations of n elements. A(0) = A(1) = 1.
func ZigzagNumber(n int) *big.Int {
	return combinZigzag(n)
}

// SecantNumber returns the n-th secant number |E(2n)| = A(2n), the number of
// alternating permutations of even length 2n.
func SecantNumber(n int) *big.Int {
	return combinZigzag(2 * n)
}

// TangentNumber returns the n-th tangent number A(2n+1), the number of
// alternating permutations of odd length 2n+1.
func TangentNumber(n int) *big.Int {
	return combinZigzag(2*n + 1)
}

// EulerNumber returns the n-th Euler (secant) number E(n). Odd-indexed values
// are zero; E(2m) = (-1)^m A(2m).
func EulerNumber(n int) *big.Int {
	if n < 0 || n%2 != 0 {
		return big.NewInt(0)
	}
	e := combinZigzag(n)
	if (n/2)%2 == 1 {
		e.Neg(e)
	}
	return e
}

// EulerNumberSequence returns the Euler numbers E(0), E(1), ..., E(n).
func EulerNumberSequence(n int) []*big.Int {
	if n < 0 {
		return nil
	}
	seq := make([]*big.Int, n+1)
	for i := 0; i <= n; i++ {
		seq[i] = EulerNumber(i)
	}
	return seq
}

// ----------------------------------------------------------------------------
// Necklaces and Lyndon words
// ----------------------------------------------------------------------------

// NecklaceNumber returns the number of distinct necklaces of length n using up
// to k colors, counted up to rotation: (1/n) sum_{d|n} phi(d) k^(n/d).
func NecklaceNumber(n, k int) *big.Int {
	if n <= 0 {
		return big.NewInt(0)
	}
	sum := big.NewInt(0)
	for _, d := range combinDivisors(n) {
		term := new(big.Int).Mul(big.NewInt(int64(combinTotient(d))), combinPowInt(k, n/d))
		sum.Add(sum, term)
	}
	sum.Div(sum, big.NewInt(int64(n)))
	return sum
}

// LyndonWordNumber returns the number of Lyndon words of length n over a k-letter
// alphabet: (1/n) sum_{d|n} mu(d) k^(n/d).
func LyndonWordNumber(n, k int) *big.Int {
	if n <= 0 {
		return big.NewInt(0)
	}
	sum := big.NewInt(0)
	for _, d := range combinDivisors(n) {
		term := new(big.Int).Mul(big.NewInt(int64(combinMobius(d))), combinPowInt(k, n/d))
		sum.Add(sum, term)
	}
	sum.Div(sum, big.NewInt(int64(n)))
	return sum
}
