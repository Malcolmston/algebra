package seq

import (
	"math/big"
	"sort"
	"strconv"
	"strings"
)

// CollatzNext returns the successor of n under the Collatz (3n+1) map: n/2 when
// n is even and 3n+1 when n is odd. n must be positive.
func CollatzNext(n uint64) uint64 {
	if n == 0 {
		panic("seq: CollatzNext requires n >= 1")
	}
	if n%2 == 0 {
		return n / 2
	}
	return 3*n + 1
}

// CollatzSequence returns the Collatz trajectory of n: the sequence starting at
// n and repeatedly applying CollatzNext until it reaches 1, inclusive of both
// endpoints. n must be positive. Intermediate values are assumed to stay within
// the range of uint64.
func CollatzSequence(n uint64) []uint64 {
	if n == 0 {
		panic("seq: CollatzSequence requires n >= 1")
	}
	out := []uint64{n}
	for n != 1 {
		n = CollatzNext(n)
		out = append(out, n)
	}
	return out
}

// CollatzSteps returns the total stopping time of n: the number of applications
// of the Collatz map required to reach 1. CollatzSteps(1) is 0. n must be
// positive.
func CollatzSteps(n uint64) int {
	if n == 0 {
		panic("seq: CollatzSteps requires n >= 1")
	}
	steps := 0
	for n != 1 {
		n = CollatzNext(n)
		steps++
	}
	return steps
}

// CollatzMax returns the largest value attained along the Collatz trajectory of
// n, including n itself. n must be positive.
func CollatzMax(n uint64) uint64 {
	if n == 0 {
		panic("seq: CollatzMax requires n >= 1")
	}
	max := n
	for n != 1 {
		n = CollatzNext(n)
		if n > max {
			max = n
		}
	}
	return max
}

// seqDigitSquareSum returns the sum of the squares of the base-10 digits of n.
func seqDigitSquareSum(n uint64) uint64 {
	var s uint64
	for n > 0 {
		d := n % 10
		s += d * d
		n /= 10
	}
	return s
}

// IsHappy reports whether n is a happy number: starting from n and repeatedly
// replacing the value by the sum of the squares of its decimal digits, the
// process reaches 1. Unhappy numbers instead enter the cycle that contains 4.
// n must be positive.
func IsHappy(n uint64) bool {
	if n == 0 {
		panic("seq: IsHappy requires n >= 1")
	}
	seen := make(map[uint64]bool)
	for n != 1 && !seen[n] {
		seen[n] = true
		n = seqDigitSquareSum(n)
	}
	return n == 1
}

// HappyNumbers returns the first count happy numbers in increasing order,
// beginning with 1, 7, 10, 13, 19, 23, 28, … count must be non-negative.
func HappyNumbers(count int) []uint64 {
	if count < 0 {
		panic("seq: HappyNumbers requires count >= 0")
	}
	out := make([]uint64, 0, count)
	for n := uint64(1); len(out) < count; n++ {
		if IsHappy(n) {
			out = append(out, n)
		}
	}
	return out
}

// IsKaprekar reports whether n is a Kaprekar number: the decimal representation
// of n² can be split into a left part and a non-zero right part that add up to
// n. For example 45² = 2025 splits as 20+25 = 45, and 4879² = 23804641 splits
// as 238+04641 = 4879. The value 1 qualifies by convention. n must be positive.
// This is the classical definition (OEIS A006886), which permits the split at
// any digit boundary.
func IsKaprekar(n uint64) bool {
	if n == 0 {
		panic("seq: IsKaprekar requires n >= 1")
	}
	if n == 1 {
		return true
	}
	nb := new(big.Int).SetUint64(n)
	sq := new(big.Int).Mul(nb, nb)
	s := sq.String()
	left := new(big.Int)
	right := new(big.Int)
	sum := new(big.Int)
	for k := 1; k < len(s); k++ {
		rightStr := s[len(s)-k:]
		right.SetString(rightStr, 10)
		if right.Sign() == 0 {
			continue // right part must be non-zero
		}
		left.SetString(s[:len(s)-k], 10)
		sum.Add(left, right)
		if sum.Cmp(nb) == 0 {
			return true
		}
	}
	return false
}

// KaprekarNumbers returns every Kaprekar number in the closed interval
// [1, limit], in increasing order: 1, 9, 45, 55, 99, 297, 703, 999, …
func KaprekarNumbers(limit uint64) []uint64 {
	var out []uint64
	for n := uint64(1); n <= limit; n++ {
		if IsKaprekar(n) {
			out = append(out, n)
		}
		if n == limit { // guard against overflow when limit == max uint64
			break
		}
	}
	return out
}

// KaprekarStep applies one step of Kaprekar's routine to a four-digit number:
// the digits of n (zero-padded to four places) are arranged into the largest
// and smallest numbers they form, and the smaller is subtracted from the
// larger. n must be in the range [0, 9999]. Applied repeatedly to any four-digit
// number that does not have all identical digits, the routine converges to
// Kaprekar's constant 6174.
func KaprekarStep(n int) int {
	if n < 0 || n > 9999 {
		panic("seq: KaprekarStep requires 0 <= n <= 9999")
	}
	digits := []int{n / 1000 % 10, n / 100 % 10, n / 10 % 10, n % 10}
	sort.Ints(digits)
	asc := digits[0]*1000 + digits[1]*100 + digits[2]*10 + digits[3]
	desc := digits[3]*1000 + digits[2]*100 + digits[1]*10 + digits[0]
	return desc - asc
}

// KaprekarSequence returns the trajectory of the four-digit Kaprekar routine
// starting from n, inclusive of n, continuing until it reaches the fixed point
// 6174, the fixed point 0 (for repdigits such as 1111), or a value that has
// already been seen. n must be in the range [0, 9999].
func KaprekarSequence(n int) []int {
	if n < 0 || n > 9999 {
		panic("seq: KaprekarSequence requires 0 <= n <= 9999")
	}
	out := []int{n}
	seen := map[int]bool{n: true}
	for n != 6174 && n != 0 {
		n = KaprekarStep(n)
		out = append(out, n)
		if seen[n] {
			break
		}
		seen[n] = true
	}
	return out
}

// RecamanSequence returns the first n terms of Recaman's sequence, defined by
// a₀ = 0 and, for k ≥ 1,
//
//	a(k) = a(k-1) − k  if that value is positive and not already in the sequence,
//	a(k) = a(k-1) + k  otherwise.
//
// The sequence begins 0, 1, 3, 6, 2, 7, 13, 20, 12, 21, 11, … n must be
// non-negative.
func RecamanSequence(n int) []int64 {
	if n < 0 {
		panic("seq: RecamanSequence requires n >= 0")
	}
	out := make([]int64, 0, n)
	seen := make(map[int64]bool)
	var cur int64
	for k := 0; k < n; k++ {
		if k == 0 {
			cur = 0
		} else {
			back := cur - int64(k)
			if back > 0 && !seen[back] {
				cur = back
			} else {
				cur = cur + int64(k)
			}
		}
		out = append(out, cur)
		seen[cur] = true
	}
	return out
}

// LookAndSayStep returns the look-and-say successor of the decimal digit string
// s: each maximal run of identical digits is replaced by its length followed by
// the digit. For example "1211" becomes "111221". s must be non-empty and
// consist only of decimal digits.
func LookAndSayStep(s string) string {
	if s == "" {
		panic("seq: LookAndSayStep requires a non-empty string")
	}
	for _, r := range s {
		if r < '0' || r > '9' {
			panic("seq: LookAndSayStep requires a decimal digit string")
		}
	}
	var b strings.Builder
	runeS := []byte(s)
	i := 0
	for i < len(runeS) {
		j := i
		for j < len(runeS) && runeS[j] == runeS[i] {
			j++
		}
		b.WriteString(strconv.Itoa(j - i))
		b.WriteByte(runeS[i])
		i = j
	}
	return b.String()
}

// LookAndSaySequence returns the first n terms of the look-and-say sequence
// generated from the seed string, inclusive of the seed. Each subsequent term
// is LookAndSayStep applied to the previous one. seed must be a non-empty
// decimal digit string and n must be non-negative. With seed "1" the sequence
// is "1", "11", "21", "1211", "111221", …
func LookAndSaySequence(seed string, n int) []string {
	if n < 0 {
		panic("seq: LookAndSaySequence requires n >= 0")
	}
	if n > 0 && seed == "" {
		panic("seq: LookAndSaySequence requires a non-empty seed")
	}
	out := make([]string, 0, n)
	cur := seed
	for i := 0; i < n; i++ {
		out = append(out, cur)
		cur = LookAndSayStep(cur)
	}
	return out
}
