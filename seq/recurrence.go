package seq

import (
	"math/big"
	"math/bits"
)

// seqFibPair returns the pair (F(n), F(n+1)) of consecutive Fibonacci numbers
// using the fast-doubling identities
//
//	F(2k)   = F(k) * (2*F(k+1) - F(k))
//	F(2k+1) = F(k)^2 + F(k+1)^2
//
// The result is exact for every n whose (n+1)-th Fibonacci number fits in a
// uint64; the caller is responsible for the overflow range.
func seqFibPair(n int) (uint64, uint64) {
	var a, b uint64 = 0, 1 // F(0), F(1)
	for i := bits.Len(uint(n)) - 1; i >= 0; i-- {
		c := a * (2*b - a) // F(2k)
		d := a*a + b*b     // F(2k+1)
		if (n>>uint(i))&1 == 0 {
			a, b = c, d
		} else {
			a, b = d, c+d
		}
	}
	return a, b
}

// Fibonacci returns the n-th Fibonacci number Fₙ, defined by F₀ = 0, F₁ = 1 and
// Fₙ = Fₙ₋₁ + Fₙ₋₂. It is computed by the fast-doubling method in O(log n)
// arithmetic operations. n must be non-negative. The result is exact for
// n ≤ 93; beyond that Fₙ exceeds the range of uint64 and the value wraps.
func Fibonacci(n int) uint64 {
	if n < 0 {
		panic("seq: Fibonacci requires n >= 0")
	}
	a, _ := seqFibPair(n)
	return a
}

// FibonacciBig returns the n-th Fibonacci number Fₙ as an arbitrary-precision
// integer, computed by the fast-doubling method in O(log n) big-integer
// multiplications. n must be non-negative.
func FibonacciBig(n int) *big.Int {
	if n < 0 {
		panic("seq: FibonacciBig requires n >= 0")
	}
	a, _ := seqFibPairBig(n)
	return a
}

// seqFibPairBig is the arbitrary-precision counterpart of seqFibPair.
func seqFibPairBig(n int) (*big.Int, *big.Int) {
	a := big.NewInt(0) // F(0)
	b := big.NewInt(1) // F(1)
	for i := bits.Len(uint(n)) - 1; i >= 0; i-- {
		// c = a*(2b - a)
		t := new(big.Int).Lsh(b, 1)
		t.Sub(t, a)
		c := new(big.Int).Mul(a, t)
		// d = a^2 + b^2
		a2 := new(big.Int).Mul(a, a)
		b2 := new(big.Int).Mul(b, b)
		d := a2.Add(a2, b2)
		if (n>>uint(i))&1 == 0 {
			a.Set(c)
			b.Set(d)
		} else {
			a.Set(d)
			b.Add(c, d)
		}
	}
	return a, b
}

// FibonacciSequence returns the first n Fibonacci numbers F₀ … Fₙ₋₁. n must be
// non-negative; a value of 0 yields an empty slice. Values are exact for
// indices up to 93.
func FibonacciSequence(n int) []uint64 {
	return seqLinear2(0, 1, n)
}

// Lucas returns the n-th Lucas number Lₙ, defined by L₀ = 2, L₁ = 1 and
// Lₙ = Lₙ₋₁ + Lₙ₋₂. It is derived from the fast-doubling Fibonacci pair via the
// identity Lₙ = 2·Fₙ₊₁ − Fₙ. n must be non-negative. The result is exact for
// n ≤ 91; beyond that Lₙ exceeds the range of uint64.
func Lucas(n int) uint64 {
	if n < 0 {
		panic("seq: Lucas requires n >= 0")
	}
	a, b := seqFibPair(n)
	return 2*b - a
}

// LucasBig returns the n-th Lucas number Lₙ as an arbitrary-precision integer.
// n must be non-negative.
func LucasBig(n int) *big.Int {
	if n < 0 {
		panic("seq: LucasBig requires n >= 0")
	}
	a, b := seqFibPairBig(n)
	res := new(big.Int).Lsh(b, 1) // 2*F(n+1)
	return res.Sub(res, a)        // 2*F(n+1) - F(n)
}

// LucasSequence returns the first n Lucas numbers L₀ … Lₙ₋₁. n must be
// non-negative.
func LucasSequence(n int) []uint64 {
	return seqLinear2(2, 1, n)
}

// Pell returns the n-th Pell number Pₙ, defined by P₀ = 0, P₁ = 1 and
// Pₙ = 2·Pₙ₋₁ + Pₙ₋₂. n must be non-negative. The result is exact for n ≤ 50;
// beyond that Pₙ exceeds the range of uint64.
func Pell(n int) uint64 {
	if n < 0 {
		panic("seq: Pell requires n >= 0")
	}
	if n == 0 {
		return 0
	}
	var a, b uint64 = 0, 1
	for i := 1; i < n; i++ {
		a, b = b, 2*b+a
	}
	return b
}

// PellBig returns the n-th Pell number Pₙ as an arbitrary-precision integer.
// n must be non-negative.
func PellBig(n int) *big.Int {
	if n < 0 {
		panic("seq: PellBig requires n >= 0")
	}
	a, b := big.NewInt(0), big.NewInt(1)
	if n == 0 {
		return a
	}
	for i := 1; i < n; i++ {
		next := new(big.Int).Lsh(b, 1) // 2b
		next.Add(next, a)              // 2b + a
		a, b = b, next
	}
	return b
}

// PellSequence returns the first n Pell numbers P₀ … Pₙ₋₁. n must be
// non-negative.
func PellSequence(n int) []uint64 {
	if n < 0 {
		panic("seq: PellSequence requires n >= 0")
	}
	out := make([]uint64, 0, n)
	var a, b uint64 = 0, 1
	for i := 0; i < n; i++ {
		out = append(out, a)
		a, b = b, 2*b+a
	}
	return out
}

// PellLucas returns the n-th Pell-Lucas (companion Pell) number Qₙ, defined by
// Q₀ = 2, Q₁ = 2 and Qₙ = 2·Qₙ₋₁ + Qₙ₋₂: 2, 2, 6, 14, 34, 82, … n must be
// non-negative. The result is exact for n ≤ 50.
func PellLucas(n int) uint64 {
	if n < 0 {
		panic("seq: PellLucas requires n >= 0")
	}
	var a, b uint64 = 2, 2
	if n == 0 {
		return a
	}
	for i := 1; i < n; i++ {
		a, b = b, 2*b+a
	}
	return b
}

// Jacobsthal returns the n-th Jacobsthal number Jₙ, defined by J₀ = 0, J₁ = 1
// and Jₙ = Jₙ₋₁ + 2·Jₙ₋₂: 0, 1, 1, 3, 5, 11, 21, … equivalently
// Jₙ = (2ⁿ − (−1)ⁿ)/3. n must be non-negative. The result is exact for n ≤ 63.
func Jacobsthal(n int) uint64 {
	if n < 0 {
		panic("seq: Jacobsthal requires n >= 0")
	}
	var a, b uint64 = 0, 1
	if n == 0 {
		return a
	}
	for i := 1; i < n; i++ {
		a, b = b, b+2*a
	}
	return b
}

// JacobsthalLucas returns the n-th Jacobsthal-Lucas number, defined by the same
// recurrence with initial values 2 and 1: 2, 1, 5, 7, 17, 31, … equivalently
// 2ⁿ + (−1)ⁿ. n must be non-negative. The result is exact for n ≤ 63.
func JacobsthalLucas(n int) uint64 {
	if n < 0 {
		panic("seq: JacobsthalLucas requires n >= 0")
	}
	var a, b uint64 = 2, 1
	if n == 0 {
		return a
	}
	for i := 1; i < n; i++ {
		a, b = b, b+2*a
	}
	return b
}

// JacobsthalSequence returns the first n Jacobsthal numbers J₀ … Jₙ₋₁. n must
// be non-negative.
func JacobsthalSequence(n int) []uint64 {
	if n < 0 {
		panic("seq: JacobsthalSequence requires n >= 0")
	}
	out := make([]uint64, 0, n)
	var a, b uint64 = 0, 1
	for i := 0; i < n; i++ {
		out = append(out, a)
		a, b = b, b+2*a
	}
	return out
}

// Tribonacci returns the n-th Tribonacci number Tₙ, defined by T₀ = 0, T₁ = 0,
// T₂ = 1 and Tₙ = Tₙ₋₁ + Tₙ₋₂ + Tₙ₋₃: 0, 0, 1, 1, 2, 4, 7, 13, 24, 44, …
// n must be non-negative. The result is exact for n ≤ 91.
func Tribonacci(n int) uint64 {
	if n < 0 {
		panic("seq: Tribonacci requires n >= 0")
	}
	var a, b, c uint64 = 0, 0, 1
	if n < 3 {
		return [...]uint64{0, 0, 1}[n]
	}
	for i := 3; i <= n; i++ {
		a, b, c = b, c, a+b+c
	}
	return c
}

// TribonacciSequence returns the first n Tribonacci numbers T₀ … Tₙ₋₁. n must
// be non-negative.
func TribonacciSequence(n int) []uint64 {
	if n < 0 {
		panic("seq: TribonacciSequence requires n >= 0")
	}
	out := make([]uint64, 0, n)
	var a, b, c uint64 = 0, 0, 1
	for i := 0; i < n; i++ {
		out = append(out, a)
		a, b, c = b, c, a+b+c
	}
	return out
}

// Tetranacci returns the n-th Tetranacci number, defined by initial values
// 0, 0, 0, 1 and the four-term sum recurrence: 0, 0, 0, 1, 1, 2, 4, 8, 15, 29,
// 56, … n must be non-negative. The result is exact for n ≤ 88.
func Tetranacci(n int) uint64 {
	if n < 0 {
		panic("seq: Tetranacci requires n >= 0")
	}
	if n < 4 {
		return [...]uint64{0, 0, 0, 1}[n]
	}
	var a, b, c, d uint64 = 0, 0, 0, 1
	for i := 4; i <= n; i++ {
		a, b, c, d = b, c, d, a+b+c+d
	}
	return d
}

// Padovan returns the n-th Padovan number, defined by P₀ = P₁ = P₂ = 1 and
// Pₙ = Pₙ₋₂ + Pₙ₋₃: 1, 1, 1, 2, 2, 3, 4, 5, 7, 9, 12, 16, … n must be
// non-negative. The result is exact for n ≤ 219.
func Padovan(n int) uint64 {
	if n < 0 {
		panic("seq: Padovan requires n >= 0")
	}
	if n < 3 {
		return 1
	}
	var a, b, c uint64 = 1, 1, 1 // P(0), P(1), P(2)
	for i := 3; i <= n; i++ {
		// P(i) = P(i-2) + P(i-3) = b + a
		a, b, c = b, c, a+b
	}
	return c
}

// PadovanSequence returns the first n Padovan numbers P₀ … Pₙ₋₁. n must be
// non-negative.
func PadovanSequence(n int) []uint64 {
	if n < 0 {
		panic("seq: PadovanSequence requires n >= 0")
	}
	out := make([]uint64, 0, n)
	var a, b, c uint64 = 1, 1, 1
	for i := 0; i < n; i++ {
		out = append(out, a)
		a, b, c = b, c, a+b
	}
	return out
}

// Perrin returns the n-th Perrin number, defined by P₀ = 3, P₁ = 0, P₂ = 2 and
// Pₙ = Pₙ₋₂ + Pₙ₋₃: 3, 0, 2, 3, 2, 5, 5, 7, 10, 12, 17, 22, … n must be
// non-negative. The result is exact for n ≤ 218.
func Perrin(n int) uint64 {
	if n < 0 {
		panic("seq: Perrin requires n >= 0")
	}
	if n < 3 {
		return [...]uint64{3, 0, 2}[n]
	}
	var a, b, c uint64 = 3, 0, 2 // P(0), P(1), P(2)
	for i := 3; i <= n; i++ {
		a, b, c = b, c, a+b
	}
	return c
}

// PerrinSequence returns the first n Perrin numbers P₀ … Pₙ₋₁. n must be
// non-negative.
func PerrinSequence(n int) []uint64 {
	if n < 0 {
		panic("seq: PerrinSequence requires n >= 0")
	}
	out := make([]uint64, 0, n)
	var a, b, c uint64 = 3, 0, 2
	for i := 0; i < n; i++ {
		out = append(out, a)
		a, b, c = b, c, a+b
	}
	return out
}

// seqLinear2 returns the first n terms of the order-2 recurrence
// x(i) = x(i-1) + x(i-2) with the given two initial values.
func seqLinear2(x0, x1 uint64, n int) []uint64 {
	if n < 0 {
		panic("seq: sequence length must be >= 0")
	}
	out := make([]uint64, 0, n)
	a, b := x0, x1
	for i := 0; i < n; i++ {
		out = append(out, a)
		a, b = b, a+b
	}
	return out
}

// LinearRecurrence evaluates the term aₙ of the constant-coefficient linear
// recurrence
//
//	a(i) = coeffs[0]·a(i-1) + coeffs[1]·a(i-2) + … + coeffs[k-1]·a(i-k)
//
// of order k = len(coeffs), whose first k terms a₀ … a₍ₖ₋₁₎ are supplied in
// init. coeffs and init must be non-empty and of equal length, and n must be
// non-negative. Arithmetic is performed in int64, so negative terms (as in the
// Perrin sequence) are handled naturally; the caller is responsible for the
// overflow range.
func LinearRecurrence(coeffs, init []int64, n int) int64 {
	k := len(coeffs)
	if k == 0 || k != len(init) {
		panic("seq: LinearRecurrence requires non-empty coeffs and init of equal length")
	}
	if n < 0 {
		panic("seq: LinearRecurrence requires n >= 0")
	}
	if n < k {
		return init[n]
	}
	window := make([]int64, k)
	copy(window, init) // window[j] holds a(i-k+j)
	for i := k; i <= n; i++ {
		var next int64
		for j := 0; j < k; j++ {
			// coeffs[j] multiplies a(i-1-j) = window[k-1-j]
			next += coeffs[j] * window[k-1-j]
		}
		copy(window, window[1:])
		window[k-1] = next
	}
	return window[k-1]
}

// LinearRecurrenceSequence returns the first m terms a₀ … a₍ₘ₋₁₎ of the linear
// recurrence described by LinearRecurrence. coeffs and init must be non-empty
// and of equal length, and m must be non-negative.
func LinearRecurrenceSequence(coeffs, init []int64, m int) []int64 {
	k := len(coeffs)
	if k == 0 || k != len(init) {
		panic("seq: LinearRecurrenceSequence requires non-empty coeffs and init of equal length")
	}
	if m < 0 {
		panic("seq: LinearRecurrenceSequence requires m >= 0")
	}
	out := make([]int64, 0, m)
	window := make([]int64, k)
	copy(window, init)
	for i := 0; i < m; i++ {
		if i < k {
			out = append(out, init[i])
			continue
		}
		var next int64
		for j := 0; j < k; j++ {
			next += coeffs[j] * window[k-1-j]
		}
		copy(window, window[1:])
		window[k-1] = next
		out = append(out, next)
	}
	return out
}

// LinearRecurrenceBig evaluates term aₙ of the constant-coefficient linear
// recurrence described by LinearRecurrence using arbitrary-precision integer
// arithmetic, so no overflow can occur. coeffs and init must be non-empty and
// of equal length, and n must be non-negative.
func LinearRecurrenceBig(coeffs, init []*big.Int, n int) *big.Int {
	k := len(coeffs)
	if k == 0 || k != len(init) {
		panic("seq: LinearRecurrenceBig requires non-empty coeffs and init of equal length")
	}
	if n < 0 {
		panic("seq: LinearRecurrenceBig requires n >= 0")
	}
	if n < k {
		return new(big.Int).Set(init[n])
	}
	window := make([]*big.Int, k)
	for j := range init {
		window[j] = new(big.Int).Set(init[j])
	}
	for i := k; i <= n; i++ {
		next := new(big.Int)
		term := new(big.Int)
		for j := 0; j < k; j++ {
			term.Mul(coeffs[j], window[k-1-j])
			next.Add(next, term)
		}
		copy(window, window[1:])
		window[k-1] = next
	}
	return new(big.Int).Set(window[k-1])
}

// IsFibonacci reports whether n is a Fibonacci number. Every non-negative
// integer that appears in the sequence 0, 1, 1, 2, 3, 5, 8, 13, … qualifies.
func IsFibonacci(n uint64) bool {
	var a, b uint64 = 0, 1
	for a < n {
		a, b = b, a+b
	}
	return a == n
}

// Zeckendorf returns the Zeckendorf representation of n: the unique set of
// non-consecutive Fibonacci numbers (each ≥ 1, drawn from 1, 2, 3, 5, 8, …)
// whose sum is n, returned in descending order. The representation of 0 is the
// empty slice.
func Zeckendorf(n uint64) []uint64 {
	if n == 0 {
		return []uint64{}
	}
	// Generate Fibonacci numbers 1, 2, 3, 5, 8, … not exceeding n.
	fibs := []uint64{1, 2}
	for {
		next := fibs[len(fibs)-1] + fibs[len(fibs)-2]
		if next > n {
			break
		}
		fibs = append(fibs, next)
	}
	var out []uint64
	for i := len(fibs) - 1; i >= 0; i-- {
		if fibs[i] <= n {
			out = append(out, fibs[i])
			n -= fibs[i]
			i-- // skip the immediately smaller Fibonacci number
		}
	}
	return out
}
