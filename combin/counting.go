// Package combin implements combinatorics: counting, arrangements and
// combinatorial numbers.
//
// The package collects the classical counting functions used throughout
// discrete mathematics. Where results grow without bound the routines return
// arbitrary-precision *big.Int values so that answers stay exact; parallel
// float64 and natural-logarithm variants are provided for the cases where an
// approximate but fast answer (or one that survives magnitudes far beyond the
// range of a machine integer) is preferable.
//
// The following families are covered:
//
//   - factorials: Factorial, DoubleFactorial, Multifactorial, Subfactorial
//     (derangements), Superfactorial, HyperFactorial and Primorial;
//   - arrangements: Permutations (nPr) and Permutations with repetition;
//   - selections: Combinations (nCr), Binomial, combinations with repetition,
//     the central binomial coefficient and the real-argument GammaBinomial;
//   - Multinomial coefficients;
//   - Rising and Falling factorials (Pochhammer symbols);
//   - the combinatorial number triangles and sequences: Catalan, Motzkin,
//     Schroeder, Bell, Stirling numbers of the first and second kind, Lah
//     numbers, Narayana, Eulerian, Delannoy and telephone numbers;
//   - Pascal's triangle rows and integer-partition counts.
//
// Every routine is deterministic and depends only on the Go standard library.
package combin

import (
	"math"
	"math/big"
)

// combinlgamma returns the natural logarithm of the absolute value of the
// Gamma function evaluated at x. It is the shared engine behind the package's
// logarithmic routines.
func combinlgamma(x float64) float64 {
	v, _ := math.Lgamma(x)
	return v
}

// combincheckNonNeg panics with the given message when n is negative. It keeps
// the exported routines terse while giving callers a clear diagnostic.
func combincheckNonNeg(n int, msg string) {
	if n < 0 {
		panic("combin: " + msg)
	}
}

// -------------------------------------------------------------------------
// Factorials
// -------------------------------------------------------------------------

// Factorial returns n! = 1·2·…·n as an exact integer. Factorial(0) is 1.
// It panics when n is negative.
func Factorial(n int) *big.Int {
	combincheckNonNeg(n, "Factorial of negative number")
	r := big.NewInt(1)
	t := new(big.Int)
	for i := 2; i <= n; i++ {
		r.Mul(r, t.SetInt64(int64(i)))
	}
	return r
}

// FactorialFloat returns n! as a float64, using the Gamma function so that the
// result stays finite up to n == 170. It panics when n is negative.
func FactorialFloat(n int) float64 {
	combincheckNonNeg(n, "FactorialFloat of negative number")
	return math.Gamma(float64(n) + 1)
}

// LogFactorial returns the natural logarithm of n!, computed from the log-Gamma
// function so that it is accurate for very large n. It panics when n is
// negative.
func LogFactorial(n int) float64 {
	combincheckNonNeg(n, "LogFactorial of negative number")
	return combinlgamma(float64(n) + 1)
}

// DoubleFactorial returns the double factorial n!! = n·(n-2)·(n-4)·…, stopping
// at 2 or 1. By convention 0!! = (-1)!! = 1. It panics when n < -1.
func DoubleFactorial(n int) *big.Int {
	if n < -1 {
		panic("combin: DoubleFactorial undefined for n < -1")
	}
	r := big.NewInt(1)
	t := new(big.Int)
	for i := n; i > 1; i -= 2 {
		r.Mul(r, t.SetInt64(int64(i)))
	}
	return r
}

// DoubleFactorialFloat returns the double factorial n!! as a float64.
// It panics when n < -1.
func DoubleFactorialFloat(n int) float64 {
	if n < -1 {
		panic("combin: DoubleFactorialFloat undefined for n < -1")
	}
	r := 1.0
	for i := n; i > 1; i -= 2 {
		r *= float64(i)
	}
	return r
}

// Multifactorial returns the k-fold factorial n·(n-k)·(n-2k)·… down to the
// smallest positive term. Multifactorial(n, 1) equals Factorial(n) and
// Multifactorial(n, 2) equals DoubleFactorial(n). It panics when k < 1 or
// n < 0.
func Multifactorial(n, k int) *big.Int {
	if k < 1 {
		panic("combin: Multifactorial requires k >= 1")
	}
	combincheckNonNeg(n, "Multifactorial of negative number")
	r := big.NewInt(1)
	t := new(big.Int)
	for i := n; i > 0; i -= k {
		r.Mul(r, t.SetInt64(int64(i)))
	}
	return r
}

// Subfactorial returns the number of derangements !n, i.e. the number of
// permutations of n elements with no fixed point. It satisfies
// !n = (n-1)·(!(n-1) + !(n-2)) with !0 = 1 and !1 = 0. It panics when n < 0.
func Subfactorial(n int) *big.Int {
	combincheckNonNeg(n, "Subfactorial of negative number")
	if n == 0 {
		return big.NewInt(1)
	}
	if n == 1 {
		return big.NewInt(0)
	}
	a := big.NewInt(1) // !(i-2)
	b := big.NewInt(0) // !(i-1)
	t := new(big.Int)
	for i := 2; i <= n; i++ {
		t.Add(a, b)
		t.Mul(t, big.NewInt(int64(i-1)))
		a.Set(b)
		b.Set(t)
	}
	return new(big.Int).Set(b)
}

// SubfactorialFloat returns the number of derangements !n as a float64, using
// the closed form round(n!/e). It panics when n is negative.
func SubfactorialFloat(n int) float64 {
	combincheckNonNeg(n, "SubfactorialFloat of negative number")
	if n == 0 {
		return 1
	}
	return math.Round(math.Gamma(float64(n)+1) / math.E)
}

// Derangement is a synonym for Subfactorial: it returns the number of
// permutations of n elements that leave no element fixed. It panics when n < 0.
func Derangement(n int) *big.Int {
	return Subfactorial(n)
}

// Superfactorial returns the product of the first n factorials,
// 1!·2!·…·n!, with Superfactorial(0) = 1. It panics when n is negative.
func Superfactorial(n int) *big.Int {
	combincheckNonNeg(n, "Superfactorial of negative number")
	result := big.NewInt(1)
	fact := big.NewInt(1)
	t := new(big.Int)
	for k := 1; k <= n; k++ {
		fact.Mul(fact, t.SetInt64(int64(k)))
		result.Mul(result, fact)
	}
	return result
}

// HyperFactorial returns the hyperfactorial H(n) = 1^1·2^2·…·n^n, with
// HyperFactorial(0) = 1. It panics when n is negative.
func HyperFactorial(n int) *big.Int {
	combincheckNonNeg(n, "HyperFactorial of negative number")
	result := big.NewInt(1)
	for k := 1; k <= n; k++ {
		kk := big.NewInt(int64(k))
		result.Mul(result, new(big.Int).Exp(kk, kk, nil))
	}
	return result
}

// Primorial returns n# — the product of all prime numbers not exceeding n.
// Primorial(0) and Primorial(1) are 1. It panics when n is negative.
func Primorial(n int) *big.Int {
	combincheckNonNeg(n, "Primorial of negative number")
	if n < 2 {
		return big.NewInt(1)
	}
	composite := make([]bool, n+1)
	r := big.NewInt(1)
	for i := 2; i <= n; i++ {
		if !composite[i] {
			r.Mul(r, big.NewInt(int64(i)))
			for j := i * 2; j <= n; j += i {
				composite[j] = true
			}
		}
	}
	return r
}

// -------------------------------------------------------------------------
// Arrangements (permutations)
// -------------------------------------------------------------------------

// Permutations returns the number of ordered arrangements of r items chosen
// from n, nPr = n!/(n-r)!. It returns 0 when r > n and panics when either
// argument is negative.
func Permutations(n, r int) *big.Int {
	combincheckNonNeg(n, "Permutations of negative n")
	combincheckNonNeg(r, "Permutations of negative r")
	if r > n {
		return big.NewInt(0)
	}
	result := big.NewInt(1)
	t := new(big.Int)
	for i := 0; i < r; i++ {
		result.Mul(result, t.SetInt64(int64(n-i)))
	}
	return result
}

// PermutationsFloat returns nPr as a float64, evaluated through logarithms so
// that it remains finite for large arguments. It returns 0 when r > n and
// panics when either argument is negative.
func PermutationsFloat(n, r int) float64 {
	combincheckNonNeg(n, "PermutationsFloat of negative n")
	combincheckNonNeg(r, "PermutationsFloat of negative r")
	if r > n {
		return 0
	}
	v := math.Exp(LogPermutations(n, r))
	if v < 9e15 {
		return math.Round(v)
	}
	return v
}

// LogPermutations returns the natural logarithm of nPr. It panics when either
// argument is negative or when r > n (for which the logarithm is undefined).
func LogPermutations(n, r int) float64 {
	combincheckNonNeg(n, "LogPermutations of negative n")
	combincheckNonNeg(r, "LogPermutations of negative r")
	if r > n {
		panic("combin: LogPermutations requires r <= n")
	}
	return LogFactorial(n) - LogFactorial(n-r)
}

// PermutationsWithRepetition returns n^r, the number of ordered sequences of
// length r drawn from n symbols when repetition is allowed. It panics when
// either argument is negative.
func PermutationsWithRepetition(n, r int) *big.Int {
	combincheckNonNeg(n, "PermutationsWithRepetition of negative n")
	combincheckNonNeg(r, "PermutationsWithRepetition of negative r")
	return new(big.Int).Exp(big.NewInt(int64(n)), big.NewInt(int64(r)), nil)
}

// -------------------------------------------------------------------------
// Selections (combinations / binomial coefficients)
// -------------------------------------------------------------------------

// Binomial returns the binomial coefficient C(n, k) = n!/(k!·(n-k)!) using a
// multiplicative scheme that keeps every intermediate value integral. It
// returns 0 when k is outside [0, n] and panics when n is negative.
func Binomial(n, k int) *big.Int {
	combincheckNonNeg(n, "Binomial of negative n")
	if k < 0 || k > n {
		return big.NewInt(0)
	}
	if k > n-k {
		k = n - k
	}
	result := big.NewInt(1)
	t := new(big.Int)
	for i := 0; i < k; i++ {
		result.Mul(result, t.SetInt64(int64(n-i)))
		result.Div(result, t.SetInt64(int64(i+1)))
	}
	return result
}

// Combinations returns the number of unordered selections of k items from n,
// nCr = C(n, k). It is a synonym for Binomial. It returns 0 when k is outside
// [0, n] and panics when n is negative.
func Combinations(n, r int) *big.Int {
	return Binomial(n, r)
}

// CombinationsFloat returns C(n, k) as a float64, evaluated via logarithms so
// that it stays finite for large arguments. It returns 0 when k is outside
// [0, n] and panics when n is negative.
func CombinationsFloat(n, r int) float64 {
	return BinomialFloat(n, r)
}

// BinomialFloat returns the binomial coefficient C(n, k) as a float64. Small
// results are rounded to the nearest integer; very large results are returned
// as their floating-point approximation. It panics when n is negative.
func BinomialFloat(n, k int) float64 {
	combincheckNonNeg(n, "BinomialFloat of negative n")
	if k < 0 || k > n {
		return 0
	}
	v := math.Exp(LogBinomial(n, k))
	if v < 9e15 {
		return math.Round(v)
	}
	return v
}

// LogBinomial returns the natural logarithm of the binomial coefficient
// C(n, k). It panics when n is negative or k is outside [0, n].
func LogBinomial(n, k int) float64 {
	combincheckNonNeg(n, "LogBinomial of negative n")
	if k < 0 || k > n {
		panic("combin: LogBinomial requires 0 <= k <= n")
	}
	return LogFactorial(n) - LogFactorial(k) - LogFactorial(n-k)
}

// LogCombinations returns the natural logarithm of C(n, k). It is a synonym for
// LogBinomial. It panics when n is negative or k is outside [0, n].
func LogCombinations(n, r int) float64 {
	return LogBinomial(n, r)
}

// CombinationsWithRepetition returns the number of multisets of size r drawn
// from n distinct items, C(n+r-1, r). It panics when either argument is
// negative.
func CombinationsWithRepetition(n, r int) *big.Int {
	combincheckNonNeg(n, "CombinationsWithRepetition of negative n")
	combincheckNonNeg(r, "CombinationsWithRepetition of negative r")
	if n == 0 {
		if r == 0 {
			return big.NewInt(1)
		}
		return big.NewInt(0)
	}
	return Binomial(n+r-1, r)
}

// CentralBinomial returns the central binomial coefficient C(2n, n).
// It panics when n is negative.
func CentralBinomial(n int) *big.Int {
	combincheckNonNeg(n, "CentralBinomial of negative n")
	return Binomial(2*n, n)
}

// GammaBinomial returns the generalized binomial coefficient for real
// arguments, C(x, y) = Γ(x+1) / (Γ(y+1)·Γ(x-y+1)). For non-negative integer
// arguments it agrees with the ordinary binomial coefficient.
func GammaBinomial(x, y float64) float64 {
	return math.Gamma(x+1) / (math.Gamma(y+1) * math.Gamma(x-y+1))
}

// -------------------------------------------------------------------------
// Multinomial coefficients
// -------------------------------------------------------------------------

// Multinomial returns the multinomial coefficient (k1+k2+…)! / (k1!·k2!·…),
// the number of ways to partition a set of that many items into labelled groups
// of the given sizes. It panics when any argument is negative.
func Multinomial(ks ...int) *big.Int {
	result := big.NewInt(1)
	n := 0
	for _, k := range ks {
		combincheckNonNeg(k, "Multinomial of negative part")
		n += k
		result.Mul(result, Binomial(n, k))
	}
	return result
}

// MultinomialFloat returns the multinomial coefficient as a float64, evaluated
// through logarithms. It panics when any argument is negative.
func MultinomialFloat(ks ...int) float64 {
	v := math.Exp(LogMultinomial(ks...))
	if v < 9e15 {
		return math.Round(v)
	}
	return v
}

// LogMultinomial returns the natural logarithm of the multinomial coefficient
// for the given part sizes. It panics when any argument is negative.
func LogMultinomial(ks ...int) float64 {
	total := 0
	for _, k := range ks {
		combincheckNonNeg(k, "LogMultinomial of negative part")
		total += k
	}
	r := LogFactorial(total)
	for _, k := range ks {
		r -= LogFactorial(k)
	}
	return r
}

// -------------------------------------------------------------------------
// Rising and falling factorials (Pochhammer symbols)
// -------------------------------------------------------------------------

// RisingFactorial returns the rising factorial (Pochhammer symbol)
// x^(n) = x·(x+1)·…·(x+n-1), with the empty product RisingFactorial(x, 0) = 1.
// It panics when n is negative.
func RisingFactorial(x, n int) *big.Int {
	combincheckNonNeg(n, "RisingFactorial of negative length")
	r := big.NewInt(1)
	t := new(big.Int)
	for i := 0; i < n; i++ {
		r.Mul(r, t.SetInt64(int64(x+i)))
	}
	return r
}

// FallingFactorial returns the falling factorial
// (x)_n = x·(x-1)·…·(x-n+1), with the empty product FallingFactorial(x, 0) = 1.
// It panics when n is negative.
func FallingFactorial(x, n int) *big.Int {
	combincheckNonNeg(n, "FallingFactorial of negative length")
	r := big.NewInt(1)
	t := new(big.Int)
	for i := 0; i < n; i++ {
		r.Mul(r, t.SetInt64(int64(x-i)))
	}
	return r
}

// RisingFactorialFloat returns the rising factorial x^(n) for a real base x.
// It panics when n is negative.
func RisingFactorialFloat(x float64, n int) float64 {
	combincheckNonNeg(n, "RisingFactorialFloat of negative length")
	r := 1.0
	for i := 0; i < n; i++ {
		r *= x + float64(i)
	}
	return r
}

// FallingFactorialFloat returns the falling factorial (x)_n for a real base x.
// It panics when n is negative.
func FallingFactorialFloat(x float64, n int) float64 {
	combincheckNonNeg(n, "FallingFactorialFloat of negative length")
	r := 1.0
	for i := 0; i < n; i++ {
		r *= x - float64(i)
	}
	return r
}

// -------------------------------------------------------------------------
// Catalan, Motzkin and Schroeder numbers
// -------------------------------------------------------------------------

// Catalan returns the nth Catalan number, C(2n, n)/(n+1). The sequence begins
// 1, 1, 2, 5, 14, 42, … It panics when n is negative.
func Catalan(n int) *big.Int {
	combincheckNonNeg(n, "Catalan of negative n")
	c := Binomial(2*n, n)
	c.Div(c, big.NewInt(int64(n+1)))
	return c
}

// CatalanFloat returns the nth Catalan number as a float64.
// It panics when n is negative.
func CatalanFloat(n int) float64 {
	combincheckNonNeg(n, "CatalanFloat of negative n")
	return BinomialFloat(2*n, n) / float64(n+1)
}

// CatalanTriangle returns the entry C(n, k) of Catalan's triangle, the number
// of lattice paths that never rise above the diagonal, equal to
// binom(n+k, k)·(n-k+1)/(n+1). Row n sums to the (n+1)th Catalan number.
// It returns 0 when k is outside [0, n] and panics when n is negative.
func CatalanTriangle(n, k int) *big.Int {
	combincheckNonNeg(n, "CatalanTriangle of negative n")
	if k < 0 || k > n {
		return big.NewInt(0)
	}
	r := Binomial(n+k, k)
	r.Mul(r, big.NewInt(int64(n-k+1)))
	r.Div(r, big.NewInt(int64(n+1)))
	return r
}

// MotzkinNumber returns the nth Motzkin number, counting the ways of drawing
// non-intersecting chords between n points on a circle. The sequence begins
// 1, 1, 2, 4, 9, 21, 51, … It panics when n is negative.
func MotzkinNumber(n int) *big.Int {
	combincheckNonNeg(n, "MotzkinNumber of negative n")
	m := make([]*big.Int, n+1)
	m[0] = big.NewInt(1)
	if n >= 1 {
		m[1] = big.NewInt(1)
	}
	for i := 2; i <= n; i++ {
		v := new(big.Int).Set(m[i-1])
		for k := 0; k <= i-2; k++ {
			v.Add(v, new(big.Int).Mul(m[k], m[i-2-k]))
		}
		m[i] = v
	}
	return new(big.Int).Set(m[n])
}

// SchroederNumber returns the nth large Schröder number, counting the lattice
// paths from (0,0) to (n,n) using steps east, north and north-east that never
// rise above the diagonal. The sequence begins 1, 2, 6, 22, 90, 394, …
// It panics when n is negative.
func SchroederNumber(n int) *big.Int {
	combincheckNonNeg(n, "SchroederNumber of negative n")
	if n == 0 {
		return big.NewInt(1)
	}
	if n == 1 {
		return big.NewInt(2)
	}
	a := big.NewInt(1)
	b := big.NewInt(2)
	for i := 2; i <= n; i++ {
		c := new(big.Int).Mul(big.NewInt(int64(6*i-3)), b)
		c.Sub(c, new(big.Int).Mul(big.NewInt(int64(i-2)), a))
		c.Div(c, big.NewInt(int64(i+1)))
		a.Set(b)
		b.Set(c)
	}
	return new(big.Int).Set(b)
}

// -------------------------------------------------------------------------
// Bell numbers
// -------------------------------------------------------------------------

// Bell returns the nth Bell number, the number of partitions of a set of n
// elements. The sequence begins 1, 1, 2, 5, 15, 52, 203, … It panics when n
// is negative.
func Bell(n int) *big.Int {
	combincheckNonNeg(n, "Bell of negative n")
	return BellTriangleRow(n)[0]
}

// BellTriangleRow returns row n of the Bell (Aitken) triangle. Its first entry
// is the nth Bell number and its last entry is the (n+1)th Bell number.
// It panics when n is negative.
func BellTriangleRow(n int) []*big.Int {
	combincheckNonNeg(n, "BellTriangleRow of negative n")
	row := []*big.Int{big.NewInt(1)}
	for i := 1; i <= n; i++ {
		next := make([]*big.Int, i+1)
		next[0] = new(big.Int).Set(row[len(row)-1])
		for j := 1; j <= i; j++ {
			next[j] = new(big.Int).Add(next[j-1], row[j-1])
		}
		row = next
	}
	return row
}

// -------------------------------------------------------------------------
// Stirling numbers
// -------------------------------------------------------------------------

// StirlingFirstUnsigned returns the unsigned Stirling number of the first kind
// c(n, k), the number of permutations of n elements with exactly k disjoint
// cycles. It satisfies c(n, k) = c(n-1, k-1) + (n-1)·c(n-1, k) with c(0, 0) = 1.
// It returns 0 for arguments outside the valid triangle.
func StirlingFirstUnsigned(n, k int) *big.Int {
	if n < 0 || k < 0 || k > n {
		return big.NewInt(0)
	}
	if n == 0 {
		if k == 0 {
			return big.NewInt(1)
		}
		return big.NewInt(0)
	}
	s := make([]*big.Int, k+1)
	for j := range s {
		s[j] = big.NewInt(0)
	}
	s[0] = big.NewInt(1)
	for i := 1; i <= n; i++ {
		ns := make([]*big.Int, k+1)
		ns[0] = big.NewInt(0)
		for j := 1; j <= k; j++ {
			term := new(big.Int).Mul(big.NewInt(int64(i-1)), s[j])
			term.Add(term, s[j-1])
			ns[j] = term
		}
		s = ns
	}
	return s[k]
}

// StirlingFirst returns the signed Stirling number of the first kind
// s(n, k) = (-1)^(n-k)·c(n, k), the coefficients relating falling factorials to
// ordinary powers. It returns 0 for arguments outside the valid triangle.
func StirlingFirst(n, k int) *big.Int {
	u := StirlingFirstUnsigned(n, k)
	if (n-k)%2 != 0 {
		u.Neg(u)
	}
	return u
}

// StirlingFirstRow returns the unsigned Stirling numbers of the first kind
// c(n, 0), c(n, 1), …, c(n, n). It panics when n is negative.
func StirlingFirstRow(n int) []*big.Int {
	combincheckNonNeg(n, "StirlingFirstRow of negative n")
	row := make([]*big.Int, n+1)
	for k := 0; k <= n; k++ {
		row[k] = StirlingFirstUnsigned(n, k)
	}
	return row
}

// StirlingSecond returns the Stirling number of the second kind S(n, k), the
// number of ways to partition n elements into k non-empty unlabelled subsets.
// It satisfies S(n, k) = k·S(n-1, k) + S(n-1, k-1) with S(0, 0) = 1. It returns
// 0 for arguments outside the valid triangle.
func StirlingSecond(n, k int) *big.Int {
	if n < 0 || k < 0 || k > n {
		return big.NewInt(0)
	}
	if n == 0 {
		if k == 0 {
			return big.NewInt(1)
		}
		return big.NewInt(0)
	}
	s := make([]*big.Int, k+1)
	for j := range s {
		s[j] = big.NewInt(0)
	}
	s[0] = big.NewInt(1)
	for i := 1; i <= n; i++ {
		ns := make([]*big.Int, k+1)
		ns[0] = big.NewInt(0)
		for j := 1; j <= k; j++ {
			term := new(big.Int).Mul(big.NewInt(int64(j)), s[j])
			term.Add(term, s[j-1])
			ns[j] = term
		}
		s = ns
	}
	return s[k]
}

// StirlingSecondRow returns the Stirling numbers of the second kind
// S(n, 0), S(n, 1), …, S(n, n). It panics when n is negative.
func StirlingSecondRow(n int) []*big.Int {
	combincheckNonNeg(n, "StirlingSecondRow of negative n")
	row := make([]*big.Int, n+1)
	for k := 0; k <= n; k++ {
		row[k] = StirlingSecond(n, k)
	}
	return row
}

// -------------------------------------------------------------------------
// Lah numbers
// -------------------------------------------------------------------------

// LahNumber returns the unsigned Lah number L(n, k) = C(n-1, k-1)·n!/k!, the
// number of ways to partition n elements into k non-empty ordered lists.
// L(0, 0) is 1. It returns 0 for arguments outside the valid triangle.
func LahNumber(n, k int) *big.Int {
	if n == 0 && k == 0 {
		return big.NewInt(1)
	}
	if n < 1 || k < 1 || k > n {
		return big.NewInt(0)
	}
	res := Binomial(n-1, k-1)
	num := big.NewInt(1)
	t := new(big.Int)
	for i := k + 1; i <= n; i++ {
		num.Mul(num, t.SetInt64(int64(i)))
	}
	res.Mul(res, num)
	return res
}

// SignedLahNumber returns the signed Lah number (-1)^n·L(n, k), the
// coefficients that express rising factorials in terms of falling factorials.
// It returns 0 for arguments outside the valid triangle.
func SignedLahNumber(n, k int) *big.Int {
	l := LahNumber(n, k)
	if n%2 != 0 {
		l.Neg(l)
	}
	return l
}

// -------------------------------------------------------------------------
// Pascal's triangle
// -------------------------------------------------------------------------

// PascalRow returns row n of Pascal's triangle, the binomial coefficients
// C(n, 0), C(n, 1), …, C(n, n). It panics when n is negative.
func PascalRow(n int) []*big.Int {
	combincheckNonNeg(n, "PascalRow of negative n")
	row := make([]*big.Int, n+1)
	row[0] = big.NewInt(1)
	for k := 1; k <= n; k++ {
		v := new(big.Int).Mul(row[k-1], big.NewInt(int64(n-k+1)))
		v.Div(v, big.NewInt(int64(k)))
		row[k] = v
	}
	return row
}

// BinomialRow is a synonym for PascalRow: it returns the binomial coefficients
// C(n, 0), C(n, 1), …, C(n, n). It panics when n is negative.
func BinomialRow(n int) []*big.Int {
	return PascalRow(n)
}

// PascalTriangle returns rows 0 through n of Pascal's triangle.
// It panics when n is negative.
func PascalTriangle(n int) [][]*big.Int {
	combincheckNonNeg(n, "PascalTriangle of negative n")
	t := make([][]*big.Int, n+1)
	for i := 0; i <= n; i++ {
		t[i] = PascalRow(i)
	}
	return t
}

// -------------------------------------------------------------------------
// Linear-recurrence sequences
// -------------------------------------------------------------------------

// Fibonacci returns the nth Fibonacci number, F(0) = 0, F(1) = 1 and
// F(n) = F(n-1) + F(n-2). It panics when n is negative.
func Fibonacci(n int) *big.Int {
	combincheckNonNeg(n, "Fibonacci of negative n")
	if n == 0 {
		return big.NewInt(0)
	}
	a := big.NewInt(0)
	b := big.NewInt(1)
	t := new(big.Int)
	for i := 2; i <= n; i++ {
		t.Add(a, b)
		a.Set(b)
		b.Set(t)
	}
	return new(big.Int).Set(b)
}

// Lucas returns the nth Lucas number, L(0) = 2, L(1) = 1 and
// L(n) = L(n-1) + L(n-2). It panics when n is negative.
func Lucas(n int) *big.Int {
	combincheckNonNeg(n, "Lucas of negative n")
	if n == 0 {
		return big.NewInt(2)
	}
	a := big.NewInt(2)
	b := big.NewInt(1)
	t := new(big.Int)
	for i := 2; i <= n; i++ {
		t.Add(a, b)
		a.Set(b)
		b.Set(t)
	}
	return new(big.Int).Set(b)
}

// Tribonacci returns the nth tribonacci number with T(0) = T(1) = 0, T(2) = 1
// and T(n) = T(n-1) + T(n-2) + T(n-3). The sequence begins
// 0, 0, 1, 1, 2, 4, 7, 13, … It panics when n is negative.
func Tribonacci(n int) *big.Int {
	combincheckNonNeg(n, "Tribonacci of negative n")
	if n < 2 {
		return big.NewInt(0)
	}
	if n == 2 {
		return big.NewInt(1)
	}
	a := big.NewInt(0)
	b := big.NewInt(0)
	c := big.NewInt(1)
	t := new(big.Int)
	for i := 3; i <= n; i++ {
		t.Add(a, b)
		t.Add(t, c)
		a.Set(b)
		b.Set(c)
		c.Set(t)
	}
	return new(big.Int).Set(c)
}

// TelephoneNumber returns the nth telephone (involution) number, the number of
// self-inverse permutations of n elements. It satisfies
// T(n) = T(n-1) + (n-1)·T(n-2) with T(0) = T(1) = 1 and begins
// 1, 1, 2, 4, 10, 26, 76, … It panics when n is negative.
func TelephoneNumber(n int) *big.Int {
	combincheckNonNeg(n, "TelephoneNumber of negative n")
	if n == 0 || n == 1 {
		return big.NewInt(1)
	}
	a := big.NewInt(1)
	b := big.NewInt(1)
	t := new(big.Int)
	for i := 2; i <= n; i++ {
		t.Mul(big.NewInt(int64(i-1)), a)
		t.Add(t, b)
		a.Set(b)
		b.Set(t)
	}
	return new(big.Int).Set(b)
}

// -------------------------------------------------------------------------
// Partitions and compositions
// -------------------------------------------------------------------------

// PartitionCount returns p(n), the number of unordered ways to write n as a sum
// of positive integers, computed with Euler's pentagonal-number recurrence. The
// sequence begins 1, 1, 2, 3, 5, 7, 11, … It panics when n is negative.
func PartitionCount(n int) *big.Int {
	combincheckNonNeg(n, "PartitionCount of negative n")
	p := make([]*big.Int, n+1)
	p[0] = big.NewInt(1)
	for i := 1; i <= n; i++ {
		sum := big.NewInt(0)
		for k := 1; ; k++ {
			g1 := k * (3*k - 1) / 2
			g2 := k * (3*k + 1) / 2
			if g1 > i && g2 > i {
				break
			}
			neg := k%2 == 0
			if g1 <= i {
				if neg {
					sum.Sub(sum, p[i-g1])
				} else {
					sum.Add(sum, p[i-g1])
				}
			}
			if g2 <= i {
				if neg {
					sum.Sub(sum, p[i-g2])
				} else {
					sum.Add(sum, p[i-g2])
				}
			}
		}
		p[i] = sum
	}
	return new(big.Int).Set(p[n])
}

// PartitionsIntoKParts returns the number of partitions of n into exactly k
// positive parts, using the recurrence p(n, k) = p(n-1, k-1) + p(n-k, k) with
// p(0, 0) = 1. It returns 0 for arguments outside the valid range and panics
// when either argument is negative.
func PartitionsIntoKParts(n, k int) *big.Int {
	combincheckNonNeg(n, "PartitionsIntoKParts of negative n")
	combincheckNonNeg(k, "PartitionsIntoKParts of negative k")
	if k > n {
		return big.NewInt(0)
	}
	if n == 0 && k == 0 {
		return big.NewInt(1)
	}
	if k == 0 {
		return big.NewInt(0)
	}
	dp := make([][]*big.Int, n+1)
	for i := range dp {
		dp[i] = make([]*big.Int, k+1)
		for j := range dp[i] {
			dp[i][j] = big.NewInt(0)
		}
	}
	dp[0][0] = big.NewInt(1)
	for i := 1; i <= n; i++ {
		for j := 1; j <= k && j <= i; j++ {
			dp[i][j] = new(big.Int).Add(dp[i-1][j-1], dp[i-j][j])
		}
	}
	return new(big.Int).Set(dp[n][k])
}

// Compositions returns the number of ordered compositions of n into positive
// parts, 2^(n-1) for n >= 1 and 1 for n = 0. It panics when n is negative.
func Compositions(n int) *big.Int {
	combincheckNonNeg(n, "Compositions of negative n")
	if n == 0 {
		return big.NewInt(1)
	}
	return new(big.Int).Lsh(big.NewInt(1), uint(n-1))
}

// WeakCompositions returns the number of ways to write n as an ordered sum of k
// non-negative integers, C(n+k-1, k-1). It panics when either argument is
// negative.
func WeakCompositions(n, k int) *big.Int {
	combincheckNonNeg(n, "WeakCompositions of negative n")
	combincheckNonNeg(k, "WeakCompositions of negative k")
	if k == 0 {
		if n == 0 {
			return big.NewInt(1)
		}
		return big.NewInt(0)
	}
	return Binomial(n+k-1, k-1)
}

// -------------------------------------------------------------------------
// Combinatorial triangles
// -------------------------------------------------------------------------

// NarayanaNumber returns the Narayana number N(n, k) = C(n, k)·C(n, k-1)/n,
// counting Dyck paths of length 2n with exactly k peaks. Row n sums to the nth
// Catalan number. It returns 0 for arguments outside [1, n].
func NarayanaNumber(n, k int) *big.Int {
	if n <= 0 || k < 1 || k > n {
		return big.NewInt(0)
	}
	r := new(big.Int).Mul(Binomial(n, k), Binomial(n, k-1))
	r.Div(r, big.NewInt(int64(n)))
	return r
}

// EulerianNumber returns the Eulerian number A(n, k), the number of
// permutations of n elements with exactly k ascents. It satisfies
// A(n, k) = (k+1)·A(n-1, k) + (n-k)·A(n-1, k-1). It returns 0 for arguments
// outside the valid triangle and panics when n is negative.
func EulerianNumber(n, k int) *big.Int {
	combincheckNonNeg(n, "EulerianNumber of negative n")
	if n == 0 {
		if k == 0 {
			return big.NewInt(1)
		}
		return big.NewInt(0)
	}
	if k < 0 || k >= n {
		return big.NewInt(0)
	}
	row := []*big.Int{big.NewInt(1)} // n == 1
	for i := 2; i <= n; i++ {
		nr := make([]*big.Int, i)
		for j := 0; j < i; j++ {
			term := big.NewInt(0)
			if j < len(row) {
				term.Add(term, new(big.Int).Mul(big.NewInt(int64(j+1)), row[j]))
			}
			if j-1 >= 0 && j-1 < len(row) {
				term.Add(term, new(big.Int).Mul(big.NewInt(int64(i-j)), row[j-1]))
			}
			nr[j] = term
		}
		row = nr
	}
	return new(big.Int).Set(row[k])
}

// EulerianSecondOrder returns the second-order Eulerian number <<n, k>>,
// satisfying <<n, k>> = (k+1)·<<n-1, k>> + (2n-1-k)·<<n-1, k-1>> with
// <<0, 0>> = 1. These count Stirling permutations by descents. It returns 0 for
// arguments outside the valid triangle and panics when n is negative.
func EulerianSecondOrder(n, k int) *big.Int {
	combincheckNonNeg(n, "EulerianSecondOrder of negative n")
	if n == 0 {
		if k == 0 {
			return big.NewInt(1)
		}
		return big.NewInt(0)
	}
	if k < 0 || k >= n {
		return big.NewInt(0)
	}
	row := []*big.Int{big.NewInt(1)} // n == 1
	for i := 2; i <= n; i++ {
		nr := make([]*big.Int, i)
		for j := 0; j < i; j++ {
			term := big.NewInt(0)
			if j < len(row) {
				term.Add(term, new(big.Int).Mul(big.NewInt(int64(j+1)), row[j]))
			}
			if j-1 >= 0 && j-1 < len(row) {
				term.Add(term, new(big.Int).Mul(big.NewInt(int64(2*i-1-j)), row[j-1]))
			}
			nr[j] = term
		}
		row = nr
	}
	return new(big.Int).Set(row[k])
}

// DelannoyNumber returns the Delannoy number D(m, n), the number of lattice
// paths from (0,0) to (m,n) using steps east, north and north-east. It
// satisfies D(m, n) = D(m-1, n) + D(m, n-1) + D(m-1, n-1). It panics when
// either argument is negative.
func DelannoyNumber(m, n int) *big.Int {
	combincheckNonNeg(m, "DelannoyNumber of negative m")
	combincheckNonNeg(n, "DelannoyNumber of negative n")
	dp := make([][]*big.Int, m+1)
	for i := 0; i <= m; i++ {
		dp[i] = make([]*big.Int, n+1)
		for j := 0; j <= n; j++ {
			if i == 0 || j == 0 {
				dp[i][j] = big.NewInt(1)
				continue
			}
			v := new(big.Int).Add(dp[i-1][j], dp[i][j-1])
			v.Add(v, dp[i-1][j-1])
			dp[i][j] = v
		}
	}
	return new(big.Int).Set(dp[m][n])
}

// RencontresNumber returns the number of permutations of n elements with
// exactly k fixed points, C(n, k)·!(n-k). Row n sums to n!. It returns 0 for
// arguments outside [0, n] and panics when either argument is negative.
func RencontresNumber(n, k int) *big.Int {
	combincheckNonNeg(n, "RencontresNumber of negative n")
	combincheckNonNeg(k, "RencontresNumber of negative k")
	if k > n {
		return big.NewInt(0)
	}
	r := Binomial(n, k)
	r.Mul(r, Subfactorial(n-k))
	return r
}

// -------------------------------------------------------------------------
// Harmonic numbers and asymptotics
// -------------------------------------------------------------------------

// HarmonicNumber returns the nth harmonic number H(n) = 1 + 1/2 + … + 1/n as a
// float64, with H(0) = 0. It panics when n is negative.
func HarmonicNumber(n int) float64 {
	combincheckNonNeg(n, "HarmonicNumber of negative n")
	h := 0.0
	for i := n; i >= 1; i-- {
		h += 1 / float64(i)
	}
	return h
}

// HarmonicNumberGeneralized returns the generalized harmonic number
// H(n, m) = 1 + 1/2^m + … + 1/n^m as a float64, with H(0, m) = 0. It panics
// when n is negative.
func HarmonicNumberGeneralized(n, m int) float64 {
	combincheckNonNeg(n, "HarmonicNumberGeneralized of negative n")
	h := 0.0
	for i := n; i >= 1; i-- {
		h += math.Pow(float64(i), -float64(m))
	}
	return h
}

// StirlingApproximation returns Stirling's asymptotic estimate of n!,
// sqrt(2·π·n)·(n/e)^n, as a float64. StirlingApproximation(0) is 1. It panics
// when n is negative.
func StirlingApproximation(n int) float64 {
	combincheckNonNeg(n, "StirlingApproximation of negative n")
	if n == 0 {
		return 1
	}
	nf := float64(n)
	return math.Sqrt(2*math.Pi*nf) * math.Pow(nf/math.E, nf)
}
