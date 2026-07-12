package ntheory

import "math/big"

// Fibonacci returns the n-th Fibonacci number Fₙ, with F₀ = 0, F₁ = 1 and
// Fₙ = Fₙ₋₁ + Fₙ₋₂, as an arbitrary-precision integer. n must be non-negative.
func Fibonacci(n int64) *big.Int {
	if n < 0 {
		panic("ntheory: Fibonacci requires n >= 0")
	}
	a, b := big.NewInt(0), big.NewInt(1)
	for i := int64(0); i < n; i++ {
		a.Add(a, b)
		a, b = b, a
	}
	return a
}

// Lucas returns the n-th Lucas number Lₙ, with L₀ = 2, L₁ = 1 and
// Lₙ = Lₙ₋₁ + Lₙ₋₂, as an arbitrary-precision integer. n must be non-negative.
func Lucas(n int64) *big.Int {
	if n < 0 {
		panic("ntheory: Lucas requires n >= 0")
	}
	a, b := big.NewInt(2), big.NewInt(1)
	for i := int64(0); i < n; i++ {
		a.Add(a, b)
		a, b = b, a
	}
	return a
}

// Tribonacci returns the n-th Tribonacci number Tₙ, with T₀ = 0, T₁ = 1,
// T₂ = 1 and Tₙ = Tₙ₋₁ + Tₙ₋₂ + Tₙ₋₃, as an arbitrary-precision integer. n must
// be non-negative.
func Tribonacci(n int64) *big.Int {
	if n < 0 {
		panic("ntheory: Tribonacci requires n >= 0")
	}
	a, b, c := big.NewInt(0), big.NewInt(1), big.NewInt(1)
	if n == 0 {
		return a
	}
	if n == 1 || n == 2 {
		return big.NewInt(1)
	}
	for i := int64(3); i <= n; i++ {
		next := new(big.Int).Add(a, b)
		next.Add(next, c)
		a, b, c = b, c, next
	}
	return c
}

// IsqrtBig returns the integer square root of n: the greatest integer r with
// r*r <= n, using math/big's exact Sqrt. n must be non-negative.
func IsqrtBig(n *big.Int) *big.Int {
	if n.Sign() < 0 {
		panic("ntheory: IsqrtBig requires n >= 0")
	}
	return new(big.Int).Sqrt(n)
}

// IsSquare reports whether n is a perfect square (n == k*k for some integer k).
// Negative numbers are never perfect squares; 0 and 1 are.
func IsSquare(n int64) bool {
	if n < 0 {
		return false
	}
	r := IsqrtBig(big.NewInt(n))
	r.Mul(r, r)
	return r.Int64() == n
}

// Bernoulli returns the n-th Bernoulli number Bₙ as an exact rational
// (math/big.Rat), using the convention B₁ = -1/2. All odd-indexed Bernoulli
// numbers beyond B₁ are zero. n must be non-negative.
//
// It is computed from the recurrence Σ_{k=0}^{n} C(n+1, k)·Bₖ = 0.
func Bernoulli(n int64) *big.Rat {
	if n < 0 {
		panic("ntheory: Bernoulli requires n >= 0")
	}
	b := make([]*big.Rat, n+1)
	for m := int64(0); m <= n; m++ {
		// B_m = -1/(m+1) * Σ_{k=0}^{m-1} C(m+1, k) B_k, with B_0 = 1.
		b[m] = new(big.Rat).SetInt64(0)
		if m == 0 {
			b[m] = new(big.Rat).SetInt64(1)
			continue
		}
		sum := new(big.Rat)
		for k := int64(0); k < m; k++ {
			c := new(big.Rat).SetInt(new(big.Int).Binomial(m+1, k))
			term := new(big.Rat).Mul(c, b[k])
			sum.Add(sum, term)
		}
		// B_m = -sum / (m+1).
		denom := new(big.Rat).SetInt64(m + 1)
		sum.Quo(sum, denom)
		sum.Neg(sum)
		b[m] = sum
	}
	return b[n]
}
