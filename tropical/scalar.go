package tropical

import "math"

// MinPlusZero returns the tropical zero of the min-plus semiring, +Inf.
func MinPlusZero() float64 { return math.Inf(1) }

// MaxPlusZero returns the tropical zero of the max-plus semiring, -Inf.
func MaxPlusZero() float64 { return math.Inf(-1) }

// MinPlusOne returns the tropical one of the min-plus semiring, 0.
func MinPlusOne() float64 { return 0 }

// MaxPlusOne returns the tropical one of the max-plus semiring, 0.
func MaxPlusOne() float64 { return 0 }

// MinPlusAdd returns the min-plus tropical sum of a and b, that is min(a, b).
func MinPlusAdd(a, b float64) float64 { return MinPlusSemiring().Add(a, b) }

// MaxPlusAdd returns the max-plus tropical sum of a and b, that is max(a, b).
func MaxPlusAdd(a, b float64) float64 { return MaxPlusSemiring().Add(a, b) }

// MinPlusMul returns the min-plus tropical product of a and b (a+b, with +Inf
// absorbing).
func MinPlusMul(a, b float64) float64 { return MinPlusSemiring().Mul(a, b) }

// MaxPlusMul returns the max-plus tropical product of a and b (a+b, with -Inf
// absorbing).
func MaxPlusMul(a, b float64) float64 { return MaxPlusSemiring().Mul(a, b) }

// MinPlusDiv returns the min-plus tropical quotient a/b (a-b on finite values).
func MinPlusDiv(a, b float64) float64 { return MinPlusSemiring().Div(a, b) }

// MaxPlusDiv returns the max-plus tropical quotient a/b (a-b on finite values).
func MaxPlusDiv(a, b float64) float64 { return MaxPlusSemiring().Div(a, b) }

// MinPlusPow returns the min-plus tropical power a^n (n*a on finite values).
func MinPlusPow(a float64, n int) float64 { return MinPlusSemiring().Pow(a, n) }

// MaxPlusPow returns the max-plus tropical power a^n (n*a on finite values).
func MaxPlusPow(a float64, n int) float64 { return MaxPlusSemiring().Pow(a, n) }

// MinPlusStar returns the min-plus scalar Kleene star: 0 for a >= 0 and -Inf
// otherwise.
func MinPlusStar(a float64) float64 { return MinPlusSemiring().Star(a) }

// MaxPlusStar returns the max-plus scalar Kleene star: 0 for a <= 0 and +Inf
// otherwise.
func MaxPlusStar(a float64) float64 { return MaxPlusSemiring().Star(a) }

// MinPlusSum returns the min-plus tropical sum (minimum) of the arguments, or
// +Inf when none are given.
func MinPlusSum(xs ...float64) float64 { return MinPlusSemiring().Sum(xs...) }

// MaxPlusSum returns the max-plus tropical sum (maximum) of the arguments, or
// -Inf when none are given.
func MaxPlusSum(xs ...float64) float64 { return MaxPlusSemiring().Sum(xs...) }

// MinPlusProd returns the min-plus tropical product (ordinary sum) of the
// arguments, or 0 when none are given.
func MinPlusProd(xs ...float64) float64 { return MinPlusSemiring().Prod(xs...) }

// MaxPlusProd returns the max-plus tropical product (ordinary sum) of the
// arguments, or 0 when none are given.
func MaxPlusProd(xs ...float64) float64 { return MaxPlusSemiring().Prod(xs...) }

// IsInf reports whether x is either positive or negative infinity, i.e. a
// tropical zero for one of the two semirings.
func IsInf(x float64) bool { return math.IsInf(x, 0) }
