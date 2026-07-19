package odesolvers

import "math"

// Field is the right-hand side of a first-order ODE system y' = f(t, y).
//
// A Field must return a newly allocated slice of the same length as y and must
// not modify y. This aliasing-free contract lets the integrators reuse the
// returned slices freely.
type Field func(t float64, y []float64) []float64

// Clone returns a copy of v. The returned slice shares no storage with v.
func Clone(v []float64) []float64 {
	out := make([]float64, len(v))
	copy(out, v)
	return out
}

// Zeros returns a newly allocated slice of n zeros.
func Zeros(n int) []float64 { return make([]float64, n) }

// Fill returns a newly allocated slice of length n with every element set to x.
func Fill(n int, x float64) []float64 {
	out := make([]float64, n)
	for i := range out {
		out[i] = x
	}
	return out
}

// Add returns the elementwise sum a+b. It panics if the lengths differ.
func Add(a, b []float64) []float64 {
	mustSameLen(a, b)
	out := make([]float64, len(a))
	for i := range a {
		out[i] = a[i] + b[i]
	}
	return out
}

// Sub returns the elementwise difference a-b. It panics if the lengths differ.
func Sub(a, b []float64) []float64 {
	mustSameLen(a, b)
	out := make([]float64, len(a))
	for i := range a {
		out[i] = a[i] - b[i]
	}
	return out
}

// Scale returns the vector s*v.
func Scale(s float64, v []float64) []float64 {
	out := make([]float64, len(v))
	for i := range v {
		out[i] = s * v[i]
	}
	return out
}

// Neg returns the elementwise negation -v.
func Neg(v []float64) []float64 { return Scale(-1, v) }

// AXPY returns a fresh vector a*x + y (the classic "a x plus y"). It panics if
// the lengths of x and y differ.
func AXPY(a float64, x, y []float64) []float64 {
	mustSameLen(x, y)
	out := make([]float64, len(x))
	for i := range x {
		out[i] = a*x[i] + y[i]
	}
	return out
}

// AXPYInPlace performs y += a*x, mutating and returning y. It panics if the
// lengths differ.
func AXPYInPlace(a float64, x, y []float64) []float64 {
	mustSameLen(x, y)
	for i := range x {
		y[i] += a * x[i]
	}
	return y
}

// LinearCombination returns sum_i coeffs[i] * vecs[i]. All vectors must have the
// same length n, which must be supplied so an empty combination is well defined.
func LinearCombination(n int, coeffs []float64, vecs [][]float64) []float64 {
	if len(coeffs) != len(vecs) {
		panic("odesolvers: LinearCombination coeff/vector count mismatch")
	}
	out := make([]float64, n)
	for k, c := range coeffs {
		v := vecs[k]
		if len(v) != n {
			panic("odesolvers: LinearCombination dimension mismatch")
		}
		for i := 0; i < n; i++ {
			out[i] += c * v[i]
		}
	}
	return out
}

// Dot returns the Euclidean inner product of a and b.
func Dot(a, b []float64) float64 {
	mustSameLen(a, b)
	var s float64
	for i := range a {
		s += a[i] * b[i]
	}
	return s
}

// Norm2 returns the Euclidean (L2) norm of v.
func Norm2(v []float64) float64 { return math.Sqrt(Dot(v, v)) }

// Norm1 returns the L1 norm (sum of absolute values) of v.
func Norm1(v []float64) float64 {
	var s float64
	for _, x := range v {
		s += math.Abs(x)
	}
	return s
}

// NormInf returns the maximum-absolute-value (infinity) norm of v.
func NormInf(v []float64) float64 {
	var m float64
	for _, x := range v {
		if a := math.Abs(x); a > m {
			m = a
		}
	}
	return m
}

// RMSNorm returns the root-mean-square norm sqrt(mean(v_i^2)) of v. For an
// empty vector it returns 0.
func RMSNorm(v []float64) float64 {
	if len(v) == 0 {
		return 0
	}
	return Norm2(v) / math.Sqrt(float64(len(v)))
}

// Distance returns the Euclidean distance between a and b.
func Distance(a, b []float64) float64 { return Norm2(Sub(a, b)) }

// MaxAbs returns the largest absolute value among the elements of v, or 0 for
// an empty slice.
func MaxAbs(v []float64) float64 { return NormInf(v) }

// Linspace returns n evenly spaced samples over the closed interval [a, b].
// For n <= 0 it returns an empty slice; for n == 1 it returns [a].
func Linspace(a, b float64, n int) []float64 {
	if n <= 0 {
		return []float64{}
	}
	if n == 1 {
		return []float64{a}
	}
	out := make([]float64, n)
	step := (b - a) / float64(n-1)
	for i := 0; i < n; i++ {
		out[i] = a + step*float64(i)
	}
	out[n-1] = b
	return out
}

// WeightedRMSNorm returns the scaled RMS norm sqrt(mean((v_i/scale_i)^2)),
// the error measure used by the adaptive integrators. It panics if the lengths
// differ.
func WeightedRMSNorm(v, scale []float64) float64 {
	mustSameLen(v, scale)
	if len(v) == 0 {
		return 0
	}
	var s float64
	for i := range v {
		r := v[i] / scale[i]
		s += r * r
	}
	return math.Sqrt(s / float64(len(v)))
}

// mustSameLen panics if a and b have different lengths.
func mustSameLen(a, b []float64) {
	if len(a) != len(b) {
		panic("odesolvers: vector length mismatch")
	}
}
