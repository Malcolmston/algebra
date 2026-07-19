package graphspectral

import (
	"math"
	"sort"
)

// Dot returns the Euclidean inner product of a and b. It returns 0 when the
// slices have different lengths.
func Dot(a, b []float64) float64 {
	if len(a) != len(b) {
		return 0
	}
	var s float64
	for i := range a {
		s += a[i] * b[i]
	}
	return s
}

// Norm2 returns the Euclidean (L2) norm of v.
func Norm2(v []float64) float64 {
	return math.Sqrt(Dot(v, v))
}

// Norm1 returns the L1 (taxicab) norm of v, the sum of absolute values.
func Norm1(v []float64) float64 {
	var s float64
	for _, x := range v {
		s += math.Abs(x)
	}
	return s
}

// NormInf returns the L-infinity (maximum absolute value) norm of v.
func NormInf(v []float64) float64 {
	var m float64
	for _, x := range v {
		if a := math.Abs(x); a > m {
			m = a
		}
	}
	return m
}

// Normalize returns a new unit-length copy of v under the Euclidean norm. A zero
// vector is returned unchanged (as a copy).
func Normalize(v []float64) []float64 {
	n := Norm2(v)
	out := make([]float64, len(v))
	if n == 0 {
		copy(out, v)
		return out
	}
	for i, x := range v {
		out[i] = x / n
	}
	return out
}

// NormalizeL1 returns a new copy of v scaled so that its L1 norm is one. A vector
// whose entries sum in absolute value to zero is returned unchanged.
func NormalizeL1(v []float64) []float64 {
	n := Norm1(v)
	out := make([]float64, len(v))
	if n == 0 {
		copy(out, v)
		return out
	}
	for i, x := range v {
		out[i] = x / n
	}
	return out
}

// VecAdd returns the element-wise sum a+b. It returns nil on a length mismatch.
func VecAdd(a, b []float64) []float64 {
	if len(a) != len(b) {
		return nil
	}
	out := make([]float64, len(a))
	for i := range a {
		out[i] = a[i] + b[i]
	}
	return out
}

// VecSub returns the element-wise difference a-b. It returns nil on a length
// mismatch.
func VecSub(a, b []float64) []float64 {
	if len(a) != len(b) {
		return nil
	}
	out := make([]float64, len(a))
	for i := range a {
		out[i] = a[i] - b[i]
	}
	return out
}

// VecScale returns a new vector equal to s*v.
func VecScale(v []float64, s float64) []float64 {
	out := make([]float64, len(v))
	for i, x := range v {
		out[i] = s * x
	}
	return out
}

// VecAXPY returns a*x + y (a scalar times x, plus y). It returns nil on a length
// mismatch.
func VecAXPY(a float64, x, y []float64) []float64 {
	if len(x) != len(y) {
		return nil
	}
	out := make([]float64, len(x))
	for i := range x {
		out[i] = a*x[i] + y[i]
	}
	return out
}

// VecSum returns the sum of the entries of v.
func VecSum(v []float64) float64 {
	var s float64
	for _, x := range v {
		s += x
	}
	return s
}

// VecMean returns the arithmetic mean of the entries of v, or 0 for an empty
// slice.
func VecMean(v []float64) float64 {
	if len(v) == 0 {
		return 0
	}
	return VecSum(v) / float64(len(v))
}

// VecMax returns the maximum entry of v, or 0 for an empty slice.
func VecMax(v []float64) float64 {
	if len(v) == 0 {
		return 0
	}
	m := v[0]
	for _, x := range v[1:] {
		if x > m {
			m = x
		}
	}
	return m
}

// VecMin returns the minimum entry of v, or 0 for an empty slice.
func VecMin(v []float64) float64 {
	if len(v) == 0 {
		return 0
	}
	m := v[0]
	for _, x := range v[1:] {
		if x < m {
			m = x
		}
	}
	return m
}

// ArgMax returns the index of the (first) maximum entry of v, or -1 if v is
// empty.
func ArgMax(v []float64) int {
	if len(v) == 0 {
		return -1
	}
	idx := 0
	for i, x := range v {
		if x > v[idx] {
			idx = i
		}
	}
	return idx
}

// ArgMin returns the index of the (first) minimum entry of v, or -1 if v is
// empty.
func ArgMin(v []float64) int {
	if len(v) == 0 {
		return -1
	}
	idx := 0
	for i, x := range v {
		if x < v[idx] {
			idx = i
		}
	}
	return idx
}

// VecClone returns a copy of v.
func VecClone(v []float64) []float64 {
	out := make([]float64, len(v))
	copy(out, v)
	return out
}

// Zeros returns a new zero vector of length n.
func Zeros(n int) []float64 { return make([]float64, n) }

// Ones returns a new vector of length n whose entries are all one.
func Ones(n int) []float64 {
	out := make([]float64, n)
	for i := range out {
		out[i] = 1
	}
	return out
}

// VecEqual reports whether a and b are element-wise identical.
func VecEqual(a, b []float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// VecApproxEqual reports whether a and b agree to within absolute tolerance tol
// in every coordinate.
func VecApproxEqual(a, b []float64, tol float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if math.Abs(a[i]-b[i]) > tol {
			return false
		}
	}
	return true
}

// VecDistance returns the Euclidean distance between a and b, or +Inf on a
// length mismatch.
func VecDistance(a, b []float64) float64 {
	if len(a) != len(b) {
		return math.Inf(1)
	}
	var s float64
	for i := range a {
		d := a[i] - b[i]
		s += d * d
	}
	return math.Sqrt(s)
}

// CosineSimilarity returns the cosine of the angle between a and b. It returns 0
// if either vector is zero or the lengths differ.
func CosineSimilarity(a, b []float64) float64 {
	if len(a) != len(b) {
		return 0
	}
	na, nb := Norm2(a), Norm2(b)
	if na == 0 || nb == 0 {
		return 0
	}
	return Dot(a, b) / (na * nb)
}

// SortedFloats returns a sorted (ascending) copy of v.
func SortedFloats(v []float64) []float64 {
	out := VecClone(v)
	sort.Float64s(out)
	return out
}
