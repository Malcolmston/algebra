package spectralpde

import "math"

// VectorZeros returns a zero vector of length n.
func VectorZeros(n int) []float64 { return make([]float64, n) }

// VectorFill returns a vector of length n with every entry set to v.
func VectorFill(n int, v float64) []float64 {
	out := make([]float64, n)
	for i := range out {
		out[i] = v
	}
	return out
}

// VectorCopy returns a copy of x.
func VectorCopy(x []float64) []float64 {
	out := make([]float64, len(x))
	copy(out, x)
	return out
}

// VectorAdd returns x+y.
func VectorAdd(x, y []float64) []float64 {
	out := make([]float64, len(x))
	for i := range x {
		out[i] = x[i] + y[i]
	}
	return out
}

// VectorSub returns x-y.
func VectorSub(x, y []float64) []float64 {
	out := make([]float64, len(x))
	for i := range x {
		out[i] = x[i] - y[i]
	}
	return out
}

// VectorScale returns s*x.
func VectorScale(x []float64, s float64) []float64 {
	out := make([]float64, len(x))
	for i := range x {
		out[i] = s * x[i]
	}
	return out
}

// VectorHadamard returns the elementwise product of x and y.
func VectorHadamard(x, y []float64) []float64 {
	out := make([]float64, len(x))
	for i := range x {
		out[i] = x[i] * y[i]
	}
	return out
}

// AXPY returns a*x + y.
func AXPY(a float64, x, y []float64) []float64 {
	out := make([]float64, len(x))
	for i := range x {
		out[i] = a*x[i] + y[i]
	}
	return out
}

// DotProduct returns the Euclidean inner product of x and y.
func DotProduct(x, y []float64) float64 {
	var s float64
	for i := range x {
		s += x[i] * y[i]
	}
	return s
}

// Norm2 returns the Euclidean (l2) norm of x.
func Norm2(x []float64) float64 {
	return math.Sqrt(DotProduct(x, x))
}

// Norm1 returns the l1 norm of x.
func Norm1(x []float64) float64 {
	var s float64
	for _, v := range x {
		s += math.Abs(v)
	}
	return s
}

// NormInf returns the l-infinity (maximum absolute) norm of x.
func NormInf(x []float64) float64 {
	var m float64
	for _, v := range x {
		if a := math.Abs(v); a > m {
			m = a
		}
	}
	return m
}

// NormP returns the discrete l^p norm of x for p >= 1.
func NormP(x []float64, p float64) float64 {
	var s float64
	for _, v := range x {
		s += math.Pow(math.Abs(v), p)
	}
	return math.Pow(s, 1/p)
}

// MaxAbs returns the largest absolute value in x.
func MaxAbs(x []float64) float64 { return NormInf(x) }

// Sum returns the sum of the entries of x.
func Sum(x []float64) float64 {
	var s float64
	for _, v := range x {
		s += v
	}
	return s
}

// Mean returns the arithmetic mean of x.
func Mean(x []float64) float64 {
	if len(x) == 0 {
		return 0
	}
	return Sum(x) / float64(len(x))
}

// LinfError returns the l-infinity norm of the difference a-b.
func LinfError(a, b []float64) float64 {
	return NormInf(VectorSub(a, b))
}

// L2Error returns the Euclidean norm of the difference a-b.
func L2Error(a, b []float64) float64 {
	return Norm2(VectorSub(a, b))
}

// RMSError returns the root-mean-square difference of a and b.
func RMSError(a, b []float64) float64 {
	if len(a) == 0 {
		return 0
	}
	return Norm2(VectorSub(a, b)) / math.Sqrt(float64(len(a)))
}

// RelativeL2Error returns ||a-b||_2 / ||b||_2, guarding against division by
// zero.
func RelativeL2Error(a, b []float64) float64 {
	den := Norm2(b)
	if den == 0 {
		return Norm2(VectorSub(a, b))
	}
	return Norm2(VectorSub(a, b)) / den
}

// ApplyFunc evaluates f at each entry of x and returns the resulting vector.
func ApplyFunc(f func(float64) float64, x []float64) []float64 {
	out := make([]float64, len(x))
	for i, v := range x {
		out[i] = f(v)
	}
	return out
}

// Reverse returns x with its entries in reverse order.
func Reverse(x []float64) []float64 {
	n := len(x)
	out := make([]float64, n)
	for i := 0; i < n; i++ {
		out[i] = x[n-1-i]
	}
	return out
}
