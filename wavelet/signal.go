package wavelet

import (
	"math"
	"sort"
)

// PeriodicExtend returns a copy of x extended by pad samples on each side using
// periodic (wrap-around) boundary handling, so the result has length
// len(x)+2*pad. It panics if x is empty and pad is positive.
func PeriodicExtend(x []float64, pad int) []float64 {
	n := len(x)
	if pad < 0 {
		pad = 0
	}
	out := make([]float64, n+2*pad)
	for i := range out {
		idx := ((i-pad)%n + n) % n
		out[i] = x[idx]
	}
	return out
}

// SymmetricExtend returns a copy of x extended by pad samples on each side using
// whole-sample symmetric (mirror) boundary handling, so the result has length
// len(x)+2*pad. The boundary samples are not repeated. It panics if x is empty
// and pad is positive.
func SymmetricExtend(x []float64, pad int) []float64 {
	n := len(x)
	if pad < 0 {
		pad = 0
	}
	out := make([]float64, n+2*pad)
	period := 2 * n
	if n == 1 {
		period = 1
	}
	for i := range out {
		k := ((i - pad) % period)
		if k < 0 {
			k += period
		}
		if k >= n {
			k = period - 1 - k
		}
		out[i] = x[k]
	}
	return out
}

// ZeroExtend returns a copy of x extended by pad zero samples on each side, so
// the result has length len(x)+2*pad.
func ZeroExtend(x []float64, pad int) []float64 {
	if pad < 0 {
		pad = 0
	}
	out := make([]float64, len(x)+2*pad)
	copy(out[pad:], x)
	return out
}

// Convolve returns the full linear convolution of x and h, a slice of length
// len(x)+len(h)-1. It returns an empty slice if either input is empty.
func Convolve(x, h []float64) []float64 {
	if len(x) == 0 || len(h) == 0 {
		return []float64{}
	}
	out := make([]float64, len(x)+len(h)-1)
	for i, xi := range x {
		for j, hj := range h {
			out[i+j] += xi * hj
		}
	}
	return out
}

// Downsample returns every factor-th sample of x starting at index 0. A factor
// of 2 keeps the even-indexed samples. It panics if factor is not positive.
func Downsample(x []float64, factor int) []float64 {
	if factor <= 0 {
		panic("wavelet: Downsample factor must be positive")
	}
	out := make([]float64, 0, (len(x)+factor-1)/factor)
	for i := 0; i < len(x); i += factor {
		out = append(out, x[i])
	}
	return out
}

// Upsample returns x with factor-1 zeros inserted after each sample, producing
// a slice of length len(x)*factor. It panics if factor is not positive.
func Upsample(x []float64, factor int) []float64 {
	if factor <= 0 {
		panic("wavelet: Upsample factor must be positive")
	}
	out := make([]float64, len(x)*factor)
	for i, v := range x {
		out[i*factor] = v
	}
	return out
}

// L2Norm returns the Euclidean (l2) norm of x, sqrt(sum x_i^2).
func L2Norm(x []float64) float64 {
	return math.Sqrt(Energy(x))
}

// Energy returns the sum of squares of the elements of x.
func Energy(x []float64) float64 {
	var s float64
	for _, v := range x {
		s += v * v
	}
	return s
}

// MaxAbs returns the largest absolute value among the elements of x, or 0 for
// an empty slice.
func MaxAbs(x []float64) float64 {
	var m float64
	for _, v := range x {
		if a := math.Abs(v); a > m {
			m = a
		}
	}
	return m
}

// Median returns the median of x. For an even number of elements it returns the
// average of the two central values. It does not modify x and returns NaN for
// an empty slice.
func Median(x []float64) float64 {
	n := len(x)
	if n == 0 {
		return math.NaN()
	}
	cp := append([]float64(nil), x...)
	sort.Float64s(cp)
	if n%2 == 1 {
		return cp[n/2]
	}
	return 0.5 * (cp[n/2-1] + cp[n/2])
}

// MeanAbsoluteDeviation returns the median absolute deviation of x about its
// own median, median(|x_i - median(x)|). It is a robust measure of spread and
// returns NaN for an empty slice.
func MeanAbsoluteDeviation(x []float64) float64 {
	if len(x) == 0 {
		return math.NaN()
	}
	med := Median(x)
	dev := make([]float64, len(x))
	for i, v := range x {
		dev[i] = math.Abs(v - med)
	}
	return Median(dev)
}

// IsPowerOfTwo reports whether n is a positive power of two.
func IsPowerOfTwo(n int) bool {
	return n > 0 && n&(n-1) == 0
}

// NextPowerOfTwo returns the smallest power of two that is greater than or equal
// to n. It returns 1 for n less than or equal to 1.
func NextPowerOfTwo(n int) int {
	if n <= 1 {
		return 1
	}
	p := 1
	for p < n {
		p <<= 1
	}
	return p
}

// ZeroPadToPowerOfTwo returns a copy of x zero-padded on the right up to the
// next power-of-two length. If len(x) is already a power of two the input is
// returned unchanged (as a copy).
func ZeroPadToPowerOfTwo(x []float64) []float64 {
	target := NextPowerOfTwo(len(x))
	out := make([]float64, target)
	copy(out, x)
	return out
}

// wavelet2Valuation returns the number of times n can be halved while remaining
// even (the 2-adic valuation), which is the maximum number of dyadic
// decomposition levels feasible for a signal of length n.
func wavelet2Valuation(n int) int {
	if n <= 0 {
		return 0
	}
	v := 0
	for n%2 == 0 {
		v++
		n /= 2
	}
	return v
}
