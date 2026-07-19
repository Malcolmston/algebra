package timeseries

import (
	"math"
	"sort"
)

// Mean returns the arithmetic mean of the series, or NaN if it is empty.
func Mean(x []float64) float64 {
	return mean(x)
}

// Sum returns the sum of all observations in the series.
func Sum(x []float64) float64 {
	return sumf(x)
}

// Variance returns the sample (unbiased, divide-by-n−1) variance of the
// series. It returns NaN if fewer than two observations are present.
func Variance(x []float64) float64 {
	n := len(x)
	if n < 2 {
		return math.NaN()
	}
	m := mean(x)
	var s float64
	for _, v := range x {
		d := v - m
		s += d * d
	}
	return s / float64(n-1)
}

// PopVariance returns the population (divide-by-n) variance of the series. It
// returns NaN for an empty series.
func PopVariance(x []float64) float64 {
	n := len(x)
	if n == 0 {
		return math.NaN()
	}
	m := mean(x)
	var s float64
	for _, v := range x {
		d := v - m
		s += d * d
	}
	return s / float64(n)
}

// StdDev returns the sample standard deviation (square root of [Variance]).
func StdDev(x []float64) float64 {
	return math.Sqrt(Variance(x))
}

// PopStdDev returns the population standard deviation (square root of
// [PopVariance]).
func PopStdDev(x []float64) float64 {
	return math.Sqrt(PopVariance(x))
}

// Min returns the smallest value in the series, or NaN if it is empty.
func Min(x []float64) float64 {
	if len(x) == 0 {
		return math.NaN()
	}
	m := x[0]
	for _, v := range x[1:] {
		if v < m {
			m = v
		}
	}
	return m
}

// Max returns the largest value in the series, or NaN if it is empty.
func Max(x []float64) float64 {
	if len(x) == 0 {
		return math.NaN()
	}
	m := x[0]
	for _, v := range x[1:] {
		if v > m {
			m = v
		}
	}
	return m
}

// Range returns the difference between the maximum and minimum values.
func Range(x []float64) float64 {
	return Max(x) - Min(x)
}

// Argmin returns the index of the smallest value, or -1 if the series is
// empty. The first index is returned on ties.
func Argmin(x []float64) int {
	if len(x) == 0 {
		return -1
	}
	idx := 0
	for i, v := range x {
		if v < x[idx] {
			idx = i
		}
	}
	return idx
}

// Argmax returns the index of the largest value, or -1 if the series is
// empty. The first index is returned on ties.
func Argmax(x []float64) int {
	if len(x) == 0 {
		return -1
	}
	idx := 0
	for i, v := range x {
		if v > x[idx] {
			idx = i
		}
	}
	return idx
}

// Median returns the median of the series, or NaN if it is empty. The input is
// not modified.
func Median(x []float64) float64 {
	n := len(x)
	if n == 0 {
		return math.NaN()
	}
	s := copyf(x)
	sort.Float64s(s)
	if n%2 == 1 {
		return s[n/2]
	}
	return 0.5 * (s[n/2-1] + s[n/2])
}

// Quantile returns the q-quantile (0 ≤ q ≤ 1) of the series using linear
// interpolation between order statistics. It returns NaN for an empty series
// or an out-of-range q. The input is not modified.
func Quantile(x []float64, q float64) float64 {
	n := len(x)
	if n == 0 || q < 0 || q > 1 || math.IsNaN(q) {
		return math.NaN()
	}
	s := copyf(x)
	sort.Float64s(s)
	if n == 1 {
		return s[0]
	}
	pos := q * float64(n-1)
	lo := int(math.Floor(pos))
	hi := lo + 1
	if hi >= n {
		return s[n-1]
	}
	frac := pos - float64(lo)
	return s[lo]*(1-frac) + s[hi]*frac
}

// Percentile returns the p-th percentile (0 ≤ p ≤ 100) of the series.
func Percentile(x []float64, p float64) float64 {
	return Quantile(x, p/100)
}

// MeanAbsoluteDeviation returns the mean of the absolute deviations from the
// series mean.
func MeanAbsoluteDeviation(x []float64) float64 {
	if len(x) == 0 {
		return math.NaN()
	}
	m := mean(x)
	var s float64
	for _, v := range x {
		s += math.Abs(v - m)
	}
	return s / float64(len(x))
}

// RootMeanSquare returns sqrt(mean(x²)), the quadratic mean of the series.
func RootMeanSquare(x []float64) float64 {
	if len(x) == 0 {
		return math.NaN()
	}
	var s float64
	for _, v := range x {
		s += v * v
	}
	return math.Sqrt(s / float64(len(x)))
}

// Energy returns the sum of squared observations of the series.
func Energy(x []float64) float64 {
	var s float64
	for _, v := range x {
		s += v * v
	}
	return s
}

// Skewness returns the sample skewness (third standardized moment) using the
// population standard deviation in the denominator. It returns NaN if the
// series has fewer than two observations or zero variance.
func Skewness(x []float64) float64 {
	n := len(x)
	if n < 2 {
		return math.NaN()
	}
	m := mean(x)
	var m2, m3 float64
	for _, v := range x {
		d := v - m
		m2 += d * d
		m3 += d * d * d
	}
	m2 /= float64(n)
	m3 /= float64(n)
	if m2 == 0 {
		return math.NaN()
	}
	return m3 / math.Pow(m2, 1.5)
}

// Kurtosis returns the excess kurtosis (fourth standardized moment minus 3) of
// the series. It returns NaN if the series has fewer than two observations or
// zero variance.
func Kurtosis(x []float64) float64 {
	n := len(x)
	if n < 2 {
		return math.NaN()
	}
	m := mean(x)
	var m2, m4 float64
	for _, v := range x {
		d := v - m
		dd := d * d
		m2 += dd
		m4 += dd * dd
	}
	m2 /= float64(n)
	m4 /= float64(n)
	if m2 == 0 {
		return math.NaN()
	}
	return m4/(m2*m2) - 3
}

// CoefficientOfVariation returns the ratio of the sample standard deviation to
// the mean, a scale-free measure of dispersion.
func CoefficientOfVariation(x []float64) float64 {
	m := mean(x)
	if m == 0 {
		return math.NaN()
	}
	return StdDev(x) / m
}

// First returns the first observation of the series, or NaN if empty.
func First(x []float64) float64 {
	if len(x) == 0 {
		return math.NaN()
	}
	return x[0]
}

// Last returns the final observation of the series, or NaN if empty.
func Last(x []float64) float64 {
	if len(x) == 0 {
		return math.NaN()
	}
	return x[len(x)-1]
}

// Reverse returns a new series with the observations in reverse order.
func Reverse(x []float64) []float64 {
	n := len(x)
	out := make([]float64, n)
	for i, v := range x {
		out[n-1-i] = v
	}
	return out
}
