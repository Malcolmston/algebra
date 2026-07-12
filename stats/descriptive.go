package stats

import (
	"math"
	"sort"
)

// Sum returns the sum of xs. The sum of an empty slice is 0.
func Sum(xs []float64) float64 {
	s := 0.0
	for _, x := range xs {
		s += x
	}
	return s
}

// Product returns the product of xs. The product of an empty slice is 1.
func Product(xs []float64) float64 {
	p := 1.0
	for _, x := range xs {
		p *= x
	}
	return p
}

// Min returns the smallest value in xs, or NaN if xs is empty.
func Min(xs []float64) float64 {
	if len(xs) == 0 {
		return math.NaN()
	}
	m := xs[0]
	for _, x := range xs[1:] {
		if x < m {
			m = x
		}
	}
	return m
}

// Max returns the largest value in xs, or NaN if xs is empty.
func Max(xs []float64) float64 {
	if len(xs) == 0 {
		return math.NaN()
	}
	m := xs[0]
	for _, x := range xs[1:] {
		if x > m {
			m = x
		}
	}
	return m
}

// Range returns the difference between the largest and smallest values in xs,
// or NaN if xs is empty.
func Range(xs []float64) float64 {
	if len(xs) == 0 {
		return math.NaN()
	}
	return Max(xs) - Min(xs)
}

// Mean returns the arithmetic mean of xs, or NaN if xs is empty.
func Mean(xs []float64) float64 {
	if len(xs) == 0 {
		return math.NaN()
	}
	return Sum(xs) / float64(len(xs))
}

// Median returns the median of xs, or NaN if xs is empty. For an even number
// of elements it returns the average of the two middle values. The input
// slice is not modified.
func Median(xs []float64) float64 {
	n := len(xs)
	if n == 0 {
		return math.NaN()
	}
	s := append([]float64(nil), xs...)
	sort.Float64s(s)
	if n%2 == 1 {
		return s[n/2]
	}
	return (s[n/2-1] + s[n/2]) / 2
}

// Mode returns the value or values that occur most frequently in xs, sorted
// in ascending order. If every value is unique (maximum frequency 1) there is
// no mode and Mode returns nil. Values are compared for exact equality, so
// Mode is most meaningful on discrete or rounded data.
func Mode(xs []float64) []float64 {
	if len(xs) == 0 {
		return nil
	}
	counts := make(map[float64]int, len(xs))
	best := 0
	for _, x := range xs {
		counts[x]++
		if counts[x] > best {
			best = counts[x]
		}
	}
	if best <= 1 {
		return nil
	}
	var modes []float64
	for v, c := range counts {
		if c == best {
			modes = append(modes, v)
		}
	}
	sort.Float64s(modes)
	return modes
}

// Variance returns the unbiased sample variance of xs using the n-1 (Bessel)
// denominator. It returns NaN if xs has fewer than two elements.
func Variance(xs []float64) float64 {
	n := len(xs)
	if n < 2 {
		return math.NaN()
	}
	m := Mean(xs)
	ss := 0.0
	for _, x := range xs {
		d := x - m
		ss += d * d
	}
	return ss / float64(n-1)
}

// PopVariance returns the population variance of xs using the n denominator.
// It returns NaN if xs is empty.
func PopVariance(xs []float64) float64 {
	n := len(xs)
	if n == 0 {
		return math.NaN()
	}
	m := Mean(xs)
	ss := 0.0
	for _, x := range xs {
		d := x - m
		ss += d * d
	}
	return ss / float64(n)
}

// StdDev returns the sample standard deviation of xs (the square root of
// [Variance]). It returns NaN if xs has fewer than two elements.
func StdDev(xs []float64) float64 {
	return math.Sqrt(Variance(xs))
}

// PopStdDev returns the population standard deviation of xs (the square root
// of [PopVariance]). It returns NaN if xs is empty.
func PopStdDev(xs []float64) float64 {
	return math.Sqrt(PopVariance(xs))
}

// Quantile returns the q-quantile of xs for q in [0, 1] using linear
// interpolation between the closest ranks (the "type 7" method used by NumPy
// and R by default). It returns NaN if xs is empty or q is outside [0, 1].
// The input slice is not modified.
func Quantile(xs []float64, q float64) float64 {
	n := len(xs)
	if n == 0 || q < 0 || q > 1 || math.IsNaN(q) {
		return math.NaN()
	}
	s := append([]float64(nil), xs...)
	sort.Float64s(s)
	if n == 1 {
		return s[0]
	}
	pos := q * float64(n-1)
	lo := int(math.Floor(pos))
	hi := int(math.Ceil(pos))
	if lo == hi {
		return s[lo]
	}
	frac := pos - float64(lo)
	return s[lo] + frac*(s[hi]-s[lo])
}

// Percentile returns the p-th percentile of xs for p in [0, 100]. It is a
// convenience wrapper around [Quantile] with q = p/100.
func Percentile(xs []float64, p float64) float64 {
	return Quantile(xs, p/100)
}

// IQR returns the interquartile range of xs, the difference between the 75th
// and 25th percentiles. It returns NaN if xs is empty.
func IQR(xs []float64) float64 {
	if len(xs) == 0 {
		return math.NaN()
	}
	return Quantile(xs, 0.75) - Quantile(xs, 0.25)
}

// Skewness returns the population coefficient of skewness of xs (the third
// standardized moment). It returns NaN if xs has fewer than two elements or
// zero variance.
func Skewness(xs []float64) float64 {
	n := len(xs)
	if n < 2 {
		return math.NaN()
	}
	m := Mean(xs)
	var m2, m3 float64
	for _, x := range xs {
		d := x - m
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

// Kurtosis returns the population excess kurtosis of xs (the fourth
// standardized moment minus 3, so a normal distribution has kurtosis 0). It
// returns NaN if xs has fewer than two elements or zero variance.
func Kurtosis(xs []float64) float64 {
	n := len(xs)
	if n < 2 {
		return math.NaN()
	}
	m := Mean(xs)
	var m2, m4 float64
	for _, x := range xs {
		d := x - m
		d2 := d * d
		m2 += d2
		m4 += d2 * d2
	}
	m2 /= float64(n)
	m4 /= float64(n)
	if m2 == 0 {
		return math.NaN()
	}
	return m4/(m2*m2) - 3
}

// GeometricMean returns the geometric mean of xs, the n-th root of the product
// of the values. It returns NaN if xs is empty or contains a non-positive
// value.
func GeometricMean(xs []float64) float64 {
	n := len(xs)
	if n == 0 {
		return math.NaN()
	}
	sumLog := 0.0
	for _, x := range xs {
		if x <= 0 {
			return math.NaN()
		}
		sumLog += math.Log(x)
	}
	return math.Exp(sumLog / float64(n))
}

// HarmonicMean returns the harmonic mean of xs, the reciprocal of the mean of
// the reciprocals. It returns NaN if xs is empty or contains a non-positive
// value.
func HarmonicMean(xs []float64) float64 {
	n := len(xs)
	if n == 0 {
		return math.NaN()
	}
	sumRecip := 0.0
	for _, x := range xs {
		if x <= 0 {
			return math.NaN()
		}
		sumRecip += 1 / x
	}
	return float64(n) / sumRecip
}

// Covariance returns the unbiased sample covariance between xs and ys using
// the n-1 denominator. It returns NaN if the slices differ in length or have
// fewer than two elements.
func Covariance(xs, ys []float64) float64 {
	n := len(xs)
	if n != len(ys) || n < 2 {
		return math.NaN()
	}
	mx, my := Mean(xs), Mean(ys)
	s := 0.0
	for i := 0; i < n; i++ {
		s += (xs[i] - mx) * (ys[i] - my)
	}
	return s / float64(n-1)
}

// Correlation returns the Pearson product-moment correlation coefficient
// between xs and ys, a value in [-1, 1]. It returns NaN if the slices differ
// in length, have fewer than two elements, or either has zero variance.
func Correlation(xs, ys []float64) float64 {
	n := len(xs)
	if n != len(ys) || n < 2 {
		return math.NaN()
	}
	mx, my := Mean(xs), Mean(ys)
	var sxy, sxx, syy float64
	for i := 0; i < n; i++ {
		dx := xs[i] - mx
		dy := ys[i] - my
		sxy += dx * dy
		sxx += dx * dx
		syy += dy * dy
	}
	den := math.Sqrt(sxx * syy)
	if den == 0 {
		return math.NaN()
	}
	return sxy / den
}

// ZScore returns the standard score of x given a distribution mean and
// standard deviation: (x - mean) / stdDev. It returns NaN if stdDev is zero.
func ZScore(x, mean, stdDev float64) float64 {
	if stdDev == 0 {
		return math.NaN()
	}
	return (x - mean) / stdDev
}

// WeightedMean returns the weighted arithmetic mean of xs with the given
// weights: sum(w*x) / sum(w). It returns NaN if the slices differ in length,
// are empty, or the weights sum to zero.
func WeightedMean(xs, weights []float64) float64 {
	n := len(xs)
	if n == 0 || n != len(weights) {
		return math.NaN()
	}
	var num, den float64
	for i := 0; i < n; i++ {
		num += weights[i] * xs[i]
		den += weights[i]
	}
	if den == 0 {
		return math.NaN()
	}
	return num / den
}
