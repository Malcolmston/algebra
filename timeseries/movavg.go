package timeseries

import (
	"math"
	"sort"
)

// MovingAverage returns the causal (trailing) simple moving average with window
// length w: out[i] is the mean of the up-to-w most recent samples ending at i.
// At the start, where fewer than w samples exist, the average uses the samples
// seen so far. The result has the same length as x. It returns nil if w < 1.
func MovingAverage(x []float64, w int) []float64 {
	if w < 1 {
		return nil
	}
	out := make([]float64, len(x))
	var sum float64
	for i := range x {
		sum += x[i]
		if i >= w {
			sum -= x[i-w]
		}
		denom := w
		if i+1 < w {
			denom = i + 1
		}
		out[i] = sum / float64(denom)
	}
	return out
}

// MovingAverageCentered returns the centered simple moving average with window
// length w, truncating the window at the boundaries so the output length equals
// the input length. For even w the window is biased one sample toward the past.
// It returns nil if w < 1.
func MovingAverageCentered(x []float64, w int) []float64 {
	if w < 1 {
		return nil
	}
	n := len(x)
	out := make([]float64, n)
	half := w / 2
	for i := range x {
		lo := i - half
		hi := i + (w - 1 - half)
		if lo < 0 {
			lo = 0
		}
		if hi >= n {
			hi = n - 1
		}
		var s float64
		for j := lo; j <= hi; j++ {
			s += x[j]
		}
		out[i] = s / float64(hi-lo+1)
	}
	return out
}

// MovingAverageValid returns the centered simple moving average over full
// windows only, so the result has length len(x)−w+1. It returns nil if w < 1 or
// w > len(x).
func MovingAverageValid(x []float64, w int) []float64 {
	n := len(x)
	if w < 1 || w > n {
		return nil
	}
	out := make([]float64, n-w+1)
	var sum float64
	for i := 0; i < w; i++ {
		sum += x[i]
	}
	out[0] = sum / float64(w)
	for i := w; i < n; i++ {
		sum += x[i] - x[i-w]
		out[i-w+1] = sum / float64(w)
	}
	return out
}

// WeightedMovingAverage returns the causal weighted moving average whose window
// length equals len(weights). out[i] is Σ weights[j]·x[i−len+1+j] divided by the
// sum of the weights used; at the boundary the trailing weights are applied to
// the available samples. It returns nil for empty weights.
func WeightedMovingAverage(x []float64, weights []float64) []float64 {
	w := len(weights)
	if w == 0 {
		return nil
	}
	n := len(x)
	out := make([]float64, n)
	for i := 0; i < n; i++ {
		var num, den float64
		for j := 0; j < w; j++ {
			idx := i - w + 1 + j
			if idx < 0 {
				continue
			}
			num += weights[j] * x[idx]
			den += weights[j]
		}
		if den == 0 {
			out[i] = math.NaN()
		} else {
			out[i] = num / den
		}
	}
	return out
}

// TriangularMovingAverage returns a trailing weighted moving average whose
// weights increase linearly from 1 up to w (giving most weight to the most
// recent sample), a smooth double-averaging filter. It returns nil if w < 1.
func TriangularMovingAverage(x []float64, w int) []float64 {
	if w < 1 {
		return nil
	}
	weights := make([]float64, w)
	for i := range weights {
		weights[i] = float64(i + 1)
	}
	return WeightedMovingAverage(x, weights)
}

// ExponentialMovingAverage returns the exponentially weighted moving average
// with smoothing factor alpha in (0,1]: out[0] = x[0] and
// out[i] = alpha·x[i] + (1−alpha)·out[i−1]. It returns nil for an out-of-range
// alpha.
func ExponentialMovingAverage(x []float64, alpha float64) []float64 {
	if alpha <= 0 || alpha > 1 {
		return nil
	}
	out := make([]float64, len(x))
	if len(x) == 0 {
		return out
	}
	out[0] = x[0]
	for i := 1; i < len(x); i++ {
		out[i] = alpha*x[i] + (1-alpha)*out[i-1]
	}
	return out
}

// ExponentialMovingAverageSpan returns the exponential moving average using the
// span convention alpha = 2/(span+1), common in technical analysis. It returns
// nil if span < 1.
func ExponentialMovingAverageSpan(x []float64, span int) []float64 {
	if span < 1 {
		return nil
	}
	return ExponentialMovingAverage(x, 2/float64(span+1))
}

// DoubleExponentialMovingAverage returns the DEMA indicator,
// 2·EMA(x) − EMA(EMA(x)), which reduces the lag of a plain EMA. It returns nil
// for an out-of-range alpha.
func DoubleExponentialMovingAverage(x []float64, alpha float64) []float64 {
	e1 := ExponentialMovingAverage(x, alpha)
	if e1 == nil {
		return nil
	}
	e2 := ExponentialMovingAverage(e1, alpha)
	out := make([]float64, len(x))
	for i := range x {
		out[i] = 2*e1[i] - e2[i]
	}
	return out
}

// TripleExponentialMovingAverage returns the TEMA indicator,
// 3·EMA − 3·EMA² + EMA³, which further reduces lag. It returns nil for an
// out-of-range alpha.
func TripleExponentialMovingAverage(x []float64, alpha float64) []float64 {
	e1 := ExponentialMovingAverage(x, alpha)
	if e1 == nil {
		return nil
	}
	e2 := ExponentialMovingAverage(e1, alpha)
	e3 := ExponentialMovingAverage(e2, alpha)
	out := make([]float64, len(x))
	for i := range x {
		out[i] = 3*e1[i] - 3*e2[i] + e3[i]
	}
	return out
}

// CumulativeMovingAverage returns the running (expanding-window) mean of the
// series: out[i] is the mean of x[0..i].
func CumulativeMovingAverage(x []float64) []float64 {
	out := make([]float64, len(x))
	var s float64
	for i, v := range x {
		s += v
		out[i] = s / float64(i+1)
	}
	return out
}

// RollingMean returns the trailing rolling mean over full windows of length w,
// producing len(x)−w+1 values. It returns nil if w < 1 or w > len(x).
func RollingMean(x []float64, w int) []float64 {
	return MovingAverageValid(x, w)
}

// RollingSum returns the trailing rolling sum over full windows of length w. It
// returns nil if w < 1 or w > len(x).
func RollingSum(x []float64, w int) []float64 {
	n := len(x)
	if w < 1 || w > n {
		return nil
	}
	out := make([]float64, n-w+1)
	var sum float64
	for i := 0; i < w; i++ {
		sum += x[i]
	}
	out[0] = sum
	for i := w; i < n; i++ {
		sum += x[i] - x[i-w]
		out[i-w+1] = sum
	}
	return out
}

// RollingVariance returns the trailing rolling sample variance over full
// windows of length w. It returns nil if w < 2 or w > len(x).
func RollingVariance(x []float64, w int) []float64 {
	n := len(x)
	if w < 2 || w > n {
		return nil
	}
	out := make([]float64, n-w+1)
	for i := 0; i+w <= n; i++ {
		out[i] = Variance(x[i : i+w])
	}
	return out
}

// RollingStdDev returns the trailing rolling sample standard deviation over
// full windows of length w. It returns nil if w < 2 or w > len(x).
func RollingStdDev(x []float64, w int) []float64 {
	v := RollingVariance(x, w)
	if v == nil {
		return nil
	}
	for i := range v {
		v[i] = math.Sqrt(v[i])
	}
	return v
}

// RollingMin returns the trailing rolling minimum over full windows of length
// w. It returns nil if w < 1 or w > len(x).
func RollingMin(x []float64, w int) []float64 {
	n := len(x)
	if w < 1 || w > n {
		return nil
	}
	out := make([]float64, n-w+1)
	for i := 0; i+w <= n; i++ {
		out[i] = Min(x[i : i+w])
	}
	return out
}

// RollingMax returns the trailing rolling maximum over full windows of length
// w. It returns nil if w < 1 or w > len(x).
func RollingMax(x []float64, w int) []float64 {
	n := len(x)
	if w < 1 || w > n {
		return nil
	}
	out := make([]float64, n-w+1)
	for i := 0; i+w <= n; i++ {
		out[i] = Max(x[i : i+w])
	}
	return out
}

// RollingMedian returns the trailing rolling median over full windows of length
// w. It returns nil if w < 1 or w > len(x).
func RollingMedian(x []float64, w int) []float64 {
	n := len(x)
	if w < 1 || w > n {
		return nil
	}
	out := make([]float64, n-w+1)
	buf := make([]float64, w)
	for i := 0; i+w <= n; i++ {
		copy(buf, x[i:i+w])
		sort.Float64s(buf)
		if w%2 == 1 {
			out[i] = buf[w/2]
		} else {
			out[i] = 0.5 * (buf[w/2-1] + buf[w/2])
		}
	}
	return out
}

// MedianFilter returns the centered running median with window length w,
// truncating the window at the boundaries so the output length equals the input
// length. It returns nil if w < 1.
func MedianFilter(x []float64, w int) []float64 {
	if w < 1 {
		return nil
	}
	n := len(x)
	out := make([]float64, n)
	half := w / 2
	for i := range x {
		lo := i - half
		hi := i + (w - 1 - half)
		if lo < 0 {
			lo = 0
		}
		if hi >= n {
			hi = n - 1
		}
		out[i] = Median(x[lo : hi+1])
	}
	return out
}

// ExpandingMean returns the expanding-window mean, identical to
// [CumulativeMovingAverage].
func ExpandingMean(x []float64) []float64 {
	return CumulativeMovingAverage(x)
}

// ExpandingSum returns the expanding-window sum, identical to [CumSum].
func ExpandingSum(x []float64) []float64 {
	return CumSum(x)
}

// ExpandingMax returns the running maximum: out[i] = max(x[0..i]).
func ExpandingMax(x []float64) []float64 {
	out := make([]float64, len(x))
	for i, v := range x {
		if i == 0 || v > out[i-1] {
			out[i] = v
		} else {
			out[i] = out[i-1]
		}
	}
	return out
}

// ExpandingMin returns the running minimum: out[i] = min(x[0..i]).
func ExpandingMin(x []float64) []float64 {
	out := make([]float64, len(x))
	for i, v := range x {
		if i == 0 || v < out[i-1] {
			out[i] = v
		} else {
			out[i] = out[i-1]
		}
	}
	return out
}
