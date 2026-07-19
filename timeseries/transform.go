package timeseries

import "math"

// LinearFit holds the coefficients of a straight line y = Intercept + Slope·t.
type LinearFit struct {
	Intercept float64
	Slope     float64
}

// At evaluates the fitted line at index t.
func (f LinearFit) At(t float64) float64 {
	return f.Intercept + f.Slope*t
}

// Diff returns the first difference of the series: out[i] = x[i+1] − x[i]. The
// result has length len(x)−1 (empty for a series shorter than two).
func Diff(x []float64) []float64 {
	if len(x) < 2 {
		return []float64{}
	}
	out := make([]float64, len(x)-1)
	for i := 1; i < len(x); i++ {
		out[i-1] = x[i] - x[i-1]
	}
	return out
}

// DiffOrder applies the first difference d times. Differencing d times shortens
// the series by d. A d of zero returns a copy of x.
func DiffOrder(x []float64, d int) []float64 {
	out := copyf(x)
	for i := 0; i < d; i++ {
		out = Diff(out)
	}
	return out
}

// SeasonalDiff returns the seasonal difference at lag s: out[i] = x[i+s] − x[i].
// The result has length len(x)−s. It returns an empty slice if s < 1 or s ≥
// len(x).
func SeasonalDiff(x []float64, s int) []float64 {
	if s < 1 || s >= len(x) {
		return []float64{}
	}
	out := make([]float64, len(x)-s)
	for i := s; i < len(x); i++ {
		out[i-s] = x[i] - x[i-s]
	}
	return out
}

// Integrate is the inverse of [Diff]: given the first-difference series d and
// an initial value x0, it reconstructs the original series by cumulative
// summation. The result has length len(d)+1.
func Integrate(d []float64, x0 float64) []float64 {
	out := make([]float64, len(d)+1)
	out[0] = x0
	for i, v := range d {
		out[i+1] = out[i] + v
	}
	return out
}

// SeasonalIntegrate is the inverse of [SeasonalDiff] at lag s: given the
// seasonal-difference series d and the s initial values seed (the first season
// of the original series), it reconstructs the original series of length
// len(d)+s. It returns nil if len(seed) != s.
func SeasonalIntegrate(d []float64, seed []float64, s int) []float64 {
	if s < 1 || len(seed) != s {
		return nil
	}
	out := make([]float64, len(d)+s)
	copy(out, seed)
	for i := 0; i < len(d); i++ {
		out[i+s] = out[i] + d[i]
	}
	return out
}

// CumSum returns the cumulative sum of the series: out[i] = x[0]+…+x[i].
func CumSum(x []float64) []float64 {
	out := make([]float64, len(x))
	var s float64
	for i, v := range x {
		s += v
		out[i] = s
	}
	return out
}

// CumProd returns the cumulative product of the series.
func CumProd(x []float64) []float64 {
	out := make([]float64, len(x))
	p := 1.0
	for i, v := range x {
		p *= v
		out[i] = p
	}
	return out
}

// Lag shifts the series forward by k steps, padding the first k positions with
// NaN. The result has the same length as x. A negative k delegates to [Lead].
func Lag(x []float64, k int) []float64 {
	if k < 0 {
		return Lead(x, -k)
	}
	out := make([]float64, len(x))
	for i := range out {
		if i < k {
			out[i] = math.NaN()
		} else {
			out[i] = x[i-k]
		}
	}
	return out
}

// Lead shifts the series backward by k steps, padding the last k positions with
// NaN. The result has the same length as x. A negative k delegates to [Lag].
func Lead(x []float64, k int) []float64 {
	if k < 0 {
		return Lag(x, -k)
	}
	out := make([]float64, len(x))
	n := len(x)
	for i := range out {
		if i+k < n {
			out[i] = x[i+k]
		} else {
			out[i] = math.NaN()
		}
	}
	return out
}

// Shift shifts the series forward by k steps, filling vacated positions with
// fill. A negative k shifts backward. The result has the same length as x.
func Shift(x []float64, k int, fill float64) []float64 {
	out := make([]float64, len(x))
	n := len(x)
	for i := range out {
		j := i - k
		if j >= 0 && j < n {
			out[i] = x[j]
		} else {
			out[i] = fill
		}
	}
	return out
}

// LogTransform returns the natural logarithm of each observation. Non-positive
// inputs map to NaN.
func LogTransform(x []float64) []float64 {
	out := make([]float64, len(x))
	for i, v := range x {
		if v <= 0 {
			out[i] = math.NaN()
		} else {
			out[i] = math.Log(v)
		}
	}
	return out
}

// ExpTransform returns the exponential of each observation (the inverse of
// [LogTransform]).
func ExpTransform(x []float64) []float64 {
	out := make([]float64, len(x))
	for i, v := range x {
		out[i] = math.Exp(v)
	}
	return out
}

// SqrtTransform returns the square root of each observation. Negative inputs
// map to NaN.
func SqrtTransform(x []float64) []float64 {
	out := make([]float64, len(x))
	for i, v := range x {
		out[i] = math.Sqrt(v)
	}
	return out
}

// BoxCox applies the Box–Cox power transform with parameter lambda. For lambda
// = 0 it returns the natural log; otherwise (x^lambda − 1)/lambda. Non-positive
// inputs map to NaN.
func BoxCox(x []float64, lambda float64) []float64 {
	out := make([]float64, len(x))
	for i, v := range x {
		if v <= 0 {
			out[i] = math.NaN()
			continue
		}
		if lambda == 0 {
			out[i] = math.Log(v)
		} else {
			out[i] = (math.Pow(v, lambda) - 1) / lambda
		}
	}
	return out
}

// InverseBoxCox inverts the Box–Cox transform with parameter lambda.
func InverseBoxCox(y []float64, lambda float64) []float64 {
	out := make([]float64, len(y))
	for i, v := range y {
		if lambda == 0 {
			out[i] = math.Exp(v)
		} else {
			out[i] = math.Pow(lambda*v+1, 1/lambda)
		}
	}
	return out
}

// Demean returns the series with its mean subtracted from every observation.
func Demean(x []float64) []float64 {
	m := mean(x)
	out := make([]float64, len(x))
	for i, v := range x {
		out[i] = v - m
	}
	return out
}

// Standardize returns the z-scored series (x − mean)/stddev using the sample
// standard deviation. If the standard deviation is zero the result is all
// zeros.
func Standardize(x []float64) []float64 {
	m := mean(x)
	sd := StdDev(x)
	out := make([]float64, len(x))
	if sd == 0 || math.IsNaN(sd) {
		return out
	}
	for i, v := range x {
		out[i] = (v - m) / sd
	}
	return out
}

// MinMaxNormalize rescales the series linearly to the unit interval [0,1]. If
// all values are equal the result is all zeros.
func MinMaxNormalize(x []float64) []float64 {
	lo := Min(x)
	hi := Max(x)
	out := make([]float64, len(x))
	if hi == lo {
		return out
	}
	for i, v := range x {
		out[i] = (v - lo) / (hi - lo)
	}
	return out
}

// Rescale linearly maps the series to the interval [a,b].
func Rescale(x []float64, a, b float64) []float64 {
	n := MinMaxNormalize(x)
	out := make([]float64, len(n))
	for i, v := range n {
		out[i] = a + v*(b-a)
	}
	return out
}

// Clip constrains every observation to the closed interval [lo,hi].
func Clip(x []float64, lo, hi float64) []float64 {
	out := make([]float64, len(x))
	for i, v := range x {
		switch {
		case v < lo:
			out[i] = lo
		case v > hi:
			out[i] = hi
		default:
			out[i] = v
		}
	}
	return out
}

// SimpleReturns returns the simple (arithmetic) returns of the series:
// out[i] = x[i+1]/x[i] − 1. The result has length len(x)−1.
func SimpleReturns(x []float64) []float64 {
	if len(x) < 2 {
		return []float64{}
	}
	out := make([]float64, len(x)-1)
	for i := 1; i < len(x); i++ {
		out[i-1] = x[i]/x[i-1] - 1
	}
	return out
}

// LogReturns returns the log returns of the series: out[i] = ln(x[i+1]/x[i]).
// The result has length len(x)−1.
func LogReturns(x []float64) []float64 {
	if len(x) < 2 {
		return []float64{}
	}
	out := make([]float64, len(x)-1)
	for i := 1; i < len(x); i++ {
		out[i-1] = math.Log(x[i] / x[i-1])
	}
	return out
}

// FracDiffWeights returns the first n binomial weights of the fractional
// differencing operator (1−B)^d. weight[0] is 1 and
// weight[k] = −weight[k−1]·(d−k+1)/k.
func FracDiffWeights(d float64, n int) []float64 {
	if n <= 0 {
		return []float64{}
	}
	w := make([]float64, n)
	w[0] = 1
	for k := 1; k < n; k++ {
		w[k] = -w[k-1] * (d - float64(k) + 1) / float64(k)
	}
	return w
}

// FractionalDifference applies the fractional differencing operator (1−B)^d to
// the series using an expanding window of binomial weights. The result has the
// same length as x; early observations use only the available past. For
// integer d this reproduces ordinary differencing (padded at the front).
func FractionalDifference(x []float64, d float64) []float64 {
	n := len(x)
	w := FracDiffWeights(d, n)
	out := make([]float64, n)
	for t := 0; t < n; t++ {
		var s float64
		for k := 0; k <= t; k++ {
			s += w[k] * x[t-k]
		}
		out[t] = s
	}
	return out
}

// FitLinearTrend fits a straight line y = a + b·t to the series against the
// index t = 0,1,…,n−1 by ordinary least squares and returns the coefficients.
func FitLinearTrend(x []float64) LinearFit {
	n := len(x)
	if n == 0 {
		return LinearFit{Intercept: math.NaN(), Slope: math.NaN()}
	}
	if n == 1 {
		return LinearFit{Intercept: x[0], Slope: 0}
	}
	nf := float64(n)
	var st, sy, stt, sty float64
	for i, v := range x {
		t := float64(i)
		st += t
		sy += v
		stt += t * t
		sty += t * v
	}
	denom := nf*stt - st*st
	if denom == 0 {
		return LinearFit{Intercept: mean(x), Slope: 0}
	}
	b := (nf*sty - st*sy) / denom
	a := (sy - b*st) / nf
	return LinearFit{Intercept: a, Slope: b}
}

// Detrend removes a least-squares linear trend from the series, returning the
// residuals x[i] − (a + b·i).
func Detrend(x []float64) []float64 {
	f := FitLinearTrend(x)
	out := make([]float64, len(x))
	for i, v := range x {
		out[i] = v - f.At(float64(i))
	}
	return out
}

// TrendLine returns the fitted trend values a + b·i for the series.
func TrendLine(x []float64) []float64 {
	f := FitLinearTrend(x)
	out := make([]float64, len(x))
	for i := range x {
		out[i] = f.At(float64(i))
	}
	return out
}
