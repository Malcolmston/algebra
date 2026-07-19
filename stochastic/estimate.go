package stochastic

import "math"

// SampleMean returns the arithmetic mean of x, or NaN for an empty slice.
func SampleMean(x []float64) float64 {
	if len(x) == 0 {
		return math.NaN()
	}
	s := 0.0
	for _, v := range x {
		s += v
	}
	return s / float64(len(x))
}

// SampleVariance returns the unbiased (Bessel-corrected) sample variance of x.
func SampleVariance(x []float64) float64 {
	n := len(x)
	if n < 2 {
		return 0
	}
	m := SampleMean(x)
	s := 0.0
	for _, v := range x {
		d := v - m
		s += d * d
	}
	return s / float64(n-1)
}

// SampleStdDev returns the unbiased sample standard deviation of x.
func SampleStdDev(x []float64) float64 { return math.Sqrt(SampleVariance(x)) }

// SampleCovariance returns the unbiased sample covariance of x and y, which must
// have the same length.
func SampleCovariance(x, y []float64) float64 {
	n := len(x)
	if n != len(y) || n < 2 {
		return 0
	}
	mx := SampleMean(x)
	my := SampleMean(y)
	s := 0.0
	for i := 0; i < n; i++ {
		s += (x[i] - mx) * (y[i] - my)
	}
	return s / float64(n-1)
}

// StandardError returns the standard error of the mean of x.
func StandardError(x []float64) float64 {
	if len(x) == 0 {
		return math.NaN()
	}
	return SampleStdDev(x) / math.Sqrt(float64(len(x)))
}

// ConfidenceInterval95 returns the approximate 95% confidence interval for the
// mean of x using the normal quantile 1.959964.
func ConfidenceInterval95(x []float64) (lo, hi float64) {
	m := SampleMean(x)
	se := StandardError(x)
	h := 1.959963984540054 * se
	return m - h, m + h
}

// Autocovariance returns the sample autocovariance of x at the given lag.
func Autocovariance(x []float64, lag int) float64 {
	n := len(x)
	if lag < 0 {
		lag = -lag
	}
	if lag >= n {
		return 0
	}
	m := SampleMean(x)
	s := 0.0
	for i := 0; i+lag < n; i++ {
		s += (x[i] - m) * (x[i+lag] - m)
	}
	return s / float64(n)
}

// Autocorrelation returns the sample autocorrelation of x at the given lag,
// normalized by the lag-0 autocovariance.
func Autocorrelation(x []float64, lag int) float64 {
	c0 := Autocovariance(x, 0)
	if c0 == 0 {
		return 0
	}
	return Autocovariance(x, lag) / c0
}

// LogReturns returns the log-returns log(v[i]/v[i-1]) of the path values. All
// values must be positive.
func LogReturns(p Path) []float64 {
	if len(p.Values) < 2 {
		return nil
	}
	out := make([]float64, len(p.Values)-1)
	for i := 1; i < len(p.Values); i++ {
		out[i-1] = math.Log(p.Values[i] / p.Values[i-1])
	}
	return out
}

// RealizedVariance returns the sum of squared log-returns of a positive-valued
// path, an estimator of the integrated variance.
func RealizedVariance(p Path) float64 {
	r := LogReturns(p)
	s := 0.0
	for _, v := range r {
		s += v * v
	}
	return s
}

// RealizedVolatility returns the square root of the realized variance,
// annualized by dividing the total horizon out. It equals sqrt(RealizedVariance
// / T) where T is the path horizon.
func RealizedVolatility(p Path) float64 {
	T := p.EndTime() - p.StartTime()
	if T <= 0 {
		return 0
	}
	return math.Sqrt(RealizedVariance(p) / T)
}

// EstimateGBMParams estimates the drift mu and volatility sigma of geometric
// Brownian motion from a positive-valued sample path with equally spaced times.
// The estimator uses the mean and variance of the log-returns.
func EstimateGBMParams(p Path) (mu, sigma float64) {
	r := LogReturns(p)
	if len(r) < 2 {
		return math.NaN(), math.NaN()
	}
	dt := (p.EndTime() - p.StartTime()) / float64(len(r))
	if dt <= 0 {
		return math.NaN(), math.NaN()
	}
	m := SampleMean(r)
	v := SampleVariance(r)
	sigma = math.Sqrt(v / dt)
	mu = m/dt + 0.5*sigma*sigma
	return mu, sigma
}

// EstimateDrift estimates the constant drift of a path with equally spaced times
// as the mean increment divided by the time step.
func EstimateDrift(p Path) float64 {
	inc := p.Increments()
	if len(inc) == 0 {
		return math.NaN()
	}
	dt := (p.EndTime() - p.StartTime()) / float64(len(inc))
	if dt <= 0 {
		return math.NaN()
	}
	return SampleMean(inc) / dt
}

// EstimateVolatility estimates the diffusion volatility of a path with equally
// spaced times as the standard deviation of increments divided by sqrt(dt).
func EstimateVolatility(p Path) float64 {
	inc := p.Increments()
	if len(inc) < 2 {
		return math.NaN()
	}
	dt := (p.EndTime() - p.StartTime()) / float64(len(inc))
	if dt <= 0 {
		return math.NaN()
	}
	return SampleStdDev(inc) / math.Sqrt(dt)
}

// EstimateOUParams estimates the mean-reversion speed theta, long-run mean mu
// and volatility sigma of an Ornstein-Uhlenbeck process from an equally spaced
// sample path, using the least-squares (AR(1)) estimator on X(t+dt) = a + b*X(t).
func EstimateOUParams(p Path) (theta, mu, sigma float64) {
	n := len(p.Values)
	if n < 3 {
		return math.NaN(), math.NaN(), math.NaN()
	}
	x := p.Values[:n-1]
	y := p.Values[1:]
	mx := SampleMean(x)
	my := SampleMean(y)
	var sxy, sxx float64
	for i := range x {
		sxy += (x[i] - mx) * (y[i] - my)
		sxx += (x[i] - mx) * (x[i] - mx)
	}
	if sxx == 0 {
		return math.NaN(), math.NaN(), math.NaN()
	}
	b := sxy / sxx
	a := my - b*mx
	dt := (p.EndTime() - p.StartTime()) / float64(n-1)
	if b <= 0 || b >= 1 || dt <= 0 {
		return math.NaN(), math.NaN(), math.NaN()
	}
	theta = -math.Log(b) / dt
	mu = a / (1 - b)
	// residual variance
	var sse float64
	for i := range x {
		e := y[i] - (a + b*x[i])
		sse += e * e
	}
	resVar := sse / float64(len(x)-2)
	sigma = math.Sqrt(2 * theta * resVar / (1 - b*b))
	return theta, mu, sigma
}

// EnsembleMean returns the pointwise mean across a set of paths sampled on the
// same grid. Paths shorter than the longest are ignored beyond their length.
func EnsembleMean(paths []Path) []float64 {
	if len(paths) == 0 {
		return nil
	}
	n := len(paths[0].Values)
	sum := make([]float64, n)
	cnt := make([]int, n)
	for _, p := range paths {
		for i := 0; i < n && i < len(p.Values); i++ {
			sum[i] += p.Values[i]
			cnt[i]++
		}
	}
	for i := range sum {
		if cnt[i] > 0 {
			sum[i] /= float64(cnt[i])
		}
	}
	return sum
}

// EnsembleVariance returns the pointwise unbiased variance across a set of
// paths sampled on the same grid.
func EnsembleVariance(paths []Path) []float64 {
	if len(paths) == 0 {
		return nil
	}
	n := len(paths[0].Values)
	m := EnsembleMean(paths)
	ss := make([]float64, n)
	cnt := make([]int, n)
	for _, p := range paths {
		for i := 0; i < n && i < len(p.Values); i++ {
			d := p.Values[i] - m[i]
			ss[i] += d * d
			cnt[i]++
		}
	}
	for i := range ss {
		if cnt[i] > 1 {
			ss[i] /= float64(cnt[i] - 1)
		} else {
			ss[i] = 0
		}
	}
	return ss
}

// EnsembleFinalValues returns the terminal value of each path.
func EnsembleFinalValues(paths []Path) []float64 {
	out := make([]float64, len(paths))
	for i, p := range paths {
		out[i] = p.Final()
	}
	return out
}
