package timeseries

import "math"

// AutoCovarianceAt returns the biased sample autocovariance of the series at
// lag k: (1/n)·Σ (x_t − x̄)(x_{t+k} − x̄). It returns NaN if k is negative or ≥
// len(x).
func AutoCovarianceAt(x []float64, k int) float64 {
	n := len(x)
	if k < 0 || k >= n {
		return math.NaN()
	}
	m := mean(x)
	var s float64
	for t := 0; t+k < n; t++ {
		s += (x[t] - m) * (x[t+k] - m)
	}
	return s / float64(n)
}

// AutoCovariance returns the biased sample autocovariances of the series for
// lags 0,1,…,maxlag. Element 0 is the population variance.
func AutoCovariance(x []float64, maxlag int) []float64 {
	if maxlag < 0 {
		return []float64{}
	}
	n := len(x)
	m := mean(x)
	out := make([]float64, maxlag+1)
	for k := 0; k <= maxlag; k++ {
		if k >= n {
			out[k] = 0
			continue
		}
		var s float64
		for t := 0; t+k < n; t++ {
			s += (x[t] - m) * (x[t+k] - m)
		}
		out[k] = s / float64(n)
	}
	return out
}

// AutoCorrelationAt returns the sample autocorrelation at lag k, i.e.
// [AutoCovarianceAt](x,k) divided by the lag-0 autocovariance.
func AutoCorrelationAt(x []float64, k int) float64 {
	g0 := AutoCovarianceAt(x, 0)
	if g0 == 0 {
		return math.NaN()
	}
	return AutoCovarianceAt(x, k) / g0
}

// AutoCorrelation returns the sample autocorrelation function (ACF) for lags
// 0,1,…,maxlag. Element 0 is always 1 (unless the series has zero variance, in
// which case NaN values are returned).
func AutoCorrelation(x []float64, maxlag int) []float64 {
	g := AutoCovariance(x, maxlag)
	if len(g) == 0 {
		return g
	}
	g0 := g[0]
	out := make([]float64, len(g))
	if g0 == 0 {
		for i := range out {
			out[i] = math.NaN()
		}
		return out
	}
	for i, v := range g {
		out[i] = v / g0
	}
	return out
}

// PartialAutoCorrelation returns the sample partial autocorrelation function
// (PACF) for lags 0,1,…,maxlag computed with the Durbin–Levinson recursion.
// Element 0 is 1 by convention and element k is the last coefficient of the
// fitted AR(k) model.
func PartialAutoCorrelation(x []float64, maxlag int) []float64 {
	rho := AutoCorrelation(x, maxlag)
	out := make([]float64, maxlag+1)
	if maxlag < 0 {
		return []float64{}
	}
	out[0] = 1
	if maxlag == 0 {
		return out
	}
	phi := make([]float64, maxlag+1)
	prev := make([]float64, maxlag+1)
	phi[1] = rho[1]
	out[1] = rho[1]
	for k := 2; k <= maxlag; k++ {
		copy(prev, phi)
		num := rho[k]
		den := 1.0
		for j := 1; j < k; j++ {
			num -= prev[j] * rho[k-j]
			den -= prev[j] * rho[j]
		}
		var pkk float64
		if den != 0 {
			pkk = num / den
		}
		phi[k] = pkk
		for j := 1; j < k; j++ {
			phi[j] = prev[j] - pkk*prev[k-j]
		}
		out[k] = pkk
	}
	return out
}

// CrossCovarianceAt returns the sample cross-covariance between x and y at lag
// k: (1/n)·Σ (x_t − x̄)(y_{t+k} − ȳ) for k ≥ 0, and the symmetric definition for
// k < 0. The two series must have equal length; otherwise NaN is returned.
func CrossCovarianceAt(x, y []float64, k int) float64 {
	n := len(x)
	if n == 0 || len(y) != n || k <= -n || k >= n {
		return math.NaN()
	}
	mx := mean(x)
	my := mean(y)
	var s float64
	if k >= 0 {
		for t := 0; t+k < n; t++ {
			s += (x[t] - mx) * (y[t+k] - my)
		}
	} else {
		for t := -k; t < n; t++ {
			s += (x[t] - mx) * (y[t+k] - my)
		}
	}
	return s / float64(n)
}

// CrossCorrelationAt returns the sample cross-correlation between x and y at
// lag k, normalized by the product of the two series' standard deviations.
func CrossCorrelationAt(x, y []float64, k int) float64 {
	cov := CrossCovarianceAt(x, y, k)
	sx := PopStdDev(x)
	sy := PopStdDev(y)
	if sx == 0 || sy == 0 {
		return math.NaN()
	}
	return cov / (sx * sy)
}

// CrossCorrelation returns the cross-correlation function of x and y for lags
// −maxlag,…,0,…,maxlag as a slice of length 2·maxlag+1, indexed so that
// out[maxlag+k] is the correlation at lag k.
func CrossCorrelation(x, y []float64, maxlag int) []float64 {
	if maxlag < 0 {
		return []float64{}
	}
	out := make([]float64, 2*maxlag+1)
	for k := -maxlag; k <= maxlag; k++ {
		out[maxlag+k] = CrossCorrelationAt(x, y, k)
	}
	return out
}

// ACFConfidenceBound returns the approximate two-sided 95% confidence bound
// 1.96/√n for the autocorrelations of white noise of length n.
func ACFConfidenceBound(n int) float64 {
	if n <= 0 {
		return math.NaN()
	}
	return 1.96 / math.Sqrt(float64(n))
}

// LjungBox returns the Ljung–Box Q statistic for the first h autocorrelation
// lags, a test for overall autocorrelation. Under the null of no
// autocorrelation Q is approximately chi-squared with h degrees of freedom.
func LjungBox(x []float64, h int) float64 {
	n := len(x)
	if n < 2 || h < 1 {
		return math.NaN()
	}
	rho := AutoCorrelation(x, h)
	nf := float64(n)
	var q float64
	for k := 1; k <= h; k++ {
		q += rho[k] * rho[k] / (nf - float64(k))
	}
	return nf * (nf + 2) * q
}

// BoxPierce returns the Box–Pierce Q statistic for the first h autocorrelation
// lags, the simpler large-sample precursor to [LjungBox].
func BoxPierce(x []float64, h int) float64 {
	n := len(x)
	if n < 2 || h < 1 {
		return math.NaN()
	}
	rho := AutoCorrelation(x, h)
	var q float64
	for k := 1; k <= h; k++ {
		q += rho[k] * rho[k]
	}
	return float64(n) * q
}

// DurbinWatson returns the Durbin–Watson statistic of a residual series,
// Σ(e_t − e_{t−1})² / Σe_t², which detects lag-1 autocorrelation. Values near 2
// indicate no autocorrelation.
func DurbinWatson(e []float64) float64 {
	if len(e) < 2 {
		return math.NaN()
	}
	var num, den float64
	for i, v := range e {
		den += v * v
		if i > 0 {
			d := v - e[i-1]
			num += d * d
		}
	}
	if den == 0 {
		return math.NaN()
	}
	return num / den
}
