package timeseries

import (
	"errors"
	"math"
)

// ADFResult holds the outcome of an augmented Dickey–Fuller test.
type ADFResult struct {
	Statistic float64 // ADF t-statistic on the lagged-level coefficient
	Gamma     float64 // estimated coefficient on y_{t-1}
	Lags      int     // number of augmenting difference lags used
	NObs      int     // number of observations in the regression
	Trend     string  // deterministic terms: "c" (constant) or "ct" (constant+trend)
}

// ADFTest performs the augmented Dickey–Fuller test for a unit root in the
// series. The regression is Δy_t = α + (β·t if trend) + γ·y_{t−1} + Σ δ_i·Δy_{t−i}
// + ε_t with the requested number of augmenting lags; the returned statistic is
// γ̂/se(γ̂). A more negative statistic is stronger evidence against a unit root
// (against non-stationarity). Set trend true to include a linear time trend. It
// returns an error if the series is too short.
func ADFTest(x []float64, lags int, trend bool) (*ADFResult, error) {
	if lags < 0 {
		return nil, errors.New("timeseries: ADF lags must be >= 0")
	}
	n := len(x)
	// Rows require y_{t-1}, lags of Δy, and Δy_t. First usable t index is
	// lags+1.
	start := lags + 1
	rows := n - start
	nParams := 2 + lags // intercept, gamma, plus lag diffs
	if trend {
		nParams++
	}
	if rows < nParams+1 {
		return nil, errors.New("timeseries: series too short for ADF with these lags")
	}
	dy := Diff(x) // dy[i] = x[i+1]-x[i]
	X := make([][]float64, rows)
	y := make([]float64, rows)
	gammaIdx := 1
	for t := start; t < n; t++ {
		r := t - start
		col := make([]float64, nParams)
		col[0] = 1 // intercept
		col[1] = x[t-1]
		idx := 2
		if trend {
			col[idx] = float64(t)
			idx++
		}
		for i := 1; i <= lags; i++ {
			// Δy_{t-i} = dy[t-i-1]
			col[idx] = dy[t-i-1]
			idx++
		}
		X[r] = col
		y[r] = dy[t-1] // Δy_t = x[t]-x[t-1] = dy[t-1]
	}
	beta, se, _, ok := olsStats(X, y)
	if !ok {
		return nil, errors.New("timeseries: ADF regression failed")
	}
	gamma := beta[gammaIdx]
	stat := math.NaN()
	if se[gammaIdx] > 0 {
		stat = gamma / se[gammaIdx]
	}
	tr := "c"
	if trend {
		tr = "ct"
	}
	return &ADFResult{
		Statistic: stat,
		Gamma:     gamma,
		Lags:      lags,
		NObs:      rows,
		Trend:     tr,
	}, nil
}

// DickeyFuller performs the simple (non-augmented) Dickey–Fuller test with a
// constant, equivalent to ADFTest with zero augmenting lags and no trend.
func DickeyFuller(x []float64) (*ADFResult, error) {
	return ADFTest(x, 0, false)
}

// ADFCriticalValue returns an approximate MacKinnon asymptotic critical value
// for the ADF t-statistic at the given significance level (0.01, 0.05 or 0.10)
// for the constant-only ("c") or constant-and-trend ("ct") specification.
// Values come from MacKinnon's response-surface constants (the asymptotic term).
func ADFCriticalValue(level float64, trend string) float64 {
	// Asymptotic critical values.
	switch trend {
	case "ct":
		switch {
		case approxEqual(level, 0.01, 1e-9):
			return -3.96
		case approxEqual(level, 0.05, 1e-9):
			return -3.41
		case approxEqual(level, 0.10, 1e-9):
			return -3.12
		}
	default: // "c"
		switch {
		case approxEqual(level, 0.01, 1e-9):
			return -3.43
		case approxEqual(level, 0.05, 1e-9):
			return -2.86
		case approxEqual(level, 0.10, 1e-9):
			return -2.57
		}
	}
	return math.NaN()
}

// IsStationaryADF reports whether the augmented Dickey–Fuller test rejects the
// unit-root null at the given significance level (typically 0.05), i.e. whether
// the series appears stationary. It uses no deterministic trend.
func IsStationaryADF(x []float64, lags int, level float64) bool {
	res, err := ADFTest(x, lags, false)
	if err != nil {
		return false
	}
	crit := ADFCriticalValue(level, "c")
	if math.IsNaN(crit) || math.IsNaN(res.Statistic) {
		return false
	}
	return res.Statistic < crit
}

// NumberOfDifferences estimates the order of differencing needed to make the
// series stationary by repeatedly applying the ADF test (constant, given lags)
// and differencing until stationarity is achieved or maxD is reached.
func NumberOfDifferences(x []float64, lags, maxD int, level float64) int {
	cur := copyf(x)
	for d := 0; d < maxD; d++ {
		if len(cur) < 3*(lags+2) {
			return d
		}
		if IsStationaryADF(cur, lags, level) {
			return d
		}
		cur = Diff(cur)
	}
	return maxD
}

// KPSSResult holds the outcome of a KPSS stationarity test.
type KPSSResult struct {
	Statistic float64 // KPSS LM statistic
	Lags      int     // truncation lag for the long-run variance
	Trend     string  // "c" for level stationarity, "ct" for trend stationarity
}

// KPSSTest performs the Kwiatkowski–Phillips–Schmidt–Shin test whose null
// hypothesis is that the series is (level- or trend-) stationary. Residuals are
// taken from a regression on a constant (trend=false) or a constant and linear
// trend (trend=true); the statistic is Σ S_t² / (n²·σ²_LR) with a Newey–West
// long-run variance estimate using the given truncation lag. A large statistic
// is evidence against stationarity. It returns an error if the series is too
// short.
func KPSSTest(x []float64, lag int, trend bool) (*KPSSResult, error) {
	n := len(x)
	if n < 4 || lag < 0 {
		return nil, errors.New("timeseries: series too short or invalid lag for KPSS")
	}
	// Regression residuals.
	var resid []float64
	if trend {
		X := make([][]float64, n)
		y := make([]float64, n)
		for i := 0; i < n; i++ {
			X[i] = []float64{1, float64(i)}
			y[i] = x[i]
		}
		beta, ok := leastSquares(X, y)
		if !ok {
			return nil, errors.New("timeseries: KPSS regression failed")
		}
		resid = make([]float64, n)
		for i := 0; i < n; i++ {
			resid[i] = x[i] - (beta[0] + beta[1]*float64(i))
		}
	} else {
		m := mean(x)
		resid = make([]float64, n)
		for i := 0; i < n; i++ {
			resid[i] = x[i] - m
		}
	}
	// Partial sums.
	S := make([]float64, n)
	var c float64
	for i := 0; i < n; i++ {
		c += resid[i]
		S[i] = c
	}
	var sumS2 float64
	for _, v := range S {
		sumS2 += v * v
	}
	// Newey–West long-run variance.
	var s0 float64
	for _, v := range resid {
		s0 += v * v
	}
	s0 /= float64(n)
	lr := s0
	for l := 1; l <= lag; l++ {
		var g float64
		for t := l; t < n; t++ {
			g += resid[t] * resid[t-l]
		}
		g /= float64(n)
		w := 1 - float64(l)/float64(lag+1)
		lr += 2 * w * g
	}
	if lr <= 0 {
		return nil, errors.New("timeseries: non-positive long-run variance in KPSS")
	}
	stat := sumS2 / (float64(n) * float64(n) * lr)
	tr := "c"
	if trend {
		tr = "ct"
	}
	return &KPSSResult{Statistic: stat, Lags: lag, Trend: tr}, nil
}

// VarianceRatio returns the Lo–MacKinlay variance ratio of the series at
// horizon q: the variance of q-period differences divided by q times the
// variance of one-period differences. A value near 1 is consistent with a
// random walk; values away from 1 indicate mean reversion (<1) or trending
// (>1). It returns NaN for q < 2 or a too-short series.
func VarianceRatio(x []float64, q int) float64 {
	n := len(x)
	if q < 2 || n < q+1 {
		return math.NaN()
	}
	d1 := Diff(x)
	var1 := PopVariance(d1)
	if var1 == 0 {
		return math.NaN()
	}
	// q-period differences.
	dq := make([]float64, n-q)
	for i := 0; i+q < n; i++ {
		dq[i] = x[i+q] - x[i]
	}
	varq := PopVariance(dq)
	return varq / (float64(q) * var1)
}
