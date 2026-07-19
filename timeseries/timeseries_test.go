package timeseries

import (
	"fmt"
	"math"
	"testing"
)

const tol = 1e-9

func approxSlice(got, want []float64, tol float64) bool {
	if len(got) != len(want) {
		return false
	}
	for i := range got {
		if math.IsNaN(want[i]) {
			if !math.IsNaN(got[i]) {
				return false
			}
			continue
		}
		if math.Abs(got[i]-want[i]) > tol {
			return false
		}
	}
	return true
}

// --- descriptive ---

func TestDescriptive(t *testing.T) {
	x := []float64{2, 4, 4, 4, 5, 5, 7, 9}
	if got := Mean(x); !approxEqual(got, 5, tol) {
		t.Errorf("Mean = %v, want 5", got)
	}
	if got := Sum(x); !approxEqual(got, 40, tol) {
		t.Errorf("Sum = %v, want 40", got)
	}
	if got := PopVariance(x); !approxEqual(got, 4, tol) {
		t.Errorf("PopVariance = %v, want 4", got)
	}
	if got := PopStdDev(x); !approxEqual(got, 2, tol) {
		t.Errorf("PopStdDev = %v, want 2", got)
	}
	if got := Variance(x); !approxEqual(got, 32.0/7, tol) {
		t.Errorf("Variance = %v, want %v", got, 32.0/7)
	}
	if got := Min(x); got != 2 {
		t.Errorf("Min = %v, want 2", got)
	}
	if got := Max(x); got != 9 {
		t.Errorf("Max = %v, want 9", got)
	}
	if got := Range(x); got != 7 {
		t.Errorf("Range = %v, want 7", got)
	}
	if got := Median([]float64{3, 1, 2, 4}); !approxEqual(got, 2.5, tol) {
		t.Errorf("Median = %v, want 2.5", got)
	}
	if got := Argmax(x); got != 7 {
		t.Errorf("Argmax = %v, want 7", got)
	}
	if got := Argmin(x); got != 0 {
		t.Errorf("Argmin = %v, want 0", got)
	}
}

func TestQuantile(t *testing.T) {
	x := []float64{1, 2, 3, 4}
	tests := []struct {
		q, want float64
	}{
		{0, 1}, {1, 4}, {0.5, 2.5}, {0.25, 1.75},
	}
	for _, tc := range tests {
		if got := Quantile(x, tc.q); !approxEqual(got, tc.want, tol) {
			t.Errorf("Quantile(%v) = %v, want %v", tc.q, got, tc.want)
		}
	}
}

func TestMomentsNaN(t *testing.T) {
	if !math.IsNaN(Mean(nil)) {
		t.Error("Mean(nil) should be NaN")
	}
	if !math.IsNaN(Variance([]float64{1})) {
		t.Error("Variance of single element should be NaN")
	}
}

func TestSkewKurtosis(t *testing.T) {
	// Symmetric data has zero skewness.
	x := []float64{1, 2, 3, 4, 5}
	if got := Skewness(x); math.Abs(got) > tol {
		t.Errorf("Skewness of symmetric = %v, want 0", got)
	}
	if k := Kurtosis(x); math.IsNaN(k) {
		t.Error("Kurtosis unexpectedly NaN")
	}
}

// --- transforms ---

func TestDiffIntegrate(t *testing.T) {
	x := []float64{1, 2, 4, 7, 11}
	d := Diff(x)
	if !approxSlice(d, []float64{1, 2, 3, 4}, tol) {
		t.Errorf("Diff = %v", d)
	}
	back := Integrate(d, x[0])
	if !approxSlice(back, x, tol) {
		t.Errorf("Integrate round-trip = %v, want %v", back, x)
	}
	if got := DiffOrder(x, 2); !approxSlice(got, []float64{1, 1, 1}, tol) {
		t.Errorf("DiffOrder 2 = %v", got)
	}
}

func TestSeasonalDiff(t *testing.T) {
	x := []float64{1, 2, 3, 5, 7, 9}
	d := SeasonalDiff(x, 3)
	if !approxSlice(d, []float64{4, 5, 6}, tol) {
		t.Errorf("SeasonalDiff = %v", d)
	}
	back := SeasonalIntegrate(d, x[:3], 3)
	if !approxSlice(back, x, tol) {
		t.Errorf("SeasonalIntegrate = %v, want %v", back, x)
	}
}

func TestCumSumProd(t *testing.T) {
	if got := CumSum([]float64{1, 2, 3, 4}); !approxSlice(got, []float64{1, 3, 6, 10}, tol) {
		t.Errorf("CumSum = %v", got)
	}
	if got := CumProd([]float64{1, 2, 3, 4}); !approxSlice(got, []float64{1, 2, 6, 24}, tol) {
		t.Errorf("CumProd = %v", got)
	}
}

func TestLagLead(t *testing.T) {
	x := []float64{1, 2, 3}
	if got := Lag(x, 1); !approxSlice(got, []float64{math.NaN(), 1, 2}, tol) {
		t.Errorf("Lag = %v", got)
	}
	if got := Lead(x, 1); !approxSlice(got, []float64{2, 3, math.NaN()}, tol) {
		t.Errorf("Lead = %v", got)
	}
	if got := Shift(x, 1, 0); !approxSlice(got, []float64{0, 1, 2}, tol) {
		t.Errorf("Shift = %v", got)
	}
}

func TestReturns(t *testing.T) {
	x := []float64{100, 110, 121}
	if got := SimpleReturns(x); !approxSlice(got, []float64{0.1, 0.1}, tol) {
		t.Errorf("SimpleReturns = %v", got)
	}
	lr := LogReturns(x)
	if !approxSlice(lr, []float64{math.Log(1.1), math.Log(1.1)}, tol) {
		t.Errorf("LogReturns = %v", lr)
	}
}

func TestBoxCox(t *testing.T) {
	x := []float64{1, 2, 4}
	y := BoxCox(x, 0.5)
	want := []float64{0, 2*math.Sqrt(2) - 2, 2}
	if !approxSlice(y, want, tol) {
		t.Errorf("BoxCox = %v, want %v", y, want)
	}
	back := InverseBoxCox(y, 0.5)
	if !approxSlice(back, x, 1e-9) {
		t.Errorf("InverseBoxCox = %v, want %v", back, x)
	}
	// lambda 0 is the log.
	if got := BoxCox([]float64{math.E}, 0); !approxEqual(got[0], 1, tol) {
		t.Errorf("BoxCox lambda 0 = %v, want 1", got[0])
	}
}

func TestFracDiffWeights(t *testing.T) {
	// (1-B)^1 weights: 1, -1, 0, 0, ...
	w := FracDiffWeights(1, 4)
	if !approxSlice(w, []float64{1, -1, 0, 0}, tol) {
		t.Errorf("FracDiffWeights(1) = %v", w)
	}
	// (1-B)^2 weights: 1, -2, 1, 0
	w2 := FracDiffWeights(2, 4)
	if !approxSlice(w2, []float64{1, -2, 1, 0}, tol) {
		t.Errorf("FracDiffWeights(2) = %v", w2)
	}
}

func TestFitLinearTrend(t *testing.T) {
	// y = 3 + 2t exactly.
	x := []float64{3, 5, 7, 9, 11}
	f := FitLinearTrend(x)
	if !approxEqual(f.Intercept, 3, tol) || !approxEqual(f.Slope, 2, tol) {
		t.Errorf("FitLinearTrend = %+v, want {3,2}", f)
	}
	d := Detrend(x)
	if !approxSlice(d, []float64{0, 0, 0, 0, 0}, 1e-9) {
		t.Errorf("Detrend of perfect line = %v", d)
	}
}

// --- ACF/PACF ---

func TestAutoCorrelation(t *testing.T) {
	x := []float64{1, 2, 3, 4, 5}
	acf := AutoCorrelation(x, 4)
	want := []float64{1, 0.4, -0.1, -0.4, -0.4}
	if !approxSlice(acf, want, tol) {
		t.Errorf("AutoCorrelation = %v, want %v", acf, want)
	}
	if got := AutoCovarianceAt(x, 0); !approxEqual(got, 2, tol) {
		t.Errorf("AutoCovarianceAt(0) = %v, want 2", got)
	}
}

func TestPACF(t *testing.T) {
	x := []float64{1, 2, 3, 4, 5}
	pacf := PartialAutoCorrelation(x, 2)
	want := []float64{1, 0.4, -0.30952380952380953}
	if !approxSlice(pacf, want, 1e-9) {
		t.Errorf("PACF = %v, want %v", pacf, want)
	}
}

func TestCrossCorrelation(t *testing.T) {
	x := []float64{1, 2, 3, 4}
	// Perfectly correlated with itself at lag 0.
	if got := CrossCorrelationAt(x, x, 0); !approxEqual(got, 1, tol) {
		t.Errorf("CrossCorrelationAt(0) = %v, want 1", got)
	}
	cc := CrossCorrelation(x, x, 1)
	if len(cc) != 3 || !approxEqual(cc[1], 1, tol) {
		t.Errorf("CrossCorrelation = %v", cc)
	}
}

func TestDurbinWatson(t *testing.T) {
	// Alternating residuals: num = 5·2² = 20, den = 6, so DW = 20/6.
	e := []float64{1, -1, 1, -1, 1, -1}
	dw := DurbinWatson(e)
	if !approxEqual(dw, 20.0/6.0, 1e-9) {
		t.Errorf("DurbinWatson alternating = %v, want %v", dw, 20.0/6.0)
	}
	// Uncorrelated-looking residuals give DW near 2 for a constant series step.
	if dw <= 2 {
		t.Errorf("DurbinWatson should exceed 2 for negative autocorrelation, got %v", dw)
	}
}

// --- moving averages ---

func TestMovingAverages(t *testing.T) {
	if got := MovingAverage([]float64{1, 2, 3, 4}, 2); !approxSlice(got, []float64{1, 1.5, 2.5, 3.5}, tol) {
		t.Errorf("MovingAverage = %v", got)
	}
	if got := MovingAverageCentered([]float64{1, 2, 3, 4, 5}, 3); !approxSlice(got, []float64{1.5, 2, 3, 4, 4.5}, tol) {
		t.Errorf("MovingAverageCentered = %v", got)
	}
	if got := MovingAverageValid([]float64{1, 2, 3, 4}, 2); !approxSlice(got, []float64{1.5, 2.5, 3.5}, tol) {
		t.Errorf("MovingAverageValid = %v", got)
	}
	if got := ExponentialMovingAverage([]float64{1, 2, 3}, 0.5); !approxSlice(got, []float64{1, 1.5, 2.25}, tol) {
		t.Errorf("EMA = %v", got)
	}
	if got := CumulativeMovingAverage([]float64{1, 2, 3, 4}); !approxSlice(got, []float64{1, 1.5, 2, 2.5}, tol) {
		t.Errorf("CumulativeMovingAverage = %v", got)
	}
}

func TestRolling(t *testing.T) {
	x := []float64{1, 2, 3, 4, 5}
	if got := RollingSum(x, 2); !approxSlice(got, []float64{3, 5, 7, 9}, tol) {
		t.Errorf("RollingSum = %v", got)
	}
	if got := RollingMean(x, 2); !approxSlice(got, []float64{1.5, 2.5, 3.5, 4.5}, tol) {
		t.Errorf("RollingMean = %v", got)
	}
	if got := RollingMin(x, 3); !approxSlice(got, []float64{1, 2, 3}, tol) {
		t.Errorf("RollingMin = %v", got)
	}
	if got := RollingMax(x, 3); !approxSlice(got, []float64{3, 4, 5}, tol) {
		t.Errorf("RollingMax = %v", got)
	}
	if got := RollingMedian(x, 3); !approxSlice(got, []float64{2, 3, 4}, tol) {
		t.Errorf("RollingMedian = %v", got)
	}
	if got := ExpandingMax(x); !approxSlice(got, []float64{1, 2, 3, 4, 5}, tol) {
		t.Errorf("ExpandingMax = %v", got)
	}
}

func TestWeightedTriangular(t *testing.T) {
	// Weighted MA with equal weights equals simple MA.
	got := WeightedMovingAverage([]float64{2, 4, 6}, []float64{1, 1})
	if !approxSlice(got, []float64{2, 3, 5}, tol) {
		t.Errorf("WeightedMovingAverage = %v", got)
	}
	tri := TriangularMovingAverage([]float64{1, 1, 1, 1}, 2)
	if !approxSlice(tri, []float64{1, 1, 1, 1}, tol) {
		t.Errorf("TriangularMovingAverage of const = %v", tri)
	}
}

// --- smoothing ---

func TestSES(t *testing.T) {
	got := SimpleExponentialSmoothing([]float64{3, 5, 9, 20}, 0.5)
	if !approxSlice(got, []float64{3, 4, 6.5, 13.25}, tol) {
		t.Errorf("SES = %v", got)
	}
	fc := SESForecast([]float64{3, 5, 9, 20}, 0.5, 3)
	if !approxSlice(fc, []float64{13.25, 13.25, 13.25}, tol) {
		t.Errorf("SESForecast = %v", fc)
	}
}

func TestHoltLinear(t *testing.T) {
	m, err := HoltLinear([]float64{1, 2, 3, 4}, 0.5, 0.5)
	if err != nil {
		t.Fatal(err)
	}
	fc := m.Forecast(2)
	if !approxSlice(fc, []float64{5, 6}, tol) {
		t.Errorf("Holt Forecast = %v, want [5 6]", fc)
	}
	fitted := m.Fitted()
	if !approxSlice(fitted, []float64{1, 2, 3, 4}, tol) {
		t.Errorf("Holt Fitted = %v", fitted)
	}
	// A perfect linear series is fit exactly, so SSE is ~0.
	if sse := m.SSE([]float64{1, 2, 3, 4}); sse > 1e-9 {
		t.Errorf("Holt SSE = %v, want ~0", sse)
	}
}

func TestHoltLinearErrors(t *testing.T) {
	if _, err := HoltLinear([]float64{1}, 0.5, 0.5); err == nil {
		t.Error("expected error for short series")
	}
	if _, err := HoltLinear([]float64{1, 2, 3}, 1.5, 0.5); err == nil {
		t.Error("expected error for bad alpha")
	}
}

func TestHoltWinters(t *testing.T) {
	// Trend + additive seasonality, period 4.
	n := 24
	x := make([]float64, n)
	season := []float64{0, 2, -2, 0}
	for i := 0; i < n; i++ {
		x[i] = float64(i) + season[i%4]
	}
	m, err := HoltWinters(x, 0.4, 0.1, 0.3, 4, false)
	if err != nil {
		t.Fatal(err)
	}
	fc := m.Forecast(4)
	if len(fc) != 4 {
		t.Fatalf("forecast length = %d", len(fc))
	}
	// The upward trend means forecasts should exceed the recent mean.
	if fc[0] < x[n-1]-5 {
		t.Errorf("HoltWinters forecast too low: %v", fc)
	}
	// Multiplicative path should also fit.
	if _, err := HoltWinters(x, 0.4, 0.1, 0.3, 4, true); err != nil {
		t.Errorf("multiplicative HoltWinters error: %v", err)
	}
}

func TestBrownDouble(t *testing.T) {
	// Perfect linear trend: Brown's method should extrapolate the slope.
	x := []float64{1, 2, 3, 4, 5, 6, 7, 8}
	fc := BrownDoubleExponential(x, 0.5, 2)
	if len(fc) != 2 {
		t.Fatalf("length = %d", len(fc))
	}
	if fc[1]-fc[0] < 0.5 {
		t.Errorf("Brown forecast slope too small: %v", fc)
	}
}

// --- AR models ---

func makeAR1(phi float64, n int) []float64 {
	seed := 20240101.0
	rnd := func() float64 {
		seed = float64(int64(seed*1103515245+12345) % 2147483648)
		return seed/2147483648 - 0.5
	}
	x := make([]float64, n)
	for i := 1; i < n; i++ {
		x[i] = phi*x[i-1] + rnd()
	}
	return x
}

func TestYuleWalker(t *testing.T) {
	x := makeAR1(0.7, 500)
	m, err := YuleWalker(x, 1)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(m.Phi[0]-0.7) > 0.1 {
		t.Errorf("YuleWalker phi = %v, want ~0.7", m.Phi[0])
	}
	if m.Sigma2 <= 0 {
		t.Errorf("Sigma2 = %v, want > 0", m.Sigma2)
	}
}

func TestARFitLeastSquares(t *testing.T) {
	x := makeAR1(0.5, 500)
	m, err := ARFitLeastSquares(x, 1)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(m.Phi[0]-0.5) > 0.1 {
		t.Errorf("LS phi = %v, want ~0.5", m.Phi[0])
	}
}

func TestBurgAR(t *testing.T) {
	x := makeAR1(0.6, 500)
	m, err := BurgAR(x, 1)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(m.Phi[0]-0.6) > 0.1 {
		t.Errorf("Burg phi = %v, want ~0.6", m.Phi[0])
	}
}

func TestARForecastValues(t *testing.T) {
	m := &ARModel{Order: 1, Phi: []float64{0.5}, Mean: 0, Intercept: 0}
	x := []float64{8}
	fc := m.Forecast(x, 3)
	if !approxSlice(fc, []float64{4, 2, 1}, tol) {
		t.Errorf("AR Forecast = %v, want [4 2 1]", fc)
	}
}

func TestLevinsonDurbin(t *testing.T) {
	// AR(1) autocovariance gamma_k = phi^k * sigma^2/(1-phi^2).
	phi := 0.6
	g0 := 1.0
	gamma := []float64{g0, phi * g0, phi * phi * g0}
	coeffs, v := LevinsonDurbin(gamma, 1)
	if math.Abs(coeffs[0]-phi) > tol {
		t.Errorf("Levinson phi = %v, want %v", coeffs[0], phi)
	}
	if v <= 0 {
		t.Errorf("innovation variance = %v", v)
	}
}

// --- ARMA/MA/ARIMA ---

func TestMAFit(t *testing.T) {
	m, err := MAFit([]float64{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, 2)
	if err != nil {
		t.Fatal(err)
	}
	if len(m.Theta) != 2 {
		t.Errorf("MA order = %d", len(m.Theta))
	}
	// Forecasts beyond the MA order revert to the mean.
	fc := m.Forecast([]float64{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, 5)
	if !approxEqual(fc[3], m.Mean, 1e-9) || !approxEqual(fc[4], m.Mean, 1e-9) {
		t.Errorf("MA forecast should revert to mean %v: %v", m.Mean, fc)
	}
}

func TestARMAFit(t *testing.T) {
	x := makeAR1(0.6, 400)
	m, err := ARMAFit(x, 1, 1)
	if err != nil {
		t.Fatal(err)
	}
	if len(m.Phi) != 1 || len(m.Theta) != 1 {
		t.Errorf("ARMA orders wrong: %+v", m)
	}
	fc := m.Forecast(x, 3)
	if len(fc) != 3 {
		t.Errorf("ARMA forecast length = %d", len(fc))
	}
	for _, v := range fc {
		if math.IsNaN(v) || math.IsInf(v, 0) {
			t.Errorf("ARMA forecast has bad value: %v", fc)
		}
	}
}

func TestARIMAFit(t *testing.T) {
	// Random walk with drift.
	seed := 7.0
	rnd := func() float64 {
		seed = float64(int64(seed*1103515245+12345) % 2147483648)
		return seed/2147483648 - 0.5
	}
	n := 200
	x := make([]float64, n)
	for i := 1; i < n; i++ {
		x[i] = x[i-1] + 0.3 + rnd()
	}
	m, err := ARIMAFit(x, 1, 1, 0)
	if err != nil {
		t.Fatal(err)
	}
	fc := m.Forecast(5)
	if len(fc) != 5 {
		t.Fatalf("ARIMA forecast length = %d", len(fc))
	}
	// A positive-drift random walk should keep rising.
	if fc[4] <= fc[0] {
		t.Errorf("ARIMA forecast not rising: %v", fc)
	}
}

func TestARMAToMA(t *testing.T) {
	// AR(1): psi_j = phi^j.
	psi := ARMAToMA([]float64{0.5}, nil, 4)
	if !approxSlice(psi, []float64{0.5, 0.25, 0.125, 0.0625}, tol) {
		t.Errorf("ARMAToMA = %v", psi)
	}
}

func TestARMAToAR(t *testing.T) {
	// MA(1) with theta=0.6: pi_k = -(-0.6)^k.
	pi := ARMAToAR(nil, []float64{0.6}, 4)
	want := []float64{0.6, -0.36, 0.216, -0.1296}
	if !approxSlice(pi, want, tol) {
		t.Errorf("ARMAToAR = %v, want %v", pi, want)
	}
}

func TestInnovationsAlgorithm(t *testing.T) {
	gamma := []float64{2, 1, 0}
	theta, v := InnovationsAlgorithm(gamma, 2)
	if len(v) != 3 || v[0] != 2 {
		t.Errorf("innovations v = %v", v)
	}
	// theta[1][1] = gamma[1]/v[0] = 0.5.
	if !approxEqual(theta[1][1], 0.5, tol) {
		t.Errorf("theta[1][1] = %v, want 0.5", theta[1][1])
	}
}

// --- spectral ---

func TestDFTRoundTrip(t *testing.T) {
	x := []float64{1, 2, 3, 4}
	X := DFT(x)
	back := InverseDFT(X)
	for i := range x {
		if math.Abs(real(back[i])-x[i]) > 1e-9 || math.Abs(imag(back[i])) > 1e-9 {
			t.Errorf("DFT round-trip[%d] = %v, want %v", i, back[i], x[i])
		}
	}
}

func TestDFTConstant(t *testing.T) {
	// DFT of a constant has all energy at frequency 0.
	x := []float64{5, 5, 5, 5}
	X := DFT(x)
	if math.Abs(real(X[0])-20) > 1e-9 {
		t.Errorf("X[0] = %v, want 20", X[0])
	}
	for k := 1; k < 4; k++ {
		if math.Abs(real(X[k])) > 1e-9 || math.Abs(imag(X[k])) > 1e-9 {
			t.Errorf("X[%d] = %v, want 0", k, X[k])
		}
	}
}

func TestPeriodogram(t *testing.T) {
	// A pure sinusoid of period 4 over 64 samples peaks at that period.
	n := 64
	x := make([]float64, n)
	for i := 0; i < n; i++ {
		x[i] = math.Sin(2 * math.Pi * float64(i) / 4)
	}
	if p := DominantPeriod(x); math.Abs(p-4) > 1e-6 {
		t.Errorf("DominantPeriod = %v, want 4", p)
	}
	freqs, power := Periodogram(x)
	if len(freqs) != n/2+1 || len(power) != n/2+1 {
		t.Errorf("periodogram lengths = %d,%d", len(freqs), len(power))
	}
}

func TestSpectralEntropy(t *testing.T) {
	// White-ish noise has high spectral entropy; a sinusoid has low.
	n := 128
	sine := make([]float64, n)
	for i := 0; i < n; i++ {
		sine[i] = math.Sin(2 * math.Pi * float64(i) / 8)
	}
	seed := 3.0
	rnd := func() float64 {
		seed = float64(int64(seed*1103515245+12345) % 2147483648)
		return seed/2147483648 - 0.5
	}
	noise := make([]float64, n)
	for i := range noise {
		noise[i] = rnd()
	}
	if SpectralEntropy(sine) >= SpectralEntropy(noise) {
		t.Errorf("expected sine entropy < noise entropy: %v vs %v",
			SpectralEntropy(sine), SpectralEntropy(noise))
	}
}

func TestARSpectralDensity(t *testing.T) {
	m := &ARModel{Order: 1, Phi: []float64{0.5}, Sigma2: 1}
	freqs, dens := ARSpectralDensity(m, 5)
	if len(freqs) != 5 || len(dens) != 5 {
		t.Fatalf("lengths = %d,%d", len(freqs), len(dens))
	}
	// For phi>0 the AR(1) spectrum is largest at frequency 0.
	if dens[0] <= dens[4] {
		t.Errorf("AR spectrum not peaked at 0: %v", dens)
	}
}

// --- decomposition ---

func TestSeasonalDecompose(t *testing.T) {
	x := []float64{1, 3, 2, 5, 2, 4, 3, 6, 3, 5, 4, 7}
	d := SeasonalDecompose(x, 4, false)
	if d == nil {
		t.Fatal("decomposition nil")
	}
	// Additive seasonal indices sum to zero.
	var s float64
	for i := 0; i < 4; i++ {
		s += d.Seasonal[i]
	}
	if math.Abs(s) > 1e-9 {
		t.Errorf("seasonal indices sum = %v, want 0", s)
	}
	// Reconstruction where trend is defined.
	for i := range x {
		if math.IsNaN(d.Trend[i]) {
			continue
		}
		recon := d.Trend[i] + d.Seasonal[i] + d.Residual[i]
		if math.Abs(recon-x[i]) > 1e-9 {
			t.Errorf("reconstruction[%d] = %v, want %v", i, recon, x[i])
		}
	}
}

func TestSeasonallyAdjust(t *testing.T) {
	x := []float64{1, 3, 2, 5, 2, 4, 3, 6, 3, 5, 4, 7}
	adj := SeasonallyAdjust(x, 4, false)
	if len(adj) != len(x) {
		t.Errorf("adjusted length = %d", len(adj))
	}
}

// --- stationarity ---

func TestADF(t *testing.T) {
	// Stationary AR(1) should reject the unit-root null (very negative stat).
	x := makeAR1(0.3, 300)
	res, err := ADFTest(x, 1, false)
	if err != nil {
		t.Fatal(err)
	}
	if res.Statistic > ADFCriticalValue(0.05, "c") {
		t.Errorf("ADF stat = %v not below 5%% critical value", res.Statistic)
	}
	if !IsStationaryADF(x, 1, 0.05) {
		t.Error("stationary AR(1) not flagged stationary")
	}
}

func TestADFRandomWalk(t *testing.T) {
	seed := 42.0
	rnd := func() float64 {
		seed = float64(int64(seed*1103515245+12345) % 2147483648)
		return seed/2147483648 - 0.5
	}
	n := 300
	rw := make([]float64, n)
	for i := 1; i < n; i++ {
		rw[i] = rw[i-1] + rnd()
	}
	if IsStationaryADF(rw, 1, 0.05) {
		t.Error("random walk should not be flagged stationary")
	}
	if NumberOfDifferences(rw, 1, 2, 0.05) < 1 {
		t.Error("random walk should need at least one difference")
	}
}

func TestKPSS(t *testing.T) {
	// Stationary noise -> small statistic; random walk -> large.
	seed := 11.0
	rnd := func() float64 {
		seed = float64(int64(seed*1103515245+12345) % 2147483648)
		return seed/2147483648 - 0.5
	}
	n := 300
	noise := make([]float64, n)
	rw := make([]float64, n)
	for i := 0; i < n; i++ {
		noise[i] = rnd()
		if i > 0 {
			rw[i] = rw[i-1] + noise[i]
		}
	}
	kn, err := KPSSTest(noise, 5, false)
	if err != nil {
		t.Fatal(err)
	}
	krw, err := KPSSTest(rw, 5, false)
	if err != nil {
		t.Fatal(err)
	}
	if kn.Statistic >= krw.Statistic {
		t.Errorf("KPSS noise (%v) should be below random walk (%v)", kn.Statistic, krw.Statistic)
	}
}

func TestVarianceRatio(t *testing.T) {
	// Random walk should have variance ratio near 1.
	seed := 5.0
	rnd := func() float64 {
		seed = float64(int64(seed*1103515245+12345) % 2147483648)
		return seed/2147483648 - 0.5
	}
	n := 2000
	rw := make([]float64, n)
	for i := 1; i < n; i++ {
		rw[i] = rw[i-1] + rnd()
	}
	vr := VarianceRatio(rw, 2)
	if math.Abs(vr-1) > 0.15 {
		t.Errorf("VarianceRatio random walk = %v, want ~1", vr)
	}
}

// --- embedding ---

func TestLagMatrix(t *testing.T) {
	x := []float64{10, 20, 30, 40}
	X, y := LagMatrix(x, 2)
	if len(X) != 2 || len(y) != 2 {
		t.Fatalf("dims = %d,%d", len(X), len(y))
	}
	if !approxSlice(X[0], []float64{20, 10}, tol) || y[0] != 30 {
		t.Errorf("row0 = %v, y0 = %v", X[0], y[0])
	}
	if !approxSlice(X[1], []float64{30, 20}, tol) || y[1] != 40 {
		t.Errorf("row1 = %v, y1 = %v", X[1], y[1])
	}
}

func TestEmbed(t *testing.T) {
	x := []float64{1, 2, 3, 4, 5}
	e := Embed(x, 2, 1)
	if len(e) != 4 {
		t.Fatalf("rows = %d", len(e))
	}
	if !approxSlice(e[0], []float64{1, 2}, tol) || !approxSlice(e[3], []float64{4, 5}, tol) {
		t.Errorf("embed = %v", e)
	}
}

func TestHankelToeplitz(t *testing.T) {
	h := HankelMatrix([]float64{1, 2, 3, 4}, 2)
	if !approxSlice(h[0], []float64{1, 2, 3}, tol) || !approxSlice(h[1], []float64{2, 3, 4}, tol) {
		t.Errorf("Hankel = %v", h)
	}
	tp := ToeplitzMatrix([]float64{1, 2, 3})
	if !approxSlice(tp[0], []float64{1, 2, 3}, tol) || !approxSlice(tp[2], []float64{3, 2, 1}, tol) {
		t.Errorf("Toeplitz = %v", tp)
	}
}

func TestSlidingWindows(t *testing.T) {
	w := SlidingWindows([]float64{1, 2, 3, 4, 5}, 2, 2)
	if len(w) != 2 {
		t.Fatalf("windows = %d", len(w))
	}
	if !approxSlice(w[0], []float64{1, 2}, tol) || !approxSlice(w[1], []float64{3, 4}, tol) {
		t.Errorf("windows = %v", w)
	}
}

// --- metrics ---

func TestMetrics(t *testing.T) {
	actual := []float64{1, 2, 3, 4}
	forecast := []float64{1.1, 1.9, 3.2, 3.8}
	mae := MeanAbsoluteError(actual, forecast)
	if !approxEqual(mae, 0.15, 1e-9) {
		t.Errorf("MAE = %v, want 0.15", mae)
	}
	mse := MeanSquaredError(actual, forecast)
	if !approxEqual(mse, 0.025, 1e-9) {
		t.Errorf("MSE = %v, want 0.025", mse)
	}
	if !approxEqual(RootMeanSquaredError(actual, forecast), math.Sqrt(0.025), 1e-9) {
		t.Errorf("RMSE mismatch")
	}
	// Perfect forecast: R^2 = 1.
	if r := RSquared(actual, actual); !approxEqual(r, 1, tol) {
		t.Errorf("RSquared perfect = %v, want 1", r)
	}
	if me := MeanError(actual, forecast); math.IsNaN(me) {
		t.Error("MeanError NaN")
	}
}

func TestMASE(t *testing.T) {
	train := []float64{1, 2, 3, 4, 5}
	actual := []float64{6, 7}
	forecast := []float64{6, 7}
	// Perfect forecast -> MASE 0.
	if m := MeanAbsoluteScaledError(actual, forecast, train, 1); !approxEqual(m, 0, tol) {
		t.Errorf("MASE perfect = %v, want 0", m)
	}
}

func TestTheilU(t *testing.T) {
	// Perfect one-step forecast -> U = 0.
	actual := []float64{1, 2, 3, 4}
	if u := TheilU(actual, actual); !approxEqual(u, 0, tol) {
		t.Errorf("TheilU perfect = %v, want 0", u)
	}
}

// --- solver sanity ---

func TestSolveLinear(t *testing.T) {
	A := [][]float64{{2, 1}, {1, 3}}
	b := []float64{3, 5}
	z, ok := solveLinear(A, b)
	if !ok {
		t.Fatal("solveLinear failed")
	}
	// Solution: x = 0.8, y = 1.4.
	if !approxEqual(z[0], 0.8, 1e-9) || !approxEqual(z[1], 1.4, 1e-9) {
		t.Errorf("solveLinear = %v", z)
	}
}

func TestInvertMatrix(t *testing.T) {
	A := [][]float64{{4, 7}, {2, 6}}
	inv, ok := invertMatrix(A)
	if !ok {
		t.Fatal("invertMatrix failed")
	}
	want := [][]float64{{0.6, -0.7}, {-0.2, 0.4}}
	for i := range want {
		if !approxSlice(inv[i], want[i], 1e-9) {
			t.Errorf("inv row %d = %v, want %v", i, inv[i], want[i])
		}
	}
}

// --- runnable example ---

func ExampleAutoCorrelation() {
	x := []float64{1, 2, 3, 4, 5}
	acf := AutoCorrelation(x, 2)
	fmt.Printf("%.1f %.1f %.1f\n", acf[0], acf[1], acf[2])
	// Output: 1.0 0.4 -0.1
}

func ExampleHoltLinear() {
	x := []float64{1, 2, 3, 4}
	m, _ := HoltLinear(x, 0.5, 0.5)
	fmt.Println(m.Forecast(2))
	// Output: [5 6]
}
