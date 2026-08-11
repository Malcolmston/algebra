# timeseries — self-contained, dependency-free time-series analysis toolkit written entirely

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/timeseries.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/timeseries)

Package timeseries is a self-contained, dependency-free time-series analysis
toolkit written entirely with the Go standard library (only math, math/cmplx,
sort and errors). It operates on plain []float64 sample vectors and does not
import the parent algebra module.

The package covers the classical time-series workflow end to end.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/timeseries"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func ACFConfidenceBound(n int) float64` | ACFConfidenceBound returns the approximate two-sided 95% confidence bound 1.96/√n for the autocorrelations of white noise of length n. |
| `func ADFCriticalValue(level float64, trend string) float64` | ADFCriticalValue returns an approximate MacKinnon asymptotic critical value for the ADF t-statistic at the given significance level (0.01, 0.05 or… |
| `func ARMAToAR(phi, theta []float64, n int) []float64` | ARMAToAR returns the first n coefficients π_1,…,π_n of the AR(∞) representation Σ π_j x_{t−j} = e_t implied by the AR coefficients phi and… |
| `func ARMAToMA(phi, theta []float64, n int) []float64` | ARMAToMA returns the first n coefficients ψ_1,…,ψ_n of the MA(∞) representation x_t = Σ ψ_j e_{t−j} implied by the AR coefficients phi and… |
| `func ARSpectralDensity(m *ARModel, nf int) (freqs, density []float64)` | ARSpectralDensity evaluates the AR power spectral density implied by the model at nf equally spaced frequencies in [0, 0.5]. |
| `func Argmax(x []float64) int` | Argmax returns the index of the largest value, or -1 if the series is empty. |
| `func Argmin(x []float64) int` | Argmin returns the index of the smallest value, or -1 if the series is empty. |
| `func AutoCorrelation(x []float64, maxlag int) []float64` | AutoCorrelation returns the sample autocorrelation function (ACF) for lags 0,1,…,maxlag. |
| `func AutoCorrelationAt(x []float64, k int) float64` | AutoCorrelationAt returns the sample autocorrelation at lag k, i.e. |
| `func AutoCovariance(x []float64, maxlag int) []float64` | AutoCovariance returns the biased sample autocovariances of the series for lags 0,1,…,maxlag. |
| `func AutoCovarianceAt(x []float64, k int) float64` | AutoCovarianceAt returns the biased sample autocovariance of the series at lag k: (1/n)·Σ (x_t − x̄)(x_{t+k} − x̄). |
| `func AutocorrelationMatrix(x []float64, p int) [][]float64` | AutocorrelationMatrix returns the p×p symmetric Toeplitz autocorrelation matrix of the series, with element (i,j) equal to the sample… |
| `func BoxCox(x []float64, lambda float64) []float64` | BoxCox applies the Box–Cox power transform with parameter lambda. |
| `func BoxPierce(x []float64, h int) float64` | BoxPierce returns the Box–Pierce Q statistic for the first h autocorrelation lags, the simpler large-sample precursor to `LjungBox`. |
| `func BrownDoubleExponential(x []float64, alpha float64, h int) []float64` | BrownDoubleExponential returns Brown's single-parameter double exponential smoothing (linear-trend) forecast for the series. |
| `func Clip(x []float64, lo, hi float64) []float64` | Clip constrains every observation to the closed interval [lo,hi]. |
| `func CoefficientOfVariation(x []float64) float64` | CoefficientOfVariation returns the ratio of the sample standard deviation to the mean, a scale-free measure of dispersion. |
| `func CrossCorrelation(x, y []float64, maxlag int) []float64` | CrossCorrelation returns the cross-correlation function of x and y for lags −maxlag,…,0,…,maxlag as a slice of length 2·maxlag+1, indexed so… |
| `func CrossCorrelationAt(x, y []float64, k int) float64` | CrossCorrelationAt returns the sample cross-correlation between x and y at lag k, normalized by the product of the two series' standard deviations. |
| `func CrossCovarianceAt(x, y []float64, k int) float64` | CrossCovarianceAt returns the sample cross-covariance between x and y at lag k: (1/n)·Σ (x_t − x̄)(y_{t+k} − ȳ) for k ≥ 0, and the… |
| `func CumProd(x []float64) []float64` | CumProd returns the cumulative product of the series. |
| `func CumSum(x []float64) []float64` | CumSum returns the cumulative sum of the series: out`i` = x[0]+…+x`i`. |
| `func CumulativeMovingAverage(x []float64) []float64` | CumulativeMovingAverage returns the running (expanding-window) mean of the series: out`i` is the mean of x[0..i]. |
| `func CumulativePeriodogram(x []float64) []float64` | CumulativePeriodogram returns the normalized cumulative periodogram of the series over the non-zero Fourier frequencies, a monotone sequence rising… |
| `func DFT(x []float64) []complex128` | DFT returns the discrete Fourier transform of the real series x, a complex spectrum of the same length computed by direct summation: X_k = Σ_t… |
| `func Demean(x []float64) []float64` | Demean returns the series with its mean subtracted from every observation. |
| `func Detrend(x []float64) []float64` | Detrend removes a least-squares linear trend from the series, returning the residuals x`i` − (a + b·i). |
| `func Diff(x []float64) []float64` | Diff returns the first difference of the series: out`i` = x[i+1] − x`i`. |
| `func DiffOrder(x []float64, d int) []float64` | DiffOrder applies the first difference d times. |
| `func DominantFrequency(x []float64) float64` | DominantFrequency returns the non-zero Fourier frequency with the largest periodogram power, i.e. |
| `func DominantPeriod(x []float64) float64` | DominantPeriod returns the reciprocal of `DominantFrequency`, i.e. |
| `func DoubleExponentialMovingAverage(x []float64, alpha float64) []float64` | DoubleExponentialMovingAverage returns the DEMA indicator, 2·EMA(x) − EMA(EMA(x)), which reduces the lag of a plain EMA. |
| `func DurbinWatson(e []float64) float64` | DurbinWatson returns the Durbin–Watson statistic of a residual series, Σ(e_t − e_{t−1})² / Σe_t², which detects lag-1 autocorrelation. |
| `func Embed(x []float64, m, tau int) [][]float64` | Embed constructs a time-delay (Takens) embedding of the series with embedding dimension m and delay tau. |
| `func Energy(x []float64) float64` | Energy returns the sum of squared observations of the series. |
| `func EstimateSeasonalPeriod(x []float64, maxLag int) int` | EstimateSeasonalPeriod returns the lag (≥ 2) of the largest peak in the autocorrelation function up to maxLag, a data-driven estimate of the… |
| `func ExpTransform(x []float64) []float64` | ExpTransform returns the exponential of each observation (the inverse of `LogTransform`). |
| `func ExpandingMax(x []float64) []float64` | ExpandingMax returns the running maximum: out`i` = max(x[0..i]). |
| `func ExpandingMean(x []float64) []float64` | ExpandingMean returns the expanding-window mean, identical to `CumulativeMovingAverage`. |
| `func ExpandingMin(x []float64) []float64` | ExpandingMin returns the running minimum: out`i` = min(x[0..i]). |
| `func ExpandingSum(x []float64) []float64` | ExpandingSum returns the expanding-window sum, identical to `CumSum`. |
| `func ExponentialMovingAverage(x []float64, alpha float64) []float64` | ExponentialMovingAverage returns the exponentially weighted moving average with smoothing factor alpha in (0,1]: out[0] = x[0] and out`i` =… |
| `func ExponentialMovingAverageSpan(x []float64, span int) []float64` | ExponentialMovingAverageSpan returns the exponential moving average using the span convention alpha = 2/(span+1), common in technical analysis. |
| `func First(x []float64) float64` | First returns the first observation of the series, or NaN if empty. |
| `func FlattenMatrix(m [][]float64) []float64` | FlattenMatrix returns a row-major flattening of a matrix into a single slice, a convenience for feeding embedded matrices to routines that expect… |
| `func FourierFrequencies(n int, d float64) []float64` | FourierFrequencies returns the n Fourier frequencies (in cycles per sample unit) associated with a length-n DFT: k/(n·d) for k = 0,…,n−1, where… |
| `func FracDiffWeights(d float64, n int) []float64` | FracDiffWeights returns the first n binomial weights of the fractional differencing operator (1−B)^d. |
| `func FractionalDifference(x []float64, d float64) []float64` | FractionalDifference applies the fractional differencing operator (1−B)^d to the series using an expanding window of binomial weights. |
| `func HankelMatrix(x []float64, rows int) [][]float64` | HankelMatrix builds the Hankel matrix of the series with the given number of rows: element (i,j) is x[i+j]. |
| `func InnovationsAlgorithm(gamma []float64, maxlag int) ([][]float64, []float64)` | InnovationsAlgorithm runs the innovations algorithm on the autocovariance sequence gamma (length ≥ maxlag+1) and returns the innovation… |
| `func Integrate(d []float64, x0 float64) []float64` | Integrate is the inverse of `Diff`: given the first-difference series d and an initial value x0, it reconstructs the original series by cumulative… |
| `func InverseBoxCox(y []float64, lambda float64) []float64` | InverseBoxCox inverts the Box–Cox transform with parameter lambda. |
| `func InverseDFT(X []complex128) []complex128` | InverseDFT returns the inverse discrete Fourier transform of the complex spectrum X, a complex series of the same length. |
| `func IsStationaryADF(x []float64, lags int, level float64) bool` | IsStationaryADF reports whether the augmented Dickey–Fuller test rejects the unit-root null at the given significance level (typically 0.05), i.e. |
| `func Kurtosis(x []float64) float64` | Kurtosis returns the excess kurtosis (fourth standardized moment minus 3) of the series. |
| `func Lag(x []float64, k int) []float64` | Lag shifts the series forward by k steps, padding the first k positions with NaN. The result has the same length as x. |
| `func LagMatrix(x []float64, p int) ([][]float64, []float64)` | LagMatrix builds a design matrix of lagged values for autoregressive modeling. |
| `func Last(x []float64) float64` | Last returns the final observation of the series, or NaN if empty. |
| `func Lead(x []float64, k int) []float64` | Lead shifts the series backward by k steps, padding the last k positions with NaN. The result has the same length as x. |
| `func LevinsonDurbin(gamma []float64, p int) ([]float64, float64)` | LevinsonDurbin solves the Yule–Walker equations for an AR(p) model from the autocovariance sequence gamma (which must have length at least p+1). |
| `func LjungBox(x []float64, h int) float64` | LjungBox returns the Ljung–Box Q statistic for the first h autocorrelation lags, a test for overall autocorrelation. |
| `func LogReturns(x []float64) []float64` | LogReturns returns the log returns of the series: out`i` = ln(x[i+1]/x`i`). |
| `func LogTransform(x []float64) []float64` | LogTransform returns the natural logarithm of each observation. |
| `func Max(x []float64) float64` | Max returns the largest value in the series, or NaN if it is empty. |
| `func Mean(x []float64) float64` | Mean returns the arithmetic mean of the series, or NaN if it is empty. |
| `func MeanAbsoluteDeviation(x []float64) float64` | MeanAbsoluteDeviation returns the mean of the absolute deviations from the series mean. |
| `func MeanAbsoluteError(actual, forecast []float64) float64` | MeanAbsoluteError returns the mean of \|actual − forecast\| over the aligned pairs. |
| `func MeanAbsolutePercentageError(actual, forecast []float64) float64` | MeanAbsolutePercentageError returns the mean of \|(actual − forecast)/actual\| as a fraction (multiply by 100 for a percentage). |
| `func MeanAbsoluteScaledError(actual, forecast, train []float64, m int) float64` | MeanAbsoluteScaledError returns the MASE of the forecast, the mean absolute error scaled by the mean absolute error of a naive one-step… |
| `func MeanError(actual, forecast []float64) float64` | MeanError returns the mean forecast error (bias) actual − forecast over the aligned pairs. |
| `func MeanSquaredError(actual, forecast []float64) float64` | MeanSquaredError returns the mean of (actual − forecast)² over the aligned pairs. |
| `func Median(x []float64) float64` | Median returns the median of the series, or NaN if it is empty. |
| `func MedianAbsoluteError(actual, forecast []float64) float64` | MedianAbsoluteError returns the median of \|actual − forecast\| over the aligned pairs, a robust accuracy measure. |
| `func MedianFilter(x []float64, w int) []float64` | MedianFilter returns the centered running median with window length w, truncating the window at the boundaries so the output length equals the input… |
| `func Min(x []float64) float64` | Min returns the smallest value in the series, or NaN if it is empty. |
| `func MinMaxNormalize(x []float64) []float64` | MinMaxNormalize rescales the series linearly to the unit interval [0,1]. |
| `func MovingAverage(x []float64, w int) []float64` | MovingAverage returns the causal (trailing) simple moving average with window length w: out`i` is the mean of the up-to-w most recent samples ending… |
| `func MovingAverageCentered(x []float64, w int) []float64` | MovingAverageCentered returns the centered simple moving average with window length w, truncating the window at the boundaries so the output length… |
| `func MovingAverageValid(x []float64, w int) []float64` | MovingAverageValid returns the centered simple moving average over full windows only, so the result has length len(x)−w+1. |
| `func NumberOfDifferences(x []float64, lags, maxD int, level float64) int` | NumberOfDifferences estimates the order of differencing needed to make the series stationary by repeatedly applying the ADF test (constant, given… |
| `func PartialAutoCorrelation(x []float64, maxlag int) []float64` | PartialAutoCorrelation returns the sample partial autocorrelation function (PACF) for lags 0,1,…,maxlag computed with the Durbin–Levinson… |
| `func Percentile(x []float64, p float64) float64` | Percentile returns the p-th percentile (0 ≤ p ≤ 100) of the series. |
| `func Periodogram(x []float64) (freqs, power []float64)` | Periodogram returns the one-sided raw periodogram of the series, an estimate of the power at each Fourier frequency. |
| `func PopStdDev(x []float64) float64` | PopStdDev returns the population standard deviation (square root of `PopVariance`). |
| `func PopVariance(x []float64) float64` | PopVariance returns the population (divide-by-n) variance of the series. |
| `func Quantile(x []float64, q float64) float64` | Quantile returns the q-quantile (0 ≤ q ≤ 1) of the series using linear interpolation between order statistics. |
| `func RSquared(actual, forecast []float64) float64` | RSquared returns the coefficient of determination R² = 1 − SS_res/SS_tot of the forecast relative to the mean of the actual series. |
| `func Range(x []float64) float64` | Range returns the difference between the maximum and minimum values. |
| `func ReflectionCoefficients(x []float64, p int) []float64` | ReflectionCoefficients returns the partial autocorrelation (reflection / PARCOR) coefficients k_1,…,k_p produced by the Levinson–Durbin recursion… |
| `func Rescale(x []float64, a, b float64) []float64` | Rescale linearly maps the series to the interval [a,b]. |
| `func Reverse(x []float64) []float64` | Reverse returns a new series with the observations in reverse order. |
| `func RollingMax(x []float64, w int) []float64` | RollingMax returns the trailing rolling maximum over full windows of length w. |
| `func RollingMean(x []float64, w int) []float64` | RollingMean returns the trailing rolling mean over full windows of length w, producing len(x)−w+1 values. |
| `func RollingMedian(x []float64, w int) []float64` | RollingMedian returns the trailing rolling median over full windows of length w. |
| `func RollingMin(x []float64, w int) []float64` | RollingMin returns the trailing rolling minimum over full windows of length w. |
| `func RollingStdDev(x []float64, w int) []float64` | RollingStdDev returns the trailing rolling sample standard deviation over full windows of length w. |
| `func RollingSum(x []float64, w int) []float64` | RollingSum returns the trailing rolling sum over full windows of length w. |
| `func RollingVariance(x []float64, w int) []float64` | RollingVariance returns the trailing rolling sample variance over full windows of length w. |
| `func RootMeanSquare(x []float64) float64` | RootMeanSquare returns sqrt(mean(x²)), the quadratic mean of the series. |
| `func RootMeanSquaredError(actual, forecast []float64) float64` | RootMeanSquaredError returns the square root of `MeanSquaredError`. |
| `func SESForecast(x []float64, alpha float64, h int) []float64` | SESForecast returns the h-step-ahead forecast of a simple exponential smoothing model fitted to x with parameter alpha. |
| `func SeasonalDiff(x []float64, s int) []float64` | SeasonalDiff returns the seasonal difference at lag s: out`i` = x[i+s] − x`i`. |
| `func SeasonalIndices(x []float64, period int, mult bool) []float64` | SeasonalIndices returns the per-phase seasonal factors of the series for the given period (length period): additive deviations that sum to zero when… |
| `func SeasonalIntegrate(d []float64, seed []float64, s int) []float64` | SeasonalIntegrate is the inverse of `SeasonalDiff` at lag s: given the seasonal-difference series d and the s initial values seed (the first season… |
| `func SeasonallyAdjust(x []float64, period int, mult bool) []float64` | SeasonallyAdjust removes the estimated seasonal component from the series, returning the seasonally adjusted series (trend plus residual for the… |
| `func Shift(x []float64, k int, fill float64) []float64` | Shift shifts the series forward by k steps, filling vacated positions with fill. |
| `func SimpleExponentialSmoothing(x []float64, alpha float64) []float64` | SimpleExponentialSmoothing returns the simple exponentially smoothed series with smoothing parameter alpha in (0,1): s[0] = x[0] and s`t` =… |
| `func SimpleReturns(x []float64) []float64` | SimpleReturns returns the simple (arithmetic) returns of the series: out`i` = x[i+1]/x`i` − 1. |
| `func Skewness(x []float64) float64` | Skewness returns the sample skewness (third standardized moment) using the population standard deviation in the denominator. |
| `func SlidingWindows(x []float64, w, step int) [][]float64` | SlidingWindows returns the list of consecutive full windows of length w over the series, each stepped by step positions. |
| `func SpectralEntropy(x []float64) float64` | SpectralEntropy returns the normalized spectral (Shannon) entropy of the series computed from its periodogram, a value in [0,1] that is near 1 for… |
| `func SqrtTransform(x []float64) []float64` | SqrtTransform returns the square root of each observation. |
| `func Standardize(x []float64) []float64` | Standardize returns the z-scored series (x − mean)/stddev using the sample standard deviation. |
| `func StdDev(x []float64) float64` | StdDev returns the sample standard deviation (square root of `Variance`). |
| `func Sum(x []float64) float64` | Sum returns the sum of all observations in the series. |
| `func SymmetricMAPE(actual, forecast []float64) float64` | SymmetricMAPE returns the symmetric mean absolute percentage error, the mean of \|actual − forecast\| / ((\|actual\| + \|forecast\|)/2) as a… |
| `func TakensThetaAutoMI(x []float64, maxLag int) int` | TakensThetaAutoMI returns the lag at which the sample autocorrelation of the series first drops to or below 1/e, a common heuristic for choosing the… |
| `func TheilU(actual, forecast []float64) float64` | TheilU returns Theil's U2 forecast-accuracy statistic, the ratio of the RMSE of the forecast to the RMSE of a naive no-change forecast. |
| `func TimeDelayEmbedding(x []float64, dimension, delay int) [][]float64` | TimeDelayEmbedding is an alias for `Embed` with the conventional argument order (dimension then delay), returning the delay-coordinate matrix. |
| `func ToeplitzMatrix(c []float64) [][]float64` | ToeplitzMatrix builds the symmetric Toeplitz matrix whose first row and column are c: element (i,j) is c[\|i−j\|]. |
| `func TrendComponent(x []float64, period int) []float64` | TrendComponent returns the classical moving-average trend estimate of the series for the given period, with NaN at boundary positions. |
| `func TrendLine(x []float64) []float64` | TrendLine returns the fitted trend values a + b·i for the series. |
| `func TriangularMovingAverage(x []float64, w int) []float64` | TriangularMovingAverage returns a trailing weighted moving average whose weights increase linearly from 1 up to w (giving most weight to the most… |
| `func TripleExponentialMovingAverage(x []float64, alpha float64) []float64` | TripleExponentialMovingAverage returns the TEMA indicator, 3·EMA − 3·EMA² + EMA³, which further reduces lag. |
| `func Variance(x []float64) float64` | Variance returns the sample (unbiased, divide-by-n−1) variance of the series. |
| `func VarianceRatio(x []float64, q int) float64` | VarianceRatio returns the Lo–MacKinlay variance ratio of the series at horizon q: the variance of q-period differences divided by q times the… |
| `func WeightedMovingAverage(x []float64, weights []float64) []float64` | WeightedMovingAverage returns the causal weighted moving average whose window length equals len(weights). |

### Types

| Type | What it is |
| --- | --- |
| `ADFResult` | ADFResult holds the outcome of an augmented Dickey–Fuller test. |
| `ARIMAModel` | ARIMAModel is a fitted ARIMA(p,d,q) model: after d-fold differencing the series is modeled as ARMA(p,q). |
| `ARMAModel` | ARMAModel is a fitted autoregressive moving-average model x_t − μ = Σ Phi`i`·(x_{t−i−1} − μ) + e_t + Σ Theta`j`·e_{t−j−1}. |
| `ARModel` | ARModel is a fitted autoregressive model of a given order. |
| `Decomposition` | Decomposition holds the additive or multiplicative decomposition of a series into trend, seasonal and residual components, each aligned with the… |
| `HoltModel` | HoltModel is a fitted Holt linear-trend (double exponential smoothing) model. |
| `HoltWintersModel` | HoltWintersModel is a fitted Holt–Winters triple exponential smoothing model with additive or multiplicative seasonality of a fixed period. |
| `KPSSResult` | KPSSResult holds the outcome of a KPSS stationarity test. |
| `LinearFit` | LinearFit holds the coefficients of a straight line y = Intercept + Slope·t. |
| `MAModel` | MAModel is a fitted moving-average model x_t = μ + e_t + Σ Theta`i`·e_{t−i−1}, with innovation variance Sigma2. |

<details>
<summary><code>ADFResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ADFTest(x []float64, lags int, trend bool) (*ADFResult, error)` | ADFTest performs the augmented Dickey–Fuller test for a unit root in the series. |
| `func DickeyFuller(x []float64) (*ADFResult, error)` | DickeyFuller performs the simple (non-augmented) Dickey–Fuller test with a constant, equivalent to ADFTest with zero augmenting lags and no trend. |

</details>

<details>
<summary><code>ARIMAModel</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ARIMAFit(x []float64, p, d, q int) (*ARIMAModel, error)` | ARIMAFit fits an ARIMA(p,d,q) model by differencing the series d times and fitting an ARMA(p,q) model to the result. |
| `func (m *ARIMAModel) Forecast(h int) []float64` | Forecast returns the h-step-ahead forecasts of the ARIMA model on the original scale, integrating the differenced ARMA forecasts back up through the… |

</details>

<details>
<summary><code>ARMAModel</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ARMAFit(x []float64, p, q int) (*ARMAModel, error)` | ARMAFit fits an ARMA(p,q) model with the two-stage Hannan–Rissanen procedure: a long autoregression estimates the innovations, then x_t is… |
| `func (m *ARMAModel) Forecast(x []float64, h int) []float64` | Forecast returns the h-step-ahead forecasts of the ARMA model given the observed history x, feeding forecasts back for the AR part and dropping… |
| `func (m *ARMAModel) Residuals(x []float64) []float64` | Residuals returns the estimated innovations of the ARMA model reconstructed from the series by inverting the recursion, using zero pre-sample values. |

</details>

<details>
<summary><code>ARModel</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ARFitLeastSquares(x []float64, p int) (*ARModel, error)` | ARFitLeastSquares fits an AR(p) model by ordinary least squares regression of x_t on its p lagged values and an intercept. |
| `func BurgAR(x []float64, p int) (*ARModel, error)` | BurgAR fits an AR(p) model using Burg's method, which minimizes the sum of forward and backward prediction errors and is well conditioned for short… |
| `func YuleWalker(x []float64, p int) (*ARModel, error)` | YuleWalker fits an AR(p) model to the series by solving the Yule–Walker equations via the Levinson–Durbin recursion. |
| `func (m *ARModel) Forecast(x []float64, h int) []float64` | Forecast returns the h-step-ahead recursive forecasts of the AR model given the observed history x, feeding predicted values back in as needed. |
| `func (m *ARModel) Predict(x []float64) []float64` | Predict returns the in-sample one-step-ahead predictions of the AR model for the series x. |
| `func (m *ARModel) Residuals(x []float64) []float64` | Residuals returns the in-sample one-step-ahead residuals x_t − x̂_t of the AR model, with NaN for the first Order positions. |

</details>

<details>
<summary><code>Decomposition</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func SeasonalDecompose(x []float64, period int, mult bool) *Decomposition` | SeasonalDecompose performs a classical seasonal decomposition of the series with the given period. |

</details>

<details>
<summary><code>HoltModel</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func HoltLinear(x []float64, alpha, beta float64) (*HoltModel, error)` | HoltLinear fits Holt's linear-trend method to the series with level and trend smoothing parameters alpha and beta, both in (0,1]. |
| `func (m *HoltModel) Fitted() []float64` | Fitted returns the in-sample one-step-ahead forecasts, ŷ_t = level_{t−1} + trend_{t−1} for t ≥ 1, with the first element set to the initial… |
| `func (m *HoltModel) Forecast(h int) []float64` | Forecast returns the h-step-ahead forecasts from the fitted model: ŷ_{n+k} = level + k·trend for k = 1,…,h. |
| `func (m *HoltModel) SSE(x []float64) float64` | SSE returns the in-sample sum of squared one-step-ahead forecast errors for the observations x the model was fitted to (indices 1…n−1). |

</details>

<details>
<summary><code>HoltWintersModel</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func HoltWinters(x []float64, alpha, beta, gamma float64, period int, mult bool) (*HoltWintersModel, error)` | HoltWinters fits the Holt–Winters seasonal method to the series with the given smoothing parameters (all in (0,1]) and seasonal period. |
| `func (m *HoltWintersModel) Fitted() []float64` | Fitted returns the in-sample one-step-ahead forecasts of the Holt–Winters model. |
| `func (m *HoltWintersModel) Forecast(h int) []float64` | Forecast returns the h-step-ahead forecasts from the fitted Holt–Winters model, recycling the most recent seasonal factors. |

</details>

<details>
<summary><code>KPSSResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func KPSSTest(x []float64, lag int, trend bool) (*KPSSResult, error)` | KPSSTest performs the Kwiatkowski–Phillips–Schmidt–Shin test whose null hypothesis is that the series is (level- or trend-) stationary. |

</details>

<details>
<summary><code>LinearFit</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func FitLinearTrend(x []float64) LinearFit` | FitLinearTrend fits a straight line y = a + b·t to the series against the index t = 0,1,…,n−1 by ordinary least squares and returns the… |
| `func (f LinearFit) At(t float64) float64` | At evaluates the fitted line at index t. |

</details>

<details>
<summary><code>MAModel</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func MAFit(x []float64, q int) (*MAModel, error)` | MAFit estimates an MA(q) model using the innovations algorithm applied to the sample autocovariances of the series. |
| `func (m *MAModel) Forecast(x []float64, h int) []float64` | Forecast returns the h-step-ahead forecasts of the MA(q) model given the observed history x. |
| `func (m *MAModel) Residuals(x []float64) []float64` | Residuals returns the estimated innovations of the MA model reconstructed from the series by inverting the moving-average recursion: e_t = (x_t −… |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/timeseries).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
