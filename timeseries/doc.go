// Package timeseries is a self-contained, dependency-free time-series analysis
// toolkit written entirely with the Go standard library (only math, math/cmplx,
// sort and errors). It operates on plain []float64 sample vectors and does not
// import the parent algebra module.
//
// The package covers the classical time-series workflow end to end.
//
// # Descriptive statistics and transforms
//
// Summaries such as [Mean], [Variance], [StdDev], [Median], [Quantile],
// [Skewness] and [Kurtosis], together with reshaping and preprocessing
// operations: [Diff], [DiffOrder], [SeasonalDiff] and their inverses
// [Integrate] and [SeasonalIntegrate], [CumSum], [Lag]/[Lead]/[Shift],
// [Demean], [Standardize], [MinMaxNormalize], [BoxCox], [SimpleReturns],
// [LogReturns], [FractionalDifference] and least-squares [Detrend] via
// [FitLinearTrend].
//
// # Correlation structure
//
// [AutoCovariance], [AutoCorrelation] (ACF) and [PartialAutoCorrelation]
// (PACF via the Durbin–Levinson recursion), cross-correlation
// ([CrossCorrelation]), and diagnostic statistics [LjungBox], [BoxPierce] and
// [DurbinWatson].
//
// # Smoothing and moving averages
//
// Simple, centered, weighted, triangular and exponential moving averages
// ([MovingAverage], [WeightedMovingAverage], [ExponentialMovingAverage],
// [DoubleExponentialMovingAverage], [TripleExponentialMovingAverage]), a family
// of rolling and expanding window statistics, exponential smoothing
// ([SimpleExponentialSmoothing]), Holt's linear-trend method ([HoltLinear]) and
// Holt–Winters seasonal smoothing ([HoltWinters]).
//
// # Parametric models
//
// Autoregressive fitting by Yule–Walker ([YuleWalker]), ordinary least squares
// ([ARFitLeastSquares]) and Burg's method ([BurgAR]); moving-average estimation
// via the innovations algorithm ([MAFit]); ARMA fitting with the
// Hannan–Rissanen procedure ([ARMAFit]); and ARIMA fitting ([ARIMAFit]) with
// automatic differencing and integration. Each model provides forecasting and
// residual methods, and [ARMAToMA]/[ARMAToAR] convert between representations.
//
// # Spectral analysis and decomposition
//
// The discrete Fourier transform ([DFT]), the raw [Periodogram], AR spectral
// density ([ARSpectralDensity]), [SpectralEntropy] and dominant-cycle helpers;
// classical additive and multiplicative seasonal decomposition
// ([SeasonalDecompose]) with [SeasonallyAdjust] and [SeasonalIndices].
//
// # Stationarity and embedding
//
// Augmented Dickey–Fuller ([ADFTest]) and KPSS ([KPSSTest]) helpers,
// [VarianceRatio], number-of-differences estimation ([NumberOfDifferences]),
// delay-coordinate embedding ([Embed]), and matrix builders [LagMatrix],
// [HankelMatrix] and [ToeplitzMatrix]. Forecast-accuracy metrics such as
// [RootMeanSquaredError], [MeanAbsolutePercentageError],
// [MeanAbsoluteScaledError], [RSquared] and [TheilU] round out the package.
//
// Scalar functions return NaN on empty or otherwise undefined input rather than
// panicking, so callers can test with math.IsNaN. Model-fitting constructors
// return an error for invalid arguments.
package timeseries
