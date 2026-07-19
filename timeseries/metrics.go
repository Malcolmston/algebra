package timeseries

import "math"

// MeanError returns the mean forecast error (bias) actual − forecast over the
// aligned pairs. It returns NaN if the lengths differ or are zero.
func MeanError(actual, forecast []float64) float64 {
	if len(actual) != len(forecast) || len(actual) == 0 {
		return math.NaN()
	}
	var s float64
	for i := range actual {
		s += actual[i] - forecast[i]
	}
	return s / float64(len(actual))
}

// MeanAbsoluteError returns the mean of |actual − forecast| over the aligned
// pairs. It returns NaN if the lengths differ or are zero.
func MeanAbsoluteError(actual, forecast []float64) float64 {
	if len(actual) != len(forecast) || len(actual) == 0 {
		return math.NaN()
	}
	var s float64
	for i := range actual {
		s += math.Abs(actual[i] - forecast[i])
	}
	return s / float64(len(actual))
}

// MeanSquaredError returns the mean of (actual − forecast)² over the aligned
// pairs. It returns NaN if the lengths differ or are zero.
func MeanSquaredError(actual, forecast []float64) float64 {
	if len(actual) != len(forecast) || len(actual) == 0 {
		return math.NaN()
	}
	var s float64
	for i := range actual {
		d := actual[i] - forecast[i]
		s += d * d
	}
	return s / float64(len(actual))
}

// RootMeanSquaredError returns the square root of [MeanSquaredError].
func RootMeanSquaredError(actual, forecast []float64) float64 {
	return math.Sqrt(MeanSquaredError(actual, forecast))
}

// MeanAbsolutePercentageError returns the mean of |(actual − forecast)/actual|
// as a fraction (multiply by 100 for a percentage). Pairs with a zero actual
// value are skipped. It returns NaN if the lengths differ or no valid pair
// exists.
func MeanAbsolutePercentageError(actual, forecast []float64) float64 {
	if len(actual) != len(forecast) || len(actual) == 0 {
		return math.NaN()
	}
	var s float64
	var cnt int
	for i := range actual {
		if actual[i] == 0 {
			continue
		}
		s += math.Abs((actual[i] - forecast[i]) / actual[i])
		cnt++
	}
	if cnt == 0 {
		return math.NaN()
	}
	return s / float64(cnt)
}

// SymmetricMAPE returns the symmetric mean absolute percentage error, the mean
// of |actual − forecast| / ((|actual| + |forecast|)/2) as a fraction. Pairs
// where both values are zero are skipped.
func SymmetricMAPE(actual, forecast []float64) float64 {
	if len(actual) != len(forecast) || len(actual) == 0 {
		return math.NaN()
	}
	var s float64
	var cnt int
	for i := range actual {
		denom := (math.Abs(actual[i]) + math.Abs(forecast[i])) / 2
		if denom == 0 {
			continue
		}
		s += math.Abs(actual[i]-forecast[i]) / denom
		cnt++
	}
	if cnt == 0 {
		return math.NaN()
	}
	return s / float64(cnt)
}

// MedianAbsoluteError returns the median of |actual − forecast| over the
// aligned pairs, a robust accuracy measure.
func MedianAbsoluteError(actual, forecast []float64) float64 {
	if len(actual) != len(forecast) || len(actual) == 0 {
		return math.NaN()
	}
	abs := make([]float64, len(actual))
	for i := range actual {
		abs[i] = math.Abs(actual[i] - forecast[i])
	}
	return Median(abs)
}

// MeanAbsoluteScaledError returns the MASE of the forecast, the mean absolute
// error scaled by the mean absolute error of a naive one-step (seasonal-period
// m) forecast computed on the training series train. A value below 1 means the
// forecast beats the naive benchmark. It returns NaN on length mismatch or an
// undefined scale.
func MeanAbsoluteScaledError(actual, forecast, train []float64, m int) float64 {
	if len(actual) != len(forecast) || len(actual) == 0 || m < 1 || len(train) <= m {
		return math.NaN()
	}
	var scale float64
	var cnt int
	for t := m; t < len(train); t++ {
		scale += math.Abs(train[t] - train[t-m])
		cnt++
	}
	if cnt == 0 {
		return math.NaN()
	}
	scale /= float64(cnt)
	if scale == 0 {
		return math.NaN()
	}
	return MeanAbsoluteError(actual, forecast) / scale
}

// RSquared returns the coefficient of determination R² = 1 − SS_res/SS_tot of
// the forecast relative to the mean of the actual series. It returns NaN on
// length mismatch or when the actual series has zero variance.
func RSquared(actual, forecast []float64) float64 {
	if len(actual) != len(forecast) || len(actual) == 0 {
		return math.NaN()
	}
	m := mean(actual)
	var ssRes, ssTot float64
	for i := range actual {
		dr := actual[i] - forecast[i]
		dt := actual[i] - m
		ssRes += dr * dr
		ssTot += dt * dt
	}
	if ssTot == 0 {
		return math.NaN()
	}
	return 1 - ssRes/ssTot
}

// TheilU returns Theil's U2 forecast-accuracy statistic, the ratio of the RMSE
// of the forecast to the RMSE of a naive no-change forecast. Values below 1
// indicate the forecast improves on the naive benchmark. The actual and
// forecast slices are aligned one-step-ahead predictions; the naive benchmark
// uses the previous actual value, so both must have length ≥ 2. It returns NaN
// on invalid input.
func TheilU(actual, forecast []float64) float64 {
	n := len(actual)
	if n < 2 || len(forecast) != n {
		return math.NaN()
	}
	var num, den float64
	for i := 1; i < n; i++ {
		d1 := forecast[i] - actual[i]
		d2 := actual[i-1] - actual[i]
		num += d1 * d1
		den += d2 * d2
	}
	if den == 0 {
		return math.NaN()
	}
	return math.Sqrt(num / den)
}
