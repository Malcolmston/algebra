package timeseries

import "math"

// Decomposition holds the additive or multiplicative decomposition of a series
// into trend, seasonal and residual components, each aligned with the original
// observations. Positions where the trend is undefined (near the boundaries)
// hold NaN.
type Decomposition struct {
	Observed       []float64
	Trend          []float64
	Seasonal       []float64
	Residual       []float64
	Period         int
	Multiplicative bool
}

// symmetricMA computes the trend by a centered moving average of length equal
// to the seasonal period. For an even period a 2×period moving average (a
// weighted average with half weights at the ends) is used so the result stays
// centered; boundary positions that lack a full window are set to NaN.
func symmetricMA(x []float64, period int) []float64 {
	n := len(x)
	trend := make([]float64, n)
	for i := range trend {
		trend[i] = math.NaN()
	}
	if period < 1 {
		return trend
	}
	if period%2 == 0 {
		half := period / 2
		for i := half; i < n-half; i++ {
			var s float64
			s += 0.5 * x[i-half]
			s += 0.5 * x[i+half]
			for j := -half + 1; j <= half-1; j++ {
				s += x[i+j]
			}
			trend[i] = s / float64(period)
		}
	} else {
		half := period / 2
		for i := half; i < n-half; i++ {
			var s float64
			for j := -half; j <= half; j++ {
				s += x[i+j]
			}
			trend[i] = s / float64(period)
		}
	}
	return trend
}

// SeasonalDecompose performs a classical seasonal decomposition of the series
// with the given period. If mult is true a multiplicative model
// (x = trend·seasonal·residual) is used, otherwise an additive model
// (x = trend + seasonal + residual). The seasonal component is the average
// detrended value per phase, normalized to sum to zero (additive) or to average
// one (multiplicative), and repeated across the series. It returns nil if the
// period is invalid or the series is shorter than two periods.
func SeasonalDecompose(x []float64, period int, mult bool) *Decomposition {
	n := len(x)
	if period < 2 || n < 2*period {
		return nil
	}
	trend := symmetricMA(x, period)
	// Detrend.
	detr := make([]float64, n)
	for i := 0; i < n; i++ {
		if math.IsNaN(trend[i]) {
			detr[i] = math.NaN()
			continue
		}
		if mult {
			if trend[i] == 0 {
				detr[i] = math.NaN()
			} else {
				detr[i] = x[i] / trend[i]
			}
		} else {
			detr[i] = x[i] - trend[i]
		}
	}
	// Average per phase.
	sums := make([]float64, period)
	counts := make([]int, period)
	for i := 0; i < n; i++ {
		if math.IsNaN(detr[i]) {
			continue
		}
		ph := i % period
		sums[ph] += detr[i]
		counts[ph]++
	}
	seasonalPhase := make([]float64, period)
	for p := 0; p < period; p++ {
		if counts[p] > 0 {
			seasonalPhase[p] = sums[p] / float64(counts[p])
		} else if mult {
			seasonalPhase[p] = 1
		}
	}
	// Normalize.
	if mult {
		var mp float64
		for _, v := range seasonalPhase {
			mp += v
		}
		mp /= float64(period)
		if mp != 0 {
			for p := range seasonalPhase {
				seasonalPhase[p] /= mp
			}
		}
	} else {
		var mp float64
		for _, v := range seasonalPhase {
			mp += v
		}
		mp /= float64(period)
		for p := range seasonalPhase {
			seasonalPhase[p] -= mp
		}
	}
	seasonal := make([]float64, n)
	residual := make([]float64, n)
	for i := 0; i < n; i++ {
		seasonal[i] = seasonalPhase[i%period]
		if math.IsNaN(trend[i]) {
			residual[i] = math.NaN()
			continue
		}
		if mult {
			if trend[i]*seasonal[i] == 0 {
				residual[i] = math.NaN()
			} else {
				residual[i] = x[i] / (trend[i] * seasonal[i])
			}
		} else {
			residual[i] = x[i] - trend[i] - seasonal[i]
		}
	}
	return &Decomposition{
		Observed:       copyf(x),
		Trend:          trend,
		Seasonal:       seasonal,
		Residual:       residual,
		Period:         period,
		Multiplicative: mult,
	}
}

// SeasonalIndices returns the per-phase seasonal factors of the series for the
// given period (length period): additive deviations that sum to zero when mult
// is false, or multiplicative factors that average to one when mult is true. It
// returns nil for an invalid period or too-short series.
func SeasonalIndices(x []float64, period int, mult bool) []float64 {
	d := SeasonalDecompose(x, period, mult)
	if d == nil {
		return nil
	}
	out := make([]float64, period)
	copy(out, d.Seasonal[:period])
	return out
}

// SeasonallyAdjust removes the estimated seasonal component from the series,
// returning the seasonally adjusted series (trend plus residual for the
// additive model, or the series divided by the seasonal factors for the
// multiplicative model). It returns nil for an invalid period or too-short
// series.
func SeasonallyAdjust(x []float64, period int, mult bool) []float64 {
	d := SeasonalDecompose(x, period, mult)
	if d == nil {
		return nil
	}
	out := make([]float64, len(x))
	for i := range x {
		if mult {
			if d.Seasonal[i] == 0 {
				out[i] = math.NaN()
			} else {
				out[i] = x[i] / d.Seasonal[i]
			}
		} else {
			out[i] = x[i] - d.Seasonal[i]
		}
	}
	return out
}

// TrendComponent returns the classical moving-average trend estimate of the
// series for the given period, with NaN at boundary positions. It returns nil
// if the period is invalid.
func TrendComponent(x []float64, period int) []float64 {
	if period < 1 {
		return nil
	}
	return symmetricMA(x, period)
}
