package timeseries

import (
	"errors"
	"math"
)

// SimpleExponentialSmoothing returns the simple exponentially smoothed series
// with smoothing parameter alpha in (0,1): s[0] = x[0] and
// s[t] = alpha·x[t] + (1−alpha)·s[t−1]. It returns nil for an out-of-range
// alpha. This is identical to [ExponentialMovingAverage].
func SimpleExponentialSmoothing(x []float64, alpha float64) []float64 {
	if alpha <= 0 || alpha >= 1 {
		if alpha == 1 {
			return copyf(x)
		}
		return nil
	}
	return ExponentialMovingAverage(x, alpha)
}

// SESForecast returns the h-step-ahead forecast of a simple exponential
// smoothing model fitted to x with parameter alpha. Because the model has no
// trend, every forecast equals the final smoothed level, so the returned slice
// of length h is constant. It returns nil for an out-of-range alpha or empty x.
func SESForecast(x []float64, alpha float64, h int) []float64 {
	s := SimpleExponentialSmoothing(x, alpha)
	if s == nil || len(s) == 0 || h < 1 {
		return nil
	}
	level := s[len(s)-1]
	out := make([]float64, h)
	for i := range out {
		out[i] = level
	}
	return out
}

// HoltModel is a fitted Holt linear-trend (double exponential smoothing) model.
type HoltModel struct {
	Alpha float64   // level smoothing parameter
	Beta  float64   // trend smoothing parameter
	Level []float64 // fitted level component, one per observation
	Trend []float64 // fitted trend component, one per observation
}

// HoltLinear fits Holt's linear-trend method to the series with level and trend
// smoothing parameters alpha and beta, both in (0,1]. The level is initialized
// to x[0] and the trend to x[1]−x[0]. It returns an error if the parameters are
// out of range or the series has fewer than two observations.
func HoltLinear(x []float64, alpha, beta float64) (*HoltModel, error) {
	if len(x) < 2 {
		return nil, errors.New("timeseries: HoltLinear needs at least 2 observations")
	}
	if alpha <= 0 || alpha > 1 || beta <= 0 || beta > 1 {
		return nil, errors.New("timeseries: HoltLinear parameters must be in (0,1]")
	}
	n := len(x)
	level := make([]float64, n)
	trend := make([]float64, n)
	level[0] = x[0]
	trend[0] = x[1] - x[0]
	for t := 1; t < n; t++ {
		prevLevel := level[t-1]
		level[t] = alpha*x[t] + (1-alpha)*(prevLevel+trend[t-1])
		trend[t] = beta*(level[t]-prevLevel) + (1-beta)*trend[t-1]
	}
	return &HoltModel{Alpha: alpha, Beta: beta, Level: level, Trend: trend}, nil
}

// Forecast returns the h-step-ahead forecasts from the fitted model:
// ŷ_{n+k} = level + k·trend for k = 1,…,h.
func (m *HoltModel) Forecast(h int) []float64 {
	if h < 1 {
		return nil
	}
	n := len(m.Level)
	l := m.Level[n-1]
	b := m.Trend[n-1]
	out := make([]float64, h)
	for k := 1; k <= h; k++ {
		out[k-1] = l + float64(k)*b
	}
	return out
}

// Fitted returns the in-sample one-step-ahead forecasts, ŷ_t = level_{t−1} +
// trend_{t−1} for t ≥ 1, with the first element set to the initial level.
func (m *HoltModel) Fitted() []float64 {
	n := len(m.Level)
	out := make([]float64, n)
	out[0] = m.Level[0]
	for t := 1; t < n; t++ {
		out[t] = m.Level[t-1] + m.Trend[t-1]
	}
	return out
}

// SSE returns the in-sample sum of squared one-step-ahead forecast errors for
// the observations x the model was fitted to (indices 1…n−1).
func (m *HoltModel) SSE(x []float64) float64 {
	fitted := m.Fitted()
	var s float64
	for t := 1; t < len(x) && t < len(fitted); t++ {
		e := x[t] - fitted[t]
		s += e * e
	}
	return s
}

// HoltWintersModel is a fitted Holt–Winters triple exponential smoothing model
// with additive or multiplicative seasonality of a fixed period.
type HoltWintersModel struct {
	Alpha          float64   // level smoothing parameter
	Beta           float64   // trend smoothing parameter
	Gamma          float64   // seasonal smoothing parameter
	Period         int       // seasonal period m
	Multiplicative bool      // seasonality type
	Level          []float64 // fitted level component
	Trend          []float64 // fitted trend component
	Season         []float64 // fitted seasonal component
}

// HoltWinters fits the Holt–Winters seasonal method to the series with the
// given smoothing parameters (all in (0,1]) and seasonal period. If mult is
// true multiplicative seasonality is used, otherwise additive. The series must
// contain at least two full seasons. It returns an error on invalid arguments.
func HoltWinters(x []float64, alpha, beta, gamma float64, period int, mult bool) (*HoltWintersModel, error) {
	if period < 2 {
		return nil, errors.New("timeseries: HoltWinters period must be >= 2")
	}
	if len(x) < 2*period {
		return nil, errors.New("timeseries: HoltWinters needs at least two full seasons")
	}
	if alpha <= 0 || alpha > 1 || beta <= 0 || beta > 1 || gamma <= 0 || gamma > 1 {
		return nil, errors.New("timeseries: HoltWinters parameters must be in (0,1]")
	}
	n := len(x)
	m := period

	// Initial level: mean of the first season.
	var l0 float64
	for i := 0; i < m; i++ {
		l0 += x[i]
	}
	l0 /= float64(m)
	// Initial trend: average per-step change between the first two seasons.
	var b0 float64
	for i := 0; i < m; i++ {
		b0 += (x[m+i] - x[i]) / float64(m)
	}
	b0 /= float64(m)
	// Initial seasonal factors from the first season.
	season := make([]float64, n)
	for i := 0; i < m; i++ {
		if mult {
			if l0 == 0 {
				season[i] = 1
			} else {
				season[i] = x[i] / l0
			}
		} else {
			season[i] = x[i] - l0
		}
	}

	level := make([]float64, n)
	trend := make([]float64, n)
	prevL, prevB := l0, b0
	for t := m; t < n; t++ {
		sPrev := season[t-m]
		var lt float64
		if mult {
			if sPrev == 0 {
				sPrev = 1e-8
			}
			lt = alpha*(x[t]/sPrev) + (1-alpha)*(prevL+prevB)
		} else {
			lt = alpha*(x[t]-sPrev) + (1-alpha)*(prevL+prevB)
		}
		bt := beta*(lt-prevL) + (1-beta)*prevB
		var st float64
		if mult {
			if lt == 0 {
				st = sPrev
			} else {
				st = gamma*(x[t]/lt) + (1-gamma)*sPrev
			}
		} else {
			st = gamma*(x[t]-lt) + (1-gamma)*sPrev
		}
		level[t] = lt
		trend[t] = bt
		season[t] = st
		prevL, prevB = lt, bt
	}
	// Fill the level and trend of the first season with the initial values so
	// forecasting and inspection see sensible numbers.
	for t := 0; t < m; t++ {
		level[t] = l0
		trend[t] = b0
	}
	return &HoltWintersModel{
		Alpha: alpha, Beta: beta, Gamma: gamma, Period: m, Multiplicative: mult,
		Level: level, Trend: trend, Season: season,
	}, nil
}

// Forecast returns the h-step-ahead forecasts from the fitted Holt–Winters
// model, recycling the most recent seasonal factors.
func (m *HoltWintersModel) Forecast(h int) []float64 {
	if h < 1 {
		return nil
	}
	n := len(m.Level)
	l := m.Level[n-1]
	b := m.Trend[n-1]
	out := make([]float64, h)
	for k := 1; k <= h; k++ {
		// Seasonal factor from the same phase one period back.
		si := n - m.Period + ((k - 1) % m.Period)
		s := m.Season[si]
		if m.Multiplicative {
			out[k-1] = (l + float64(k)*b) * s
		} else {
			out[k-1] = l + float64(k)*b + s
		}
	}
	return out
}

// Fitted returns the in-sample one-step-ahead forecasts of the Holt–Winters
// model. The first Period elements (used for initialization) are set to the
// observed values via the initial seasonal factors.
func (m *HoltWintersModel) Fitted() []float64 {
	n := len(m.Level)
	out := make([]float64, n)
	for t := 0; t < n; t++ {
		if t < m.Period {
			out[t] = math.NaN()
			continue
		}
		l := m.Level[t-1]
		b := m.Trend[t-1]
		s := m.Season[t-m.Period]
		if m.Multiplicative {
			out[t] = (l + b) * s
		} else {
			out[t] = l + b + s
		}
	}
	return out
}

// BrownDoubleExponential returns Brown's single-parameter double exponential
// smoothing (linear-trend) forecast for the series. It computes the singly and
// doubly smoothed statistics with parameter alpha, forms the level and trend at
// the final point, and returns the h-step-ahead forecasts. It returns nil for
// an out-of-range alpha or empty x.
func BrownDoubleExponential(x []float64, alpha float64, h int) []float64 {
	if alpha <= 0 || alpha >= 1 || len(x) == 0 || h < 1 {
		return nil
	}
	s1 := ExponentialMovingAverage(x, alpha)
	s2 := ExponentialMovingAverage(s1, alpha)
	n := len(x)
	a := 2*s1[n-1] - s2[n-1]
	b := (alpha / (1 - alpha)) * (s1[n-1] - s2[n-1])
	out := make([]float64, h)
	for k := 1; k <= h; k++ {
		out[k-1] = a + float64(k)*b
	}
	return out
}
