package fin

import "math"

// finMean returns the arithmetic mean of xs, or NaN for an empty slice.
func finMean(xs []float64) float64 {
	if len(xs) == 0 {
		return math.NaN()
	}
	var s float64
	for _, x := range xs {
		s += x
	}
	return s / float64(len(xs))
}

// finStdDev returns the sample (n−1) standard deviation of xs, or NaN when
// fewer than two observations are present.
func finStdDev(xs []float64) float64 {
	n := len(xs)
	if n < 2 {
		return math.NaN()
	}
	m := finMean(xs)
	var ss float64
	for _, x := range xs {
		d := x - m
		ss += d * d
	}
	return math.Sqrt(ss / float64(n-1))
}

// Volatility returns the sample standard deviation of a series of periodic
// returns, the standard measure of investment risk. It returns NaN for fewer
// than two observations. Use [Annualize] to scale it to another horizon.
func Volatility(returns []float64) float64 {
	return finStdDev(returns)
}

// Annualize scales a per-period volatility to an annual figure given the number
// of periods per year, multiplying by √periodsPerYear (the square-root-of-time
// rule). It returns NaN for periodsPerYear ≤ 0.
func Annualize(periodVol float64, periodsPerYear float64) float64 {
	if periodsPerYear <= 0 {
		return math.NaN()
	}
	return periodVol * math.Sqrt(periodsPerYear)
}

// SharpeRatio returns the Sharpe ratio of a series of periodic returns given a
// per-period riskFree rate: the mean excess return divided by the standard
// deviation of returns. It measures risk-adjusted performance and returns NaN
// for fewer than two observations or zero volatility.
func SharpeRatio(returns []float64, riskFree float64) float64 {
	sd := finStdDev(returns)
	if math.IsNaN(sd) || sd == 0 {
		return math.NaN()
	}
	return (finMean(returns) - riskFree) / sd
}

// SortinoRatio returns the Sortino ratio of a series of periodic returns given
// a per-period riskFree rate: the mean excess return divided by the downside
// deviation (the root-mean-square of shortfalls below riskFree). Unlike
// [SharpeRatio] it penalises only downside volatility. It returns NaN when
// there is no downside deviation or fewer than one observation.
func SortinoRatio(returns []float64, riskFree float64) float64 {
	if len(returns) == 0 {
		return math.NaN()
	}
	var ss float64
	for _, r := range returns {
		if r < riskFree {
			d := r - riskFree
			ss += d * d
		}
	}
	dd := math.Sqrt(ss / float64(len(returns)))
	if dd == 0 {
		return math.NaN()
	}
	return (finMean(returns) - riskFree) / dd
}

// HoldingPeriodReturn returns the total return of chaining a series of periodic
// returns (each expressed as a decimal, e.g. 0.05 for 5%): the product of
// (1+r) minus one. For an empty series it returns 0.
func HoldingPeriodReturn(returns []float64) float64 {
	prod := 1.0
	for _, r := range returns {
		prod *= 1 + r
	}
	return prod - 1
}

// AnnualizedReturn returns the geometric mean per-period return of a series
// compounded to a single period figure: (∏(1+rᵢ))^(1/n) − 1. This is the
// constant per-period rate that reproduces the series' [HoldingPeriodReturn].
// It returns NaN for an empty series.
func AnnualizedReturn(returns []float64) float64 {
	n := len(returns)
	if n == 0 {
		return math.NaN()
	}
	prod := 1.0
	for _, r := range returns {
		prod *= 1 + r
	}
	if prod <= 0 {
		return math.NaN()
	}
	return math.Pow(prod, 1/float64(n)) - 1
}

// SimpleReturn returns the simple (arithmetic) return from an opening value
// begin to a closing value end: (end − begin)/begin. It returns NaN when begin
// is zero.
func SimpleReturn(begin, end float64) float64 {
	if begin == 0 {
		return math.NaN()
	}
	return (end - begin) / begin
}

// LogReturn returns the continuously compounded (log) return from an opening
// value begin to a closing value end: ln(end/begin). It returns NaN when either
// value is non-positive.
func LogReturn(begin, end float64) float64 {
	if begin <= 0 || end <= 0 {
		return math.NaN()
	}
	return math.Log(end / begin)
}
