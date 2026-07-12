package stats

import "math"

// LinearRegression fits an ordinary-least-squares straight line y = slope·x +
// intercept to the paired data (xs, ys) and returns its slope, intercept and
// coefficient of determination R² (the fraction of variance in ys explained by
// the fit; 1 for a perfect fit).
//
// It returns NaN for all three results if the slices differ in length, have
// fewer than two elements, or xs has zero variance (a vertical line, for which
// no finite slope exists).
func LinearRegression(xs, ys []float64) (slope, intercept, r2 float64) {
	n := len(xs)
	if n != len(ys) || n < 2 {
		return math.NaN(), math.NaN(), math.NaN()
	}
	mx, my := Mean(xs), Mean(ys)
	var sxx, sxy, syy float64
	for i := 0; i < n; i++ {
		dx := xs[i] - mx
		dy := ys[i] - my
		sxx += dx * dx
		sxy += dx * dy
		syy += dy * dy
	}
	if sxx == 0 {
		return math.NaN(), math.NaN(), math.NaN()
	}
	slope = sxy / sxx
	intercept = my - slope*mx
	if syy == 0 {
		// ys is constant: the horizontal fit is exact.
		r2 = 1
	} else {
		r := sxy / math.Sqrt(sxx*syy)
		r2 = r * r
	}
	return slope, intercept, r2
}
