package spectralpde

import "math"

// ChebyshevL2Projection returns the first n+1 Chebyshev coefficients of the
// L2-orthogonal (Galerkin) projection of f with respect to the weight
// 1/sqrt(1-x^2), computed with an accurate Gauss-Chebyshev rule. The result
// is the best weighted-L2 polynomial approximation of degree n.
func ChebyshevL2Projection(f func(float64) float64, n int) []float64 {
	q := n + 8
	nodes := ChebyshevGaussNodes(q)
	fv := ApplyFunc(f, nodes)
	c := make([]float64, n+1)
	for k := 0; k <= n; k++ {
		var s float64
		for i := range nodes {
			s += fv[i] * ChebyshevT(k, nodes[i])
		}
		s *= math.Pi / float64(q) // Gauss-Chebyshev weight
		if k == 0 {
			c[k] = s / math.Pi
		} else {
			c[k] = 2 * s / math.Pi
		}
	}
	return c
}

// GalerkinProjectChebyshev is a convenience wrapper for ChebyshevL2Projection.
func GalerkinProjectChebyshev(f func(float64) float64, n int) []float64 {
	return ChebyshevL2Projection(f, n)
}

// GalerkinProjectLegendre returns the first n+1 Legendre coefficients of the
// L2-orthogonal projection of f (see LegendreProjection).
func GalerkinProjectLegendre(f func(float64) float64, n int) []float64 {
	return LegendreProjection(f, n)
}

// ChebyshevInterpolantError returns the maximum absolute difference between f
// and its degree-N Chebyshev interpolant, sampled at the given evaluation
// points.
func ChebyshevInterpolantError(f func(float64) float64, N int, eval []float64) float64 {
	coeffs := ChebyshevFit(f, N)
	var maxErr float64
	for _, x := range eval {
		e := math.Abs(f(x) - ClenshawEval(coeffs, x))
		if e > maxErr {
			maxErr = e
		}
	}
	return maxErr
}

// SpectralConvergenceRate estimates the exponential decay rate r in
// err ~ C*exp(-r*n) by a least-squares fit of log(err) against n. A larger
// positive r indicates faster spectral convergence. Non-positive errors are
// ignored.
func SpectralConvergenceRate(ns []int, errs []float64) float64 {
	var xs, ys []float64
	for i := range errs {
		if errs[i] > 0 {
			xs = append(xs, float64(ns[i]))
			ys = append(ys, math.Log(errs[i]))
		}
	}
	slope, _ := leastSquaresLine(xs, ys)
	return -slope
}

// ConvergenceOrder estimates the algebraic order p in err ~ C*h^p by a
// least-squares fit of log(err) against log(h).
func ConvergenceOrder(hs, errs []float64) float64 {
	var xs, ys []float64
	for i := range errs {
		if errs[i] > 0 && hs[i] > 0 {
			xs = append(xs, math.Log(hs[i]))
			ys = append(ys, math.Log(errs[i]))
		}
	}
	slope, _ := leastSquaresLine(xs, ys)
	return slope
}

// leastSquaresLine fits y = slope*x + intercept and returns (slope, intercept).
func leastSquaresLine(x, y []float64) (slope, intercept float64) {
	n := float64(len(x))
	if n < 2 {
		return 0, 0
	}
	var sx, sy, sxx, sxy float64
	for i := range x {
		sx += x[i]
		sy += y[i]
		sxx += x[i] * x[i]
		sxy += x[i] * y[i]
	}
	den := n*sxx - sx*sx
	if den == 0 {
		return 0, sy / n
	}
	slope = (n*sxy - sx*sy) / den
	intercept = (sy - slope*sx) / n
	return slope, intercept
}

// IsSpectrallyConverging reports whether the errors decay at least as fast as
// exp(-rate*n) over the tested resolutions, using SpectralConvergenceRate.
func IsSpectrallyConverging(ns []int, errs []float64, rate float64) bool {
	return SpectralConvergenceRate(ns, errs) >= rate
}

// CoefficientDecayRate estimates the exponential decay rate of the magnitudes
// of Chebyshev (or other) spectral coefficients, |c_k| ~ C*exp(-r*k).
func CoefficientDecayRate(coeffs []float64) float64 {
	var xs, ys []float64
	for k := range coeffs {
		a := math.Abs(coeffs[k])
		if a > 0 {
			xs = append(xs, float64(k))
			ys = append(ys, math.Log(a))
		}
	}
	slope, _ := leastSquaresLine(xs, ys)
	return -slope
}
