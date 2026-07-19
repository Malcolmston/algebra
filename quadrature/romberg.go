package quadrature

import "math"

// RombergTable computes the full Romberg extrapolation table for the integral
// of f over [a, b] using k levels of trapezoidal refinement (rows 0..k-1). The
// returned lower-triangular table has T[i][j] with T[i][0] the composite
// trapezoid estimate on 2^i subintervals and T[i][j] the j-th Richardson
// extrapolation. The best estimate is T[k-1][k-1].
func RombergTable(f Func, a, b float64, k int) [][]float64 {
	if k < 1 {
		k = 1
	}
	T := make([][]float64, k)
	for i := range T {
		T[i] = make([]float64, k)
	}
	h := b - a
	T[0][0] = 0.5 * h * (f(a) + f(b))
	for i := 1; i < k; i++ {
		h *= 0.5
		// refine trapezoid: add the new midpoints
		var sum float64
		n := 1 << (i - 1)
		for j := 0; j < n; j++ {
			sum += f(a + (2*float64(j)+1)*h)
		}
		T[i][0] = 0.5*T[i-1][0] + h*sum
		pow4 := 1.0
		for j := 1; j <= i; j++ {
			pow4 *= 4
			T[i][j] = (pow4*T[i][j-1] - T[i-1][j-1]) / (pow4 - 1)
		}
	}
	return T
}

// Romberg approximates the integral of f over [a, b] using Romberg
// extrapolation, stopping when two successive diagonal estimates agree to the
// given absolute-or-relative tolerance or after maxLevels refinements. It
// returns the best estimate and an error estimate.
func Romberg(f Func, a, b, tol float64, maxLevels int) Result {
	if maxLevels < 1 {
		maxLevels = 1
	}
	if tol <= 0 {
		tol = 1e-12
	}
	prevRow := make([]float64, maxLevels+1)
	curRow := make([]float64, maxLevels+1)
	h := b - a
	prevRow[0] = 0.5 * h * (f(a) + f(b))
	evals := 2
	best := prevRow[0]
	var errEst float64
	for i := 1; i <= maxLevels; i++ {
		h *= 0.5
		var sum float64
		n := 1 << (i - 1)
		for j := 0; j < n; j++ {
			sum += f(a + (2*float64(j)+1)*h)
		}
		evals += n
		curRow[0] = 0.5*prevRow[0] + h*sum
		pow4 := 1.0
		for j := 1; j <= i; j++ {
			pow4 *= 4
			curRow[j] = (pow4*curRow[j-1] - prevRow[j-1]) / (pow4 - 1)
		}
		errEst = math.Abs(curRow[i] - prevRow[i-1])
		best = curRow[i]
		if errEst <= tol*math.Abs(best)+tol {
			return Result{Value: best, AbsErr: errEst, Evals: evals, Success: true}
		}
		prevRow, curRow = curRow, prevRow
	}
	return Result{Value: best, AbsErr: errEst, Evals: evals, Success: false}
}

// RichardsonExtrapolate performs one step of Richardson extrapolation on two
// estimates A(h) and A(h/ratio) of the same quantity whose leading error term
// is proportional to h^order, returning the higher-order combination
// (ratio^order * fine - coarse) / (ratio^order - 1).
func RichardsonExtrapolate(coarse, fine, ratio, order float64) float64 {
	p := math.Pow(ratio, order)
	return (p*fine - coarse) / (p - 1)
}

// adaptiveSimpson performs the recursive step of adaptive Simpson integration.
func adaptiveSimpson(f Func, a, b, fa, fb, fm, whole, tol float64, depth int, evals *int) float64 {
	m := 0.5 * (a + b)
	lm := 0.5 * (a + m)
	rm := 0.5 * (m + b)
	flm := f(lm)
	frm := f(rm)
	*evals += 2
	left := (m - a) / 6 * (fa + 4*flm + fm)
	right := (b - m) / 6 * (fm + 4*frm + fb)
	if depth <= 0 || math.Abs(left+right-whole) <= 15*tol {
		return left + right + (left+right-whole)/15
	}
	return adaptiveSimpson(f, a, m, fa, fm, flm, left, tol/2, depth-1, evals) +
		adaptiveSimpson(f, m, b, fm, fb, frm, right, tol/2, depth-1, evals)
}

// AdaptiveSimpson approximates the integral of f over [a, b] to the requested
// absolute tolerance using recursive adaptive Simpson quadrature with error
// control based on comparing a whole-panel Simpson estimate against the sum of
// its two halves.
func AdaptiveSimpson(f Func, a, b, tol float64) float64 {
	return AdaptiveSimpsonResult(f, a, b, tol).Value
}

// AdaptiveSimpsonResult is like AdaptiveSimpson but returns a Result carrying
// the error estimate and evaluation count.
func AdaptiveSimpsonResult(f Func, a, b, tol float64) Result {
	if a == b {
		return Result{Success: true}
	}
	if tol <= 0 {
		tol = 1e-10
	}
	fa := f(a)
	fb := f(b)
	m := 0.5 * (a + b)
	fm := f(m)
	evals := 3
	whole := (b - a) / 6 * (fa + 4*fm + fb)
	v := adaptiveSimpson(f, a, b, fa, fb, fm, whole, tol, 50, &evals)
	return Result{Value: v, AbsErr: math.Abs(v - whole), Evals: evals, Success: true}
}

// AdaptiveTrapezoid approximates the integral of f over [a, b] to the
// requested absolute tolerance using recursive adaptive trapezoidal
// refinement.
func AdaptiveTrapezoid(f Func, a, b, tol float64) float64 {
	if a == b {
		return 0
	}
	if tol <= 0 {
		tol = 1e-10
	}
	var rec func(a, b, fa, fb, whole, tol float64, depth int) float64
	rec = func(a, b, fa, fb, whole, tol float64, depth int) float64 {
		m := 0.5 * (a + b)
		fm := f(m)
		left := 0.25 * (b - a) * (fa + fm)
		right := 0.25 * (b - a) * (fm + fb)
		if depth <= 0 || math.Abs(left+right-whole) <= 3*tol {
			return left + right + (left+right-whole)/3
		}
		return rec(a, m, fa, fm, left, tol/2, depth-1) +
			rec(m, b, fm, fb, right, tol/2, depth-1)
	}
	fa := f(a)
	fb := f(b)
	whole := 0.5 * (b - a) * (fa + fb)
	return rec(a, b, fa, fb, whole, tol, 50)
}
