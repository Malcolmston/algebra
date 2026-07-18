// Package numint implements numerical integration (quadrature) of functions
// of one and several variables using only the Go standard library.
//
// The package collects the classical fixed-order Newton-Cotes rules
// (trapezoid, midpoint, Simpson 1/3, Simpson 3/8, Boole, Weddle), their
// composite forms, and higher-accuracy schemes: Romberg extrapolation,
// adaptive Simpson and adaptive Gauss-Kronrod, fixed Gauss-Legendre rules
// from 2 to 10 points plus an n-point routine for arbitrary order,
// Gauss-Lobatto rules, Clenshaw-Curtis and Gauss-Chebyshev quadrature, and a
// tanh-sinh (double-exponential) rule for integrands with endpoint
// singularities. Sample-based helpers integrate tabulated data, including
// cumulative integrals, and tensor-product routines extend the fixed rules to
// two and three dimensions. A small deterministic Monte-Carlo family covers
// higher-dimensional volumes.
//
// All routines are deterministic and depend only on the standard library.
package numint

import (
	"math"
	"math/rand"
)

// Func is a real-valued function of a single real variable.
type Func func(x float64) float64

// Func2 is a real-valued function of two real variables.
type Func2 func(x, y float64) float64

// Func3 is a real-valued function of three real variables.
type Func3 func(x, y, z float64) float64

// FuncND is a real-valued function of a vector of real variables.
type FuncND func(x []float64) float64

// QuadResult bundles the estimated value of an integral with an estimate of
// the absolute error and the number of integrand evaluations performed.
type QuadResult struct {
	// Value is the estimated value of the integral.
	Value float64
	// ErrEst is an estimate of the absolute error of Value.
	ErrEst float64
	// Evals is the number of integrand evaluations performed.
	Evals int
}

// --- small unexported helpers ---------------------------------------------

// numintabs returns the absolute value of x.
func numintabs(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}

// numintevenUp returns n unchanged when it is a positive even number, and n+1
// when n is odd; it returns 2 for n <= 0. Several rules require an even number
// of panels.
func numintevenUp(n int) int {
	if n <= 0 {
		return 2
	}
	if n%2 != 0 {
		return n + 1
	}
	return n
}

// numintmultipleUp rounds n up to the smallest positive multiple of m.
func numintmultipleUp(n, m int) int {
	if n <= 0 {
		return m
	}
	if r := n % m; r != 0 {
		return n + (m - r)
	}
	return n
}

// numinttrapWeights returns the abscissae and weights of the composite
// trapezoid rule with n panels on [a, b].
func numinttrapWeights(a, b float64, n int) (xs, ws []float64) {
	if n < 1 {
		n = 1
	}
	h := (b - a) / float64(n)
	xs = make([]float64, n+1)
	ws = make([]float64, n+1)
	for i := 0; i <= n; i++ {
		xs[i] = a + float64(i)*h
		if i == 0 || i == n {
			ws[i] = h / 2
		} else {
			ws[i] = h
		}
	}
	return xs, ws
}

// numintsimpWeights returns the abscissae and weights of the composite
// Simpson 1/3 rule with n panels (n even) on [a, b].
func numintsimpWeights(a, b float64, n int) (xs, ws []float64) {
	n = numintevenUp(n)
	h := (b - a) / float64(n)
	xs = make([]float64, n+1)
	ws = make([]float64, n+1)
	for i := 0; i <= n; i++ {
		xs[i] = a + float64(i)*h
		switch {
		case i == 0 || i == n:
			ws[i] = h / 3
		case i%2 == 1:
			ws[i] = 4 * h / 3
		default:
			ws[i] = 2 * h / 3
		}
	}
	return xs, ws
}

// numintmidWeights returns the abscissae and weights of the composite
// midpoint rule with n panels on [a, b].
func numintmidWeights(a, b float64, n int) (xs, ws []float64) {
	if n < 1 {
		n = 1
	}
	h := (b - a) / float64(n)
	xs = make([]float64, n)
	ws = make([]float64, n)
	for i := 0; i < n; i++ {
		xs[i] = a + (float64(i)+0.5)*h
		ws[i] = h
	}
	return xs, ws
}

// numintsimpsonNonuniform integrates the quadratic through the three points
// (x0,y0), (x1,y1), (x2,y2) over [x0, x2] with possibly unequal spacing.
func numintsimpsonNonuniform(x0, x1, x2, y0, y1, y2 float64) float64 {
	h0 := x1 - x0
	h1 := x2 - x1
	if h0 == 0 || h1 == 0 {
		return 0.5 * (x2 - x0) * (y0 + y2)
	}
	return (h0 + h1) / 6 *
		((2-h1/h0)*y0 + (h0+h1)*(h0+h1)/(h0*h1)*y1 + (2-h0/h1)*y2)
}

// --- single-panel Newton-Cotes rules --------------------------------------

// TrapezoidRule applies the two-point trapezoid rule to f over [a, b].
func TrapezoidRule(f Func, a, b float64) float64 {
	return 0.5 * (b - a) * (f(a) + f(b))
}

// MidpointRule applies the one-point open midpoint rule to f over [a, b].
func MidpointRule(f Func, a, b float64) float64 {
	return (b - a) * f(0.5*(a+b))
}

// SimpsonRule applies the three-point Simpson 1/3 rule to f over [a, b].
func SimpsonRule(f Func, a, b float64) float64 {
	m := 0.5 * (a + b)
	return (b - a) / 6 * (f(a) + 4*f(m) + f(b))
}

// Simpson38Rule applies the four-point Simpson 3/8 rule to f over [a, b].
func Simpson38Rule(f Func, a, b float64) float64 {
	h := (b - a) / 3
	return 3 * h / 8 * (f(a) + 3*f(a+h) + 3*f(a+2*h) + f(b))
}

// BooleRule applies the five-point Boole rule to f over [a, b].
func BooleRule(f Func, a, b float64) float64 {
	h := (b - a) / 4
	return 2 * h / 45 * (7*f(a) + 32*f(a+h) + 12*f(a+2*h) + 32*f(a+3*h) + 7*f(b))
}

// WeddleRule applies the seven-point Weddle (Newton-Cotes degree six) rule to
// f over [a, b].
func WeddleRule(f Func, a, b float64) float64 {
	h := (b - a) / 6
	return 3 * h / 10 * (f(a) + 5*f(a+h) + f(a+2*h) + 6*f(a+3*h) +
		f(a+4*h) + 5*f(a+5*h) + f(b))
}

// --- composite Newton-Cotes rules -----------------------------------------

// Trapezoid approximates the integral of f over [a, b] using the composite
// trapezoid rule with n panels.
func Trapezoid(f Func, a, b float64, n int) float64 {
	if n < 1 {
		n = 1
	}
	h := (b - a) / float64(n)
	sum := 0.5 * (f(a) + f(b))
	for i := 1; i < n; i++ {
		sum += f(a + float64(i)*h)
	}
	return sum * h
}

// Midpoint approximates the integral of f over [a, b] using the composite
// midpoint rule with n panels.
func Midpoint(f Func, a, b float64, n int) float64 {
	if n < 1 {
		n = 1
	}
	h := (b - a) / float64(n)
	sum := 0.0
	for i := 0; i < n; i++ {
		sum += f(a + (float64(i)+0.5)*h)
	}
	return sum * h
}

// Simpson approximates the integral of f over [a, b] using the composite
// Simpson 1/3 rule; n is the number of panels and is rounded up to an even
// number when odd.
func Simpson(f Func, a, b float64, n int) float64 {
	n = numintevenUp(n)
	h := (b - a) / float64(n)
	sum := f(a) + f(b)
	for i := 1; i < n; i++ {
		x := a + float64(i)*h
		if i%2 == 1 {
			sum += 4 * f(x)
		} else {
			sum += 2 * f(x)
		}
	}
	return sum * h / 3
}

// Simpson38 approximates the integral of f over [a, b] using the composite
// Simpson 3/8 rule; n is the number of panels and is rounded up to a multiple
// of three when necessary.
func Simpson38(f Func, a, b float64, n int) float64 {
	n = numintmultipleUp(n, 3)
	h := (b - a) / float64(n)
	sum := f(a) + f(b)
	for i := 1; i < n; i++ {
		x := a + float64(i)*h
		if i%3 == 0 {
			sum += 2 * f(x)
		} else {
			sum += 3 * f(x)
		}
	}
	return sum * 3 * h / 8
}

// Boole approximates the integral of f over [a, b] using the composite Boole
// rule; n is the number of panels and is rounded up to a multiple of four.
func Boole(f Func, a, b float64, n int) float64 {
	n = numintmultipleUp(n, 4)
	h := (b - a) / float64(n)
	sum := 0.0
	for i := 0; i < n; i += 4 {
		x := a + float64(i)*h
		sum += 7*f(x) + 32*f(x+h) + 12*f(x+2*h) + 32*f(x+3*h) + 7*f(x+4*h)
	}
	return sum * 2 * h / 45
}

// Weddle approximates the integral of f over [a, b] using the composite Weddle
// rule; n is the number of panels and is rounded up to a multiple of six.
func Weddle(f Func, a, b float64, n int) float64 {
	n = numintmultipleUp(n, 6)
	h := (b - a) / float64(n)
	sum := 0.0
	for i := 0; i < n; i += 6 {
		x := a + float64(i)*h
		sum += f(x) + 5*f(x+h) + f(x+2*h) + 6*f(x+3*h) +
			f(x+4*h) + 5*f(x+5*h) + f(x+6*h)
	}
	return sum * 3 * h / 10
}

// --- sample- and data-based integration -----------------------------------

// TrapezoidSamples integrates equally spaced samples ys with spacing h using
// the composite trapezoid rule.
func TrapezoidSamples(ys []float64, h float64) float64 {
	n := len(ys)
	if n < 2 {
		return 0
	}
	sum := 0.5 * (ys[0] + ys[n-1])
	for i := 1; i < n-1; i++ {
		sum += ys[i]
	}
	return sum * h
}

// SimpsonSamples integrates equally spaced samples ys with spacing h using the
// composite Simpson 1/3 rule. It requires an odd number of samples (an even
// number of panels); when given an even number of samples it integrates the
// leading odd-length block with Simpson and closes the final panel with the
// trapezoid rule.
func SimpsonSamples(ys []float64, h float64) float64 {
	n := len(ys)
	if n < 2 {
		return 0
	}
	if n == 2 {
		return 0.5 * h * (ys[0] + ys[1])
	}
	last := n - 1
	extra := 0.0
	if (n-1)%2 != 0 {
		// Even number of panels is required; peel off the final panel.
		extra = 0.5 * h * (ys[n-2] + ys[n-1])
		last = n - 2
	}
	sum := ys[0] + ys[last]
	for i := 1; i < last; i++ {
		if i%2 == 1 {
			sum += 4 * ys[i]
		} else {
			sum += 2 * ys[i]
		}
	}
	return sum*h/3 + extra
}

// MidpointSamples integrates samples ys taken at the midpoints of n panels of
// width h using the composite midpoint rule.
func MidpointSamples(ys []float64, h float64) float64 {
	sum := 0.0
	for _, y := range ys {
		sum += y
	}
	return sum * h
}

// TrapezoidData integrates tabulated data (xs, ys) with arbitrary, not
// necessarily uniform, spacing using the trapezoid rule. xs must be sorted and
// the same length as ys.
func TrapezoidData(xs, ys []float64) float64 {
	n := len(xs)
	if n < 2 || len(ys) != n {
		return 0
	}
	sum := 0.0
	for i := 1; i < n; i++ {
		sum += 0.5 * (xs[i] - xs[i-1]) * (ys[i] + ys[i-1])
	}
	return sum
}

// CumulativeTrapezoid returns the cumulative integral of tabulated data
// (xs, ys) evaluated at each xs. The returned slice has the same length as xs;
// its first element is zero and element i holds the integral from xs[0] to
// xs[i].
func CumulativeTrapezoid(xs, ys []float64) []float64 {
	n := len(xs)
	res := make([]float64, n)
	if n < 2 || len(ys) != n {
		return res
	}
	for i := 1; i < n; i++ {
		res[i] = res[i-1] + 0.5*(xs[i]-xs[i-1])*(ys[i]+ys[i-1])
	}
	return res
}

// CumulativeTrapezoidUniform returns the cumulative integral of equally spaced
// samples ys with spacing h. The returned slice has the same length as ys and
// its first element is zero.
func CumulativeTrapezoidUniform(ys []float64, h float64) []float64 {
	n := len(ys)
	res := make([]float64, n)
	for i := 1; i < n; i++ {
		res[i] = res[i-1] + 0.5*h*(ys[i]+ys[i-1])
	}
	return res
}

// CumulativeSimpson returns the cumulative integral of tabulated data
// (xs, ys). The returned slice has the same length as xs and its first element
// is zero. Interior points use a local (possibly non-uniform) Simpson rule;
// the first panel uses the trapezoid rule.
func CumulativeSimpson(xs, ys []float64) []float64 {
	n := len(xs)
	res := make([]float64, n)
	if n < 2 || len(ys) != n {
		return res
	}
	res[1] = res[0] + 0.5*(xs[1]-xs[0])*(ys[0]+ys[1])
	for i := 2; i < n; i++ {
		res[i] = res[i-2] +
			numintsimpsonNonuniform(xs[i-2], xs[i-1], xs[i], ys[i-2], ys[i-1], ys[i])
	}
	return res
}

// --- Romberg integration --------------------------------------------------

// Romberg approximates the integral of f over [a, b] using Romberg
// integration with maxLevel levels of Richardson extrapolation on the
// composite trapezoid rule. A typical value for maxLevel is 5 to 8.
func Romberg(f Func, a, b float64, maxLevel int) float64 {
	v, _ := RombergTable(f, a, b, maxLevel)
	return v
}

// RombergTol approximates the integral of f over [a, b] using Romberg
// integration, stopping once two successive diagonal estimates differ by less
// than tol in absolute value or maxLevel levels have been used. It returns the
// best available estimate.
func RombergTol(f Func, a, b, tol float64, maxLevel int) float64 {
	if maxLevel < 1 {
		maxLevel = 1
	}
	r := make([]float64, maxLevel)
	prevDiag := 0.0
	h := b - a
	r[0] = 0.5 * h * (f(a) + f(b))
	for i := 1; i < maxLevel; i++ {
		h *= 0.5
		sum := 0.0
		n := 1 << (i - 1)
		for k := 0; k < n; k++ {
			sum += f(a + (float64(2*k)+1)*h)
		}
		// In-place Richardson extrapolation of the current trapezoid row.
		prev := r[0] // R(i-1,0)
		r[0] = 0.5*prev + h*sum
		pow := 1.0
		for j := 1; j <= i; j++ {
			pow *= 4
			temp := r[j] // R(i-1,j)
			r[j] = (pow*r[j-1] - prev) / (pow - 1)
			prev = temp
		}
		if numintabs(r[i]-prevDiag) < tol {
			return r[i]
		}
		prevDiag = r[i]
	}
	return r[maxLevel-1]
}

// RombergTable performs Romberg integration of f over [a, b] with maxLevel
// levels and returns the best estimate (the last diagonal entry) together with
// the full lower-triangular Romberg table, indexed table[i][j].
func RombergTable(f Func, a, b float64, maxLevel int) (float64, [][]float64) {
	if maxLevel < 1 {
		maxLevel = 1
	}
	table := make([][]float64, maxLevel)
	for i := range table {
		table[i] = make([]float64, i+1)
	}
	h := b - a
	table[0][0] = 0.5 * h * (f(a) + f(b))
	for i := 1; i < maxLevel; i++ {
		h *= 0.5
		sum := 0.0
		n := 1 << (i - 1)
		for k := 0; k < n; k++ {
			sum += f(a + (float64(2*k)+1)*h)
		}
		table[i][0] = 0.5*table[i-1][0] + h*sum
		pow := 1.0
		for j := 1; j <= i; j++ {
			pow *= 4
			table[i][j] = (pow*table[i][j-1] - table[i-1][j-1]) / (pow - 1)
		}
	}
	return table[maxLevel-1][maxLevel-1], table
}

// --- adaptive rules -------------------------------------------------------

// AdaptiveTrapezoid approximates the integral of f over [a, b] by recursively
// bisecting subintervals until the trapezoid estimate on each meets the
// requested absolute tolerance tol.
func AdaptiveTrapezoid(f Func, a, b, tol float64) float64 {
	fa, fb := f(a), f(b)
	whole := 0.5 * (b - a) * (fa + fb)
	return numintadaptTrap(f, a, b, fa, fb, tol, whole, 50)
}

// numintadaptTrap is the recursive core of AdaptiveTrapezoid.
func numintadaptTrap(f Func, a, b, fa, fb, tol, whole float64, depth int) float64 {
	m := 0.5 * (a + b)
	fm := f(m)
	left := 0.25 * (b - a) * (fa + fm)
	right := 0.25 * (b - a) * (fm + fb)
	if depth <= 0 || numintabs(left+right-whole) <= 3*tol {
		return left + right + (left+right-whole)/3
	}
	return numintadaptTrap(f, a, m, fa, fm, tol/2, left, depth-1) +
		numintadaptTrap(f, m, b, fm, fb, tol/2, right, depth-1)
}

// AdaptiveSimpson approximates the integral of f over [a, b] using adaptive
// Simpson quadrature with an absolute tolerance tol.
func AdaptiveSimpson(f Func, a, b, tol float64) float64 {
	return AdaptiveSimpsonMaxDepth(f, a, b, tol, 50)
}

// AdaptiveSimpsonMaxDepth approximates the integral of f over [a, b] using
// adaptive Simpson quadrature with an absolute tolerance tol and a bound
// maxDepth on the recursion depth.
func AdaptiveSimpsonMaxDepth(f Func, a, b, tol float64, maxDepth int) float64 {
	fa := f(a)
	fb := f(b)
	m := 0.5 * (a + b)
	fm := f(m)
	whole := (b - a) / 6 * (fa + 4*fm + fb)
	return numintadaptSimpson(f, a, b, fa, fb, fm, whole, tol, maxDepth)
}

// numintadaptSimpson is the recursive core of the adaptive Simpson rule.
func numintadaptSimpson(f Func, a, b, fa, fb, fm, whole, tol float64, depth int) float64 {
	m := 0.5 * (a + b)
	lm := 0.5 * (a + m)
	rm := 0.5 * (m + b)
	flm := f(lm)
	frm := f(rm)
	left := (m - a) / 6 * (fa + 4*flm + fm)
	right := (b - m) / 6 * (fm + 4*frm + fb)
	if depth <= 0 || numintabs(left+right-whole) <= 15*tol {
		return left + right + (left+right-whole)/15
	}
	return numintadaptSimpson(f, a, m, fa, fm, flm, left, tol/2, depth-1) +
		numintadaptSimpson(f, m, b, fm, fb, frm, right, tol/2, depth-1)
}

// AdaptiveGaussKronrod approximates the integral of f over [a, b] using
// globally adaptive 15-point Gauss-Kronrod quadrature with an absolute
// tolerance tol. The returned QuadResult carries the value, an error estimate,
// and the number of integrand evaluations.
func AdaptiveGaussKronrod(f Func, a, b, tol float64) QuadResult {
	v, e := GaussKronrod15(f, a, b)
	val, err, ev := numintadaptGK(f, a, b, tol, v, e, 50)
	return QuadResult{Value: val, ErrEst: err, Evals: ev + 15}
}

// numintadaptGK is the recursive core of AdaptiveGaussKronrod.
func numintadaptGK(f Func, a, b, tol, whole, wholeErr float64, depth int) (float64, float64, int) {
	if depth <= 0 || wholeErr <= tol {
		return whole, wholeErr, 0
	}
	m := 0.5 * (a + b)
	lv, le := GaussKronrod15(f, a, m)
	rv, re := GaussKronrod15(f, m, b)
	if le+re <= tol {
		return lv + rv, le + re, 30
	}
	l, lerr, ln := numintadaptGK(f, a, m, tol/2, lv, le, depth-1)
	r, rerr, rn := numintadaptGK(f, m, b, tol/2, rv, re, depth-1)
	return l + r, lerr + rerr, ln + rn + 30
}

// --- fixed Gauss-Legendre rules -------------------------------------------

// numintglNodes holds the Gauss-Legendre abscissae on [-1, 1] for orders 2..10.
var numintglNodes = map[int][]float64{
	2:  {-0.5773502691896257, 0.5773502691896257},
	3:  {-0.7745966692414834, 0, 0.7745966692414834},
	4:  {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526},
	5:  {-0.9061798459386640, -0.5384693101056831, 0, 0.5384693101056831, 0.9061798459386640},
	6:  {-0.9324695142031521, -0.6612093864662645, -0.2386191860831969, 0.2386191860831969, 0.6612093864662645, 0.9324695142031521},
	7:  {-0.9491079123427585, -0.7415311855993945, -0.4058451513773972, 0, 0.4058451513773972, 0.7415311855993945, 0.9491079123427585},
	8:  {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498, 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363},
	9:  {-0.9681602395076261, -0.8360311073266358, -0.6133714327005904, -0.3242534234038089, 0, 0.3242534234038089, 0.6133714327005904, 0.8360311073266358, 0.9681602395076261},
	10: {-0.9739065285171717, -0.8650633666889845, -0.6794095682990244, -0.4333953941292472, -0.1488743389816312, 0.1488743389816312, 0.4333953941292472, 0.6794095682990244, 0.8650633666889845, 0.9739065285171717},
}

// numintglWeights holds the Gauss-Legendre weights on [-1, 1] for orders 2..10.
var numintglWeights = map[int][]float64{
	2:  {1, 1},
	3:  {0.5555555555555556, 0.8888888888888888, 0.5555555555555556},
	4:  {0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538},
	5:  {0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891},
	6:  {0.1713244923791704, 0.3607615730481386, 0.4679139345726910, 0.4679139345726910, 0.3607615730481386, 0.1713244923791704},
	7:  {0.1294849661688697, 0.2797053914892766, 0.3818300505051189, 0.4179591836734694, 0.3818300505051189, 0.2797053914892766, 0.1294849661688697},
	8:  {0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620, 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763},
	9:  {0.0812743883615744, 0.1806481606948574, 0.2606106964029354, 0.3123470770400029, 0.3302393550012598, 0.3123470770400029, 0.2606106964029354, 0.1806481606948574, 0.0812743883615744},
	10: {0.0666713443086881, 0.1494513491505806, 0.2190863625159820, 0.2692667193099963, 0.2955242247147529, 0.2955242247147529, 0.2692667193099963, 0.2190863625159820, 0.1494513491505806, 0.0666713443086881},
}

// numintfixedGL evaluates the fixed Gauss-Legendre rule of the given order on
// [a, b] using the tabulated nodes and weights.
func numintfixedGL(f Func, a, b float64, order int) float64 {
	nodes := numintglNodes[order]
	weights := numintglWeights[order]
	c := 0.5 * (b - a)
	d := 0.5 * (b + a)
	sum := 0.0
	for i := range nodes {
		sum += weights[i] * f(c*nodes[i]+d)
	}
	return sum * c
}

// GaussLegendre2 approximates the integral of f over [a, b] using the
// two-point Gauss-Legendre rule (exact for polynomials up to degree 3).
func GaussLegendre2(f Func, a, b float64) float64 { return numintfixedGL(f, a, b, 2) }

// GaussLegendre3 approximates the integral of f over [a, b] using the
// three-point Gauss-Legendre rule (exact for polynomials up to degree 5).
func GaussLegendre3(f Func, a, b float64) float64 { return numintfixedGL(f, a, b, 3) }

// GaussLegendre4 approximates the integral of f over [a, b] using the
// four-point Gauss-Legendre rule (exact for polynomials up to degree 7).
func GaussLegendre4(f Func, a, b float64) float64 { return numintfixedGL(f, a, b, 4) }

// GaussLegendre5 approximates the integral of f over [a, b] using the
// five-point Gauss-Legendre rule (exact for polynomials up to degree 9).
func GaussLegendre5(f Func, a, b float64) float64 { return numintfixedGL(f, a, b, 5) }

// GaussLegendre6 approximates the integral of f over [a, b] using the
// six-point Gauss-Legendre rule (exact for polynomials up to degree 11).
func GaussLegendre6(f Func, a, b float64) float64 { return numintfixedGL(f, a, b, 6) }

// GaussLegendre7 approximates the integral of f over [a, b] using the
// seven-point Gauss-Legendre rule (exact for polynomials up to degree 13).
func GaussLegendre7(f Func, a, b float64) float64 { return numintfixedGL(f, a, b, 7) }

// GaussLegendre8 approximates the integral of f over [a, b] using the
// eight-point Gauss-Legendre rule (exact for polynomials up to degree 15).
func GaussLegendre8(f Func, a, b float64) float64 { return numintfixedGL(f, a, b, 8) }

// GaussLegendre9 approximates the integral of f over [a, b] using the
// nine-point Gauss-Legendre rule (exact for polynomials up to degree 17).
func GaussLegendre9(f Func, a, b float64) float64 { return numintfixedGL(f, a, b, 9) }

// GaussLegendre10 approximates the integral of f over [a, b] using the
// ten-point Gauss-Legendre rule (exact for polynomials up to degree 19).
func GaussLegendre10(f Func, a, b float64) float64 { return numintfixedGL(f, a, b, 10) }

// GaussLegendreNodes returns the abscissae and weights of the n-point
// Gauss-Legendre rule on the reference interval [-1, 1]. Tabulated values are
// used for 2 <= n <= 10; higher orders are computed by Newton iteration on the
// Legendre polynomials.
func GaussLegendreNodes(n int) (nodes, weights []float64) {
	if n < 1 {
		n = 1
	}
	if nd, ok := numintglNodes[n]; ok {
		nodes = append([]float64(nil), nd...)
		weights = append([]float64(nil), numintglWeights[n]...)
		return nodes, weights
	}
	return numintlegendre(n)
}

// numintlegendre computes the n-point Gauss-Legendre nodes and weights on
// [-1, 1] via Newton iteration.
func numintlegendre(n int) (nodes, weights []float64) {
	nodes = make([]float64, n)
	weights = make([]float64, n)
	if n == 1 {
		nodes[0] = 0
		weights[0] = 2
		return nodes, weights
	}
	m := (n + 1) / 2
	for i := 0; i < m; i++ {
		z := math.Cos(math.Pi * (float64(i) + 0.75) / (float64(n) + 0.5))
		var pp float64
		for iter := 0; iter < 100; iter++ {
			p1 := 1.0
			p2 := 0.0
			for j := 0; j < n; j++ {
				p3 := p2
				p2 = p1
				p1 = ((2*float64(j)+1)*z*p2 - float64(j)*p3) / (float64(j) + 1)
			}
			pp = float64(n) * (z*p1 - p2) / (z*z - 1)
			z1 := z
			z = z1 - p1/pp
			if numintabs(z-z1) < 1e-15 {
				break
			}
		}
		nodes[i] = -z
		nodes[n-1-i] = z
		w := 2 / ((1 - z*z) * pp * pp)
		weights[i] = w
		weights[n-1-i] = w
	}
	return nodes, weights
}

// GaussLegendreN approximates the integral of f over [a, b] using the n-point
// Gauss-Legendre rule for any n >= 1.
func GaussLegendreN(f Func, a, b float64, n int) float64 {
	nodes, weights := GaussLegendreNodes(n)
	c := 0.5 * (b - a)
	d := 0.5 * (b + a)
	sum := 0.0
	for i := range nodes {
		sum += weights[i] * f(c*nodes[i]+d)
	}
	return sum * c
}

// CompositeGaussLegendre approximates the integral of f over [a, b] by
// splitting the interval into panels equal subintervals and applying the
// n-point Gauss-Legendre rule to each.
func CompositeGaussLegendre(f Func, a, b float64, panels, n int) float64 {
	if panels < 1 {
		panels = 1
	}
	nodes, weights := GaussLegendreNodes(n)
	h := (b - a) / float64(panels)
	total := 0.0
	for p := 0; p < panels; p++ {
		lo := a + float64(p)*h
		c := 0.5 * h
		d := lo + c
		sum := 0.0
		for i := range nodes {
			sum += weights[i] * f(c*nodes[i]+d)
		}
		total += sum * c
	}
	return total
}

// --- Gauss-Lobatto rules --------------------------------------------------

// GaussLobatto3 approximates the integral of f over [a, b] using the
// three-point Gauss-Lobatto rule (Simpson's rule), which includes both
// endpoints.
func GaussLobatto3(f Func, a, b float64) float64 {
	m := 0.5 * (a + b)
	return (b - a) / 6 * (f(a) + 4*f(m) + f(b))
}

// GaussLobatto4 approximates the integral of f over [a, b] using the
// four-point Gauss-Lobatto rule, which includes both endpoints and the
// interior nodes at +/- 1/sqrt(5).
func GaussLobatto4(f Func, a, b float64) float64 {
	const node = 0.4472135954999579 // 1/sqrt(5)
	c := 0.5 * (b - a)
	d := 0.5 * (b + a)
	return c * (1.0/6.0*(f(a)+f(b)) +
		5.0/6.0*(f(-node*c+d)+f(node*c+d)))
}

// GaussLobatto5 approximates the integral of f over [a, b] using the
// five-point Gauss-Lobatto rule, which includes both endpoints, the centre,
// and the interior nodes at +/- sqrt(3/7).
func GaussLobatto5(f Func, a, b float64) float64 {
	const node = 0.6546536707079771 // sqrt(3/7)
	c := 0.5 * (b - a)
	d := 0.5 * (b + a)
	return c * (0.1*(f(a)+f(b)) +
		49.0/90.0*(f(-node*c+d)+f(node*c+d)) +
		32.0/45.0*f(d))
}

// --- Gauss-Kronrod 15 -----------------------------------------------------

// numintxgk holds the abscissae of the 15-point Gauss-Kronrod rule on [-1, 1],
// listed from the outermost positive node to zero.
var numintxgk = []float64{
	0.991455371120813,
	0.949107912342759,
	0.864864423359769,
	0.741531185599394,
	0.586087235467691,
	0.405845151377397,
	0.207784955007898,
	0.000000000000000,
}

// numintwgk holds the 15-point Kronrod weights corresponding to numintxgk.
var numintwgk = []float64{
	0.022935322010529,
	0.063092092629979,
	0.104790010322250,
	0.140653259715525,
	0.169004726639267,
	0.190350578064785,
	0.204432940075298,
	0.209482141084728,
}

// numintwg holds the 7-point Gauss weights for the nodes at the odd indices of
// numintxgk.
var numintwg = []float64{
	0.129484966168870,
	0.279705391489277,
	0.381830050505119,
	0.417959183673469,
}

// GaussKronrod15Nodes returns the positive-half abscissae of the 15-point
// Gauss-Kronrod rule on [-1, 1], the corresponding Kronrod weights, and the
// four embedded 7-point Gauss weights. The full rule is symmetric about zero.
func GaussKronrod15Nodes() (abscissae, kronrodWeights, gaussWeights []float64) {
	abscissae = append([]float64(nil), numintxgk...)
	kronrodWeights = append([]float64(nil), numintwgk...)
	gaussWeights = append([]float64(nil), numintwg...)
	return abscissae, kronrodWeights, gaussWeights
}

// GaussKronrod15 approximates the integral of f over [a, b] using the 15-point
// Gauss-Kronrod rule. It returns the Kronrod estimate together with a
// conservative estimate of the absolute error derived from the difference
// between the Kronrod and embedded 7-point Gauss estimates.
func GaussKronrod15(f Func, a, b float64) (value, errEst float64) {
	centre := 0.5 * (a + b)
	hlgth := 0.5 * (b - a)
	dhlgth := numintabs(hlgth)

	fc := f(centre)
	resg := numintwg[3] * fc
	resk := numintwgk[7] * fc
	resabs := numintabs(resk)

	fv1 := make([]float64, 8)
	fv2 := make([]float64, 8)

	for j := 0; j < 3; j++ {
		jtw := 2*j + 1
		absc := hlgth * numintxgk[jtw]
		f1 := f(centre - absc)
		f2 := f(centre + absc)
		fv1[jtw] = f1
		fv2[jtw] = f2
		fsum := f1 + f2
		resg += numintwg[j] * fsum
		resk += numintwgk[jtw] * fsum
		resabs += numintwgk[jtw] * (numintabs(f1) + numintabs(f2))
	}
	for j := 0; j < 4; j++ {
		jtwm1 := 2 * j
		absc := hlgth * numintxgk[jtwm1]
		f1 := f(centre - absc)
		f2 := f(centre + absc)
		fv1[jtwm1] = f1
		fv2[jtwm1] = f2
		fsum := f1 + f2
		resk += numintwgk[jtwm1] * fsum
		resabs += numintwgk[jtwm1] * (numintabs(f1) + numintabs(f2))
	}

	reskh := resk * 0.5
	resasc := numintwgk[7] * numintabs(fc-reskh)
	for j := 0; j < 8; j++ {
		resasc += numintwgk[j] * (numintabs(fv1[j]-reskh) + numintabs(fv2[j]-reskh))
	}

	value = resk * hlgth
	resabs *= dhlgth
	resasc *= dhlgth
	errEst = numintabs((resk - resg) * hlgth)
	if resasc != 0 && errEst != 0 {
		scale := math.Pow(200*errEst/resasc, 1.5)
		if scale < 1 {
			errEst = resasc * scale
		} else {
			errEst = resasc
		}
	}
	const uflow = 2.2250738585072014e-308
	const eps = 2.220446049250313e-16
	if resabs > uflow/(50*eps) {
		if min := 50 * eps * resabs; min > errEst {
			errEst = min
		}
	}
	return value, errEst
}

// --- Clenshaw-Curtis and Gauss-Chebyshev ----------------------------------

// ClenshawCurtis approximates the integral of f over [a, b] using
// Clenshaw-Curtis quadrature with n+1 Chebyshev extreme points; n must be a
// positive even number and is rounded up when odd.
func ClenshawCurtis(f Func, a, b float64, n int) float64 {
	n = numintevenUp(n)
	c := 0.5 * (b - a)
	d := 0.5 * (b + a)
	sum := 0.0
	for k := 0; k <= n; k++ {
		theta := float64(k) * math.Pi / float64(n)
		x := math.Cos(theta)
		inner := 0.0
		for j := 1; j <= n/2; j++ {
			bj := 2.0
			if j == n/2 {
				bj = 1.0
			}
			inner += bj / float64(4*j*j-1) * math.Cos(2*float64(j)*theta)
		}
		ck := 2.0
		if k == 0 || k == n {
			ck = 1.0
		}
		w := ck / float64(n) * (1 - inner)
		sum += w * f(c*x+d)
	}
	return sum * c
}

// GaussChebyshev approximates the weighted integral of f over [-1, 1] with the
// Chebyshev weight 1/sqrt(1-x^2), that is the integral of f(x)/sqrt(1-x^2),
// using the n-point Gauss-Chebyshev rule of the first kind.
func GaussChebyshev(f Func, n int) float64 {
	if n < 1 {
		n = 1
	}
	w := math.Pi / float64(n)
	sum := 0.0
	for i := 1; i <= n; i++ {
		x := math.Cos(math.Pi * (float64(i) - 0.5) / float64(n))
		sum += f(x)
	}
	return w * sum
}

// --- tanh-sinh (double-exponential) ---------------------------------------

// TanhSinh approximates the integral of f over [a, b] using the tanh-sinh
// (double-exponential) rule with 2n+1 abscissae per level. The rule places
// nodes that cluster near the endpoints, making it well suited to integrands
// with integrable endpoint singularities. Non-finite integrand samples in the
// negligible tail are skipped.
func TanhSinh(f Func, a, b float64, n int) float64 {
	if n < 1 {
		n = 1
	}
	const tmax = 3.0
	h := tmax / float64(n)
	c := 0.5 * (b - a)
	d := 0.5 * (b + a)
	sum := 0.0
	for k := -n; k <= n; k++ {
		t := float64(k) * h
		u := 0.5 * math.Pi * math.Sinh(t)
		ch := math.Cosh(u)
		x := math.Tanh(u)
		w := 0.5 * math.Pi * math.Cosh(t) / (ch * ch)
		xa := c*x + d
		fx := f(xa)
		if math.IsNaN(fx) || math.IsInf(fx, 0) {
			continue
		}
		sum += w * fx
	}
	return sum * c * h
}

// TanhSinhTol approximates the integral of f over [a, b] using the tanh-sinh
// rule, refining the abscissa spacing until two successive level estimates
// differ by less than tol in absolute value or a level cap is reached.
func TanhSinhTol(f Func, a, b, tol float64) float64 {
	prev := TanhSinh(f, a, b, 8)
	n := 16
	for level := 0; level < 8; level++ {
		cur := TanhSinh(f, a, b, n)
		if numintabs(cur-prev) < tol {
			return cur
		}
		prev = cur
		n *= 2
	}
	return prev
}

// DoubleExponential is a synonym for TanhSinh, using the alternative name for
// the double-exponential quadrature rule.
func DoubleExponential(f Func, a, b float64, n int) float64 {
	return TanhSinh(f, a, b, n)
}

// --- high-level convenience -----------------------------------------------

// Integrate approximates the integral of f over [a, b] to roughly full
// double-precision accuracy using globally adaptive Gauss-Kronrod quadrature.
func Integrate(f Func, a, b float64) float64 {
	return AdaptiveGaussKronrod(f, a, b, 1e-12).Value
}

// RichardsonExtrapolate combines a coarse and a fine estimate of the same
// quantity, computed with step sizes differing by a factor of two, into a
// higher-order estimate, given the leading error order p. It returns
// (2^p * fine - coarse) / (2^p - 1).
func RichardsonExtrapolate(coarse, fine, order float64) float64 {
	pow := math.Pow(2, order)
	return (pow*fine - coarse) / (pow - 1)
}

// --- multidimensional tensor-product rules --------------------------------

// numinttensor2 combines one-dimensional abscissae and weights into a
// tensor-product estimate of the double integral of f.
func numinttensor2(f Func2, xs1, ws1, xs2, ws2 []float64) float64 {
	total := 0.0
	for i := range xs1 {
		row := 0.0
		for j := range xs2 {
			row += ws2[j] * f(xs1[i], xs2[j])
		}
		total += ws1[i] * row
	}
	return total
}

// Trapezoid2D approximates the double integral of f over the rectangle
// [ax, bx] x [ay, by] using the composite trapezoid rule with nx panels in x
// and ny panels in y.
func Trapezoid2D(f Func2, ax, bx, ay, by float64, nx, ny int) float64 {
	xs1, ws1 := numinttrapWeights(ax, bx, nx)
	xs2, ws2 := numinttrapWeights(ay, by, ny)
	return numinttensor2(f, xs1, ws1, xs2, ws2)
}

// Simpson2D approximates the double integral of f over the rectangle
// [ax, bx] x [ay, by] using the composite Simpson rule with nx panels in x and
// ny panels in y (each rounded up to an even number).
func Simpson2D(f Func2, ax, bx, ay, by float64, nx, ny int) float64 {
	xs1, ws1 := numintsimpWeights(ax, bx, nx)
	xs2, ws2 := numintsimpWeights(ay, by, ny)
	return numinttensor2(f, xs1, ws1, xs2, ws2)
}

// Midpoint2D approximates the double integral of f over the rectangle
// [ax, bx] x [ay, by] using the composite midpoint rule with nx panels in x
// and ny panels in y.
func Midpoint2D(f Func2, ax, bx, ay, by float64, nx, ny int) float64 {
	xs1, ws1 := numintmidWeights(ax, bx, nx)
	xs2, ws2 := numintmidWeights(ay, by, ny)
	return numinttensor2(f, xs1, ws1, xs2, ws2)
}

// GaussLegendre2D approximates the double integral of f over the rectangle
// [ax, bx] x [ay, by] using the tensor product of the n-point Gauss-Legendre
// rule in each direction.
func GaussLegendre2D(f Func2, ax, bx, ay, by float64, n int) float64 {
	nodes, weights := GaussLegendreNodes(n)
	cx := 0.5 * (bx - ax)
	dx := 0.5 * (bx + ax)
	cy := 0.5 * (by - ay)
	dy := 0.5 * (by + ay)
	total := 0.0
	for i := range nodes {
		xi := cx*nodes[i] + dx
		row := 0.0
		for j := range nodes {
			yj := cy*nodes[j] + dy
			row += weights[j] * f(xi, yj)
		}
		total += weights[i] * row
	}
	return total * cx * cy
}

// Trapezoid3D approximates the triple integral of f over the box
// [ax, bx] x [ay, by] x [az, bz] using the composite trapezoid rule with the
// given number of panels in each direction.
func Trapezoid3D(f Func3, ax, bx, ay, by, az, bz float64, nx, ny, nz int) float64 {
	xs1, ws1 := numinttrapWeights(ax, bx, nx)
	xs2, ws2 := numinttrapWeights(ay, by, ny)
	xs3, ws3 := numinttrapWeights(az, bz, nz)
	return numinttensor3(f, xs1, ws1, xs2, ws2, xs3, ws3)
}

// Simpson3D approximates the triple integral of f over the box
// [ax, bx] x [ay, by] x [az, bz] using the composite Simpson rule with the
// given number of panels in each direction (each rounded up to an even number).
func Simpson3D(f Func3, ax, bx, ay, by, az, bz float64, nx, ny, nz int) float64 {
	xs1, ws1 := numintsimpWeights(ax, bx, nx)
	xs2, ws2 := numintsimpWeights(ay, by, ny)
	xs3, ws3 := numintsimpWeights(az, bz, nz)
	return numinttensor3(f, xs1, ws1, xs2, ws2, xs3, ws3)
}

// numinttensor3 combines one-dimensional abscissae and weights into a
// tensor-product estimate of the triple integral of f.
func numinttensor3(f Func3, xs1, ws1, xs2, ws2, xs3, ws3 []float64) float64 {
	total := 0.0
	for i := range xs1 {
		for j := range xs2 {
			wij := ws1[i] * ws2[j]
			for k := range xs3 {
				total += wij * ws3[k] * f(xs1[i], xs2[j], xs3[k])
			}
		}
	}
	return total
}

// GaussLegendre3D approximates the triple integral of f over the box
// [ax, bx] x [ay, by] x [az, bz] using the tensor product of the n-point
// Gauss-Legendre rule in each direction.
func GaussLegendre3D(f Func3, ax, bx, ay, by, az, bz float64, n int) float64 {
	nodes, weights := GaussLegendreNodes(n)
	cx := 0.5 * (bx - ax)
	dx := 0.5 * (bx + ax)
	cy := 0.5 * (by - ay)
	dy := 0.5 * (by + ay)
	cz := 0.5 * (bz - az)
	dz := 0.5 * (bz + az)
	total := 0.0
	for i := range nodes {
		xi := cx*nodes[i] + dx
		for j := range nodes {
			yj := cy*nodes[j] + dy
			wij := weights[i] * weights[j]
			for k := range nodes {
				zk := cz*nodes[k] + dz
				total += wij * weights[k] * f(xi, yj, zk)
			}
		}
	}
	return total * cx * cy * cz
}

// Integrate2D approximates the double integral of f over the rectangle
// [ax, bx] x [ay, by] to high accuracy using an 8-point Gauss-Legendre tensor
// product refined over a modest grid of panels.
func Integrate2D(f Func2, ax, bx, ay, by float64) float64 {
	const panels = 8
	hx := (bx - ax) / panels
	hy := (by - ay) / panels
	total := 0.0
	for px := 0; px < panels; px++ {
		lox := ax + float64(px)*hx
		for py := 0; py < panels; py++ {
			loy := ay + float64(py)*hy
			total += GaussLegendre2D(f, lox, lox+hx, loy, loy+hy, 8)
		}
	}
	return total
}

// IntegrateND approximates the integral of f over the axis-aligned box defined
// by the lower and upper bounds using an n-point Gauss-Legendre tensor product
// in each dimension. The two bound slices must have equal length equal to the
// dimensionality.
func IntegrateND(f FuncND, lower, upper []float64, n int) float64 {
	d := len(lower)
	if d == 0 || len(upper) != d {
		return 0
	}
	nodes, weights := GaussLegendreNodes(n)
	point := make([]float64, d)
	scale := 1.0
	for i := 0; i < d; i++ {
		scale *= 0.5 * (upper[i] - lower[i])
	}
	var rec func(dim int, w float64) float64
	rec = func(dim int, w float64) float64 {
		if dim == d {
			return w * f(point)
		}
		c := 0.5 * (upper[dim] - lower[dim])
		m := 0.5 * (upper[dim] + lower[dim])
		sum := 0.0
		for i := range nodes {
			point[dim] = c*nodes[i] + m
			sum += rec(dim+1, w*weights[i])
		}
		return sum
	}
	return rec(0, 1) * scale
}

// --- deterministic Monte-Carlo integration --------------------------------

// MonteCarlo approximates the integral of f over [a, b] using n uniform random
// samples drawn from a generator seeded with seed, making the result
// deterministic for a given seed.
func MonteCarlo(f Func, a, b float64, n int, seed int64) float64 {
	if n < 1 {
		n = 1
	}
	rng := rand.New(rand.NewSource(seed))
	sum := 0.0
	for i := 0; i < n; i++ {
		x := a + (b-a)*rng.Float64()
		sum += f(x)
	}
	return (b - a) * sum / float64(n)
}

// MonteCarlo2D approximates the double integral of f over the rectangle
// [ax, bx] x [ay, by] using n uniform random samples from a generator seeded
// with seed.
func MonteCarlo2D(f Func2, ax, bx, ay, by float64, n int, seed int64) float64 {
	if n < 1 {
		n = 1
	}
	rng := rand.New(rand.NewSource(seed))
	sum := 0.0
	for i := 0; i < n; i++ {
		x := ax + (bx-ax)*rng.Float64()
		y := ay + (by-ay)*rng.Float64()
		sum += f(x, y)
	}
	area := (bx - ax) * (by - ay)
	return area * sum / float64(n)
}

// MonteCarloND approximates the integral of f over the axis-aligned box
// defined by lower and upper using n uniform random samples from a generator
// seeded with seed. The two bound slices must have equal length.
func MonteCarloND(f FuncND, lower, upper []float64, n int, seed int64) float64 {
	d := len(lower)
	if d == 0 || len(upper) != d || n < 1 {
		return 0
	}
	rng := rand.New(rand.NewSource(seed))
	vol := 1.0
	for i := 0; i < d; i++ {
		vol *= upper[i] - lower[i]
	}
	point := make([]float64, d)
	sum := 0.0
	for i := 0; i < n; i++ {
		for k := 0; k < d; k++ {
			point[k] = lower[k] + (upper[k]-lower[k])*rng.Float64()
		}
		sum += f(point)
	}
	return vol * sum / float64(n)
}
