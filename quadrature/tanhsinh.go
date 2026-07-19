package quadrature

import "math"

// tanhSinhCore integrates g over the whole real line using the trapezoidal
// rule with successive step halving, where g already includes the change of
// variables and its Jacobian. It stops when successive refinements agree to
// tol or after a fixed number of levels, returning the estimate and an error
// estimate.
func tanhSinhCore(g func(t float64) float64, tol float64) (float64, float64) {
	const tmax = 6.5
	if tol <= 0 {
		tol = 1e-12
	}
	h := 1.0
	s := g(0)
	for t := h; t <= tmax; t += h {
		s += g(t) + g(-t)
	}
	old := s * h
	errEst := math.Abs(old)
	for level := 1; level <= 14; level++ {
		h *= 0.5
		var add float64
		for t := h; t <= tmax; t += 2 * h {
			add += g(t) + g(-t)
		}
		s += add
		cur := s * h
		errEst = math.Abs(cur - old)
		if errEst <= tol*math.Abs(cur)+tol {
			return cur, errEst
		}
		old = cur
	}
	return old, errEst
}

// safeVal returns v unless it is non-finite, in which case it returns 0. This
// lets the tanh-sinh weight suppress integrable endpoint singularities whose
// bare function value overflows.
func safeVal(v float64) float64 {
	if math.IsInf(v, 0) || math.IsNaN(v) {
		return 0
	}
	return v
}

// TanhSinh approximates the integral of f over the finite interval [a, b]
// using the tanh-sinh (double-exponential) rule. The rule is robust against
// integrable singularities at the endpoints because its weights decay
// doubly-exponentially there. It returns the estimate to the requested
// tolerance.
func TanhSinh(f Func, a, b, tol float64) float64 {
	if a == b {
		return 0
	}
	c := 0.5 * (b - a)
	mid := 0.5 * (a + b)
	const piHalf = math.Pi / 2
	g := func(t float64) float64 {
		s := math.Sinh(t)
		u := piHalf * s
		cu := math.Cosh(u)
		x := math.Tanh(u)
		w := piHalf * math.Cosh(t) / (cu * cu)
		return safeVal(f(mid+c*x)) * w
	}
	v, _ := tanhSinhCore(g, tol)
	return c * v
}

// TanhSinhResult is like TanhSinh but returns a Result with the error
// estimate.
func TanhSinhResult(f Func, a, b, tol float64) Result {
	if a == b {
		return Result{Success: true}
	}
	c := 0.5 * (b - a)
	mid := 0.5 * (a + b)
	const piHalf = math.Pi / 2
	g := func(t float64) float64 {
		s := math.Sinh(t)
		u := piHalf * s
		cu := math.Cosh(u)
		x := math.Tanh(u)
		w := piHalf * math.Cosh(t) / (cu * cu)
		return safeVal(f(mid+c*x)) * w
	}
	v, e := tanhSinhCore(g, tol)
	return Result{Value: c * v, AbsErr: c * e, Success: true}
}

// DoubleExponential is an alias for TanhSinh, the name by which the method is
// also known.
func DoubleExponential(f Func, a, b, tol float64) float64 {
	return TanhSinh(f, a, b, tol)
}

// SinhSinh approximates the integral of f over the entire real line
// (-inf, inf) using the sinh-sinh double-exponential transformation
// x = sinh((pi/2) sinh t). It is suited to functions that decay
// exponentially or faster at both infinities.
func SinhSinh(f Func, tol float64) float64 {
	const piHalf = math.Pi / 2
	g := func(t float64) float64 {
		s := math.Sinh(t)
		u := piHalf * s
		x := math.Sinh(u)
		w := piHalf * math.Cosh(t) * math.Cosh(u)
		return safeVal(f(x)) * w
	}
	v, _ := tanhSinhCore(g, tol)
	return v
}

// ExpSinh approximates the integral of f over the semi-infinite interval
// [a, inf) using the exp-sinh double-exponential transformation
// x = a + exp((pi/2) sinh t). It is suited to functions that decay
// exponentially or faster as x -> inf.
func ExpSinh(f Func, a, tol float64) float64 {
	const piHalf = math.Pi / 2
	g := func(t float64) float64 {
		u := piHalf * math.Sinh(t)
		ex := math.Exp(u)
		x := a + ex
		w := ex * piHalf * math.Cosh(t)
		return safeVal(f(x)) * w
	}
	v, _ := tanhSinhCore(g, tol)
	return v
}

// ExpSinhLeft approximates the integral of f over the semi-infinite interval
// (-inf, b] using the exp-sinh transformation x = b - exp((pi/2) sinh t).
func ExpSinhLeft(f Func, b, tol float64) float64 {
	const piHalf = math.Pi / 2
	g := func(t float64) float64 {
		u := piHalf * math.Sinh(t)
		ex := math.Exp(u)
		x := b - ex
		w := ex * piHalf * math.Cosh(t)
		return safeVal(f(x)) * w
	}
	v, _ := tanhSinhCore(g, tol)
	return v
}

// IntegrateInfinite approximates the integral of f over the entire real line
// using the substitution x = t/(1-t^2), t in (-1, 1), evaluated with the
// tanh-sinh rule. It is a convenient default when f is smooth and decays
// polynomially.
func IntegrateInfinite(f Func, tol float64) float64 {
	g := func(t float64) float64 {
		d := 1 - t*t
		x := t / d
		jac := (1 + t*t) / (d * d)
		return safeVal(f(x)) * jac
	}
	return TanhSinh(func(t float64) float64 { return g(t) }, -1, 1, tol)
}
