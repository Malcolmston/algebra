package approxtheory

import "math"

// EvalDerivative evaluates the derivative of the Chebyshev series at x.
func (s *ChebSeries) EvalDerivative(x float64) float64 {
	return s.Derivative().Eval(x)
}

// Coefficient returns the k-th Chebyshev coefficient, or 0 when k is out of
// range.
func (s *ChebSeries) Coefficient(k int) float64 {
	if k < 0 || k >= len(s.Coeffs) {
		return 0
	}
	return s.Coeffs[k]
}

// Remap returns a new series with the same coefficients reinterpreted on the
// interval [a, b]; the underlying polynomial shape in the unit variable is
// preserved but its domain changes.
func (s *ChebSeries) Remap(a, b float64) *ChebSeries {
	return NewChebSeries(s.Coeffs, a, b)
}

// ToMonomial converts the Chebyshev series into monomial coefficients
// (ascending) of the polynomial in the original variable x on [A, B].
func (s *ChebSeries) ToMonomial() []float64 {
	inT := ChebyshevToMonomial(s.Coeffs) // polynomial in t
	// t = alpha*x + beta
	alpha := 2.0 / (s.B - s.A)
	beta := -(s.A + s.B) / (s.B - s.A)
	return polyComposeLinear(inT, alpha, beta)
}

// L2Norm returns the L2 norm of the series over its domain, sqrt(int f^2 dx),
// computed exactly from the coefficients via the Chebyshev product.
func (s *ChebSeries) L2Norm() float64 {
	sq := ChebProduct(s, s).Integral()
	if sq < 0 {
		sq = 0
	}
	return math.Sqrt(sq)
}

// MaxAbs returns the maximum absolute value of the series over its domain,
// estimated on a fine grid.
func (s *ChebSeries) MaxAbs() float64 {
	grid := Linspace(s.A, s.B, 2000)
	var m float64
	for _, x := range grid {
		if v := math.Abs(s.Eval(x)); v > m {
			m = v
		}
	}
	return m
}

// polyComposeLinear returns the coefficients of p(alpha*x + beta) given the
// coefficients of p (ascending).
func polyComposeLinear(p []float64, alpha, beta float64) []float64 {
	// Horner: build up using the linear polynomial (beta + alpha*x).
	lin := []float64{beta, alpha}
	var out []float64
	for i := len(p) - 1; i >= 0; i-- {
		out = PolyMul(out, lin)
		out = PolyAdd(out, []float64{p[i]})
	}
	if out == nil {
		out = []float64{0}
	}
	return out
}

// PolyCompose returns the coefficients of the composition p(q(x)) for two
// monomial polynomials.
func PolyCompose(p, q []float64) []float64 {
	var out []float64
	for i := len(p) - 1; i >= 0; i-- {
		out = PolyMul(out, q)
		out = PolyAdd(out, []float64{p[i]})
	}
	if out == nil {
		out = []float64{0}
	}
	return out
}

// PolyShift returns the coefficients of p(x + h).
func PolyShift(p []float64, h float64) []float64 {
	return polyComposeLinear(p, 1, h)
}

// GaussChebyshev2Nodes returns the n Gauss-Chebyshev nodes of the second kind,
// x_k = cos(k*pi/(n+1)), mapped to [a, b] in ascending order.
func GaussChebyshev2Nodes(n int, a, b float64) []float64 {
	out := make([]float64, n)
	for k := 1; k <= n; k++ {
		t := math.Cos(float64(k) * math.Pi / float64(n+1))
		out[n-k] = 0.5*(a+b) + 0.5*(b-a)*t
	}
	return out
}

// GaussChebyshev2Weights returns the n weights of Gauss-Chebyshev quadrature of
// the second kind on [a, b], w_k = pi/(n+1) sin^2(k*pi/(n+1)) scaled by
// (b-a)/2. These integrate f(x)*sqrt(1-t^2).
func GaussChebyshev2Weights(n int, a, b float64) []float64 {
	out := make([]float64, n)
	scale := (b - a) / 2
	for k := 1; k <= n; k++ {
		s := math.Sin(float64(k) * math.Pi / float64(n+1))
		out[n-k] = math.Pi / float64(n+1) * s * s * scale
	}
	return out
}

// GaussChebyshev2Quadrature approximates the plain integral of f over [a, b]
// using n-point Gauss-Chebyshev quadrature of the second kind, folding the
// weight sqrt(1-t^2) back in.
func GaussChebyshev2Quadrature(f func(float64) float64, n int, a, b float64) float64 {
	nodes := GaussChebyshev2Nodes(n, a, b)
	var sum float64
	for k := 1; k <= n; k++ {
		theta := float64(k) * math.Pi / float64(n+1)
		s := math.Sin(theta)
		// weight/sqrt(1-t^2) = pi/(n+1) * sin(theta)
		sum += math.Pi / float64(n+1) * s * f(nodes[n-k])
	}
	return sum * (b - a) / 2
}

// TaylorSinh returns the first n+1 Taylor coefficients of sinh about 0.
func TaylorSinh(n int) []float64 {
	out := make([]float64, n+1)
	fact := 1.0
	for k := 0; k <= n; k++ {
		if k > 0 {
			fact *= float64(k)
		}
		if k%2 == 1 {
			out[k] = 1.0 / fact
		}
	}
	return out
}

// TaylorCosh returns the first n+1 Taylor coefficients of cosh about 0.
func TaylorCosh(n int) []float64 {
	out := make([]float64, n+1)
	fact := 1.0
	for k := 0; k <= n; k++ {
		if k > 0 {
			fact *= float64(k)
		}
		if k%2 == 0 {
			out[k] = 1.0 / fact
		}
	}
	return out
}

// TaylorExpm1 returns the first n+1 Taylor coefficients of exp(x)-1 about 0.
func TaylorExpm1(n int) []float64 {
	c := TaylorExp(n)
	c[0] = 0
	return c
}

// TaylorAtan returns the first n+1 Taylor coefficients of atan about 0.
func TaylorAtan(n int) []float64 {
	out := make([]float64, n+1)
	for k := 1; k <= n; k += 2 {
		sign := 1.0
		if (k-1)/2%2 == 1 {
			sign = -1.0
		}
		out[k] = sign / float64(k)
	}
	return out
}

// TaylorGeometric returns the first n+1 Taylor coefficients of 1/(1-x), all
// equal to one.
func TaylorGeometric(n int) []float64 {
	out := make([]float64, n+1)
	for k := range out {
		out[k] = 1
	}
	return out
}

// PadeSin returns the [m/n] Pade approximant of sin about 0.
func PadeSin(m, n int) *PadeResult {
	res, _ := PadeApprox(TaylorSin(m+n), m, n)
	return res
}

// PadeCos returns the [m/n] Pade approximant of cos about 0.
func PadeCos(m, n int) *PadeResult {
	res, _ := PadeApprox(TaylorCos(m+n), m, n)
	return res
}

// PadeLog1p returns the [m/n] Pade approximant of log(1+x) about 0.
func PadeLog1p(m, n int) *PadeResult {
	res, _ := PadeApprox(TaylorLog1p(m+n), m, n)
	return res
}

// ContinuedFractionEval evaluates the continued fraction
//
//	a[0] + b[1]/(a[1] + b[2]/(a[2] + ... ))
//
// with the given partial numerators b (b[0] is ignored) and denominators a,
// using the backward recurrence. len(a) must equal len(b).
func ContinuedFractionEval(a, b []float64) float64 {
	n := len(a)
	if n == 0 {
		return math.NaN()
	}
	acc := 0.0
	for i := n - 1; i >= 1; i-- {
		denom := a[i] + acc
		if denom == 0 {
			denom = 1e-300
		}
		acc = b[i] / denom
	}
	return a[0] + acc
}

// BernsteinFromMonomial returns the Bernstein control values of degree n for a
// monomial polynomial defined on [0, 1]. It requires n >= degree of the input.
func BernsteinFromMonomial(mono []float64, n int) []float64 {
	out := make([]float64, n+1)
	for k := 0; k <= n; k++ {
		var s float64
		for i := 0; i <= k && i < len(mono); i++ {
			s += Binomial(k, i) / Binomial(n, i) * mono[i]
		}
		out[k] = s
	}
	return out
}

// BernsteinControlPoints returns the control values f(a + k/n (b-a)) used by
// the Bernstein approximation operator.
func BernsteinControlPoints(f func(float64) float64, n int, a, b float64) []float64 {
	out := make([]float64, n+1)
	for k := 0; k <= n; k++ {
		x := a
		if n > 0 {
			x = a + float64(k)/float64(n)*(b-a)
		}
		out[k] = f(x)
	}
	return out
}

// EconomizeChebyshev reduces the degree of a monomial polynomial on [a, b] by
// converting to the Chebyshev basis, discarding the highest-order coefficients
// whose cumulative absolute value does not exceed tol, and converting back. It
// returns the economized monomial polynomial on [a, b] and the discarded
// coefficient mass (an upper bound on the introduced uniform error).
func EconomizeChebyshev(mono []float64, a, b, tol float64) ([]float64, float64) {
	s := ChebFromMonomial(mono, a, b)
	c := s.Coeffs
	dropped := 0.0
	d := len(c)
	for d > 1 {
		next := dropped + math.Abs(c[d-1])
		if next > tol {
			break
		}
		dropped = next
		d--
	}
	reduced := &ChebSeries{Coeffs: c[:d], A: a, B: b}
	return reduced.ToMonomial(), dropped
}

// LebesgueConstantEquispacedAsymptotic returns the classical asymptotic growth
// estimate 2^{n+1} / (e * n * ln n) of the Lebesgue constant for n+1
// equispaced nodes (valid for large n).
func LebesgueConstantEquispacedAsymptotic(n int) float64 {
	if n < 2 {
		return 1
	}
	return math.Pow(2, float64(n+1)) / (math.E * float64(n) * math.Log(float64(n)))
}

// ConvergenceOrder estimates the observed order of convergence from two step
// sizes and their associated errors: log(e1/e2) / log(h1/h2).
func ConvergenceOrder(h1, e1, h2, e2 float64) float64 {
	return math.Log(e1/e2) / math.Log(h1/h2)
}

// NevilleEval interpolates (xs, ys) at x using Neville's algorithm and returns
// the interpolated value.
func NevilleEval(xs, ys []float64, x float64) float64 {
	n := len(xs)
	if n == 0 {
		return math.NaN()
	}
	p := make([]float64, n)
	copy(p, ys)
	for k := 1; k < n; k++ {
		for i := 0; i < n-k; i++ {
			p[i] = ((x-xs[i+k])*p[i] - (x-xs[i])*p[i+1]) / (xs[i] - xs[i+k])
		}
	}
	return p[0]
}

// RungeFunction evaluates Runge's function 1/(1+25 x^2), the classic example
// where high-degree equispaced interpolation diverges.
func RungeFunction(x float64) float64 {
	return 1 / (1 + 25*x*x)
}
