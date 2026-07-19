package quadrature

import (
	"math/big"
)

// Midpoint approximates the integral of f over [a, b] with the single-panel
// midpoint rule (b-a)*f((a+b)/2), exact for linear functions.
func Midpoint(f Func, a, b float64) float64 {
	return (b - a) * f(0.5*(a+b))
}

// Trapezoid approximates the integral of f over [a, b] with the single-panel
// trapezoidal rule (b-a)/2*(f(a)+f(b)), exact for linear functions.
func Trapezoid(f Func, a, b float64) float64 {
	return 0.5 * (b - a) * (f(a) + f(b))
}

// Simpson approximates the integral of f over [a, b] with Simpson's 1/3 rule
// using the two endpoints and the midpoint, exact for cubics.
func Simpson(f Func, a, b float64) float64 {
	m := 0.5 * (a + b)
	return (b - a) / 6 * (f(a) + 4*f(m) + f(b))
}

// Simpson38 approximates the integral of f over [a, b] with Simpson's 3/8
// rule using four equally spaced points, exact for cubics.
func Simpson38(f Func, a, b float64) float64 {
	h := (b - a) / 3
	return (b - a) / 8 * (f(a) + 3*f(a+h) + 3*f(a+2*h) + f(b))
}

// Boole approximates the integral of f over [a, b] with Boole's rule using
// five equally spaced points, exact for polynomials up to degree five.
func Boole(f Func, a, b float64) float64 {
	h := (b - a) / 4
	return (b - a) / 90 * (7*f(a) + 32*f(a+h) + 12*f(a+2*h) + 32*f(a+3*h) + 7*f(b))
}

// CompositeMidpoint approximates the integral of f over [a, b] using the
// midpoint rule on n equal subintervals.
func CompositeMidpoint(f Func, a, b float64, n int) float64 {
	if n < 1 {
		n = 1
	}
	h := (b - a) / float64(n)
	var s float64
	for i := 0; i < n; i++ {
		s += f(a + (float64(i)+0.5)*h)
	}
	return h * s
}

// CompositeTrapezoid approximates the integral of f over [a, b] using the
// trapezoidal rule on n equal subintervals.
func CompositeTrapezoid(f Func, a, b float64, n int) float64 {
	if n < 1 {
		n = 1
	}
	h := (b - a) / float64(n)
	s := 0.5 * (f(a) + f(b))
	for i := 1; i < n; i++ {
		s += f(a + float64(i)*h)
	}
	return h * s
}

// CompositeSimpson approximates the integral of f over [a, b] using the
// composite Simpson's 1/3 rule. n is the number of subintervals and is rounded
// up to the next even value.
func CompositeSimpson(f Func, a, b float64, n int) float64 {
	if n < 2 {
		n = 2
	}
	if n%2 == 1 {
		n++
	}
	h := (b - a) / float64(n)
	s := f(a) + f(b)
	for i := 1; i < n; i++ {
		x := a + float64(i)*h
		if i%2 == 1 {
			s += 4 * f(x)
		} else {
			s += 2 * f(x)
		}
	}
	return h / 3 * s
}

// CompositeSimpson38 approximates the integral of f over [a, b] using the
// composite Simpson's 3/8 rule. n is the number of subintervals and is rounded
// up to the next multiple of three.
func CompositeSimpson38(f Func, a, b float64, n int) float64 {
	if n < 3 {
		n = 3
	}
	if r := n % 3; r != 0 {
		n += 3 - r
	}
	h := (b - a) / float64(n)
	s := f(a) + f(b)
	for i := 1; i < n; i++ {
		x := a + float64(i)*h
		if i%3 == 0 {
			s += 2 * f(x)
		} else {
			s += 3 * f(x)
		}
	}
	return 3 * h / 8 * s
}

// CompositeBoole approximates the integral of f over [a, b] using the
// composite Boole rule. n is the number of subintervals and is rounded up to
// the next multiple of four.
func CompositeBoole(f Func, a, b float64, n int) float64 {
	if n < 4 {
		n = 4
	}
	if r := n % 4; r != 0 {
		n += 4 - r
	}
	h := (b - a) / float64(n)
	var s float64
	for i := 0; i < n; i += 4 {
		x := a + float64(i)*h
		s += 7*f(x) + 32*f(x+h) + 12*f(x+2*h) + 32*f(x+3*h) + 7*f(x+4*h)
	}
	return 2 * h / 45 * s
}

// TrapezoidSamples integrates uniformly spaced samples ys with spacing h using
// the composite trapezoidal rule. It requires at least two samples.
func TrapezoidSamples(ys []float64, h float64) float64 {
	n := len(ys)
	if n < 2 {
		return 0
	}
	s := 0.5 * (ys[0] + ys[n-1])
	for i := 1; i < n-1; i++ {
		s += ys[i]
	}
	return h * s
}

// SimpsonSamples integrates uniformly spaced samples ys with spacing h using
// the composite Simpson rule. It requires an odd number of samples (an even
// number of subintervals); when given an even count it falls back to the
// trapezoidal correction on the final panel.
func SimpsonSamples(ys []float64, h float64) float64 {
	n := len(ys)
	if n < 2 {
		return 0
	}
	if n == 2 {
		return 0.5 * h * (ys[0] + ys[1])
	}
	// number of panels
	panels := n - 1
	var s float64
	limit := panels
	odd := false
	if panels%2 == 1 {
		odd = true
		limit = panels - 1 // integrate first even number of panels with Simpson
	}
	s = ys[0] + ys[limit]
	for i := 1; i < limit; i++ {
		if i%2 == 1 {
			s += 4 * ys[i]
		} else {
			s += 2 * ys[i]
		}
	}
	result := h / 3 * s
	if odd {
		// add the last panel with a trapezoid
		result += 0.5 * h * (ys[limit] + ys[limit+1])
	}
	return result
}

// TrapezoidNonUniform integrates samples ys taken at the (not necessarily
// equally spaced) abscissae xs using the trapezoidal rule. The slices must
// have equal length of at least two.
func TrapezoidNonUniform(xs, ys []float64) float64 {
	n := len(xs)
	if n != len(ys) || n < 2 {
		return 0
	}
	var s float64
	for i := 1; i < n; i++ {
		s += 0.5 * (xs[i] - xs[i-1]) * (ys[i] + ys[i-1])
	}
	return s
}

// CumulativeTrapezoid returns the running integral of samples ys at abscissae
// xs computed with the trapezoidal rule. The returned slice has the same
// length as the inputs, with the first element equal to zero.
func CumulativeTrapezoid(xs, ys []float64) []float64 {
	n := len(xs)
	out := make([]float64, n)
	if n != len(ys) || n < 2 {
		return out
	}
	for i := 1; i < n; i++ {
		out[i] = out[i-1] + 0.5*(xs[i]-xs[i-1])*(ys[i]+ys[i-1])
	}
	return out
}

// NewtonCotesClosedExact returns the exact rational weights of the closed
// (n+1)-point Newton-Cotes rule on the reference interval [0, n], scaled so
// that the integral over [a, b] is (b-a)/n times the weighted sum of the
// samples. The returned weights sum to n. Panics for n < 1.
func NewtonCotesClosedExact(n int) []*big.Rat {
	if n < 1 {
		panic("quadrature: NewtonCotesClosedExact requires n >= 1")
	}
	// Nodes at 0,1,...,n. Weight_i = integral_0^n prod_{j!=i}(x-j)/(i-j) dx.
	weights := make([]*big.Rat, n+1)
	for i := 0; i <= n; i++ {
		weights[i] = lagrangeIntegral(n, i, closedNodes(n))
	}
	return weights
}

// NewtonCotesOpenExact returns the exact rational weights of the open n-point
// Newton-Cotes rule with interior nodes at 1,2,...,n on the reference interval
// [0, n+1], scaled so that the integral over [a, b] is (b-a)/(n+1) times the
// weighted sum. Panics for n < 1.
func NewtonCotesOpenExact(n int) []*big.Rat {
	if n < 1 {
		panic("quadrature: NewtonCotesOpenExact requires n >= 1")
	}
	nodes := make([]int64, n)
	for i := 0; i < n; i++ {
		nodes[i] = int64(i + 1)
	}
	weights := make([]*big.Rat, n)
	for i := 0; i < n; i++ {
		weights[i] = lagrangeIntegralOpen(n+1, i, nodes)
	}
	return weights
}

// NewtonCotesClosed returns the closed (n+1)-point Newton-Cotes weights as
// float64 values (see NewtonCotesClosedExact).
func NewtonCotesClosed(n int) []float64 {
	rats := NewtonCotesClosedExact(n)
	out := make([]float64, len(rats))
	for i, r := range rats {
		out[i], _ = r.Float64()
	}
	return out
}

// NewtonCotesOpen returns the open n-point Newton-Cotes weights as float64
// values (see NewtonCotesOpenExact).
func NewtonCotesOpen(n int) []float64 {
	rats := NewtonCotesOpenExact(n)
	out := make([]float64, len(rats))
	for i, r := range rats {
		out[i], _ = r.Float64()
	}
	return out
}

// IntegrateNewtonCotesClosed approximates the integral of f over [a, b] with a
// single closed (n+1)-point Newton-Cotes rule of order n.
func IntegrateNewtonCotesClosed(f Func, a, b float64, n int) float64 {
	w := NewtonCotesClosed(n)
	h := (b - a) / float64(n)
	var s float64
	for i := 0; i <= n; i++ {
		s += w[i] * f(a+float64(i)*h)
	}
	return h * s
}

// IntegrateNewtonCotesOpen approximates the integral of f over [a, b] with a
// single open n-point Newton-Cotes rule.
func IntegrateNewtonCotesOpen(f Func, a, b float64, n int) float64 {
	w := NewtonCotesOpen(n)
	h := (b - a) / float64(n+1)
	var s float64
	for i := 0; i < n; i++ {
		s += w[i] * f(a+float64(i+1)*h)
	}
	return h * s
}

// closedNodes returns the integer nodes 0..n.
func closedNodes(n int) []int64 {
	nodes := make([]int64, n+1)
	for i := 0; i <= n; i++ {
		nodes[i] = int64(i)
	}
	return nodes
}

// lagrangeIntegral returns the exact integral over [0, span] of the i-th
// Lagrange basis polynomial for the closed rule (span = n), divided by the
// panel width so that the weights sum to n. Concretely it computes
// integral_0^n L_i(x) dx where L_i has nodes 0..n.
func lagrangeIntegral(span, i int, nodes []int64) *big.Rat {
	// Build numerator polynomial prod_{j!=i}(x - nodes[j]) with rational
	// coefficients, integrate over [0, span], divide by prod_{j!=i}(i-j).
	coeffs := polyFromRoots(nodes, i)
	integral := integratePoly(coeffs, int64(span))
	denom := big.NewInt(1)
	for j := range nodes {
		if j == i {
			continue
		}
		denom.Mul(denom, big.NewInt(nodes[i]-nodes[j]))
	}
	res := new(big.Rat).SetFrac(integral.Num(), integral.Denom())
	res.Quo(res, new(big.Rat).SetInt(denom))
	return res
}

// lagrangeIntegralOpen behaves like lagrangeIntegral but integrates over
// [0, span] for the open rule whose nodes are given explicitly.
func lagrangeIntegralOpen(span, i int, nodes []int64) *big.Rat {
	coeffs := polyFromRoots(nodes, i)
	integral := integratePoly(coeffs, int64(span))
	denom := big.NewInt(1)
	for j := range nodes {
		if j == i {
			continue
		}
		denom.Mul(denom, big.NewInt(nodes[i]-nodes[j]))
	}
	res := new(big.Rat).Set(integral)
	res.Quo(res, new(big.Rat).SetInt(denom))
	return res
}

// polyFromRoots builds the polynomial prod_{j != skip}(x - nodes[j]) as a slice
// of big.Rat coefficients in ascending degree order.
func polyFromRoots(nodes []int64, skip int) []*big.Rat {
	poly := []*big.Rat{big.NewRat(1, 1)}
	for j, r := range nodes {
		if j == skip {
			continue
		}
		// multiply poly by (x - r)
		next := make([]*big.Rat, len(poly)+1)
		for k := range next {
			next[k] = new(big.Rat)
		}
		for k, c := range poly {
			// c * x term goes to degree k+1
			next[k+1].Add(next[k+1], c)
			// c * (-r) goes to degree k
			tmp := new(big.Rat).SetInt64(-r)
			tmp.Mul(tmp, c)
			next[k].Add(next[k], tmp)
		}
		poly = next
	}
	return poly
}

// integratePoly returns the exact integral over [0, span] of the polynomial
// with the given ascending-degree rational coefficients.
func integratePoly(coeffs []*big.Rat, span int64) *big.Rat {
	result := new(big.Rat)
	spanPow := big.NewInt(span)
	for k, c := range coeffs {
		// integral of c*x^k over [0,span] = c*span^{k+1}/(k+1)
		p := new(big.Int).Exp(spanPow, big.NewInt(int64(k+1)), nil)
		term := new(big.Rat).SetInt(p)
		term.Quo(term, new(big.Rat).SetInt64(int64(k+1)))
		term.Mul(term, c)
		result.Add(result, term)
	}
	return result
}
