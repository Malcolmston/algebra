package quadrature

import "math"

// IntegrateGaussLegendre approximates the integral of f over [a, b] with the
// n-point Gauss-Legendre rule, mapping the canonical nodes onto [a, b].
func IntegrateGaussLegendre(f Func, a, b float64, n int) float64 {
	nodes, weights := GaussLegendre(n)
	half := 0.5 * (b - a)
	mid := 0.5 * (a + b)
	var s float64
	for i, t := range nodes {
		s += weights[i] * f(mid+half*t)
	}
	return half * s
}

// IntegrateGaussLobatto approximates the integral of f over [a, b] with the
// n-point Gauss-Lobatto rule.
func IntegrateGaussLobatto(f Func, a, b float64, n int) float64 {
	nodes, weights := GaussLobatto(n)
	half := 0.5 * (b - a)
	mid := 0.5 * (a + b)
	var s float64
	for i, t := range nodes {
		s += weights[i] * f(mid+half*t)
	}
	return half * s
}

// IntegrateGaussRadau approximates the integral of f over [a, b] with the
// n-point left Gauss-Radau rule.
func IntegrateGaussRadau(f Func, a, b float64, n int) float64 {
	nodes, weights := GaussRadau(n)
	half := 0.5 * (b - a)
	mid := 0.5 * (a + b)
	var s float64
	for i, t := range nodes {
		s += weights[i] * f(mid+half*t)
	}
	return half * s
}

// IntegrateGaussHermite approximates the integral of f(x)*e^{-x^2} over the
// whole real line with the n-point Gauss-Hermite rule (physicists'
// convention).
func IntegrateGaussHermite(f Func, n int) float64 {
	nodes, weights := GaussHermite(n)
	var s float64
	for i, x := range nodes {
		s += weights[i] * f(x)
	}
	return s
}

// IntegrateGaussHermiteFunc approximates the integral of f (without any weight
// factor) over the whole real line by folding the Gaussian weight back into
// the summand: the summand is f(x_i)*e^{x_i^2}*w_i. It is accurate when f
// decays at least as fast as a Gaussian.
func IntegrateGaussHermiteFunc(f Func, n int) float64 {
	nodes, weights := GaussHermite(n)
	var s float64
	for i, x := range nodes {
		s += weights[i] * math.Exp(x*x) * f(x)
	}
	return s
}

// IntegrateGaussLaguerre approximates the integral of f(x)*e^{-x} over
// [0, inf) with the n-point Gauss-Laguerre rule.
func IntegrateGaussLaguerre(f Func, n int) float64 {
	nodes, weights := GaussLaguerre(n)
	var s float64
	for i, x := range nodes {
		s += weights[i] * f(x)
	}
	return s
}

// IntegrateGaussLaguerreFunc approximates the integral of f (without any
// weight factor) over [0, inf) by folding the e^{-x} weight back in: the
// summand is f(x_i)*e^{x_i}*w_i.
func IntegrateGaussLaguerreFunc(f Func, n int) float64 {
	nodes, weights := GaussLaguerre(n)
	var s float64
	for i, x := range nodes {
		s += weights[i] * math.Exp(x) * f(x)
	}
	return s
}

// IntegrateGaussLaguerreGen approximates the integral of f(x)*x^a*e^{-x} over
// [0, inf) with the n-point generalized Gauss-Laguerre rule.
func IntegrateGaussLaguerreGen(f Func, n int, a float64) float64 {
	nodes, weights := GaussLaguerreGen(n, a)
	var s float64
	for i, x := range nodes {
		s += weights[i] * f(x)
	}
	return s
}

// IntegrateGaussChebyshev1 approximates the integral of f(x)/sqrt(1-x^2) over
// [-1, 1] with the n-point first-kind Gauss-Chebyshev rule.
func IntegrateGaussChebyshev1(f Func, n int) float64 {
	nodes, weights := GaussChebyshev1(n)
	var s float64
	for i, x := range nodes {
		s += weights[i] * f(x)
	}
	return s
}

// IntegrateGaussChebyshev2 approximates the integral of f(x)*sqrt(1-x^2) over
// [-1, 1] with the n-point second-kind Gauss-Chebyshev rule.
func IntegrateGaussChebyshev2(f Func, n int) float64 {
	nodes, weights := GaussChebyshev2(n)
	var s float64
	for i, x := range nodes {
		s += weights[i] * f(x)
	}
	return s
}

// IntegrateGaussJacobi approximates the integral of f(x)*(1-x)^a*(1+x)^b over
// [-1, 1] with the n-point Gauss-Jacobi rule.
func IntegrateGaussJacobi(f Func, n int, a, b float64) float64 {
	nodes, weights := GaussJacobi(n, a, b)
	var s float64
	for i, x := range nodes {
		s += weights[i] * f(x)
	}
	return s
}

// ExpectationGaussHermite approximates E[f(X)] for a standard normal random
// variable X ~ N(0, 1) using the n-point probabilists' Gauss-Hermite rule.
// The result is (1/sqrt(2*pi)) times the weighted sum.
func ExpectationGaussHermite(f Func, n int) float64 {
	nodes, weights := GaussHermiteProb(n)
	var s float64
	for i, x := range nodes {
		s += weights[i] * f(x)
	}
	return s / math.Sqrt(2*math.Pi)
}

// ExpectationGaussHermiteNormal approximates E[f(X)] for a normal random
// variable X ~ N(mean, sd^2) using the n-point probabilists' Gauss-Hermite
// rule.
func ExpectationGaussHermiteNormal(f Func, mean, sd float64, n int) float64 {
	nodes, weights := GaussHermiteProb(n)
	var s float64
	for i, x := range nodes {
		s += weights[i] * f(mean+sd*x)
	}
	return s / math.Sqrt(2*math.Pi)
}
