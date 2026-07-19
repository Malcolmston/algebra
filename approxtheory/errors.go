package approxtheory

import "math"

// ChebyshevNodes returns the n+1 Chebyshev-Gauss-Lobatto nodes on [a, b],
// ascending. It is an alias for ChebPoints kept for readability in error
// analysis code.
func ChebyshevNodes(n int, a, b float64) []float64 {
	return ChebPoints(n, a, b)
}

// EquispacedNodes returns n+1 equispaced nodes spanning [a, b].
func EquispacedNodes(n int, a, b float64) []float64 {
	return Linspace(a, b, n+1)
}

// NodePolynomial evaluates the node polynomial w(x) = prod_i (x - nodes[i]) at
// x. Its size governs the interpolation error.
func NodePolynomial(nodes []float64, x float64) float64 {
	p := 1.0
	for _, xi := range nodes {
		p *= x - xi
	}
	return p
}

// MaxNodePolynomial returns the maximum of |w(x)| = |prod_i (x-nodes[i])| over
// [a, b], estimated on a fine grid.
func MaxNodePolynomial(nodes []float64, a, b float64) float64 {
	grid := Linspace(a, b, 2000)
	var m float64
	for _, x := range grid {
		if v := math.Abs(NodePolynomial(nodes, x)); v > m {
			m = v
		}
	}
	return m
}

// LebesgueFunction evaluates the Lebesgue function lambda(x) = sum_j
// |L_j(x)| for the interpolation nodes at x, where L_j are the Lagrange
// cardinal polynomials.
func LebesgueFunction(nodes []float64, x float64) float64 {
	var sum float64
	for j := range nodes {
		sum += math.Abs(LagrangeBasis(nodes, j, x))
	}
	return sum
}

// LebesgueConstant returns the Lebesgue constant, the maximum of the Lebesgue
// function over [a, b], estimated on a fine grid. It bounds how much larger the
// interpolation error can be than the best-approximation error.
func LebesgueConstant(nodes []float64, a, b float64) float64 {
	grid := Linspace(a, b, 4000)
	var m float64
	for _, x := range grid {
		if v := LebesgueFunction(nodes, x); v > m {
			m = v
		}
	}
	return m
}

// LebesgueConstantChebyshevAsymptotic returns the classical asymptotic estimate
// (2/pi) log(n+1) + 0.9625 for the Lebesgue constant of n+1 Chebyshev points.
func LebesgueConstantChebyshevAsymptotic(n int) float64 {
	return 2/math.Pi*math.Log(float64(n)+1) + 0.9625
}

// MaxError returns the maximum of |f(x) - g(x)| over [a, b], estimated on a
// fine grid. It is the standard uniform (sup-norm) error estimate.
func MaxError(f, g func(float64) float64, a, b float64) float64 {
	grid := Linspace(a, b, 4000)
	var m float64
	for _, x := range grid {
		if v := math.Abs(f(x) - g(x)); v > m {
			m = v
		}
	}
	return m
}

// MaxAbsError returns the maximum absolute difference between paired samples.
func MaxAbsError(a, b []float64) float64 {
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	var m float64
	for i := 0; i < n; i++ {
		if v := math.Abs(a[i] - b[i]); v > m {
			m = v
		}
	}
	return m
}

// L2Error returns the discrete root-integral L2 error of f-g over [a, b],
// estimated by the composite trapezoidal rule on a fine grid.
func L2Error(f, g func(float64) float64, a, b float64) float64 {
	const N = 4000
	grid := Linspace(a, b, N+1)
	h := (b - a) / float64(N)
	var sum float64
	for i, x := range grid {
		d := f(x) - g(x)
		w := 1.0
		if i == 0 || i == N {
			w = 0.5
		}
		sum += w * d * d
	}
	return math.Sqrt(sum * h)
}

// RMSError returns the root-mean-square difference between paired samples.
func RMSError(a, b []float64) float64 {
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	if n == 0 {
		return 0
	}
	var sum float64
	for i := 0; i < n; i++ {
		d := a[i] - b[i]
		sum += d * d
	}
	return math.Sqrt(sum / float64(n))
}

// MeanAbsError returns the mean absolute difference between paired samples.
func MeanAbsError(a, b []float64) float64 {
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	if n == 0 {
		return 0
	}
	var sum float64
	for i := 0; i < n; i++ {
		sum += math.Abs(a[i] - b[i])
	}
	return sum / float64(n)
}

// RelativeError returns |approx-exact| / max(|exact|, eps) guarding against
// division by zero with the small floor eps.
func RelativeError(approx, exact, eps float64) float64 {
	d := math.Abs(exact)
	if d < eps {
		d = eps
	}
	return math.Abs(approx-exact) / d
}

// InterpErrorBound returns the standard bound on the polynomial interpolation
// error, |f(x)-p(x)| <= (M / (n+1)!) * max|w(x)|, where n+1 is the number of
// nodes, derivBound is a bound on |f^{(n+1)}| over the interval and nodes are
// the interpolation abscissae.
func InterpErrorBound(nodes []float64, derivBound, a, b float64) float64 {
	n := len(nodes) // this is the order of the derivative
	fact := 1.0
	for k := 2; k <= n; k++ {
		fact *= float64(k)
	}
	return derivBound / fact * MaxNodePolynomial(nodes, a, b)
}

// Factorial returns k! as a float64.
func Factorial(k int) float64 {
	if k < 0 {
		return math.NaN()
	}
	res := 1.0
	for i := 2; i <= k; i++ {
		res *= float64(i)
	}
	return res
}

// ChebyshevTruncationError estimates the uniform error of truncating a
// Chebyshev series after index k as the sum of the absolute values of the
// discarded coefficients (a standard a-posteriori bound).
func ChebyshevTruncationError(s *ChebSeries, k int) float64 {
	return s.TailNorm(k)
}
