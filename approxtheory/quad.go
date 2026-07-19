package approxtheory

import "math"

// GaussChebyshevNodes returns the n Gauss-Chebyshev nodes (roots of T_n) on
// [a, b] in ascending order. These are the nodes of Gauss-Chebyshev quadrature
// of the first kind.
func GaussChebyshevNodes(n int, a, b float64) []float64 {
	return ChebGaussPoints(n, a, b)
}

// GaussChebyshevWeights returns the n weights of Gauss-Chebyshev quadrature of
// the first kind on [a, b]. Each weight equals pi/n scaled by (b-a)/2. The
// quadrature integrates f(x) against the weight 1/sqrt(1-t^2), so it estimates
// the integral of f divided by that weight; see GaussChebyshevQuadrature for
// the plain integral.
func GaussChebyshevWeights(n int, a, b float64) []float64 {
	w := make([]float64, n)
	val := math.Pi / float64(n)
	for i := range w {
		w[i] = val
	}
	return w
}

// GaussChebyshevQuadrature approximates the integral of f over [a, b] using
// n-point Gauss-Chebyshev quadrature of the first kind, folding the Chebyshev
// weight back in so the plain (unweighted) integral is returned.
func GaussChebyshevQuadrature(f func(float64) float64, n int, a, b float64) float64 {
	nodes := ChebGaussPoints(n, a, b)
	var sum float64
	for i := 0; i < n; i++ {
		theta := math.Pi * (float64(i) + 0.5) / float64(n)
		// weight 1/sqrt(1-t^2) is folded out by multiplying by sin(theta).
		s := math.Sin(theta)
		sum += f(nodes[i]) * s
	}
	return sum * math.Pi / float64(n) * (b - a) / 2
}

// ClenshawCurtisWeights returns the n+1 Clenshaw-Curtis quadrature weights on
// [a, b] associated with the Chebyshev-Gauss-Lobatto nodes (ordered a..b as in
// ChebPoints). The weights integrate polynomials up to degree n exactly.
func ClenshawCurtisWeights(n int, a, b float64) []float64 {
	if n == 0 {
		return []float64{b - a}
	}
	w := make([]float64, n+1)
	// Standard Clenshaw-Curtis weights on [-1,1], symmetric in node index.
	for k := 0; k <= n; k++ {
		var sum float64
		// c_0 term
		sum = 1.0
		jmax := n / 2
		for j := 1; j <= jmax; j++ {
			var bj float64
			if 2*j == n {
				bj = 1
			} else {
				bj = 2
			}
			sum += bj / (1 - float64(4*j*j)) * math.Cos(2*float64(j)*float64(k)*math.Pi/float64(n))
		}
		ck := 2.0
		if k == 0 || k == n {
			ck = 1.0
		}
		w[k] = ck / float64(n) * sum
	}
	// scale to [a,b]
	scale := (b - a) / 2
	for i := range w {
		w[i] *= scale
	}
	return w
}

// ClenshawCurtisQuadrature approximates the integral of f over [a, b] using
// n+1 point Clenshaw-Curtis quadrature at the Chebyshev-Gauss-Lobatto nodes.
func ClenshawCurtisQuadrature(f func(float64) float64, n int, a, b float64) float64 {
	nodes := ChebPoints(n, a, b)
	w := ClenshawCurtisWeights(n, a, b)
	var sum float64
	for i := 0; i <= n; i++ {
		sum += w[i] * f(nodes[i])
	}
	return sum
}

// FejerNodes returns the n Fejer nodes (Chebyshev points of the first kind,
// the roots of T_n) on [a, b], ascending.
func FejerNodes(n int, a, b float64) []float64 {
	return ChebGaussPoints(n, a, b)
}

// FejerWeights returns the n weights of Fejer's first quadrature rule on
// [a, b], associated with the roots of T_n.
func FejerWeights(n int, a, b float64) []float64 {
	w := make([]float64, n)
	for k := 0; k < n; k++ {
		theta := math.Pi * (float64(n-1-k) + 0.5) / float64(n)
		var sum float64
		jmax := n / 2
		for j := 1; j <= jmax; j++ {
			sum += math.Cos(2*float64(j)*theta) / (4*float64(j)*float64(j) - 1)
		}
		w[k] = (2.0 / float64(n)) * (1 - 2*sum)
	}
	scale := (b - a) / 2
	for i := range w {
		w[i] *= scale
	}
	return w
}

// FejerQuadrature approximates the integral of f over [a, b] using Fejer's
// first rule with n nodes.
func FejerQuadrature(f func(float64) float64, n int, a, b float64) float64 {
	nodes := FejerNodes(n, a, b)
	w := FejerWeights(n, a, b)
	var sum float64
	for i := 0; i < n; i++ {
		sum += w[i] * f(nodes[i])
	}
	return sum
}

// ChebyshevQuadrature integrates f over the domain of a Chebyshev fit by first
// forming the degree-n interpolant and integrating it exactly. This is a
// spectrally accurate alternative to the explicit quadrature rules.
func ChebyshevQuadrature(f func(float64) float64, n int, a, b float64) float64 {
	return ChebFit(f, n, a, b).Integral()
}
