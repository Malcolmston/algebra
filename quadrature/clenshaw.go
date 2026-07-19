package quadrature

import "math"

// ClenshawCurtisWeights returns the n+1 Clenshaw-Curtis weights on [-1, 1] for
// the nodes cos(j*pi/n), j = 0..n (Trefethen's algorithm). The nodes are the
// Chebyshev extrema and the weights sum to 2. n must be at least 1.
func ClenshawCurtisWeights(n int) (nodes, weights []float64) {
	nodes = make([]float64, n+1)
	weights = make([]float64, n+1)
	theta := make([]float64, n+1)
	for j := 0; j <= n; j++ {
		theta[j] = math.Pi * float64(j) / float64(n)
		nodes[j] = math.Cos(theta[j])
	}
	for j := 1; j < n; j++ {
		v := 1.0
		if n%2 == 0 {
			for k := 1; k <= n/2-1; k++ {
				v -= 2 * math.Cos(2*float64(k)*theta[j]) / float64(4*k*k-1)
			}
			v -= math.Cos(float64(n)*theta[j]) / float64(n*n-1)
		} else {
			for k := 1; k <= (n-1)/2; k++ {
				v -= 2 * math.Cos(2*float64(k)*theta[j]) / float64(4*k*k-1)
			}
		}
		weights[j] = 2 * v / float64(n)
	}
	if n%2 == 0 {
		weights[0] = 1 / float64(n*n-1)
	} else {
		weights[0] = 1 / float64(n*n)
	}
	weights[n] = weights[0]
	// Reverse so that nodes come out ascending (cos is decreasing in j).
	for i, j := 0, n; i < j; i, j = i+1, j-1 {
		nodes[i], nodes[j] = nodes[j], nodes[i]
		weights[i], weights[j] = weights[j], weights[i]
	}
	return nodes, weights
}

// ClenshawCurtisRule returns the (n+1)-point Clenshaw-Curtis rule on [-1, 1].
func ClenshawCurtisRule(n int) Rule {
	nodes, weights := ClenshawCurtisWeights(n)
	return Rule{Nodes: nodes, Weights: weights}
}

// ClenshawCurtis approximates the integral of f over [a, b] with the
// (n+1)-point Clenshaw-Curtis rule.
func ClenshawCurtis(f Func, a, b float64, n int) float64 {
	nodes, weights := ClenshawCurtisWeights(n)
	half := 0.5 * (b - a)
	mid := 0.5 * (a + b)
	var s float64
	for i, t := range nodes {
		s += weights[i] * f(mid+half*t)
	}
	return half * s
}

// Fejer1Weights returns the n nodes and weights of Fejer's first rule on
// [-1, 1], an open Chebyshev rule with nodes cos((2i-1)pi/(2n)). The weights
// sum to 2.
func Fejer1Weights(n int) (nodes, weights []float64) {
	nodes = make([]float64, n)
	weights = make([]float64, n)
	for i := 1; i <= n; i++ {
		theta := float64(2*i-1) * math.Pi / float64(2*n)
		v := 1.0
		for m := 1; m <= n/2; m++ {
			v -= 2 * math.Cos(2*float64(m)*theta) / float64(4*m*m-1)
		}
		w := 2 * v / float64(n)
		nodes[n-i] = math.Cos(theta)
		weights[n-i] = w
	}
	return nodes, weights
}

// Fejer1Rule returns Fejer's first rule of n points on [-1, 1].
func Fejer1Rule(n int) Rule {
	nodes, weights := Fejer1Weights(n)
	return Rule{Nodes: nodes, Weights: weights}
}

// IntegrateFejer1 approximates the integral of f over [a, b] with Fejer's
// first rule of n points.
func IntegrateFejer1(f Func, a, b float64, n int) float64 {
	nodes, weights := Fejer1Weights(n)
	half := 0.5 * (b - a)
	mid := 0.5 * (a + b)
	var s float64
	for i, t := range nodes {
		s += weights[i] * f(mid+half*t)
	}
	return half * s
}

// Fejer2Weights returns the n nodes and weights of Fejer's second rule on
// [-1, 1], an open Chebyshev rule with nodes cos(i*pi/(n+1)), i = 1..n. The
// weights sum to 2.
func Fejer2Weights(n int) (nodes, weights []float64) {
	nodes = make([]float64, n)
	weights = make([]float64, n)
	np1 := n + 1
	for i := 1; i <= n; i++ {
		theta := float64(i) * math.Pi / float64(np1)
		s := math.Sin(theta)
		var sum float64
		for m := 1; m <= (np1)/2; m++ {
			sum += math.Sin(float64(2*m-1)*theta) / float64(2*m-1)
		}
		w := 4 * s / float64(np1) * sum
		nodes[n-i] = math.Cos(theta)
		weights[n-i] = w
	}
	return nodes, weights
}

// Fejer2Rule returns Fejer's second rule of n points on [-1, 1].
func Fejer2Rule(n int) Rule {
	nodes, weights := Fejer2Weights(n)
	return Rule{Nodes: nodes, Weights: weights}
}

// IntegrateFejer2 approximates the integral of f over [a, b] with Fejer's
// second rule of n points.
func IntegrateFejer2(f Func, a, b float64, n int) float64 {
	nodes, weights := Fejer2Weights(n)
	half := 0.5 * (b - a)
	mid := 0.5 * (a + b)
	var s float64
	for i, t := range nodes {
		s += weights[i] * f(mid+half*t)
	}
	return half * s
}
