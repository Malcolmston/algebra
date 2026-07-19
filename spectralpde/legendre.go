package spectralpde

import "math"

// LegendreP evaluates the Legendre polynomial P_n at x via the three-term
// recurrence.
func LegendreP(n int, x float64) float64 {
	if n == 0 {
		return 1
	}
	if n == 1 {
		return x
	}
	pm1, p := 1.0, x
	for k := 2; k <= n; k++ {
		pm1, p = p, ((2*float64(k)-1)*x*p-(float64(k)-1)*pm1)/float64(k)
	}
	return p
}

// LegendrePDerivative evaluates the derivative P_n'(x).
func LegendrePDerivative(n int, x float64) float64 {
	if n == 0 {
		return 0
	}
	if x == 1 {
		return float64(n) * float64(n+1) / 2
	}
	if x == -1 {
		s := 1.0
		if (n+1)%2 == 1 {
			s = -1.0
		}
		return s * float64(n) * float64(n+1) / 2
	}
	pn := LegendreP(n, x)
	pnm1 := LegendreP(n-1, x)
	return float64(n) * (pnm1 - x*pn) / (1 - x*x)
}

// LegendrePSecondDerivative evaluates P_n”(x) using the Legendre differential
// equation.
func LegendrePSecondDerivative(n int, x float64) float64 {
	if x == 1 || x == -1 {
		if n < 2 {
			return 0
		}
		// Closed form P_n^{(2)}(1) = (n+2)(n+1)n(n-1)/8, with the value at -1
		// carrying an extra factor (-1)^n.
		nf := float64(n)
		val := (nf + 2) * (nf + 1) * nf * (nf - 1) / 8
		if x == -1 && n%2 == 1 {
			val = -val
		}
		return val
	}
	pn := LegendreP(n, x)
	pp := LegendrePDerivative(n, x)
	return (2*x*pp - float64(n)*float64(n+1)*pn) / (1 - x*x)
}

// LegendrePValues evaluates P_n at every point of xs.
func LegendrePValues(n int, xs []float64) []float64 {
	out := make([]float64, len(xs))
	for i, x := range xs {
		out[i] = LegendreP(n, x)
	}
	return out
}

// LegendreVandermonde returns the len(nodes) x (degree+1) matrix with entries
// P_j(nodes[i]).
func LegendreVandermonde(nodes []float64, degree int) [][]float64 {
	m := len(nodes)
	v := Zeros(m, degree+1)
	for i, x := range nodes {
		v[i][0] = 1
		if degree >= 1 {
			v[i][1] = x
		}
		for j := 2; j <= degree; j++ {
			v[i][j] = ((2*float64(j)-1)*x*v[i][j-1] - (float64(j)-1)*v[i][j-2]) / float64(j)
		}
	}
	return v
}

// LegendreGaussNodesWeights returns the n Gauss-Legendre nodes and weights on
// [-1, 1], computed by Newton iteration on P_n. Nodes are returned in
// ascending order.
func LegendreGaussNodesWeights(n int) ([]float64, []float64) {
	nodes := make([]float64, n)
	weights := make([]float64, n)
	for k := 0; k < n; k++ {
		x := math.Cos(math.Pi * (float64(k) + 0.75) / (float64(n) + 0.5))
		for it := 0; it < 100; it++ {
			p := LegendreP(n, x)
			dp := LegendrePDerivative(n, x)
			dx := -p / dp
			x += dx
			if math.Abs(dx) < 1e-15 {
				break
			}
		}
		dp := LegendrePDerivative(n, x)
		// Store ascending: index from the end since initial guesses descend.
		idx := n - 1 - k
		nodes[idx] = x
		weights[idx] = 2 / ((1 - x*x) * dp * dp)
	}
	return nodes, weights
}

// LegendreGaussNodes returns just the n Gauss-Legendre nodes on [-1, 1].
func LegendreGaussNodes(n int) []float64 {
	nodes, _ := LegendreGaussNodesWeights(n)
	return nodes
}

// LegendreGaussWeights returns just the n Gauss-Legendre weights on [-1, 1].
func LegendreGaussWeights(n int) []float64 {
	_, w := LegendreGaussNodesWeights(n)
	return w
}

// LegendreGaussLobattoNodesWeights returns the N+1 Legendre-Gauss-Lobatto
// nodes and weights on [-1, 1] (endpoints included). Nodes are in ascending
// order. N must be >= 1.
func LegendreGaussLobattoNodesWeights(N int) ([]float64, []float64) {
	nodes := make([]float64, N+1)
	weights := make([]float64, N+1)
	nodes[0] = -1
	nodes[N] = 1
	// Interior nodes are the roots of P_N'.
	for k := 1; k < N; k++ {
		x := math.Cos(math.Pi * float64(k) / float64(N))
		for it := 0; it < 100; it++ {
			dp := LegendrePDerivative(N, x)
			ddp := LegendrePSecondDerivative(N, x)
			dx := -dp / ddp
			x += dx
			if math.Abs(dx) < 1e-15 {
				break
			}
		}
		nodes[N-k] = x
	}
	for i := 0; i <= N; i++ {
		pn := LegendreP(N, nodes[i])
		weights[i] = 2 / (float64(N) * float64(N+1) * pn * pn)
	}
	return nodes, weights
}

// LegendreGaussLobattoNodes returns just the N+1 Legendre-Gauss-Lobatto nodes.
func LegendreGaussLobattoNodes(N int) []float64 {
	nodes, _ := LegendreGaussLobattoNodesWeights(N)
	return nodes
}

// LegendreGaussLobattoWeights returns just the N+1 Legendre-Gauss-Lobatto
// weights.
func LegendreGaussLobattoWeights(N int) []float64 {
	_, w := LegendreGaussLobattoNodesWeights(N)
	return w
}

// LegendreSeriesEval evaluates the Legendre series sum_k coeffs[k]*P_k at x.
func LegendreSeriesEval(coeffs []float64, x float64) float64 {
	var s float64
	for k := 0; k < len(coeffs); k++ {
		s += coeffs[k] * LegendreP(k, x)
	}
	return s
}

// LegendreProjection returns the first n+1 Legendre coefficients of f, i.e.
// the L2([-1,1]) projection onto span{P_0,...,P_n}, computed with an accurate
// Gauss-Legendre rule. c_k = (2k+1)/2 * integral_{-1}^{1} f(x) P_k(x) dx.
func LegendreProjection(f func(float64) float64, n int) []float64 {
	q := n + 2
	nodes, weights := LegendreGaussNodesWeights(q)
	coeffs := make([]float64, n+1)
	for k := 0; k <= n; k++ {
		var s float64
		for i := range nodes {
			s += weights[i] * f(nodes[i]) * LegendreP(k, nodes[i])
		}
		coeffs[k] = (2*float64(k) + 1) / 2 * s
	}
	return coeffs
}

// LegendreNormSquared returns the squared L2([-1,1]) norm of P_n, 2/(2n+1).
func LegendreNormSquared(n int) float64 {
	return 2 / (2*float64(n) + 1)
}
