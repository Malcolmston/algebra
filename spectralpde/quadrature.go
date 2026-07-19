package spectralpde

import "math"

// ClenshawCurtisNodes returns the N+1 Clenshaw-Curtis nodes, i.e. the
// Chebyshev-Gauss-Lobatto nodes on [-1, 1].
func ClenshawCurtisNodes(N int) []float64 {
	return ChebyshevGaussLobattoNodes(N)
}

// ClenshawCurtisWeights returns the N+1 Clenshaw-Curtis quadrature weights on
// [-1, 1] for the weight function 1 (Trefethen's clencurt algorithm). N must
// be >= 1.
func ClenshawCurtisWeights(N int) []float64 {
	w := make([]float64, N+1)
	theta := ChebyshevGaussLobattoAngles(N)
	if N%2 == 0 {
		w[0] = 1.0 / float64(N*N-1)
		w[N] = w[0]
		for i := 1; i < N; i++ {
			v := 1.0
			for k := 1; k <= N/2-1; k++ {
				v -= 2 * math.Cos(2*float64(k)*theta[i]) / float64(4*k*k-1)
			}
			v -= math.Cos(float64(N)*theta[i]) / float64(N*N-1)
			w[i] = 2 * v / float64(N)
		}
	} else {
		w[0] = 1.0 / float64(N*N)
		w[N] = w[0]
		for i := 1; i < N; i++ {
			v := 1.0
			for k := 1; k <= (N-1)/2; k++ {
				v -= 2 * math.Cos(2*float64(k)*theta[i]) / float64(4*k*k-1)
			}
			w[i] = 2 * v / float64(N)
		}
	}
	return w
}

// ClenshawCurtisIntegrate approximates the integral of f over [a, b] using an
// (N+1)-point Clenshaw-Curtis rule.
func ClenshawCurtisIntegrate(f func(float64) float64, N int, a, b float64) float64 {
	nodes := MapToInterval(ClenshawCurtisNodes(N), a, b)
	w := ClenshawCurtisWeights(N)
	scale := 0.5 * (b - a)
	var s float64
	for i := range nodes {
		s += w[i] * f(nodes[i])
	}
	return scale * s
}

// FejerNodes returns the n nodes of Fejér's first quadrature rule, the
// Chebyshev-Gauss nodes x_k = cos((k+1/2)*pi/n).
func FejerNodes(n int) []float64 {
	return ChebyshevGaussNodes(n)
}

// FejerWeights returns the n weights of Fejér's first quadrature rule on
// [-1, 1] for the weight function 1.
func FejerWeights(n int) []float64 {
	w := make([]float64, n)
	half := n / 2
	for k := 0; k < n; k++ {
		theta := (float64(k) + 0.5) * math.Pi / float64(n)
		s := 1.0
		for j := 1; j <= half; j++ {
			s -= 2 * math.Cos(2*float64(j)*theta) / float64(4*j*j-1)
		}
		w[k] = 2 * s / float64(n)
	}
	return w
}

// FejerIntegrate approximates the integral of f over [a, b] with an n-point
// Fejér first rule.
func FejerIntegrate(f func(float64) float64, n int, a, b float64) float64 {
	nodes := MapToInterval(FejerNodes(n), a, b)
	w := FejerWeights(n)
	scale := 0.5 * (b - a)
	var s float64
	for i := range nodes {
		s += w[i] * f(nodes[i])
	}
	return scale * s
}

// GaussLegendreIntegrate approximates the integral of f over [a, b] using an
// n-point Gauss-Legendre rule.
func GaussLegendreIntegrate(f func(float64) float64, n int, a, b float64) float64 {
	nodes, weights := LegendreGaussNodesWeights(n)
	mapped := MapToInterval(nodes, a, b)
	scale := 0.5 * (b - a)
	var s float64
	for i := range mapped {
		s += weights[i] * f(mapped[i])
	}
	return scale * s
}

// GaussChebyshevIntegrate approximates the integral over [-1, 1] of f(x) with
// respect to the Chebyshev weight 1/sqrt(1-x^2), using an n-point
// Gauss-Chebyshev rule.
func GaussChebyshevIntegrate(f func(float64) float64, n int) float64 {
	nodes := ChebyshevGaussNodes(n)
	w := math.Pi / float64(n)
	var s float64
	for _, x := range nodes {
		s += f(x)
	}
	return w * s
}

// IntegrateWeighted returns the dot product of quadrature weights and nodal
// function values.
func IntegrateWeighted(weights, values []float64) float64 {
	return DotProduct(weights, values)
}

// IntegrateFunction evaluates f at the given nodes and contracts with the
// weights.
func IntegrateFunction(f func(float64) float64, nodes, weights []float64) float64 {
	var s float64
	for i := range nodes {
		s += weights[i] * f(nodes[i])
	}
	return s
}

// TrapezoidalWeights returns the composite trapezoidal weights for n+1 equally
// spaced points on [a, b].
func TrapezoidalWeights(n int, a, b float64) []float64 {
	h := (b - a) / float64(n)
	w := make([]float64, n+1)
	for i := 0; i <= n; i++ {
		if i == 0 || i == n {
			w[i] = h / 2
		} else {
			w[i] = h
		}
	}
	return w
}

// TrapezoidalIntegrate approximates the integral of f over [a, b] with the
// composite trapezoidal rule on n subintervals.
func TrapezoidalIntegrate(f func(float64) float64, n int, a, b float64) float64 {
	x := LinSpace(a, b, n+1)
	w := TrapezoidalWeights(n, a, b)
	var s float64
	for i := range x {
		s += w[i] * f(x[i])
	}
	return s
}

// SimpsonIntegrate approximates the integral of f over [a, b] with the
// composite Simpson rule on n subintervals (n must be even and >= 2).
func SimpsonIntegrate(f func(float64) float64, n int, a, b float64) float64 {
	if n < 2 || n%2 != 0 {
		n = n + (n % 2)
		if n < 2 {
			n = 2
		}
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
	return s * h / 3
}
