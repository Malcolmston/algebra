package spectralpde

import "math"

// BarycentricWeights returns the barycentric interpolation weights
// w_j = 1 / prod_{k != j} (nodes[j] - nodes[k]) for an arbitrary set of
// distinct nodes.
func BarycentricWeights(nodes []float64) []float64 {
	n := len(nodes)
	w := make([]float64, n)
	for j := 0; j < n; j++ {
		p := 1.0
		for k := 0; k < n; k++ {
			if k != j {
				p *= nodes[j] - nodes[k]
			}
		}
		w[j] = 1 / p
	}
	return w
}

// BarycentricWeightsChebyshev returns the closed-form barycentric weights for
// the N+1 Chebyshev-Gauss-Lobatto nodes: w_j = (-1)^j * delta_j with
// delta_j = 1/2 at the endpoints and 1 in the interior.
func BarycentricWeightsChebyshev(N int) []float64 {
	w := make([]float64, N+1)
	for j := 0; j <= N; j++ {
		s := 1.0
		if j%2 == 1 {
			s = -1.0
		}
		if j == 0 || j == N {
			s *= 0.5
		}
		w[j] = s
	}
	return w
}

// BarycentricInterpolate evaluates the interpolating polynomial through
// (nodes, values) at x, using the barycentric weights w (see
// BarycentricWeights).
func BarycentricInterpolate(nodes, values, w []float64, x float64) float64 {
	var num, den float64
	for j := 0; j < len(nodes); j++ {
		d := x - nodes[j]
		if d == 0 {
			return values[j]
		}
		t := w[j] / d
		num += t * values[j]
		den += t
	}
	return num / den
}

// BarycentricInterpolateVec evaluates the barycentric interpolant at every
// point of xs.
func BarycentricInterpolateVec(nodes, values, w, xs []float64) []float64 {
	out := make([]float64, len(xs))
	for i, x := range xs {
		out[i] = BarycentricInterpolate(nodes, values, w, x)
	}
	return out
}

// LagrangeBasis evaluates the j-th Lagrange cardinal polynomial for the given
// nodes at x.
func LagrangeBasis(nodes []float64, j int, x float64) float64 {
	p := 1.0
	for k := 0; k < len(nodes); k++ {
		if k != j {
			p *= (x - nodes[k]) / (nodes[j] - nodes[k])
		}
	}
	return p
}

// LagrangeInterpolate evaluates the Lagrange interpolating polynomial through
// (nodes, values) at x directly from the cardinal-polynomial definition.
func LagrangeInterpolate(nodes, values []float64, x float64) float64 {
	var s float64
	for j := 0; j < len(nodes); j++ {
		s += values[j] * LagrangeBasis(nodes, j, x)
	}
	return s
}

// PolynomialInterpolate builds barycentric weights for the given nodes and
// evaluates the interpolant at x. It is convenient for one-off evaluations.
func PolynomialInterpolate(nodes, values []float64, x float64) float64 {
	w := BarycentricWeights(nodes)
	return BarycentricInterpolate(nodes, values, w, x)
}

// ChebyshevBarycentricInterpolate evaluates the interpolant of values sampled
// at the Chebyshev-Gauss-Lobatto nodes at x, using the closed-form Chebyshev
// barycentric weights.
func ChebyshevBarycentricInterpolate(values []float64, x float64) float64 {
	N := len(values) - 1
	nodes := ChebyshevGaussLobattoNodes(N)
	w := BarycentricWeightsChebyshev(N)
	return BarycentricInterpolate(nodes, values, w, x)
}

// NewtonDividedDifferences returns the Newton divided-difference coefficients
// for the data (nodes, values).
func NewtonDividedDifferences(nodes, values []float64) []float64 {
	n := len(values)
	coef := make([]float64, n)
	copy(coef, values)
	for j := 1; j < n; j++ {
		for i := n - 1; i >= j; i-- {
			coef[i] = (coef[i] - coef[i-1]) / (nodes[i] - nodes[i-j])
		}
	}
	return coef
}

// NewtonEval evaluates the Newton form of the interpolating polynomial with the
// given divided-difference coefficients and nodes at x.
func NewtonEval(coef, nodes []float64, x float64) float64 {
	n := len(coef)
	if n == 0 {
		return 0
	}
	res := coef[n-1]
	for i := n - 2; i >= 0; i-- {
		res = res*(x-nodes[i]) + coef[i]
	}
	return res
}

// InterpolationError returns the maximum absolute difference between f and its
// interpolant (through (nodes, f(nodes))) sampled at the given evaluation
// points.
func InterpolationError(f func(float64) float64, nodes []float64, eval []float64) float64 {
	values := ApplyFunc(f, nodes)
	w := BarycentricWeights(nodes)
	var maxErr float64
	for _, x := range eval {
		e := math.Abs(f(x) - BarycentricInterpolate(nodes, values, w, x))
		if e > maxErr {
			maxErr = e
		}
	}
	return maxErr
}
