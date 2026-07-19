package approxtheory

import "math"

// Barycentric holds a polynomial interpolant in the second (true) barycentric
// form, storing the nodes, the values at those nodes and the barycentric
// weights.
type Barycentric struct {
	Xs      []float64
	Ys      []float64
	Weights []float64
}

// BarycentricWeights computes the barycentric weights w_j = 1 / prod_{k!=j}
// (x_j - x_k) for an arbitrary set of distinct nodes.
func BarycentricWeights(xs []float64) []float64 {
	n := len(xs)
	w := make([]float64, n)
	for j := 0; j < n; j++ {
		prod := 1.0
		for k := 0; k < n; k++ {
			if k == j {
				continue
			}
			prod *= xs[j] - xs[k]
		}
		w[j] = 1 / prod
	}
	return w
}

// ChebyshevBarycentricWeights returns the closed-form barycentric weights for
// the n+1 Chebyshev-Gauss-Lobatto nodes (ordered a..b as in ChebPoints):
// w_j = (-1)^j delta_j with the endpoints halved. These weights are equal in
// magnitude and independent of the interval.
func ChebyshevBarycentricWeights(n int) []float64 {
	w := make([]float64, n+1)
	for j := 0; j <= n; j++ {
		sign := 1.0
		if j%2 == 1 {
			sign = -1.0
		}
		d := 1.0
		if j == 0 || j == n {
			d = 0.5
		}
		w[j] = sign * d
	}
	return w
}

// NewBarycentric builds a barycentric interpolant through (xs, ys) with weights
// derived from the nodes. The xs must be distinct.
func NewBarycentric(xs, ys []float64) (*Barycentric, error) {
	if len(xs) == 0 {
		return nil, ErrEmptyInput
	}
	if len(xs) != len(ys) {
		return nil, ErrDimensionMismatch
	}
	xc := make([]float64, len(xs))
	yc := make([]float64, len(ys))
	copy(xc, xs)
	copy(yc, ys)
	return &Barycentric{Xs: xc, Ys: yc, Weights: BarycentricWeights(xs)}, nil
}

// NewBarycentricWeighted builds a barycentric interpolant with caller-supplied
// weights (for example the closed-form Chebyshev weights). The three slices
// must have equal length.
func NewBarycentricWeighted(xs, ys, weights []float64) (*Barycentric, error) {
	if len(xs) == 0 {
		return nil, ErrEmptyInput
	}
	if len(xs) != len(ys) || len(xs) != len(weights) {
		return nil, ErrDimensionMismatch
	}
	xc := make([]float64, len(xs))
	yc := make([]float64, len(ys))
	wc := make([]float64, len(weights))
	copy(xc, xs)
	copy(yc, ys)
	copy(wc, weights)
	return &Barycentric{Xs: xc, Ys: yc, Weights: wc}, nil
}

// Eval evaluates the barycentric interpolant at x. When x coincides with a
// node the exact node value is returned.
func (b *Barycentric) Eval(x float64) float64 {
	var num, den float64
	for j := range b.Xs {
		d := x - b.Xs[j]
		if d == 0 {
			return b.Ys[j]
		}
		t := b.Weights[j] / d
		num += t * b.Ys[j]
		den += t
	}
	return num / den
}

// EvalSlice evaluates the interpolant at every point in xs.
func (b *Barycentric) EvalSlice(xs []float64) []float64 {
	out := make([]float64, len(xs))
	for i, x := range xs {
		out[i] = b.Eval(x)
	}
	return out
}

// ChebyshevInterpolant builds a barycentric polynomial interpolant of f at the
// n+1 Chebyshev-Gauss-Lobatto nodes on [a, b], using the stable closed-form
// weights. This is the recommended way to interpolate at Chebyshev points.
func ChebyshevInterpolant(f func(float64) float64, n int, a, b float64) *Barycentric {
	xs := ChebPoints(n, a, b)
	ys := make([]float64, len(xs))
	for i, x := range xs {
		ys[i] = f(x)
	}
	return &Barycentric{Xs: xs, Ys: ys, Weights: ChebyshevBarycentricWeights(n)}
}

// LagrangeEval evaluates, at x, the Lagrange interpolating polynomial through
// (xs, ys) directly (O(n^2) per point). It is a simple reference companion to
// the barycentric form.
func LagrangeEval(xs, ys []float64, x float64) float64 {
	n := len(xs)
	var sum float64
	for j := 0; j < n; j++ {
		term := ys[j]
		for k := 0; k < n; k++ {
			if k == j {
				continue
			}
			term *= (x - xs[k]) / (xs[j] - xs[k])
		}
		sum += term
	}
	return sum
}

// LagrangeBasis evaluates the j-th Lagrange cardinal polynomial for the nodes
// xs at x.
func LagrangeBasis(xs []float64, j int, x float64) float64 {
	term := 1.0
	for k := range xs {
		if k == j {
			continue
		}
		term *= (x - xs[k]) / (xs[j] - xs[k])
	}
	return term
}

// NewtonDividedDifferences returns the divided-difference coefficients of the
// Newton interpolating polynomial through (xs, ys).
func NewtonDividedDifferences(xs, ys []float64) []float64 {
	n := len(xs)
	c := make([]float64, n)
	copy(c, ys)
	for j := 1; j < n; j++ {
		for i := n - 1; i >= j; i-- {
			c[i] = (c[i] - c[i-1]) / (xs[i] - xs[i-j])
		}
	}
	return c
}

// NewtonEval evaluates the Newton form with the given divided-difference
// coefficients and nodes at x using Horner-like nesting.
func NewtonEval(coeffs, xs []float64, x float64) float64 {
	n := len(coeffs)
	if n == 0 {
		return math.NaN()
	}
	acc := coeffs[n-1]
	for i := n - 2; i >= 0; i-- {
		acc = acc*(x-xs[i]) + coeffs[i]
	}
	return acc
}
