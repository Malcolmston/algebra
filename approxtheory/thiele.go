package approxtheory

import "math"

// Thiele holds a rational interpolant built from Thiele's continued fraction.
// Xs are the abscissae and Rho are the continued-fraction coefficients
// (the leading reciprocal differences).
type Thiele struct {
	Xs  []float64
	Rho []float64
}

// NewThiele builds Thiele's continued-fraction rational interpolant through
// the points (xs, ys). It computes the reciprocal-difference table and returns
// an interpolant that reproduces every data value. The xs must be distinct.
func NewThiele(xs, ys []float64) (*Thiele, error) {
	n := len(xs)
	if n == 0 {
		return nil, ErrEmptyInput
	}
	if len(ys) != n {
		return nil, ErrDimensionMismatch
	}
	// y[i][k] is the inverse difference conditioned on the path x_0..x_{k-1}:
	//   y[i][0] = f(x_i)
	//   y[i][k] = (x_i - x_{k-1}) / (y[i][k-1] - y[k-1][k-1]).
	// The Thiele continued-fraction coefficients are a_k = y[k][k].
	y := make([][]float64, n)
	for i := range y {
		y[i] = make([]float64, n)
		y[i][0] = ys[i]
	}
	for k := 1; k < n; k++ {
		for i := k; i < n; i++ {
			denom := y[i][k-1] - y[k-1][k-1]
			if denom == 0 {
				denom = 1e-300 // guard against exact repeats / poles
			}
			y[i][k] = (xs[i] - xs[k-1]) / denom
		}
	}
	coeffs := make([]float64, n)
	for k := 0; k < n; k++ {
		coeffs[k] = y[k][k]
	}
	xcopy := make([]float64, n)
	copy(xcopy, xs)
	return &Thiele{Xs: xcopy, Rho: coeffs}, nil
}

// Eval evaluates the Thiele continued-fraction interpolant at x.
func (t *Thiele) Eval(x float64) float64 {
	n := len(t.Rho)
	if n == 0 {
		return math.NaN()
	}
	var acc float64
	for i := n - 1; i >= 1; i-- {
		denom := t.Rho[i] + acc
		if denom == 0 {
			denom = 1e-300
		}
		acc = (x - t.Xs[i-1]) / denom
	}
	return t.Rho[0] + acc
}

// EvalSlice evaluates the interpolant at every point in xs.
func (t *Thiele) EvalSlice(xs []float64) []float64 {
	out := make([]float64, len(xs))
	for i, x := range xs {
		out[i] = t.Eval(x)
	}
	return out
}

// Degree returns the number of interpolation nodes, one more than the depth of
// the continued fraction.
func (t *Thiele) Degree() int { return len(t.Xs) }

// ThieleInterp is a convenience wrapper that builds a Thiele interpolant from
// (xs, ys) and evaluates it at a single point x.
func ThieleInterp(xs, ys []float64, x float64) (float64, error) {
	t, err := NewThiele(xs, ys)
	if err != nil {
		return 0, err
	}
	return t.Eval(x), nil
}
