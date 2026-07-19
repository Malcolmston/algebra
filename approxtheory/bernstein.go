package approxtheory

import "math"

// Binomial returns the binomial coefficient C(n, k) as a float64, computed
// multiplicatively to avoid overflow for moderate arguments.
func Binomial(n, k int) float64 {
	if k < 0 || k > n {
		return 0
	}
	if k > n-k {
		k = n - k
	}
	res := 1.0
	for i := 0; i < k; i++ {
		res = res * float64(n-i) / float64(i+1)
	}
	return res
}

// BernsteinBasis evaluates the k-th Bernstein basis polynomial of degree n,
// b_{k,n}(t) = C(n,k) t^k (1-t)^{n-k}, at t in [0, 1].
func BernsteinBasis(k, n int, t float64) float64 {
	if k < 0 || k > n {
		return 0
	}
	return Binomial(n, k) * math.Pow(t, float64(k)) * math.Pow(1-t, float64(n-k))
}

// BernsteinBasisInterval evaluates the k-th Bernstein basis polynomial of
// degree n on the interval [a, b] at x, by mapping x to [0, 1].
func BernsteinBasisInterval(k, n int, x, a, b float64) float64 {
	t := (x - a) / (b - a)
	return BernsteinBasis(k, n, t)
}

// BernsteinPoly represents a polynomial in the Bernstein (Bezier) form with
// control values Coeffs on the interval [A, B]. The degree equals
// len(Coeffs)-1.
type BernsteinPoly struct {
	Coeffs []float64 // control points b_0..b_n
	A, B   float64
}

// NewBernsteinPoly builds a Bernstein-form polynomial with the given control
// values on [a, b].
func NewBernsteinPoly(coeffs []float64, a, b float64) *BernsteinPoly {
	c := make([]float64, len(coeffs))
	copy(c, coeffs)
	return &BernsteinPoly{Coeffs: c, A: a, B: b}
}

// Degree returns the polynomial degree of the Bernstein form.
func (p *BernsteinPoly) Degree() int {
	if len(p.Coeffs) == 0 {
		return 0
	}
	return len(p.Coeffs) - 1
}

// Eval evaluates the Bernstein polynomial at x in [A, B] using the numerically
// stable de Casteljau algorithm.
func (p *BernsteinPoly) Eval(x float64) float64 {
	n := len(p.Coeffs)
	if n == 0 {
		return 0
	}
	t := (x - p.A) / (p.B - p.A)
	beta := make([]float64, n)
	copy(beta, p.Coeffs)
	for r := 1; r < n; r++ {
		for i := 0; i < n-r; i++ {
			beta[i] = beta[i]*(1-t) + beta[i+1]*t
		}
	}
	return beta[0]
}

// EvalSlice evaluates the Bernstein polynomial at every point in xs.
func (p *BernsteinPoly) EvalSlice(xs []float64) []float64 {
	out := make([]float64, len(xs))
	for i, x := range xs {
		out[i] = p.Eval(x)
	}
	return out
}

// DeCasteljau evaluates a Bezier curve with the given control values at
// parameter t in [0, 1] using the de Casteljau algorithm. It is the
// interval-free core used by BernsteinPoly.Eval.
func DeCasteljau(control []float64, t float64) float64 {
	n := len(control)
	if n == 0 {
		return 0
	}
	beta := make([]float64, n)
	copy(beta, control)
	for r := 1; r < n; r++ {
		for i := 0; i < n-r; i++ {
			beta[i] = beta[i]*(1-t) + beta[i+1]*t
		}
	}
	return beta[0]
}

// BezierEval evaluates a Bezier curve (Bernstein polynomial on [0,1]) with the
// given control values at parameter t. It is an alias for DeCasteljau.
func BezierEval(control []float64, t float64) float64 {
	return DeCasteljau(control, t)
}

// BernsteinApprox returns the degree-n Bernstein polynomial approximation of f
// on [a, b], whose control values are the samples f(a + k/n (b-a)). The
// Bernstein operator converges uniformly to any continuous f.
func BernsteinApprox(f func(float64) float64, n int, a, b float64) *BernsteinPoly {
	if n < 0 {
		n = 0
	}
	coeffs := make([]float64, n+1)
	for k := 0; k <= n; k++ {
		x := a + float64(k)/float64(n)*(b-a)
		if n == 0 {
			x = a
		}
		coeffs[k] = f(x)
	}
	return &BernsteinPoly{Coeffs: coeffs, A: a, B: b}
}

// Elevate returns an equivalent Bernstein polynomial of degree n+1 (degree
// elevation), leaving the represented polynomial unchanged.
func (p *BernsteinPoly) Elevate() *BernsteinPoly {
	n := len(p.Coeffs) - 1
	if n < 0 {
		return p.Clone()
	}
	out := make([]float64, n+2)
	out[0] = p.Coeffs[0]
	out[n+1] = p.Coeffs[n]
	for i := 1; i <= n; i++ {
		alpha := float64(i) / float64(n+1)
		out[i] = alpha*p.Coeffs[i-1] + (1-alpha)*p.Coeffs[i]
	}
	return &BernsteinPoly{Coeffs: out, A: p.A, B: p.B}
}

// Derivative returns the Bernstein form of the derivative of p, correctly
// scaled for the domain [A, B].
func (p *BernsteinPoly) Derivative() *BernsteinPoly {
	n := len(p.Coeffs) - 1
	if n <= 0 {
		return &BernsteinPoly{Coeffs: []float64{0}, A: p.A, B: p.B}
	}
	out := make([]float64, n)
	scale := float64(n) / (p.B - p.A)
	for i := 0; i < n; i++ {
		out[i] = scale * (p.Coeffs[i+1] - p.Coeffs[i])
	}
	return &BernsteinPoly{Coeffs: out, A: p.A, B: p.B}
}

// Integral returns the definite integral of the Bernstein polynomial over its
// whole domain [A, B]. Every Bernstein basis of degree n integrates to
// (B-A)/(n+1), so the integral is the mean of the control values times (B-A).
func (p *BernsteinPoly) Integral() float64 {
	n := len(p.Coeffs)
	if n == 0 {
		return 0
	}
	var sum float64
	for _, c := range p.Coeffs {
		sum += c
	}
	return sum / float64(n) * (p.B - p.A)
}

// Clone returns a deep copy of the Bernstein polynomial.
func (p *BernsteinPoly) Clone() *BernsteinPoly {
	return NewBernsteinPoly(p.Coeffs, p.A, p.B)
}

// ToMonomial converts the Bernstein form (assumed on [0, 1]) into monomial
// coefficients in ascending order.
func (p *BernsteinPoly) ToMonomial() []float64 {
	n := len(p.Coeffs) - 1
	if n < 0 {
		return nil
	}
	// c_j = sum_{i=0}^{j} (-1)^{j-i} C(n,j) C(j,i) b_i
	out := make([]float64, n+1)
	for j := 0; j <= n; j++ {
		cnj := Binomial(n, j)
		var s float64
		for i := 0; i <= j; i++ {
			sign := 1.0
			if (j-i)%2 == 1 {
				sign = -1.0
			}
			s += sign * Binomial(j, i) * p.Coeffs[i]
		}
		out[j] = cnj * s
	}
	return out
}
