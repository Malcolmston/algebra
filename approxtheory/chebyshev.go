package approxtheory

import (
	"math"
	"sort"
)

// ChebSeries represents a finite Chebyshev series
//
//	f(x) = sum_{k=0}^{n} Coeffs[k] * T_k(t),   t = (2x-(A+B))/(B-A),
//
// on the interval [A, B]. T_k is the Chebyshev polynomial of the first kind.
type ChebSeries struct {
	Coeffs []float64 // Chebyshev coefficients, index k multiplies T_k
	A, B   float64   // domain endpoints, A < B
}

// NewChebSeries builds a ChebSeries from raw coefficients on the interval
// [a, b]. The coefficient slice is copied.
func NewChebSeries(coeffs []float64, a, b float64) *ChebSeries {
	c := make([]float64, len(coeffs))
	copy(c, coeffs)
	return &ChebSeries{Coeffs: c, A: a, B: b}
}

// mapToUnit maps x in [A, B] to t in [-1, 1].
func (s *ChebSeries) mapToUnit(x float64) float64 {
	return (2*x - (s.A + s.B)) / (s.B - s.A)
}

// mapFromUnit maps t in [-1, 1] to x in [A, B].
func (s *ChebSeries) mapFromUnit(t float64) float64 {
	return 0.5*(s.A+s.B) + 0.5*(s.B-s.A)*t
}

// ChebT evaluates the Chebyshev polynomial of the first kind T_n at x in
// [-1, 1] using the stable three-term recurrence. For |x| > 1 the polynomial
// recurrence is still applied.
func ChebT(n int, x float64) float64 {
	if n < 0 {
		return 0
	}
	if n == 0 {
		return 1
	}
	tm1, t := 1.0, x
	for k := 2; k <= n; k++ {
		tm1, t = t, 2*x*t-tm1
	}
	return t
}

// ChebU evaluates the Chebyshev polynomial of the second kind U_n at x.
func ChebU(n int, x float64) float64 {
	if n < 0 {
		return 0
	}
	if n == 0 {
		return 1
	}
	um1, u := 1.0, 2*x
	for k := 2; k <= n; k++ {
		um1, u = u, 2*x*u-um1
	}
	return u
}

// ChebTValues returns the vector [T_0(x), T_1(x), ..., T_n(x)].
func ChebTValues(n int, x float64) []float64 {
	out := make([]float64, n+1)
	if n < 0 {
		return nil
	}
	out[0] = 1
	if n >= 1 {
		out[1] = x
	}
	for k := 2; k <= n; k++ {
		out[k] = 2*x*out[k-1] - out[k-2]
	}
	return out
}

// ChebClenshaw evaluates the Chebyshev series with coefficients c (index k
// multiplies T_k) at t in [-1, 1] using Clenshaw's recurrence.
func ChebClenshaw(c []float64, t float64) float64 {
	n := len(c)
	if n == 0 {
		return 0
	}
	if n == 1 {
		return c[0]
	}
	var bk1, bk2 float64
	for k := n - 1; k >= 1; k-- {
		bk1, bk2 = 2*t*bk1-bk2+c[k], bk1
	}
	return t*bk1 - bk2 + c[0]
}

// Eval returns the value of the series at x in [A, B].
func (s *ChebSeries) Eval(x float64) float64 {
	return ChebClenshaw(s.Coeffs, s.mapToUnit(x))
}

// EvalSlice evaluates the series at every point in xs.
func (s *ChebSeries) EvalSlice(xs []float64) []float64 {
	out := make([]float64, len(xs))
	for i, x := range xs {
		out[i] = s.Eval(x)
	}
	return out
}

// Degree returns the polynomial degree of the series (len(Coeffs)-1), or 0
// when empty.
func (s *ChebSeries) Degree() int {
	if len(s.Coeffs) == 0 {
		return 0
	}
	return len(s.Coeffs) - 1
}

// Len returns the number of coefficients in the series.
func (s *ChebSeries) Len() int { return len(s.Coeffs) }

// Coefficients returns a copy of the series coefficients.
func (s *ChebSeries) Coefficients() []float64 {
	out := make([]float64, len(s.Coeffs))
	copy(out, s.Coeffs)
	return out
}

// Domain returns the interval endpoints [A, B] of the series.
func (s *ChebSeries) Domain() (float64, float64) { return s.A, s.B }

// Clone returns a deep copy of the series.
func (s *ChebSeries) Clone() *ChebSeries {
	return NewChebSeries(s.Coeffs, s.A, s.B)
}

// ChebPoints returns the n+1 Chebyshev-Gauss-Lobatto points on [a, b], ordered
// from a to b. These are the extrema of T_n mapped to the interval and are the
// natural sampling grid for ChebFit.
func ChebPoints(n int, a, b float64) []float64 {
	out := make([]float64, n+1)
	for j := 0; j <= n; j++ {
		t := math.Cos(math.Pi * float64(j) / float64(n))
		x := 0.5*(a+b) + 0.5*(b-a)*t
		out[n-j] = x // reverse so the slice runs a..b
	}
	if n == 0 {
		out[0] = 0.5 * (a + b)
	}
	return out
}

// ChebGaussPoints returns the n Chebyshev-Gauss points (roots of T_n) on
// [a, b], ordered from a to b.
func ChebGaussPoints(n int, a, b float64) []float64 {
	out := make([]float64, n)
	for j := 0; j < n; j++ {
		t := math.Cos(math.Pi * (float64(j) + 0.5) / float64(n))
		out[n-1-j] = 0.5*(a+b) + 0.5*(b-a)*t
	}
	return out
}

// ChebCoeffsFromValues computes the Chebyshev interpolation coefficients from
// samples of a function taken at the n+1 Chebyshev-Gauss-Lobatto points in the
// order returned by ChebPoints (that is, values[i] is the sample at the i-th
// point running from a to b). It performs a discrete cosine transform.
func ChebCoeffsFromValues(values []float64) []float64 {
	N := len(values) - 1
	if N < 0 {
		return nil
	}
	if N == 0 {
		return []float64{values[0]}
	}
	// f indexed by the standard node ordering x_j = cos(pi j/N), j=0..N.
	// values runs a..b which corresponds to j=N..0, so reverse it.
	f := make([]float64, N+1)
	for j := 0; j <= N; j++ {
		f[j] = values[N-j]
	}
	coeffs := make([]float64, N+1)
	for k := 0; k <= N; k++ {
		var sum float64
		for j := 0; j <= N; j++ {
			gj := 1.0
			if j == 0 || j == N {
				gj = 2
			}
			sum += f[j] / gj * math.Cos(math.Pi*float64(k*j)/float64(N))
		}
		gk := 1.0
		if k == 0 || k == N {
			gk = 2
		}
		coeffs[k] = 2.0 / float64(N) / gk * sum
	}
	return coeffs
}

// ChebValuesFromCoeffs is the inverse of ChebCoeffsFromValues: it evaluates the
// Chebyshev series at the n+1 Chebyshev-Gauss-Lobatto nodes (ordered a..b via
// ChebPoints on [-1,1]) and returns those sample values.
func ChebValuesFromCoeffs(coeffs []float64) []float64 {
	N := len(coeffs) - 1
	if N < 0 {
		return nil
	}
	out := make([]float64, N+1)
	if N == 0 {
		out[0] = coeffs[0]
		return out
	}
	for i := 0; i <= N; i++ {
		// node ordering a..b corresponds to t running -1..1.
		t := math.Cos(math.Pi * float64(N-i) / float64(N))
		out[i] = ChebClenshaw(coeffs, t)
	}
	return out
}

// ChebFit returns the degree-n Chebyshev interpolant of f on [a, b], sampling f
// at the n+1 Chebyshev-Gauss-Lobatto points and transforming to coefficients.
func ChebFit(f func(float64) float64, n int, a, b float64) *ChebSeries {
	if n < 0 {
		n = 0
	}
	pts := ChebPoints(n, a, b)
	vals := make([]float64, len(pts))
	for i, x := range pts {
		vals[i] = f(x)
	}
	coeffs := ChebCoeffsFromValues(vals)
	return &ChebSeries{Coeffs: coeffs, A: a, B: b}
}

// ChebFitGauss returns the degree n-1 Chebyshev approximation of f built from
// samples at the n Chebyshev-Gauss points (roots of T_n). The coefficients are
// the discrete Chebyshev transform associated with Gauss quadrature.
func ChebFitGauss(f func(float64) float64, n int, a, b float64) *ChebSeries {
	if n < 1 {
		n = 1
	}
	pts := ChebGaussPoints(n, a, b)
	fv := make([]float64, n)
	for i, x := range pts {
		fv[i] = f(x)
	}
	coeffs := make([]float64, n)
	for k := 0; k < n; k++ {
		var sum float64
		for j := 0; j < n; j++ {
			// pts ordered a..b maps to t = cos(theta) with theta decreasing;
			// recompute theta for node j.
			theta := math.Pi * (float64(n-1-j) + 0.5) / float64(n)
			sum += fv[j] * math.Cos(float64(k)*theta)
		}
		factor := 2.0 / float64(n)
		if k == 0 {
			factor = 1.0 / float64(n)
		}
		coeffs[k] = factor * sum
	}
	return &ChebSeries{Coeffs: coeffs, A: a, B: b}
}

// Add returns the sum of two series. They must share the same domain.
func (s *ChebSeries) Add(o *ChebSeries) *ChebSeries {
	return &ChebSeries{Coeffs: PolyAdd(s.Coeffs, o.Coeffs), A: s.A, B: s.B}
}

// Sub returns the difference s-o of two series sharing the same domain.
func (s *ChebSeries) Sub(o *ChebSeries) *ChebSeries {
	return &ChebSeries{Coeffs: PolySub(s.Coeffs, o.Coeffs), A: s.A, B: s.B}
}

// Scale returns the series multiplied by the constant c.
func (s *ChebSeries) Scale(c float64) *ChebSeries {
	return &ChebSeries{Coeffs: PolyScale(s.Coeffs, c), A: s.A, B: s.B}
}

// Derivative returns the Chebyshev series of the derivative of s, correctly
// scaled for the domain [A, B].
func (s *ChebSeries) Derivative() *ChebSeries {
	c := make([]float64, len(s.Coeffs))
	copy(c, s.Coeffs)
	n := len(c)
	if n <= 1 {
		return &ChebSeries{Coeffs: []float64{0}, A: s.A, B: s.B}
	}
	der := make([]float64, n-1)
	for j := n - 1; j > 2; j-- {
		der[j-1] = 2 * float64(j) * c[j]
		c[j-2] += float64(j) * c[j] / float64(j-2)
	}
	if n > 2 {
		der[1] = 4 * c[2]
	}
	der[0] = c[1]
	// chain rule for the linear map t = (2x-(A+B))/(B-A): dt/dx = 2/(B-A).
	scale := 2.0 / (s.B - s.A)
	for i := range der {
		der[i] *= scale
	}
	return &ChebSeries{Coeffs: der, A: s.A, B: s.B}
}

// Antiderivative returns a Chebyshev series whose derivative is s and whose
// value at A is zero (an indefinite integral pinned to the left endpoint),
// correctly scaled for the domain [A, B].
func (s *ChebSeries) Antiderivative() *ChebSeries {
	c := s.Coeffs
	n := len(c)
	if n == 0 {
		return &ChebSeries{Coeffs: []float64{0}, A: s.A, B: s.B}
	}
	tmp := make([]float64, n+1)
	if n >= 1 {
		tmp[1] = c[0]
	}
	if n >= 2 {
		tmp[2] = c[1] / 4
	}
	for j := 2; j < n; j++ {
		tmp[j+1] = c[j] / float64(2*(j+1))
		tmp[j-1] -= c[j] / float64(2*(j-1))
	}
	// scale for the interval: dx = (B-A)/2 dt.
	scale := (s.B - s.A) / 2
	for i := range tmp {
		tmp[i] *= scale
	}
	res := &ChebSeries{Coeffs: tmp, A: s.A, B: s.B}
	res.Coeffs[0] -= res.Eval(s.A) // pin F(A) = 0
	return res
}

// DefiniteIntegral returns the integral of the series from x0 to x1.
func (s *ChebSeries) DefiniteIntegral(x0, x1 float64) float64 {
	F := s.Antiderivative()
	return F.Eval(x1) - F.Eval(x0)
}

// Integral returns the integral of the series over its whole domain [A, B].
func (s *ChebSeries) Integral() float64 {
	// int_{-1}^{1} T_k dt = 0 for odd k and 2/(1-k^2) for even k.
	var sum float64
	for k, ck := range s.Coeffs {
		if k%2 == 1 {
			continue
		}
		sum += ck * 2.0 / (1.0 - float64(k*k))
	}
	return sum * (s.B - s.A) / 2
}

// Truncate returns a copy of the series with trailing coefficients whose
// absolute value is at most tol removed. At least the constant term is kept.
func (s *ChebSeries) Truncate(tol float64) *ChebSeries {
	d := len(s.Coeffs)
	for d > 1 && math.Abs(s.Coeffs[d-1]) <= tol {
		d--
	}
	return NewChebSeries(s.Coeffs[:d], s.A, s.B)
}

// MaxCoeff returns the largest absolute coefficient of the series.
func (s *ChebSeries) MaxCoeff() float64 {
	var m float64
	for _, c := range s.Coeffs {
		if a := math.Abs(c); a > m {
			m = a
		}
	}
	return m
}

// TailNorm returns the sum of absolute values of the coefficients beyond index
// k, a convenient a-posteriori estimate of the truncation error.
func (s *ChebSeries) TailNorm(k int) float64 {
	var sum float64
	for i := k + 1; i < len(s.Coeffs); i++ {
		sum += math.Abs(s.Coeffs[i])
	}
	return sum
}

// Roots returns the real roots of the Chebyshev series lying in [A, B], sorted
// ascending. Roots are located by scanning a fine grid for sign changes and
// refining each bracket with a hybrid bisection/Newton iteration.
func (s *ChebSeries) Roots() []float64 {
	deg := s.Degree()
	if deg <= 0 {
		return nil
	}
	m := 8 * (deg + 1)
	if m < 64 {
		m = 64
	}
	der := s.Derivative()
	grid := Linspace(s.A, s.B, m+1)
	var roots []float64
	prevX := grid[0]
	prevV := s.Eval(prevX)
	if prevV == 0 {
		roots = append(roots, prevX)
	}
	for i := 1; i <= m; i++ {
		x := grid[i]
		v := s.Eval(x)
		if v == 0 {
			roots = append(roots, x)
			prevX, prevV = x, v
			continue
		}
		if prevV*v < 0 {
			r := refineRoot(s, der, prevX, x)
			roots = append(roots, r)
		}
		prevX, prevV = x, v
	}
	// Deduplicate roots that coincide (e.g. an exact grid hit next to a
	// bracketed crossing).
	sort.Float64s(roots)
	span := s.B - s.A
	uniq := roots[:0]
	for i, r := range roots {
		if i == 0 || r-uniq[len(uniq)-1] > 1e-10*span {
			uniq = append(uniq, r)
		}
	}
	return uniq
}

// refineRoot refines a bracketed root of s in [lo, hi] using Newton's method
// with bisection fallback.
func refineRoot(s, der *ChebSeries, lo, hi float64) float64 {
	flo := s.Eval(lo)
	x := 0.5 * (lo + hi)
	for it := 0; it < 60; it++ {
		fx := s.Eval(x)
		if fx == 0 {
			return x
		}
		if flo*fx < 0 {
			hi = x
		} else {
			lo = x
			flo = fx
		}
		d := der.Eval(x)
		nx := x
		if d != 0 {
			nx = x - fx/d
		}
		if nx <= lo || nx >= hi {
			nx = 0.5 * (lo + hi)
		}
		if math.Abs(nx-x) <= 1e-15*(1+math.Abs(x)) {
			return nx
		}
		x = nx
	}
	return x
}
