package quadrature

import "math"

// ChebyshevT evaluates the Chebyshev polynomial of the first kind of degree n
// at x using the recurrence T_{k+1} = 2x T_k - T_{k-1}. These polynomials
// underlie the Clenshaw-Curtis and Fejer rules.
func ChebyshevT(n int, x float64) float64 {
	if n == 0 {
		return 1
	}
	tkm1, tk := 1.0, x
	for k := 1; k < n; k++ {
		tkm1, tk = tk, 2*x*tk-tkm1
	}
	return tk
}

// ChebyshevU evaluates the Chebyshev polynomial of the second kind of degree n
// at x using the recurrence U_{k+1} = 2x U_k - U_{k-1}.
func ChebyshevU(n int, x float64) float64 {
	if n == 0 {
		return 1
	}
	ukm1, uk := 1.0, 2*x
	for k := 1; k < n; k++ {
		ukm1, uk = uk, 2*x*uk-ukm1
	}
	return uk
}

// HermiteH evaluates the physicists' Hermite polynomial H_n at x using the
// recurrence H_{k+1} = 2x H_k - 2k H_{k-1}.
func HermiteH(n int, x float64) float64 {
	if n == 0 {
		return 1
	}
	hkm1, hk := 1.0, 2*x
	for k := 1; k < n; k++ {
		hkm1, hk = hk, 2*x*hk-2*float64(k)*hkm1
	}
	return hk
}

// HermiteHe evaluates the probabilists' Hermite polynomial He_n at x using the
// recurrence He_{k+1} = x He_k - k He_{k-1}.
func HermiteHe(n int, x float64) float64 {
	if n == 0 {
		return 1
	}
	hkm1, hk := 1.0, x
	for k := 1; k < n; k++ {
		hkm1, hk = hk, x*hk-float64(k)*hkm1
	}
	return hk
}

// LaguerreL evaluates the Laguerre polynomial L_n at x using the recurrence
// (k+1) L_{k+1} = (2k+1-x) L_k - k L_{k-1}.
func LaguerreL(n int, x float64) float64 {
	return LaguerreLGen(n, 0, x)
}

// LaguerreLGen evaluates the generalized Laguerre polynomial L_n^{(a)} at x
// using the recurrence (k+1) L_{k+1} = (2k+1+a-x) L_k - (k+a) L_{k-1}.
func LaguerreLGen(n int, a, x float64) float64 {
	if n == 0 {
		return 1
	}
	lkm1, lk := 1.0, 1+a-x
	for k := 1; k < n; k++ {
		kk := float64(k)
		next := ((2*kk+1+a-x)*lk - (kk+a)*lkm1) / (kk + 1)
		lkm1, lk = lk, next
	}
	return lk
}

// JacobiP evaluates the Jacobi polynomial P_n^{(a,b)} at x using the standard
// three-term recurrence. Legendre is a = b = 0, Chebyshev of the first kind is
// a = b = -1/2, and so on.
func JacobiP(n int, a, b, x float64) float64 {
	if n == 0 {
		return 1
	}
	p0 := 1.0
	p1 := 0.5 * (a - b + (a+b+2)*x)
	if n == 1 {
		return p1
	}
	for k := 2; k <= n; k++ {
		kk := float64(k)
		c1 := 2 * kk * (kk + a + b) * (2*kk + a + b - 2)
		c2 := (2*kk + a + b - 1) * (a*a - b*b)
		c3 := (2*kk + a + b - 2) * (2*kk + a + b - 1) * (2*kk + a + b)
		c4 := 2 * (kk + a - 1) * (kk + b - 1) * (2*kk + a + b)
		p2 := ((c2+c3*x)*p1 - c4*p0) / c1
		p0, p1 = p1, p2
	}
	return p1
}

// GaussChebyshev3Rule returns the third-kind Gauss-Chebyshev rule of n points.
func GaussChebyshev3Rule(n int) Rule {
	nodes, weights := GaussChebyshev3(n)
	return Rule{Nodes: nodes, Weights: weights}
}

// GaussChebyshev4Rule returns the fourth-kind Gauss-Chebyshev rule of n points.
func GaussChebyshev4Rule(n int) Rule {
	nodes, weights := GaussChebyshev4(n)
	return Rule{Nodes: nodes, Weights: weights}
}

// IntegrateGaussChebyshev3 approximates the integral of
// f(x)*sqrt((1+x)/(1-x)) over [-1, 1] with the n-point third-kind rule.
func IntegrateGaussChebyshev3(f Func, n int) float64 {
	nodes, weights := GaussChebyshev3(n)
	var s float64
	for i, x := range nodes {
		s += weights[i] * f(x)
	}
	return s
}

// IntegrateGaussChebyshev4 approximates the integral of
// f(x)*sqrt((1-x)/(1+x)) over [-1, 1] with the n-point fourth-kind rule.
func IntegrateGaussChebyshev4(f Func, n int) float64 {
	nodes, weights := GaussChebyshev4(n)
	var s float64
	for i, x := range nodes {
		s += weights[i] * f(x)
	}
	return s
}

// IntegrateGaussHermiteProb approximates the integral of f(x)*e^{-x^2/2} over
// the whole real line with the n-point probabilists' Gauss-Hermite rule.
func IntegrateGaussHermiteProb(f Func, n int) float64 {
	nodes, weights := GaussHermiteProb(n)
	var s float64
	for i, x := range nodes {
		s += weights[i] * f(x)
	}
	return s
}

// ExpectationExponential approximates E[f(X)] for an exponential random
// variable X with rate 1 (density e^{-x} on [0, inf)) using the n-point
// Gauss-Laguerre rule.
func ExpectationExponential(f Func, n int) float64 {
	return IntegrateGaussLaguerre(f, n)
}

// Average returns the mean value of f over [a, b], namely the integral divided
// by the interval length, computed with adaptive Simpson quadrature.
func Average(f Func, a, b, tol float64) float64 {
	if a == b {
		return f(a)
	}
	return AdaptiveSimpson(f, a, b, tol) / (b - a)
}

// DoubleAverage returns the mean value of f over the rectangle
// [ax, bx] x [ay, by], computed with the n-point tensor Gauss-Legendre rule.
func DoubleAverage(f Func2, ax, bx, ay, by float64, n int) float64 {
	area := (bx - ax) * (by - ay)
	return DoubleGaussLegendre(f, ax, bx, ay, by, n) / area
}

// CompositeGaussLegendre approximates the integral of f over [a, b] by
// splitting the interval into panels equal subintervals and applying the
// n-point Gauss-Legendre rule on each. This combines the high polynomial
// accuracy of Gauss quadrature with the robustness of subdivision.
func CompositeGaussLegendre(f Func, a, b float64, n, panels int) float64 {
	if panels < 1 {
		panels = 1
	}
	nodes, weights := GaussLegendre(n)
	h := (b - a) / float64(panels)
	var total float64
	for p := 0; p < panels; p++ {
		pa := a + float64(p)*h
		half := 0.5 * h
		mid := pa + half
		var s float64
		for i, t := range nodes {
			s += weights[i] * f(mid+half*t)
		}
		total += half * s
	}
	return total
}

// AdaptiveLobatto approximates the integral of f over [a, b] to the requested
// absolute tolerance by recursively subdividing where a 3-point Gauss-Lobatto
// estimate disagrees with a 5-point one.
func AdaptiveLobatto(f Func, a, b, tol float64) float64 {
	if a == b {
		return 0
	}
	if tol <= 0 {
		tol = 1e-10
	}
	n3, w3 := GaussLobatto(3)
	n5, w5 := GaussLobatto(5)
	eval := func(nodes, weights []float64, a, b float64) float64 {
		half := 0.5 * (b - a)
		mid := 0.5 * (a + b)
		var s float64
		for i, t := range nodes {
			s += weights[i] * f(mid+half*t)
		}
		return half * s
	}
	var rec func(a, b, tol float64, depth int) float64
	rec = func(a, b, tol float64, depth int) float64 {
		coarse := eval(n3, w3, a, b)
		fine := eval(n5, w5, a, b)
		if depth <= 0 || math.Abs(fine-coarse) <= tol {
			return fine
		}
		m := 0.5 * (a + b)
		return rec(a, m, tol/2, depth-1) + rec(m, b, tol/2, depth-1)
	}
	return rec(a, b, tol, 50)
}

// MidpointSamples integrates samples ys, interpreted as the values of f at the
// midpoints of n equal panels of width h, using the midpoint rule.
func MidpointSamples(ys []float64, h float64) float64 {
	var s float64
	for _, y := range ys {
		s += y
	}
	return h * s
}

// BooleSamples integrates uniformly spaced samples ys with spacing h using the
// composite Boole rule. The number of samples must be one more than a multiple
// of four; otherwise it falls back to SimpsonSamples.
func BooleSamples(ys []float64, h float64) float64 {
	n := len(ys)
	if n < 5 || (n-1)%4 != 0 {
		return SimpsonSamples(ys, h)
	}
	var s float64
	for i := 0; i+4 < n; i += 4 {
		s += 7*ys[i] + 32*ys[i+1] + 12*ys[i+2] + 32*ys[i+3] + 7*ys[i+4]
	}
	return 2 * h / 45 * s
}

// IntegrateSamples integrates uniformly spaced samples ys with spacing h,
// selecting the highest-order composite rule the sample count permits: Boole
// when the number of panels is a multiple of four, otherwise Simpson when even,
// otherwise the trapezoidal rule.
func IntegrateSamples(ys []float64, h float64) float64 {
	n := len(ys)
	if n < 2 {
		return 0
	}
	panels := n - 1
	switch {
	case panels%4 == 0 && panels >= 4:
		return BooleSamples(ys, h)
	case panels%2 == 0:
		return SimpsonSamples(ys, h)
	default:
		return TrapezoidSamples(ys, h)
	}
}

// RombergValue approximates the integral of f over [a, b] with the k-level
// Romberg table, returning only the best estimate T[k-1][k-1].
func RombergValue(f Func, a, b float64, k int) float64 {
	T := RombergTable(f, a, b, k)
	return T[len(T)-1][len(T)-1]
}

// NewtonCotesClosedRule returns the closed (n+1)-point Newton-Cotes rule on
// the canonical interval [-1, 1] with equally spaced nodes. Use Rule.Scale to
// move it to another interval.
func NewtonCotesClosedRule(n int) Rule {
	w := NewtonCotesClosed(n)
	nodes := make([]float64, n+1)
	weights := make([]float64, n+1)
	for i := 0; i <= n; i++ {
		nodes[i] = -1 + 2*float64(i)/float64(n)
		weights[i] = 2 / float64(n) * w[i]
	}
	return Rule{Nodes: nodes, Weights: weights}
}

// NewtonCotesOpenRule returns the open n-point Newton-Cotes rule on the
// canonical interval [-1, 1] with equally spaced interior nodes.
func NewtonCotesOpenRule(n int) Rule {
	w := NewtonCotesOpen(n)
	nodes := make([]float64, n)
	weights := make([]float64, n)
	for i := 0; i < n; i++ {
		nodes[i] = -1 + 2*float64(i+1)/float64(n+1)
		weights[i] = 2 / float64(n+1) * w[i]
	}
	return Rule{Nodes: nodes, Weights: weights}
}

// TripleTrapezoid approximates the triple integral of f over the box using the
// composite trapezoidal rule with nx, ny and nz subintervals.
func TripleTrapezoid(f Func3, ax, bx, ay, by, az, bz float64, nx, ny, nz int) float64 {
	return CompositeTrapezoid(func(x float64) float64 {
		return CompositeTrapezoid(func(y float64) float64 {
			return CompositeTrapezoid(func(z float64) float64 { return f(x, y, z) }, az, bz, nz)
		}, ay, by, ny)
	}, ax, bx, nx)
}

// TripleMidpoint approximates the triple integral of f over the box using the
// composite midpoint rule with nx, ny and nz subintervals.
func TripleMidpoint(f Func3, ax, bx, ay, by, az, bz float64, nx, ny, nz int) float64 {
	return CompositeMidpoint(func(x float64) float64 {
		return CompositeMidpoint(func(y float64) float64 {
			return CompositeMidpoint(func(z float64) float64 { return f(x, y, z) }, az, bz, nz)
		}, ay, by, ny)
	}, ax, bx, nx)
}

// Transform returns a new rule obtained from r by the change of variables
// u = g(x): each node x_i becomes g(x_i) and each weight w_i is multiplied by
// the absolute derivative |dg(x_i)|. It lets a canonical rule be pushed
// through a smooth substitution.
func (r Rule) Transform(g, dg Func) Rule {
	n := len(r.Nodes)
	nodes := make([]float64, n)
	weights := make([]float64, n)
	for i, x := range r.Nodes {
		nodes[i] = g(x)
		weights[i] = r.Weights[i] * math.Abs(dg(x))
	}
	return Rule{Nodes: nodes, Weights: weights}
}
