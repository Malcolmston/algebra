package quadrature

import "math"

// RecurrenceLegendre returns the monic three-term recurrence coefficients for
// the Legendre polynomials (weight 1 on [-1, 1]) for an n-point rule. The
// returned alpha has length n and beta has length n with beta[0] equal to the
// zeroth moment mu0 = 2.
func RecurrenceLegendre(n int) (alpha, beta []float64) {
	alpha = make([]float64, n)
	beta = make([]float64, n)
	if n == 0 {
		return
	}
	beta[0] = 2
	for k := 1; k < n; k++ {
		kk := float64(k)
		beta[k] = kk * kk / (4*kk*kk - 1)
	}
	return
}

// RecurrenceChebyshev1 returns the monic recurrence coefficients for the
// Chebyshev polynomials of the first kind (weight 1/sqrt(1-x^2) on [-1, 1]).
// beta[0] = pi, beta[1] = 1/2 and beta[k] = 1/4 for k >= 2.
func RecurrenceChebyshev1(n int) (alpha, beta []float64) {
	alpha = make([]float64, n)
	beta = make([]float64, n)
	if n == 0 {
		return
	}
	beta[0] = math.Pi
	if n > 1 {
		beta[1] = 0.5
	}
	for k := 2; k < n; k++ {
		beta[k] = 0.25
	}
	return
}

// RecurrenceChebyshev2 returns the monic recurrence coefficients for the
// Chebyshev polynomials of the second kind (weight sqrt(1-x^2) on [-1, 1]).
// beta[0] = pi/2 and beta[k] = 1/4 for k >= 1.
func RecurrenceChebyshev2(n int) (alpha, beta []float64) {
	alpha = make([]float64, n)
	beta = make([]float64, n)
	if n == 0 {
		return
	}
	beta[0] = math.Pi / 2
	for k := 1; k < n; k++ {
		beta[k] = 0.25
	}
	return
}

// RecurrenceHermite returns the monic recurrence coefficients for the
// physicists' Hermite polynomials (weight e^{-x^2} on the whole real line).
// beta[0] = sqrt(pi) and beta[k] = k/2.
func RecurrenceHermite(n int) (alpha, beta []float64) {
	alpha = make([]float64, n)
	beta = make([]float64, n)
	if n == 0 {
		return
	}
	beta[0] = math.Sqrt(math.Pi)
	for k := 1; k < n; k++ {
		beta[k] = float64(k) / 2
	}
	return
}

// RecurrenceHermiteProb returns the monic recurrence coefficients for the
// probabilists' Hermite polynomials (weight e^{-x^2/2} on the whole real
// line). beta[0] = sqrt(2*pi) and beta[k] = k.
func RecurrenceHermiteProb(n int) (alpha, beta []float64) {
	alpha = make([]float64, n)
	beta = make([]float64, n)
	if n == 0 {
		return
	}
	beta[0] = math.Sqrt(2 * math.Pi)
	for k := 1; k < n; k++ {
		beta[k] = float64(k)
	}
	return
}

// RecurrenceLaguerre returns the monic recurrence coefficients for the
// generalized Laguerre polynomials with parameter a (weight x^a * e^{-x} on
// [0, inf)). alpha[k] = 2k+a+1, beta[0] = Gamma(a+1) and beta[k] = k*(k+a).
// Standard Laguerre corresponds to a = 0.
func RecurrenceLaguerre(n int, a float64) (alpha, beta []float64) {
	alpha = make([]float64, n)
	beta = make([]float64, n)
	if n == 0 {
		return
	}
	beta[0] = math.Gamma(a + 1)
	for k := 0; k < n; k++ {
		alpha[k] = 2*float64(k) + a + 1
	}
	for k := 1; k < n; k++ {
		kk := float64(k)
		beta[k] = kk * (kk + a)
	}
	return
}

// RecurrenceJacobi returns the monic recurrence coefficients for the Jacobi
// polynomials with parameters a and b (weight (1-x)^a * (1+x)^b on [-1, 1]).
// It unifies the Legendre (a=b=0), Chebyshev and Gegenbauer families. Both
// parameters must be greater than -1.
func RecurrenceJacobi(n int, a, b float64) (alpha, beta []float64) {
	alpha = make([]float64, n)
	beta = make([]float64, n)
	if n == 0 {
		return
	}
	ab := a + b
	beta[0] = math.Pow(2, ab+1) * math.Gamma(a+1) * math.Gamma(b+1) / math.Gamma(ab+2)
	alpha[0] = (b - a) / (ab + 2)
	if n > 1 {
		beta[1] = 4 * (a + 1) * (b + 1) / ((ab + 2) * (ab + 2) * (ab + 3))
	}
	for k := 1; k < n; k++ {
		kk := float64(k)
		den := (2*kk + ab) * (2*kk + ab + 2)
		alpha[k] = (b*b - a*a) / den
	}
	for k := 2; k < n; k++ {
		kk := float64(k)
		num := 4 * kk * (kk + a) * (kk + b) * (kk + ab)
		den := (2*kk + ab) * (2*kk + ab) * (2*kk + ab + 1) * (2*kk + ab - 1)
		beta[k] = num / den
	}
	return
}

// GaussLegendre returns the n nodes and weights of the Gauss-Legendre rule on
// [-1, 1], which integrates polynomials up to degree 2n-1 exactly. Nodes are
// sorted ascending and are symmetric about the origin.
func GaussLegendre(n int) (nodes, weights []float64) {
	alpha, beta := RecurrenceLegendre(n)
	return GolubWelsch(alpha, beta)
}

// GaussLegendreRule returns the Gauss-Legendre rule of n points on [-1, 1] as
// a Rule value.
func GaussLegendreRule(n int) Rule {
	nodes, weights := GaussLegendre(n)
	return Rule{Nodes: nodes, Weights: weights}
}

// GaussChebyshev1 returns the n nodes and weights of the Gauss-Chebyshev rule
// of the first kind, for the integral of f(x)/sqrt(1-x^2) over [-1, 1]. The
// nodes are cos((2i-1)pi/(2n)) and every weight equals pi/n.
func GaussChebyshev1(n int) (nodes, weights []float64) {
	nodes = make([]float64, n)
	weights = make([]float64, n)
	w := math.Pi / float64(n)
	for i := 0; i < n; i++ {
		// index from the top so that nodes come out ascending
		k := n - i
		nodes[i] = math.Cos((2*float64(k) - 1) * math.Pi / (2 * float64(n)))
		weights[i] = w
	}
	return
}

// GaussChebyshev1Rule returns the first-kind Gauss-Chebyshev rule of n points.
func GaussChebyshev1Rule(n int) Rule {
	nodes, weights := GaussChebyshev1(n)
	return Rule{Nodes: nodes, Weights: weights}
}

// GaussChebyshev2 returns the n nodes and weights of the Gauss-Chebyshev rule
// of the second kind, for the integral of f(x)*sqrt(1-x^2) over [-1, 1]. The
// nodes are cos(i*pi/(n+1)) and the weights are pi/(n+1)*sin^2(i*pi/(n+1)).
func GaussChebyshev2(n int) (nodes, weights []float64) {
	nodes = make([]float64, n)
	weights = make([]float64, n)
	for i := 0; i < n; i++ {
		k := n - i
		th := float64(k) * math.Pi / float64(n+1)
		s := math.Sin(th)
		nodes[i] = math.Cos(th)
		weights[i] = math.Pi / float64(n+1) * s * s
	}
	return
}

// GaussChebyshev2Rule returns the second-kind Gauss-Chebyshev rule of n points.
func GaussChebyshev2Rule(n int) Rule {
	nodes, weights := GaussChebyshev2(n)
	return Rule{Nodes: nodes, Weights: weights}
}

// GaussChebyshev3 returns the n nodes and weights of the Gauss-Chebyshev rule
// of the third kind, for the integral of f(x)*sqrt((1+x)/(1-x)) over [-1, 1].
// It is the Gauss-Jacobi rule with a = -1/2, b = +1/2.
func GaussChebyshev3(n int) (nodes, weights []float64) {
	return GaussJacobi(n, -0.5, 0.5)
}

// GaussChebyshev4 returns the n nodes and weights of the Gauss-Chebyshev rule
// of the fourth kind, for the integral of f(x)*sqrt((1-x)/(1+x)) over [-1, 1].
// It is the Gauss-Jacobi rule with a = +1/2, b = -1/2.
func GaussChebyshev4(n int) (nodes, weights []float64) {
	return GaussJacobi(n, 0.5, -0.5)
}

// GaussHermite returns the n nodes and weights of the Gauss-Hermite rule for
// the integral of f(x)*e^{-x^2} over the whole real line (physicists'
// convention). The sum of the weights is sqrt(pi).
func GaussHermite(n int) (nodes, weights []float64) {
	alpha, beta := RecurrenceHermite(n)
	return GolubWelsch(alpha, beta)
}

// GaussHermiteRule returns the physicists' Gauss-Hermite rule of n points.
func GaussHermiteRule(n int) Rule {
	nodes, weights := GaussHermite(n)
	return Rule{Nodes: nodes, Weights: weights}
}

// GaussHermiteProb returns the n nodes and weights of the Gauss-Hermite rule
// for the integral of f(x)*e^{-x^2/2} over the whole real line (probabilists'
// convention). The sum of the weights is sqrt(2*pi).
func GaussHermiteProb(n int) (nodes, weights []float64) {
	alpha, beta := RecurrenceHermiteProb(n)
	return GolubWelsch(alpha, beta)
}

// GaussHermiteProbRule returns the probabilists' Gauss-Hermite rule of n
// points.
func GaussHermiteProbRule(n int) Rule {
	nodes, weights := GaussHermiteProb(n)
	return Rule{Nodes: nodes, Weights: weights}
}

// GaussLaguerre returns the n nodes and weights of the Gauss-Laguerre rule for
// the integral of f(x)*e^{-x} over [0, inf). The sum of the weights is 1.
func GaussLaguerre(n int) (nodes, weights []float64) {
	return GaussLaguerreGen(n, 0)
}

// GaussLaguerreRule returns the Gauss-Laguerre rule of n points.
func GaussLaguerreRule(n int) Rule {
	nodes, weights := GaussLaguerre(n)
	return Rule{Nodes: nodes, Weights: weights}
}

// GaussLaguerreGen returns the n nodes and weights of the generalized
// Gauss-Laguerre rule for the integral of f(x)*x^a*e^{-x} over [0, inf). The
// sum of the weights is Gamma(a+1). The parameter a must be greater than -1.
func GaussLaguerreGen(n int, a float64) (nodes, weights []float64) {
	alpha, beta := RecurrenceLaguerre(n, a)
	return GolubWelsch(alpha, beta)
}

// GaussLaguerreGenRule returns the generalized Gauss-Laguerre rule of n points
// with parameter a.
func GaussLaguerreGenRule(n int, a float64) Rule {
	nodes, weights := GaussLaguerreGen(n, a)
	return Rule{Nodes: nodes, Weights: weights}
}

// GaussJacobi returns the n nodes and weights of the Gauss-Jacobi rule for the
// integral of f(x)*(1-x)^a*(1+x)^b over [-1, 1]. Both parameters must be
// greater than -1.
func GaussJacobi(n int, a, b float64) (nodes, weights []float64) {
	alpha, beta := RecurrenceJacobi(n, a, b)
	return GolubWelsch(alpha, beta)
}

// GaussJacobiRule returns the Gauss-Jacobi rule of n points with parameters a
// and b.
func GaussJacobiRule(n int, a, b float64) Rule {
	nodes, weights := GaussJacobi(n, a, b)
	return Rule{Nodes: nodes, Weights: weights}
}

// GaussGegenbauer returns the n nodes and weights of the Gauss-Gegenbauer
// (ultraspherical) rule for the integral of f(x)*(1-x^2)^(lambda-1/2) over
// [-1, 1]. It is the Gauss-Jacobi rule with a = b = lambda - 1/2, and lambda
// must be greater than -1/2.
func GaussGegenbauer(n int, lambda float64) (nodes, weights []float64) {
	p := lambda - 0.5
	return GaussJacobi(n, p, p)
}

// GaussGegenbauerRule returns the Gauss-Gegenbauer rule of n points.
func GaussGegenbauerRule(n int, lambda float64) Rule {
	nodes, weights := GaussGegenbauer(n, lambda)
	return Rule{Nodes: nodes, Weights: weights}
}

// legendrePD evaluates the Legendre polynomial of degree n and its derivative
// at x using the three-term recurrence.
func legendrePD(n int, x float64) (p, dp float64) {
	if n == 0 {
		return 1, 0
	}
	if n == 1 {
		return x, 1
	}
	p0, p1 := 1.0, x
	for k := 2; k <= n; k++ {
		kk := float64(k)
		p2 := ((2*kk-1)*x*p1 - (kk-1)*p0) / kk
		p0 = p1
		p1 = p2
	}
	if x*x == 1 {
		// Derivative at the endpoints: P_n'(1) = n(n+1)/2 and
		// P_n'(-1) = (-1)^{n-1} n(n+1)/2.
		dp = float64(n) * float64(n+1) / 2
		if x < 0 && n%2 == 0 {
			dp = -dp
		}
		return p1, dp
	}
	dp = float64(n) * (x*p1 - p0) / (x*x - 1)
	return p1, dp
}

// LegendreP evaluates the Legendre polynomial of degree n at x using the
// stable three-term recurrence. It is exposed because the Legendre polynomials
// underlie the Lobatto and Radau rules.
func LegendreP(n int, x float64) float64 {
	p, _ := legendrePD(n, x)
	return p
}

// LegendrePDeriv evaluates the derivative of the Legendre polynomial of degree
// n at x.
func LegendrePDeriv(n int, x float64) float64 {
	_, dp := legendrePD(n, x)
	return dp
}

// GaussLobatto returns the n nodes and weights of the Gauss-Lobatto-Legendre
// rule on [-1, 1]. The rule fixes both endpoints -1 and +1 as nodes and is
// exact for polynomials up to degree 2n-3. n must be at least 2.
func GaussLobatto(n int) (nodes, weights []float64) {
	alpha, beta := RecurrenceLegendre(n)
	return gaussLobattoGolub(alpha, beta, -1, 1)
}

// GaussLobattoRule returns the Gauss-Lobatto-Legendre rule of n points.
func GaussLobattoRule(n int) Rule {
	nodes, weights := GaussLobatto(n)
	return Rule{Nodes: nodes, Weights: weights}
}

// gaussLobattoGolub builds a Gauss-Lobatto rule that fixes the endpoints aEnd
// and bEnd by modifying the last diagonal and off-diagonal recurrence
// coefficients (Golub 1973) and then diagonalizing the resulting Jacobi
// matrix.
func gaussLobattoGolub(alpha, beta []float64, aEnd, bEnd float64) (nodes, weights []float64) {
	n := len(alpha)
	pa1, pa2 := monicEvalPair(alpha, beta, aEnd)
	pb1, pb2 := monicEvalPair(alpha, beta, bEnd)
	// Solve [[pa1, pa2],[pb1, pb2]] [al, be]^T = [aEnd*pa1, bEnd*pb1]^T.
	det := pa1*pb2 - pa2*pb1
	rhs0 := aEnd * pa1
	rhs1 := bEnd * pb1
	al := (rhs0*pb2 - pa2*rhs1) / det
	be := (pa1*rhs1 - rhs0*pb1) / det
	mod := make([]float64, n)
	copy(mod, alpha)
	mod[n-1] = al
	modBeta := make([]float64, n)
	copy(modBeta, beta)
	modBeta[n-1] = be
	return GolubWelsch(mod, modBeta)
}

// GaussRadau returns the n nodes and weights of the Gauss-Radau-Legendre rule
// on [-1, 1] that fixes the left endpoint -1 as a node. It is exact for
// polynomials up to degree 2n-2. n must be at least 1.
func GaussRadau(n int) (nodes, weights []float64) {
	alpha, beta := RecurrenceLegendre(n)
	return gaussRadauGolub(alpha, beta, -1)
}

// GaussRadauRight returns the n nodes and weights of the Gauss-Radau-Legendre
// rule on [-1, 1] that fixes the right endpoint +1 as a node.
func GaussRadauRight(n int) (nodes, weights []float64) {
	alpha, beta := RecurrenceLegendre(n)
	return gaussRadauGolub(alpha, beta, 1)
}

// GaussRadauRule returns the left Gauss-Radau-Legendre rule of n points.
func GaussRadauRule(n int) Rule {
	nodes, weights := GaussRadau(n)
	return Rule{Nodes: nodes, Weights: weights}
}

// gaussRadauGolub builds a Gauss-Radau rule that fixes the node end by
// modifying the last diagonal recurrence coefficient (Golub 1973).
func gaussRadauGolub(alpha, beta []float64, end float64) (nodes, weights []float64) {
	n := len(alpha)
	if n == 1 {
		// Single node forced to the endpoint with the full mass.
		return []float64{end}, []float64{beta[0]}
	}
	p1, p2 := monicEvalPair(alpha, beta, end)
	mod := make([]float64, n)
	copy(mod, alpha)
	mod[n-1] = end - beta[n-1]*p2/p1
	return GolubWelsch(mod, beta)
}
