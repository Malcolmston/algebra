package infogeom

import "math"

// FisherRaoBernoulli returns the Fisher-Rao geodesic distance between two
// Bernoulli distributions with success probabilities p and q. Under the
// square-root embedding the Bernoulli manifold is an arc of a circle of radius
// two, giving 2 arccos( sqrt(pq) + sqrt((1-p)(1-q)) ). It returns ErrDomain
// when p or q lies outside [0,1].
func FisherRaoBernoulli(p, q float64) (float64, error) {
	if p < 0 || p > 1 || q < 0 || q > 1 {
		return 0, ErrDomain
	}
	bc := math.Sqrt(p*q) + math.Sqrt((1-p)*(1-q))
	if bc > 1 {
		bc = 1
	}
	return 2 * math.Acos(bc), nil
}

// FisherRaoCategorical returns the Fisher-Rao geodesic distance between two
// categorical distributions p and q on the same support. The square-root map
// sends the probability simplex isometrically to the positive orthant of a
// sphere of radius two, so the distance is 2 arccos( sum sqrt(p_i q_i) ), twice
// the angle whose cosine is the Bhattacharyya coefficient. It returns
// ErrNotProb / ErrDim on malformed input.
func FisherRaoCategorical(p, q []float64) (float64, error) {
	bc, err := BhattacharyyaCoefficient(p, q)
	if err != nil {
		return 0, err
	}
	if bc > 1 {
		bc = 1
	}
	return 2 * math.Acos(bc), nil
}

// FisherRaoGaussian returns the Fisher-Rao geodesic distance between two
// univariate Gaussians. The location-scale family carries the hyperbolic
// metric ds^2 = (d mu^2 + 2 d sigma^2)/sigma^2, whose geodesic distance is
//
//	sqrt(2) arccosh( 1 + ((mu1-mu2)^2 + 2(sigma1-sigma2)^2) / (4 sigma1 sigma2) ).
//
// It returns ErrDomain when either standard deviation is not positive.
func FisherRaoGaussian(g, h Gaussian) (float64, error) {
	if !g.Valid() || !h.Valid() {
		return 0, ErrDomain
	}
	dm := g.Mu - h.Mu
	ds := g.Sigma - h.Sigma
	arg := 1 + (dm*dm+2*ds*ds)/(4*g.Sigma*h.Sigma)
	if arg < 1 {
		arg = 1
	}
	return math.Sqrt2 * math.Acosh(arg), nil
}

// FisherRaoGaussianFixedVariance returns the Fisher-Rao distance between two
// Gaussians that share the standard deviation sigma but differ in mean, which
// reduces to the Euclidean distance |mu1-mu2|/sigma. It returns ErrDomain when
// sigma is not positive.
func FisherRaoGaussianFixedVariance(mu1, mu2, sigma float64) (float64, error) {
	if sigma <= 0 {
		return 0, ErrDomain
	}
	return math.Abs(mu1-mu2) / sigma, nil
}

// FisherRaoGaussianFixedMean returns the Fisher-Rao distance between two
// Gaussians that share the mean but differ in standard deviation, which along
// the scale axis reduces to sqrt(2) |ln(sigma1/sigma2)|. It returns ErrDomain
// when either standard deviation is not positive.
func FisherRaoGaussianFixedMean(sigma1, sigma2 float64) (float64, error) {
	if sigma1 <= 0 || sigma2 <= 0 {
		return 0, ErrDomain
	}
	return math.Sqrt2 * math.Abs(math.Log(sigma1/sigma2)), nil
}

// FisherRaoPoisson returns the Fisher-Rao geodesic distance between two Poisson
// distributions with rates lambda1 and lambda2. With scalar metric 1/lambda the
// distance integrates to 2 |sqrt(lambda1) - sqrt(lambda2)|. It returns
// ErrDomain when either rate is negative.
func FisherRaoPoisson(lambda1, lambda2 float64) (float64, error) {
	if lambda1 < 0 || lambda2 < 0 {
		return 0, ErrDomain
	}
	return 2 * math.Abs(math.Sqrt(lambda1)-math.Sqrt(lambda2)), nil
}

// FisherRaoExponential returns the Fisher-Rao geodesic distance between two
// exponential distributions with rates rate1 and rate2. With scalar metric
// 1/rate^2 the distance integrates to |ln(rate1/rate2)|. It returns ErrDomain
// when either rate is not positive.
func FisherRaoExponential(rate1, rate2 float64) (float64, error) {
	if rate1 <= 0 || rate2 <= 0 {
		return 0, ErrDomain
	}
	return math.Abs(math.Log(rate1 / rate2)), nil
}

// FisherRaoScalar returns the Fisher-Rao geodesic distance between two points a
// and b of a one-dimensional statistical model whose scalar Fisher information
// is given by fisher. The distance is the arc length integral of sqrt(fisher)
// from a to b, evaluated by composite Simpson's rule with n subintervals
// (rounded up to an even number). It returns ErrDomain when fisher is negative
// anywhere on the sampled grid, and ErrDim when n < 1.
func FisherRaoScalar(fisher func(x float64) float64, a, b float64, n int) (float64, error) {
	if n < 1 {
		return 0, ErrDim
	}
	if n%2 == 1 {
		n++
	}
	if a == b {
		return 0, nil
	}
	h := (b - a) / float64(n)
	g := func(x float64) (float64, error) {
		v := fisher(x)
		if v < 0 {
			return 0, ErrDomain
		}
		return math.Sqrt(v), nil
	}
	fa, err := g(a)
	if err != nil {
		return 0, err
	}
	fb, err := g(b)
	if err != nil {
		return 0, err
	}
	sum := fa + fb
	for i := 1; i < n; i++ {
		x := a + float64(i)*h
		fi, err := g(x)
		if err != nil {
			return 0, err
		}
		if i%2 == 1 {
			sum += 4 * fi
		} else {
			sum += 2 * fi
		}
	}
	return math.Abs(h) / 3 * sum, nil
}

// FisherRaoPathLength returns the Riemannian length of a discretised path
// through parameter space under a position-dependent metric tensor. The metric
// function returns the Fisher information matrix at a point; the length is the
// sum over segments of sqrt( d^T G(midpoint) d ) with d the segment vector. It
// returns ErrDim when the path has fewer than two points or ragged entries.
func FisherRaoPathLength(metric func(x []float64) [][]float64, path [][]float64) (float64, error) {
	if len(path) < 2 {
		return 0, ErrDim
	}
	n := len(path[0])
	var length float64
	for k := 0; k+1 < len(path); k++ {
		if len(path[k]) != n || len(path[k+1]) != n {
			return 0, ErrDim
		}
		d, err := SubVectors(path[k+1], path[k])
		if err != nil {
			return 0, err
		}
		mid := make([]float64, n)
		for i := range mid {
			mid[i] = 0.5 * (path[k][i] + path[k+1][i])
		}
		g := metric(mid)
		q, err := QuadraticForm(g, d)
		if err != nil {
			return 0, err
		}
		if q < 0 {
			q = 0
		}
		length += math.Sqrt(q)
	}
	return length, nil
}
