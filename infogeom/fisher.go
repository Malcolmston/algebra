package infogeom

import "math"

// FisherInformationBernoulli returns the scalar Fisher information
// 1/(p(1-p)) of the Bernoulli family at success probability p. It returns
// ErrDomain when p is not strictly in (0,1).
func FisherInformationBernoulli(p float64) (float64, error) {
	return Bernoulli{P: p}.FisherInformation()
}

// FisherInformationPoisson returns the scalar Fisher information 1/lambda of
// the Poisson family. It returns ErrDomain when lambda is not positive.
func FisherInformationPoisson(lambda float64) (float64, error) {
	return Poisson{Lambda: lambda}.FisherInformation()
}

// FisherInformationExponential returns the scalar Fisher information 1/rate^2
// of the exponential family. It returns ErrDomain when rate is not positive.
func FisherInformationExponential(rate float64) (float64, error) {
	return Exponential{Rate: rate}.FisherInformation()
}

// FisherInformationGaussianMean returns the Fisher information 1/sigma^2 of the
// Gaussian mean parameter with known standard deviation sigma. It returns
// ErrDomain when sigma is not positive.
func FisherInformationGaussianMean(sigma float64) (float64, error) {
	if sigma <= 0 {
		return 0, ErrDomain
	}
	return 1 / (sigma * sigma), nil
}

// FisherInformationSimplex returns the Fisher information matrix of the
// categorical family expressed in the full probability coordinates, the
// diagonal matrix diag(1/p_i). It is the metric induced on the probability
// simplex by the Fisher-Rao structure and is singular in the ambient space
// (the physical tangent space is the subspace sum dp_i = 0). It returns
// ErrNotProb when p is not a probability vector.
func FisherInformationSimplex(p []float64) ([][]float64, error) {
	if !IsProbabilityVector(p, probTol) {
		return nil, ErrNotProb
	}
	d := make([]float64, len(p))
	for i, pi := range p {
		if pi <= 0 {
			return nil, ErrDomain
		}
		d[i] = 1 / pi
	}
	return Diagonal(d), nil
}

// FisherInformationCategorical returns the (k-1)-by-(k-1) Fisher information
// matrix of the categorical family in the reduced mean coordinates
// (p_1,...,p_{k-1}) with the last outcome eliminated, namely
// g_ij = delta_ij/p_i + 1/p_k. It returns ErrNotProb when p is not a
// probability vector and ErrDomain when any used probability vanishes.
func FisherInformationCategorical(p []float64) ([][]float64, error) {
	if !IsProbabilityVector(p, probTol) {
		return nil, ErrNotProb
	}
	k := len(p)
	if k < 2 {
		return nil, ErrDim
	}
	pk := p[k-1]
	if pk <= 0 {
		return nil, ErrDomain
	}
	g := make([][]float64, k-1)
	for i := 0; i < k-1; i++ {
		if p[i] <= 0 {
			return nil, ErrDomain
		}
		g[i] = make([]float64, k-1)
		for j := 0; j < k-1; j++ {
			if i == j {
				g[i][j] = 1/p[i] + 1/pk
			} else {
				g[i][j] = 1 / pk
			}
		}
	}
	return g, nil
}

// FisherInformationMultivariateGaussianMean returns the Fisher information
// matrix of the mean parameter of a multivariate Gaussian with known
// covariance cov, which is the inverse covariance (precision) matrix. It
// returns ErrSingular when cov is not invertible.
func FisherInformationMultivariateGaussianMean(cov [][]float64) ([][]float64, error) {
	return Inverse(cov)
}

// FisherInformationFromLogPartition returns the Fisher information matrix of an
// exponential family in natural coordinates as the Hessian of the log-partition
// function A evaluated at theta, computed by finite differences with step h.
func FisherInformationFromLogPartition(a func(theta []float64) float64, theta []float64, h float64) [][]float64 {
	return NumericalHessian(a, theta, h)
}

// ObservedInformation returns the observed information matrix, the negative
// Hessian of a log-likelihood function logLik evaluated at the parameter
// theta, computed by finite differences with step h.
func ObservedInformation(logLik func(theta []float64) float64, theta []float64, h float64) [][]float64 {
	return ScaleMatrix(-1, NumericalHessian(logLik, theta, h))
}

// FisherInformationFromScores returns an empirical estimate of the Fisher
// information matrix as the average outer product of per-sample score vectors,
// (1/m) sum_j s_j s_j^T. It returns ErrDim when the score set is empty or
// ragged.
func FisherInformationFromScores(scores [][]float64) ([][]float64, error) {
	m := len(scores)
	if m == 0 {
		return nil, ErrDim
	}
	n := len(scores[0])
	if n == 0 {
		return nil, ErrDim
	}
	g := make([][]float64, n)
	for i := range g {
		g[i] = make([]float64, n)
	}
	for _, s := range scores {
		if len(s) != n {
			return nil, ErrDim
		}
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				g[i][j] += s[i] * s[j]
			}
		}
	}
	inv := 1 / float64(m)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			g[i][j] *= inv
		}
	}
	return g, nil
}

// MetricInnerProduct returns the inner product u^T G v of two tangent vectors
// under the metric tensor G. It returns ErrDim on a shape mismatch.
func MetricInnerProduct(g [][]float64, u, v []float64) (float64, error) {
	return BilinearForm(g, u, v)
}

// MetricNorm returns the length sqrt(v^T G v) of a tangent vector v under the
// metric tensor G. It returns ErrDim on a shape mismatch and ErrDomain when
// the quadratic form is negative (G not positive semidefinite at v).
func MetricNorm(g [][]float64, v []float64) (float64, error) {
	q, err := QuadraticForm(g, v)
	if err != nil {
		return 0, err
	}
	if q < 0 {
		return 0, ErrDomain
	}
	return math.Sqrt(q), nil
}

// MetricAngle returns the angle in radians between tangent vectors u and v
// measured with the metric tensor G. It returns ErrDomain when either vector
// has zero length.
func MetricAngle(g [][]float64, u, v []float64) (float64, error) {
	iuv, err := MetricInnerProduct(g, u, v)
	if err != nil {
		return 0, err
	}
	nu, err := MetricNorm(g, u)
	if err != nil {
		return 0, err
	}
	nv, err := MetricNorm(g, v)
	if err != nil {
		return 0, err
	}
	if nu == 0 || nv == 0 {
		return 0, ErrDomain
	}
	c := iuv / (nu * nv)
	if c > 1 {
		c = 1
	} else if c < -1 {
		c = -1
	}
	return math.Acos(c), nil
}
