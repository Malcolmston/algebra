package infogeom

import "math"

// AlphaEmbedding returns the alpha-representation (alpha-embedding) of the
// probability vector p, applied elementwise. For alpha != 1 the map is
// l_alpha(u) = (2/(1-alpha)) u^{(1-alpha)/2}; the limit alpha = 1 gives the
// logarithm ln u. The e-representation (alpha = 1) and m-representation
// (alpha = -1) are the two flat coordinate systems of the dually flat
// structure. It returns ErrNotProb when p is not a probability vector.
func AlphaEmbedding(p []float64, alpha float64) ([]float64, error) {
	if !IsProbabilityVector(p, probTol) {
		return nil, ErrNotProb
	}
	out := make([]float64, len(p))
	if math.Abs(alpha-1) < 1e-12 {
		for i, pi := range p {
			out[i] = math.Log(pi)
		}
		return out, nil
	}
	c := 2 / (1 - alpha)
	e := (1 - alpha) / 2
	for i, pi := range p {
		out[i] = c * math.Pow(pi, e)
	}
	return out, nil
}

// EscortDistribution returns the alpha-escort distribution of p, the normalised
// power p_i^alpha / sum_j p_j^alpha. The escort at alpha = 1 is p itself; large
// alpha concentrates on the modes and small alpha flattens the distribution. It
// returns ErrNotProb when p is not a probability vector and ErrDomain when
// alpha is negative and p has a zero entry.
func EscortDistribution(p []float64, alpha float64) ([]float64, error) {
	if !IsProbabilityVector(p, probTol) {
		return nil, ErrNotProb
	}
	w := make([]float64, len(p))
	for i, pi := range p {
		if pi <= 0 {
			if alpha <= 0 {
				return nil, ErrDomain
			}
			w[i] = 0
			continue
		}
		w[i] = math.Pow(pi, alpha)
	}
	return Normalize(w)
}

// MixtureGeodesic returns the point at parameter t in [0,1] on the m-geodesic
// (mixture geodesic) joining the distributions p and q, the convex combination
// (1-t)p + t q. Mixture geodesics are the straight lines of the mixture
// (m-) affine connection. It returns ErrDim / ErrNotProb on malformed input.
func MixtureGeodesic(p, q []float64, t float64) ([]float64, error) {
	if err := checkProbPair(p, q); err != nil {
		return nil, err
	}
	out := make([]float64, len(p))
	for i := range p {
		out[i] = (1-t)*p[i] + t*q[i]
	}
	return out, nil
}

// ExponentialGeodesic returns the point at parameter t in [0,1] on the
// e-geodesic (exponential geodesic) joining the distributions p and q, the
// normalised geometric interpolation proportional to p_i^{1-t} q_i^t. Exponential
// geodesics are the straight lines of the exponential (e-) affine connection.
// It returns ErrDim / ErrNotProb on malformed input.
func ExponentialGeodesic(p, q []float64, t float64) ([]float64, error) {
	if err := checkProbPair(p, q); err != nil {
		return nil, err
	}
	logw := make([]float64, len(p))
	for i := range p {
		switch {
		case p[i] <= 0 || q[i] <= 0:
			logw[i] = math.Inf(-1)
		default:
			logw[i] = (1-t)*math.Log(p[i]) + t*math.Log(q[i])
		}
	}
	return Softmax(logw), nil
}

// AlphaMixture returns the point at parameter t in [0,1] on the alpha-geodesic
// joining the distributions p and q. For alpha != 1 it interpolates the
// alpha-representations and maps back, giving an unnormalised weight
// ((1-t) p_i^{(1-alpha)/2} + t q_i^{(1-alpha)/2})^{2/(1-alpha)}, which is then
// renormalised; alpha = 1 recovers the exponential geodesic. The cases
// alpha = -1 and alpha = 1 reproduce the mixture and exponential geodesics
// respectively. It returns ErrDim / ErrNotProb on malformed input.
func AlphaMixture(p, q []float64, t, alpha float64) ([]float64, error) {
	if err := checkProbPair(p, q); err != nil {
		return nil, err
	}
	if math.Abs(alpha-1) < 1e-12 {
		return ExponentialGeodesic(p, q, t)
	}
	e := (1 - alpha) / 2
	inv := 1 / e
	w := make([]float64, len(p))
	for i := range p {
		m := (1-t)*math.Pow(p[i], e) + t*math.Pow(q[i], e)
		if m <= 0 {
			w[i] = 0
		} else {
			w[i] = math.Pow(m, inv)
		}
	}
	return Normalize(w)
}

// CanonicalDivergence returns the canonical divergence of a dually flat space,
//
//	D(p, q) = F(theta_p) + Fstar(eta_q) - <theta_p, eta_q>,
//
// where F is the primal potential (log-partition), Fstar its Legendre conjugate,
// theta_p the natural coordinates of p and eta_q the expectation coordinates of
// q. It coincides with the Bregman divergence of F and, for exponential
// families, with the Kullback-Leibler divergence D(p || q). It returns ErrDim
// on a length mismatch.
func CanonicalDivergence(fTheta, fStarEta float64, thetaP, etaQ []float64) (float64, error) {
	inner, err := Dot(thetaP, etaQ)
	if err != nil {
		return 0, err
	}
	return fTheta + fStarEta - inner, nil
}

// NaturalFromExpectation inverts the expectation map of an exponential family,
// returning the natural parameters theta with grad A(theta) = eta. It runs
// Newton's method from the initial guess theta0 using the family's Fisher
// information (the Hessian of A) and stops when the update is smaller than tol
// or after maxIter iterations. It returns ErrDim on a length mismatch,
// ErrSingular when the Fisher matrix is singular, and ErrDomain when the
// iteration fails to converge.
func NaturalFromExpectation(f ExponentialFamily, eta, theta0 []float64, tol float64, maxIter int) ([]float64, error) {
	if len(eta) != len(theta0) || len(eta) == 0 {
		return nil, ErrDim
	}
	theta := CloneVector(theta0)
	for it := 0; it < maxIter; it++ {
		cur := f.MeanParameters(theta)
		if len(cur) != len(eta) {
			return nil, ErrDim
		}
		resid, err := SubVectors(cur, eta)
		if err != nil {
			return nil, err
		}
		g := f.FisherInformationNatural(theta)
		delta, err := Solve(g, resid)
		if err != nil {
			return nil, err
		}
		var maxStep float64
		for i := range theta {
			theta[i] -= delta[i]
			if a := math.Abs(delta[i]); a > maxStep {
				maxStep = a
			}
		}
		if maxStep < tol {
			return theta, nil
		}
	}
	return nil, ErrDomain
}

// CumulantTensor3 returns the third-order cumulant tensor of an exponential
// family in natural coordinates, T_{ijk} = d^3 A / d theta_i d theta_j d theta_k,
// computed by differentiating the supplied Fisher (Hessian) function once with a
// central finite difference of step h. This tensor is the skewness of the
// sufficient statistics and generates the family's alpha-connections. The
// result is indexed as T[i][j][k].
func CumulantTensor3(fisher func(theta []float64) [][]float64, theta []float64, h float64) [][][]float64 {
	n := len(theta)
	t := make([][][]float64, n)
	for i := range t {
		t[i] = make([][]float64, n)
		for j := range t[i] {
			t[i][j] = make([]float64, n)
		}
	}
	xp := CloneVector(theta)
	for k := 0; k < n; k++ {
		orig := xp[k]
		xp[k] = orig + h
		hp := fisher(xp)
		xp[k] = orig - h
		hm := fisher(xp)
		xp[k] = orig
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				t[i][j][k] = (hp[i][j] - hm[i][j]) / (2 * h)
			}
		}
	}
	return t
}

// AlphaConnectionCoefficients returns the coefficients of the Amari
// alpha-connection of an exponential family in natural coordinates, the lowered
// tensor Gamma^{(alpha)}_{ij,k} = ((1-alpha)/2) T_{ijk}, where T is the third
// cumulant tensor. The e-connection (alpha = 1) has vanishing coefficients,
// reflecting the flatness of natural coordinates, while the m-connection
// (alpha = -1) gives the full skewness tensor. The result is indexed [i][j][k].
func AlphaConnectionCoefficients(fisher func(theta []float64) [][]float64, theta []float64, alpha, h float64) [][][]float64 {
	t := CumulantTensor3(fisher, theta, h)
	s := (1 - alpha) / 2
	for i := range t {
		for j := range t[i] {
			for k := range t[i][j] {
				t[i][j][k] *= s
			}
		}
	}
	return t
}

// DualAlpha returns the dual connection order -alpha. The alpha- and
// (-alpha)-connections are mutually dual with respect to the Fisher metric, the
// e- and m-connections (alpha = +1 and -1) being the canonical dual pair.
func DualAlpha(alpha float64) float64 { return -alpha }
