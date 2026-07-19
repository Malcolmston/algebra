package infogeom

import "math"

// Bernoulli is the Bernoulli distribution with success probability P in (0,1).
// It is the exponential family with sufficient statistic x, natural parameter
// the log-odds, and log-partition ln(1+e^theta).
type Bernoulli struct {
	// P is the probability of the outcome 1.
	P float64
}

// Valid reports whether the success probability lies strictly in (0,1).
func (b Bernoulli) Valid() bool { return b.P > 0 && b.P < 1 }

// Mean returns the expectation P of the Bernoulli distribution.
func (b Bernoulli) Mean() float64 { return b.P }

// Variance returns the variance P(1-P) of the Bernoulli distribution.
func (b Bernoulli) Variance() float64 { return b.P * (1 - b.P) }

// PMF returns the probability mass at x in {0,1}; it is zero for any other x.
func (b Bernoulli) PMF(x int) float64 {
	switch x {
	case 1:
		return b.P
	case 0:
		return 1 - b.P
	default:
		return 0
	}
}

// NaturalParameter returns the natural parameter (log-odds) ln(P/(1-P)). It
// returns ErrDomain when P is not in (0,1).
func (b Bernoulli) NaturalParameter() (float64, error) {
	if !b.Valid() {
		return 0, ErrDomain
	}
	return math.Log(b.P / (1 - b.P)), nil
}

// LogPartition returns the log-partition A(theta) = ln(1+e^theta) at the
// distribution's natural parameter.
func (b Bernoulli) LogPartition() (float64, error) {
	theta, err := b.NaturalParameter()
	if err != nil {
		return 0, err
	}
	return logOnePlusExp(theta), nil
}

// FisherInformation returns the scalar Fisher information 1/(P(1-P)) of the
// Bernoulli family in the mean parameterisation. It returns ErrDomain when P
// is not in (0,1).
func (b Bernoulli) FisherInformation() (float64, error) {
	if !b.Valid() {
		return 0, ErrDomain
	}
	return 1 / (b.P * (1 - b.P)), nil
}

// BernoulliFromNatural returns the Bernoulli distribution with natural
// parameter (log-odds) theta, i.e. P = sigmoid(theta).
func BernoulliFromNatural(theta float64) Bernoulli {
	return Bernoulli{P: Sigmoid(theta)}
}

// Categorical is a categorical distribution over k outcomes, represented by
// its probability vector P whose entries are non-negative and sum to one.
type Categorical struct {
	// P holds the outcome probabilities.
	P []float64
}

// Valid reports whether P is a probability vector.
func (c Categorical) Valid() bool { return IsProbabilityVector(c.P, probTol) }

// Mean returns the probability vector itself, the expectation of the one-hot
// sufficient statistic.
func (c Categorical) Mean() []float64 { return CloneVector(c.P) }

// Entropy returns the Shannon entropy of the categorical distribution in nats.
func (c Categorical) Entropy() (float64, error) { return Entropy(c.P) }

// NaturalParameters returns the k-1 natural parameters theta_i =
// ln(P_i/P_{k-1}) of the categorical family taken with the last outcome as the
// reference category. It returns ErrNotProb when P is invalid and ErrDomain
// when the reference probability is zero.
func (c Categorical) NaturalParameters() ([]float64, error) {
	if !c.Valid() {
		return nil, ErrNotProb
	}
	k := len(c.P)
	ref := c.P[k-1]
	if ref <= 0 {
		return nil, ErrDomain
	}
	theta := make([]float64, k-1)
	for i := 0; i < k-1; i++ {
		if c.P[i] <= 0 {
			return nil, ErrDomain
		}
		theta[i] = math.Log(c.P[i] / ref)
	}
	return theta, nil
}

// CategoricalFromNatural returns the categorical distribution whose k-1 natural
// parameters are theta (with the last outcome as reference), recovering P via
// the softmax over (theta, 0).
func CategoricalFromNatural(theta []float64) Categorical {
	ext := make([]float64, len(theta)+1)
	copy(ext, theta)
	// last entry stays 0 (reference category).
	return Categorical{P: Softmax(ext)}
}

// Poisson is the Poisson distribution with rate Lambda > 0.
type Poisson struct {
	// Lambda is the mean (and variance) of the distribution.
	Lambda float64
}

// Valid reports whether the rate is positive.
func (p Poisson) Valid() bool { return p.Lambda > 0 }

// Mean returns the rate Lambda.
func (p Poisson) Mean() float64 { return p.Lambda }

// Variance returns the variance Lambda of the Poisson distribution.
func (p Poisson) Variance() float64 { return p.Lambda }

// PMF returns the probability mass e^{-Lambda} Lambda^k / k! at the count k.
func (p Poisson) PMF(k int) float64 {
	if k < 0 || !p.Valid() {
		return 0
	}
	return math.Exp(float64(k)*math.Log(p.Lambda) - p.Lambda - logFactorial(k))
}

// NaturalParameter returns the natural parameter ln(Lambda). It returns
// ErrDomain when Lambda is not positive.
func (p Poisson) NaturalParameter() (float64, error) {
	if !p.Valid() {
		return 0, ErrDomain
	}
	return math.Log(p.Lambda), nil
}

// FisherInformation returns the scalar Fisher information 1/Lambda of the
// Poisson family in the rate parameterisation. It returns ErrDomain when
// Lambda is not positive.
func (p Poisson) FisherInformation() (float64, error) {
	if !p.Valid() {
		return 0, ErrDomain
	}
	return 1 / p.Lambda, nil
}

// Exponential is the exponential distribution with rate parameter Rate > 0 and
// density Rate * e^{-Rate x} on x >= 0.
type Exponential struct {
	// Rate is the inverse of the mean.
	Rate float64
}

// Valid reports whether the rate is positive.
func (e Exponential) Valid() bool { return e.Rate > 0 }

// Mean returns the mean 1/Rate.
func (e Exponential) Mean() float64 { return 1 / e.Rate }

// Variance returns the variance 1/Rate^2.
func (e Exponential) Variance() float64 { return 1 / (e.Rate * e.Rate) }

// PDF returns the density Rate * e^{-Rate x} for x >= 0 and zero otherwise.
func (e Exponential) PDF(x float64) float64 {
	if x < 0 || !e.Valid() {
		return 0
	}
	return e.Rate * math.Exp(-e.Rate*x)
}

// NaturalParameter returns the natural parameter -Rate of the exponential
// family. It returns ErrDomain when Rate is not positive.
func (e Exponential) NaturalParameter() (float64, error) {
	if !e.Valid() {
		return 0, ErrDomain
	}
	return -e.Rate, nil
}

// FisherInformation returns the scalar Fisher information 1/Rate^2 of the
// exponential family in the rate parameterisation. It returns ErrDomain when
// Rate is not positive.
func (e Exponential) FisherInformation() (float64, error) {
	if !e.Valid() {
		return 0, ErrDomain
	}
	return 1 / (e.Rate * e.Rate), nil
}

// Gaussian is the univariate normal distribution with mean Mu and standard
// deviation Sigma > 0.
type Gaussian struct {
	// Mu is the mean.
	Mu float64
	// Sigma is the standard deviation.
	Sigma float64
}

// Valid reports whether the standard deviation is positive.
func (g Gaussian) Valid() bool { return g.Sigma > 0 }

// Mean returns the mean Mu.
func (g Gaussian) Mean() float64 { return g.Mu }

// Variance returns the variance Sigma^2.
func (g Gaussian) Variance() float64 { return g.Sigma * g.Sigma }

// PDF returns the normal density at x.
func (g Gaussian) PDF(x float64) float64 {
	if !g.Valid() {
		return 0
	}
	z := (x - g.Mu) / g.Sigma
	return math.Exp(-0.5*z*z) / (g.Sigma * math.Sqrt(2*math.Pi))
}

// DifferentialEntropy returns the differential entropy 1/2 ln(2 pi e Sigma^2)
// of the Gaussian in nats. It returns ErrDomain when Sigma is not positive.
func (g Gaussian) DifferentialEntropy() (float64, error) {
	if !g.Valid() {
		return 0, ErrDomain
	}
	return 0.5 * math.Log(2*math.Pi*math.E*g.Sigma*g.Sigma), nil
}

// FisherInformationMuSigma returns the 2x2 Fisher information matrix of the
// Gaussian family in the (Mu, Sigma) parameterisation, diag(1/Sigma^2,
// 2/Sigma^2). It returns ErrDomain when Sigma is not positive.
func (g Gaussian) FisherInformationMuSigma() ([][]float64, error) {
	if !g.Valid() {
		return nil, ErrDomain
	}
	s2 := g.Sigma * g.Sigma
	return [][]float64{{1 / s2, 0}, {0, 2 / s2}}, nil
}

// FisherInformationMuVar returns the 2x2 Fisher information matrix of the
// Gaussian family in the (Mu, Sigma^2) parameterisation, diag(1/Sigma^2,
// 1/(2 Sigma^4)). It returns ErrDomain when Sigma is not positive.
func (g Gaussian) FisherInformationMuVar() ([][]float64, error) {
	if !g.Valid() {
		return nil, ErrDomain
	}
	s2 := g.Sigma * g.Sigma
	return [][]float64{{1 / s2, 0}, {0, 1 / (2 * s2 * s2)}}, nil
}

// KLDivergenceGaussian returns the Kullback-Leibler divergence D(g||h) in nats
// between two univariate Gaussians,
//
//	ln(sigma_h/sigma_g) + (sigma_g^2 + (mu_g-mu_h)^2)/(2 sigma_h^2) - 1/2.
//
// It returns ErrDomain when either standard deviation is not positive.
func KLDivergenceGaussian(g, h Gaussian) (float64, error) {
	if !g.Valid() || !h.Valid() {
		return 0, ErrDomain
	}
	dm := g.Mu - h.Mu
	return math.Log(h.Sigma/g.Sigma) +
		(g.Sigma*g.Sigma+dm*dm)/(2*h.Sigma*h.Sigma) - 0.5, nil
}
