package stats

import "math"

// Normal is a normal (Gaussian) distribution with mean Mu and standard
// deviation Sigma. Sigma must be positive.
type Normal struct {
	Mu    float64 // Mu is the mean.
	Sigma float64 // Sigma is the standard deviation (> 0).
}

// PDF returns the probability density at x.
func (n Normal) PDF(x float64) float64 {
	z := (x - n.Mu) / n.Sigma
	return math.Exp(-0.5*z*z) / (n.Sigma * math.Sqrt(2*math.Pi))
}

// CDF returns the cumulative probability P(X <= x), evaluated with math.Erf.
func (n Normal) CDF(x float64) float64 {
	return 0.5 * math.Erfc(-(x-n.Mu)/(n.Sigma*math.Sqrt2))
}

// Quantile returns the inverse CDF: the value x such that CDF(x) = p, for p in
// [0, 1]. It uses a rational approximation refined by one Halley step.
func (n Normal) Quantile(p float64) float64 {
	return n.Mu + n.Sigma*normQuantile(p)
}

// Mean returns the mean of the distribution, Mu.
func (n Normal) Mean() float64 { return n.Mu }

// Variance returns the variance of the distribution, Sigma².
func (n Normal) Variance() float64 { return n.Sigma * n.Sigma }

// NormalPDF returns the density at x of a normal distribution with the given
// mean mu and standard deviation sigma. It is a convenience wrapper around
// Normal.PDF.
func NormalPDF(x, mu, sigma float64) float64 {
	return Normal{Mu: mu, Sigma: sigma}.PDF(x)
}

// NormalCDF returns the cumulative probability P(X <= x) for a normal
// distribution with the given mean mu and standard deviation sigma. It is a
// convenience wrapper around Normal.CDF.
func NormalCDF(x, mu, sigma float64) float64 {
	return Normal{Mu: mu, Sigma: sigma}.CDF(x)
}

// Binomial is a binomial distribution: the number of successes in N
// independent trials, each succeeding with probability P. N must be >= 0 and
// P must lie in [0, 1].
type Binomial struct {
	N int     // N is the number of independent trials (>= 0).
	P float64 // P is the per-trial success probability in [0, 1].
}

// PMF returns the probability mass P(X = k), the probability of exactly k
// successes. It is 0 for k outside [0, N].
func (b Binomial) PMF(k int) float64 {
	if k < 0 || k > b.N {
		return 0
	}
	if b.P <= 0 {
		if k == 0 {
			return 1
		}
		return 0
	}
	if b.P >= 1 {
		if k == b.N {
			return 1
		}
		return 0
	}
	logC := gammaLn(float64(b.N)+1) - gammaLn(float64(k)+1) - gammaLn(float64(b.N-k)+1)
	logP := logC + float64(k)*math.Log(b.P) + float64(b.N-k)*math.Log(1-b.P)
	return math.Exp(logP)
}

// CDF returns the cumulative probability P(X <= k).
func (b Binomial) CDF(k int) float64 {
	if k < 0 {
		return 0
	}
	if k >= b.N {
		return 1
	}
	sum := 0.0
	for i := 0; i <= k; i++ {
		sum += b.PMF(i)
	}
	return sum
}

// Mean returns the mean of the distribution, N·P.
func (b Binomial) Mean() float64 { return float64(b.N) * b.P }

// Variance returns the variance of the distribution, N·P·(1-P).
func (b Binomial) Variance() float64 { return float64(b.N) * b.P * (1 - b.P) }

// Poisson is a Poisson distribution with rate parameter Lambda > 0: the number
// of events in a fixed interval when events occur independently at a constant
// average rate.
type Poisson struct {
	Lambda float64 // Lambda is the rate parameter (> 0).
}

// PMF returns the probability mass P(X = k). It is 0 for negative k.
func (p Poisson) PMF(k int) float64 {
	if k < 0 {
		return 0
	}
	logP := float64(k)*math.Log(p.Lambda) - p.Lambda - gammaLn(float64(k)+1)
	return math.Exp(logP)
}

// CDF returns the cumulative probability P(X <= k).
func (p Poisson) CDF(k int) float64 {
	if k < 0 {
		return 0
	}
	sum := 0.0
	for i := 0; i <= k; i++ {
		sum += p.PMF(i)
	}
	return sum
}

// Mean returns the mean of the distribution, Lambda.
func (p Poisson) Mean() float64 { return p.Lambda }

// Variance returns the variance of the distribution, Lambda.
func (p Poisson) Variance() float64 { return p.Lambda }

// Uniform is a continuous uniform distribution on the closed interval [A, B]
// with A < B.
type Uniform struct {
	A float64 // A is the lower bound of the support.
	B float64 // B is the upper bound of the support (A < B).
}

// PDF returns the probability density at x: 1/(B-A) inside [A, B] and 0
// outside.
func (u Uniform) PDF(x float64) float64 {
	if x < u.A || x > u.B {
		return 0
	}
	return 1 / (u.B - u.A)
}

// CDF returns the cumulative probability P(X <= x).
func (u Uniform) CDF(x float64) float64 {
	switch {
	case x < u.A:
		return 0
	case x > u.B:
		return 1
	default:
		return (x - u.A) / (u.B - u.A)
	}
}

// Quantile returns the inverse CDF for p in [0, 1].
func (u Uniform) Quantile(p float64) float64 {
	if p < 0 || p > 1 {
		return math.NaN()
	}
	return u.A + p*(u.B-u.A)
}

// Mean returns the mean of the distribution, (A+B)/2.
func (u Uniform) Mean() float64 { return (u.A + u.B) / 2 }

// Variance returns the variance of the distribution, (B-A)²/12.
func (u Uniform) Variance() float64 {
	d := u.B - u.A
	return d * d / 12
}

// Exponential is an exponential distribution with rate parameter Lambda > 0,
// modelling the waiting time between events in a Poisson process.
type Exponential struct {
	Lambda float64 // Lambda is the rate parameter (> 0).
}

// PDF returns the probability density at x. It is 0 for x < 0.
func (e Exponential) PDF(x float64) float64 {
	if x < 0 {
		return 0
	}
	return e.Lambda * math.Exp(-e.Lambda*x)
}

// CDF returns the cumulative probability P(X <= x). It is 0 for x < 0.
func (e Exponential) CDF(x float64) float64 {
	if x < 0 {
		return 0
	}
	return 1 - math.Exp(-e.Lambda*x)
}

// Quantile returns the inverse CDF for p in [0, 1].
func (e Exponential) Quantile(p float64) float64 {
	if p < 0 || p > 1 {
		return math.NaN()
	}
	return -math.Log(1-p) / e.Lambda
}

// Mean returns the mean of the distribution, 1/Lambda.
func (e Exponential) Mean() float64 { return 1 / e.Lambda }

// Variance returns the variance of the distribution, 1/Lambda².
func (e Exponential) Variance() float64 { return 1 / (e.Lambda * e.Lambda) }

// StudentT is a Student's t distribution with Nu degrees of freedom (Nu > 0).
type StudentT struct {
	Nu float64 // Nu is the degrees of freedom (> 0).
}

// PDF returns the probability density at x.
func (t StudentT) PDF(x float64) float64 {
	nu := t.Nu
	lnC := gammaLn((nu+1)/2) - gammaLn(nu/2) - 0.5*math.Log(nu*math.Pi)
	return math.Exp(lnC) * math.Pow(1+x*x/nu, -(nu+1)/2)
}

// CDF returns the cumulative probability P(X <= x), evaluated with the
// regularized incomplete beta function.
func (t StudentT) CDF(x float64) float64 {
	nu := t.Nu
	ib := 0.5 * regularizedIncompleteBeta(nu/2, 0.5, nu/(nu+x*x))
	if x > 0 {
		return 1 - ib
	}
	return ib
}

// Mean returns the mean of the distribution, 0 for Nu > 1 and NaN otherwise.
func (t StudentT) Mean() float64 {
	if t.Nu > 1 {
		return 0
	}
	return math.NaN()
}

// Variance returns the variance of the distribution: Nu/(Nu-2) for Nu > 2,
// +Inf for 1 < Nu <= 2, and NaN otherwise.
func (t StudentT) Variance() float64 {
	switch {
	case t.Nu > 2:
		return t.Nu / (t.Nu - 2)
	case t.Nu > 1:
		return math.Inf(1)
	default:
		return math.NaN()
	}
}

// ChiSquared is a chi-squared distribution with K degrees of freedom (K > 0).
type ChiSquared struct {
	K float64 // K is the degrees of freedom (> 0).
}

// PDF returns the probability density at x. It is 0 for x < 0.
func (c ChiSquared) PDF(x float64) float64 {
	if x < 0 {
		return 0
	}
	if x == 0 {
		if c.K < 2 {
			return math.Inf(1)
		}
		if c.K == 2 {
			return 0.5
		}
		return 0
	}
	k := c.K
	lnP := (k/2-1)*math.Log(x) - x/2 - (k/2)*math.Log(2) - gammaLn(k/2)
	return math.Exp(lnP)
}

// CDF returns the cumulative probability P(X <= x), evaluated with the
// regularized lower incomplete gamma function.
func (c ChiSquared) CDF(x float64) float64 {
	if x < 0 {
		return 0
	}
	return regularizedGammaP(c.K/2, x/2)
}

// Mean returns the mean of the distribution, K.
func (c ChiSquared) Mean() float64 { return c.K }

// Variance returns the variance of the distribution, 2K.
func (c ChiSquared) Variance() float64 { return 2 * c.K }

// Gamma is a gamma distribution parameterized by Shape (k, > 0) and Scale
// (θ, > 0). With this parameterization the mean is Shape·Scale. A rate
// parameter λ corresponds to Scale = 1/λ.
type Gamma struct {
	Shape float64 // Shape is the shape parameter k (> 0).
	Scale float64 // Scale is the scale parameter theta (> 0).
}

// PDF returns the probability density at x. It is 0 for x < 0.
func (g Gamma) PDF(x float64) float64 {
	if x < 0 {
		return 0
	}
	if x == 0 {
		switch {
		case g.Shape < 1:
			return math.Inf(1)
		case g.Shape == 1:
			return 1 / g.Scale
		default:
			return 0
		}
	}
	k, th := g.Shape, g.Scale
	lnP := (k-1)*math.Log(x) - x/th - k*math.Log(th) - gammaLn(k)
	return math.Exp(lnP)
}

// CDF returns the cumulative probability P(X <= x), evaluated with the
// regularized lower incomplete gamma function.
func (g Gamma) CDF(x float64) float64 {
	if x < 0 {
		return 0
	}
	return regularizedGammaP(g.Shape, x/g.Scale)
}

// Mean returns the mean of the distribution, Shape·Scale.
func (g Gamma) Mean() float64 { return g.Shape * g.Scale }

// Variance returns the variance of the distribution, Shape·Scale².
func (g Gamma) Variance() float64 { return g.Shape * g.Scale * g.Scale }
