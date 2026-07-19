package bayesian

import (
	"math"
)

// ------------------------------------------------------------------
// Poisson distribution
// ------------------------------------------------------------------

// Poisson is a Poisson distribution with rate Lambda ≥ 0.
type Poisson struct {
	Lambda float64
}

// NewPoisson constructs a Poisson distribution, returning ErrParam for a
// negative rate.
func NewPoisson(lambda float64) (Poisson, error) {
	if lambda < 0 {
		return Poisson{}, ErrParam
	}
	return Poisson{Lambda: lambda}, nil
}

// Mean returns Lambda.
func (d Poisson) Mean() float64 { return d.Lambda }

// Variance returns Lambda.
func (d Poisson) Variance() float64 { return d.Lambda }

// PMF returns the probability mass P(X = k) for integer k ≥ 0.
func (d Poisson) PMF(k int) float64 {
	if k < 0 {
		return 0
	}
	return math.Exp(d.LogPMF(k))
}

// LogPMF returns the natural logarithm of the mass at k ≥ 0.
func (d Poisson) LogPMF(k int) float64 {
	if k < 0 {
		return math.Inf(-1)
	}
	if d.Lambda == 0 {
		if k == 0 {
			return 0
		}
		return math.Inf(-1)
	}
	fk := float64(k)
	return fk*math.Log(d.Lambda) - d.Lambda - LogFactorial(k)
}

// CDF returns P(X ≤ k).
func (d Poisson) CDF(k int) float64 {
	if k < 0 {
		return 0
	}
	return RegularizedGammaQ(float64(k+1), d.Lambda)
}

// ------------------------------------------------------------------
// Binomial distribution
// ------------------------------------------------------------------

// Binomial is a Binomial distribution with N ≥ 0 trials and success
// probability P in [0,1].
type Binomial struct {
	N int
	P float64
}

// NewBinomial constructs a Binomial distribution, returning ErrParam for an
// out-of-range probability or negative N.
func NewBinomial(n int, p float64) (Binomial, error) {
	if n < 0 || p < 0 || p > 1 {
		return Binomial{}, ErrParam
	}
	return Binomial{N: n, P: p}, nil
}

// Mean returns N·P.
func (d Binomial) Mean() float64 { return float64(d.N) * d.P }

// Variance returns N·P·(1−P).
func (d Binomial) Variance() float64 { return float64(d.N) * d.P * (1 - d.P) }

// PMF returns the probability mass P(X = k).
func (d Binomial) PMF(k int) float64 {
	if k < 0 || k > d.N {
		return 0
	}
	return math.Exp(d.LogPMF(k))
}

// LogPMF returns the natural logarithm of the mass at k.
func (d Binomial) LogPMF(k int) float64 {
	if k < 0 || k > d.N {
		return math.Inf(-1)
	}
	fk := float64(k)
	fn := float64(d.N)
	lp := math.Log(d.P)
	lq := math.Log(1 - d.P)
	if d.P == 0 {
		if k == 0 {
			return 0
		}
		return math.Inf(-1)
	}
	if d.P == 1 {
		if k == d.N {
			return 0
		}
		return math.Inf(-1)
	}
	return LogChoose(d.N, k) + fk*lp + (fn-fk)*lq
}

// CDF returns P(X ≤ k).
func (d Binomial) CDF(k int) float64 {
	if k < 0 {
		return 0
	}
	if k >= d.N {
		return 1
	}
	return RegularizedIncompleteBeta(1-d.P, float64(d.N-k), float64(k+1))
}

// ------------------------------------------------------------------
// Beta-Binomial distribution (Binomial marginalized over a Beta prior)
// ------------------------------------------------------------------

// BetaBinomial is the Beta-Binomial distribution: the number of successes in N
// trials where the success probability is drawn from Beta(Alpha, Beta). It is
// the posterior predictive of a Beta-Bernoulli model.
type BetaBinomial struct {
	N     int
	Alpha float64
	Beta  float64
}

// NewBetaBinomial constructs a Beta-Binomial distribution, returning ErrParam
// for invalid parameters.
func NewBetaBinomial(n int, alpha, beta float64) (BetaBinomial, error) {
	if n < 0 || alpha <= 0 || beta <= 0 {
		return BetaBinomial{}, ErrParam
	}
	return BetaBinomial{N: n, Alpha: alpha, Beta: beta}, nil
}

// Mean returns N·α/(α+β).
func (d BetaBinomial) Mean() float64 {
	return float64(d.N) * d.Alpha / (d.Alpha + d.Beta)
}

// Variance returns the variance of the Beta-Binomial distribution.
func (d BetaBinomial) Variance() float64 {
	n := float64(d.N)
	a, b := d.Alpha, d.Beta
	s := a + b
	return n * a * b * (s + n) / (s * s * (s + 1))
}

// PMF returns the probability mass P(X = k) for k in [0,N].
func (d BetaBinomial) PMF(k int) float64 {
	if k < 0 || k > d.N {
		return 0
	}
	return math.Exp(d.LogPMF(k))
}

// LogPMF returns the natural logarithm of the mass at k.
func (d BetaBinomial) LogPMF(k int) float64 {
	if k < 0 || k > d.N {
		return math.Inf(-1)
	}
	return LogChoose(d.N, k) +
		LogBeta(float64(k)+d.Alpha, float64(d.N-k)+d.Beta) -
		LogBeta(d.Alpha, d.Beta)
}

// CDF returns P(X ≤ k) computed by summation over the support.
func (d BetaBinomial) CDF(k int) float64 {
	if k < 0 {
		return 0
	}
	if k >= d.N {
		return 1
	}
	var sum float64
	for i := 0; i <= k; i++ {
		sum += d.PMF(i)
	}
	return sum
}

// ------------------------------------------------------------------
// Negative-Binomial distribution (Gamma-Poisson mixture)
// ------------------------------------------------------------------

// NegativeBinomial is parameterized by a real number R > 0 of "failures" and a
// success probability P in (0,1); the mass is over the number of counts k ≥ 0.
// It is the posterior predictive of a Gamma-Poisson model, where R equals the
// Gamma shape and P = Rate/(Rate+1).
type NegativeBinomial struct {
	R float64
	P float64
}

// NewNegativeBinomial constructs a Negative-Binomial distribution, returning
// ErrParam for invalid parameters.
func NewNegativeBinomial(r, p float64) (NegativeBinomial, error) {
	if r <= 0 || p <= 0 || p >= 1 {
		return NegativeBinomial{}, ErrParam
	}
	return NegativeBinomial{R: r, P: p}, nil
}

// Mean returns R·(1−P)/P.
func (d NegativeBinomial) Mean() float64 {
	return d.R * (1 - d.P) / d.P
}

// Variance returns R·(1−P)/P².
func (d NegativeBinomial) Variance() float64 {
	return d.R * (1 - d.P) / (d.P * d.P)
}

// PMF returns the probability mass P(X = k) for k ≥ 0.
func (d NegativeBinomial) PMF(k int) float64 {
	if k < 0 {
		return 0
	}
	return math.Exp(d.LogPMF(k))
}

// LogPMF returns the natural logarithm of the mass at k ≥ 0 using the
// generalized (real-valued R) binomial coefficient.
func (d NegativeBinomial) LogPMF(k int) float64 {
	if k < 0 {
		return math.Inf(-1)
	}
	fk := float64(k)
	// C(k+r-1, k) = Γ(k+r)/(k! Γ(r))
	logCoef := LogGamma(fk+d.R) - LogFactorial(k) - LogGamma(d.R)
	return logCoef + d.R*math.Log(d.P) + fk*math.Log(1-d.P)
}

// CDF returns P(X ≤ k) via the regularized incomplete beta function.
func (d NegativeBinomial) CDF(k int) float64 {
	if k < 0 {
		return 0
	}
	return RegularizedIncompleteBeta(d.P, d.R, float64(k)+1)
}
