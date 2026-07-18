package probability

import (
	"math"
	"math/cmplx"
)

// Expectation returns E[g(X)] = Σ_i g(Outcomes[i]) · Probs[i], the expected
// value of the arbitrary function g applied to the random variable.
func (d Distribution) Expectation(g func(float64) float64) float64 {
	sum := 0.0
	for i, o := range d.Outcomes {
		sum += g(o) * d.Probs[i]
	}
	return sum
}

// Mean returns the expected value E[X] of the distribution.
func (d Distribution) Mean() float64 {
	sum := 0.0
	for i, o := range d.Outcomes {
		sum += o * d.Probs[i]
	}
	return sum
}

// Moment returns the k-th raw moment E[X^k]. Moment(0) is one and Moment(1) is
// the mean. It returns NaN for negative k.
func (d Distribution) Moment(k int) float64 {
	if k < 0 {
		return math.NaN()
	}
	sum := 0.0
	for i, o := range d.Outcomes {
		sum += math.Pow(o, float64(k)) * d.Probs[i]
	}
	return sum
}

// CentralMoment returns the k-th central moment E[(X - E[X])^k]. CentralMoment(2)
// is the variance. It returns NaN for negative k.
func (d Distribution) CentralMoment(k int) float64 {
	if k < 0 {
		return math.NaN()
	}
	mean := d.Mean()
	sum := 0.0
	for i, o := range d.Outcomes {
		sum += math.Pow(o-mean, float64(k)) * d.Probs[i]
	}
	return sum
}

// Variance returns Var(X) = E[(X - E[X])^2], the population variance of the
// distribution.
func (d Distribution) Variance() float64 { return d.CentralMoment(2) }

// StdDev returns the standard deviation, the non-negative square root of the
// variance.
func (d Distribution) StdDev() float64 { return math.Sqrt(d.Variance()) }

// Skewness returns the standardized third central moment,
// E[(X - E[X])^3] / σ^3, a dimensionless measure of asymmetry. It returns NaN
// when the standard deviation is zero.
func (d Distribution) Skewness() float64 {
	sd := d.StdDev()
	if sd == 0 {
		return math.NaN()
	}
	return d.CentralMoment(3) / (sd * sd * sd)
}

// Kurtosis returns the excess kurtosis, E[(X - E[X])^4] / σ^4 - 3, which is zero
// for a normal distribution. It returns NaN when the standard deviation is zero.
func (d Distribution) Kurtosis() float64 {
	v := d.Variance()
	if v == 0 {
		return math.NaN()
	}
	return d.CentralMoment(4)/(v*v) - 3
}

// MGF returns the moment-generating function M(t) = E[e^{tX}] evaluated
// numerically at t. Derivatives of M at t = 0 give the raw moments of X.
func (d Distribution) MGF(t float64) float64 {
	sum := 0.0
	for i, o := range d.Outcomes {
		sum += math.Exp(t*o) * d.Probs[i]
	}
	return sum
}

// CGF returns the cumulant-generating function K(t) = ln E[e^{tX}] = ln M(t)
// evaluated at t. Its derivatives at zero give the cumulants of X.
func (d Distribution) CGF(t float64) float64 {
	return math.Log(d.MGF(t))
}

// PGF returns the probability-generating function G(z) = E[z^X] = Σ_i Probs[i]
// z^{Outcomes[i]} evaluated at z. It is most meaningful for distributions on the
// non-negative integers, where the k-th derivative at zero recovers k!·P(X=k).
func (d Distribution) PGF(z float64) float64 {
	sum := 0.0
	for i, o := range d.Outcomes {
		sum += d.Probs[i] * math.Pow(z, o)
	}
	return sum
}

// CharacteristicFunction returns φ(t) = E[e^{itX}], the characteristic function
// of X evaluated at real argument t, as a complex value.
func (d Distribution) CharacteristicFunction(t float64) complex128 {
	var sum complex128
	for i, o := range d.Outcomes {
		sum += complex(d.Probs[i], 0) * cmplx.Exp(complex(0, t*o))
	}
	return sum
}

// MGFDerivativeMoment returns the k-th raw moment estimated by numerically
// differentiating the moment-generating function k times at t = 0 using a
// central finite-difference stencil with step h. It is provided as a numerical
// cross-check of [Distribution.Moment]; for exact moments prefer that method.
// It returns NaN for negative k or non-positive h.
func (d Distribution) MGFDerivativeMoment(k int, h float64) float64 {
	if k < 0 || h <= 0 {
		return math.NaN()
	}
	if k == 0 {
		return 1
	}
	// k-th derivative via the central finite-difference formula:
	// f^{(k)}(0) ≈ h^{-k} Σ_{j=0}^{k} (-1)^j C(k,j) f((k/2 - j)h).
	sum := 0.0
	for j := 0; j <= k; j++ {
		c := probabilityBinomCoeff(k, j)
		sign := 1.0
		if j%2 == 1 {
			sign = -1.0
		}
		x := (float64(k)/2 - float64(j)) * h
		sum += sign * c * d.MGF(x)
	}
	return sum / math.Pow(h, float64(k))
}

// probabilityBinomCoeff returns the binomial coefficient C(n, k) as a float64.
func probabilityBinomCoeff(n, k int) float64 {
	if k < 0 || k > n {
		return 0
	}
	if k > n-k {
		k = n - k
	}
	res := 1.0
	for i := 0; i < k; i++ {
		res *= float64(n - i)
		res /= float64(i + 1)
	}
	return res
}

// Entropy returns the Shannon entropy H(X) = -Σ_i Probs[i] ln Probs[i] measured
// in nats. Outcomes with zero probability contribute nothing.
func (d Distribution) Entropy() float64 {
	sum := 0.0
	for _, p := range d.Probs {
		if p > 0 {
			sum -= p * math.Log(p)
		}
	}
	return sum
}

// EntropyBits returns the Shannon entropy measured in bits, i.e. the entropy in
// nats divided by ln 2.
func (d Distribution) EntropyBits() float64 {
	return d.Entropy() / math.Ln2
}
