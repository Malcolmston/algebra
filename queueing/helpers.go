package queueing

import (
	"errors"
	"math"
)

// Common sentinel errors returned by the constructors in this package.
var (
	// ErrNonPositiveRate indicates a non-positive arrival or service rate.
	ErrNonPositiveRate = errors.New("queueing: rate must be positive")
	// ErrServers indicates an invalid number of servers.
	ErrServers = errors.New("queueing: number of servers must be positive")
	// ErrCapacity indicates an invalid system capacity.
	ErrCapacity = errors.New("queueing: capacity is invalid")
	// ErrUnstable indicates the model does not admit a steady state.
	ErrUnstable = errors.New("queueing: system is unstable")
	// ErrDimension indicates mismatched slice dimensions.
	ErrDimension = errors.New("queueing: dimension mismatch")
	// ErrNegative indicates an unexpected negative parameter.
	ErrNegative = errors.New("queueing: parameter must be non-negative")
)

// Factorial returns n! as a float64. It returns NaN for negative n and +Inf
// once the exact value overflows float64 (for n greater than 170).
func Factorial(n int) float64 {
	if n < 0 {
		return math.NaN()
	}
	f := 1.0
	for i := 2; i <= n; i++ {
		f *= float64(i)
	}
	return f
}

// LogFactorial returns the natural logarithm of n! using the log-gamma
// function, remaining finite for large n. It returns NaN for negative n.
func LogFactorial(n int) float64 {
	if n < 0 {
		return math.NaN()
	}
	lg, _ := math.Lgamma(float64(n + 1))
	return lg
}

// OfferedLoad returns the offered load (traffic intensity in erlangs)
// a = lambda/mu, the mean number of busy servers in an unrestricted system.
// It returns NaN when mu is not positive or lambda is negative.
func OfferedLoad(lambda, mu float64) float64 {
	if mu <= 0 || lambda < 0 {
		return math.NaN()
	}
	return lambda / mu
}

// Utilization returns the per-server utilization rho = lambda/(c*mu) of a
// c-server queue. It returns NaN for invalid arguments.
func Utilization(lambda, mu float64, c int) float64 {
	if mu <= 0 || lambda < 0 || c <= 0 {
		return math.NaN()
	}
	return lambda / (float64(c) * mu)
}

// SquaredCoefficientOfVariation returns the squared coefficient of variation
// (SCV) variance/mean^2 of a distribution given its mean and variance. It
// returns NaN when the mean is zero or an argument is negative.
func SquaredCoefficientOfVariation(mean, variance float64) float64 {
	if mean == 0 || variance < 0 {
		return math.NaN()
	}
	return variance / (mean * mean)
}

// PoissonPMF returns the Poisson probability mass P(N=k) with mean lambda,
// computed in log space for numerical stability. It returns NaN for lambda<0
// and 0 for negative k.
func PoissonPMF(k int, lambda float64) float64 {
	if lambda < 0 {
		return math.NaN()
	}
	if k < 0 {
		return 0
	}
	if lambda == 0 {
		if k == 0 {
			return 1
		}
		return 0
	}
	logp := float64(k)*math.Log(lambda) - lambda - LogFactorial(k)
	return math.Exp(logp)
}

// PoissonCDF returns the Poisson cumulative probability P(N<=k) with mean
// lambda. It returns NaN for lambda<0 and 0 for k<0.
func PoissonCDF(k int, lambda float64) float64 {
	if lambda < 0 {
		return math.NaN()
	}
	if k < 0 {
		return 0
	}
	sum := 0.0
	for i := 0; i <= k; i++ {
		sum += PoissonPMF(i, lambda)
	}
	if sum > 1 {
		sum = 1
	}
	return sum
}

// ExponentialTail returns the exponential survival probability
// P(X>t)=exp(-rate*t) for t>=0. It returns NaN for a non-positive rate and 1
// for t<=0.
func ExponentialTail(rate, t float64) float64 {
	if rate <= 0 {
		return math.NaN()
	}
	if t <= 0 {
		return 1
	}
	return math.Exp(-rate * t)
}

// powFactRatio returns a^n / n! computed iteratively to avoid intermediate
// overflow. a must be non-negative and n non-negative.
func powFactRatio(a float64, n int) float64 {
	term := 1.0
	for i := 1; i <= n; i++ {
		term *= a / float64(i)
	}
	return term
}
