package stochastic

import (
	"math"
	"math/rand"
)

// NewRNG returns a new deterministic pseudo-random source seeded with the given
// value. Two generators created with the same seed produce identical streams,
// which makes every simulation in this package reproducible.
func NewRNG(seed int64) *rand.Rand {
	return rand.New(rand.NewSource(seed))
}

// UniformSample returns a sample from the continuous uniform distribution on
// the half-open interval [a, b). If a > b the endpoints are swapped.
func UniformSample(rng *rand.Rand, a, b float64) float64 {
	if a > b {
		a, b = b, a
	}
	return a + (b-a)*rng.Float64()
}

// UniformIntSample returns a uniformly distributed integer in the closed
// interval [a, b]. If a > b the endpoints are swapped.
func UniformIntSample(rng *rand.Rand, a, b int) int {
	if a > b {
		a, b = b, a
	}
	return a + rng.Intn(b-a+1)
}

// StandardNormalSample returns a sample from the standard normal distribution
// N(0, 1).
func StandardNormalSample(rng *rand.Rand) float64 {
	return rng.NormFloat64()
}

// NormalSample returns a sample from the normal distribution N(mu, sigma^2)
// with mean mu and standard deviation sigma. Negative sigma is treated as its
// absolute value.
func NormalSample(rng *rand.Rand, mu, sigma float64) float64 {
	return mu + math.Abs(sigma)*rng.NormFloat64()
}

// BoxMullerPair returns two independent standard normal variates generated from
// a single pair of uniforms using the Box-Muller transform.
func BoxMullerPair(rng *rand.Rand) (float64, float64) {
	u1 := 1 - rng.Float64()
	u2 := rng.Float64()
	r := math.Sqrt(-2 * math.Log(u1))
	return r * math.Cos(2*math.Pi*u2), r * math.Sin(2*math.Pi*u2)
}

// ExponentialSample returns a sample from the exponential distribution with the
// given rate (inverse mean). The mean of the distribution is 1/rate.
func ExponentialSample(rng *rand.Rand, rate float64) float64 {
	if rate <= 0 {
		return math.Inf(1)
	}
	return -math.Log(1-rng.Float64()) / rate
}

// BernoulliSample returns true with probability p and false otherwise. Values
// of p outside [0, 1] are clamped.
func BernoulliSample(rng *rand.Rand, p float64) bool {
	if p <= 0 {
		return false
	}
	if p >= 1 {
		return true
	}
	return rng.Float64() < p
}

// RademacherSample returns +1 or -1, each with probability 1/2.
func RademacherSample(rng *rand.Rand) int {
	if rng.Float64() < 0.5 {
		return -1
	}
	return 1
}

// GeometricSample returns the number of Bernoulli(p) trials up to and including
// the first success. The result is at least 1 and has mean 1/p.
func GeometricSample(rng *rand.Rand, p float64) int {
	if p >= 1 {
		return 1
	}
	if p <= 0 {
		return math.MaxInt
	}
	return int(math.Floor(math.Log(1-rng.Float64())/math.Log(1-p))) + 1
}

// BinomialSample returns a sample from the binomial distribution Binomial(n, p),
// the number of successes in n independent Bernoulli(p) trials.
func BinomialSample(rng *rand.Rand, n int, p float64) int {
	if n <= 0 {
		return 0
	}
	if p <= 0 {
		return 0
	}
	if p >= 1 {
		return n
	}
	k := 0
	for i := 0; i < n; i++ {
		if rng.Float64() < p {
			k++
		}
	}
	return k
}

// PoissonSample returns a sample from the Poisson distribution with mean lambda.
// Knuth's multiplication method is used for small means and Hoermann's
// transformed-rejection method for large means, so the routine is exact and
// fast across the whole range.
func PoissonSample(rng *rand.Rand, lambda float64) int {
	if lambda <= 0 {
		return 0
	}
	if lambda < 30 {
		l := math.Exp(-lambda)
		k := 0
		p := 1.0
		for {
			k++
			p *= rng.Float64()
			if p <= l {
				return k - 1
			}
		}
	}
	return poissonLarge(rng, lambda)
}

// poissonLarge implements the PTRS transformed-rejection algorithm of Hoermann
// (1993) for lambda >= 10.
func poissonLarge(rng *rand.Rand, lambda float64) int {
	b := 0.931 + 2.53*math.Sqrt(lambda)
	a := -0.059 + 0.02483*b
	invAlpha := 1.1239 + 1.1328/(b-3.4)
	vr := 0.9277 - 3.6224/(b-2)
	for {
		u := rng.Float64() - 0.5
		v := rng.Float64()
		us := 0.5 - math.Abs(u)
		k := math.Floor((2*a/us+b)*u + lambda + 0.43)
		if us >= 0.07 && v <= vr {
			return int(k)
		}
		if k < 0 || (us < 0.013 && v > us) {
			continue
		}
		lg, _ := math.Lgamma(k + 1)
		if math.Log(v)+math.Log(invAlpha)-math.Log(a/(us*us)+b) <= -lambda+k*math.Log(lambda)-lg {
			return int(k)
		}
	}
}

// GammaSample returns a sample from the gamma distribution with the given shape
// (k) and scale (theta) parameters, using the Marsaglia-Tsang method. The mean
// is shape*scale.
func GammaSample(rng *rand.Rand, shape, scale float64) float64 {
	if shape <= 0 || scale <= 0 {
		return 0
	}
	if shape < 1 {
		u := rng.Float64()
		return GammaSample(rng, shape+1, scale) * math.Pow(u, 1/shape)
	}
	d := shape - 1.0/3.0
	c := 1 / math.Sqrt(9*d)
	for {
		x := rng.NormFloat64()
		v := 1 + c*x
		if v <= 0 {
			continue
		}
		v = v * v * v
		u := rng.Float64()
		if u < 1-0.0331*x*x*x*x {
			return d * v * scale
		}
		if math.Log(u) < 0.5*x*x+d*(1-v+math.Log(v)) {
			return d * v * scale
		}
	}
}

// BetaSample returns a sample from the beta distribution with shape parameters
// alpha and beta on the interval (0, 1).
func BetaSample(rng *rand.Rand, alpha, beta float64) float64 {
	x := GammaSample(rng, alpha, 1)
	y := GammaSample(rng, beta, 1)
	if x+y == 0 {
		return 0
	}
	return x / (x + y)
}

// ChiSquaredSample returns a sample from the chi-squared distribution with k
// degrees of freedom.
func ChiSquaredSample(rng *rand.Rand, k float64) float64 {
	return GammaSample(rng, k/2, 2)
}

// StudentTSample returns a sample from Student's t distribution with nu degrees
// of freedom.
func StudentTSample(rng *rand.Rand, nu float64) float64 {
	z := rng.NormFloat64()
	v := ChiSquaredSample(rng, nu)
	return z / math.Sqrt(v/nu)
}

// CauchySample returns a sample from the Cauchy distribution with the given
// location and scale.
func CauchySample(rng *rand.Rand, location, scale float64) float64 {
	return location + scale*math.Tan(math.Pi*(rng.Float64()-0.5))
}

// LaplaceSample returns a sample from the Laplace (double-exponential)
// distribution with the given location and scale (b).
func LaplaceSample(rng *rand.Rand, location, scale float64) float64 {
	u := rng.Float64() - 0.5
	return location - scale*sign(u)*math.Log(1-2*math.Abs(u))
}

// LogNormalSample returns a sample from the log-normal distribution whose
// natural logarithm is N(mu, sigma^2).
func LogNormalSample(rng *rand.Rand, mu, sigma float64) float64 {
	return math.Exp(mu + math.Abs(sigma)*rng.NormFloat64())
}

// WeibullSample returns a sample from the Weibull distribution with shape k and
// scale lambda.
func WeibullSample(rng *rand.Rand, k, lambda float64) float64 {
	if k <= 0 || lambda <= 0 {
		return 0
	}
	return lambda * math.Pow(-math.Log(1-rng.Float64()), 1/k)
}

// RayleighSample returns a sample from the Rayleigh distribution with the given
// scale parameter sigma.
func RayleighSample(rng *rand.Rand, sigma float64) float64 {
	return sigma * math.Sqrt(-2*math.Log(1-rng.Float64()))
}

// ParetoSample returns a sample from the Pareto distribution with minimum value
// xm and tail index alpha.
func ParetoSample(rng *rand.Rand, xm, alpha float64) float64 {
	if xm <= 0 || alpha <= 0 {
		return 0
	}
	return xm * math.Pow(1-rng.Float64(), -1/alpha)
}

// DiscreteSample returns an index in [0, len(weights)) chosen with probability
// proportional to the corresponding weight. Negative weights are treated as
// zero. It returns -1 only when every weight is non-positive.
func DiscreteSample(rng *rand.Rand, weights []float64) int {
	total := 0.0
	for _, w := range weights {
		if w > 0 {
			total += w
		}
	}
	if total <= 0 {
		return -1
	}
	u := rng.Float64() * total
	acc := 0.0
	for i, w := range weights {
		if w > 0 {
			acc += w
			if u < acc {
				return i
			}
		}
	}
	return len(weights) - 1
}

// CategoricalSample is a convenience wrapper around DiscreteSample for
// probability vectors that already sum to one.
func CategoricalSample(rng *rand.Rand, probs []float64) int {
	return DiscreteSample(rng, probs)
}

// NormalVector returns a slice of n independent standard normal variates.
func NormalVector(rng *rand.Rand, n int) []float64 {
	if n < 0 {
		n = 0
	}
	v := make([]float64, n)
	for i := range v {
		v[i] = rng.NormFloat64()
	}
	return v
}

// UniformVector returns a slice of n independent uniform variates on [a, b).
func UniformVector(rng *rand.Rand, n int, a, b float64) []float64 {
	if n < 0 {
		n = 0
	}
	v := make([]float64, n)
	for i := range v {
		v[i] = UniformSample(rng, a, b)
	}
	return v
}

// ExponentialVector returns a slice of n independent exponential variates with
// the given rate.
func ExponentialVector(rng *rand.Rand, n int, rate float64) []float64 {
	if n < 0 {
		n = 0
	}
	v := make([]float64, n)
	for i := range v {
		v[i] = ExponentialSample(rng, rate)
	}
	return v
}

// Permutation returns a uniformly random permutation of the integers
// [0, n) using the Fisher-Yates shuffle.
func Permutation(rng *rand.Rand, n int) []int {
	if n < 0 {
		n = 0
	}
	p := make([]int, n)
	for i := range p {
		p[i] = i
	}
	for i := n - 1; i > 0; i-- {
		j := rng.Intn(i + 1)
		p[i], p[j] = p[j], p[i]
	}
	return p
}

// ShuffleFloat randomly permutes the elements of x in place using the
// Fisher-Yates shuffle.
func ShuffleFloat(rng *rand.Rand, x []float64) {
	for i := len(x) - 1; i > 0; i-- {
		j := rng.Intn(i + 1)
		x[i], x[j] = x[j], x[i]
	}
}

func sign(x float64) float64 {
	switch {
	case x > 0:
		return 1
	case x < 0:
		return -1
	default:
		return 0
	}
}
