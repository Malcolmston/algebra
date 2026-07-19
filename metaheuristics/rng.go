package metaheuristics

import (
	"math"
	"math/rand"
)

// RNG is a deterministic pseudo-random number generator used by every
// stochastic routine in the package. It wraps math/rand's source so that a
// given seed always produces the same sequence, independent of the wall clock.
type RNG struct {
	src *rand.Rand
}

// NewRNG returns a deterministic [RNG] seeded with the supplied value. Two
// generators created with the same seed produce identical sequences.
func NewRNG(seed int64) *RNG {
	return &RNG{src: rand.New(rand.NewSource(seed))}
}

// Float64 returns a pseudo-random float64 in the half-open interval [0, 1).
func (g *RNG) Float64() float64 {
	return g.src.Float64()
}

// Float64Range returns a pseudo-random float64 uniformly distributed in the
// half-open interval [lo, hi). If hi <= lo it returns lo.
func (g *RNG) Float64Range(lo, hi float64) float64 {
	if hi <= lo {
		return lo
	}
	return lo + g.src.Float64()*(hi-lo)
}

// Intn returns a pseudo-random int in [0, n). It panics if n <= 0.
func (g *RNG) Intn(n int) int {
	return g.src.Intn(n)
}

// IntRange returns a pseudo-random int in [lo, hi). If hi <= lo it returns lo.
func (g *RNG) IntRange(lo, hi int) int {
	if hi <= lo {
		return lo
	}
	return lo + g.src.Intn(hi-lo)
}

// NormFloat64 returns a normally distributed float64 with mean 0 and standard
// deviation 1.
func (g *RNG) NormFloat64() float64 {
	return g.src.NormFloat64()
}

// Gaussian returns a normally distributed float64 with the given mean and
// standard deviation.
func (g *RNG) Gaussian(mean, std float64) float64 {
	return mean + std*g.src.NormFloat64()
}

// Perm returns a pseudo-random permutation of the integers [0, n).
func (g *RNG) Perm(n int) []int {
	return g.src.Perm(n)
}

// Shuffle pseudo-randomly permutes the elements of s in place using the
// Fisher-Yates algorithm.
func (g *RNG) Shuffle(s []int) {
	for i := len(s) - 1; i > 0; i-- {
		j := g.src.Intn(i + 1)
		s[i], s[j] = s[j], s[i]
	}
}

// ShuffleFloat pseudo-randomly permutes the elements of s in place.
func (g *RNG) ShuffleFloat(s []float64) {
	for i := len(s) - 1; i > 0; i-- {
		j := g.src.Intn(i + 1)
		s[i], s[j] = s[j], s[i]
	}
}

// UniformVec returns a new vector of length b.Dim() with each coordinate drawn
// uniformly from the corresponding [lo, hi) interval of b.
func (g *RNG) UniformVec(b Bounds) []float64 {
	x := make([]float64, b.Dim())
	for i := range x {
		x[i] = g.Float64Range(b.Lower[i], b.Upper[i])
	}
	return x
}

// GaussianVec returns a new vector of length n with each coordinate drawn from
// a normal distribution of the given mean and standard deviation.
func (g *RNG) GaussianVec(n int, mean, std float64) []float64 {
	x := make([]float64, n)
	for i := range x {
		x[i] = g.Gaussian(mean, std)
	}
	return x
}

// Choice returns a pseudo-random element index in [0, n) selected with
// probability proportional to the supplied non-negative weights. If all
// weights are zero (or weights is empty) it falls back to a uniform choice.
// It panics if len(weights) != n is not required; n is taken from weights.
func (g *RNG) Choice(weights []float64) int {
	n := len(weights)
	if n == 0 {
		return 0
	}
	total := 0.0
	for _, w := range weights {
		if w > 0 {
			total += w
		}
	}
	if total <= 0 || math.IsInf(total, 0) || math.IsNaN(total) {
		return g.src.Intn(n)
	}
	t := g.src.Float64() * total
	acc := 0.0
	for i, w := range weights {
		if w > 0 {
			acc += w
			if t < acc {
				return i
			}
		}
	}
	return n - 1
}

// Sample returns k distinct indices drawn without replacement from [0, n).
// If k >= n it returns a full permutation of [0, n). It panics if k < 0.
func (g *RNG) Sample(n, k int) []int {
	if k < 0 {
		panic("metaheuristics: negative sample size")
	}
	p := g.src.Perm(n)
	if k >= n {
		return p
	}
	return p[:k]
}
