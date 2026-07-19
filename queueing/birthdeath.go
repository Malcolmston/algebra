package queueing

import "math"

// BirthDeath models a finite continuous-time birth–death chain on states
// 0,1,...,N. Birth[i] is the transition rate from state i to i+1 (for
// i=0..N-1) and Death[i] is the rate from state i+1 to i (for i=0..N-1), so
// both slices have length N. Many elementary queues are special cases of this
// chain.
type BirthDeath struct {
	Birth []float64 // up rates: Birth[i] is rate i -> i+1
	Death []float64 // down rates: Death[i] is rate i+1 -> i
}

// NewBirthDeath constructs a [BirthDeath] chain from equal-length birth and
// death rate slices. It returns an error when the slices differ in length, are
// empty, or contain a negative rate or a non-positive death rate.
func NewBirthDeath(birth, death []float64) (BirthDeath, error) {
	if len(birth) != len(death) {
		return BirthDeath{}, ErrDimension
	}
	if len(birth) == 0 {
		return BirthDeath{}, ErrDimension
	}
	for i := range birth {
		if birth[i] < 0 {
			return BirthDeath{}, ErrNegative
		}
		if death[i] <= 0 {
			return BirthDeath{}, ErrNonPositiveRate
		}
	}
	b := make([]float64, len(birth))
	d := make([]float64, len(death))
	copy(b, birth)
	copy(d, death)
	return BirthDeath{Birth: b, Death: d}, nil
}

// States returns the number of states N+1 in the chain.
func (c BirthDeath) States() int { return len(c.Birth) + 1 }

// StationaryDistribution returns the normalized stationary distribution over
// states 0..N obtained from the detailed-balance product
// p_n = p_0 * prod_{i<n} Birth[i]/Death[i].
func (c BirthDeath) StationaryDistribution() []float64 {
	n := c.States()
	p := make([]float64, n)
	p[0] = 1
	total := 1.0
	for i := 1; i < n; i++ {
		p[i] = p[i-1] * c.Birth[i-1] / c.Death[i-1]
		total += p[i]
	}
	for i := range p {
		p[i] /= total
	}
	return p
}

// ProbN returns the stationary probability of state n, or 0 when n is outside
// the range 0..N.
func (c BirthDeath) ProbN(n int) float64 {
	if n < 0 || n >= c.States() {
		return 0
	}
	return c.StationaryDistribution()[n]
}

// MeanNumber returns the stationary mean state value sum_n n*p_n, i.e. the mean
// number in system.
func (c BirthDeath) MeanNumber() float64 {
	p := c.StationaryDistribution()
	sum := 0.0
	for n, pn := range p {
		sum += float64(n) * pn
	}
	return sum
}

// VarianceNumber returns the variance of the stationary state value.
func (c BirthDeath) VarianceNumber() float64 {
	p := c.StationaryDistribution()
	mean := 0.0
	for n, pn := range p {
		mean += float64(n) * pn
	}
	second := 0.0
	for n, pn := range p {
		d := float64(n) - mean
		second += d * d * pn
	}
	return second
}

// Throughput returns the stationary birth rate sum_n Birth[n]*p_n, which equals
// the long-run rate of upward transitions.
func (c BirthDeath) Throughput() float64 {
	p := c.StationaryDistribution()
	sum := 0.0
	for i, up := range c.Birth {
		sum += up * p[i]
	}
	return sum
}

// BirthDeathStationary returns the stationary distribution of a birth–death
// chain given standalone birth and death rate slices, or nil when the inputs
// are invalid (mismatched lengths, empty, or a non-positive death rate).
func BirthDeathStationary(birth, death []float64) []float64 {
	c, err := NewBirthDeath(birth, death)
	if err != nil {
		return nil
	}
	return c.StationaryDistribution()
}

// DetailedBalance reports whether a candidate distribution p satisfies the
// detailed-balance equations p[i]*Birth[i] == p[i+1]*Death[i] within the given
// absolute tolerance.
func (c BirthDeath) DetailedBalance(p []float64, tol float64) bool {
	if len(p) != c.States() {
		return false
	}
	for i := 0; i < len(c.Birth); i++ {
		if math.Abs(p[i]*c.Birth[i]-p[i+1]*c.Death[i]) > tol {
			return false
		}
	}
	return true
}
