package gametheory

import (
	"errors"
	"math"
)

// MixedStrategy is a probability distribution over a player's pure strategies:
// entry i is the probability of playing pure strategy i. A valid mixed strategy
// has non-negative entries summing to one.
type MixedStrategy []float64

// ErrDistribution is returned when a slice is not a valid probability
// distribution.
var ErrDistribution = errors.New("gametheory: probabilities must be non-negative and sum to one")

// NewMixedStrategy validates probs as a probability distribution (non-negative
// entries summing to one within a small tolerance) and returns a copy of it as
// a MixedStrategy.
func NewMixedStrategy(probs []float64) (MixedStrategy, error) {
	if len(probs) == 0 {
		return nil, ErrDistribution
	}
	var sum float64
	for _, p := range probs {
		if p < -1e-12 {
			return nil, ErrDistribution
		}
		sum += p
	}
	if math.Abs(sum-1) > 1e-9 {
		return nil, ErrDistribution
	}
	return append(MixedStrategy(nil), probs...), nil
}

// UniformStrategy returns the uniform mixed strategy over n pure strategies,
// assigning probability 1/n to each. It panics if n is not positive.
func UniformStrategy(n int) MixedStrategy {
	if n <= 0 {
		panic("gametheory: UniformStrategy requires n > 0")
	}
	m := make(MixedStrategy, n)
	for i := range m {
		m[i] = 1 / float64(n)
	}
	return m
}

// PureStrategy returns the degenerate mixed strategy over n pure strategies
// that plays pure strategy i with probability one. It panics if i is out of
// range.
func PureStrategy(n, i int) MixedStrategy {
	if i < 0 || i >= n {
		panic("gametheory: PureStrategy index out of range")
	}
	m := make(MixedStrategy, n)
	m[i] = 1
	return m
}

// IsValid reports whether the receiver is a valid probability distribution:
// every entry is at least -tol and the entries sum to one within tol.
func (m MixedStrategy) IsValid(tol float64) bool {
	if len(m) == 0 {
		return false
	}
	var sum float64
	for _, p := range m {
		if p < -tol {
			return false
		}
		sum += p
	}
	return math.Abs(sum-1) <= tol
}

// Support returns the sorted indices of pure strategies played with probability
// greater than tol.
func (m MixedStrategy) Support(tol float64) []int {
	var s []int
	for i, p := range m {
		if p > tol {
			s = append(s, i)
		}
	}
	return s
}

// Entropy returns the Shannon entropy of the distribution in bits (base-2 log).
// Zero-probability outcomes contribute nothing.
func (m MixedStrategy) Entropy() float64 {
	var h float64
	for _, p := range m {
		if p > 0 {
			h -= p * math.Log2(p)
		}
	}
	return h
}

// RowBestResponses returns the sorted row-player pure strategies that maximize
// the row player's expected payoff against the column player's mixed strategy q,
// treating strategies within tol of the maximum as tied best responses.
func (g Game) RowBestResponses(q MixedStrategy, tol float64) []int {
	m := g.NumRowStrategies()
	vals := make([]float64, m)
	best := math.Inf(-1)
	for i := 0; i < m; i++ {
		var v float64
		for j := range g.Row[i] {
			v += g.Row[i][j] * q[j]
		}
		vals[i] = v
		if v > best {
			best = v
		}
	}
	var br []int
	for i := 0; i < m; i++ {
		if vals[i] >= best-tol {
			br = append(br, i)
		}
	}
	return br
}

// ColBestResponses returns the sorted column-player pure strategies that
// maximize the column player's expected payoff against the row player's mixed
// strategy p, treating strategies within tol of the maximum as tied best
// responses.
func (g Game) ColBestResponses(p MixedStrategy, tol float64) []int {
	n := g.NumColStrategies()
	vals := make([]float64, n)
	best := math.Inf(-1)
	for j := 0; j < n; j++ {
		var v float64
		for i := range g.Col {
			v += g.Col[i][j] * p[i]
		}
		vals[j] = v
		if v > best {
			best = v
		}
	}
	var br []int
	for j := 0; j < n; j++ {
		if vals[j] >= best-tol {
			br = append(br, j)
		}
	}
	return br
}

// BestResponseRowValue returns the maximum expected payoff the row player can
// obtain against the column player's mixed strategy q.
func (g Game) BestResponseRowValue(q MixedStrategy) float64 {
	best := math.Inf(-1)
	for i := range g.Row {
		var v float64
		for j := range g.Row[i] {
			v += g.Row[i][j] * q[j]
		}
		if v > best {
			best = v
		}
	}
	return best
}

// BestResponseColValue returns the maximum expected payoff the column player can
// obtain against the row player's mixed strategy p.
func (g Game) BestResponseColValue(p MixedStrategy) float64 {
	n := g.NumColStrategies()
	best := math.Inf(-1)
	for j := 0; j < n; j++ {
		var v float64
		for i := range g.Col {
			v += g.Col[i][j] * p[i]
		}
		if v > best {
			best = v
		}
	}
	return best
}
