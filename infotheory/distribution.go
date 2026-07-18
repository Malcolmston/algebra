package infotheory

import (
	"errors"
	"math"
	"sort"
)

// ErrNotNormalizable is returned when a set of weights or counts cannot be
// turned into a probability distribution because their total is not positive.
var ErrNotNormalizable = errors.New("infotheory: weights sum to a non-positive value")

// UniformDistribution returns the uniform probability distribution over n
// outcomes, a slice of length n whose every entry is 1/n. It returns nil for
// n <= 0.
func UniformDistribution(n int) []float64 {
	if n <= 0 {
		return nil
	}
	p := make([]float64, n)
	v := 1.0 / float64(n)
	for i := range p {
		p[i] = v
	}
	return p
}

// NormalizeWeights returns a probability distribution proportional to the
// supplied non-negative weights, that is each weight divided by their total. It
// returns ErrNotNormalizable if the weights do not sum to a positive value.
func NormalizeWeights(weights []float64) ([]float64, error) {
	var sum float64
	for _, w := range weights {
		sum += w
	}
	if sum <= 0 {
		return nil, ErrNotNormalizable
	}
	p := make([]float64, len(weights))
	for i, w := range weights {
		p[i] = w / sum
	}
	return p, nil
}

// NormalizeCounts returns a probability distribution proportional to the
// supplied integer counts, dividing each count by their total. It returns
// ErrNotNormalizable if the counts do not sum to a positive value.
func NormalizeCounts(counts []int) ([]float64, error) {
	var total int
	for _, c := range counts {
		total += c
	}
	if total <= 0 {
		return nil, ErrNotNormalizable
	}
	p := make([]float64, len(counts))
	for i, c := range counts {
		p[i] = float64(c) / float64(total)
	}
	return p, nil
}

// IsProbabilityDistribution reports whether p consists of non-negative entries
// summing to one within the absolute tolerance tol. A negative entry always
// makes the report false.
func IsProbabilityDistribution(p []float64, tol float64) bool {
	var sum float64
	for _, pi := range p {
		if pi < 0 {
			return false
		}
		sum += pi
	}
	return math.Abs(sum-1) <= tol
}

// EmpiricalDistribution returns the empirical probability distribution of the
// samples together with the slice of distinct values, sorted ascending, that
// index it. Entry i of the returned distribution is the fraction of samples
// equal to the i-th distinct value. It returns ErrNotNormalizable when samples
// is empty.
func EmpiricalDistribution(samples []float64) (values []float64, dist []float64, err error) {
	if len(samples) == 0 {
		return nil, nil, ErrNotNormalizable
	}
	counts := make(map[float64]int)
	for _, s := range samples {
		counts[s]++
	}
	values = make([]float64, 0, len(counts))
	for v := range counts {
		values = append(values, v)
	}
	sort.Float64s(values)
	dist = make([]float64, len(values))
	n := float64(len(samples))
	for i, v := range values {
		dist[i] = float64(counts[v]) / n
	}
	return values, dist, nil
}
