package markov

import "math/rand"

// SampleCategorical draws an index from the categorical distribution p using
// rng. It assumes p is a non-negative vector; it is normalized implicitly by
// scaling the uniform draw by the total mass. It returns -1 for an empty or
// all-zero p.
func SampleCategorical(p []float64, rng *rand.Rand) int {
	if len(p) == 0 || rng == nil {
		return -1
	}
	var total float64
	for _, x := range p {
		if x > 0 {
			total += x
		}
	}
	if total <= 0 {
		return -1
	}
	u := rng.Float64() * total
	var cum float64
	for i, x := range p {
		if x > 0 {
			cum += x
			if u < cum {
				return i
			}
		}
	}
	// Fall back to the last positive index.
	for i := len(p) - 1; i >= 0; i-- {
		if p[i] > 0 {
			return i
		}
	}
	return -1
}

// SampleFromCDF returns the smallest index i such that cdf[i] >= u, given a
// non-decreasing cumulative distribution cdf and a uniform draw u in [0,1). It
// performs a binary search. It returns -1 for an empty cdf.
func SampleFromCDF(cdf []float64, u float64) int {
	n := len(cdf)
	if n == 0 {
		return -1
	}
	lo, hi := 0, n-1
	for lo < hi {
		mid := (lo + hi) / 2
		if cdf[mid] < u {
			lo = mid + 1
		} else {
			hi = mid
		}
	}
	return lo
}

// CumulativeDistribution returns the cumulative sums of the (non-negative)
// probability vector p, normalized so the final entry is 1.
func CumulativeDistribution(p []float64) []float64 {
	cdf := make([]float64, len(p))
	var cum float64
	for i, x := range p {
		cum += x
		cdf[i] = cum
	}
	if cum > 0 {
		for i := range cdf {
			cdf[i] /= cum
		}
	}
	return cdf
}

// AliasTable implements Walker's alias method for O(1) sampling from a fixed
// discrete distribution.
type AliasTable struct {
	prob  []float64
	alias []int
	n     int
}

// NewAliasTable builds an alias table for the discrete distribution weights
// (which need not be normalized). It returns nil if weights is empty or has no
// positive mass.
func NewAliasTable(weights []float64) *AliasTable {
	n := len(weights)
	if n == 0 {
		return nil
	}
	var total float64
	for _, w := range weights {
		if w < 0 {
			return nil
		}
		total += w
	}
	if total <= 0 {
		return nil
	}
	prob := make([]float64, n)
	alias := make([]int, n)
	scaled := make([]float64, n)
	for i, w := range weights {
		scaled[i] = w / total * float64(n)
	}
	var small, large []int
	for i, s := range scaled {
		if s < 1 {
			small = append(small, i)
		} else {
			large = append(large, i)
		}
	}
	for len(small) > 0 && len(large) > 0 {
		s := small[len(small)-1]
		small = small[:len(small)-1]
		l := large[len(large)-1]
		large = large[:len(large)-1]
		prob[s] = scaled[s]
		alias[s] = l
		scaled[l] = scaled[l] - (1 - scaled[s])
		if scaled[l] < 1 {
			small = append(small, l)
		} else {
			large = append(large, l)
		}
	}
	for _, l := range large {
		prob[l] = 1
	}
	for _, s := range small {
		prob[s] = 1
	}
	return &AliasTable{prob: prob, alias: alias, n: n}
}

// Sample draws an index from the alias table using rng.
func (t *AliasTable) Sample(rng *rand.Rand) int {
	if t == nil || t.n == 0 || rng == nil {
		return -1
	}
	i := rng.Intn(t.n)
	if rng.Float64() < t.prob[i] {
		return i
	}
	return t.alias[i]
}

// Len returns the number of outcomes in the alias table.
func (t *AliasTable) Len() int {
	if t == nil {
		return 0
	}
	return t.n
}
