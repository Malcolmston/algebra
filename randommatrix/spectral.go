package randommatrix

import (
	"math"
	"sort"
)

// SortedCopy returns a sorted (ascending) copy of xs, leaving the input intact.
func SortedCopy(xs []float64) []float64 {
	out := make([]float64, len(xs))
	copy(out, xs)
	sort.Float64s(out)
	return out
}

// Min returns the minimum of a non-empty slice; it returns NaN for the empty
// slice.
func Min(xs []float64) float64 {
	if len(xs) == 0 {
		return math.NaN()
	}
	m := xs[0]
	for _, v := range xs[1:] {
		if v < m {
			m = v
		}
	}
	return m
}

// Max returns the maximum of a non-empty slice; it returns NaN for the empty
// slice.
func Max(xs []float64) float64 {
	if len(xs) == 0 {
		return math.NaN()
	}
	m := xs[0]
	for _, v := range xs[1:] {
		if v > m {
			m = v
		}
	}
	return m
}

// Mean returns the arithmetic mean of xs, or NaN for the empty slice.
func Mean(xs []float64) float64 {
	if len(xs) == 0 {
		return math.NaN()
	}
	var s float64
	for _, v := range xs {
		s += v
	}
	return s / float64(len(xs))
}

// Histogram is a fixed-width histogram of a real sample. Edges has NumBins+1
// entries; Counts[i] counts samples in [Edges[i], Edges[i+1]) except the last
// bin, which is closed on the right.
type Histogram struct {
	Edges  []float64
	Counts []int
	Total  int
}

// NewHistogram bins the data into the given number of equal-width bins spanning
// [lo, hi]. Values outside the range are ignored. It panics if bins < 1 or
// hi <= lo.
func NewHistogram(data []float64, bins int, lo, hi float64) *Histogram {
	if bins < 1 {
		panic("randommatrix: histogram needs at least one bin")
	}
	if hi <= lo {
		panic("randommatrix: histogram needs hi > lo")
	}
	edges := make([]float64, bins+1)
	w := (hi - lo) / float64(bins)
	for i := 0; i <= bins; i++ {
		edges[i] = lo + float64(i)*w
	}
	counts := make([]int, bins)
	total := 0
	for _, v := range data {
		if v < lo || v > hi {
			continue
		}
		idx := int((v - lo) / w)
		if idx >= bins {
			idx = bins - 1
		}
		if idx < 0 {
			idx = 0
		}
		counts[idx]++
		total++
	}
	return &Histogram{Edges: edges, Counts: counts, Total: total}
}

// NewHistogramAuto bins the data into the given number of bins spanning its
// observed range [min, max]. It panics if the data is empty or degenerate.
func NewHistogramAuto(data []float64, bins int) *Histogram {
	lo, hi := Min(data), Max(data)
	if hi == lo {
		hi = lo + 1
	}
	return NewHistogram(data, bins, lo, hi)
}

// NumBins returns the number of bins.
func (h *Histogram) NumBins() int { return len(h.Counts) }

// BinWidth returns the common width of the bins.
func (h *Histogram) BinWidth() float64 {
	if len(h.Edges) < 2 {
		return 0
	}
	return h.Edges[1] - h.Edges[0]
}

// BinCenters returns the midpoint of each bin.
func (h *Histogram) BinCenters() []float64 {
	c := make([]float64, len(h.Counts))
	for i := range c {
		c[i] = 0.5 * (h.Edges[i] + h.Edges[i+1])
	}
	return c
}

// Density returns the empirical probability density: Counts[i] divided by
// (Total * BinWidth). The returned values integrate to one over the range when
// no samples were discarded.
func (h *Histogram) Density() []float64 {
	d := make([]float64, len(h.Counts))
	w := h.BinWidth()
	if h.Total == 0 || w == 0 {
		return d
	}
	for i, c := range h.Counts {
		d[i] = float64(c) / (float64(h.Total) * w)
	}
	return d
}

// Probabilities returns the fraction of samples in each bin.
func (h *Histogram) Probabilities() []float64 {
	p := make([]float64, len(h.Counts))
	if h.Total == 0 {
		return p
	}
	for i, c := range h.Counts {
		p[i] = float64(c) / float64(h.Total)
	}
	return p
}

// ModeBin returns the index of the bin with the most samples.
func (h *Histogram) ModeBin() int {
	best, bi := -1, 0
	for i, c := range h.Counts {
		if c > best {
			best = c
			bi = i
		}
	}
	return bi
}

// EmpiricalSpectralDensity returns a histogram of the eigenvalue sample using
// the given number of bins over the observed spectral range. It is a synonym
// for NewHistogramAuto tailored to spectra.
func EmpiricalSpectralDensity(eigs []float64, bins int) *Histogram {
	return NewHistogramAuto(eigs, bins)
}

// EmpiricalCDF returns the value of the empirical cumulative distribution
// function of the sample at x, i.e. the fraction of samples not exceeding x.
func EmpiricalCDF(data []float64, x float64) float64 {
	if len(data) == 0 {
		return math.NaN()
	}
	var c int
	for _, v := range data {
		if v <= x {
			c++
		}
	}
	return float64(c) / float64(len(data))
}

// KolmogorovSmirnovStatistic returns the maximum absolute difference between the
// empirical CDF of the sample and the reference CDF, evaluated at the sample
// points.
func KolmogorovSmirnovStatistic(data []float64, cdf func(float64) float64) float64 {
	if len(data) == 0 {
		return math.NaN()
	}
	s := SortedCopy(data)
	n := float64(len(s))
	var d float64
	for i, x := range s {
		fx := cdf(x)
		lo := math.Abs(float64(i)/n - fx)
		hi := math.Abs(float64(i+1)/n - fx)
		if lo > d {
			d = lo
		}
		if hi > d {
			d = hi
		}
	}
	return d
}

// SpectralGaps returns the consecutive differences of the sorted eigenvalues.
func SpectralGaps(eigs []float64) []float64 {
	s := SortedCopy(eigs)
	if len(s) < 2 {
		return []float64{}
	}
	g := make([]float64, len(s)-1)
	for i := 1; i < len(s); i++ {
		g[i-1] = s[i] - s[i-1]
	}
	return g
}

// SpectralQuantile returns the q-quantile (0 <= q <= 1) of the eigenvalue sample
// using linear interpolation between order statistics.
func SpectralQuantile(eigs []float64, q float64) float64 {
	if len(eigs) == 0 {
		return math.NaN()
	}
	s := SortedCopy(eigs)
	if q <= 0 {
		return s[0]
	}
	if q >= 1 {
		return s[len(s)-1]
	}
	pos := q * float64(len(s)-1)
	lo := int(math.Floor(pos))
	frac := pos - float64(lo)
	if lo+1 >= len(s) {
		return s[lo]
	}
	return s[lo]*(1-frac) + s[lo+1]*frac
}
