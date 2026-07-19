package clustering

import (
	"math"
	"sort"
)

// VectorAdd returns the element-wise sum a+b. The shorter length governs the
// result.
func VectorAdd(a, b []float64) []float64 {
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	out := make([]float64, n)
	for i := 0; i < n; i++ {
		out[i] = a[i] + b[i]
	}
	return out
}

// VectorSub returns the element-wise difference a-b.
func VectorSub(a, b []float64) []float64 {
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	out := make([]float64, n)
	for i := 0; i < n; i++ {
		out[i] = a[i] - b[i]
	}
	return out
}

// VectorScale returns v scaled by the scalar s.
func VectorScale(v []float64, s float64) []float64 {
	out := make([]float64, len(v))
	for i, x := range v {
		out[i] = x * s
	}
	return out
}

// VectorAddInPlace adds b into a element-wise, modifying and returning a.
func VectorAddInPlace(a, b []float64) []float64 {
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	for i := 0; i < n; i++ {
		a[i] += b[i]
	}
	return a
}

// Centroid returns the coordinate-wise mean of the given points. It returns nil
// for an empty set.
func Centroid(points [][]float64) []float64 {
	if len(points) == 0 {
		return nil
	}
	dim := len(points[0])
	c := make([]float64, dim)
	for _, p := range points {
		for j := 0; j < dim && j < len(p); j++ {
			c[j] += p[j]
		}
	}
	inv := 1 / float64(len(points))
	for j := range c {
		c[j] *= inv
	}
	return c
}

// WeightedCentroid returns the weighted coordinate-wise mean of points using the
// given nonnegative weights. It returns nil if the total weight is zero.
func WeightedCentroid(points [][]float64, weights []float64) []float64 {
	if len(points) == 0 {
		return nil
	}
	dim := len(points[0])
	c := make([]float64, dim)
	var total float64
	for i, p := range points {
		w := 1.0
		if i < len(weights) {
			w = weights[i]
		}
		total += w
		for j := 0; j < dim && j < len(p); j++ {
			c[j] += w * p[j]
		}
	}
	if total == 0 {
		return nil
	}
	for j := range c {
		c[j] /= total
	}
	return c
}

// Medoid returns the index of the point in data that minimises the total
// distance to all other points in the given index subset, using the supplied
// metric. If indices is nil, all points are considered.
func Medoid(data [][]float64, indices []int, metric Metric) int {
	if metric == nil {
		metric = Euclidean
	}
	if indices == nil {
		indices = make([]int, len(data))
		for i := range indices {
			indices[i] = i
		}
	}
	best := -1
	bestCost := math.Inf(1)
	for _, i := range indices {
		var cost float64
		for _, j := range indices {
			if i != j {
				cost += metric(data[i], data[j])
			}
		}
		if cost < bestCost {
			bestCost = cost
			best = i
		}
	}
	return best
}

// ColumnMeans returns the mean of each column (feature) of data.
func ColumnMeans(data [][]float64) []float64 {
	if len(data) == 0 {
		return nil
	}
	dim := len(data[0])
	m := make([]float64, dim)
	for _, row := range data {
		for j := 0; j < dim && j < len(row); j++ {
			m[j] += row[j]
		}
	}
	inv := 1 / float64(len(data))
	for j := range m {
		m[j] *= inv
	}
	return m
}

// ColumnStdDevs returns the population standard deviation of each column of
// data.
func ColumnStdDevs(data [][]float64) []float64 {
	means := ColumnMeans(data)
	if means == nil {
		return nil
	}
	dim := len(means)
	v := make([]float64, dim)
	for _, row := range data {
		for j := 0; j < dim && j < len(row); j++ {
			d := row[j] - means[j]
			v[j] += d * d
		}
	}
	inv := 1 / float64(len(data))
	for j := range v {
		v[j] = math.Sqrt(v[j] * inv)
	}
	return v
}

// ColumnMinMax returns the per-column minimum and maximum of data.
func ColumnMinMax(data [][]float64) (min, max []float64) {
	if len(data) == 0 {
		return nil, nil
	}
	dim := len(data[0])
	min = make([]float64, dim)
	max = make([]float64, dim)
	copy(min, data[0])
	copy(max, data[0])
	for _, row := range data {
		for j := 0; j < dim && j < len(row); j++ {
			if row[j] < min[j] {
				min[j] = row[j]
			}
			if row[j] > max[j] {
				max[j] = row[j]
			}
		}
	}
	return min, max
}

// Standardize returns a copy of data with each column transformed to zero mean
// and unit variance (z-score normalisation). Columns with zero variance are left
// centred but unscaled.
func Standardize(data [][]float64) [][]float64 {
	means := ColumnMeans(data)
	sds := ColumnStdDevs(data)
	out := make([][]float64, len(data))
	for i, row := range data {
		out[i] = make([]float64, len(row))
		for j := range row {
			if j < len(sds) && sds[j] != 0 {
				out[i][j] = (row[j] - means[j]) / sds[j]
			} else if j < len(means) {
				out[i][j] = row[j] - means[j]
			} else {
				out[i][j] = row[j]
			}
		}
	}
	return out
}

// MinMaxScale returns a copy of data with each column linearly rescaled to the
// range [0, 1]. Columns with zero range are set to zero.
func MinMaxScale(data [][]float64) [][]float64 {
	min, max := ColumnMinMax(data)
	out := make([][]float64, len(data))
	for i, row := range data {
		out[i] = make([]float64, len(row))
		for j := range row {
			if j < len(min) {
				rng := max[j] - min[j]
				if rng == 0 {
					out[i][j] = 0
				} else {
					out[i][j] = (row[j] - min[j]) / rng
				}
			} else {
				out[i][j] = row[j]
			}
		}
	}
	return out
}

// CovarianceMatrix returns the sample covariance matrix of data (dividing by
// n-1). It returns nil if there are fewer than two samples.
func CovarianceMatrix(data [][]float64) [][]float64 {
	n := len(data)
	if n < 2 {
		return nil
	}
	dim := len(data[0])
	means := ColumnMeans(data)
	cov := Zeros(dim, dim)
	for _, row := range data {
		for i := 0; i < dim; i++ {
			di := row[i] - means[i]
			for j := i; j < dim; j++ {
				cov[i][j] += di * (row[j] - means[j])
			}
		}
	}
	inv := 1 / float64(n-1)
	for i := 0; i < dim; i++ {
		for j := i; j < dim; j++ {
			cov[i][j] *= inv
			cov[j][i] = cov[i][j]
		}
	}
	return cov
}

// CovarianceMatrixMLE returns the maximum-likelihood covariance matrix of data
// (dividing by n).
func CovarianceMatrixMLE(data [][]float64) [][]float64 {
	n := len(data)
	if n < 1 {
		return nil
	}
	cov := CovarianceMatrix(data)
	if cov == nil {
		// n == 1: covariance is zero.
		dim := len(data[0])
		return Zeros(dim, dim)
	}
	scale := float64(n-1) / float64(n)
	for i := range cov {
		for j := range cov[i] {
			cov[i][j] *= scale
		}
	}
	return cov
}

// Quantile returns the q-quantile (0 <= q <= 1) of the values using linear
// interpolation. The input is copied and sorted internally.
func Quantile(values []float64, q float64) float64 {
	if len(values) == 0 {
		return math.NaN()
	}
	if q <= 0 {
		return minSlice(values)
	}
	if q >= 1 {
		return maxSlice(values)
	}
	v := append([]float64(nil), values...)
	sort.Float64s(v)
	pos := q * float64(len(v)-1)
	lo := int(math.Floor(pos))
	hi := int(math.Ceil(pos))
	if lo == hi {
		return v[lo]
	}
	frac := pos - float64(lo)
	return v[lo]*(1-frac) + v[hi]*frac
}

func minSlice(v []float64) float64 {
	m := math.Inf(1)
	for _, x := range v {
		if x < m {
			m = x
		}
	}
	return m
}

func maxSlice(v []float64) float64 {
	m := math.Inf(-1)
	for _, x := range v {
		if x > m {
			m = x
		}
	}
	return m
}

// LabelCounts returns a map from cluster label to the number of samples with
// that label.
func LabelCounts(labels []int) map[int]int {
	counts := make(map[int]int)
	for _, l := range labels {
		counts[l]++
	}
	return counts
}

// ClusterIndices groups sample indices by their cluster label, returning a map
// from label to the sorted list of sample indices assigned to it.
func ClusterIndices(labels []int) map[int][]int {
	groups := make(map[int][]int)
	for i, l := range labels {
		groups[l] = append(groups[l], i)
	}
	return groups
}

// UniqueLabels returns the sorted distinct labels present in labels.
func UniqueLabels(labels []int) []int {
	seen := make(map[int]struct{})
	for _, l := range labels {
		seen[l] = struct{}{}
	}
	out := make([]int, 0, len(seen))
	for l := range seen {
		out = append(out, l)
	}
	sort.Ints(out)
	return out
}

// CountClusters returns the number of distinct non-negative cluster labels,
// ignoring the noise label -1 used by density-based methods.
func CountClusters(labels []int) int {
	seen := make(map[int]struct{})
	for _, l := range labels {
		if l >= 0 {
			seen[l] = struct{}{}
		}
	}
	return len(seen)
}
