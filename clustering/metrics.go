package clustering

import (
	"errors"
	"math"
	"sort"
)

// Metric is a distance function between two equal-length real vectors. A metric
// should return a nonnegative value that is zero when the two points coincide.
type Metric func(a, b []float64) float64

// checkDim panics with ErrDimensionMismatch-style behaviour avoided; metrics
// instead return NaN on mismatched lengths so they compose safely. Internally we
// use the shorter length where sensible, but exported metrics assume equal
// length as documented.

// Euclidean returns the Euclidean (L2) distance between a and b.
func Euclidean(a, b []float64) float64 {
	return math.Sqrt(SquaredEuclidean(a, b))
}

// SquaredEuclidean returns the squared Euclidean distance between a and b. It
// avoids the square root and is preferred inside inner loops where only
// relative distances matter.
func SquaredEuclidean(a, b []float64) float64 {
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	var s float64
	for i := 0; i < n; i++ {
		d := a[i] - b[i]
		s += d * d
	}
	return s
}

// Manhattan returns the Manhattan (L1, taxicab) distance between a and b.
func Manhattan(a, b []float64) float64 {
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	var s float64
	for i := 0; i < n; i++ {
		s += math.Abs(a[i] - b[i])
	}
	return s
}

// Chebyshev returns the Chebyshev (L-infinity, maximum coordinate) distance
// between a and b.
func Chebyshev(a, b []float64) float64 {
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	var max float64
	for i := 0; i < n; i++ {
		if d := math.Abs(a[i] - b[i]); d > max {
			max = d
		}
	}
	return max
}

// Minkowski returns the Minkowski distance of order p between a and b. For p=1
// it reduces to the Manhattan distance and for p=2 to the Euclidean distance.
func Minkowski(a, b []float64, p float64) float64 {
	if p <= 0 {
		return math.NaN()
	}
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	if math.IsInf(p, 1) {
		return Chebyshev(a, b)
	}
	var s float64
	for i := 0; i < n; i++ {
		s += math.Pow(math.Abs(a[i]-b[i]), p)
	}
	return math.Pow(s, 1/p)
}

// MinkowskiMetric returns a Metric computing the Minkowski distance of order p.
func MinkowskiMetric(p float64) Metric {
	return func(a, b []float64) float64 { return Minkowski(a, b, p) }
}

// Dot returns the dot (inner) product of a and b.
func Dot(a, b []float64) float64 {
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	var s float64
	for i := 0; i < n; i++ {
		s += a[i] * b[i]
	}
	return s
}

// Norm returns the Euclidean (L2) norm of the vector v.
func Norm(v []float64) float64 {
	return math.Sqrt(Dot(v, v))
}

// NormP returns the Lp norm of the vector v for p > 0.
func NormP(v []float64, p float64) float64 {
	if p <= 0 {
		return math.NaN()
	}
	if math.IsInf(p, 1) {
		var max float64
		for _, x := range v {
			if a := math.Abs(x); a > max {
				max = a
			}
		}
		return max
	}
	var s float64
	for _, x := range v {
		s += math.Pow(math.Abs(x), p)
	}
	return math.Pow(s, 1/p)
}

// CosineSimilarity returns the cosine of the angle between a and b, in the range
// [-1, 1]. If either vector has zero norm it returns 0.
func CosineSimilarity(a, b []float64) float64 {
	na := Norm(a)
	nb := Norm(b)
	if na == 0 || nb == 0 {
		return 0
	}
	return Dot(a, b) / (na * nb)
}

// Cosine returns the cosine distance 1 - CosineSimilarity(a, b), in the range
// [0, 2].
func Cosine(a, b []float64) float64 {
	return 1 - CosineSimilarity(a, b)
}

// Hamming returns the number of coordinates in which a and b differ.
func Hamming(a, b []float64) float64 {
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	var c float64
	for i := 0; i < n; i++ {
		if a[i] != b[i] {
			c++
		}
	}
	return c
}

// NormalizedHamming returns the fraction of coordinates in which a and b differ,
// in the range [0, 1].
func NormalizedHamming(a, b []float64) float64 {
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	if n == 0 {
		return 0
	}
	return Hamming(a, b) / float64(n)
}

// Canberra returns the Canberra distance between a and b, a weighted version of
// the Manhattan distance that is sensitive to small changes near zero.
func Canberra(a, b []float64) float64 {
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	var s float64
	for i := 0; i < n; i++ {
		den := math.Abs(a[i]) + math.Abs(b[i])
		if den == 0 {
			continue
		}
		s += math.Abs(a[i]-b[i]) / den
	}
	return s
}

// BrayCurtis returns the Bray-Curtis dissimilarity between a and b. It is
// typically used for nonnegative abundance data and lies in [0, 1] for such
// data.
func BrayCurtis(a, b []float64) float64 {
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	var num, den float64
	for i := 0; i < n; i++ {
		num += math.Abs(a[i] - b[i])
		den += math.Abs(a[i] + b[i])
	}
	if den == 0 {
		return 0
	}
	return num / den
}

// PearsonCorrelation returns the Pearson product-moment correlation coefficient
// between a and b, in the range [-1, 1].
func PearsonCorrelation(a, b []float64) float64 {
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	if n == 0 {
		return 0
	}
	var ma, mb float64
	for i := 0; i < n; i++ {
		ma += a[i]
		mb += b[i]
	}
	ma /= float64(n)
	mb /= float64(n)
	var num, va, vb float64
	for i := 0; i < n; i++ {
		da := a[i] - ma
		db := b[i] - mb
		num += da * db
		va += da * da
		vb += db * db
	}
	if va == 0 || vb == 0 {
		return 0
	}
	return num / math.Sqrt(va*vb)
}

// CorrelationDistance returns 1 - PearsonCorrelation(a, b), in the range [0, 2].
func CorrelationDistance(a, b []float64) float64 {
	return 1 - PearsonCorrelation(a, b)
}

// JaccardBinary returns the Jaccard distance between two binary vectors, treating
// any nonzero coordinate as set. The result is 1 minus the ratio of the size of
// the intersection to the size of the union.
func JaccardBinary(a, b []float64) float64 {
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	var inter, union float64
	for i := 0; i < n; i++ {
		ai := a[i] != 0
		bi := b[i] != 0
		if ai && bi {
			inter++
		}
		if ai || bi {
			union++
		}
	}
	if union == 0 {
		return 0
	}
	return 1 - inter/union
}

// WeightedEuclidean returns the Euclidean distance between a and b with each
// coordinate weighted by the corresponding entry of w.
func WeightedEuclidean(a, b, w []float64) float64 {
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	if len(w) < n {
		n = len(w)
	}
	var s float64
	for i := 0; i < n; i++ {
		d := a[i] - b[i]
		s += w[i] * d * d
	}
	return math.Sqrt(s)
}

// WeightedEuclideanMetric returns a Metric computing the weighted Euclidean
// distance with the given per-coordinate weights.
func WeightedEuclideanMetric(w []float64) Metric {
	ww := CloneVector(w)
	return func(a, b []float64) float64 { return WeightedEuclidean(a, b, ww) }
}

// Mahalanobis returns the Mahalanobis distance between a and b given the inverse
// covariance matrix invCov, defined as sqrt((a-b)^T invCov (a-b)).
func Mahalanobis(a, b []float64, invCov [][]float64) float64 {
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	diff := make([]float64, n)
	for i := 0; i < n; i++ {
		diff[i] = a[i] - b[i]
	}
	var s float64
	for i := 0; i < n; i++ {
		var row float64
		for j := 0; j < n; j++ {
			row += invCov[i][j] * diff[j]
		}
		s += diff[i] * row
	}
	if s < 0 {
		s = 0
	}
	return math.Sqrt(s)
}

// MahalanobisMetric returns a Metric computing the Mahalanobis distance using
// the given inverse covariance matrix.
func MahalanobisMetric(invCov [][]float64) Metric {
	c := CloneMatrix(invCov)
	return func(a, b []float64) float64 { return Mahalanobis(a, b, c) }
}

// NewMahalanobisMetric builds a Mahalanobis Metric from a covariance matrix by
// inverting it. It returns an error if the covariance matrix is singular.
func NewMahalanobisMetric(cov [][]float64) (Metric, error) {
	inv, err := Inverse(cov)
	if err != nil {
		return nil, err
	}
	return MahalanobisMetric(inv), nil
}

// DistanceMatrix returns the full symmetric n x n matrix of pairwise distances
// between the rows of data using the given metric. If metric is nil, Euclidean
// is used.
func DistanceMatrix(data [][]float64, metric Metric) [][]float64 {
	if metric == nil {
		metric = Euclidean
	}
	n := len(data)
	d := make([][]float64, n)
	for i := range d {
		d[i] = make([]float64, n)
	}
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			dist := metric(data[i], data[j])
			d[i][j] = dist
			d[j][i] = dist
		}
	}
	return d
}

// CondensedDistanceMatrix returns the pairwise distances between rows of data in
// condensed upper-triangular form: entries for pairs (i, j) with i < j in
// row-major order. The length of the result is n*(n-1)/2.
func CondensedDistanceMatrix(data [][]float64, metric Metric) []float64 {
	if metric == nil {
		metric = Euclidean
	}
	n := len(data)
	out := make([]float64, 0, n*(n-1)/2)
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			out = append(out, metric(data[i], data[j]))
		}
	}
	return out
}

// CondensedIndex returns the position in a condensed distance matrix (as
// produced by CondensedDistanceMatrix) of the pair (i, j) for n points. It
// returns -1 if i == j or an index is out of range.
func CondensedIndex(n, i, j int) int {
	if i == j || i < 0 || j < 0 || i >= n || j >= n {
		return -1
	}
	if i > j {
		i, j = j, i
	}
	return i*n - (i*(i+1))/2 + (j - i - 1)
}

// NearestNeighbor returns the index of the row of data closest to point under
// the given metric, together with that distance. It returns (-1, +Inf) for empty
// data.
func NearestNeighbor(point []float64, data [][]float64, metric Metric) (int, float64) {
	if metric == nil {
		metric = Euclidean
	}
	best := -1
	bestDist := math.Inf(1)
	for i, row := range data {
		if d := metric(point, row); d < bestDist {
			bestDist = d
			best = i
		}
	}
	return best, bestDist
}

// KNearestNeighbors returns the indices of the k rows of data closest to point
// under the given metric, ordered from nearest to farthest.
func KNearestNeighbors(point []float64, data [][]float64, k int, metric Metric) []int {
	if metric == nil {
		metric = Euclidean
	}
	all := make([]neighborDist, len(data))
	for i, row := range data {
		all[i] = neighborDist{i, metric(point, row)}
	}
	sort.SliceStable(all, func(a, b int) bool { return all[a].d < all[b].d })
	if k > len(all) {
		k = len(all)
	}
	out := make([]int, k)
	for i := 0; i < k; i++ {
		out[i] = all[i].idx
	}
	return out
}

type neighborDist struct {
	idx int
	d   float64
}

// MetricByName returns a Metric for a named distance function. Recognised names
// (case-insensitive) are "euclidean", "sqeuclidean", "manhattan"/"cityblock"/"l1",
// "chebyshev"/"linf", "cosine", "correlation", "canberra", "braycurtis",
// "hamming" and "jaccard". It returns an error for unknown names.
func MetricByName(name string) (Metric, error) {
	switch normalizeName(name) {
	case "euclidean", "l2":
		return Euclidean, nil
	case "sqeuclidean", "squaredeuclidean":
		return SquaredEuclidean, nil
	case "manhattan", "cityblock", "l1", "taxicab":
		return Manhattan, nil
	case "chebyshev", "linf", "chebychev", "maximum":
		return Chebyshev, nil
	case "cosine":
		return Cosine, nil
	case "correlation":
		return CorrelationDistance, nil
	case "canberra":
		return Canberra, nil
	case "braycurtis":
		return BrayCurtis, nil
	case "hamming":
		return Hamming, nil
	case "jaccard":
		return JaccardBinary, nil
	default:
		return nil, errors.New("clustering: unknown metric " + name)
	}
}

func normalizeName(s string) string {
	out := make([]rune, 0, len(s))
	for _, r := range s {
		switch {
		case r >= 'A' && r <= 'Z':
			out = append(out, r+('a'-'A'))
		case r == ' ' || r == '_' || r == '-':
			// skip separators
		default:
			out = append(out, r)
		}
	}
	return string(out)
}
