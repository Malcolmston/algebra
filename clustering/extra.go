package clustering

import (
	"math"
	"math/rand"
)

// BisectingKMeans clusters data into k clusters top-down: it starts with all
// points in one cluster and repeatedly splits the cluster with the largest
// within-cluster sum of squares into two using 2-means, until k clusters exist.
func BisectingKMeans(data [][]float64, k int, rng *rand.Rand) (*KMeansResult, error) {
	if len(data) == 0 {
		return nil, ErrEmptyData
	}
	if k <= 0 || k > len(data) {
		return nil, ErrInvalidK
	}
	if rng == nil {
		rng = rand.New(rand.NewSource(1))
	}
	// Each cluster is a slice of original indices.
	clusters := [][]int{make([]int, len(data))}
	for i := range data {
		clusters[0][i] = i
	}
	for len(clusters) < k {
		// Pick the cluster with the largest SSE (and at least 2 points).
		bestC := -1
		bestSSE := -1.0
		for ci, members := range clusters {
			if len(members) < 2 {
				continue
			}
			pts := subset(data, members)
			sse := Inertia(pts, zeroLabels(len(pts)))
			if sse > bestSSE {
				bestSSE = sse
				bestC = ci
			}
		}
		if bestC < 0 {
			break
		}
		members := clusters[bestC]
		pts := subset(data, members)
		res, err := KMeansWithOptions(pts, 2, KMeansOptions{NInit: 5, Rand: rng})
		if err != nil {
			return nil, err
		}
		var left, right []int
		for i, lbl := range res.Labels {
			if lbl == 0 {
				left = append(left, members[i])
			} else {
				right = append(right, members[i])
			}
		}
		if len(left) == 0 || len(right) == 0 {
			// Degenerate split; leave the cluster intact to avoid an infinite
			// loop.
			break
		}
		clusters[bestC] = left
		clusters = append(clusters, right)
	}
	labels := make([]int, len(data))
	centroids := make([][]float64, len(clusters))
	for ci, members := range clusters {
		for _, m := range members {
			labels[m] = ci
		}
		centroids[ci] = Centroid(subset(data, members))
	}
	inertia := ComputeInertia(data, centroids, labels, SquaredEuclidean)
	return &KMeansResult{
		K:         len(clusters),
		Centroids: centroids,
		Labels:    labels,
		Inertia:   inertia,
		Converged: true,
	}, nil
}

func subset(data [][]float64, idx []int) [][]float64 {
	out := make([][]float64, len(idx))
	for i, j := range idx {
		out[i] = data[j]
	}
	return out
}

func zeroLabels(n int) []int {
	return make([]int, n)
}

// FuzzyCMeansResult holds the outcome of fuzzy c-means clustering.
type FuzzyCMeansResult struct {
	// Centers holds the c cluster centres.
	Centers [][]float64
	// Membership[i][j] is the degree (in [0,1]) to which sample i belongs to
	// cluster j; each row sums to 1.
	Membership [][]float64
	// Labels holds the hard assignment (arg-max membership) for each sample.
	Labels []int
	// Objective is the final value of the fuzzy c-means objective function.
	Objective float64
	// Iterations is the number of iterations performed.
	Iterations int
}

// FuzzyCMeans runs the fuzzy c-means (soft k-means) algorithm on data with c
// clusters and fuzziness exponent m > 1 (a common value is 2). Each sample is
// given a graded membership to every cluster.
func FuzzyCMeans(data [][]float64, c int, m float64, maxIter int, tol float64, rng *rand.Rand) (*FuzzyCMeansResult, error) {
	if len(data) == 0 {
		return nil, ErrEmptyData
	}
	if c <= 0 || c > len(data) {
		return nil, ErrInvalidK
	}
	if m <= 1 {
		m = 2
	}
	if maxIter <= 0 {
		maxIter = 300
	}
	if tol <= 0 {
		tol = 1e-5
	}
	if rng == nil {
		rng = rand.New(rand.NewSource(1))
	}
	n := len(data)
	dim := len(data[0])
	// Initialise membership randomly, normalised per row.
	u := make([][]float64, n)
	for i := 0; i < n; i++ {
		u[i] = make([]float64, c)
		var s float64
		for j := 0; j < c; j++ {
			u[i][j] = rng.Float64() + 1e-6
			s += u[i][j]
		}
		for j := 0; j < c; j++ {
			u[i][j] /= s
		}
	}
	centers := Zeros(c, dim)
	exp := 2 / (m - 1)
	iter := 0
	prevObj := math.Inf(1)
	for iter = 0; iter < maxIter; iter++ {
		// Update centres.
		for j := 0; j < c; j++ {
			num := make([]float64, dim)
			var den float64
			for i := 0; i < n; i++ {
				w := math.Pow(u[i][j], m)
				den += w
				for d := 0; d < dim; d++ {
					num[d] += w * data[i][d]
				}
			}
			if den == 0 {
				den = 1e-300
			}
			for d := 0; d < dim; d++ {
				centers[j][d] = num[d] / den
			}
		}
		// Update membership.
		for i := 0; i < n; i++ {
			dists := make([]float64, c)
			zero := -1
			for j := 0; j < c; j++ {
				dists[j] = Euclidean(data[i], centers[j])
				if dists[j] == 0 {
					zero = j
				}
			}
			if zero >= 0 {
				for j := 0; j < c; j++ {
					if j == zero {
						u[i][j] = 1
					} else {
						u[i][j] = 0
					}
				}
				continue
			}
			for j := 0; j < c; j++ {
				var sum float64
				for l := 0; l < c; l++ {
					sum += math.Pow(dists[j]/dists[l], exp)
				}
				u[i][j] = 1 / sum
			}
		}
		// Objective.
		var obj float64
		for i := 0; i < n; i++ {
			for j := 0; j < c; j++ {
				obj += math.Pow(u[i][j], m) * SquaredEuclidean(data[i], centers[j])
			}
		}
		if math.Abs(prevObj-obj) < tol {
			prevObj = obj
			iter++
			break
		}
		prevObj = obj
	}
	labels := make([]int, n)
	for i := 0; i < n; i++ {
		best := 0
		for j := 1; j < c; j++ {
			if u[i][j] > u[i][best] {
				best = j
			}
		}
		labels[i] = best
	}
	return &FuzzyCMeansResult{
		Centers:    centers,
		Membership: u,
		Labels:     labels,
		Objective:  prevObj,
		Iterations: iter,
	}, nil
}

// PartitionCoefficient returns the fuzzy partition coefficient of a membership
// matrix, in the range [1/c, 1]; values near 1 indicate crisp (well separated)
// clustering.
func PartitionCoefficient(membership [][]float64) float64 {
	n := len(membership)
	if n == 0 {
		return 0
	}
	var s float64
	for i := range membership {
		for j := range membership[i] {
			s += membership[i][j] * membership[i][j]
		}
	}
	return s / float64(n)
}

// PartitionEntropy returns the fuzzy partition entropy of a membership matrix;
// lower values indicate crisper clustering.
func PartitionEntropy(membership [][]float64) float64 {
	n := len(membership)
	if n == 0 {
		return 0
	}
	var s float64
	for i := range membership {
		for j := range membership[i] {
			u := membership[i][j]
			if u > 0 {
				s -= u * math.Log(u)
			}
		}
	}
	return s / float64(n)
}
