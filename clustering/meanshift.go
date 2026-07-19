package clustering

import (
	"math"
)

// MeanShiftResult holds the outcome of a mean-shift clustering run.
type MeanShiftResult struct {
	// Labels holds the cluster assignment of each sample.
	Labels []int
	// ClusterCenters holds the discovered modes, one per cluster.
	ClusterCenters [][]float64
	// Converged holds the shifted position of every input sample.
	Converged [][]float64
	// Iterations is the maximum number of iterations any point took to
	// converge.
	Iterations int
}

// FlatKernel is a uniform (flat) kernel: it returns 1 within the bandwidth and 0
// outside, giving each in-window point equal weight.
func FlatKernel(distance, bandwidth float64) float64 {
	if distance <= bandwidth {
		return 1
	}
	return 0
}

// GaussianKernel returns the (unnormalised) Gaussian kernel weight
// exp(-distance^2 / (2*bandwidth^2)).
func GaussianKernel(distance, bandwidth float64) float64 {
	if bandwidth <= 0 {
		return 0
	}
	x := distance / bandwidth
	return math.Exp(-0.5 * x * x)
}

// MeanShift runs the mean-shift clustering algorithm on data using a flat
// kernel of the given bandwidth and the given metric (nil means Euclidean).
// Each point is iteratively shifted toward the local mean of points within the
// bandwidth until convergence, and points whose modes are within bandwidth are
// merged into the same cluster.
func MeanShift(data [][]float64, bandwidth float64, metric Metric) (*MeanShiftResult, error) {
	return MeanShiftWithKernel(data, bandwidth, metric, FlatKernel, 300, 1e-4)
}

// MeanShiftWithKernel runs mean-shift with a caller-supplied kernel, iteration
// cap and convergence tolerance.
func MeanShiftWithKernel(data [][]float64, bandwidth float64, metric Metric, kernel func(dist, bw float64) float64, maxIter int, tol float64) (*MeanShiftResult, error) {
	if len(data) == 0 {
		return nil, ErrEmptyData
	}
	if bandwidth <= 0 {
		return nil, ErrInvalidK
	}
	if metric == nil {
		metric = Euclidean
	}
	if kernel == nil {
		kernel = FlatKernel
	}
	if maxIter <= 0 {
		maxIter = 300
	}
	if tol <= 0 {
		tol = 1e-4
	}
	n := len(data)
	dim := len(data[0])
	shifted := CloneMatrix(data)
	maxIters := 0
	for i := 0; i < n; i++ {
		point := shifted[i]
		it := 0
		for ; it < maxIter; it++ {
			mean := make([]float64, dim)
			var wsum float64
			for j := 0; j < n; j++ {
				d := metric(point, data[j])
				w := kernel(d, bandwidth)
				if w == 0 {
					continue
				}
				wsum += w
				for c := 0; c < dim; c++ {
					mean[c] += w * data[j][c]
				}
			}
			if wsum == 0 {
				break
			}
			for c := 0; c < dim; c++ {
				mean[c] /= wsum
			}
			shift := Euclidean(mean, point)
			point = mean
			if shift < tol {
				break
			}
		}
		shifted[i] = point
		if it > maxIters {
			maxIters = it
		}
	}
	// Merge modes that are within bandwidth of each other.
	centers, labels := mergeModes(shifted, bandwidth, metric)
	return &MeanShiftResult{
		Labels:         labels,
		ClusterCenters: centers,
		Converged:      shifted,
		Iterations:     maxIters,
	}, nil
}

// mergeModes groups converged points whose pairwise distance is below bandwidth,
// assigning each a cluster label and returning the representative centres.
func mergeModes(modes [][]float64, bandwidth float64, metric Metric) ([][]float64, []int) {
	var centers [][]float64
	labels := make([]int, len(modes))
	for i := range labels {
		labels[i] = -1
	}
	for i, m := range modes {
		assigned := -1
		for c, center := range centers {
			if metric(m, center) < bandwidth {
				assigned = c
				break
			}
		}
		if assigned < 0 {
			assigned = len(centers)
			centers = append(centers, CloneVector(m))
		}
		labels[i] = assigned
	}
	// Recompute centres as the mean of their members for stability.
	counts := make([]int, len(centers))
	sums := make([][]float64, len(centers))
	for c := range sums {
		sums[c] = make([]float64, len(modes[0]))
	}
	for i, l := range labels {
		counts[l]++
		for j := range sums[l] {
			sums[l][j] += modes[i][j]
		}
	}
	for c := range centers {
		if counts[c] > 0 {
			for j := range centers[c] {
				centers[c][j] = sums[c][j] / float64(counts[c])
			}
		}
	}
	return centers, labels
}

// EstimateBandwidth returns a bandwidth estimate for mean-shift based on the
// average distance to the k-th nearest neighbour, using the given quantile of
// the pairwise distances. quantile is in (0, 1]; a common choice is 0.3.
func EstimateBandwidth(data [][]float64, quantile float64, metric Metric) float64 {
	if metric == nil {
		metric = Euclidean
	}
	if quantile <= 0 || quantile > 1 {
		quantile = 0.3
	}
	n := len(data)
	if n < 2 {
		return 1
	}
	k := int(float64(n) * quantile)
	if k < 1 {
		k = 1
	}
	if k >= n {
		k = n - 1
	}
	// For each point, distance to its k-th nearest neighbour; average them.
	kd := KDistances(data, k, metric)
	var sum float64
	for _, d := range kd {
		sum += d
	}
	return sum / float64(n)
}
