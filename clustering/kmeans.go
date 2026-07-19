package clustering

import (
	"math"
	"math/rand"
)

// KMeansResult holds the outcome of a k-means clustering run.
type KMeansResult struct {
	// K is the number of clusters.
	K int
	// Centroids holds the K cluster centres, one per row.
	Centroids [][]float64
	// Labels holds the cluster assignment for each input sample.
	Labels []int
	// Inertia is the sum of squared distances of samples to their assigned
	// centroid.
	Inertia float64
	// Iterations is the number of Lloyd iterations performed.
	Iterations int
	// Converged reports whether the run stopped because assignments stopped
	// changing (rather than hitting the iteration cap).
	Converged bool
}

// KMeansOptions configures a k-means run.
type KMeansOptions struct {
	// MaxIter is the maximum number of Lloyd iterations per run (default 300).
	MaxIter int
	// Tol is the centroid-shift tolerance below which the run is considered
	// converged (default 1e-4).
	Tol float64
	// NInit is the number of independent restarts; the best (lowest-inertia)
	// result is kept (default 10).
	NInit int
	// Init selects the initialisation strategy: "kmeans++" (default) or
	// "random".
	Init string
	// Rand is the random source. If nil, a deterministic source seeded with 1 is
	// used so results are reproducible.
	Rand *rand.Rand
}

func (o KMeansOptions) withDefaults() KMeansOptions {
	if o.MaxIter <= 0 {
		o.MaxIter = 300
	}
	if o.Tol <= 0 {
		o.Tol = 1e-4
	}
	if o.NInit <= 0 {
		o.NInit = 10
	}
	if o.Init == "" {
		o.Init = "kmeans++"
	}
	if o.Rand == nil {
		o.Rand = rand.New(rand.NewSource(1))
	}
	return o
}

// KMeans runs k-means clustering on data with the default options, returning the
// best result over the default number of restarts.
func KMeans(data [][]float64, k int) (*KMeansResult, error) {
	return KMeansWithOptions(data, k, KMeansOptions{})
}

// KMeansWithOptions runs k-means clustering on data using the supplied options.
// It performs opts.NInit restarts and returns the result with the lowest
// inertia.
func KMeansWithOptions(data [][]float64, k int, opts KMeansOptions) (*KMeansResult, error) {
	if len(data) == 0 {
		return nil, ErrEmptyData
	}
	if k <= 0 || k > len(data) {
		return nil, ErrInvalidK
	}
	opts = opts.withDefaults()
	var best *KMeansResult
	for run := 0; run < opts.NInit; run++ {
		var init [][]float64
		switch normalizeName(opts.Init) {
		case "random":
			init = RandomInit(data, k, opts.Rand)
		case "forgy":
			init = ForgyInit(data, k, opts.Rand)
		default:
			init = KMeansPlusPlusInit(data, k, opts.Rand, Euclidean)
		}
		res := lloyd(data, init, opts.MaxIter, opts.Tol)
		if best == nil || res.Inertia < best.Inertia {
			best = res
		}
	}
	return best, nil
}

// lloyd runs Lloyd's algorithm from the given initial centroids.
func lloyd(data, initCentroids [][]float64, maxIter int, tol float64) *KMeansResult {
	k := len(initCentroids)
	centroids := CloneMatrix(initCentroids)
	labels := make([]int, len(data))
	converged := false
	iter := 0
	for iter = 0; iter < maxIter; iter++ {
		changed := AssignClustersInto(data, centroids, SquaredEuclidean, labels)
		newCentroids := UpdateCentroids(data, labels, k, centroids)
		shift := maxCentroidShift(centroids, newCentroids)
		centroids = newCentroids
		if !changed || shift <= tol {
			converged = true
			iter++
			break
		}
	}
	// Final assignment against the last centroids.
	AssignClustersInto(data, centroids, SquaredEuclidean, labels)
	inertia := ComputeInertia(data, centroids, labels, SquaredEuclidean)
	return &KMeansResult{
		K:          k,
		Centroids:  centroids,
		Labels:     labels,
		Inertia:    inertia,
		Iterations: iter,
		Converged:  converged,
	}
}

// KMeansPlusPlusInit selects k initial centroids from data using the k-means++
// seeding strategy, which spreads the initial centres according to squared
// distance. metric may be nil (Euclidean).
func KMeansPlusPlusInit(data [][]float64, k int, rng *rand.Rand, metric Metric) [][]float64 {
	if metric == nil {
		metric = Euclidean
	}
	if rng == nil {
		rng = rand.New(rand.NewSource(1))
	}
	n := len(data)
	centroids := make([][]float64, 0, k)
	first := rng.Intn(n)
	centroids = append(centroids, CloneVector(data[first]))
	dist2 := make([]float64, n)
	for i := range dist2 {
		d := metric(data[i], centroids[0])
		dist2[i] = d * d
	}
	for len(centroids) < k {
		var total float64
		for _, d := range dist2 {
			total += d
		}
		var chosen int
		if total == 0 {
			chosen = rng.Intn(n)
		} else {
			target := rng.Float64() * total
			var cum float64
			chosen = n - 1
			for i, d := range dist2 {
				cum += d
				if cum >= target {
					chosen = i
					break
				}
			}
		}
		centroids = append(centroids, CloneVector(data[chosen]))
		newC := centroids[len(centroids)-1]
		for i := range dist2 {
			d := metric(data[i], newC)
			if d2 := d * d; d2 < dist2[i] {
				dist2[i] = d2
			}
		}
	}
	return centroids
}

// RandomInit selects k distinct data points uniformly at random as initial
// centroids.
func RandomInit(data [][]float64, k int, rng *rand.Rand) [][]float64 {
	if rng == nil {
		rng = rand.New(rand.NewSource(1))
	}
	perm := rng.Perm(len(data))
	centroids := make([][]float64, k)
	for i := 0; i < k; i++ {
		centroids[i] = CloneVector(data[perm[i]])
	}
	return centroids
}

// ForgyInit is an alias for RandomInit: it picks k random data points as the
// initial centroids (the Forgy method).
func ForgyInit(data [][]float64, k int, rng *rand.Rand) [][]float64 {
	return RandomInit(data, k, rng)
}

// AssignClusters returns, for each sample, the index of the nearest centroid
// under the given metric. metric may be nil (Euclidean).
func AssignClusters(data, centroids [][]float64, metric Metric) []int {
	labels := make([]int, len(data))
	AssignClustersInto(data, centroids, metric, labels)
	return labels
}

// AssignClustersInto assigns each sample to its nearest centroid, writing the
// result into labels (which must have length len(data)). It returns true if any
// assignment changed from the previous contents of labels.
func AssignClustersInto(data, centroids [][]float64, metric Metric, labels []int) bool {
	if metric == nil {
		metric = SquaredEuclidean
	}
	changed := false
	for i, row := range data {
		best := 0
		bestDist := math.Inf(1)
		for c, cen := range centroids {
			if d := metric(row, cen); d < bestDist {
				bestDist = d
				best = c
			}
		}
		if labels[i] != best {
			changed = true
			labels[i] = best
		}
	}
	return changed
}

// UpdateCentroids recomputes each centroid as the mean of the samples assigned
// to it. Empty clusters retain their previous centroid taken from prev (which
// may be nil, in which case an empty cluster keeps a zero centroid).
func UpdateCentroids(data [][]float64, labels []int, k int, prev [][]float64) [][]float64 {
	dim := len(data[0])
	centroids := Zeros(k, dim)
	counts := make([]int, k)
	for i, row := range data {
		c := labels[i]
		counts[c]++
		for j := 0; j < dim && j < len(row); j++ {
			centroids[c][j] += row[j]
		}
	}
	for c := 0; c < k; c++ {
		if counts[c] == 0 {
			if prev != nil && c < len(prev) {
				copy(centroids[c], prev[c])
			}
			continue
		}
		inv := 1 / float64(counts[c])
		for j := 0; j < dim; j++ {
			centroids[c][j] *= inv
		}
	}
	return centroids
}

// ComputeInertia returns the sum over all samples of the distance from each
// sample to its assigned centroid, using the given metric. For the standard
// k-means objective pass SquaredEuclidean.
func ComputeInertia(data, centroids [][]float64, labels []int, metric Metric) float64 {
	if metric == nil {
		metric = SquaredEuclidean
	}
	var s float64
	for i, row := range data {
		s += metric(row, centroids[labels[i]])
	}
	return s
}

// KMeansPredict assigns each row of newData to the nearest centroid of a fitted
// result, returning the cluster labels.
func (r *KMeansResult) Predict(newData [][]float64) []int {
	return AssignClusters(newData, r.Centroids, SquaredEuclidean)
}

// Transform returns, for each row of newData, the Euclidean distance to every
// centroid of the fitted result.
func (r *KMeansResult) Transform(newData [][]float64) [][]float64 {
	out := make([][]float64, len(newData))
	for i, row := range newData {
		out[i] = make([]float64, r.K)
		for c, cen := range r.Centroids {
			out[i][c] = Euclidean(row, cen)
		}
	}
	return out
}

func maxCentroidShift(a, b [][]float64) float64 {
	var max float64
	for i := range a {
		if d := SquaredEuclidean(a[i], b[i]); d > max {
			max = d
		}
	}
	return math.Sqrt(max)
}

// MiniBatchKMeans runs a mini-batch variant of k-means, updating centroids from
// random subsets of size batchSize for the given number of iterations. It is
// faster on large datasets at the cost of some accuracy. It returns a result
// with final full-data labels and inertia.
func MiniBatchKMeans(data [][]float64, k, batchSize, iterations int, rng *rand.Rand) (*KMeansResult, error) {
	if len(data) == 0 {
		return nil, ErrEmptyData
	}
	if k <= 0 || k > len(data) {
		return nil, ErrInvalidK
	}
	if rng == nil {
		rng = rand.New(rand.NewSource(1))
	}
	if batchSize <= 0 {
		batchSize = 100
	}
	if batchSize > len(data) {
		batchSize = len(data)
	}
	if iterations <= 0 {
		iterations = 100
	}
	centroids := KMeansPlusPlusInit(data, k, rng, Euclidean)
	counts := make([]int, k)
	for it := 0; it < iterations; it++ {
		for b := 0; b < batchSize; b++ {
			idx := rng.Intn(len(data))
			row := data[idx]
			// nearest centroid
			best := 0
			bestDist := math.Inf(1)
			for c, cen := range centroids {
				if d := SquaredEuclidean(row, cen); d < bestDist {
					bestDist = d
					best = c
				}
			}
			counts[best]++
			eta := 1 / float64(counts[best])
			for j := range centroids[best] {
				centroids[best][j] = (1-eta)*centroids[best][j] + eta*row[j]
			}
		}
	}
	labels := AssignClusters(data, centroids, SquaredEuclidean)
	inertia := ComputeInertia(data, centroids, labels, SquaredEuclidean)
	return &KMeansResult{
		K:          k,
		Centroids:  centroids,
		Labels:     labels,
		Inertia:    inertia,
		Iterations: iterations,
		Converged:  true,
	}, nil
}
