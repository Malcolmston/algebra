package clustering

import (
	"math"
	"math/rand"
)

// KMedoidsResult holds the outcome of a k-medoids (PAM) clustering run.
type KMedoidsResult struct {
	// K is the number of clusters.
	K int
	// MedoidIndices holds the indices into the original data of the chosen
	// medoids.
	MedoidIndices []int
	// Medoids holds copies of the medoid points themselves.
	Medoids [][]float64
	// Labels holds the cluster assignment (0..K-1) of every sample.
	Labels []int
	// Cost is the total distance of samples to their assigned medoid.
	Cost float64
	// Iterations is the number of swap-improvement iterations performed.
	Iterations int
}

// KMedoids runs the Partitioning Around Medoids (PAM) algorithm on data using
// the given metric (nil means Euclidean). It uses the BUILD phase for
// initialisation followed by the SWAP phase, and returns the resulting medoids
// and labels.
func KMedoids(data [][]float64, k int, metric Metric, rng *rand.Rand) (*KMedoidsResult, error) {
	if len(data) == 0 {
		return nil, ErrEmptyData
	}
	if k <= 0 || k > len(data) {
		return nil, ErrInvalidK
	}
	if metric == nil {
		metric = Euclidean
	}
	if rng == nil {
		rng = rand.New(rand.NewSource(1))
	}
	dm := DistanceMatrix(data, metric)
	medoids := pamBuild(dm, k)
	medoids, iters := pamSwap(dm, medoids)
	labels, cost := assignToMedoids(dm, medoids)
	medPoints := make([][]float64, k)
	for i, m := range medoids {
		medPoints[i] = CloneVector(data[m])
	}
	return &KMedoidsResult{
		K:             k,
		MedoidIndices: medoids,
		Medoids:       medPoints,
		Labels:        labels,
		Cost:          cost,
		Iterations:    iters,
	}, nil
}

// KMedoidsFromDistances runs PAM directly on a precomputed symmetric distance
// matrix, avoiding recomputation of pairwise distances.
func KMedoidsFromDistances(dm [][]float64, k int) (*KMedoidsResult, error) {
	n := len(dm)
	if n == 0 {
		return nil, ErrEmptyData
	}
	if k <= 0 || k > n {
		return nil, ErrInvalidK
	}
	medoids := pamBuild(dm, k)
	medoids, iters := pamSwap(dm, medoids)
	labels, cost := assignToMedoids(dm, medoids)
	return &KMedoidsResult{
		K:             k,
		MedoidIndices: medoids,
		Labels:        labels,
		Cost:          cost,
		Iterations:    iters,
	}, nil
}

// pamBuild selects k initial medoids greedily so as to minimise total distance.
func pamBuild(dm [][]float64, k int) []int {
	n := len(dm)
	chosen := make([]bool, n)
	medoids := make([]int, 0, k)

	// First medoid: the point minimising total distance to all others.
	best := 0
	bestSum := math.Inf(1)
	for i := 0; i < n; i++ {
		var s float64
		for j := 0; j < n; j++ {
			s += dm[i][j]
		}
		if s < bestSum {
			bestSum = s
			best = i
		}
	}
	medoids = append(medoids, best)
	chosen[best] = true

	// nearest[j] = distance from j to the closest chosen medoid so far.
	nearest := make([]float64, n)
	for j := 0; j < n; j++ {
		nearest[j] = dm[best][j]
	}

	for len(medoids) < k {
		bestCand := -1
		bestGain := -math.Inf(1)
		for c := 0; c < n; c++ {
			if chosen[c] {
				continue
			}
			var gain float64
			for j := 0; j < n; j++ {
				if d := nearest[j] - dm[c][j]; d > 0 {
					gain += d
				}
			}
			if gain > bestGain {
				bestGain = gain
				bestCand = c
			}
		}
		if bestCand < 0 {
			break
		}
		medoids = append(medoids, bestCand)
		chosen[bestCand] = true
		for j := 0; j < n; j++ {
			if dm[bestCand][j] < nearest[j] {
				nearest[j] = dm[bestCand][j]
			}
		}
	}
	return medoids
}

// pamSwap performs the SWAP phase, repeatedly swapping a medoid with a
// non-medoid whenever this reduces the total cost.
func pamSwap(dm [][]float64, medoids []int) ([]int, int) {
	n := len(dm)
	current := append([]int(nil), medoids...)
	isMedoid := make([]bool, n)
	for _, m := range current {
		isMedoid[m] = true
	}
	_, curCost := assignToMedoids(dm, current)
	iters := 0
	for {
		bestDelta := 0.0
		bestI, bestH := -1, -1
		for mi := range current {
			for h := 0; h < n; h++ {
				if isMedoid[h] {
					continue
				}
				trial := append([]int(nil), current...)
				trial[mi] = h
				_, cost := assignToMedoids(dm, trial)
				if delta := cost - curCost; delta < bestDelta-1e-12 {
					bestDelta = delta
					bestI = mi
					bestH = h
				}
			}
		}
		if bestI < 0 {
			break
		}
		isMedoid[current[bestI]] = false
		current[bestI] = bestH
		isMedoid[bestH] = true
		curCost += bestDelta
		iters++
	}
	return current, iters
}

// assignToMedoids assigns each point to its nearest medoid and returns the
// labels (0..k-1 indexing into medoids) and the total cost.
func assignToMedoids(dm [][]float64, medoids []int) ([]int, float64) {
	n := len(dm)
	labels := make([]int, n)
	var cost float64
	for i := 0; i < n; i++ {
		best := 0
		bestDist := math.Inf(1)
		for mi, m := range medoids {
			if d := dm[i][m]; d < bestDist {
				bestDist = d
				best = mi
			}
		}
		labels[i] = best
		cost += bestDist
	}
	return labels, cost
}

// Predict assigns each row of newData to the nearest medoid of a fitted result
// using the given metric.
func (r *KMedoidsResult) Predict(newData [][]float64, metric Metric) []int {
	if metric == nil {
		metric = Euclidean
	}
	labels := make([]int, len(newData))
	for i, row := range newData {
		best := 0
		bestDist := math.Inf(1)
		for mi, med := range r.Medoids {
			if d := metric(row, med); d < bestDist {
				bestDist = d
				best = mi
			}
		}
		labels[i] = best
	}
	return labels
}
