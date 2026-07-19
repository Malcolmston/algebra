package clustering

import (
	"math"
	"sort"
)

// NoiseLabel is the label assigned to points that DBSCAN classifies as noise
// (not belonging to any cluster).
const NoiseLabel = -1

// DBSCANResult holds the outcome of a DBSCAN run.
type DBSCANResult struct {
	// Labels holds the cluster label of each sample; noise points have label
	// NoiseLabel (-1). Cluster labels are 0-based.
	Labels []int
	// Core marks whether each sample is a core point.
	Core []bool
	// NClusters is the number of clusters found (excluding noise).
	NClusters int
}

// DBSCAN runs density-based spatial clustering (DBSCAN) on data. A point is a
// core point if at least minPts points (including itself) lie within distance
// eps under the given metric (nil means Euclidean). Clusters are grown from core
// points; points not reachable from any core point are labelled as noise
// (NoiseLabel).
func DBSCAN(data [][]float64, eps float64, minPts int, metric Metric) (*DBSCANResult, error) {
	if len(data) == 0 {
		return nil, ErrEmptyData
	}
	if metric == nil {
		metric = Euclidean
	}
	if minPts < 1 {
		minPts = 1
	}
	n := len(data)
	neighbors := RegionQueryAll(data, eps, metric)
	labels := make([]int, n)
	for i := range labels {
		labels[i] = -2 // undefined
	}
	core := make([]bool, n)
	for i := 0; i < n; i++ {
		core[i] = len(neighbors[i]) >= minPts
	}
	clusterID := 0
	for i := 0; i < n; i++ {
		if labels[i] != -2 {
			continue
		}
		if !core[i] {
			labels[i] = NoiseLabel
			continue
		}
		// Start a new cluster and expand via BFS over the seed set.
		labels[i] = clusterID
		queue := append([]int(nil), neighbors[i]...)
		for qi := 0; qi < len(queue); qi++ {
			p := queue[qi]
			if labels[p] == NoiseLabel {
				labels[p] = clusterID // border point
			}
			if labels[p] != -2 {
				continue
			}
			labels[p] = clusterID
			if core[p] {
				queue = append(queue, neighbors[p]...)
			}
		}
		clusterID++
	}
	return &DBSCANResult{Labels: labels, Core: core, NClusters: clusterID}, nil
}

// RegionQueryAll returns, for each sample, the indices of all samples within
// distance eps (including the point itself), using the given metric.
func RegionQueryAll(data [][]float64, eps float64, metric Metric) [][]int {
	if metric == nil {
		metric = Euclidean
	}
	n := len(data)
	neighbors := make([][]int, n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if metric(data[i], data[j]) <= eps {
				neighbors[i] = append(neighbors[i], j)
			}
		}
	}
	return neighbors
}

// RegionQuery returns the indices of samples within distance eps of the query
// point (including any coincident points), using the given metric.
func RegionQuery(data [][]float64, point []float64, eps float64, metric Metric) []int {
	if metric == nil {
		metric = Euclidean
	}
	var out []int
	for j, row := range data {
		if metric(point, row) <= eps {
			out = append(out, j)
		}
	}
	return out
}

// KDistances returns, for each sample, the distance to its k-th nearest
// neighbour (excluding the point itself). Sorting these values produces the
// classic "k-distance graph" used to choose eps for DBSCAN.
func KDistances(data [][]float64, k int, metric Metric) []float64 {
	if metric == nil {
		metric = Euclidean
	}
	n := len(data)
	out := make([]float64, n)
	for i := 0; i < n; i++ {
		dists := make([]float64, 0, n-1)
		for j := 0; j < n; j++ {
			if i != j {
				dists = append(dists, metric(data[i], data[j]))
			}
		}
		sort.Float64s(dists)
		if k-1 < len(dists) {
			out[i] = dists[k-1]
		} else if len(dists) > 0 {
			out[i] = dists[len(dists)-1]
		}
	}
	return out
}

// SortedKDistances returns the k-distances of all points sorted in ascending
// order, ready to plot as a k-distance graph for eps selection.
func SortedKDistances(data [][]float64, k int, metric Metric) []float64 {
	d := KDistances(data, k, metric)
	sort.Float64s(d)
	return d
}

// OPTICSResult holds the ordering and reachability produced by the OPTICS
// algorithm.
type OPTICSResult struct {
	// Order is the cluster ordering of point indices.
	Order []int
	// Reachability holds the reachability distance of each point in the order
	// produced; the first point's reachability is +Inf.
	Reachability []float64
	// CoreDistance holds the core distance of each point (indexed by original
	// sample index), or +Inf if the point is not a core point.
	CoreDistance []float64
}

// OPTICS computes the OPTICS ordering of data, generalising DBSCAN across a
// range of eps values. maxEps bounds the neighbourhood radius (use math.Inf(1)
// for unbounded) and minPts sets the density threshold. Metric may be nil.
func OPTICS(data [][]float64, maxEps float64, minPts int, metric Metric) (*OPTICSResult, error) {
	if len(data) == 0 {
		return nil, ErrEmptyData
	}
	if metric == nil {
		metric = Euclidean
	}
	if minPts < 1 {
		minPts = 1
	}
	n := len(data)
	processed := make([]bool, n)
	reachOut := make([]float64, 0, n)
	order := make([]int, 0, n)
	reachability := make([]float64, n)
	coreDist := make([]float64, n)
	for i := range reachability {
		reachability[i] = math.Inf(1)
	}

	coreDistance := func(i int) float64 {
		dists := make([]float64, 0, n)
		for j := 0; j < n; j++ {
			d := metric(data[i], data[j])
			if d <= maxEps {
				dists = append(dists, d)
			}
		}
		if len(dists) < minPts {
			return math.Inf(1)
		}
		sort.Float64s(dists)
		return dists[minPts-1]
	}
	for i := 0; i < n; i++ {
		coreDist[i] = coreDistance(i)
	}

	for i := 0; i < n; i++ {
		if processed[i] {
			continue
		}
		// Process point i and expand its seed list.
		seeds := map[int]float64{}
		process := func(p int) {
			processed[p] = true
			order = append(order, p)
			reachOut = append(reachOut, reachability[p])
			cd := coreDist[p]
			if math.IsInf(cd, 1) {
				return
			}
			for j := 0; j < n; j++ {
				if processed[j] {
					continue
				}
				d := metric(data[p], data[j])
				if d > maxEps {
					continue
				}
				newReach := math.Max(cd, d)
				if newReach < reachability[j] {
					reachability[j] = newReach
					seeds[j] = newReach
				}
			}
		}
		process(i)
		for len(seeds) > 0 {
			// Pop the seed with the smallest reachability.
			bestP := -1
			bestR := math.Inf(1)
			for p, r := range seeds {
				if r < bestR || (r == bestR && p < bestP) {
					bestR = r
					bestP = p
				}
			}
			delete(seeds, bestP)
			if processed[bestP] {
				continue
			}
			process(bestP)
		}
	}
	return &OPTICSResult{Order: order, Reachability: reachOut, CoreDistance: coreDist}, nil
}

// ExtractDBSCAN extracts a flat DBSCAN-style clustering from an OPTICS result at
// the given eps threshold, returning a label per original sample (noise =
// NoiseLabel).
func (o *OPTICSResult) ExtractDBSCAN(eps float64) []int {
	n := len(o.Order)
	labels := make([]int, n)
	for i := range labels {
		labels[i] = NoiseLabel
	}
	clusterID := -1
	for k, p := range o.Order {
		r := o.Reachability[k]
		if r > eps {
			if o.CoreDistance[p] <= eps {
				clusterID++
				labels[p] = clusterID
			} else {
				labels[p] = NoiseLabel
			}
		} else {
			labels[p] = clusterID
		}
	}
	return labels
}
