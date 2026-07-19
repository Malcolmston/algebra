package clustering

import (
	"math"
	"sort"
)

// Linkage identifies an inter-cluster linkage criterion for agglomerative
// clustering.
type Linkage int

const (
	// SingleLinkage merges clusters by the minimum pairwise distance (nearest
	// neighbour).
	SingleLinkage Linkage = iota
	// CompleteLinkage merges clusters by the maximum pairwise distance (farthest
	// neighbour).
	CompleteLinkage
	// AverageLinkage merges clusters by the mean pairwise distance (UPGMA).
	AverageLinkage
	// WardLinkage merges clusters so as to minimise the increase in total
	// within-cluster variance.
	WardLinkage
	// CentroidLinkage merges clusters by the distance between their centroids
	// (UPGMC).
	CentroidLinkage
	// MedianLinkage merges clusters using the WPGMC (median) update.
	MedianLinkage
	// WeightedLinkage merges clusters using the WPGMA (weighted average) update.
	WeightedLinkage
)

// String returns the canonical name of the linkage.
func (l Linkage) String() string {
	switch l {
	case SingleLinkage:
		return "single"
	case CompleteLinkage:
		return "complete"
	case AverageLinkage:
		return "average"
	case WardLinkage:
		return "ward"
	case CentroidLinkage:
		return "centroid"
	case MedianLinkage:
		return "median"
	case WeightedLinkage:
		return "weighted"
	default:
		return "unknown"
	}
}

// LinkageByName returns the Linkage for a name such as "single", "complete",
// "average"/"upgma", "ward", "centroid", "median" or "weighted".
func LinkageByName(name string) (Linkage, bool) {
	switch normalizeName(name) {
	case "single", "nearest":
		return SingleLinkage, true
	case "complete", "farthest", "maximum":
		return CompleteLinkage, true
	case "average", "upgma":
		return AverageLinkage, true
	case "ward":
		return WardLinkage, true
	case "centroid", "upgmc":
		return CentroidLinkage, true
	case "median", "wpgmc":
		return MedianLinkage, true
	case "weighted", "wpgma", "mcquitty":
		return WeightedLinkage, true
	default:
		return 0, false
	}
}

// MergeStep records a single agglomeration in a dendrogram. A and B are the ids
// of the merged clusters, Distance is the linkage distance at which they merged
// and Size is the number of original samples in the resulting cluster. Original
// samples have ids 0..n-1; the i-th merge creates cluster with id n+i.
type MergeStep struct {
	A        int
	B        int
	Distance float64
	Size     int
}

// Dendrogram is the result of agglomerative clustering: the ordered list of
// merges and the number of original observations.
type Dendrogram struct {
	// N is the number of original observations (leaves).
	N int
	// Merges holds the n-1 merge steps in increasing merge order.
	Merges []MergeStep
	// Linkage is the linkage criterion used.
	Linkage Linkage
}

// Linkage computes an agglomerative hierarchical clustering of data using the
// Lance-Williams update formula for the given linkage, and returns the
// dendrogram. The base distances use the supplied metric (nil means Euclidean);
// for Ward, centroid and median linkage the squared Euclidean distance is used
// internally regardless of metric, matching the usual definitions.
func LinkageCluster(data [][]float64, linkage Linkage, metric Metric) (*Dendrogram, error) {
	n := len(data)
	if n == 0 {
		return nil, ErrEmptyData
	}
	if metric == nil {
		metric = Euclidean
	}
	// Build initial distance matrix. Geometric linkages use squared Euclidean.
	geometric := linkage == WardLinkage || linkage == CentroidLinkage || linkage == MedianLinkage
	d := make([][]float64, n)
	for i := range d {
		d[i] = make([]float64, n)
	}
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			var dist float64
			if geometric {
				dist = SquaredEuclidean(data[i], data[j])
			} else {
				dist = metric(data[i], data[j])
			}
			d[i][j] = dist
			d[j][i] = dist
		}
	}
	return linkageFromMatrix(d, n, linkage)
}

// LinkageFromDistanceMatrix computes agglomerative clustering from a precomputed
// symmetric distance matrix. For Ward, centroid and median linkage the entries
// are treated as squared Euclidean distances.
func LinkageFromDistanceMatrix(dm [][]float64, linkage Linkage) (*Dendrogram, error) {
	n := len(dm)
	if n == 0 {
		return nil, ErrEmptyData
	}
	d := CloneMatrix(dm)
	return linkageFromMatrix(d, n, linkage)
}

func linkageFromMatrix(d [][]float64, n int, linkage Linkage) (*Dendrogram, error) {
	// active holds the currently live cluster ids; sizes maps id->size; the row
	// index in d corresponds to slot. We keep a slice of slots each holding a
	// cluster id.
	slotID := make([]int, n)
	sizes := make([]int, n)
	active := make([]bool, n)
	for i := 0; i < n; i++ {
		slotID[i] = i
		sizes[i] = 1
		active[i] = true
	}
	sizeByID := make(map[int]int, 2*n)
	for i := 0; i < n; i++ {
		sizeByID[i] = 1
	}
	merges := make([]MergeStep, 0, n-1)
	nextID := n

	for step := 0; step < n-1; step++ {
		// Find closest pair of active slots.
		bi, bj := -1, -1
		best := math.Inf(1)
		for i := 0; i < n; i++ {
			if !active[i] {
				continue
			}
			for j := i + 1; j < n; j++ {
				if !active[j] {
					continue
				}
				if d[i][j] < best {
					best = d[i][j]
					bi, bj = i, j
				}
			}
		}
		if bi < 0 {
			break
		}
		si := sizeByID[slotID[bi]]
		sj := sizeByID[slotID[bj]]
		reported := best
		if linkage == WardLinkage || linkage == CentroidLinkage || linkage == MedianLinkage {
			// Stored values are squared distances; report actual distance.
			reported = math.Sqrt(math.Max(0, best))
		}
		merges = append(merges, MergeStep{
			A:        slotID[bi],
			B:        slotID[bj],
			Distance: reported,
			Size:     si + sj,
		})

		// Update distances from the merged cluster (kept in slot bi) to all
		// other active slots using Lance-Williams.
		for m := 0; m < n; m++ {
			if !active[m] || m == bi || m == bj {
				continue
			}
			sk := sizeByID[slotID[m]]
			dik := d[bi][m]
			djk := d[bj][m]
			dij := best
			nd := lanceWilliams(linkage, dik, djk, dij, si, sj, sk)
			d[bi][m] = nd
			d[m][bi] = nd
		}
		// Deactivate slot bj; slot bi now represents the new cluster.
		active[bj] = false
		newID := nextID
		nextID++
		slotID[bi] = newID
		sizeByID[newID] = si + sj
	}

	return &Dendrogram{N: n, Merges: merges, Linkage: linkage}, nil
}

// lanceWilliams computes the updated distance between a merged cluster (i∪j) and
// another cluster k using the Lance-Williams recurrence for the given linkage.
func lanceWilliams(linkage Linkage, dik, djk, dij float64, si, sj, sk int) float64 {
	fi, fj, fk := float64(si), float64(sj), float64(sk)
	switch linkage {
	case SingleLinkage:
		return math.Min(dik, djk)
	case CompleteLinkage:
		return math.Max(dik, djk)
	case AverageLinkage:
		return (fi*dik + fj*djk) / (fi + fj)
	case WeightedLinkage:
		return 0.5*dik + 0.5*djk
	case WardLinkage:
		t := fi + fj + fk
		return ((fi+fk)*dik + (fj+fk)*djk - fk*dij) / t
	case CentroidLinkage:
		s := fi + fj
		return (fi*dik+fj*djk)/s - (fi*fj*dij)/(s*s)
	case MedianLinkage:
		return 0.5*dik + 0.5*djk - 0.25*dij
	default:
		return math.Min(dik, djk)
	}
}

// MergeDistances returns the merge distances in the order the merges occurred.
func (dg *Dendrogram) MergeDistances() []float64 {
	out := make([]float64, len(dg.Merges))
	for i, m := range dg.Merges {
		out[i] = m.Distance
	}
	return out
}

// CutTree cuts the dendrogram to produce exactly k flat clusters and returns a
// label per original observation in 0..k-1. k must be between 1 and N.
func (dg *Dendrogram) CutTree(k int) ([]int, error) {
	n := dg.N
	if k <= 0 || k > n {
		return nil, ErrInvalidK
	}
	// Union-find over ids 0..n-1+len(merges). Apply the first n-k merges.
	parent := make([]int, n+len(dg.Merges))
	for i := range parent {
		parent[i] = i
	}
	var find func(int) int
	find = func(x int) int {
		for parent[x] != x {
			parent[x] = parent[parent[x]]
			x = parent[x]
		}
		return x
	}
	nMerges := n - k
	for i := 0; i < nMerges; i++ {
		m := dg.Merges[i]
		newID := n + i
		ra := find(m.A)
		rb := find(m.B)
		parent[ra] = newID
		parent[rb] = newID
	}
	// Map roots (over the leaves) to compact labels.
	rootLabel := make(map[int]int)
	labels := make([]int, n)
	next := 0
	for i := 0; i < n; i++ {
		r := find(i)
		if _, ok := rootLabel[r]; !ok {
			rootLabel[r] = next
			next++
		}
		labels[i] = rootLabel[r]
	}
	return labels, nil
}

// CutTreeByHeight cuts the dendrogram at the given linkage height and returns a
// flat clustering: merges with distance strictly less than height are applied,
// so every remaining cluster is separated by at least height.
func (dg *Dendrogram) CutTreeByHeight(height float64) []int {
	n := dg.N
	parent := make([]int, n+len(dg.Merges))
	for i := range parent {
		parent[i] = i
	}
	var find func(int) int
	find = func(x int) int {
		for parent[x] != x {
			parent[x] = parent[parent[x]]
			x = parent[x]
		}
		return x
	}
	for i, m := range dg.Merges {
		if m.Distance >= height {
			break
		}
		newID := n + i
		parent[find(m.A)] = newID
		parent[find(m.B)] = newID
	}
	rootLabel := make(map[int]int)
	labels := make([]int, n)
	next := 0
	for i := 0; i < n; i++ {
		r := find(i)
		if _, ok := rootLabel[r]; !ok {
			rootLabel[r] = next
			next++
		}
		labels[i] = rootLabel[r]
	}
	return labels
}

// Leaves returns the original observation indices belonging to the cluster with
// the given dendrogram id (0..N-1 for leaves, N..2N-2 for internal nodes).
func (dg *Dendrogram) Leaves(id int) []int {
	n := dg.N
	if id < n {
		return []int{id}
	}
	m := dg.Merges[id-n]
	out := append(dg.Leaves(m.A), dg.Leaves(m.B)...)
	sort.Ints(out)
	return out
}

// CopheneticDistances returns the cophenetic distance matrix: the height at
// which each pair of original observations is first joined in the dendrogram.
func (dg *Dendrogram) CopheneticDistances() [][]float64 {
	n := dg.N
	c := make([][]float64, n)
	for i := range c {
		c[i] = make([]float64, n)
	}
	for i, m := range dg.Merges {
		left := dg.Leaves(m.A)
		right := dg.Leaves(m.B)
		_ = i
		for _, a := range left {
			for _, b := range right {
				c[a][b] = m.Distance
				c[b][a] = m.Distance
			}
		}
	}
	return c
}

// CopheneticCorrelation returns the Pearson correlation between the original
// pairwise distances (condensed) and the cophenetic distances from the
// dendrogram, a common measure of how faithfully a dendrogram preserves the
// input distances.
func (dg *Dendrogram) CopheneticCorrelation(originalCondensed []float64) float64 {
	n := dg.N
	coph := dg.CopheneticDistances()
	cond := make([]float64, 0, len(originalCondensed))
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			cond = append(cond, coph[i][j])
		}
	}
	return PearsonCorrelation(originalCondensed, cond)
}

// AgglomerativeClustering is a convenience wrapper that clusters data with the
// given linkage and metric and directly returns a flat clustering into k
// clusters.
func AgglomerativeClustering(data [][]float64, k int, linkage Linkage, metric Metric) ([]int, error) {
	dg, err := LinkageCluster(data, linkage, metric)
	if err != nil {
		return nil, err
	}
	return dg.CutTree(k)
}
