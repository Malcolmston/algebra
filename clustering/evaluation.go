package clustering

import (
	"math"
	"math/rand"
)

// SilhouetteSamples returns the silhouette coefficient of every sample. For a
// sample i with mean intra-cluster distance a(i) and the smallest mean distance
// b(i) to samples of any other cluster, the coefficient is
// (b(i)-a(i))/max(a(i),b(i)). Singleton clusters give a coefficient of 0.
// metric may be nil (Euclidean).
func SilhouetteSamples(data [][]float64, labels []int, metric Metric) []float64 {
	if metric == nil {
		metric = Euclidean
	}
	n := len(data)
	dm := DistanceMatrix(data, metric)
	clusters := ClusterIndices(labels)
	out := make([]float64, n)
	for i := 0; i < n; i++ {
		own := clusters[labels[i]]
		// a(i): mean distance to other points in the same cluster.
		var a float64
		if len(own) > 1 {
			var sum float64
			for _, j := range own {
				if j != i {
					sum += dm[i][j]
				}
			}
			a = sum / float64(len(own)-1)
		}
		// b(i): min over other clusters of mean distance.
		b := math.Inf(1)
		for lbl, members := range clusters {
			if lbl == labels[i] {
				continue
			}
			var sum float64
			for _, j := range members {
				sum += dm[i][j]
			}
			mean := sum / float64(len(members))
			if mean < b {
				b = mean
			}
		}
		if len(own) <= 1 || math.IsInf(b, 1) {
			out[i] = 0
			continue
		}
		denom := math.Max(a, b)
		if denom == 0 {
			out[i] = 0
		} else {
			out[i] = (b - a) / denom
		}
	}
	return out
}

// SilhouetteScore returns the mean silhouette coefficient over all samples, a
// value in [-1, 1] where higher is better.
func SilhouetteScore(data [][]float64, labels []int, metric Metric) float64 {
	s := SilhouetteSamples(data, labels, metric)
	if len(s) == 0 {
		return 0
	}
	var sum float64
	for _, v := range s {
		sum += v
	}
	return sum / float64(len(s))
}

// DaviesBouldinIndex returns the Davies-Bouldin index of a clustering, the mean
// over clusters of the worst-case ratio of within-cluster scatter to
// between-cluster separation. Lower values indicate better clustering.
func DaviesBouldinIndex(data [][]float64, labels []int) float64 {
	clusters := ClusterIndices(labels)
	ids := UniqueLabels(labels)
	k := len(ids)
	if k < 2 {
		return 0
	}
	centroids := make(map[int][]float64, k)
	scatter := make(map[int]float64, k)
	for _, id := range ids {
		members := clusters[id]
		pts := make([][]float64, len(members))
		for i, m := range members {
			pts[i] = data[m]
		}
		c := Centroid(pts)
		centroids[id] = c
		var s float64
		for _, m := range members {
			s += Euclidean(data[m], c)
		}
		scatter[id] = s / float64(len(members))
	}
	var total float64
	for _, i := range ids {
		worst := math.Inf(-1)
		for _, j := range ids {
			if i == j {
				continue
			}
			sep := Euclidean(centroids[i], centroids[j])
			if sep == 0 {
				continue
			}
			r := (scatter[i] + scatter[j]) / sep
			if r > worst {
				worst = r
			}
		}
		if !math.IsInf(worst, -1) {
			total += worst
		}
	}
	return total / float64(k)
}

// CalinskiHarabaszIndex returns the Calinski-Harabasz index (variance ratio
// criterion): the ratio of between-cluster dispersion to within-cluster
// dispersion, scaled by the degrees of freedom. Higher values indicate better
// clustering.
func CalinskiHarabaszIndex(data [][]float64, labels []int) float64 {
	n := len(data)
	clusters := ClusterIndices(labels)
	ids := UniqueLabels(labels)
	k := len(ids)
	if k < 2 || n <= k {
		return 0
	}
	overall := Centroid(data)
	var between, within float64
	for _, id := range ids {
		members := clusters[id]
		pts := make([][]float64, len(members))
		for i, m := range members {
			pts[i] = data[m]
		}
		c := Centroid(pts)
		between += float64(len(members)) * SquaredEuclidean(c, overall)
		for _, m := range members {
			within += SquaredEuclidean(data[m], c)
		}
	}
	if within == 0 {
		return math.Inf(1)
	}
	return (between / within) * float64(n-k) / float64(k-1)
}

// DunnIndex returns the Dunn index of a clustering: the ratio of the minimum
// inter-cluster distance to the maximum intra-cluster diameter. Higher is
// better. metric may be nil (Euclidean).
func DunnIndex(data [][]float64, labels []int, metric Metric) float64 {
	if metric == nil {
		metric = Euclidean
	}
	clusters := ClusterIndices(labels)
	ids := UniqueLabels(labels)
	if len(ids) < 2 {
		return 0
	}
	// Max intra-cluster diameter.
	maxDiam := 0.0
	for _, id := range ids {
		members := clusters[id]
		for a := 0; a < len(members); a++ {
			for b := a + 1; b < len(members); b++ {
				if d := metric(data[members[a]], data[members[b]]); d > maxDiam {
					maxDiam = d
				}
			}
		}
	}
	if maxDiam == 0 {
		return math.Inf(1)
	}
	// Min inter-cluster distance.
	minSep := math.Inf(1)
	for ii := 0; ii < len(ids); ii++ {
		for jj := ii + 1; jj < len(ids); jj++ {
			mi := clusters[ids[ii]]
			mj := clusters[ids[jj]]
			for _, a := range mi {
				for _, b := range mj {
					if d := metric(data[a], data[b]); d < minSep {
						minSep = d
					}
				}
			}
		}
	}
	return minSep / maxDiam
}

// Inertia returns the total within-cluster sum of squared distances of samples
// to their cluster centroid (the k-means objective).
func Inertia(data [][]float64, labels []int) float64 {
	clusters := ClusterIndices(labels)
	var total float64
	for _, members := range clusters {
		pts := make([][]float64, len(members))
		for i, m := range members {
			pts[i] = data[m]
		}
		c := Centroid(pts)
		for _, m := range members {
			total += SquaredEuclidean(data[m], c)
		}
	}
	return total
}

// WithinClusterSumOfSquares returns the sum of squared distances to the centroid
// for each cluster, keyed by label.
func WithinClusterSumOfSquares(data [][]float64, labels []int) map[int]float64 {
	clusters := ClusterIndices(labels)
	out := make(map[int]float64, len(clusters))
	for lbl, members := range clusters {
		pts := make([][]float64, len(members))
		for i, m := range members {
			pts[i] = data[m]
		}
		c := Centroid(pts)
		var s float64
		for _, m := range members {
			s += SquaredEuclidean(data[m], c)
		}
		out[lbl] = s
	}
	return out
}

// BetweenClusterSumOfSquares returns the between-cluster sum of squares: the
// weighted squared distance of each cluster centroid to the overall mean.
func BetweenClusterSumOfSquares(data [][]float64, labels []int) float64 {
	overall := Centroid(data)
	clusters := ClusterIndices(labels)
	var total float64
	for _, members := range clusters {
		pts := make([][]float64, len(members))
		for i, m := range members {
			pts[i] = data[m]
		}
		c := Centroid(pts)
		total += float64(len(members)) * SquaredEuclidean(c, overall)
	}
	return total
}

// TotalSumOfSquares returns the total sum of squared distances of all samples to
// the overall mean.
func TotalSumOfSquares(data [][]float64) float64 {
	overall := Centroid(data)
	var s float64
	for _, row := range data {
		s += SquaredEuclidean(row, overall)
	}
	return s
}

// ElbowInertias runs k-means for each k in kValues and returns the resulting
// inertia for each, the sequence used to locate an "elbow" in the inertia curve.
func ElbowInertias(data [][]float64, kValues []int, rng *rand.Rand) ([]float64, error) {
	opts := KMeansOptions{NInit: 5, Rand: rng}
	out := make([]float64, len(kValues))
	for i, k := range kValues {
		res, err := KMeansWithOptions(data, k, opts)
		if err != nil {
			return nil, err
		}
		out[i] = res.Inertia
	}
	return out, nil
}

// ElbowPoint returns the index into kValues of the elbow of the inertia curve,
// chosen as the point of maximum distance to the straight line joining the first
// and last points (the "knee"/kneedle heuristic). Requires at least three
// points.
func ElbowPoint(kValues []int, inertias []float64) int {
	n := len(inertias)
	if n < 3 {
		return 0
	}
	x0, y0 := float64(kValues[0]), inertias[0]
	x1, y1 := float64(kValues[n-1]), inertias[n-1]
	dx := x1 - x0
	dy := y1 - y0
	norm := math.Hypot(dx, dy)
	best := 1
	bestDist := -1.0
	for i := 1; i < n-1; i++ {
		x := float64(kValues[i])
		y := inertias[i]
		// Perpendicular distance from point to the line.
		d := math.Abs(dy*x-dx*y+x1*y0-y1*x0) / norm
		if d > bestDist {
			bestDist = d
			best = i
		}
	}
	return best
}

// GapStatistic computes the gap statistic for a given k: the difference between
// the expected log within-cluster dispersion under a uniform reference
// distribution and the observed log dispersion. Larger gaps favour k. nRefs
// reference datasets are drawn from the bounding box of the data. It returns the
// gap value and its standard error.
func GapStatistic(data [][]float64, k, nRefs int, rng *rand.Rand) (gap, stderr float64, err error) {
	if rng == nil {
		rng = rand.New(rand.NewSource(1))
	}
	if nRefs <= 0 {
		nRefs = 10
	}
	opts := KMeansOptions{NInit: 5, Rand: rng}
	res, e := KMeansWithOptions(data, k, opts)
	if e != nil {
		return 0, 0, e
	}
	logWk := math.Log(res.Inertia + 1e-300)

	min, max := ColumnMinMax(data)
	dim := len(min)
	refLogs := make([]float64, nRefs)
	for r := 0; r < nRefs; r++ {
		ref := make([][]float64, len(data))
		for i := range ref {
			ref[i] = make([]float64, dim)
			for j := 0; j < dim; j++ {
				ref[i][j] = min[j] + rng.Float64()*(max[j]-min[j])
			}
		}
		rres, e := KMeansWithOptions(ref, k, opts)
		if e != nil {
			return 0, 0, e
		}
		refLogs[r] = math.Log(rres.Inertia + 1e-300)
	}
	var meanRef float64
	for _, v := range refLogs {
		meanRef += v
	}
	meanRef /= float64(nRefs)
	gap = meanRef - logWk
	var variance float64
	for _, v := range refLogs {
		d := v - meanRef
		variance += d * d
	}
	variance /= float64(nRefs)
	stderr = math.Sqrt(variance) * math.Sqrt(1+1/float64(nRefs))
	return gap, stderr, nil
}
