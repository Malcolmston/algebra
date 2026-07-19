package graphspectral

import "math"

// SpectralBisection partitions the vertices into two blocks using the sign of
// the Fiedler vector: vertices with a negative Fiedler coordinate are labelled 0
// and the rest are labelled 1. This is the classical spectral graph bisection.
func SpectralBisection(g *Graph) ([]int, error) {
	f, err := FiedlerVector(g)
	if err != nil {
		return nil, err
	}
	part := make([]int, len(f))
	for i, x := range f {
		if x >= 0 {
			part[i] = 1
		}
	}
	return part, nil
}

// SpectralBisectionMedian partitions the vertices into two balanced blocks by
// thresholding the Fiedler vector at its median: the vertices whose Fiedler
// coordinate is above the median form block 1, the rest block 0.
func SpectralBisectionMedian(g *Graph) ([]int, error) {
	f, err := FiedlerVector(g)
	if err != nil {
		return nil, err
	}
	med := median(f)
	part := make([]int, len(f))
	for i, x := range f {
		if x > med {
			part[i] = 1
		}
	}
	return part, nil
}

func median(v []float64) float64 {
	s := SortedFloats(v)
	n := len(s)
	if n == 0 {
		return 0
	}
	if n%2 == 1 {
		return s[n/2]
	}
	return 0.5 * (s[n/2-1] + s[n/2])
}

// SpectralClustering partitions the vertices into k clusters by embedding each
// vertex into the space spanned by the k eigenvectors of the normalized
// Laplacian with the smallest eigenvalues, then running deterministic k-means
// (Lloyd's algorithm with farthest-first seeding) on the embedded points. It
// returns a label in [0, k) for every vertex. It returns ErrInvalidArgument if k
// is not in [1, n].
func SpectralClustering(g *Graph, k int) ([]int, error) {
	if k < 1 || k > g.n {
		return nil, ErrInvalidArgument
	}
	if k == 1 {
		return make([]int, g.n), nil
	}
	e, err := EigenSymmetric(g.NormalizedLaplacian())
	if err != nil {
		return nil, err
	}
	e.SortAscending()
	// build the n-by-k embedding from the k smallest eigenvectors
	points := make([][]float64, g.n)
	for i := 0; i < g.n; i++ {
		row := make([]float64, k)
		for c := 0; c < k; c++ {
			row[c] = e.Vectors.At(i, c)
		}
		points[i] = row
	}
	return KMeans(points, k, 100), nil
}

// KMeans runs Lloyd's algorithm on the given points (each a fixed-length
// coordinate slice) with k clusters and at most maxIter iterations, returning a
// cluster label in [0, k) for every point. Initial centres are chosen
// deterministically by farthest-first traversal starting from the first point,
// so the result is reproducible. It returns nil if k is not in [1, len(points)].
func KMeans(points [][]float64, k, maxIter int) []int {
	n := len(points)
	if k < 1 || k > n {
		return nil
	}
	if maxIter <= 0 {
		maxIter = 100
	}
	// farthest-first seeding
	centers := make([][]float64, k)
	centers[0] = VecClone(points[0])
	for c := 1; c < k; c++ {
		bestIdx, bestDist := -1, -1.0
		for i := 0; i < n; i++ {
			d := math.Inf(1)
			for j := 0; j < c; j++ {
				if dd := VecDistance(points[i], centers[j]); dd < d {
					d = dd
				}
			}
			if d > bestDist {
				bestDist = d
				bestIdx = i
			}
		}
		centers[c] = VecClone(points[bestIdx])
	}
	labels := make([]int, n)
	for it := 0; it < maxIter; it++ {
		changed := false
		for i := 0; i < n; i++ {
			best, bd := 0, math.Inf(1)
			for c := 0; c < k; c++ {
				if d := VecDistance(points[i], centers[c]); d < bd {
					bd = d
					best = c
				}
			}
			if labels[i] != best {
				labels[i] = best
				changed = true
			}
		}
		// recompute centres
		dim := len(points[0])
		sums := make([][]float64, k)
		counts := make([]int, k)
		for c := range sums {
			sums[c] = make([]float64, dim)
		}
		for i := 0; i < n; i++ {
			c := labels[i]
			counts[c]++
			for d := 0; d < dim; d++ {
				sums[c][d] += points[i][d]
			}
		}
		for c := 0; c < k; c++ {
			if counts[c] == 0 {
				continue
			}
			for d := 0; d < dim; d++ {
				sums[c][d] /= float64(counts[c])
			}
			centers[c] = sums[c]
		}
		if !changed && it > 0 {
			break
		}
	}
	return labels
}
