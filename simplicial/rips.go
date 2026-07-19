package simplicial

// buildAdjacency returns the boolean adjacency matrix in which i and j are
// adjacent when i≠j and their distance does not exceed threshold.
func buildAdjacency(dist [][]float64, threshold float64) [][]bool {
	n := len(dist)
	adj := make([][]bool, n)
	for i := range adj {
		adj[i] = make([]bool, n)
		for j := 0; j < n; j++ {
			if i != j && dist[i][j] <= threshold {
				adj[i][j] = true
			}
		}
	}
	return adj
}

// VietorisRipsFromDistances builds the Vietoris–Rips complex at scale epsilon
// directly from a symmetric distance matrix, up to dimension maxDim. A set of
// vertices spans a simplex exactly when all of its pairwise distances are at
// most epsilon (the clique complex of the neighbourhood graph). Every vertex is
// always included.
func VietorisRipsFromDistances(dist [][]float64, epsilon float64, maxDim int) *Complex {
	n := len(dist)
	c := NewComplex()
	for i := 0; i < n; i++ {
		c.AddVertex(i)
	}
	if maxDim < 1 {
		return c
	}
	adj := buildAdjacency(dist, epsilon)

	var rec func(clique, cand []int)
	rec = func(clique, cand []int) {
		if len(clique) == maxDim+1 {
			return
		}
		for idx, w := range cand {
			next := append(append([]int(nil), clique...), w)
			c.AddSimplex(NewSimplex(next...))
			var newcand []int
			for _, x := range cand[idx+1:] {
				if adj[w][x] {
					newcand = append(newcand, x)
				}
			}
			rec(next, newcand)
		}
	}
	for i := 0; i < n; i++ {
		var cand []int
		for j := i + 1; j < n; j++ {
			if adj[i][j] {
				cand = append(cand, j)
			}
		}
		rec([]int{i}, cand)
	}
	return c
}

// VietorisRips builds the Vietoris–Rips complex of a point cloud at scale
// epsilon under the Euclidean metric, up to dimension maxDim.
func VietorisRips(pc *PointCloud, epsilon float64, maxDim int) *Complex {
	return VietorisRipsFromDistances(pc.DistanceMatrix(), epsilon, maxDim)
}

// VietorisRipsMetric builds the Vietoris–Rips complex of a point cloud at scale
// epsilon under an arbitrary metric, up to dimension maxDim.
func VietorisRipsMetric(pc *PointCloud, epsilon float64, maxDim int, metric Metric) *Complex {
	return VietorisRipsFromDistances(pc.DistanceMatrixWith(metric), epsilon, maxDim)
}

// NeighborhoodGraph returns the 1-skeleton (graph) of the Vietoris–Rips complex
// of a point cloud at scale epsilon: vertices for every point and an edge
// between points at Euclidean distance at most epsilon.
func NeighborhoodGraph(pc *PointCloud, epsilon float64) *Complex {
	return VietorisRips(pc, epsilon, 1)
}

// Cech builds the Čech complex of a point cloud at radius r up to dimension
// maxDim. A set of points spans a simplex exactly when the closed balls of
// radius r around them share a common point — equivalently when the minimal
// enclosing ball of the points has radius at most r. The Čech complex is a
// subcomplex of the Vietoris–Rips complex at scale 2r and, unlike it, is
// homotopy equivalent to the union of the balls (the nerve lemma).
func Cech(pc *PointCloud, r float64, maxDim int) *Complex {
	n := pc.Len()
	c := NewComplex()
	for i := 0; i < n; i++ {
		c.AddVertex(i)
	}
	if maxDim < 1 || r < 0 {
		return c
	}
	dist := pc.DistanceMatrix()
	adj := buildAdjacency(dist, 2*r)

	gather := func(clique []int) [][]float64 {
		out := make([][]float64, len(clique))
		for i, v := range clique {
			out[i] = pc.points[v]
		}
		return out
	}

	var rec func(clique, cand []int)
	rec = func(clique, cand []int) {
		if len(clique) == maxDim+1 {
			return
		}
		for idx, w := range cand {
			next := append(append([]int(nil), clique...), w)
			if _, rr := MinimalEnclosingBall(gather(next)); rr > r+mebEps {
				continue
			}
			c.AddSimplex(NewSimplex(next...))
			var newcand []int
			for _, x := range cand[idx+1:] {
				if adj[w][x] {
					newcand = append(newcand, x)
				}
			}
			rec(next, newcand)
		}
	}
	for i := 0; i < n; i++ {
		var cand []int
		for j := i + 1; j < n; j++ {
			if adj[i][j] {
				cand = append(cand, j)
			}
		}
		rec([]int{i}, cand)
	}
	return c
}
