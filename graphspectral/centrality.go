package graphspectral

import "math"

// DegreeCentrality returns the degree centrality of every vertex: its weighted
// degree divided by n-1, so that in an unweighted graph a vertex adjacent to all
// others scores 1. For n < 2 the raw degrees are returned.
func DegreeCentrality(g *Graph) []float64 {
	out := g.Degrees()
	if g.n < 2 {
		return out
	}
	den := float64(g.n - 1)
	for i := range out {
		out[i] /= den
	}
	return out
}

// EigenvectorCentrality returns the eigenvector centrality: the (non-negative)
// dominant eigenvector of the adjacency matrix, scaled to unit Euclidean norm.
// It uses the power method with the given iteration budget and tolerance;
// defaults of 1000 and 1e-12 apply when non-positive values are passed. For a
// connected graph the Perron-Frobenius theorem guarantees a positive vector.
func EigenvectorCentrality(g *Graph, maxIter int, tol float64) ([]float64, error) {
	if g.n == 0 {
		return nil, ErrEmpty
	}
	if maxIter <= 0 {
		maxIter = 1000
	}
	if tol <= 0 {
		tol = 1e-12
	}
	_, vec, err := PowerIteration(g.Adjacency(), maxIter, tol)
	if err != nil && err != ErrNoConvergence {
		return nil, err
	}
	// make the vector non-negative (flip sign if predominantly negative)
	var s float64
	for _, x := range vec {
		s += x
	}
	if s < 0 {
		vec = VecScale(vec, -1)
	}
	// clamp tiny negatives from round-off and renormalize
	for i := range vec {
		if vec[i] < 0 && vec[i] > -tol*10 {
			vec[i] = 0
		}
	}
	return Normalize(vec), nil
}

// KatzCentrality returns the Katz centrality x = (I - α·A)^{-1}·β·1, solving the
// corresponding linear system. For convergence of the underlying series α must
// be smaller than the reciprocal of the spectral radius; the linear solve
// itself only requires I - α·A to be nonsingular. It returns ErrInvalidArgument
// for a non-positive beta.
func KatzCentrality(g *Graph, alpha, beta float64) ([]float64, error) {
	if g.n == 0 {
		return nil, ErrEmpty
	}
	if beta <= 0 {
		return nil, ErrInvalidArgument
	}
	a := g.Adjacency()
	m := IdentityMatrix(g.n)
	for i := 0; i < g.n; i++ {
		for j := 0; j < g.n; j++ {
			m.Add(i, j, -alpha*a.At(i, j))
		}
	}
	b := make([]float64, g.n)
	for i := range b {
		b[i] = beta
	}
	return SolveLinear(m, b)
}

// ClosenessCentrality returns the closeness centrality of every vertex using
// weighted shortest-path distances. For a vertex that reaches r-1 others the
// value is (r-1)/Σd · (r-1)/(n-1), the Wasserman-Faust normalization that keeps
// disconnected graphs meaningful. Isolated vertices score 0.
func ClosenessCentrality(g *Graph) []float64 {
	out := make([]float64, g.n)
	if g.n < 2 {
		return out
	}
	for v := 0; v < g.n; v++ {
		dist := DijkstraDistances(g, v)
		var sum float64
		reach := 0
		for u := 0; u < g.n; u++ {
			if u == v || math.IsInf(dist[u], 1) {
				continue
			}
			sum += dist[u]
			reach++
		}
		if sum > 0 {
			out[v] = (float64(reach) / sum) * (float64(reach) / float64(g.n-1))
		}
	}
	return out
}

// HarmonicCentrality returns the harmonic centrality of every vertex,
// Σ_{u≠v} 1/d(v,u) divided by n-1, using weighted shortest-path distances.
// Unreachable vertices contribute nothing, so the measure is well defined on
// disconnected graphs.
func HarmonicCentrality(g *Graph) []float64 {
	out := make([]float64, g.n)
	if g.n < 2 {
		return out
	}
	for v := 0; v < g.n; v++ {
		dist := DijkstraDistances(g, v)
		var s float64
		for u := 0; u < g.n; u++ {
			if u == v || math.IsInf(dist[u], 1) {
				continue
			}
			s += 1 / dist[u]
		}
		out[v] = s / float64(g.n-1)
	}
	return out
}

// PageRank returns the PageRank vector of the graph by power iteration on the
// Google matrix with the given damping factor (typically 0.85). Dangling
// vertices (degree 0) redistribute their mass uniformly. Defaults of maxIter =
// 1000 and tol = 1e-12 apply when non-positive values are given. The result is a
// probability distribution summing to 1. It returns ErrInvalidArgument if
// damping is not in [0, 1).
func PageRank(g *Graph, damping float64, maxIter int, tol float64) ([]float64, error) {
	uniform := make([]float64, g.n)
	for i := range uniform {
		uniform[i] = 1 / float64(g.n)
	}
	return PersonalizedPageRank(g, damping, uniform, maxIter, tol)
}

// PersonalizedPageRank returns PageRank with a custom teleportation
// distribution: at each step the random surfer restarts according to restart
// (which need not be normalized; it is normalized internally) instead of
// uniformly. It returns ErrDimensionMismatch if len(restart) != n and
// ErrInvalidArgument if damping is not in [0, 1) or restart sums to zero.
func PersonalizedPageRank(g *Graph, damping float64, restart []float64, maxIter int, tol float64) ([]float64, error) {
	if g.n == 0 {
		return nil, ErrEmpty
	}
	if damping < 0 || damping >= 1 {
		return nil, ErrInvalidArgument
	}
	if len(restart) != g.n {
		return nil, ErrDimensionMismatch
	}
	if maxIter <= 0 {
		maxIter = 1000
	}
	if tol <= 0 {
		tol = 1e-12
	}
	tele := NormalizeL1(restart)
	if VecSum(tele) == 0 {
		return nil, ErrInvalidArgument
	}
	deg := g.Degrees()
	r := VecClone(tele)
	for it := 0; it < maxIter; it++ {
		next := make([]float64, g.n)
		var dangling float64
		for i := 0; i < g.n; i++ {
			if deg[i] == 0 {
				dangling += r[i]
			}
		}
		for j := 0; j < g.n; j++ {
			var inflow float64
			for i := 0; i < g.n; i++ {
				if deg[i] == 0 {
					continue
				}
				w := g.adj.At(i, j)
				if w != 0 {
					inflow += r[i] * w / deg[i]
				}
			}
			next[j] = (1-damping)*tele[j] + damping*(inflow+dangling*tele[j])
		}
		var diff float64
		for i := 0; i < g.n; i++ {
			diff += math.Abs(next[i] - r[i])
		}
		r = next
		if diff < tol {
			break
		}
	}
	// normalize against accumulated round-off
	return NormalizeL1(r), nil
}
