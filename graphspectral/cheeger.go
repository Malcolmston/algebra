package graphspectral

import "math"

// EdgeBoundary returns the total weight of edges with exactly one endpoint in
// the vertex set S (given as a membership map). This is the size of the cut
// (S, V\S).
func EdgeBoundary(g *Graph, s map[int]bool) float64 {
	var w float64
	for i := 0; i < g.n; i++ {
		if !s[i] {
			continue
		}
		for j := 0; j < g.n; j++ {
			if i == j || s[j] {
				continue
			}
			w += g.adj.At(i, j)
		}
	}
	return w
}

// Volume returns the sum of the weighted degrees of the vertices in S.
func Volume(g *Graph, s map[int]bool) float64 {
	var v float64
	for i := 0; i < g.n; i++ {
		if s[i] {
			v += g.WeightedDegree(i)
		}
	}
	return v
}

// Conductance returns the conductance of the cut (S, V\S):
// |∂S| / min(vol(S), vol(V\S)). It returns 0 when S is empty or all of V, or
// when the smaller side has zero volume.
func Conductance(g *Graph, s map[int]bool) float64 {
	count := 0
	for i := 0; i < g.n; i++ {
		if s[i] {
			count++
		}
	}
	if count == 0 || count == g.n {
		return 0
	}
	comp := make(map[int]bool, g.n-count)
	for i := 0; i < g.n; i++ {
		if !s[i] {
			comp[i] = true
		}
	}
	volS := Volume(g, s)
	volC := Volume(g, comp)
	den := math.Min(volS, volC)
	if den == 0 {
		return 0
	}
	return EdgeBoundary(g, s) / den
}

// CutSize returns the total weight of edges crossing between the two blocks of a
// bipartition given by part (entries 0 or 1).
func CutSize(g *Graph, part []int) float64 {
	var w float64
	for i := 0; i < g.n; i++ {
		for j := i + 1; j < g.n; j++ {
			if part[i] != part[j] {
				w += g.adj.At(i, j)
			}
		}
	}
	return w
}

// NormalizedCut returns the normalized-cut value of a bipartition (part entries
// 0 or 1): cut(A,B)/vol(A) + cut(A,B)/vol(B).
func NormalizedCut(g *Graph, part []int) float64 {
	a := make(map[int]bool)
	b := make(map[int]bool)
	for i := 0; i < g.n; i++ {
		if part[i] == 0 {
			a[i] = true
		} else {
			b[i] = true
		}
	}
	cut := EdgeBoundary(g, a)
	volA := Volume(g, a)
	volB := Volume(g, b)
	var nc float64
	if volA > 0 {
		nc += cut / volA
	}
	if volB > 0 {
		nc += cut / volB
	}
	return nc
}

// CheegerConstant returns the exact Cheeger (conductance) constant, the minimum
// conductance over all nonempty proper vertex subsets, computed by brute force.
// Because it enumerates all 2^n subsets it is intended for small graphs; it
// returns ErrInvalidArgument when n exceeds 22.
func CheegerConstant(g *Graph) (float64, error) {
	if g.n > 22 {
		return 0, ErrInvalidArgument
	}
	if g.n < 2 {
		return 0, ErrEmpty
	}
	best := math.Inf(1)
	total := 1 << uint(g.n)
	for mask := 1; mask < total-1; mask++ {
		s := make(map[int]bool)
		for i := 0; i < g.n; i++ {
			if mask&(1<<uint(i)) != 0 {
				s[i] = true
			}
		}
		c := Conductance(g, s)
		if c > 0 && c < best {
			best = c
		}
	}
	if math.IsInf(best, 1) {
		return 0, nil
	}
	return best, nil
}

// CheegerLowerBound returns λ/2, the lower Cheeger bound on the conductance
// constant, where λ is the second-smallest eigenvalue of the normalized
// Laplacian. It satisfies λ/2 <= h(G).
func CheegerLowerBound(g *Graph) (float64, error) {
	lam, err := spectralGapNormalized(g)
	if err != nil {
		return 0, err
	}
	return lam / 2, nil
}

// CheegerUpperBound returns sqrt(2λ), the upper Cheeger bound on the conductance
// constant, where λ is the second-smallest eigenvalue of the normalized
// Laplacian. It satisfies h(G) <= sqrt(2λ).
func CheegerUpperBound(g *Graph) (float64, error) {
	lam, err := spectralGapNormalized(g)
	if err != nil {
		return 0, err
	}
	return math.Sqrt(2 * lam), nil
}

// spectralGapNormalized returns the second-smallest normalized Laplacian
// eigenvalue.
func spectralGapNormalized(g *Graph) (float64, error) {
	if g.n < 2 {
		return 0, ErrEmpty
	}
	vals, err := NormalizedLaplacianSpectrum(g)
	if err != nil {
		return 0, err
	}
	return vals[1], nil
}

// NormalizedAlgebraicConnectivity returns the second-smallest eigenvalue of the
// normalized Laplacian, the spectral quantity appearing in the Cheeger
// inequality.
func NormalizedAlgebraicConnectivity(g *Graph) (float64, error) {
	return spectralGapNormalized(g)
}
