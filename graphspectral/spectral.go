package graphspectral

import "math"

// AdjacencySpectrum returns the eigenvalues of the adjacency matrix in ascending
// order.
func AdjacencySpectrum(g *Graph) ([]float64, error) {
	return Eigenvalues(g.Adjacency())
}

// LaplacianSpectrum returns the eigenvalues of the combinatorial Laplacian in
// ascending order. The smallest is always (numerically) zero.
func LaplacianSpectrum(g *Graph) ([]float64, error) {
	return Eigenvalues(g.Laplacian())
}

// NormalizedLaplacianSpectrum returns the eigenvalues of the symmetric
// normalized Laplacian in ascending order. All lie in [0, 2].
func NormalizedLaplacianSpectrum(g *Graph) ([]float64, error) {
	return Eigenvalues(g.NormalizedLaplacian())
}

// SignlessLaplacianSpectrum returns the eigenvalues of the signless Laplacian
// Q = D + A in ascending order.
func SignlessLaplacianSpectrum(g *Graph) ([]float64, error) {
	return Eigenvalues(g.SignlessLaplacian())
}

// AlgebraicConnectivity returns the second-smallest eigenvalue of the Laplacian
// (the Fiedler value). It is positive exactly when the graph is connected. It
// returns ErrEmpty for graphs with fewer than two vertices.
func AlgebraicConnectivity(g *Graph) (float64, error) {
	if g.n < 2 {
		return 0, ErrEmpty
	}
	vals, err := LaplacianSpectrum(g)
	if err != nil {
		return 0, err
	}
	return vals[1], nil
}

// FiedlerValue is a synonym for [AlgebraicConnectivity].
func FiedlerValue(g *Graph) (float64, error) { return AlgebraicConnectivity(g) }

// FiedlerVector returns the eigenvector of the Laplacian associated with the
// second-smallest eigenvalue. It returns ErrEmpty for graphs with fewer than two
// vertices.
func FiedlerVector(g *Graph) ([]float64, error) {
	if g.n < 2 {
		return nil, ErrEmpty
	}
	e, err := EigenSymmetric(g.Laplacian())
	if err != nil {
		return nil, err
	}
	e.SortAscending()
	return e.Vector(1), nil
}

// SpectralGap returns the adjacency spectral gap: the difference between the
// largest and second-largest adjacency eigenvalues. It returns ErrEmpty for
// graphs with fewer than two vertices.
func SpectralGap(g *Graph) (float64, error) {
	if g.n < 2 {
		return 0, ErrEmpty
	}
	vals, err := AdjacencySpectrum(g)
	if err != nil {
		return 0, err
	}
	n := len(vals)
	return vals[n-1] - vals[n-2], nil
}

// LaplacianSpectralGap returns the difference between the second-smallest and
// smallest Laplacian eigenvalues. For a connected graph this equals the
// algebraic connectivity.
func LaplacianSpectralGap(g *Graph) (float64, error) {
	if g.n < 2 {
		return 0, ErrEmpty
	}
	vals, err := LaplacianSpectrum(g)
	if err != nil {
		return 0, err
	}
	return vals[1] - vals[0], nil
}

// SpectralRadiusAdjacency returns the largest absolute adjacency eigenvalue,
// also known as the spectral radius of the graph.
func SpectralRadiusAdjacency(g *Graph) (float64, error) {
	return SpectralRadius(g.Adjacency())
}

// GraphEnergy returns the graph energy: the sum of the absolute values of the
// adjacency eigenvalues.
func GraphEnergy(g *Graph) (float64, error) {
	vals, err := AdjacencySpectrum(g)
	if err != nil {
		return 0, err
	}
	var s float64
	for _, x := range vals {
		s += math.Abs(x)
	}
	return s, nil
}

// LaplacianEnergy returns the Laplacian energy Σ |μ_i - 2m/n|, where the μ_i are
// the Laplacian eigenvalues and 2m/n is the average (weighted) degree.
func LaplacianEnergy(g *Graph) (float64, error) {
	vals, err := LaplacianSpectrum(g)
	if err != nil {
		return 0, err
	}
	if g.n == 0 {
		return 0, nil
	}
	avg := VecSum(vals) / float64(g.n)
	var s float64
	for _, x := range vals {
		s += math.Abs(x - avg)
	}
	return s, nil
}

// EstradaIndex returns the Estrada index Σ exp(λ_i) over the adjacency
// eigenvalues.
func EstradaIndex(g *Graph) (float64, error) {
	vals, err := AdjacencySpectrum(g)
	if err != nil {
		return 0, err
	}
	var s float64
	for _, x := range vals {
		s += math.Exp(x)
	}
	return s, nil
}

// SpanningTreeCount returns the (weighted) number of spanning trees via the
// Matrix-Tree theorem: the determinant of the Laplacian with its last row and
// column deleted. For an unweighted connected graph this is an integer. It
// returns 0 for a disconnected graph and ErrEmpty for a graph with no vertices.
func SpanningTreeCount(g *Graph) (float64, error) {
	if g.n == 0 {
		return 0, ErrEmpty
	}
	if g.n == 1 {
		return 1, nil
	}
	l := g.Laplacian()
	drop := map[int]bool{g.n - 1: true}
	minor := l.SubMatrix(drop, drop)
	return Determinant(minor)
}

// NumberOfSpanningTrees returns [SpanningTreeCount] rounded to the nearest
// integer, appropriate for unweighted graphs.
func NumberOfSpanningTrees(g *Graph) (int, error) {
	c, err := SpanningTreeCount(g)
	if err != nil {
		return 0, err
	}
	return int(math.Round(c)), nil
}

// LaplacianPseudoinverse returns the Moore-Penrose pseudoinverse L⁺ of the
// combinatorial Laplacian, computed from its eigendecomposition by inverting the
// eigenvalues that exceed tol and zeroing the rest. If tol <= 0 a default of
// 1e-9 is used.
func LaplacianPseudoinverse(g *Graph, tol float64) (*Matrix, error) {
	if tol <= 0 {
		tol = 1e-9
	}
	e, err := EigenSymmetric(g.Laplacian())
	if err != nil {
		return nil, err
	}
	n := g.n
	out := NewMatrix(n, n)
	for k := 0; k < e.Len(); k++ {
		lam := e.Values[k]
		if lam <= tol {
			continue
		}
		u := e.Vector(k)
		inv := 1 / lam
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				out.Add(i, j, inv*u[i]*u[j])
			}
		}
	}
	return out, nil
}

// EffectiveResistance returns the effective (electrical) resistance between
// vertices s and t, treating each edge as a unit conductance times its weight.
// It is computed as L⁺_ss + L⁺_tt - 2·L⁺_st.
func EffectiveResistance(g *Graph, s, t int) (float64, error) {
	if s < 0 || s >= g.n || t < 0 || t >= g.n {
		return 0, ErrOutOfRange
	}
	lp, err := LaplacianPseudoinverse(g, 1e-9)
	if err != nil {
		return 0, err
	}
	return lp.At(s, s) + lp.At(t, t) - 2*lp.At(s, t), nil
}

// KirchhoffIndex returns the Kirchhoff index: the sum of effective resistances
// over all unordered vertex pairs, equal to n·Σ 1/μ_i over the nonzero Laplacian
// eigenvalues.
func KirchhoffIndex(g *Graph) (float64, error) {
	vals, err := LaplacianSpectrum(g)
	if err != nil {
		return 0, err
	}
	var s float64
	for _, x := range vals {
		if x > 1e-9 {
			s += 1 / x
		}
	}
	return float64(g.n) * s, nil
}

// NumberOfComponentsSpectral returns the number of connected components inferred
// from the spectrum: the multiplicity of the zero eigenvalue of the Laplacian
// (eigenvalues below tol are treated as zero). If tol <= 0 a default of 1e-8 is
// used.
func NumberOfComponentsSpectral(g *Graph, tol float64) (int, error) {
	if tol <= 0 {
		tol = 1e-8
	}
	vals, err := LaplacianSpectrum(g)
	if err != nil {
		return 0, err
	}
	c := 0
	for _, x := range vals {
		if x < tol {
			c++
		}
	}
	return c, nil
}

// SpectralOrdering returns a permutation of the vertices sorted by their Fiedler
// vector value in ascending order. Reordering a graph this way tends to place
// well-connected vertices near one another.
func SpectralOrdering(g *Graph) ([]int, error) {
	f, err := FiedlerVector(g)
	if err != nil {
		return nil, err
	}
	idx := make([]int, len(f))
	for i := range idx {
		idx[i] = i
	}
	// simple insertion sort by f value keeps it stable and dependency-free
	for i := 1; i < len(idx); i++ {
		for j := i; j > 0 && f[idx[j-1]] > f[idx[j]]; j-- {
			idx[j-1], idx[j] = idx[j], idx[j-1]
		}
	}
	return idx, nil
}
