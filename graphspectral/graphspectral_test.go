package graphspectral

import (
	"fmt"
	"math"
	"testing"
)

const tol = 1e-6

func approx(a, b float64) bool { return math.Abs(a-b) <= tol }

// sortedApprox reports whether two slices agree entrywise after both are sorted.
func sortedApprox(a, b []float64, eps float64) bool {
	if len(a) != len(b) {
		return false
	}
	sa := SortedFloats(a)
	sb := SortedFloats(b)
	for i := range sa {
		if math.Abs(sa[i]-sb[i]) > eps {
			return false
		}
	}
	return true
}

// ---------------------------------------------------------------------------
// Linear algebra
// ---------------------------------------------------------------------------

func TestMatrixMul(t *testing.T) {
	a, _ := NewMatrixFromRows([][]float64{{1, 2}, {3, 4}})
	b, _ := NewMatrixFromRows([][]float64{{5, 6}, {7, 8}})
	want, _ := NewMatrixFromRows([][]float64{{19, 22}, {43, 50}})
	got, err := a.Mul(b)
	if err != nil {
		t.Fatal(err)
	}
	if !got.ApproxEqual(want, tol) {
		t.Fatalf("Mul = %v, want %v", got, want)
	}
}

func TestDeterminantSolveInverse(t *testing.T) {
	a, _ := NewMatrixFromRows([][]float64{{2, 1, 1}, {1, 3, 2}, {1, 0, 0}})
	det, err := Determinant(a)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(det, -1) {
		t.Fatalf("det = %v, want -1", det)
	}
	inv, err := Inverse(a)
	if err != nil {
		t.Fatal(err)
	}
	prod, _ := a.Mul(inv)
	if !prod.ApproxEqual(IdentityMatrix(3), 1e-9) {
		t.Fatalf("A·A^{-1} != I: %v", prod)
	}
	x, err := SolveLinear(a, []float64{4, 5, 6})
	if err != nil {
		t.Fatal(err)
	}
	chk, _ := a.MulVec(x)
	if !VecApproxEqual(chk, []float64{4, 5, 6}, 1e-9) {
		t.Fatalf("A·x = %v, want [4 5 6]", chk)
	}
}

func TestPowAndRank(t *testing.T) {
	a, _ := NewMatrixFromRows([][]float64{{2, 0}, {0, 3}})
	p, err := a.Pow(3)
	if err != nil {
		t.Fatal(err)
	}
	want, _ := NewMatrixFromRows([][]float64{{8, 0}, {0, 27}})
	if !p.ApproxEqual(want, tol) {
		t.Fatalf("A^3 = %v", p)
	}
	sing, _ := NewMatrixFromRows([][]float64{{1, 2}, {2, 4}})
	if r := Rank(sing, 1e-9); r != 1 {
		t.Fatalf("rank = %d, want 1", r)
	}
}

// ---------------------------------------------------------------------------
// Eigenvalues
// ---------------------------------------------------------------------------

func TestEigenSymmetric2x2(t *testing.T) {
	m, _ := NewMatrixFromRows([][]float64{{2, 1}, {1, 2}})
	e, err := EigenSymmetric(m)
	if err != nil {
		t.Fatal(err)
	}
	e.SortAscending()
	if !VecApproxEqual(e.Values, []float64{1, 3}, tol) {
		t.Fatalf("eigenvalues = %v, want [1 3]", e.Values)
	}
	// verify residual A·v - λ·v ≈ 0 for each pair
	for i := 0; i < e.Len(); i++ {
		lam, v := e.Pair(i)
		av, _ := m.MulVec(v)
		res := VecSub(av, VecScale(v, lam))
		if Norm2(res) > 1e-9 {
			t.Fatalf("residual for pair %d = %v", i, res)
		}
	}
}

func TestPowerIteration(t *testing.T) {
	m, _ := NewMatrixFromRows([][]float64{{2, 1}, {1, 2}})
	lam, v, err := PowerIteration(m, 1000, 1e-12)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(lam, 3) {
		t.Fatalf("dominant eigenvalue = %v, want 3", lam)
	}
	want := Normalize([]float64{1, 1})
	if !VecApproxEqual(v, want, 1e-6) {
		t.Fatalf("dominant eigenvector = %v, want %v", v, want)
	}
}

// ---------------------------------------------------------------------------
// Spectra of standard graphs
// ---------------------------------------------------------------------------

func TestLaplacianSpectra(t *testing.T) {
	s3 := math.Sqrt(3)
	tests := []struct {
		name    string
		g       *Graph
		lap     []float64
		adj     []float64
		trees   int
		algConn float64
		energy  float64
		specRad float64
	}{
		{"K4", CompleteGraph(4), []float64{0, 4, 4, 4}, []float64{3, -1, -1, -1}, 16, 4, 6, 3},
		{"C4", CycleGraph(4), []float64{0, 2, 2, 4}, []float64{2, 0, 0, -2}, 4, 2, 4, 2},
		{"P3", PathGraph(3), []float64{0, 1, 3}, []float64{-math.Sqrt2, 0, math.Sqrt2}, 1, 1, 2 * math.Sqrt2, math.Sqrt2},
		{"Star4", StarGraph(4), []float64{0, 1, 1, 4}, []float64{-s3, 0, 0, s3}, 1, 1, 2 * s3, s3},
		{"K5", CompleteGraph(5), []float64{0, 5, 5, 5, 5}, []float64{4, -1, -1, -1, -1}, 125, 5, 8, 4},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			lap, err := LaplacianSpectrum(tt.g)
			if err != nil {
				t.Fatal(err)
			}
			if !sortedApprox(lap, tt.lap, tol) {
				t.Errorf("Laplacian spectrum = %v, want %v", lap, tt.lap)
			}
			adj, err := AdjacencySpectrum(tt.g)
			if err != nil {
				t.Fatal(err)
			}
			if !sortedApprox(adj, tt.adj, tol) {
				t.Errorf("adjacency spectrum = %v, want %v", adj, tt.adj)
			}
			trees, err := NumberOfSpanningTrees(tt.g)
			if err != nil {
				t.Fatal(err)
			}
			if trees != tt.trees {
				t.Errorf("spanning trees = %d, want %d", trees, tt.trees)
			}
			ac, err := AlgebraicConnectivity(tt.g)
			if err != nil {
				t.Fatal(err)
			}
			if !approx(ac, tt.algConn) {
				t.Errorf("algebraic connectivity = %v, want %v", ac, tt.algConn)
			}
			en, err := GraphEnergy(tt.g)
			if err != nil {
				t.Fatal(err)
			}
			if !approx(en, tt.energy) {
				t.Errorf("graph energy = %v, want %v", en, tt.energy)
			}
			sr, err := SpectralRadiusAdjacency(tt.g)
			if err != nil {
				t.Fatal(err)
			}
			if !approx(sr, tt.specRad) {
				t.Errorf("spectral radius = %v, want %v", sr, tt.specRad)
			}
		})
	}
}

func TestNormalizedLaplacianBounds(t *testing.T) {
	// All normalized Laplacian eigenvalues lie in [0, 2].
	for _, g := range []*Graph{CompleteGraph(5), CycleGraph(6), PathGraph(4), StarGraph(5)} {
		vals, err := NormalizedLaplacianSpectrum(g)
		if err != nil {
			t.Fatal(err)
		}
		for _, x := range vals {
			if x < -tol || x > 2+tol {
				t.Fatalf("normalized eigenvalue %v out of [0,2]", x)
			}
		}
	}
}

// ---------------------------------------------------------------------------
// Structural queries
// ---------------------------------------------------------------------------

func TestConnectivityAndBipartite(t *testing.T) {
	if !CompleteGraph(4).IsConnected() {
		t.Error("K4 should be connected")
	}
	twoTri := twoTriangles()
	if twoTri.IsConnected() {
		t.Error("two triangles should be disconnected")
	}
	if got := twoTri.NumConnectedComponents(); got != 2 {
		t.Errorf("components = %d, want 2", got)
	}
	// spectral component count agrees
	c, err := NumberOfComponentsSpectral(twoTri, 0)
	if err != nil {
		t.Fatal(err)
	}
	if c != 2 {
		t.Errorf("spectral components = %d, want 2", c)
	}
	if ok, _ := CycleGraph(4).IsBipartite(); !ok {
		t.Error("C4 should be bipartite")
	}
	if ok, _ := CycleGraph(5).IsBipartite(); ok {
		t.Error("C5 should not be bipartite")
	}
	if ok, _ := CompleteBipartiteGraph(2, 3).IsBipartite(); !ok {
		t.Error("K2,3 should be bipartite")
	}
}

func TestDegreesAndRegular(t *testing.T) {
	g := CompleteGraph(4)
	if reg, d := g.IsRegular(); !reg || !approx(d, 3) {
		t.Errorf("K4 regular=%v deg=%v, want true 3", reg, d)
	}
	if !approx(g.Density(), 1) {
		t.Errorf("K4 density = %v, want 1", g.Density())
	}
	if StarGraph(4).Density() >= 1 {
		t.Error("star density should be < 1")
	}
}

// ---------------------------------------------------------------------------
// Effective resistance / Kirchhoff
// ---------------------------------------------------------------------------

func TestEffectiveResistance(t *testing.T) {
	p := PathGraph(3) // 0-1-2, two unit resistors in series
	r, err := EffectiveResistance(p, 0, 2)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(r, 2) {
		t.Errorf("R(0,2) on P3 = %v, want 2", r)
	}
	c := CycleGraph(4) // opposite vertices: two length-2 paths in parallel
	r02, _ := EffectiveResistance(c, 0, 2)
	if !approx(r02, 1) {
		t.Errorf("R(0,2) on C4 = %v, want 1", r02)
	}
	kf, err := KirchhoffIndex(c)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(kf, 5) {
		t.Errorf("Kirchhoff index of C4 = %v, want 5", kf)
	}
}

// ---------------------------------------------------------------------------
// Cheeger
// ---------------------------------------------------------------------------

func TestCheegerBounds(t *testing.T) {
	for _, g := range []*Graph{CompleteGraph(5), CycleGraph(6), PathGraph(5)} {
		lo, err := CheegerLowerBound(g)
		if err != nil {
			t.Fatal(err)
		}
		hi, err := CheegerUpperBound(g)
		if err != nil {
			t.Fatal(err)
		}
		h, err := CheegerConstant(g)
		if err != nil {
			t.Fatal(err)
		}
		if h < lo-tol || h > hi+tol {
			t.Errorf("Cheeger constant %v not in [%v, %v]", h, lo, hi)
		}
	}
}

// ---------------------------------------------------------------------------
// Clustering / bisection
// ---------------------------------------------------------------------------

func TestSpectralBisectionPath(t *testing.T) {
	part, err := SpectralBisection(PathGraph(4))
	if err != nil {
		t.Fatal(err)
	}
	if part[0] != part[1] || part[2] != part[3] || part[0] == part[2] {
		t.Errorf("bisection of P4 = %v, want {0,1} vs {2,3}", part)
	}
}

func TestSpectralClusteringTwoTriangles(t *testing.T) {
	g := twoTriangles()
	labels, err := SpectralClustering(g, 2)
	if err != nil {
		t.Fatal(err)
	}
	if labels[0] != labels[1] || labels[1] != labels[2] {
		t.Errorf("triangle A not grouped: %v", labels)
	}
	if labels[3] != labels[4] || labels[4] != labels[5] {
		t.Errorf("triangle B not grouped: %v", labels)
	}
	if labels[0] == labels[3] {
		t.Errorf("triangles not separated: %v", labels)
	}
}

// ---------------------------------------------------------------------------
// Centralities and PageRank
// ---------------------------------------------------------------------------

func TestDegreeCentrality(t *testing.T) {
	dc := DegreeCentrality(CompleteGraph(4))
	for i, x := range dc {
		if !approx(x, 1) {
			t.Errorf("degree centrality[%d] = %v, want 1", i, x)
		}
	}
}

func TestEigenvectorCentralityRegular(t *testing.T) {
	// On a regular graph every vertex has equal eigenvector centrality.
	ec, err := EigenvectorCentrality(CycleGraph(5), 0, 0)
	if err != nil {
		t.Fatal(err)
	}
	want := ec[0]
	for i, x := range ec {
		if !approx(x, want) {
			t.Errorf("ec[%d] = %v, want %v", i, x, want)
		}
		if x <= 0 {
			t.Errorf("ec[%d] = %v, want positive", i, x)
		}
	}
}

func TestKatzCentralityStar(t *testing.T) {
	// On a star the centre should dominate the leaves.
	g := StarGraph(5)
	x, err := KatzCentrality(g, 0.1, 1)
	if err != nil {
		t.Fatal(err)
	}
	for i := 1; i < g.Order(); i++ {
		if x[0] <= x[i] {
			t.Errorf("centre %v not greater than leaf %v", x[0], x[i])
		}
	}
}

func TestPageRankUniform(t *testing.T) {
	// PageRank on a vertex-transitive graph is uniform.
	pr, err := PageRank(CompleteGraph(4), 0.85, 0, 0)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(VecSum(pr), 1) {
		t.Errorf("PageRank sum = %v, want 1", VecSum(pr))
	}
	for i, x := range pr {
		if !approx(x, 0.25) {
			t.Errorf("pr[%d] = %v, want 0.25", i, x)
		}
	}
}

func TestPageRankStar(t *testing.T) {
	// On a star the hub gets the most PageRank.
	pr, err := PageRank(StarGraph(5), 0.85, 0, 0)
	if err != nil {
		t.Fatal(err)
	}
	for i := 1; i < len(pr); i++ {
		if pr[0] <= pr[i] {
			t.Errorf("hub %v not greater than leaf %v", pr[0], pr[i])
		}
	}
}

func TestClosenessPath(t *testing.T) {
	// On P3, the middle vertex is most central.
	cc := ClosenessCentrality(PathGraph(3))
	if cc[1] <= cc[0] || cc[1] <= cc[2] {
		t.Errorf("closeness = %v, middle should be largest", cc)
	}
}

func TestShortestPaths(t *testing.T) {
	d := DijkstraDistances(PathGraph(4), 0)
	if !VecApproxEqual(d, []float64{0, 1, 2, 3}, tol) {
		t.Errorf("distances = %v, want [0 1 2 3]", d)
	}
	if !approx(Diameter(PathGraph(4)), 3) {
		t.Errorf("diameter of P4 = %v, want 3", Diameter(PathGraph(4)))
	}
}

// ---------------------------------------------------------------------------
// Matrix identities
// ---------------------------------------------------------------------------

func TestIncidenceLaplacian(t *testing.T) {
	// B·Bᵀ = L for the unweighted Laplacian.
	g := CycleGraph(5)
	b := g.IncidenceMatrix()
	bt := b.Transpose()
	bbt, _ := b.Mul(bt)
	if !bbt.ApproxEqual(g.Laplacian(), tol) {
		t.Errorf("B·Bᵀ != L")
	}
}

func TestLaplacianRowSumsZero(t *testing.T) {
	for _, g := range []*Graph{CompleteGraph(4), CycleGraph(5), StarGraph(6)} {
		for _, s := range g.Laplacian().RowSums() {
			if math.Abs(s) > tol {
				t.Errorf("Laplacian row sum = %v, want 0", s)
			}
		}
	}
}

// twoTriangles builds two disjoint triangles on 6 vertices.
func twoTriangles() *Graph {
	g := NewGraph(6)
	g.AddEdge(0, 1)
	g.AddEdge(1, 2)
	g.AddEdge(0, 2)
	g.AddEdge(3, 4)
	g.AddEdge(4, 5)
	g.AddEdge(3, 5)
	return g
}

// ExampleNumberOfSpanningTrees demonstrates the Matrix-Tree theorem: the
// complete graph K4 has 4^(4-2) = 16 spanning trees by Cayley's formula.
func ExampleNumberOfSpanningTrees() {
	g := CompleteGraph(4)
	trees, _ := NumberOfSpanningTrees(g)
	ac, _ := AlgebraicConnectivity(g)
	fmt.Printf("spanning trees: %d\n", trees)
	fmt.Printf("algebraic connectivity: %.1f\n", ac)
	// Output:
	// spanning trees: 16
	// algebraic connectivity: 4.0
}
