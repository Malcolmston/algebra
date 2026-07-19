package simplicial

import (
	"fmt"
	"math"
	"math/big"
	"testing"
)

func eqIntSlice(a, b []int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func bigsToInt64(bs []*big.Int) []int64 {
	out := make([]int64, len(bs))
	for i, b := range bs {
		out[i] = b.Int64()
	}
	return out
}

func eqInt64Slice(a, b []int64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func TestSimplexBasics(t *testing.T) {
	s := NewSimplex(2, 0, 1, 2) // duplicate 2 collapses
	if s.Dim() != 2 {
		t.Fatalf("dim = %d, want 2", s.Dim())
	}
	if !eqIntSlice(s.Vertices(), []int{0, 1, 2}) {
		t.Fatalf("vertices = %v", s.Vertices())
	}
	if s.Key() != "0,1,2" {
		t.Fatalf("key = %q", s.Key())
	}
	if got := len(s.Faces()); got != 3 {
		t.Fatalf("faces = %d, want 3", got)
	}
	if got := len(s.Closure()); got != 7 { // 2^3-1
		t.Fatalf("closure = %d, want 7", got)
	}
	if !Edge(1, 3).IsFaceOf(NewSimplex(1, 2, 3)) {
		t.Fatalf("edge should be a face")
	}
	if NewSimplex(0, 4).IsFaceOf(NewSimplex(0, 1, 2)) {
		t.Fatalf("{0,4} is not a face of {0,1,2}")
	}
}

func TestSimplexBoundarySigns(t *testing.T) {
	b := Triangle(0, 1, 2).Boundary()
	// ∂[0,1,2] = [1,2] - [0,2] + [0,1]
	want := []struct {
		sign int
		key  string
	}{
		{+1, "1,2"},
		{-1, "0,2"},
		{+1, "0,1"},
	}
	if len(b) != 3 {
		t.Fatalf("boundary len = %d", len(b))
	}
	for i, term := range b {
		if term.Sign != want[i].sign || term.Face.Key() != want[i].key {
			t.Fatalf("term %d = (%d,%s), want (%d,%s)", i, term.Sign, term.Face.Key(), want[i].sign, want[i].key)
		}
	}
}

func TestComplexFVectorAndEuler(t *testing.T) {
	tests := []struct {
		name    string
		c       *Complex
		fvec    []int
		euler   int
		dim     int
		numVert int
	}{
		{"triangle-filled", StandardSimplex(2), []int{3, 3, 1}, 1, 2, 3},
		{"circle", CycleComplex(5), []int{5, 5}, 0, 1, 5},
		{"S1-as-triangle", SphereComplex(1), []int{3, 3}, 0, 1, 3},
		{"S2", SphereComplex(2), []int{4, 6, 4}, 2, 2, 4},
		{"torus", TorusComplex(), []int{7, 21, 14}, 0, 2, 7},
		{"RP2", ProjectivePlaneComplex(), []int{6, 15, 10}, 1, 2, 6},
		{"3points", DiscreteComplex(3), []int{3}, 3, 0, 3},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := tt.c.FVector(); !eqIntSlice(got, tt.fvec) {
				t.Errorf("fvector = %v, want %v", got, tt.fvec)
			}
			if got := tt.c.EulerCharacteristic(); got != tt.euler {
				t.Errorf("euler = %d, want %d", got, tt.euler)
			}
			if got := tt.c.Dimension(); got != tt.dim {
				t.Errorf("dim = %d, want %d", got, tt.dim)
			}
			if got := tt.c.NumVertices(); got != tt.numVert {
				t.Errorf("numVert = %d, want %d", got, tt.numVert)
			}
			if !tt.c.IsValid() {
				t.Errorf("complex is not downward closed")
			}
			// Euler-Poincaré: chi from f-vector equals chi from Betti numbers.
			if a, b := tt.c.EulerCharacteristic(), tt.c.EulerCharacteristicFromBetti(); a != b {
				t.Errorf("euler mismatch: f=%d betti=%d", a, b)
			}
		})
	}
}

func TestBettiNumbers(t *testing.T) {
	tests := []struct {
		name   string
		c      *Complex
		bettiQ []int
		betti2 []int
	}{
		{"point", PointComplex(), []int{1}, []int{1}},
		{"3points", DiscreteComplex(3), []int{3}, []int{3}},
		{"path", PathComplex(4), []int{1, 0}, []int{1, 0}},
		{"circle", CycleComplex(6), []int{1, 1}, []int{1, 1}},
		{"S1", SphereComplex(1), []int{1, 1}, []int{1, 1}},
		{"S2", SphereComplex(2), []int{1, 0, 1}, []int{1, 0, 1}},
		{"S3", SphereComplex(3), []int{1, 0, 0, 1}, []int{1, 0, 0, 1}},
		{"disk", StandardSimplex(2), []int{1, 0, 0}, []int{1, 0, 0}},
		{"torus", TorusComplex(), []int{1, 2, 1}, []int{1, 2, 1}},
		{"RP2", ProjectivePlaneComplex(), []int{1, 0, 0}, []int{1, 1, 1}},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := tt.c.BettiNumbers(); !eqIntSlice(got, tt.bettiQ) {
				t.Errorf("betti Q = %v, want %v", got, tt.bettiQ)
			}
			if got := tt.c.BettiNumbersGF2(); !eqIntSlice(got, tt.betti2) {
				t.Errorf("betti GF2 = %v, want %v", got, tt.betti2)
			}
			// b_0 must equal the number of connected components.
			if got, want := tt.c.BettiNumber(0), tt.c.NumConnectedComponents(); got != want {
				t.Errorf("b0=%d but components=%d", got, want)
			}
		})
	}
}

func TestTorsionAndHomologyGroups(t *testing.T) {
	rp2 := ProjectivePlaneComplex()
	tors := rp2.TorsionCoefficients(1)
	if !eqInt64Slice(bigsToInt64(tors), []int64{2}) {
		t.Fatalf("RP2 H1 torsion = %v, want [2]", bigsToInt64(tors))
	}
	h1 := rp2.HomologyZ(1)
	if h1.FreeRank != 0 || h1.String() != "Z/2" {
		t.Fatalf("RP2 H1 = %s (rank %d), want Z/2", h1.String(), h1.FreeRank)
	}
	if s := rp2.HomologyZ(0).String(); s != "Z" {
		t.Fatalf("RP2 H0 = %s, want Z", s)
	}
	// Torus is torsion free with H1 = Z^2.
	tor := TorusComplex()
	if len(tor.TorsionCoefficients(1)) != 0 {
		t.Fatalf("torus should be torsion free")
	}
	if s := tor.HomologyZ(1).String(); s != "Z^2" {
		t.Fatalf("torus H1 = %s, want Z^2", s)
	}
	if s := tor.HomologyZ(2).String(); s != "Z" {
		t.Fatalf("torus H2 = %s, want Z", s)
	}
}

func TestBoundaryOfBoundaryIsZero(t *testing.T) {
	for _, c := range []*Complex{SphereComplex(2), TorusComplex(), ProjectivePlaneComplex(), StandardSimplex(3)} {
		for k := 2; k <= c.Dimension(); k++ {
			d1 := c.BoundaryMatrixZ(k)
			d2 := c.BoundaryMatrixZ(k - 1)
			prod, err := d2.Mul(d1) // ∂_{k-1} ∘ ∂_k
			if err != nil {
				t.Fatalf("mul error: %v", err)
			}
			if !prod.IsZero() {
				t.Fatalf("∂∂ != 0 at k=%d", k)
			}
		}
	}
}

func TestSmithNormalForm(t *testing.T) {
	tests := []struct {
		name string
		r, c int
		vals []int64
		diag []int64 // expected non-zero invariant factors
		rank int
	}{
		{"diag26", 2, 2, []int64{2, 0, 0, 6}, []int64{2, 6}, 2},
		{"rank-def", 2, 2, []int64{1, 2, 2, 4}, []int64{1}, 1},
		{"classic", 3, 3, []int64{2, 4, 4, -6, 6, 12, 10, -4, -16}, []int64{2, 6, 12}, 3},
		{"identity", 3, 3, []int64{1, 0, 0, 0, 1, 0, 0, 0, 1}, []int64{1, 1, 1}, 3},
		{"zero", 2, 3, []int64{0, 0, 0, 0, 0, 0}, nil, 0},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			m := IntMatrixFromInt64(tt.r, tt.c, tt.vals)
			snf := m.SmithNormalForm()
			if got := bigsToInt64(snf.InvariantFactors()); !eqInt64Slice(got, tt.diag) {
				t.Errorf("invariant factors = %v, want %v", got, tt.diag)
			}
			if snf.Rank() != tt.rank {
				t.Errorf("rank = %d, want %d", snf.Rank(), tt.rank)
			}
			// verify U*A*V == D
			ua, _ := snf.U.Mul(m)
			uav, _ := ua.Mul(snf.V)
			if !uav.Equal(snf.D) {
				t.Errorf("U*A*V != D\nUAV=\n%sD=\n%s", uav, snf.D)
			}
			// U and V must be unimodular (det = ±1)
			for _, tr := range []*IntMatrix{snf.U, snf.V} {
				det, err := tr.Determinant()
				if err != nil {
					t.Fatalf("det error: %v", err)
				}
				if a := new(big.Int).Abs(det); a.Cmp(big.NewInt(1)) != 0 {
					t.Errorf("transform not unimodular, det = %s", det)
				}
			}
			// divisibility chain d_i | d_{i+1}
			inv := snf.InvariantFactors()
			for i := 0; i+1 < len(inv); i++ {
				rem := new(big.Int).Rem(inv[i+1], inv[i])
				if rem.Sign() != 0 {
					t.Errorf("divisibility fails: %s does not divide %s", inv[i], inv[i+1])
				}
			}
		})
	}
}

func TestGF2Matrix(t *testing.T) {
	m := NewGF2Matrix(3, 3)
	vals := [][]int{{1, 0, 1}, {0, 1, 1}, {1, 1, 0}}
	for i := range vals {
		for j := range vals[i] {
			m.Set(i, j, vals[i][j])
		}
	}
	if m.Rank() != 2 { // rows sum to zero over GF2
		t.Fatalf("rank = %d, want 2", m.Rank())
	}
	if m.Nullity() != 1 {
		t.Fatalf("nullity = %d, want 1", m.Nullity())
	}
	kb := m.KernelBasis()
	if len(kb) != 1 {
		t.Fatalf("kernel basis size = %d, want 1", len(kb))
	}
	// verify m * kernelvector = 0
	for _, v := range kb {
		for i := 0; i < m.Rows(); i++ {
			s := 0
			for j := 0; j < m.Cols(); j++ {
				s ^= m.At(i, j) & v[j]
			}
			if s != 0 {
				t.Fatalf("kernel vector not in kernel")
			}
		}
	}
}

func TestRatMatrixDeterminantAndRank(t *testing.T) {
	m := NewRatMatrix(2, 2)
	m.SetInt(0, 0, 1)
	m.SetInt(0, 1, 2)
	m.SetInt(1, 0, 3)
	m.SetInt(1, 1, 4)
	det, err := m.Determinant()
	if err != nil {
		t.Fatal(err)
	}
	if det.Cmp(big.NewRat(-2, 1)) != 0 {
		t.Fatalf("det = %s, want -2", det.RatString())
	}
	if m.Rank() != 2 {
		t.Fatalf("rank = %d, want 2", m.Rank())
	}
	// singular matrix
	s := NewRatMatrix(2, 2)
	s.SetInt(0, 0, 1)
	s.SetInt(0, 1, 2)
	s.SetInt(1, 0, 2)
	s.SetInt(1, 1, 4)
	if d, _ := s.Determinant(); d.Sign() != 0 {
		t.Fatalf("singular det should be 0, got %s", d.RatString())
	}
	if s.Rank() != 1 {
		t.Fatalf("rank = %d, want 1", s.Rank())
	}
}

func TestMinimalEnclosingBall(t *testing.T) {
	tests := []struct {
		name   string
		pts    [][]float64
		radius float64
	}{
		{"two-points", [][]float64{{0, 0}, {2, 0}}, 1},
		{"unit-square", [][]float64{{0, 0}, {1, 0}, {1, 1}, {0, 1}}, math.Sqrt2 / 2},
		{"triangle", [][]float64{{0, 0}, {4, 0}, {2, 0}}, 2},
		{"single", [][]float64{{3, 7}}, 0},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			c, r := MinimalEnclosingBall(tt.pts)
			if math.Abs(r-tt.radius) > 1e-9 {
				t.Errorf("radius = %v, want %v", r, tt.radius)
			}
			// every point must lie in the ball
			for _, p := range tt.pts {
				if EuclideanDistance(p, c) > r+1e-9 {
					t.Errorf("point %v outside ball (c=%v r=%v)", p, c, r)
				}
			}
		})
	}
}

func TestVietorisRipsCircle(t *testing.T) {
	// 6 points on a circle of radius 1; nearest-neighbour spacing = 1 (hexagon).
	pc := CirclePoints(6, 1)
	// epsilon just above the edge length (1.0) but below the diagonal: recovers
	// a hexagonal circle -> b0=1, b1=1.
	vr := VietorisRips(pc, 1.01, 2)
	if b := vr.BettiNumbers(); !eqIntSlice(b, []int{1, 1}) {
		t.Fatalf("Rips hexagon betti = %v, want [1 1]", b)
	}
	// large epsilon fills everything into a single blob (contractible).
	vrBig := VietorisRips(pc, 10, 2)
	if b0 := vrBig.BettiNumber(0); b0 != 1 {
		t.Fatalf("b0 = %d, want 1", b0)
	}
	if b1 := vrBig.BettiNumber(1); b1 != 0 {
		t.Fatalf("b1 = %d, want 0 for large epsilon", b1)
	}
}

func TestCechVsRips(t *testing.T) {
	pc := NewPointCloud([]float64{0, 0}, []float64{2, 0}, []float64{1, 1.7})
	// At r=1 the three balls pairwise touch (edges present) but do not share a
	// common point, so no triangle: b1 = 1 (a loop).
	cech := Cech(pc, 1.0, 2)
	if cech.NumSimplicesOfDim(1) != 3 {
		t.Fatalf("cech edges = %d, want 3", cech.NumSimplicesOfDim(1))
	}
	if cech.NumSimplicesOfDim(2) != 0 {
		t.Fatalf("cech triangles = %d, want 0", cech.NumSimplicesOfDim(2))
	}
	if b := cech.BettiNumbers(); !eqIntSlice(b, []int{1, 1}) {
		t.Fatalf("cech betti = %v, want [1 1]", b)
	}
	// Rips at epsilon = 2 (max pairwise distance) fills the triangle.
	rips := VietorisRips(pc, 2.0, 2)
	if rips.NumSimplicesOfDim(2) != 1 {
		t.Fatalf("rips triangles = %d, want 1", rips.NumSimplicesOfDim(2))
	}
}

func TestPersistentHomologyCircle(t *testing.T) {
	pc := CirclePoints(8, 1)
	f := RipsFiltration(pc, 3.0, 2)
	p := PersistentHomology(f)

	// Exactly one essential H0 class (connected) and one essential-or-long H1
	// class (the circle).
	ess0 := 0
	for _, pr := range p.PairsOfDim(0) {
		if pr.IsEssential() {
			ess0++
		}
	}
	if ess0 != 1 {
		t.Fatalf("essential H0 classes = %d, want 1", ess0)
	}
	// There must be a persistent H1 feature (the hole) born before it dies.
	var best float64
	for _, pr := range p.PairsOfDim(1) {
		if l := pr.Persistence(); l > best {
			best = l
		}
	}
	if best <= 0 {
		t.Fatalf("expected a positive-persistence H1 feature, got %v", best)
	}
	// number of H0 pairs (finite + essential) equals number of points
	if len(p.PairsOfDim(0)) != pc.Len() {
		t.Fatalf("H0 pairs = %d, want %d", len(p.PairsOfDim(0)), pc.Len())
	}
	// Betti at t=0 (only vertices present) is the number of points.
	if b := p.BettiAt(0, 0); b != pc.Len() {
		t.Fatalf("BettiAt(0,0) = %d, want %d", b, pc.Len())
	}
}

func TestPersistenceMatchesStaticHomology(t *testing.T) {
	// The essential classes of a Rips filtration at threshold T must reproduce
	// the Betti numbers of the Rips complex at that scale.
	pc := CirclePoints(6, 1)
	T := 1.01
	f := RipsFiltration(pc, T, 2)
	p := PersistentHomology(f)
	static := VietorisRips(pc, T, 2)
	for k := 0; k <= 1; k++ {
		if got, want := p.BettiAt(k, T-1e-6), static.BettiNumber(k); got != want {
			t.Errorf("dim %d: persistence betti %d != static %d", k, got, want)
		}
	}
}

func TestConeAndSuspension(t *testing.T) {
	// Cone on any complex is contractible.
	cone := Cone(CycleComplex(5), 99)
	if b := cone.BettiNumbers(); !eqIntSlice(b, []int{1, 0, 0}) {
		t.Fatalf("cone betti = %v, want [1 0 0]", b)
	}
	// Suspension of S^1 is S^2.
	susp := Suspension(SphereComplex(1), 100, 101)
	if b := susp.BettiNumbers(); !eqIntSlice(b, []int{1, 0, 1}) {
		t.Fatalf("suspension(S1) betti = %v, want [1 0 1] (S2)", b)
	}
}

func TestLinkAndStar(t *testing.T) {
	// In S^2 (boundary of tetrahedron), the link of any vertex is a triangle
	// (a circle): b0=1, b1=1.
	s2 := SphereComplex(2)
	link := s2.Link(Vertex(0))
	if b := link.BettiNumbers(); !eqIntSlice(b, []int{1, 1}) {
		t.Fatalf("link betti = %v, want [1 1]", b)
	}
	star := s2.Star(Vertex(0))
	if len(star) == 0 {
		t.Fatalf("star should be non-empty")
	}
}

func TestDisjointUnionComponents(t *testing.T) {
	u := DisjointUnion(SphereComplex(1), SphereComplex(1), 100)
	if got := u.NumConnectedComponents(); got != 2 {
		t.Fatalf("components = %d, want 2", got)
	}
	if b := u.BettiNumbers(); !eqIntSlice(b, []int{2, 2}) {
		t.Fatalf("betti = %v, want [2 2]", b)
	}
}

func ExampleComplex_BettiNumbers() {
	// The torus T^2 has Betti numbers 1, 2, 1.
	t := TorusComplex()
	fmt.Println(t.BettiNumbers())
	fmt.Println("Euler:", t.EulerCharacteristic())
	// Output:
	// [1 2 1]
	// Euler: 0
}

func ExampleComplex_HomologyZ() {
	// The real projective plane has 2-torsion in H_1.
	rp2 := ProjectivePlaneComplex()
	fmt.Println("H0 =", rp2.HomologyZ(0))
	fmt.Println("H1 =", rp2.HomologyZ(1))
	fmt.Println("H2 =", rp2.HomologyZ(2))
	// Output:
	// H0 = Z
	// H1 = Z/2
	// H2 = 0
}

func ExampleVietorisRips() {
	// Six points on a circle, joined when closer than 1.01, recover a loop.
	pc := CirclePoints(6, 1)
	vr := VietorisRips(pc, 1.01, 2)
	fmt.Println(vr.BettiNumbers())
	// Output:
	// [1 1]
}
