package designs

import (
	"fmt"
	"math/rand"
	"testing"
)

func TestNumberTheory(t *testing.T) {
	if g := Gcd(24, 36); g != 12 {
		t.Errorf("Gcd(24,36)=%d want 12", g)
	}
	if l := Lcm(4, 6); l != 12 {
		t.Errorf("Lcm(4,6)=%d want 12", l)
	}
	g, x, y := ExtendedGcd(240, 46)
	if g != 2 || 240*x+46*y != g {
		t.Errorf("ExtendedGcd wrong: g=%d x=%d y=%d", g, x, y)
	}
	inv, err := ModInverse(3, 11)
	if err != nil || inv != 4 {
		t.Errorf("ModInverse(3,11)=%d,%v want 4", inv, err)
	}
	if v := ModExp(2, 10, 1000); v != 24 {
		t.Errorf("ModExp(2,10,1000)=%d want 24", v)
	}
	if EulerPhi(36) != 12 {
		t.Errorf("EulerPhi(36) wrong")
	}
	if !IsPrimePower(27) || IsPrimePower(12) {
		t.Errorf("IsPrimePower wrong")
	}
	if got := Binomial(10, 3); got != 120 {
		t.Errorf("Binomial(10,3)=%d want 120", got)
	}
	pr, err := PrimitiveRoot(7)
	if err != nil || pr != 3 {
		t.Errorf("PrimitiveRoot(7)=%d want 3", pr)
	}
}

func TestLegendreAndResidues(t *testing.T) {
	tests := []struct {
		a, p, want int
	}{
		{2, 7, 1}, {3, 7, -1}, {7, 7, 0}, {1, 11, 1}, {2, 11, -1},
	}
	for _, tt := range tests {
		got, err := LegendreSymbol(tt.a, tt.p)
		if err != nil || got != tt.want {
			t.Errorf("LegendreSymbol(%d,%d)=%d,%v want %d", tt.a, tt.p, got, err, tt.want)
		}
	}
	qr := QuadraticResidues(7)
	if len(qr) != 3 || qr[0] != 1 || qr[1] != 2 || qr[2] != 4 {
		t.Errorf("QuadraticResidues(7)=%v want [1 2 4]", qr)
	}
}

func TestGaloisField(t *testing.T) {
	for _, q := range []int{2, 3, 4, 5, 7, 8, 9, 16, 25, 27, 32} {
		f, err := NewGaloisField(q)
		if err != nil {
			t.Fatalf("NewGaloisField(%d): %v", q, err)
		}
		if f.Order() != q {
			t.Errorf("Order=%d want %d", f.Order(), q)
		}
		// Field axioms across a sample of elements.
		for a := 0; a < q; a++ {
			if f.Add(a, f.Neg(a)) != 0 {
				t.Fatalf("q=%d additive inverse of %d", q, a)
			}
			if a != 0 {
				inv, err := f.Inv(a)
				if err != nil || f.Mul(a, inv) != 1 {
					t.Fatalf("q=%d inverse of %d", q, a)
				}
			}
		}
		// Distributivity.
		for _, tr := range [][3]int{{1, 2, 3}, {2, 3, 4}} {
			a, b, c := tr[0]%q, tr[1]%q, tr[2]%q
			if f.Mul(a, f.Add(b, c)) != f.Add(f.Mul(a, b), f.Mul(a, c)) {
				t.Fatalf("q=%d distributivity", q)
			}
		}
		p, err := f.PrimitiveElement()
		if err != nil || !f.IsPrimitiveElement(p) {
			t.Fatalf("q=%d primitive element", q)
		}
		ord, _ := f.ElementOrder(p)
		if ord != q-1 {
			t.Fatalf("q=%d primitive order=%d want %d", q, ord, q-1)
		}
	}
}

func TestBIBDParams(t *testing.T) {
	tests := []struct {
		name         string
		v, k, lambda int
		wantR, wantB int
		symmetric    bool
		fisher       bool
	}{
		{"Fano", 7, 3, 1, 3, 7, true, true},
		{"AG23", 9, 3, 1, 4, 12, false, true},
		{"PG(2,3)", 13, 4, 1, 4, 13, true, true},
		{"biplane", 11, 5, 2, 5, 11, true, true},
	}
	for _, tt := range tests {
		p, err := NewBIBDParams(tt.v, tt.k, tt.lambda)
		if err != nil {
			t.Fatalf("%s: %v", tt.name, err)
		}
		if p.R != tt.wantR || p.B != tt.wantB {
			t.Errorf("%s: r=%d b=%d want r=%d b=%d", tt.name, p.R, p.B, tt.wantR, tt.wantB)
		}
		if !p.Valid() {
			t.Errorf("%s: not valid", tt.name)
		}
		if p.IsSymmetric() != tt.symmetric {
			t.Errorf("%s: symmetric=%v want %v", tt.name, p.IsSymmetric(), tt.symmetric)
		}
		if p.FisherInequalityHolds() != tt.fisher {
			t.Errorf("%s: fisher mismatch", tt.name)
		}
	}
	// Invalid parameter combination (no 2-(8,3,1) design: 3 does not divide 8*7/... )
	if _, err := NewBIBDParams(8, 3, 1); err == nil {
		t.Errorf("expected error for non-integral (8,3,1)")
	}
}

func TestSymmetricComplementDerivedResidual(t *testing.T) {
	p, _ := NewBIBDParams(7, 3, 1) // Fano, symmetric
	c, err := p.Complement()
	if err != nil {
		t.Fatalf("complement: %v", err)
	}
	if c.V != 7 || c.K != 4 || c.Lambda != 2 {
		t.Errorf("complement of Fano = 2-(%d,%d,%d) want 2-(7,4,2)", c.V, c.K, c.Lambda)
	}
	// PG(2,3) is 2-(13,4,1) symmetric.
	pg, _ := NewBIBDParams(13, 4, 1)
	res, err := pg.Residual()
	if err != nil {
		t.Fatalf("residual: %v", err)
	}
	if res.V != 9 || res.K != 3 || res.Lambda != 1 {
		t.Errorf("residual of PG(2,3) = 2-(%d,%d,%d) want 2-(9,3,1)", res.V, res.K, res.Lambda)
	}
}

func TestLatinAndMOLS(t *testing.T) {
	l, err := CyclicLatinSquare(5)
	if err != nil || !l.IsLatin() {
		t.Fatalf("cyclic latin square 5")
	}
	if l.At(2, 3) != 0 {
		t.Errorf("cyclic L[2][3]=%d want 0", l.At(2, 3))
	}
	for _, q := range []int{3, 4, 5, 7, 8, 9} {
		squares, err := MOLS(q)
		if err != nil {
			t.Fatalf("MOLS(%d): %v", q, err)
		}
		if len(squares) != q-1 {
			t.Errorf("MOLS(%d) count=%d want %d", q, len(squares), q-1)
		}
		if !IsMOLS(squares) {
			t.Errorf("MOLS(%d) not mutually orthogonal", q)
		}
	}
	// Order 6 has no pair of orthogonal Latin squares (Euler), but two cyclic
	// squares of order 5 are orthogonal.
	a, _ := MOLS(5)
	ok, err := AreOrthogonal(a[0], a[1])
	if err != nil || !ok {
		t.Errorf("expected orthogonal pair of order 5")
	}
}

func TestHadamard(t *testing.T) {
	for k := 0; k <= 5; k++ {
		h, err := SylvesterHadamard(k)
		if err != nil {
			t.Fatalf("Sylvester(%d): %v", k, err)
		}
		if h.Order() != 1<<k {
			t.Errorf("order=%d want %d", h.Order(), 1<<k)
		}
		if !h.IsHadamard() {
			t.Errorf("Sylvester(%d) not Hadamard", k)
		}
	}
	// Paley I on q=7,11 -> orders 8,12.
	for _, q := range []int{3, 7, 11, 19, 23} {
		h, err := PaleyConstructionI(q)
		if err != nil {
			t.Fatalf("PaleyI(%d): %v", q, err)
		}
		if h.Order() != q+1 || !h.IsHadamard() {
			t.Errorf("PaleyI(%d) order=%d not Hadamard", q, h.Order())
		}
	}
	// Paley II on q=5,13 -> orders 12,28.
	for _, q := range []int{5, 13, 17} {
		h, err := PaleyConstructionII(q)
		if err != nil {
			t.Fatalf("PaleyII(%d): %v", q, err)
		}
		if h.Order() != 2*(q+1) || !h.IsHadamard() {
			t.Errorf("PaleyII(%d) order=%d not Hadamard", q, h.Order())
		}
	}
	// Hadamard design from order 8 -> 2-(7,3,1) (Fano).
	h, _ := SylvesterHadamard(3)
	d, err := HadamardToDesign(h)
	if err != nil {
		t.Fatalf("HadamardToDesign: %v", err)
	}
	p, ok := d.IsBIBD()
	if !ok || p.V != 7 || p.K != 3 || p.Lambda != 1 {
		t.Errorf("Hadamard design = %v ok=%v want 2-(7,3,1)", p, ok)
	}
}

func TestDifferenceSets(t *testing.T) {
	// Paley difference set of QRs mod 7 -> (7,3,1).
	d, err := QuadraticResidueDifferenceSet(7)
	if err != nil {
		t.Fatalf("QR diffset: %v", err)
	}
	n, k, lambda, err := d.Parameters()
	if err != nil || n != 7 || k != 3 || lambda != 1 {
		t.Errorf("QR(7) diffset = (%d,%d,%d) want (7,3,1)", n, k, lambda)
	}
	if !d.IsPlanar() {
		t.Errorf("QR(7) should be planar")
	}
	// Develop into symmetric 2-(7,3,1) design.
	dd := d.Develop()
	bp, ok := dd.IsBIBD()
	if !ok || bp.V != 7 || bp.K != 3 || bp.Lambda != 1 {
		t.Errorf("developed design = %v want 2-(7,3,1)", bp)
	}
	// Singer planar difference set of order 2,3,4 -> (7,3,1),(13,4,1),(21,5,1).
	for _, q := range []int{2, 3, 4, 5} {
		s, err := SingerDifferenceSet(q)
		if err != nil {
			t.Fatalf("Singer(%d): %v", q, err)
		}
		n, k, lambda, err := s.Parameters()
		wantN := q*q + q + 1
		if err != nil || n != wantN || k != q+1 || lambda != 1 {
			t.Errorf("Singer(%d) = (%d,%d,%d) want (%d,%d,1)", q, n, k, lambda, wantN, q+1)
		}
	}
}

func TestProjectivePlane(t *testing.T) {
	for _, q := range []int{2, 3, 4, 5, 7, 8, 9} {
		pl, err := NewProjectivePlane(q)
		if err != nil {
			t.Fatalf("PG(2,%d): %v", q, err)
		}
		if pl.NumPoints() != q*q+q+1 || pl.NumLines() != q*q+q+1 {
			t.Errorf("PG(2,%d) size wrong", q)
		}
		if !pl.IsProjectivePlane() {
			t.Errorf("PG(2,%d) axioms fail", q)
		}
		// Every line has q+1 points.
		for l := 0; l < pl.NumLines(); l++ {
			if len(pl.PointsOnLine(l)) != q+1 {
				t.Fatalf("PG(2,%d) line %d has %d points", q, l, len(pl.PointsOnLine(l)))
			}
		}
		// Join then meet consistency.
		li, err := pl.LineThroughPoints(0, 1)
		if err != nil {
			t.Fatalf("join: %v", err)
		}
		if !pl.IsIncident(0, li) || !pl.IsIncident(1, li) {
			t.Errorf("join line not incident to its points")
		}
	}
}

func TestAffinePlane(t *testing.T) {
	for _, q := range []int{2, 3, 4, 5} {
		ap, err := NewAffinePlane(q)
		if err != nil {
			t.Fatalf("AG(2,%d): %v", q, err)
		}
		if ap.NumPoints() != q*q || ap.NumLines() != q*q+q {
			t.Errorf("AG(2,%d) size wrong: %d points %d lines", q, ap.NumPoints(), ap.NumLines())
		}
		if !ap.IsAffinePlane() {
			t.Errorf("AG(2,%d) axioms fail", q)
		}
		if ap.NumParallelClasses() != q+1 {
			t.Errorf("AG(2,%d) parallel classes=%d want %d", q, ap.NumParallelClasses(), q+1)
		}
	}
}

func TestSteinerTriple(t *testing.T) {
	// Bose construction for n = 3 (mod 6).
	for _, n := range []int{3, 9, 15, 21, 27} {
		d, err := BoseSteinerTripleSystem(n)
		if err != nil {
			t.Fatalf("Bose(%d): %v", n, err)
		}
		if !IsSteinerTripleSystem(d) {
			t.Errorf("Bose(%d) is not an STS", n)
		}
		want := n * (n - 1) / 6
		if d.NumBlocks() != want {
			t.Errorf("Bose(%d) blocks=%d want %d", n, d.NumBlocks(), want)
		}
	}
	// Hill-climbing for all admissible orders up to 25.
	rng := rand.New(rand.NewSource(12345))
	for _, n := range []int{7, 9, 13, 15, 19, 21, 25} {
		d, err := SteinerTripleSystem(n, rng)
		if err != nil {
			t.Fatalf("STS(%d): %v", n, err)
		}
		if !IsSteinerTripleSystem(d) {
			t.Errorf("STS(%d) hill-climb invalid", n)
		}
	}
	if _, err := SteinerTripleSystem(8, rng); err == nil {
		t.Errorf("STS(8) should be inadmissible")
	}
}

func TestSteinerQuadruple(t *testing.T) {
	for _, k := range []int{2, 3, 4} {
		d, err := BooleanQuadrupleSystem(k)
		if err != nil {
			t.Fatalf("Boolean SQS(2^%d): %v", k, err)
		}
		if !IsSteinerQuadrupleSystem(d) {
			t.Errorf("Boolean SQS(2^%d) invalid", k)
		}
		v := 1 << k
		want := v * (v - 1) * (v - 2) / 24
		if d.NumBlocks() != want {
			t.Errorf("SQS(%d) blocks=%d want %d", v, d.NumBlocks(), want)
		}
	}
	if _, err := SteinerQuadrupleParameters(10); err != nil {
		t.Errorf("SQS params 10 should be admissible")
	}
}

func TestResolvableAndOneFactorization(t *testing.T) {
	// Affine plane AG(2,3) gives a resolvable 2-(9,3,1) design (Kirkman STS(9)).
	r, err := AffineResolution(3)
	if err != nil {
		t.Fatalf("AffineResolution(3): %v", err)
	}
	if r.NumClasses() != 4 {
		t.Errorf("AG(2,3) resolution classes=%d want 4", r.NumClasses())
	}
	if !r.IsValid() {
		t.Errorf("AG(2,3) resolution invalid")
	}
	// One-factorization of K_8.
	for _, m := range []int{2, 4, 6, 8, 10} {
		f, err := OneFactorization(m)
		if err != nil {
			t.Fatalf("OneFactorization(%d): %v", m, err)
		}
		if !IsOneFactorization(m, f) {
			t.Errorf("OneFactorization(%d) invalid", m)
		}
	}
	// Round robin for 5 teams: 5 rounds, each pair meets once.
	sched, err := RoundRobinSchedule(5)
	if err != nil {
		t.Fatalf("RoundRobin(5): %v", err)
	}
	met := make(map[[2]int]int)
	for _, round := range sched {
		for _, e := range round {
			met[edgeSorted(e[0], e[1])]++
		}
	}
	if len(met) != 10 {
		t.Errorf("RoundRobin(5) distinct games=%d want 10", len(met))
	}
	for pair, c := range met {
		if c != 1 {
			t.Errorf("pair %v met %d times", pair, c)
		}
	}
}

func TestDesignRecovery(t *testing.T) {
	// The unique 2-(7,3,1) design given explicitly (Fano).
	blocks := [][]int{
		{0, 1, 2}, {0, 3, 4}, {0, 5, 6}, {1, 3, 5},
		{1, 4, 6}, {2, 3, 6}, {2, 4, 5},
	}
	d, err := NewDesign(7, blocks)
	if err != nil {
		t.Fatalf("NewDesign: %v", err)
	}
	p, err := d.Parameters()
	if err != nil {
		t.Fatalf("Parameters: %v", err)
	}
	if p.V != 7 || p.B != 7 || p.R != 3 || p.K != 3 || p.Lambda != 1 {
		t.Errorf("Fano params = %+v", p)
	}
	if !p.IsSymmetric() {
		t.Errorf("Fano should be symmetric")
	}
	// The dual of a symmetric design is also a 2-design with the same params.
	dual := d.Dual()
	pd, ok := dual.IsBIBD()
	if !ok || pd.V != 7 || pd.K != 3 || pd.Lambda != 1 {
		t.Errorf("dual params = %+v ok=%v", pd, ok)
	}
	if lam, ok := d.IsPairBalanced(); !ok || lam != 1 {
		t.Errorf("pair balance lambda=%d ok=%v", lam, ok)
	}
}

func ExampleFanoPlane() {
	pl := FanoPlane()
	fmt.Println("points:", pl.NumPoints())
	fmt.Println("lines:", pl.NumLines())
	d := pl.IncidenceDesign()
	p, _ := d.Parameters()
	fmt.Printf("design: 2-(%d,%d,%d)\n", p.V, p.K, p.Lambda)
	// Output:
	// points: 7
	// lines: 7
	// design: 2-(7,3,1)
}

func ExampleMOLS() {
	squares, _ := MOLS(3)
	fmt.Println("count:", len(squares))
	fmt.Println("orthogonal:", IsMOLS(squares))
	// Output:
	// count: 2
	// orthogonal: true
}
