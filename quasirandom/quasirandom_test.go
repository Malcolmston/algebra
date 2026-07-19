package quasirandom

import (
	"fmt"
	"math"
	"math/big"
	"testing"
)

const tol = 1e-12

func close(a, b, t float64) bool {
	if math.IsNaN(a) || math.IsNaN(b) {
		return false
	}
	return math.Abs(a-b) <= t
}

func vecClose(a, b []float64, t float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !close(a[i], b[i], t) {
			return false
		}
	}
	return true
}

func inUnit(p []float64) bool {
	for _, x := range p {
		if x < 0 || x >= 1 {
			return false
		}
	}
	return true
}

func TestIsPrimeAndPrime(t *testing.T) {
	cases := []struct {
		n    int
		want bool
	}{
		{-1, false}, {0, false}, {1, false}, {2, true}, {3, true},
		{4, false}, {17, true}, {21, false}, {97, true}, {100, false},
	}
	for _, c := range cases {
		if got := IsPrime(c.n); got != c.want {
			t.Errorf("IsPrime(%d)=%v want %v", c.n, got, c.want)
		}
	}
	wantPrimes := []int{2, 3, 5, 7, 11, 13, 17, 19, 23, 29}
	for i, w := range wantPrimes {
		p, err := Prime(i + 1)
		if err != nil || p != w {
			t.Errorf("Prime(%d)=%d,%v want %d", i+1, p, err, w)
		}
	}
	if NextPrime(13) != 17 || PrevPrime(17) != 13 || NextPrimeGE(17) != 17 {
		t.Errorf("prime navigation wrong")
	}
	if PrimeIndex(13) != 6 {
		t.Errorf("PrimeIndex(13)=%d want 6", PrimeIndex(13))
	}
	if !AreCoprimeBases([]int{2, 3, 5}) || AreCoprimeBases([]int{2, 4}) {
		t.Errorf("AreCoprimeBases wrong")
	}
}

func TestRadicalInverse(t *testing.T) {
	cases := []struct {
		base int
		n    uint64
		want float64
	}{
		{2, 0, 0}, {2, 1, 0.5}, {2, 2, 0.25}, {2, 3, 0.75}, {2, 4, 0.125},
		{2, 5, 0.625}, {2, 6, 0.375}, {2, 7, 0.875},
		{3, 1, 1.0 / 3}, {3, 2, 2.0 / 3}, {3, 3, 1.0 / 9}, {3, 4, 4.0 / 9},
		{5, 1, 0.2}, {10, 12, 0.21},
	}
	for _, c := range cases {
		got, err := RadicalInverse(c.base, c.n)
		if err != nil || !close(got, c.want, tol) {
			t.Errorf("RadicalInverse(%d,%d)=%v,%v want %v", c.base, c.n, got, err, c.want)
		}
		if c.base == 2 {
			if b2 := RadicalInverseBase2(c.n); !close(b2, c.want, tol) {
				t.Errorf("RadicalInverseBase2(%d)=%v want %v", c.n, b2, c.want)
			}
		}
	}
	if _, err := RadicalInverse(1, 3); err != ErrBadBase {
		t.Errorf("expected ErrBadBase")
	}
}

func TestScrambledRadicalInverseIdentity(t *testing.T) {
	perm := IdentityPermutation(5)
	for n := uint64(0); n < 50; n++ {
		a, _ := RadicalInverse(5, n)
		b, err := ScrambledRadicalInverse(5, n, perm)
		if err != nil || !close(a, b, tol) {
			t.Errorf("scrambled identity mismatch at %d: %v vs %v", n, a, b)
		}
	}
	if _, err := ScrambledRadicalInverse(5, 1, []int{0, 1}); err != ErrDimension {
		t.Errorf("expected ErrDimension for bad perm length")
	}
}

func TestVanDerCorputSequence(t *testing.T) {
	seq, err := VanDerCorputSequence(2, 8)
	want := []float64{0, 0.5, 0.25, 0.75, 0.125, 0.625, 0.375, 0.875}
	if err != nil || !vecClose(seq, want, tol) {
		t.Errorf("VanDerCorputSequence base 2 = %v", seq)
	}
}

func TestDigits(t *testing.T) {
	d, _ := Digits(13, 2)
	if !intsEqual(d, []int{1, 0, 1, 1}) {
		t.Errorf("Digits(13,2)=%v", d)
	}
	n, _ := FromDigits([]int{1, 0, 1, 1}, 2)
	if n != 13 {
		t.Errorf("FromDigits=%d", n)
	}
	if s, _ := DigitSum(255, 2); s != 8 {
		t.Errorf("DigitSum(255,2)=%d", s)
	}
	if g := GrayCode(6); g != 5 || InverseGrayCode(5) != 6 {
		t.Errorf("gray code roundtrip failed")
	}
	if r := ReverseBits(0b1011, 4); r != 0b1101 {
		t.Errorf("ReverseBits=%b", r)
	}
	if !IsPowerOfTwo(16) || IsPowerOfTwo(17) || NextPowerOfTwo(17) != 32 {
		t.Errorf("power-of-two helpers wrong")
	}
	if PopCount(255) != 8 || BitLength(255) != 8 {
		t.Errorf("bit helpers wrong")
	}
}

func intsEqual(a, b []int) bool {
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

func TestHalton(t *testing.T) {
	p, err := Halton(2, 1)
	if err != nil || !vecClose(p, []float64{0.5, 1.0 / 3}, tol) {
		t.Errorf("Halton(2,1)=%v", p)
	}
	p3, _ := Halton(3, 4)
	// bases 2,3,5: phi2(4)=0.125, phi3(4)=4/9, phi5(4)=0.8
	if !vecClose(p3, []float64{0.125, 4.0 / 9, 0.8}, tol) {
		t.Errorf("Halton(3,4)=%v", p3)
	}
	seq, _ := HaltonSequence(2, 20)
	for i, q := range seq {
		if !inUnit(q) {
			t.Errorf("Halton point %d out of unit cube: %v", i, q)
		}
	}
	c, _ := HaltonCoordinate(2, 1, 4)
	if !close(c, 4.0/9, tol) {
		t.Errorf("HaltonCoordinate=%v", c)
	}
}

func TestHammersley(t *testing.T) {
	set, err := HammersleySet(2, 4)
	want := [][]float64{{0, 0}, {0.25, 0.5}, {0.5, 0.25}, {0.75, 0.75}}
	if err != nil {
		t.Fatal(err)
	}
	for i := range want {
		if !vecClose(set[i], want[i], tol) {
			t.Errorf("Hammersley point %d = %v want %v", i, set[i], want[i])
		}
	}
	if _, err := Hammersley(2, 4, 4); err != ErrDimension {
		t.Errorf("expected out-of-range error")
	}
}

func TestFaure(t *testing.T) {
	if b, _ := FaureBase(2); b != 2 {
		t.Errorf("FaureBase(2)=%d", b)
	}
	if b, _ := FaureBase(4); b != 5 {
		t.Errorf("FaureBase(4)=%d", b)
	}
	f, _ := Faure(3, 3)
	if !vecClose(f, []float64{1.0 / 9, 4.0 / 9, 7.0 / 9}, tol) {
		t.Errorf("Faure(3,3)=%v", f)
	}
	f2, _ := Faure(2, 2)
	if !vecClose(f2, []float64{0.25, 0.75}, tol) {
		t.Errorf("Faure(2,2)=%v", f2)
	}
	// First coordinate of Faure equals the radical inverse in the base.
	for n := uint64(0); n < 40; n++ {
		fp, _ := Faure(4, n)
		ri, _ := RadicalInverse(5, n)
		if !close(fp[0], ri, tol) {
			t.Errorf("Faure first coord mismatch at %d", n)
		}
		if !inUnit(fp) {
			t.Errorf("Faure point out of cube at %d: %v", n, fp)
		}
	}
}

func TestFaureGeneratorMatrix(t *testing.T) {
	// Power-zero matrix is the identity Pascal matrix truncation.
	m, err := FaureGeneratorMatrix(5, 0, 3)
	if err != nil {
		t.Fatal(err)
	}
	want := [][]int{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}
	for i := range want {
		if !intsEqual(m[i], want[i]) {
			t.Errorf("row %d = %v want %v", i, m[i], want[i])
		}
	}
}

func TestSobol(t *testing.T) {
	want := [][]float64{
		{0.5, 0.5}, {0.75, 0.25}, {0.25, 0.75}, {0.375, 0.375}, {0.875, 0.875},
	}
	for i, w := range want {
		p, err := SobolPoint(2, uint64(i+1))
		if err != nil || !vecClose(p, w, tol) {
			t.Errorf("SobolPoint(2,%d)=%v want %v", i+1, p, w)
		}
	}
	// The first Sobol dimension is the base-2 radical inverse read in
	// Gray-code order: x_n = phi_2(gray(n)).
	for n := uint64(1); n < 64; n++ {
		p, _ := SobolPoint(1, n)
		ri := RadicalInverseBase2(GrayCode(n))
		if !close(p[0], ri, tol) {
			t.Errorf("Sobol dim1 mismatch at %d: %v vs %v", n, p[0], ri)
		}
	}
}

func TestSobolGeneratorStreamAndPoint(t *testing.T) {
	s, err := NewSobol(3)
	if err != nil {
		t.Fatal(err)
	}
	first := s.Next()
	if !vecClose(first, []float64{0, 0, 0}, tol) {
		t.Errorf("first Sobol point should be origin, got %v", first)
	}
	// Streaming order must match direct addressing.
	s.Reset()
	for i := uint64(0); i < 100; i++ {
		got := s.Next()
		want := s.Point(i)
		if !vecClose(got, want, tol) {
			t.Errorf("stream vs Point mismatch at %d: %v vs %v", i, got, want)
		}
	}
	// Skip is equivalent to advancing.
	a, _ := NewSobol(2)
	b, _ := NewSobol(2)
	a.Skip(10)
	for i := 0; i < 10; i++ {
		b.Next()
	}
	if a.Index() != b.Index() {
		t.Errorf("Skip index mismatch")
	}
	if !vecClose(a.Next(), b.Next(), tol) {
		t.Errorf("Skip point mismatch")
	}
	if _, err := NewSobol(MaxSobolDimension() + 1); err != ErrDimension {
		t.Errorf("expected dimension error above MaxSobolDimension")
	}
}

func TestSobolDirectionData(t *testing.T) {
	// Dimension 1 has all direction numbers m_i = 1, so v_i = 2^-i.
	dn, err := DirectionNumbers(1)
	if err != nil {
		t.Fatal(err)
	}
	for i, v := range dn {
		if !close(v, math.Ldexp(1, -(i+1)), 1e-15) {
			t.Errorf("DirectionNumbers(1)[%d]=%v", i, v)
		}
	}
	// Polynomial x^3 + x + 1 for dimension 4 (encoding 1, degree 3).
	coeffs, _ := PrimitivePolynomialCoefficients(4)
	if !intsEqual(coeffs, []int{1, 0, 1, 1}) {
		t.Errorf("dim4 primitive polynomial coeffs=%v", coeffs)
	}
}

func TestFaurePermutation(t *testing.T) {
	cases := []struct {
		base int
		want []int
	}{
		{2, []int{0, 1}},
		{3, []int{0, 1, 2}},
		{4, []int{0, 2, 1, 3}},
		{5, []int{0, 3, 2, 1, 4}},
	}
	for _, c := range cases {
		got := FaurePermutation(c.base)
		if !intsEqual(got, c.want) {
			t.Errorf("FaurePermutation(%d)=%v want %v", c.base, got, c.want)
		}
		if !IsPermutation(got) {
			t.Errorf("FaurePermutation(%d) is not a permutation", c.base)
		}
	}
}

func TestPermutations(t *testing.T) {
	p := []int{2, 0, 1}
	inv := InvertPermutation(p)
	comp := ComposePermutations(p, inv)
	if !intsEqual(comp, IdentityPermutation(3)) {
		t.Errorf("p∘p^-1 = %v", comp)
	}
	if IsPermutation([]int{0, 0}) {
		t.Errorf("non-permutation reported as permutation")
	}
	aff := AffinePermutation(5, 2, 1)
	if !IsPermutation(aff) {
		t.Errorf("AffinePermutation(5,2,1) not a permutation: %v", aff)
	}
}

func TestL2StarDiscrepancy(t *testing.T) {
	// Single centered point in 1D: sqrt(1/12).
	d, err := L2StarDiscrepancy([][]float64{{0.5}})
	if err != nil || !close(d, math.Sqrt(1.0/12), 1e-12) {
		t.Errorf("L2StarDiscrepancy single point = %v", d)
	}
	// Regular grid at (2i-1)/(2N) has known discrepancy; monotone improvement
	// for a low-discrepancy set versus a clustered one.
	good, _ := HammersleySet(2, 64)
	bad := make([][]float64, 64)
	for i := range bad {
		bad[i] = []float64{0.01, 0.01}
	}
	dg, _ := L2StarDiscrepancy(good)
	db, _ := L2StarDiscrepancy(bad)
	if !(dg < db) {
		t.Errorf("expected Hammersley discrepancy %v < clustered %v", dg, db)
	}
}

func TestStarDiscrepancyExactVs1D(t *testing.T) {
	pts := [][]float64{{0.1}, {0.3}, {0.55}, {0.8}, {0.95}}
	exact, _ := StarDiscrepancy(pts)
	formula, _ := StarDiscrepancy1D(pts)
	if !close(exact, formula, 1e-12) {
		t.Errorf("exact %v vs 1D closed form %v", exact, formula)
	}
	// Two points example with a hand value.
	if sd, _ := StarDiscrepancy1D([][]float64{{0.25}, {0.75}}); !close(sd, 0.25, tol) {
		t.Errorf("StarDiscrepancy1D=%v want 0.25", sd)
	}
}

func TestStarDiscrepancy2D(t *testing.T) {
	pts := [][]float64{{0.25, 0.25}, {0.75, 0.75}}
	exact, _ := StarDiscrepancy(pts)
	grid, _ := StarDiscrepancyGrid(pts, 200)
	// The grid estimate is a lower bound on the exact value.
	if grid > exact+1e-9 {
		t.Errorf("grid estimate %v exceeds exact %v", grid, exact)
	}
	if exact <= 0 || exact >= 1 {
		t.Errorf("exact discrepancy out of range: %v", exact)
	}
}

func TestOtherDiscrepancies(t *testing.T) {
	single := [][]float64{{0.5}}
	if d, _ := CenteredL2Discrepancy(single); !close(d, math.Sqrt(1.0/12), 1e-12) {
		t.Errorf("CenteredL2 single point = %v", d)
	}
	if d, _ := Diaphony(single); !close(d, math.Pi/math.Sqrt(3), 1e-12) {
		t.Errorf("Diaphony single point = %v want %v", d, math.Pi/math.Sqrt(3))
	}
	// All discrepancy variants are finite and non-negative for a real set.
	pts, _ := HammersleySet(2, 16)
	for name, fn := range map[string]func([][]float64) (float64, error){
		"wrap":     WrapAroundL2Discrepancy,
		"sym":      SymmetricL2Discrepancy,
		"mod":      ModifiedL2StarDiscrepancy,
		"centered": CenteredL2Discrepancy,
		"diaphony": Diaphony,
	} {
		v, err := fn(pts)
		if err != nil || v < 0 || math.IsNaN(v) {
			t.Errorf("%s discrepancy invalid: %v,%v", name, v, err)
		}
	}
	if _, err := L2StarDiscrepancy(nil); err != ErrEmptyPointSet {
		t.Errorf("expected empty point set error")
	}
}

func TestLocalDiscrepancy(t *testing.T) {
	pts := [][]float64{{0.2, 0.2}, {0.8, 0.8}}
	// Box [0,0.5)^2 contains exactly the first point: 1/2 - 0.25 = 0.25.
	d, err := LocalDiscrepancy(pts, []float64{0.5, 0.5})
	if err != nil || !close(d, 0.25, tol) {
		t.Errorf("LocalDiscrepancy=%v", d)
	}
	c, _ := CountInBox(pts, []float64{0.5, 0.5})
	if c != 1 {
		t.Errorf("CountInBox=%d", c)
	}
}

func TestKroneckerAndRoberts(t *testing.T) {
	if !close(WeylSequence(GoldenRatio, 1), Frac(GoldenRatio), tol) {
		t.Errorf("WeylSequence golden ratio wrong")
	}
	// Plastic constant satisfies x^(d+1) = x + 1.
	g1, _ := PlasticConstant(1)
	if !close(g1, GoldenRatio, 1e-12) {
		t.Errorf("PlasticConstant(1)=%v want golden ratio", g1)
	}
	g2, _ := PlasticConstant(2)
	if !close(g2*g2*g2, g2+1, 1e-12) {
		t.Errorf("PlasticConstant(2) does not satisfy its equation")
	}
	seq, _ := RobertsSequence(2, 50)
	for i, p := range seq {
		if !inUnit(p) {
			t.Errorf("Roberts point %d out of cube: %v", i, p)
		}
	}
	kp, _ := KroneckerPoint([]float64{GoldenRatio, math.Sqrt2}, 3)
	if !inUnit(kp) {
		t.Errorf("Kronecker point out of cube: %v", kp)
	}
}

func TestLattices(t *testing.T) {
	if Fibonacci(10) != 55 || Fibonacci(7) != 13 {
		t.Errorf("Fibonacci wrong")
	}
	// Fibonacci lattice with F_7 = 13 points, generator (1, F_6=8).
	lat, err := FibonacciLattice(7)
	if err != nil || len(lat) != 13 {
		t.Fatalf("FibonacciLattice(7) len=%d err=%v", len(lat), err)
	}
	for _, p := range lat {
		if !inUnit(p) {
			t.Errorf("lattice point out of cube: %v", p)
		}
	}
	gen, _ := KorobovGenerator(3, 3, 100)
	if !intsEqual(gen, []int{1, 3, 9}) {
		t.Errorf("KorobovGenerator=%v want [1 3 9]", gen)
	}
	lp, _ := RankOneLatticePoint([]int{1, 8}, 13, 5)
	if !vecClose(lp, []float64{5.0 / 13, Frac(40.0 / 13)}, tol) {
		t.Errorf("RankOneLatticePoint=%v", lp)
	}
}

func TestIntegration(t *testing.T) {
	// ∫_0^1 x^2 dx = 1/3.
	f1 := func(x []float64) float64 { return x[0] * x[0] }
	v, err := QMCIntegrateHalton(f1, 1, 4000)
	if err != nil || !close(v, 1.0/3, 5e-3) {
		t.Errorf("QMC ∫x^2 = %v", v)
	}
	// ∫_[0,1]^2 (x+y) = 1.
	f2 := func(x []float64) float64 { return x[0] + x[1] }
	v2, _ := QMCIntegrateSobol(f2, 2, 4096)
	if !close(v2, 1.0, 5e-3) {
		t.Errorf("QMC ∫(x+y) = %v", v2)
	}
	v3, _ := QMCIntegrateFaure(f2, 2, 4000)
	if !close(v3, 1.0, 5e-3) {
		t.Errorf("Faure QMC ∫(x+y) = %v", v3)
	}
	// Integrate over a box: ∫_[1,3] x dx = 4.
	fb := func(x []float64) float64 { return x[0] }
	pts, _ := HaltonSequence(1, 4000)
	vb, _ := IntegrateBox(fb, []float64{1}, []float64{3}, pts)
	if !close(vb, 4.0, 1e-2) {
		t.Errorf("IntegrateBox = %v want 4", vb)
	}
}

func TestStatisticsHelpers(t *testing.T) {
	vals := []float64{1, 2, 3, 4, 5}
	if !close(SampleMean(vals), 3, tol) {
		t.Errorf("SampleMean=%v", SampleMean(vals))
	}
	if !close(SampleVariance(vals), 2.5, tol) {
		t.Errorf("SampleVariance=%v", SampleVariance(vals))
	}
	if KoksmaHlawkaBound(2, 0.1) != 0.2 {
		t.Errorf("KoksmaHlawkaBound wrong")
	}
}

func TestPointSetUtilities(t *testing.T) {
	pts := [][]float64{{0.2, 0.4}, {0.6, 0.8}}
	lo, hi, _ := BoundingBox(pts)
	if !vecClose(lo, []float64{0.2, 0.4}, tol) || !vecClose(hi, []float64{0.6, 0.8}, tol) {
		t.Errorf("BoundingBox=%v,%v", lo, hi)
	}
	c, _ := Centroid(pts)
	if !vecClose(c, []float64{0.4, 0.6}, tol) {
		t.Errorf("Centroid=%v", c)
	}
	d, _ := MinimumDistance(pts)
	if !close(d, math.Hypot(0.4, 0.4), tol) {
		t.Errorf("MinimumDistance=%v", d)
	}
	scaled, _ := ScaleToBox(pts, []float64{0, 0}, []float64{2, 2})
	back, _ := BoxToUnit(scaled, []float64{0, 0}, []float64{2, 2})
	if !vecClose(back[0], pts[0], tol) {
		t.Errorf("ScaleToBox/BoxToUnit roundtrip failed: %v", back[0])
	}
	refl := ReflectPoint([]float64{0.25, 0.75})
	if !vecClose(refl, []float64{0.75, 0.25}, tol) {
		t.Errorf("ReflectPoint=%v", refl)
	}
	proj, _ := Project(pts, []int{1})
	if !vecClose(proj[0], []float64{0.4}, tol) {
		t.Errorf("Project=%v", proj[0])
	}
}

func TestRadicalInverseBig(t *testing.T) {
	for n := int64(0); n < 200; n++ {
		f, _ := RadicalInverseBigFloat(2, big.NewInt(n))
		g, _ := RadicalInverse(2, uint64(n))
		if !close(f, g, 1e-12) {
			t.Errorf("big vs uint radical inverse mismatch at %d: %v vs %v", n, f, g)
		}
	}
	r, _ := RadicalInverseBig(3, big.NewInt(4))
	if r.Cmp(big.NewRat(4, 9)) != 0 {
		t.Errorf("RadicalInverseBig(3,4)=%v want 4/9", r)
	}
	b, _ := BigBinomial(5, 2)
	if b.Int64() != 10 {
		t.Errorf("BigBinomial(5,2)=%v", b)
	}
}

func TestDiscrepancyDecreasesWithN(t *testing.T) {
	// The Sobol L2 star discrepancy should shrink as N grows.
	var prev float64 = math.Inf(1)
	for _, n := range []int{16, 64, 256, 1024} {
		pts, _ := SobolSequenceSkip(2, 1, n)
		d, _ := L2StarDiscrepancy(pts)
		if d > prev {
			t.Errorf("Sobol discrepancy increased at N=%d: %v > %v", n, d, prev)
		}
		prev = d
	}
}

// ExampleRadicalInverse prints the first four terms of the base-two van der
// Corput sequence.
func ExampleRadicalInverse() {
	for n := uint64(1); n <= 4; n++ {
		v, _ := RadicalInverse(2, n)
		fmt.Printf("%.3f\n", v)
	}
	// Output:
	// 0.500
	// 0.250
	// 0.750
	// 0.125
}

// ExampleSobol prints the first few two-dimensional Sobol points.
func ExampleSobol() {
	s, _ := NewSobol(2)
	for i := 0; i < 5; i++ {
		p := s.Next()
		fmt.Printf("%.3f %.3f\n", p[0], p[1])
	}
	// Output:
	// 0.000 0.000
	// 0.500 0.500
	// 0.750 0.250
	// 0.250 0.750
	// 0.375 0.375
}

// ExampleHalton prints the first three two-dimensional Halton points.
func ExampleHalton() {
	for n := uint64(0); n < 3; n++ {
		p, _ := Halton(2, n)
		fmt.Printf("%.4f %.4f\n", p[0], p[1])
	}
	// Output:
	// 0.0000 0.0000
	// 0.5000 0.3333
	// 0.2500 0.6667
}
