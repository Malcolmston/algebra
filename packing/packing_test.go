package packing

import (
	"errors"
	"fmt"
	"math"
	"testing"
)

const tol = 1e-9

func approx(a, b, eps float64) bool { return math.Abs(a-b) <= eps }

// ---------------------------------------------------------------------------
// Bin packing.
// ---------------------------------------------------------------------------

func TestBinPackingHeuristics(t *testing.T) {
	sizes := []float64{4, 8, 1, 4, 2, 1}
	cap := 10.0
	algos := []struct {
		name string
		fn   func([]float64, float64) (Packing, error)
	}{
		{"NextFit", NextFit},
		{"FirstFit", FirstFit},
		{"BestFit", BestFit},
		{"WorstFit", WorstFit},
		{"AlmostWorstFit", AlmostWorstFit},
		{"NextFitDecreasing", NextFitDecreasing},
		{"FirstFitDecreasing", FirstFitDecreasing},
		{"BestFitDecreasing", BestFitDecreasing},
		{"WorstFitDecreasing", WorstFitDecreasing},
		{"AlmostWorstFitDecreasing", AlmostWorstFitDecreasing},
	}
	for _, a := range algos {
		p, err := a.fn(sizes, cap)
		if err != nil {
			t.Fatalf("%s: unexpected error %v", a.name, err)
		}
		if !p.Valid() {
			t.Errorf("%s: produced invalid packing %+v", a.name, p.Bins)
		}
		if p.NumItems() != len(sizes) {
			t.Errorf("%s: packed %d items, want %d", a.name, p.NumItems(), len(sizes))
		}
		if p.NumBins() < LowerBoundL1(sizes, cap) {
			t.Errorf("%s: used %d bins, below L1 lower bound %d", a.name, p.NumBins(), LowerBoundL1(sizes, cap))
		}
		if got := p.TotalSize(); !approx(got, 20, tol) {
			t.Errorf("%s: total size %v, want 20", a.name, got)
		}
	}
}

func TestBinPackingBinCounts(t *testing.T) {
	// Classic instance where next-fit is worse than first/best-fit.
	sizes := []float64{0.5, 0.5, 0.5, 0.5, 0.5, 0.5}
	cap := 1.0
	nf, _ := NextFit(sizes, cap)
	ff, _ := FirstFit(sizes, cap)
	if nf.NumBins() != 3 {
		t.Errorf("NextFit bins = %d, want 3", nf.NumBins())
	}
	if ff.NumBins() != 3 {
		t.Errorf("FirstFit bins = %d, want 3", ff.NumBins())
	}

	// An instance where decreasing order helps.
	s2 := []float64{6, 5, 4, 3, 2, 1}
	ff2, _ := FirstFit(s2, 10)
	ffd2, _ := FirstFitDecreasing(s2, 10)
	if ffd2.NumBins() > ff2.NumBins() {
		t.Errorf("FFD (%d) should not use more bins than FF (%d)", ffd2.NumBins(), ff2.NumBins())
	}
	if ffd2.NumBins() != 3 {
		t.Errorf("FFD bins = %d, want 3", ffd2.NumBins())
	}
}

func TestBinPackingErrors(t *testing.T) {
	if _, err := FirstFit([]float64{1, 2}, 0); !errors.Is(err, ErrCapacity) {
		t.Errorf("expected ErrCapacity, got %v", err)
	}
	if _, err := FirstFit([]float64{1, -2}, 10); !errors.Is(err, ErrNegativeSize) {
		t.Errorf("expected ErrNegativeSize, got %v", err)
	}
	if _, err := FirstFit([]float64{11}, 10); !errors.Is(err, ErrItemTooLarge) {
		t.Errorf("expected ErrItemTooLarge, got %v", err)
	}
}

func TestLowerBounds(t *testing.T) {
	tests := []struct {
		name     string
		sizes    []float64
		cap      float64
		wantL1   int
		wantL2   int
		wantBest int
	}{
		{"exact fill", []float64{5, 5, 5, 5}, 10, 2, 2, 2},
		{"triple 0.6", []float64{0.6, 0.6, 0.6}, 1, 2, 3, 3},
		{"quad 0.7", []float64{0.7, 0.7, 0.7, 0.7}, 1, 3, 4, 4},
		{"empty", nil, 10, 0, 0, 0},
		{"mixed", []float64{4, 8, 1, 4, 2, 1}, 10, 2, 2, 2},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := LowerBoundL1(tt.sizes, tt.cap); got != tt.wantL1 {
				t.Errorf("L1 = %d, want %d", got, tt.wantL1)
			}
			if got := LowerBoundL2(tt.sizes, tt.cap); got != tt.wantL2 {
				t.Errorf("L2 = %d, want %d", got, tt.wantL2)
			}
			if got := BinPackingLowerBound(tt.sizes, tt.cap); got != tt.wantBest {
				t.Errorf("BinPackingLowerBound = %d, want %d", got, tt.wantBest)
			}
		})
	}
}

func TestLowerBoundIsValid(t *testing.T) {
	// L2 must never exceed the number of bins any real packing uses.
	sizes := []float64{0.7, 0.7, 0.7, 0.6, 0.6, 0.3, 0.3, 0.4}
	cap := 1.0
	lb := BinPackingLowerBound(sizes, cap)
	ffd, _ := FirstFitDecreasing(sizes, cap)
	if lb > ffd.NumBins() {
		t.Errorf("lower bound %d exceeds a feasible packing %d", lb, ffd.NumBins())
	}
}

func TestPackingMetrics(t *testing.T) {
	sizes := []float64{4, 6, 3, 7}
	p, _ := FirstFit(sizes, 10)
	if got := p.TotalSize(); !approx(got, 20, tol) {
		t.Errorf("TotalSize = %v, want 20", got)
	}
	waste := float64(p.NumBins())*10 - 20
	if got := p.Waste(); !approx(got, waste, tol) {
		t.Errorf("Waste = %v, want %v", got, waste)
	}
	if p.MaxLoad() > 10+1e-9 {
		t.Errorf("MaxLoad %v exceeds capacity", p.MaxLoad())
	}
	if p.Fullness() < 0 || p.Fullness() > 1 {
		t.Errorf("Fullness %v out of range", p.Fullness())
	}
	if got := SumSizes(sizes); !approx(got, 20, tol) {
		t.Errorf("SumSizes = %v, want 20", got)
	}
	if got := MaxSize(sizes); got != 7 {
		t.Errorf("MaxSize = %v, want 7", got)
	}
}

func TestApproxRatios(t *testing.T) {
	if NextFitApproxRatio() != 2 || WorstFitApproxRatio() != 2 {
		t.Error("next/worst-fit ratio should be 2")
	}
	if FirstFitApproxRatio() != 1.7 || BestFitApproxRatio() != 1.7 {
		t.Error("first/best-fit ratio should be 1.7")
	}
	if !approx(FirstFitDecreasingApproxRatio(), 11.0/9.0, tol) {
		t.Error("FFD ratio should be 11/9")
	}
	if !approx(NextFitDecreasingApproxRatio(), 1.6910302, 1e-6) {
		t.Error("NFD ratio should be ~1.691")
	}
}

// ---------------------------------------------------------------------------
// Sphere packing.
// ---------------------------------------------------------------------------

func TestBallVolumes(t *testing.T) {
	tests := []struct {
		n    int
		want float64
	}{
		{0, 1},
		{1, 2},
		{2, math.Pi},
		{3, 4 * math.Pi / 3},
		{4, math.Pi * math.Pi / 2},
		{8, math.Pow(math.Pi, 4) / 24},
	}
	for _, tt := range tests {
		if got := UnitBallVolume(tt.n); !approx(got, tt.want, 1e-12) {
			t.Errorf("UnitBallVolume(%d) = %v, want %v", tt.n, got, tt.want)
		}
	}
	if got := BallVolume(3, 2); !approx(got, 8*4*math.Pi/3, 1e-9) {
		t.Errorf("BallVolume(3,2) = %v, want %v", got, 8*4*math.Pi/3)
	}
	if got := BallSurface(3, 1); !approx(got, 4*math.Pi, 1e-9) {
		t.Errorf("BallSurface(3,1) = %v, want 4pi", got)
	}
	if got := UnitSphereSurface(3); !approx(got, 4*math.Pi, 1e-9) {
		t.Errorf("UnitSphereSurface(3) = %v, want 4pi", got)
	}
}

func TestLatticeCenterDensities(t *testing.T) {
	tests := []struct {
		name string
		got  float64
		want float64
	}{
		{"Z^1", ZnCenterDensity(1), 0.5},
		{"Z^3", ZnCenterDensity(3), 0.125},
		{"A_2", AnCenterDensity(2), 1 / (2 * math.Sqrt(3))},
		{"A_3", AnCenterDensity(3), 1 / (4 * math.Sqrt2)},
		{"D_3", DnCenterDensity(3), 1 / (4 * math.Sqrt2)},
		{"E6", E6CenterDensity(), 1 / (8 * math.Sqrt(3))},
		{"E7", E7CenterDensity(), 1.0 / 16},
		{"E8", E8CenterDensity(), 1.0 / 16},
		{"Leech", LeechCenterDensity(), 1},
	}
	for _, tt := range tests {
		if !approx(tt.got, tt.want, 1e-12) {
			t.Errorf("%s center density = %v, want %v", tt.name, tt.got, tt.want)
		}
	}
}

func TestLatticeDensities(t *testing.T) {
	// A_2 achieves the optimal planar density pi/sqrt(12).
	if !approx(AnDensity(2), HexagonalPackingDensity(), 1e-12) {
		t.Errorf("A_2 density %v != hexagonal %v", AnDensity(2), HexagonalPackingDensity())
	}
	// E8 density = pi^4 / 384.
	if !approx(E8Density(), math.Pow(math.Pi, 4)/384, 1e-12) {
		t.Errorf("E8 density = %v, want pi^4/384", E8Density())
	}
	// Leech density = V_24.
	if !approx(LeechDensity(), UnitBallVolume(24), 1e-15) {
		t.Errorf("Leech density = %v, want V_24 = %v", LeechDensity(), UnitBallVolume(24))
	}
	// Every density lies in (0,1].
	for n := 1; n <= 8; n++ {
		if d := BestLatticeDensity(n); d <= 0 || d > 1+1e-12 {
			t.Errorf("BestLatticeDensity(%d) = %v out of range", n, d)
		}
	}
}

func TestKissingNumbers(t *testing.T) {
	tests := []struct {
		n    int
		want int
	}{{1, 2}, {2, 6}, {3, 12}, {4, 24}, {8, 240}, {24, 196560}}
	for _, tt := range tests {
		if got := KissingNumberOptimal(tt.n); got != tt.want {
			t.Errorf("KissingNumberOptimal(%d) = %d, want %d", tt.n, got, tt.want)
		}
		if !KissingNumberKnown(tt.n) {
			t.Errorf("KissingNumberKnown(%d) = false", tt.n)
		}
	}
	if KissingNumberOptimal(5) != -1 || KissingNumberKnown(5) {
		t.Error("dimension 5 kissing number should be unknown")
	}
	// Lattice lower bounds attain the optimum in the special dimensions.
	if KissingNumberLowerBound(8) != 240 {
		t.Errorf("KissingNumberLowerBound(8) = %d, want 240", KissingNumberLowerBound(8))
	}
	if KissingNumberLowerBound(24) != 196560 {
		t.Errorf("KissingNumberLowerBound(24) = %d, want 196560", KissingNumberLowerBound(24))
	}
	if AnKissingNumber(2) != 6 {
		t.Errorf("A_2 kissing = %d, want 6", AnKissingNumber(2))
	}
	if DnKissingNumber(4) != 24 {
		t.Errorf("D_4 kissing = %d, want 24", DnKissingNumber(4))
	}
}

func TestCoveringRadii(t *testing.T) {
	if !approx(ZnCoveringRadius(4), 1, tol) {
		t.Errorf("Z^4 covering radius = %v, want 1", ZnCoveringRadius(4))
	}
	if !approx(AnCoveringRadius(2), math.Sqrt(2.0/3.0), tol) {
		t.Errorf("A_2 covering radius = %v, want sqrt(2/3)", AnCoveringRadius(2))
	}
	if !approx(DnCoveringRadius(3), 1, tol) {
		t.Errorf("D_3 covering radius = %v, want 1", DnCoveringRadius(3))
	}
	if !approx(DnCoveringRadius(9), 1.5, tol) {
		t.Errorf("D_9 covering radius = %v, want 1.5", DnCoveringRadius(9))
	}
	if !approx(E8CoveringRadius(), 1, tol) {
		t.Errorf("E8 covering radius = %v, want 1", E8CoveringRadius())
	}
	if !approx(LeechCoveringRadius(), math.Sqrt2, tol) {
		t.Errorf("Leech covering radius = %v, want sqrt(2)", LeechCoveringRadius())
	}
	// A_n^* covering radius squared = n(n+2)/(12(n+1)).
	if !approx(AnStarCoveringRadius(3)*AnStarCoveringRadius(3), 15.0/48.0, tol) {
		t.Errorf("A_3^* covering radius^2 = %v, want 15/48", AnStarCoveringRadius(3)*AnStarCoveringRadius(3))
	}
	// Every thickness is at least 1.
	if ZnThickness(3) < 1 {
		t.Errorf("Z^3 thickness %v < 1", ZnThickness(3))
	}
	if AnStarThickness(3) < 1 {
		t.Errorf("A_3^* thickness %v < 1", AnStarThickness(3))
	}
}

func TestHermiteConstant(t *testing.T) {
	tests := []struct {
		n    int
		want float64
	}{
		{1, 1},
		{2, 2 / math.Sqrt(3)},
		{3, math.Cbrt(2)},
		{4, math.Sqrt2},
		{8, 2},
		{24, 4},
	}
	for _, tt := range tests {
		if got := HermiteConstant(tt.n); !approx(got, tt.want, 1e-12) {
			t.Errorf("HermiteConstant(%d) = %v, want %v", tt.n, got, tt.want)
		}
	}
	if HermiteConstant(9) != -1 || HermiteConstantKnown(9) {
		t.Error("Hermite constant should be unknown for n=9")
	}
	// The E8 center density equals (gamma_8/4)^(8/2) = (1/2)^4 = 1/16.
	if !approx(math.Pow(HermiteConstant(8)/4, 4), E8CenterDensity(), 1e-12) {
		t.Error("Hermite relation for E8 failed")
	}
}

func TestMinkowskiHlawka(t *testing.T) {
	// zeta(2)/2 = pi^2/12.
	if !approx(MinkowskiHlawkaBound(2), math.Pi*math.Pi/12, 1e-9) {
		t.Errorf("MinkowskiHlawkaBound(2) = %v, want pi^2/12", MinkowskiHlawkaBound(2))
	}
	// The bound is a lower bound: the best lattice we know is at least this dense.
	for n := 2; n <= 8; n++ {
		if BestLatticeDensity(n) < MinkowskiHlawkaBound(n)-1e-9 {
			t.Errorf("dim %d: best density %v below MH bound %v", n, BestLatticeDensity(n), MinkowskiHlawkaBound(n))
		}
	}
	if !math.IsNaN(MinkowskiHlawkaBound(1)) {
		t.Error("MinkowskiHlawkaBound(1) should be NaN")
	}
}

func TestLatticeInfo(t *testing.T) {
	e8 := E8Info()
	if e8.Dimension != 8 || e8.KissingNumber != 240 {
		t.Errorf("E8Info = %+v", e8)
	}
	if !approx(e8.CenterDensity, 1.0/16, 1e-12) {
		t.Errorf("E8Info center density = %v", e8.CenterDensity)
	}
	leech := LeechInfo()
	if leech.Dimension != 24 || leech.KissingNumber != 196560 {
		t.Errorf("LeechInfo = %+v", leech)
	}
	infos := []LatticeInfo{ZnInfo(3), AnInfo(4), DnInfo(5), E6Info(), E7Info(), E8Info(), LeechInfo()}
	for _, info := range infos {
		// Center density must equal rho^n / covolume.
		want := math.Pow(info.PackingRadius, float64(info.Dimension)) / info.Covolume
		if !approx(info.CenterDensity, want, 1e-12) {
			t.Errorf("%s: center density %v inconsistent with %v", info.Name, info.CenterDensity, want)
		}
		if info.Thickness < 1-1e-9 {
			t.Errorf("%s: thickness %v < 1", info.Name, info.Thickness)
		}
	}
}

// ---------------------------------------------------------------------------
// Circle packing.
// ---------------------------------------------------------------------------

func TestCircleInSquare(t *testing.T) {
	tests := []struct {
		n         int
		wantM     float64
		wantKnown bool
	}{
		{2, math.Sqrt2, true},
		{4, 1, true},
		{5, 1 / math.Sqrt2, true},
		{9, 0.5, true},
		{16, 1.0 / 3, true}, // 4x4 grid
		{25, 1.0 / 4, true}, // 5x5 grid
		{11, 0, false},      // not tabulated
	}
	for _, tt := range tests {
		m, ok := CircleInSquareSpread(tt.n)
		if ok != tt.wantKnown {
			t.Errorf("n=%d known = %v, want %v", tt.n, ok, tt.wantKnown)
			continue
		}
		if ok && !approx(m, tt.wantM, 1e-9) {
			t.Errorf("n=%d spread = %v, want %v", tt.n, m, tt.wantM)
		}
	}
	// n=1 special case.
	r, ok := CircleInSquareRadius(1)
	if !ok || !approx(r, 0.5, tol) {
		t.Errorf("CircleInSquareRadius(1) = %v,%v want 0.5,true", r, ok)
	}
	// Radius / density consistency for n=4: m=1, r=1/4, density=4*pi/16=pi/4.
	r4, _ := CircleInSquareRadius(4)
	if !approx(r4, 0.25, 1e-12) {
		t.Errorf("CircleInSquareRadius(4) = %v, want 0.25", r4)
	}
	d4, _ := CircleInSquareDensity(4)
	if !approx(d4, math.Pi/4, 1e-12) {
		t.Errorf("CircleInSquareDensity(4) = %v, want pi/4", d4)
	}
	if !CircleInSquareKnown(16) || CircleInSquareKnown(11) {
		t.Error("CircleInSquareKnown mismatch")
	}
	if !approx(GridSpread(3), 0.5, tol) {
		t.Errorf("GridSpread(3) = %v, want 0.5", GridSpread(3))
	}
}

func TestCircleInCircle(t *testing.T) {
	tests := []struct {
		n    int
		want float64
	}{
		{1, 1},
		{2, 2},
		{3, 1 + 2/math.Sqrt(3)},
		{4, 1 + math.Sqrt2},
		{6, 3},
		{7, 3},
	}
	for _, tt := range tests {
		got, ok := CircleInCircleRatio(tt.n)
		if !ok || !approx(got, tt.want, 1e-9) {
			t.Errorf("CircleInCircleRatio(%d) = %v,%v want %v", tt.n, got, ok, tt.want)
		}
	}
	if _, ok := CircleInCircleRatio(50); ok {
		t.Error("n=50 should not be tabulated")
	}
	// Density for n=2 is 2/4 = 0.5.
	d2, _ := CircleInCircleDensity(2)
	if !approx(d2, 0.5, 1e-12) {
		t.Errorf("CircleInCircleDensity(2) = %v, want 0.5", d2)
	}
	// Ring formula: 4 circles => 1 + sqrt(2).
	if !approx(RingRatio(4), 1+math.Sqrt2, 1e-12) {
		t.Errorf("RingRatio(4) = %v, want 1+sqrt2", RingRatio(4))
	}
	if !approx(RingRatio(6), 3, 1e-9) {
		t.Errorf("RingRatio(6) = %v, want 3", RingRatio(6))
	}
	if !approx(RingDensity(6), 6.0/9, 1e-9) {
		t.Errorf("RingDensity(6) = %v, want 2/3", RingDensity(6))
	}
	if !CircleInCircleKnown(6) || CircleInCircleKnown(99) {
		t.Error("CircleInCircleKnown mismatch")
	}
}

func TestPlanarConstants(t *testing.T) {
	if !approx(HexagonalPackingDensity(), math.Pi/math.Sqrt(12), 1e-12) {
		t.Errorf("hex packing density = %v", HexagonalPackingDensity())
	}
	if HexagonalPackingDensity() <= SquarePackingDensity() {
		t.Error("hexagonal packing should beat square packing")
	}
	if !approx(SquarePackingDensity(), math.Pi/4, 1e-12) {
		t.Errorf("square packing density = %v", SquarePackingDensity())
	}
	if HexagonalCoveringThickness() >= SquareCoveringThickness() {
		t.Error("hexagonal covering should be thinner than square covering")
	}
	if HexagonalCoveringThickness() < 1 {
		t.Error("covering thickness must be >= 1")
	}
	if HexagonalKissingNumber() != 6 || SquareKissingNumber() != 4 {
		t.Error("planar kissing numbers wrong")
	}
}

// ---------------------------------------------------------------------------
// Knapsack and set cover.
// ---------------------------------------------------------------------------

func TestFractionalKnapsack(t *testing.T) {
	values := []float64{60, 100, 120}
	weights := []float64{10, 20, 30}
	v, frac, err := FractionalKnapsack(values, weights, 50)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(v, 240, 1e-9) {
		t.Errorf("fractional value = %v, want 240", v)
	}
	if !approx(frac[0], 1, tol) || !approx(frac[1], 1, tol) || !approx(frac[2], 2.0/3, 1e-9) {
		t.Errorf("fractions = %v, want [1 1 0.6667]", frac)
	}
	// The LP relaxation upper-bounds the integer optimum.
	ub, _ := KnapsackUpperBound(values, weights, 50)
	gv, _, _, _ := GreedyKnapsack(values, weights, 50)
	if gv > ub+1e-9 {
		t.Errorf("greedy value %v exceeds LP bound %v", gv, ub)
	}
}

func TestGreedyKnapsack(t *testing.T) {
	values := []float64{60, 100, 120}
	weights := []float64{10, 20, 30}
	v, chosen, w, err := GreedyKnapsack(values, weights, 50)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(v, 160, tol) {
		t.Errorf("greedy value = %v, want 160", v)
	}
	if !approx(w, 30, tol) || len(chosen) != 2 {
		t.Errorf("greedy chose %v (weight %v)", chosen, w)
	}
	if !KnapsackFeasible(weights, chosen, 50) {
		t.Error("greedy selection should be feasible")
	}
	if !approx(KnapsackValue(values, chosen), v, tol) {
		t.Error("KnapsackValue mismatch")
	}
	if !approx(KnapsackWeight(weights, chosen), w, tol) {
		t.Error("KnapsackWeight mismatch")
	}

	// Case where the modified greedy beats plain greedy: one heavy high-value
	// item versus several tiny low-value items.
	v2 := []float64{100, 1, 1, 1}
	w2 := []float64{10, 1, 1, 1}
	// Plain greedy prefers ratio-1 small items first; capacity 10 fits the big
	// item alone for value 100.
	best, _, _, _ := GreedyKnapsackBest(v2, w2, 10)
	if best < 100-1e-9 {
		t.Errorf("modified greedy value = %v, want >= 100", best)
	}
	if ModifiedGreedyKnapsackApproxRatio() != 2 {
		t.Error("modified greedy ratio should be 2")
	}
}

func TestKnapsackErrors(t *testing.T) {
	if _, _, err := FractionalKnapsack([]float64{1}, []float64{1, 2}, 5); !errors.Is(err, ErrDimensionMismatch) {
		t.Errorf("expected ErrDimensionMismatch, got %v", err)
	}
	if _, _, err := FractionalKnapsack([]float64{1}, []float64{-1}, 5); !errors.Is(err, ErrNegativeWeight) {
		t.Errorf("expected ErrNegativeWeight, got %v", err)
	}
}

func TestSetCover(t *testing.T) {
	sets := [][]int{{0, 1, 2}, {2, 3}, {3, 4}, {0, 1, 4}}
	chosen, err := SetCoverGreedy(5, sets)
	if err != nil {
		t.Fatal(err)
	}
	if !IsSetCover(5, sets, chosen) {
		t.Errorf("greedy result %v does not cover the universe", chosen)
	}
	// Greedy first picks the largest set {0,1,2}.
	if chosen[0] != 0 {
		t.Errorf("greedy first pick = %d, want 0", chosen[0])
	}

	// Uncoverable universe.
	if _, err := SetCoverGreedy(6, sets); !errors.Is(err, ErrNotCoverable) {
		t.Errorf("expected ErrNotCoverable, got %v", err)
	}

	// Weighted: a cheap pair of sets should beat one expensive covering set.
	wsets := [][]int{{0, 1, 2, 3, 4}, {0, 1, 2}, {3, 4}}
	costs := []float64{100, 1, 1}
	wc, cost, err := WeightedSetCoverGreedy(5, wsets, costs)
	if err != nil {
		t.Fatal(err)
	}
	if !IsSetCover(5, wsets, wc) {
		t.Errorf("weighted result %v does not cover", wc)
	}
	if cost > 2+1e-9 {
		t.Errorf("weighted greedy cost = %v, want <= 2", cost)
	}
}

func TestHarmonicAndBounds(t *testing.T) {
	if !approx(HarmonicNumber(0), 0, tol) {
		t.Error("H(0) should be 0")
	}
	if !approx(HarmonicNumber(4), 25.0/12, 1e-12) {
		t.Errorf("H(4) = %v, want 25/12", HarmonicNumber(4))
	}
	sets := [][]int{{0, 1, 2}, {2, 3}}
	if MaxSetSize(sets) != 3 {
		t.Errorf("MaxSetSize = %d, want 3", MaxSetSize(sets))
	}
	if !approx(SetCoverGreedyBound(3), HarmonicNumber(3), tol) {
		t.Error("SetCoverGreedyBound(3) should equal H(3)")
	}
}

// ---------------------------------------------------------------------------
// Runnable examples.
// ---------------------------------------------------------------------------

func ExampleFirstFitDecreasing() {
	sizes := []float64{6, 5, 4, 3, 2, 1}
	p, _ := FirstFitDecreasing(sizes, 10)
	fmt.Printf("bins=%d lowerbound=%d\n", p.NumBins(), BinPackingLowerBound(sizes, 10))
	// Output: bins=3 lowerbound=3
}

func ExampleE8Density() {
	fmt.Printf("E8 kissing number: %d\n", E8KissingNumber())
	fmt.Printf("E8 center density: %.4f\n", E8CenterDensity())
	// Output:
	// E8 kissing number: 240
	// E8 center density: 0.0625
}

func ExampleSetCoverGreedy() {
	sets := [][]int{{0, 1, 2}, {2, 3}, {3, 4}, {0, 1, 4}}
	chosen, _ := SetCoverGreedy(5, sets)
	fmt.Printf("sets chosen: %d, covers all: %v\n", len(chosen), IsSetCover(5, sets, chosen))
	// Output: sets chosen: 2, covers all: true
}

func ExampleFractionalKnapsack() {
	value, _, _ := FractionalKnapsack([]float64{60, 100, 120}, []float64{10, 20, 30}, 50)
	fmt.Printf("optimal fractional value: %.0f\n", value)
	// Output: optimal fractional value: 240
}
