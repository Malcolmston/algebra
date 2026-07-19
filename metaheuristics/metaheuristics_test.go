package metaheuristics

import (
	"fmt"
	"math"
	"testing"
)

const eps = 1e-9

func almostEqual(a, b, tol float64) bool {
	if math.IsInf(a, 0) || math.IsInf(b, 0) {
		return a == b
	}
	return math.Abs(a-b) <= tol
}

func TestBenchmarkOptima(t *testing.T) {
	tests := []struct {
		name string
		f    ObjectiveFunc
		x    []float64
		want float64
		tol  float64
	}{
		{"Sphere@0", Sphere, []float64{0, 0, 0}, 0, eps},
		{"Sphere@1", Sphere, []float64{1, 2, 2}, 9, eps},
		{"Rastrigin@0", Rastrigin, []float64{0, 0, 0, 0}, 0, eps},
		{"Ackley@0", Ackley, []float64{0, 0}, 0, 1e-9},
		{"Rosenbrock@1", Rosenbrock, []float64{1, 1, 1}, 0, eps},
		{"Rosenbrock@0", Rosenbrock, []float64{0, 0}, 1, eps},
		{"Griewank@0", Griewank, []float64{0, 0, 0}, 0, eps},
		{"Schwefel@opt", Schwefel, []float64{420.9687, 420.9687}, 0, 1e-3},
		{"Zakharov@0", Zakharov, []float64{0, 0}, 0, eps},
		{"StyblinskiTang@opt", StyblinskiTang, []float64{-2.903534, -2.903534}, -78.33234, 1e-3},
		{"Booth@opt", Booth, []float64{1, 3}, 0, eps},
		{"Matyas@0", Matyas, []float64{0, 0}, 0, eps},
		{"Beale@opt", Beale, []float64{3, 0.5}, 0, eps},
		{"Himmelblau@opt", Himmelblau, []float64{3, 2}, 0, eps},
		{"ThreeHumpCamel@0", ThreeHumpCamel, []float64{0, 0}, 0, eps},
		{"SixHumpCamel@opt", SixHumpCamel, []float64{0.0898, -0.7126}, -1.031628, 1e-4},
		{"Easom@opt", Easom, []float64{math.Pi, math.Pi}, -1, eps},
		{"LeviN13@opt", LeviN13, []float64{1, 1}, 0, 1e-9},
		{"GoldsteinPrice@opt", GoldsteinPrice, []float64{0, -1}, 3, 1e-9},
		{"McCormick@opt", McCormick, []float64{-0.54719, -1.54719}, -1.9133, 1e-4},
		{"DixonPrice@0-first", DixonPrice, []float64{1, 1 / math.Sqrt2}, 0, 1e-9},
		{"SumSquares@0", SumSquares, []float64{0, 0}, 0, eps},
		{"Michalewicz@0-plateau", Michalewicz, []float64{0, 0}, 0, eps},
		{"Powell@0", Powell, []float64{0, 0, 0, 0}, 0, eps},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := tt.f(tt.x)
			if !almostEqual(got, tt.want, tt.tol) {
				t.Errorf("%s = %v, want %v (tol %v)", tt.name, got, tt.want, tt.tol)
			}
		})
	}
}

func TestTridOptimum(t *testing.T) {
	// For n=2 the Trid global minimum is -2 at (2,2).
	got := Trid([]float64{2, 2})
	if !almostEqual(got, -2, eps) {
		t.Errorf("Trid(2,2) = %v, want -2", got)
	}
}

func TestBoundsMethods(t *testing.T) {
	b, err := NewBounds([]float64{-1, -2}, []float64{1, 2})
	if err != nil {
		t.Fatal(err)
	}
	if b.Dim() != 2 {
		t.Errorf("Dim = %d, want 2", b.Dim())
	}
	if !b.Valid() {
		t.Error("expected valid bounds")
	}
	c := b.Center()
	if !almostEqual(c[0], 0, eps) || !almostEqual(c[1], 0, eps) {
		t.Errorf("Center = %v, want [0 0]", c)
	}
	w := b.Width()
	if !almostEqual(w[0], 2, eps) || !almostEqual(w[1], 4, eps) {
		t.Errorf("Width = %v, want [2 4]", w)
	}
	if !b.Contains([]float64{0, 0}) {
		t.Error("expected origin inside box")
	}
	if b.Contains([]float64{5, 5}) {
		t.Error("expected (5,5) outside box")
	}
	clipped := b.Clip([]float64{5, -5})
	if !almostEqual(clipped[0], 1, eps) || !almostEqual(clipped[1], -2, eps) {
		t.Errorf("Clip = %v, want [1 -2]", clipped)
	}
	refl := b.Reflect([]float64{1.5, 0})
	if refl[0] < -1 || refl[0] > 1 {
		t.Errorf("Reflect out of box: %v", refl)
	}
}

func TestBoundsErrors(t *testing.T) {
	if _, err := NewBounds([]float64{0}, []float64{1, 2}); err != ErrDimMismatch {
		t.Errorf("want ErrDimMismatch, got %v", err)
	}
	if _, err := NewBounds(nil, nil); err != ErrEmptyBounds {
		t.Errorf("want ErrEmptyBounds, got %v", err)
	}
}

func TestVectorHelpers(t *testing.T) {
	a := []float64{1, 2, 3}
	b := []float64{4, 5, 6}
	if got := VecDot(a, b); !almostEqual(got, 32, eps) {
		t.Errorf("VecDot = %v, want 32", got)
	}
	if got := VecNorm([]float64{3, 4}); !almostEqual(got, 5, eps) {
		t.Errorf("VecNorm = %v, want 5", got)
	}
	if got := VecNorm1([]float64{-1, 2, -3}); !almostEqual(got, 6, eps) {
		t.Errorf("VecNorm1 = %v, want 6", got)
	}
	if got := VecNormInf([]float64{-1, 2, -3}); !almostEqual(got, 3, eps) {
		t.Errorf("VecNormInf = %v, want 3", got)
	}
	if got := VecDist([]float64{0, 0}, []float64{3, 4}); !almostEqual(got, 5, eps) {
		t.Errorf("VecDist = %v, want 5", got)
	}
	sum := VecAdd(a, b)
	if !almostEqual(sum[0], 5, eps) || !almostEqual(sum[2], 9, eps) {
		t.Errorf("VecAdd = %v", sum)
	}
	sc := VecScale(a, 2)
	if !almostEqual(sc[1], 4, eps) {
		t.Errorf("VecScale = %v", sc)
	}
	axpy := VecAXPY(a, 2, b)
	if !almostEqual(axpy[0], 9, eps) {
		t.Errorf("VecAXPY = %v", axpy)
	}
	mean := VecMean([][]float64{{0, 0}, {2, 4}})
	if !almostEqual(mean[0], 1, eps) || !almostEqual(mean[1], 2, eps) {
		t.Errorf("VecMean = %v", mean)
	}
	if got := Lerp(0, 10, 0.25); !almostEqual(got, 2.5, eps) {
		t.Errorf("Lerp = %v, want 2.5", got)
	}
	if got := Clamp(5, 0, 3); !almostEqual(got, 3, eps) {
		t.Errorf("Clamp = %v, want 3", got)
	}
}

func TestRNGDeterminism(t *testing.T) {
	a := NewRNG(42)
	b := NewRNG(42)
	for i := 0; i < 100; i++ {
		if a.Float64() != b.Float64() {
			t.Fatalf("RNG not deterministic at %d", i)
		}
	}
	// different seeds diverge
	c := NewRNG(1)
	d := NewRNG(2)
	same := true
	for i := 0; i < 10; i++ {
		if c.Float64() != d.Float64() {
			same = false
			break
		}
	}
	if same {
		t.Error("different seeds produced identical stream")
	}
}

func TestRNGChoice(t *testing.T) {
	rng := NewRNG(7)
	// weight all on index 2
	counts := make([]int, 3)
	for i := 0; i < 1000; i++ {
		counts[rng.Choice([]float64{0, 0, 1})]++
	}
	if counts[2] != 1000 {
		t.Errorf("Choice with single weight: counts=%v", counts)
	}
}

func TestAcceptanceProbability(t *testing.T) {
	if p := AcceptanceProbability(1, 0, 1); p != 1 {
		t.Errorf("downhill move probability = %v, want 1", p)
	}
	if p := AcceptanceProbability(0, 1, 0); p != 0 {
		t.Errorf("uphill at T=0 probability = %v, want 0", p)
	}
	// exp(-1) at deltaE=1, T=1
	if p := AcceptanceProbability(0, 1, 1); !almostEqual(p, math.Exp(-1), 1e-12) {
		t.Errorf("uphill probability = %v, want %v", p, math.Exp(-1))
	}
}

func TestCoolingSchedules(t *testing.T) {
	exp := ExponentialCooling(0.5)
	if !almostEqual(exp(100, 0, 10), 100, eps) {
		t.Error("exp cooling at k=0 should be T0")
	}
	if !almostEqual(exp(100, 2, 10), 25, eps) {
		t.Errorf("exp cooling k=2 = %v, want 25", exp(100, 2, 10))
	}
	lin := LinearCooling()
	if !almostEqual(lin(100, 5, 10), 50, eps) {
		t.Errorf("linear cooling = %v, want 50", lin(100, 5, 10))
	}
	if !almostEqual(lin(100, 10, 10), 0, eps) {
		t.Errorf("linear cooling at end = %v, want 0", lin(100, 10, 10))
	}
	cau := CauchyCooling()
	if !almostEqual(cau(100, 0, 10), 100, eps) {
		t.Error("cauchy at k=0 should be T0")
	}
	if !almostEqual(cau(100, 1, 10), 50, eps) {
		t.Errorf("cauchy k=1 = %v, want 50", cau(100, 1, 10))
	}
	// schedules are monotonically non-increasing
	for _, s := range []CoolingSchedule{ExponentialCooling(0.9), LinearCooling(), LogarithmicCooling(), BoltzmannCooling(), CauchyCooling(), QuadraticCooling()} {
		prev := math.Inf(1)
		for k := 0; k < 50; k++ {
			v := s(10, k, 50)
			if v > prev+1e-9 {
				t.Errorf("schedule not non-increasing at k=%d: %v > %v", k, v, prev)
			}
			prev = v
		}
	}
}

// optimizerConverges is a helper asserting an optimizer gets close to 0 on the
// 2-D Sphere function.
func runSphere2D(t *testing.T, name string, f func(b Bounds, rng *RNG) (Result, error), tol float64) {
	t.Helper()
	b := UniformBounds(2, -5.12, 5.12)
	rng := NewRNG(20240101)
	res, err := f(b, rng)
	if err != nil {
		t.Fatalf("%s: %v", name, err)
	}
	if res.F > tol {
		t.Errorf("%s: F = %v, want <= %v (X=%v)", name, res.F, tol, res.X)
	}
	if len(res.X) != 2 {
		t.Errorf("%s: result dim = %d, want 2", name, len(res.X))
	}
}

func TestOptimizersOnSphere(t *testing.T) {
	runSphere2D(t, "GA", func(b Bounds, rng *RNG) (Result, error) {
		cfg := DefaultGAConfig(b)
		cfg.Generations = 150
		return RunGA(Sphere, cfg, rng)
	}, 1e-2)

	runSphere2D(t, "PSO", func(b Bounds, rng *RNG) (Result, error) {
		return RunPSO(Sphere, DefaultPSOConfig(b), rng)
	}, 1e-6)

	runSphere2D(t, "DE", func(b Bounds, rng *RNG) (Result, error) {
		cfg := DefaultDEConfig(b)
		cfg.Generations = 200
		return RunDE(Sphere, cfg, rng)
	}, 1e-8)

	runSphere2D(t, "CMAES", func(b Bounds, rng *RNG) (Result, error) {
		cfg := DefaultCMAESConfig(b)
		cfg.MaxIterations = 200
		return RunCMAES(Sphere, cfg, rng)
	}, 1e-8)

	runSphere2D(t, "SA", func(b Bounds, rng *RNG) (Result, error) {
		cfg := DefaultAnnealConfig(b)
		cfg.MaxIterations = 8000
		return SimulatedAnnealing(Sphere, cfg, nil, rng)
	}, 1e-2)

	runSphere2D(t, "Harmony", func(b Bounds, rng *RNG) (Result, error) {
		return RunHarmonySearch(Sphere, DefaultHarmonyConfig(b), rng)
	}, 1e-2)

	runSphere2D(t, "HillClimbRestarts", func(b Bounds, rng *RNG) (Result, error) {
		cfg := DefaultHillClimbConfig(b)
		cfg.MaxIterations = 500
		return HillClimbRestarts(Sphere, cfg, 10, rng)
	}, 1e-3)

	runSphere2D(t, "ContinuousTabu", func(b Bounds, rng *RNG) (Result, error) {
		cfg := DefaultContinuousTabuConfig(b)
		cfg.Iterations = 2000
		return ContinuousTabuSearch(Sphere, cfg, nil, rng)
	}, 1e-1)

	runSphere2D(t, "CoordinateDescent", func(b Bounds, rng *RNG) (Result, error) {
		cfg := DefaultHillClimbConfig(b)
		cfg.MaxIterations = 200
		return CoordinateDescent(Sphere, cfg, []float64{4, -3})
	}, 1e-6)
}

func TestDEStrategiesOnRastrigin(t *testing.T) {
	b := UniformBounds(2, -5.12, 5.12)
	strategies := []DEStrategy{DERand1, DEBest1, DECurrentToBest1, DERand2, DEBest2}
	for _, s := range strategies {
		t.Run(s.String(), func(t *testing.T) {
			rng := NewRNG(99)
			cfg := DefaultDEConfig(b)
			cfg.Strategy = s
			cfg.Generations = 400
			cfg.PopSize = 40
			res, err := RunDE(Rastrigin, cfg, rng)
			if err != nil {
				t.Fatal(err)
			}
			if res.F > 1e-4 {
				t.Errorf("%s on Rastrigin: F = %v", s, res.F)
			}
		})
	}
}

func TestCMAESOnRosenbrock(t *testing.T) {
	b := UniformBounds(2, -5, 10)
	rng := NewRNG(2024)
	cfg := DefaultCMAESConfig(b)
	cfg.MaxIterations = 500
	res, err := RunCMAES(Rosenbrock, cfg, rng)
	if err != nil {
		t.Fatal(err)
	}
	if res.F > 1e-6 {
		t.Errorf("CMA-ES on Rosenbrock: F = %v, X = %v", res.F, res.X)
	}
}

func TestJacobiEigen(t *testing.T) {
	// symmetric matrix with known eigenvalues 2 and 4 (diag after rotation)
	a := [][]float64{{3, 1}, {1, 3}}
	vals, vecs := JacobiEigen(a)
	// eigenvalues of [[3,1],[1,3]] are 2 and 4
	got := []float64{vals[0], vals[1]}
	found2, found4 := false, false
	for _, v := range got {
		if almostEqual(v, 2, 1e-9) {
			found2 = true
		}
		if almostEqual(v, 4, 1e-9) {
			found4 = true
		}
	}
	if !found2 || !found4 {
		t.Errorf("eigenvalues = %v, want {2,4}", got)
	}
	// eigenvectors orthonormal: columns dot to 0
	dot := vecs[0][0]*vecs[0][1] + vecs[1][0]*vecs[1][1]
	if !almostEqual(dot, 0, 1e-9) {
		t.Errorf("eigenvectors not orthogonal, dot=%v", dot)
	}
}

func TestTSPHelpers(t *testing.T) {
	// unit square
	cities := []City{{0, 0}, {0, 1}, {1, 1}, {1, 0}}
	dist := DistanceMatrix(cities)
	if !almostEqual(dist[0][2], math.Sqrt2, eps) {
		t.Errorf("diagonal distance = %v, want sqrt2", dist[0][2])
	}
	// perimeter tour has length 4
	if l := TourLength([]int{0, 1, 2, 3}, dist); !almostEqual(l, 4, eps) {
		t.Errorf("square perimeter = %v, want 4", l)
	}
	// a crossing tour is longer
	crossing := TourLength([]int{0, 2, 1, 3}, dist)
	if crossing <= 4 {
		t.Errorf("crossing tour %v should exceed 4", crossing)
	}
	// 2-opt should recover the perimeter tour
	improved, l := TwoOpt([]int{0, 2, 1, 3}, dist)
	if !almostEqual(l, 4, eps) {
		t.Errorf("2-opt length = %v, want 4 (tour %v)", l, improved)
	}
	nn := NearestNeighborTour(dist, 0)
	if len(nn) != 4 {
		t.Errorf("NN tour length = %d, want 4", len(nn))
	}
}

func TestACOFindsOptimalSquare(t *testing.T) {
	cities := []City{{0, 0}, {0, 1}, {1, 1}, {1, 0}}
	dist := DistanceMatrix(cities)
	rng := NewRNG(123)
	cfg := DefaultACOConfig()
	cfg.Ants = 10
	cfg.Iterations = 50
	res, err := RunACO(dist, cfg, rng)
	if err != nil {
		t.Fatal(err)
	}
	if !almostEqual(res.Length, 4, 1e-9) {
		t.Errorf("ACO best length = %v, want 4 (tour %v)", res.Length, res.Tour)
	}
}

func TestACOOnRandomInstance(t *testing.T) {
	rng := NewRNG(555)
	n := 12
	cities := make([]City, n)
	for i := range cities {
		cities[i] = City{X: rng.Float64() * 100, Y: rng.Float64() * 100}
	}
	dist := DistanceMatrix(cities)
	nn := TourLength(NearestNeighborTour(dist, 0), dist)
	cfg := DefaultACOConfig()
	cfg.Ants = 20
	cfg.Iterations = 100
	res, err := RunACO(dist, cfg, NewRNG(1))
	if err != nil {
		t.Fatal(err)
	}
	// ACO should do at least as well as nearest-neighbour construction.
	if res.Length > nn+1e-9 {
		t.Errorf("ACO length %v worse than NN %v", res.Length, nn)
	}
}

func TestTabuCombinatorial(t *testing.T) {
	// Minimize a simple quadratic over integer states in [0,10]^1 encoded as a
	// single element; neighbours are +/-1.
	target := 7
	energy := func(s []int) float64 {
		d := float64(s[0] - target)
		return d * d
	}
	neighbors := func(s []int) []TabuMove {
		var moves []TabuMove
		for _, delta := range []int{-1, 1} {
			v := s[0] + delta
			if v < 0 || v > 10 {
				continue
			}
			ns := []int{v}
			moves = append(moves, TabuMove{State: ns, Value: energy(ns), Key: fmt.Sprint(v)})
		}
		return moves
	}
	best, val, _, err := TabuSearchCombinatorial([]int{0}, neighbors, TabuSearchConfig{Tenure: 3, Iterations: 50})
	if err != nil {
		t.Fatal(err)
	}
	if best[0] != target || !almostEqual(val, 0, eps) {
		t.Errorf("tabu best = %v (val %v), want %d", best, val, target)
	}
}

func TestTabuListEviction(t *testing.T) {
	tl := NewTabuList(2)
	tl.Add("a")
	tl.Add("b")
	if !tl.Contains("a") || !tl.Contains("b") {
		t.Error("expected a,b in list")
	}
	tl.Add("c") // evicts a
	if tl.Contains("a") {
		t.Error("expected a evicted")
	}
	if !tl.Contains("c") {
		t.Error("expected c present")
	}
	if tl.Len() != 2 {
		t.Errorf("Len = %d, want 2", tl.Len())
	}
	tl.Clear()
	if tl.Len() != 0 {
		t.Errorf("after Clear Len = %d, want 0", tl.Len())
	}
}

func TestGACrossoverOperators(t *testing.T) {
	rng := NewRNG(3)
	a := []float64{0, 0, 0, 0}
	b := []float64{1, 1, 1, 1}
	// arithmetic children are convex combinations -> in [0,1]
	c1, c2 := ArithmeticCrossover(a, b, rng)
	for i := range c1 {
		if c1[i] < -eps || c1[i] > 1+eps || c2[i] < -eps || c2[i] > 1+eps {
			t.Errorf("arithmetic child out of range: %v %v", c1, c2)
		}
	}
	// one-point crossover preserves gene values from parents
	d1, d2 := OnePointCrossover(a, b, rng)
	for i := range d1 {
		if d1[i] != 0 && d1[i] != 1 {
			t.Errorf("one-point child gene %v not from a parent", d1[i])
		}
		if d1[i]+d2[i] != 1 {
			t.Errorf("one-point children not complementary at %d", i)
		}
	}
}

func TestNegateAndShift(t *testing.T) {
	neg := Negate(Sphere)
	if !almostEqual(neg([]float64{2}), -4, eps) {
		t.Errorf("Negate = %v, want -4", neg([]float64{2}))
	}
	// shift moves optimum: Sphere shifted by [1] has minimum at x=1
	sh := Shift(Sphere, []float64{1})
	if !almostEqual(sh([]float64{1}), 0, eps) {
		t.Errorf("Shifted sphere at 1 = %v, want 0", sh([]float64{1}))
	}
	sc := Scale(Sphere, 2, 3)
	if !almostEqual(sc([]float64{1}), 5, eps) {
		t.Errorf("Scale = %v, want 5", sc([]float64{1}))
	}
	pen := Penalized(Sphere, UniformBounds(1, -1, 1), 10)
	// x=2 violates by 1 -> penalty 10, plus f=4 -> 14
	if !almostEqual(pen([]float64{2}), 14, eps) {
		t.Errorf("Penalized = %v, want 14", pen([]float64{2}))
	}
}

func TestStandardBenchmarks(t *testing.T) {
	bs := StandardBenchmarks()
	if len(bs) < 5 {
		t.Fatalf("expected several benchmarks, got %d", len(bs))
	}
	for _, bm := range bs {
		b := bm.BoundsFor(3)
		if b.Dim() != 3 {
			t.Errorf("%s BoundsFor(3) dim = %d", bm.Name, b.Dim())
		}
		gm := bm.GlobalMin(3)
		if math.IsNaN(gm) {
			t.Errorf("%s global min is NaN", bm.Name)
		}
	}
}

func TestConstrictionFactor(t *testing.T) {
	// classic phi = 4.1 gives chi ~ 0.7298
	chi := ConstrictionFactor(2.05, 2.05)
	if !almostEqual(chi, 0.729843788, 1e-6) {
		t.Errorf("ConstrictionFactor = %v, want ~0.7298", chi)
	}
	// phi <= 4 returns 1
	if c := ConstrictionFactor(1, 1); c != 1 {
		t.Errorf("ConstrictionFactor(1,1) = %v, want 1", c)
	}
}

func TestReproducibility(t *testing.T) {
	b := UniformBounds(2, -5, 5)
	cfg := DefaultDEConfig(b)
	cfg.Generations = 50
	r1, _ := RunDE(Sphere, cfg, NewRNG(77))
	r2, _ := RunDE(Sphere, cfg, NewRNG(77))
	if r1.F != r2.F {
		t.Errorf("DE not reproducible: %v vs %v", r1.F, r2.F)
	}
	if r1.X[0] != r2.X[0] || r1.X[1] != r2.X[1] {
		t.Errorf("DE X not reproducible: %v vs %v", r1.X, r2.X)
	}
}

func TestInvalidConfigs(t *testing.T) {
	b := UniformBounds(2, -1, 1)
	if _, err := RunGA(Sphere, GAConfig{Bounds: b, PopSize: 1, Generations: 1}, NewRNG(1)); err != ErrInvalidConfig {
		t.Errorf("want ErrInvalidConfig for tiny pop, got %v", err)
	}
	if _, err := RunPSO(Sphere, PSOConfig{Bounds: Bounds{}}, NewRNG(1)); err != ErrEmptyBounds {
		t.Errorf("want ErrEmptyBounds, got %v", err)
	}
}

// ExampleRunDE demonstrates minimizing the 2-D sphere function with
// differential evolution using a fixed seed for reproducibility.
func ExampleRunDE() {
	bounds := UniformBounds(2, -5.12, 5.12)
	cfg := DefaultDEConfig(bounds)
	cfg.Generations = 200
	res, err := RunDE(Sphere, cfg, NewRNG(1))
	if err != nil {
		fmt.Println("error:", err)
		return
	}
	fmt.Printf("minimum ~ %.6f\n", res.F)
	// Output: minimum ~ 0.000000
}

// ExampleSimulatedAnnealing demonstrates finding the minimum of the Rastrigin
// function with simulated annealing.
func ExampleSimulatedAnnealing() {
	bounds := UniformBounds(2, -5.12, 5.12)
	cfg := DefaultAnnealConfig(bounds)
	cfg.MaxIterations = 20000
	res, _ := SimulatedAnnealing(Rastrigin, cfg, nil, NewRNG(2025))
	fmt.Printf("f < 1e-2: %v\n", res.F < 1e-2)
	// Output: f < 1e-2: true
}
