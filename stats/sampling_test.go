package stats

import (
	"math"
	"math/rand"
	"testing"
)

// statsSampleMoments draws n variates via draw and returns their sample mean
// and variance.
func statsSampleMoments(n int, draw func() float64) (mean, variance float64) {
	xs := make([]float64, n)
	for i := range xs {
		xs[i] = draw()
	}
	return Mean(xs), Variance(xs)
}

// TestSampleFormulaKnownAnswer checks the closed-form samplers against their
// defining formula. A second Rand seeded identically reproduces the exact
// uniform draws the sampler consumes, so equality is exact (no tolerance).
func TestSampleFormulaKnownAnswer(t *testing.T) {
	const seed = 20260718
	cases := []struct {
		name string
		got  func(r *rand.Rand) float64
		want func(r *rand.Rand) float64
	}{
		{
			"Uniform",
			func(r *rand.Rand) float64 { return Uniform{A: -2, B: 5}.Sample(r) },
			func(r *rand.Rand) float64 { return -2 + 7*r.Float64() },
		},
		{
			"Exponential",
			func(r *rand.Rand) float64 { return Exponential{Lambda: 1.5}.Sample(r) },
			func(r *rand.Rand) float64 { return Exponential{Lambda: 1.5}.Quantile(r.Float64()) },
		},
		{
			"Weibull",
			func(r *rand.Rand) float64 { return Weibull{Shape: 2, Scale: 3}.Sample(r) },
			func(r *rand.Rand) float64 { return Weibull{Shape: 2, Scale: 3}.Quantile(r.Float64()) },
		},
		{
			"LogNormal",
			func(r *rand.Rand) float64 { return LogNormal{Mu: 0.2, Sigma: 0.7}.Sample(r) },
			func(r *rand.Rand) float64 { return LogNormal{Mu: 0.2, Sigma: 0.7}.Quantile(r.Float64()) },
		},
	}
	for _, c := range cases {
		rg := rand.New(rand.NewSource(seed))
		rw := rand.New(rand.NewSource(seed))
		for i := 0; i < 50; i++ {
			got, want := c.got(rg), c.want(rw)
			if got != want {
				t.Fatalf("%s[%d]: got %v want %v", c.name, i, got, want)
			}
		}
	}
}

// TestGeometricKnownAnswer checks Geometric.SampleInt against its exact
// inverse-transform formula floor(ln(1-U)/ln(1-P)).
func TestGeometricKnownAnswer(t *testing.T) {
	const seed = 7
	const p = 0.25
	rg := rand.New(rand.NewSource(seed))
	rw := rand.New(rand.NewSource(seed))
	for i := 0; i < 50; i++ {
		got := Geometric{P: p}.SampleInt(rg)
		u := rw.Float64()
		want := int(math.Floor(math.Log(1-u) / math.Log(1-p)))
		if got != want {
			t.Fatalf("Geometric[%d]: got %d want %d", i, got, want)
		}
	}
}

// TestSampleContinuousMoments checks that the continuous samplers reproduce the
// theoretical mean and variance of each distribution within tolerance.
func TestSampleContinuousMoments(t *testing.T) {
	r := rand.New(rand.NewSource(42))
	const n = 300000
	cases := []struct {
		name       string
		draw       func() float64
		mean, vari float64
		mtol, vtol float64
	}{
		{"Uniform", func() float64 { return Uniform{A: 0, B: 4}.Sample(r) }, 2, 16.0 / 12, 0.02, 0.05},
		{"Exponential", func() float64 { return Exponential{Lambda: 2}.Sample(r) }, 0.5, 0.25, 0.01, 0.02},
		{"Normal", func() float64 { return Normal{Mu: 2, Sigma: 3}.Sample(r) }, 2, 9, 0.05, 0.2},
		{"Gamma", func() float64 { return Gamma{Shape: 3, Scale: 2}.Sample(r) }, 6, 12, 0.1, 0.5},
		{"GammaSmallShape", func() float64 { return Gamma{Shape: 0.5, Scale: 2}.Sample(r) }, 1, 2, 0.05, 0.2},
		{"ChiSquared", func() float64 { return ChiSquared{K: 4}.Sample(r) }, 4, 8, 0.1, 0.4},
		{"Beta", func() float64 { return Beta{Alpha: 2, Beta: 3}.Sample(r) }, 0.4, 0.04, 0.01, 0.01},
		{"StudentT", func() float64 { return StudentT{Nu: 10}.Sample(r) }, 0, 10.0 / 8, 0.02, 0.1},
		{"FDist", func() float64 { return FDist{D1: 10, D2: 10}.Sample(r) }, 10.0 / 8, 0, 0.05, math.Inf(1)},
		{"Weibull", func() float64 { return Weibull{Shape: 2, Scale: 1}.Sample(r) }, math.Sqrt(math.Pi) / 2, 1 - math.Pi/4, 0.01, 0.02},
		{"LogNormal", func() float64 { return LogNormal{Mu: 0, Sigma: 0.5}.Sample(r) }, math.Exp(0.125), (math.Exp(0.25) - 1) * math.Exp(0.25), 0.02, 0.05},
	}
	for _, c := range cases {
		m, v := statsSampleMoments(n, c.draw)
		if math.Abs(m-c.mean) > c.mtol {
			t.Errorf("%s: mean=%.4f want %.4f", c.name, m, c.mean)
		}
		if !math.IsInf(c.vtol, 1) && math.Abs(v-c.vari) > c.vtol {
			t.Errorf("%s: var=%.4f want %.4f", c.name, v, c.vari)
		}
	}
}

// TestSampleDiscreteMoments checks that the discrete samplers reproduce the
// theoretical mean (and variance where cheap to state) within tolerance.
func TestSampleDiscreteMoments(t *testing.T) {
	r := rand.New(rand.NewSource(99))
	const n = 300000

	pm, pv := statsSampleMoments(n, func() float64 { return float64(Poisson{Lambda: 4}.SampleInt(r)) })
	if math.Abs(pm-4) > 0.05 || math.Abs(pv-4) > 0.15 {
		t.Errorf("Poisson(4): mean=%.3f var=%.3f want 4,4", pm, pv)
	}

	lpm, _ := statsSampleMoments(n, func() float64 { return float64(Poisson{Lambda: 60}.SampleInt(r)) })
	if math.Abs(lpm-60) > 0.3 {
		t.Errorf("Poisson(60): mean=%.3f want 60", lpm)
	}

	bm, bv := statsSampleMoments(n, func() float64 { return float64(Binomial{N: 10, P: 0.3}.SampleInt(r)) })
	if math.Abs(bm-3) > 0.05 || math.Abs(bv-2.1) > 0.1 {
		t.Errorf("Binomial(10,0.3): mean=%.3f var=%.3f want 3,2.1", bm, bv)
	}

	// Large N*P exercises the per-trial branch together with P>0.5 reflection.
	lbm, _ := statsSampleMoments(n, func() float64 { return float64(Binomial{N: 200, P: 0.8}.SampleInt(r)) })
	if math.Abs(lbm-160) > 0.5 {
		t.Errorf("Binomial(200,0.8): mean=%.3f want 160", lbm)
	}

	// Geometric with the failures parameterization has mean (1-P)/P.
	gm, _ := statsSampleMoments(n, func() float64 { return float64(Geometric{P: 0.25}.SampleInt(r)) })
	if math.Abs(gm-3) > 0.1 {
		t.Errorf("Geometric(0.25): mean=%.3f want 3", gm)
	}

	// NegativeBinomial failures count has mean R*(1-P)/P.
	nbm, _ := statsSampleMoments(n, func() float64 { return float64(NegativeBinomial{R: 5, P: 0.5}.SampleInt(r)) })
	if math.Abs(nbm-5) > 0.1 {
		t.Errorf("NegativeBinomial(5,0.5): mean=%.3f want 5", nbm)
	}

	// Hypergeometric: N=50, K=20, Draws=10 -> mean = Draws*K/N = 4.
	hm, _ := statsSampleMoments(n, func() float64 { return float64(Hypergeometric{N: 50, K: 20, Draws: 10}.SampleInt(r)) })
	if math.Abs(hm-4) > 0.05 {
		t.Errorf("Hypergeometric: mean=%.3f want 4", hm)
	}
}

// TestHypergeometricSupport verifies every draw lands inside the exact support
// max(0, Draws-(N-K)) <= k <= min(K, Draws).
func TestHypergeometricSupport(t *testing.T) {
	r := rand.New(rand.NewSource(3))
	h := Hypergeometric{N: 20, K: 12, Draws: 15}
	lo := 15 - (20 - 12) // 7
	hi := 12             // min(K, Draws)
	for i := 0; i < 10000; i++ {
		k := h.SampleInt(r)
		if k < lo || k > hi {
			t.Fatalf("draw %d out of support [%d,%d]", k, lo, hi)
		}
	}
}

// TestSampleDeterministic verifies that two Rand sources with the same seed
// yield byte-identical sequences from every sampler.
func TestSampleDeterministic(t *testing.T) {
	const seed = 123
	conts := []ContinuousSampler{
		Uniform{A: 0, B: 1}, Exponential{Lambda: 2}, Weibull{Shape: 1.5, Scale: 2},
		LogNormal{Mu: 0, Sigma: 1}, Normal{Mu: 0, Sigma: 1}, Gamma{Shape: 2.5, Scale: 1},
		ChiSquared{K: 3}, Beta{Alpha: 2, Beta: 5}, StudentT{Nu: 7}, FDist{D1: 5, D2: 9},
	}
	for _, d := range conts {
		r1 := rand.New(rand.NewSource(seed))
		r2 := rand.New(rand.NewSource(seed))
		for i := 0; i < 200; i++ {
			if a, b := d.Sample(r1), d.Sample(r2); a != b {
				t.Fatalf("%T not deterministic at %d: %v != %v", d, i, a, b)
			}
		}
	}
	discs := []DiscreteSampler{
		Poisson{Lambda: 4}, Poisson{Lambda: 55}, Binomial{N: 10, P: 0.3},
		Binomial{N: 150, P: 0.7}, Geometric{P: 0.2}, NegativeBinomial{R: 4, P: 0.4},
		Hypergeometric{N: 40, K: 15, Draws: 12},
	}
	for _, d := range discs {
		r1 := rand.New(rand.NewSource(seed))
		r2 := rand.New(rand.NewSource(seed))
		for i := 0; i < 200; i++ {
			if a, b := d.SampleInt(r1), d.SampleInt(r2); a != b {
				t.Fatalf("%T not deterministic at %d: %d != %d", d, i, a, b)
			}
		}
	}
}

// TestSampleNAndSampleIntN checks the batch helpers return the requested length
// and reproduce the distribution mean.
func TestSampleNAndSampleIntN(t *testing.T) {
	r := rand.New(rand.NewSource(5))
	xs := SampleN(Normal{Mu: 10, Sigma: 2}, 50000, r)
	if len(xs) != 50000 {
		t.Fatalf("SampleN length=%d", len(xs))
	}
	if math.Abs(Mean(xs)-10) > 0.05 {
		t.Errorf("SampleN mean=%.3f want 10", Mean(xs))
	}

	ks := SampleIntN(Poisson{Lambda: 7}, 50000, r)
	if len(ks) != 50000 {
		t.Fatalf("SampleIntN length=%d", len(ks))
	}
	sum := 0
	for _, k := range ks {
		sum += k
	}
	if m := float64(sum) / float64(len(ks)); math.Abs(m-7) > 0.1 {
		t.Errorf("SampleIntN mean=%.3f want 7", m)
	}
}

// TestShuffle verifies Shuffle is an in-place permutation (preserves the
// multiset) and is deterministic for a fixed seed.
func TestShuffle(t *testing.T) {
	orig := []float64{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}
	a := append([]float64(nil), orig...)
	Shuffle(a, rand.New(rand.NewSource(1)))

	sum, prod := 0.0, 1.0
	for _, v := range a {
		sum += v
		prod *= v
	}
	if sum != 55 || prod != 3628800 {
		t.Fatalf("Shuffle changed multiset: sum=%v prod=%v", sum, prod)
	}

	b := append([]float64(nil), orig...)
	Shuffle(b, rand.New(rand.NewSource(1)))
	for i := range a {
		if a[i] != b[i] {
			t.Fatalf("Shuffle not deterministic at %d: %v != %v", i, a[i], b[i])
		}
	}

	// A single-element (and empty) slice must be left untouched.
	one := []float64{42}
	Shuffle(one, rand.New(rand.NewSource(2)))
	if len(one) != 1 || one[0] != 42 {
		t.Fatalf("Shuffle mishandled singleton: %v", one)
	}
	Shuffle(nil, rand.New(rand.NewSource(2)))
}

// TestResample verifies bootstrap draws keep the length, draw only from the
// input, and are deterministic for a fixed seed.
func TestResample(t *testing.T) {
	xs := []float64{2, 4, 6, 8}
	set := map[float64]bool{2: true, 4: true, 6: true, 8: true}

	out := Resample(xs, rand.New(rand.NewSource(9)))
	if len(out) != len(xs) {
		t.Fatalf("Resample length=%d want %d", len(out), len(xs))
	}
	for _, v := range out {
		if !set[v] {
			t.Fatalf("Resample produced foreign value %v", v)
		}
	}

	out2 := Resample(xs, rand.New(rand.NewSource(9)))
	for i := range out {
		if out[i] != out2[i] {
			t.Fatalf("Resample not deterministic at %d", i)
		}
	}

	if empty := Resample(nil, rand.New(rand.NewSource(1))); empty == nil || len(empty) != 0 {
		t.Fatalf("Resample(nil) = %v, want empty non-nil slice", empty)
	}
}

// TestSamplerInterfaces confirms the concrete types satisfy the sampler
// interfaces (compile-time and via the batch helpers).
func TestSamplerInterfaces(t *testing.T) {
	var _ ContinuousSampler = Uniform{}
	var _ ContinuousSampler = Normal{}
	var _ ContinuousSampler = FDist{}
	var _ DiscreteSampler = Poisson{}
	var _ DiscreteSampler = Hypergeometric{}

	r := rand.New(rand.NewSource(1))
	if got := SampleN(Uniform{A: 0, B: 1}, 3, r); len(got) != 3 {
		t.Fatalf("interface SampleN length=%d", len(got))
	}
	if got := SampleIntN(Binomial{N: 5, P: 0.5}, 3, r); len(got) != 3 {
		t.Fatalf("interface SampleIntN length=%d", len(got))
	}
}

func BenchmarkNormalSample(b *testing.B) {
	r := rand.New(rand.NewSource(1))
	d := Normal{Mu: 0, Sigma: 1}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = d.Sample(r)
	}
}

func BenchmarkGammaSample(b *testing.B) {
	r := rand.New(rand.NewSource(1))
	d := Gamma{Shape: 2.5, Scale: 1}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = d.Sample(r)
	}
}

func BenchmarkPoissonSampleSmall(b *testing.B) {
	r := rand.New(rand.NewSource(1))
	d := Poisson{Lambda: 4}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = d.SampleInt(r)
	}
}

func BenchmarkPoissonSampleLarge(b *testing.B) {
	r := rand.New(rand.NewSource(1))
	d := Poisson{Lambda: 200}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = d.SampleInt(r)
	}
}

func BenchmarkBinomialSample(b *testing.B) {
	r := rand.New(rand.NewSource(1))
	d := Binomial{N: 40, P: 0.3}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = d.SampleInt(r)
	}
}

func BenchmarkSampleN(b *testing.B) {
	r := rand.New(rand.NewSource(1))
	d := Normal{Mu: 0, Sigma: 1}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = SampleN(d, 256, r)
	}
}

func BenchmarkResample(b *testing.B) {
	r := rand.New(rand.NewSource(1))
	xs := make([]float64, 1024)
	for i := range xs {
		xs[i] = float64(i)
	}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = Resample(xs, r)
	}
}

func BenchmarkShuffle(b *testing.B) {
	r := rand.New(rand.NewSource(1))
	xs := make([]float64, 1024)
	for i := range xs {
		xs[i] = float64(i)
	}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Shuffle(xs, r)
	}
}
