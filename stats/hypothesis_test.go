package stats

import (
	"math"
	"testing"
)

// pTol is the tolerance used when comparing p-values against reference values
// computed with an independent implementation (SciPy). It is loose enough to
// absorb small differences between special-function implementations while still
// pinning each answer to four significant figures.
const pTol = 1e-4

func TestOneSampleTTest(t *testing.T) {
	// xs = 1..5: mean 3, sample variance 2.5, s = sqrt(2.5), se = s/sqrt(5).
	// t = (3-0)/se = 4.242640687..., df = 4.
	xs := []float64{1, 2, 3, 4, 5}
	r := OneSampleTTest(xs, 0)
	approx(t, r.Statistic, 4.242640687119285, 1e-9, "onesample t")
	approx(t, r.DF, 4, 0, "onesample df")
	approx(t, r.PValue, 0.01323559956368269, pTol, "onesample p")

	// mu0 equal to the mean gives t = 0 and p = 1.
	r0 := OneSampleTTest(xs, 3)
	approx(t, r0.Statistic, 0, 1e-12, "onesample t zero")
	approx(t, r0.PValue, 1, pTol, "onesample p one")

	// Degenerate: too few points and zero variance.
	if !math.IsNaN(OneSampleTTest([]float64{1}, 0).Statistic) {
		t.Fatalf("onesample n<2: want NaN")
	}
	if !math.IsNaN(OneSampleTTest([]float64{2, 2, 2}, 0).PValue) {
		t.Fatalf("onesample zero variance: want NaN")
	}
}

func TestTwoSampleTTest(t *testing.T) {
	// Two well-separated groups with identical variance 2.5 and equal size 5.
	// Pooled se = 1, t = (3-8)/1 = -5, df = 8.
	xs := []float64{1, 2, 3, 4, 5}
	ys := []float64{6, 7, 8, 9, 10}

	pooled := TwoSampleTTest(xs, ys, true)
	approx(t, pooled.Statistic, -5, 1e-9, "pooled t")
	approx(t, pooled.DF, 8, 0, "pooled df")
	approx(t, pooled.PValue, 0.0010528257933665399, pTol, "pooled p")

	// Welch: with equal sizes and equal variances df collapses to 8 as well.
	welch := TwoSampleTTest(xs, ys, false)
	approx(t, welch.Statistic, -5, 1e-9, "welch t")
	approx(t, welch.DF, 8, 1e-9, "welch df")
	approx(t, welch.PValue, 0.0010528257933665399, pTol, "welch p")

	if !math.IsNaN(TwoSampleTTest([]float64{1}, ys, true).Statistic) {
		t.Fatalf("twosample n<2: want NaN")
	}
}

func TestPairedTTest(t *testing.T) {
	// Differences are 1..5, reproducing the one-sample case above.
	xs := []float64{2, 4, 6, 8, 10}
	ys := []float64{1, 2, 3, 4, 5}
	r := PairedTTest(xs, ys)
	approx(t, r.Statistic, 4.242640687119285, 1e-9, "paired t")
	approx(t, r.DF, 4, 0, "paired df")
	approx(t, r.PValue, 0.01323559956368269, pTol, "paired p")

	if !math.IsNaN(PairedTTest([]float64{1, 2}, []float64{1}).Statistic) {
		t.Fatalf("paired mismatched length: want NaN")
	}
}

func TestChiSquareGoodnessOfFit(t *testing.T) {
	// Observed 10,20,30,40 vs uniform expected 25 each:
	// X2 = (225+25+25+225)/25 = 20, df = 3.
	obs := []float64{10, 20, 30, 40}
	exp := []float64{25, 25, 25, 25}
	r := ChiSquareGoodnessOfFit(obs, exp)
	approx(t, r.Statistic, 20, 1e-12, "gof X2")
	approx(t, r.DF, 3, 0, "gof df")
	approx(t, r.PValue, 0.00016974243555282632, pTol, "gof p")

	if !math.IsNaN(ChiSquareGoodnessOfFit(obs, []float64{25, 0, 25, 25}).Statistic) {
		t.Fatalf("gof expected<=0: want NaN")
	}
	if !math.IsNaN(ChiSquareGoodnessOfFit(obs, []float64{25, 25}).Statistic) {
		t.Fatalf("gof length mismatch: want NaN")
	}
}

func TestChiSquareIndependence(t *testing.T) {
	// 2x2 table with row totals 30,70 and column totals 40,60, total 100.
	// X2 = 4/12 + 4/18 + 4/28 + 4/42 = 0.79365..., df = 1.
	table := [][]float64{{10, 20}, {30, 40}}
	r := ChiSquareIndependence(table)
	approx(t, r.Statistic, 0.7936507936507936, 1e-12, "indep X2")
	approx(t, r.DF, 1, 0, "indep df")
	approx(t, r.PValue, 0.37299848361348686, pTol, "indep p")

	if !math.IsNaN(ChiSquareIndependence([][]float64{{1, 2}, {3}}).Statistic) {
		t.Fatalf("indep ragged: want NaN")
	}
	if !math.IsNaN(ChiSquareIndependence([][]float64{{0, 0}, {0, 0}}).Statistic) {
		t.Fatalf("indep zero total: want NaN")
	}
}

func TestOneWayANOVA(t *testing.T) {
	// Three groups with equal spread; SSB = 54, SSW = 6, F = 27, df 2 and 6.
	g1 := []float64{1, 2, 3}
	g2 := []float64{4, 5, 6}
	g3 := []float64{7, 8, 9}
	r := OneWayANOVA(g1, g2, g3)
	approx(t, r.Statistic, 27, 1e-9, "anova F")
	approx(t, r.DF, 2, 0, "anova dfB")
	approx(t, r.PValue, 0.0010000000000000002, pTol, "anova p")

	if !math.IsNaN(OneWayANOVA(g1).Statistic) {
		t.Fatalf("anova single group: want NaN")
	}
	if !math.IsNaN(OneWayANOVA([]float64{5, 5}, []float64{5, 5}).Statistic) {
		t.Fatalf("anova zero within variance: want NaN")
	}
}

func TestMannWhitneyU(t *testing.T) {
	// Complete separation: xs all below ys, no ties.
	// R1 = 10, U1 = 0, U = 0, z = (8-0.5)/sqrt(12) = 2.16506..., p = 0.03038...
	xs := []float64{1, 2, 3, 4}
	ys := []float64{5, 6, 7, 8}
	r := MannWhitneyU(xs, ys)
	approx(t, r.Statistic, 0, 1e-12, "mwu U")
	if !math.IsNaN(r.DF) {
		t.Fatalf("mwu df: want NaN, got %v", r.DF)
	}
	approx(t, r.PValue, 0.03038282197657749, pTol, "mwu p")

	// Symmetry: swapping the samples yields the same U and p-value.
	rs := MannWhitneyU(ys, xs)
	approx(t, rs.Statistic, r.Statistic, 1e-12, "mwu U symmetric")
	approx(t, rs.PValue, r.PValue, 1e-12, "mwu p symmetric")

	// Determinism with ties: repeated calls give identical results and the
	// tie-corrected path stays finite.
	xt := []float64{1, 2, 2, 3, 5}
	yt := []float64{2, 3, 4, 4, 6}
	a := MannWhitneyU(xt, yt)
	b := MannWhitneyU(xt, yt)
	// DF is always NaN for this test, so compare the finite fields directly
	// rather than the whole struct (NaN != NaN would make == always fail).
	if a.Statistic != b.Statistic || a.PValue != b.PValue {
		t.Fatalf("mwu not deterministic: %+v vs %+v", a, b)
	}
	if math.IsNaN(a.PValue) {
		t.Fatalf("mwu tie path: want finite p, got NaN")
	}

	if !math.IsNaN(MannWhitneyU(nil, ys).Statistic) {
		t.Fatalf("mwu empty sample: want NaN")
	}
}

func BenchmarkMannWhitneyU(b *testing.B) {
	n := 1000
	xs := make([]float64, n)
	ys := make([]float64, n)
	for i := 0; i < n; i++ {
		xs[i] = float64(i % 97)
		ys[i] = float64((i*7)%101) + 0.5
	}
	b.ReportAllocs()
	b.ResetTimer()
	var sink float64
	for i := 0; i < b.N; i++ {
		sink += MannWhitneyU(xs, ys).PValue
	}
	_ = sink
}

func BenchmarkOneWayANOVA(b *testing.B) {
	groups := make([][]float64, 20)
	for g := range groups {
		col := make([]float64, 500)
		for i := range col {
			col[i] = float64((i*3+g)%89) + float64(g)
		}
		groups[g] = col
	}
	b.ReportAllocs()
	b.ResetTimer()
	var sink float64
	for i := 0; i < b.N; i++ {
		sink += OneWayANOVA(groups...).Statistic
	}
	_ = sink
}
