package stats

import (
	"math"
	"math/rand"
	"testing"
)

func TestMeanConfidenceInterval(t *testing.T) {
	xs := []float64{1, 2, 3, 4, 5}
	// mean=3, sd=1.581139, se=0.707107, t_{0.975,4}=2.776445, margin=1.963243.
	lo, hi := MeanConfidenceInterval(xs, 0.95)
	approx(t, (lo+hi)/2, 3, tol, "ci center")
	approx(t, lo, 1.0367568385, 1e-4, "ci lower")
	approx(t, hi, 4.9632431615, 1e-4, "ci upper")
}

func TestStudentTQuantile(t *testing.T) {
	// Known t critical values.
	approx(t, StudentT{Nu: 4}.Quantile(0.975), 2.7764451, 1e-5, "t 0.975,4")
	approx(t, StudentT{Nu: 10}.Quantile(0.95), 1.8124611, 1e-5, "t 0.95,10")
	// Large df approaches the normal quantile.
	approx(t, StudentT{Nu: 100000}.Quantile(0.975), 1.9599640, 1e-4, "t large df")
	// Boundary behaviour.
	if !math.IsInf(StudentT{Nu: 4}.Quantile(1), 1) {
		t.Errorf("Quantile(1) should be +Inf")
	}
	if !math.IsInf(StudentT{Nu: 4}.Quantile(0), -1) {
		t.Errorf("Quantile(0) should be -Inf")
	}
}

func TestMeanConfidenceIntervalZ(t *testing.T) {
	xs := []float64{1, 2, 3, 4, 5}
	// Known sigma=1.581139, mean=3, se=0.707107, z=1.959964, margin=1.385904.
	lo, hi := MeanConfidenceIntervalZ(xs, 1.5811388, 0.95)
	approx(t, (lo+hi)/2, 3, tol, "z ci center")
	approx(t, lo, 3-1.3859038, 1e-4, "z ci lower")
	approx(t, hi, 3+1.3859038, 1e-4, "z ci upper")

	// Guards: empty sample, non-positive sigma.
	if l, _ := MeanConfidenceIntervalZ(nil, 1, 0.95); !math.IsNaN(l) {
		t.Error("empty sample should be NaN")
	}
	if l, _ := MeanConfidenceIntervalZ(xs, 0, 0.95); !math.IsNaN(l) {
		t.Error("sigma<=0 should be NaN")
	}
}

func TestProportionConfidenceInterval(t *testing.T) {
	// 40 successes in 100 trials, 95% Wilson score interval.
	lo, hi := ProportionConfidenceInterval(40, 100, 0.95)
	approx(t, lo, 0.3094012864, 1e-6, "prop lower")
	approx(t, hi, 0.4979974132, 1e-6, "prop upper")
	// The interval must bracket the point estimate 0.4.
	if !(lo < 0.4 && 0.4 < hi) {
		t.Errorf("interval (%.4f, %.4f) does not bracket 0.4", lo, hi)
	}

	// Guards.
	if l, _ := ProportionConfidenceInterval(5, 0, 0.95); !math.IsNaN(l) {
		t.Error("n=0 should be NaN")
	}
	if l, _ := ProportionConfidenceInterval(-1, 100, 0.95); !math.IsNaN(l) {
		t.Error("negative successes should be NaN")
	}
	// Bounds are clamped to [0,1] at the extremes.
	lo0, hi0 := ProportionConfidenceInterval(0, 10, 0.95)
	if lo0 < 0 || hi0 > 1 {
		t.Errorf("bounds not clamped: (%.4f, %.4f)", lo0, hi0)
	}
}

func TestBootstrapDeterministic(t *testing.T) {
	xs := []float64{2, 4, 4, 4, 5, 5, 7, 9}
	a := Bootstrap(xs, Mean, 1000, 0.95, rand.New(rand.NewSource(2024)))
	b := Bootstrap(xs, Mean, 1000, 0.95, rand.New(rand.NewSource(2024)))
	// Identical seeds must give identical results.
	if a != b {
		t.Fatalf("bootstrap not deterministic: %+v vs %+v", a, b)
	}
	// The point estimate is the statistic on the original sample.
	approx(t, a.Estimate, Mean(xs), tol, "bootstrap estimate")
	// The percentile interval brackets the estimate.
	if !(a.Lo <= a.Estimate && a.Estimate <= a.Hi) {
		t.Errorf("estimate %.3f not in interval (%.3f, %.3f)", a.Estimate, a.Lo, a.Hi)
	}
}

func TestBootstrapCoverage(t *testing.T) {
	r := rand.New(rand.NewSource(7))
	xs := make([]float64, 500)
	src := rand.New(rand.NewSource(1))
	for i := range xs {
		xs[i] = Normal{Mu: 10, Sigma: 2}.Sample(src)
	}
	res := Bootstrap(xs, Mean, 2000, 0.95, r)
	approx(t, res.Estimate, Mean(xs), tol, "bootci estimate")
	if !(res.Lo < res.Estimate && res.Estimate < res.Hi) {
		t.Errorf("estimate %.3f not inside (%.3f, %.3f)", res.Estimate, res.Lo, res.Hi)
	}
	// The 95% bootstrap CI for the mean should contain the true mean 10.
	if !(res.Lo < 10 && 10 < res.Hi) {
		t.Errorf("true mean 10 not in CI (%.3f, %.3f)", res.Lo, res.Hi)
	}
	// StdErr is positive for a non-degenerate sample.
	if !(res.StdErr > 0) {
		t.Errorf("StdErr = %v, want > 0", res.StdErr)
	}

	// BootstrapMean is a convenience wrapper equivalent to Bootstrap with Mean.
	rw := rand.New(rand.NewSource(7))
	rm := rand.New(rand.NewSource(7))
	want := Bootstrap(xs, Mean, 2000, 0.95, rw)
	got := BootstrapMean(xs, 2000, 0.95, rm)
	if want != got {
		t.Errorf("BootstrapMean %+v != Bootstrap %+v", got, want)
	}
}

func TestInferenceGuards(t *testing.T) {
	if lo, _ := MeanConfidenceInterval([]float64{1}, 0.95); !math.IsNaN(lo) {
		t.Error("CI with n<2 should be NaN")
	}
	if lo, _ := MeanConfidenceInterval([]float64{1, 2}, 1.5); !math.IsNaN(lo) {
		t.Error("CI with level>=1 should be NaN")
	}
	nan := Bootstrap(nil, Mean, 10, 0.95, rand.New(rand.NewSource(1)))
	if !math.IsNaN(nan.Estimate) {
		t.Error("bootstrap of empty should be NaN")
	}
	if r := Bootstrap([]float64{1, 2, 3}, nil, 10, 0.95, rand.New(rand.NewSource(1))); !math.IsNaN(r.Estimate) {
		t.Error("bootstrap with nil stat should be NaN")
	}
}

func BenchmarkBootstrap(b *testing.B) {
	xs := make([]float64, 200)
	for i := range xs {
		xs[i] = float64(i)
	}
	r := rand.New(rand.NewSource(1))
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = Bootstrap(xs, Mean, 500, 0.95, r)
	}
}
