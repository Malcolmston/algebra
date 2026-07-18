package stats

import (
	"math"
	"testing"
)

// statsApproxEqualOnline reports whether a and b are equal within tol, treating two
// NaNs as equal so undefined-result cases can be asserted directly.
func statsApproxEqualOnline(a, b, tol float64) bool {
	if math.IsNaN(a) || math.IsNaN(b) {
		return math.IsNaN(a) && math.IsNaN(b)
	}
	return math.Abs(a-b) <= tol
}

func TestAccumulatorMatchesBatch(t *testing.T) {
	const tol = 1e-9
	cases := [][]float64{
		{},
		{42},
		{2, 4},
		{1, 2, 3, 4, 5},
		{-3, -1, 0, 1, 3, 100},
		{2.5, 2.5, 2.5, 2.5},
		{1000000, 1000001, 1000002, 1000003},
	}
	for _, xs := range cases {
		var a Accumulator
		a.PushAll(xs)

		if got, want := a.Count(), len(xs); got != want {
			t.Errorf("Count(%v) = %d, want %d", xs, got, want)
		}
		if got, want := a.Sum(), Sum(xs); !statsApproxEqualOnline(got, want, tol) {
			t.Errorf("Sum(%v) = %v, want %v", xs, got, want)
		}
		if got, want := a.Mean(), Mean(xs); !statsApproxEqualOnline(got, want, tol) {
			t.Errorf("Mean(%v) = %v, want %v", xs, got, want)
		}
		if got, want := a.Variance(), Variance(xs); !statsApproxEqualOnline(got, want, tol) {
			t.Errorf("Variance(%v) = %v, want %v", xs, got, want)
		}
		if got, want := a.PopVariance(), PopVariance(xs); !statsApproxEqualOnline(got, want, tol) {
			t.Errorf("PopVariance(%v) = %v, want %v", xs, got, want)
		}
		if got, want := a.StdDev(), StdDev(xs); !statsApproxEqualOnline(got, want, tol) {
			t.Errorf("StdDev(%v) = %v, want %v", xs, got, want)
		}
		if got, want := a.PopStdDev(), PopStdDev(xs); !statsApproxEqualOnline(got, want, tol) {
			t.Errorf("PopStdDev(%v) = %v, want %v", xs, got, want)
		}
		if got, want := a.Skewness(), Skewness(xs); !statsApproxEqualOnline(got, want, tol) {
			t.Errorf("Skewness(%v) = %v, want %v", xs, got, want)
		}
		if got, want := a.Kurtosis(), Kurtosis(xs); !statsApproxEqualOnline(got, want, tol) {
			t.Errorf("Kurtosis(%v) = %v, want %v", xs, got, want)
		}
		if got, want := a.Min(), Min(xs); !statsApproxEqualOnline(got, want, tol) {
			t.Errorf("Min(%v) = %v, want %v", xs, got, want)
		}
		if got, want := a.Max(), Max(xs); !statsApproxEqualOnline(got, want, tol) {
			t.Errorf("Max(%v) = %v, want %v", xs, got, want)
		}
	}
}

func TestAccumulatorKnownAnswers(t *testing.T) {
	const tol = 1e-12
	var a Accumulator
	a.PushAll([]float64{2, 4, 4, 4, 5, 5, 7, 9})

	// Hand-computed reference values for this classic dataset.
	checks := []struct {
		name string
		got  float64
		want float64
	}{
		{"Count", float64(a.Count()), 8},
		{"Sum", a.Sum(), 40},
		{"Mean", a.Mean(), 5},
		{"PopVariance", a.PopVariance(), 4},  // M2=32, /8
		{"Variance", a.Variance(), 32.0 / 7}, // /7
		{"PopStdDev", a.PopStdDev(), 2},      // sqrt(4)
		{"Min", a.Min(), 2},
		{"Max", a.Max(), 9},
	}
	for _, c := range checks {
		if !statsApproxEqualOnline(c.got, c.want, tol) {
			t.Errorf("%s = %v, want %v", c.name, c.got, c.want)
		}
	}
}

func TestAccumulatorEmptyAndSingle(t *testing.T) {
	var a Accumulator
	if a.Count() != 0 {
		t.Errorf("empty Count = %d, want 0", a.Count())
	}
	if a.Sum() != 0 {
		t.Errorf("empty Sum = %v, want 0", a.Sum())
	}
	for name, got := range map[string]float64{
		"Mean":        a.Mean(),
		"Variance":    a.Variance(),
		"PopVariance": a.PopVariance(),
		"Skewness":    a.Skewness(),
		"Kurtosis":    a.Kurtosis(),
		"Min":         a.Min(),
		"Max":         a.Max(),
	} {
		if !math.IsNaN(got) {
			t.Errorf("empty %s = %v, want NaN", name, got)
		}
	}

	a.Push(7)
	if got := a.Mean(); got != 7 {
		t.Errorf("single Mean = %v, want 7", got)
	}
	if got := a.PopVariance(); got != 0 {
		t.Errorf("single PopVariance = %v, want 0", got)
	}
	if got := a.Variance(); !math.IsNaN(got) {
		t.Errorf("single Variance = %v, want NaN", got)
	}
}

func TestAccumulatorMerge(t *testing.T) {
	const tol = 1e-9
	full := []float64{-3, -1, 0, 1, 3, 100, 2, 4, 4, 4, 5, 5, 7, 9, 2.5, 8.1}

	// Split into three uneven shards, summarize each, then merge.
	shards := [][]float64{full[:3], full[3:9], full[9:]}
	var merged Accumulator
	for _, s := range shards {
		var part Accumulator
		part.PushAll(s)
		merged.Merge(part)
	}

	var single Accumulator
	single.PushAll(full)

	if merged.Count() != single.Count() {
		t.Fatalf("merged Count = %d, want %d", merged.Count(), single.Count())
	}
	for name, pair := range map[string][2]float64{
		"Mean":        {merged.Mean(), single.Mean()},
		"Variance":    {merged.Variance(), single.Variance()},
		"PopVariance": {merged.PopVariance(), single.PopVariance()},
		"Skewness":    {merged.Skewness(), single.Skewness()},
		"Kurtosis":    {merged.Kurtosis(), single.Kurtosis()},
		"Min":         {merged.Min(), single.Min()},
		"Max":         {merged.Max(), single.Max()},
	} {
		if !statsApproxEqualOnline(pair[0], pair[1], tol) {
			t.Errorf("merged %s = %v, want %v", name, pair[0], pair[1])
		}
	}
}

func TestAccumulatorMergeEmptyIdentity(t *testing.T) {
	var a Accumulator
	a.PushAll([]float64{1, 2, 3})
	before := a

	var empty Accumulator
	a.Merge(empty)
	if a != before {
		t.Errorf("Merge(empty) changed accumulator: %+v vs %+v", a, before)
	}

	var fresh Accumulator
	fresh.Merge(a)
	if fresh != a {
		t.Errorf("empty.Merge(a) = %+v, want %+v", fresh, a)
	}
}

func TestAccumulatorReset(t *testing.T) {
	var a Accumulator
	a.PushAll([]float64{1, 2, 3, 4})
	a.Reset()
	if a.Count() != 0 {
		t.Errorf("after Reset Count = %d, want 0", a.Count())
	}
	if !math.IsNaN(a.Mean()) {
		t.Errorf("after Reset Mean = %v, want NaN", a.Mean())
	}
}

func TestCovAccumulatorMatchesBatch(t *testing.T) {
	const tol = 1e-9
	type pair struct{ xs, ys []float64 }
	cases := []pair{
		{[]float64{1, 2, 3, 4, 5}, []float64{2, 4, 6, 8, 10}}, // perfect positive
		{[]float64{1, 2, 3, 4, 5}, []float64{5, 4, 3, 2, 1}},  // perfect negative
		{[]float64{1, 2, 3, 4, 5}, []float64{2, 1, 4, 3, 7}},  // noisy
		{[]float64{10, 20, 30}, []float64{7, 7, 7}},           // constant y
	}
	for _, c := range cases {
		var a CovAccumulator
		a.PushAll(c.xs, c.ys)

		if got, want := a.Count(), len(c.xs); got != want {
			t.Errorf("Count = %d, want %d", got, want)
		}
		if got, want := a.Covariance(), Covariance(c.xs, c.ys); !statsApproxEqualOnline(got, want, tol) {
			t.Errorf("Covariance(%v,%v) = %v, want %v", c.xs, c.ys, got, want)
		}
		if got, want := a.Correlation(), Correlation(c.xs, c.ys); !statsApproxEqualOnline(got, want, tol) {
			t.Errorf("Correlation(%v,%v) = %v, want %v", c.xs, c.ys, got, want)
		}
		slope, intercept, _ := LinearRegression(c.xs, c.ys)
		if got := a.Slope(); !statsApproxEqualOnline(got, slope, tol) {
			t.Errorf("Slope(%v,%v) = %v, want %v", c.xs, c.ys, got, slope)
		}
		if got := a.Intercept(); !statsApproxEqualOnline(got, intercept, tol) {
			t.Errorf("Intercept(%v,%v) = %v, want %v", c.xs, c.ys, got, intercept)
		}
	}
}

func TestCovAccumulatorKnownAnswer(t *testing.T) {
	const tol = 1e-12
	var a CovAccumulator
	a.PushAll([]float64{1, 2, 3, 4, 5}, []float64{2, 4, 6, 8, 10})
	if got := a.Slope(); !statsApproxEqualOnline(got, 2, tol) {
		t.Errorf("Slope = %v, want 2", got)
	}
	if got := a.Intercept(); !statsApproxEqualOnline(got, 0, tol) {
		t.Errorf("Intercept = %v, want 0", got)
	}
	if got := a.Correlation(); !statsApproxEqualOnline(got, 1, tol) {
		t.Errorf("Correlation = %v, want 1", got)
	}
	if got := a.Covariance(); !statsApproxEqualOnline(got, 5, tol) {
		t.Errorf("Covariance = %v, want 5", got)
	}
}

func TestCovAccumulatorEmpty(t *testing.T) {
	var a CovAccumulator
	for name, got := range map[string]float64{
		"Covariance":  a.Covariance(),
		"Correlation": a.Correlation(),
		"Slope":       a.Slope(),
		"Intercept":   a.Intercept(),
	} {
		if !math.IsNaN(got) {
			t.Errorf("empty %s = %v, want NaN", name, got)
		}
	}
}

func TestCovAccumulatorMerge(t *testing.T) {
	const tol = 1e-9
	xs := []float64{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}
	ys := []float64{2, 1, 4, 3, 7, 5, 9, 8, 12, 10}

	var merged CovAccumulator
	splits := [][2]int{{0, 3}, {3, 7}, {7, 10}}
	for _, s := range splits {
		var part CovAccumulator
		part.PushAll(xs[s[0]:s[1]], ys[s[0]:s[1]])
		merged.Merge(part)
	}

	var single CovAccumulator
	single.PushAll(xs, ys)

	for name, pair := range map[string][2]float64{
		"Covariance":  {merged.Covariance(), single.Covariance()},
		"Correlation": {merged.Correlation(), single.Correlation()},
		"Slope":       {merged.Slope(), single.Slope()},
		"Intercept":   {merged.Intercept(), single.Intercept()},
	} {
		if !statsApproxEqualOnline(pair[0], pair[1], tol) {
			t.Errorf("merged %s = %v, want %v", name, pair[0], pair[1])
		}
	}
}

func TestCovAccumulatorReset(t *testing.T) {
	var a CovAccumulator
	a.PushAll([]float64{1, 2, 3}, []float64{4, 5, 6})
	a.Reset()
	if a.Count() != 0 {
		t.Errorf("after Reset Count = %d, want 0", a.Count())
	}
	if !math.IsNaN(a.Covariance()) {
		t.Errorf("after Reset Covariance = %v, want NaN", a.Covariance())
	}
}

// statsBenchData is a fixed deterministic dataset for the benchmarks so runs are
// comparable.
var statsBenchData = func() []float64 {
	xs := make([]float64, 10000)
	for i := range xs {
		xs[i] = math.Sin(float64(i)) * float64(i%97)
	}
	return xs
}()

func BenchmarkAccumulatorPushAll(b *testing.B) {
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		var a Accumulator
		a.PushAll(statsBenchData)
		_ = a.Variance()
	}
}

func BenchmarkAccumulatorMoments(b *testing.B) {
	var a Accumulator
	a.PushAll(statsBenchData)
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = a.Skewness()
		_ = a.Kurtosis()
	}
}

func BenchmarkCovAccumulatorPushAll(b *testing.B) {
	ys := make([]float64, len(statsBenchData))
	for i, x := range statsBenchData {
		ys[i] = x*1.5 + 3
	}
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		var a CovAccumulator
		a.PushAll(statsBenchData, ys)
		_ = a.Slope()
	}
}
