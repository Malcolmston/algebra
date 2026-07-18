package stats

import (
	"math"
	"testing"
)

func TestGeometricFailuresDist(t *testing.T) {
	g := GeometricFailures{P: 0.5}

	// PMF known answers: (1-P)^k * P.
	pmf := []struct {
		k    int
		want float64
	}{
		{-1, 0},
		{0, 0.5},
		{1, 0.25},
		{2, 0.125},
		{3, 0.0625},
	}
	for _, c := range pmf {
		approx(t, g.PMF(c.k), c.want, tol, "geomfail pmf")
	}

	// CDF known answers: 1 - (1-P)^(k+1).
	cdf := []struct {
		k    int
		want float64
	}{
		{-1, 0},
		{0, 0.5},
		{1, 0.75},
		{2, 0.875},
		{3, 0.9375},
	}
	for _, c := range cdf {
		approx(t, g.CDF(c.k), c.want, tol, "geomfail cdf")
	}

	// Quantile is the smallest k with CDF(k) >= p.
	q := []struct {
		p    float64
		want int
	}{
		{0, 0},
		{0.5, 0},
		{0.6, 1},
		{0.75, 1},
		{0.8, 2},
		{0.9, 3},
	}
	for _, c := range q {
		if got := g.Quantile(c.p); got != c.want {
			t.Errorf("geomfail quantile(%v) = %d, want %d", c.p, got, c.want)
		}
	}

	approx(t, g.Mean(), 1, tol, "geomfail mean")
	approx(t, g.Variance(), 2, tol, "geomfail var")

	// A different P: mean (1-P)/P, variance (1-P)/P^2.
	g2 := GeometricFailures{P: 0.25}
	approx(t, g2.Mean(), 3, tol, "geomfail mean p=.25")
	approx(t, g2.Variance(), 12, tol, "geomfail var p=.25")
	approx(t, g2.PMF(0), 0.25, tol, "geomfail pmf0 p=.25")
	approx(t, g2.CDF(0), 0.25, tol, "geomfail cdf0 p=.25")

	// Invalid parameters and out-of-range p.
	approx(t, GeometricFailures{P: 0}.PMF(0), math.NaN(), 0, "geomfail invalid pmf")
	approx(t, GeometricFailures{P: 1.5}.Mean(), math.NaN(), 0, "geomfail invalid mean")
	if got := g.Quantile(-0.1); got != 0 {
		t.Errorf("geomfail quantile(-0.1) = %d, want 0", got)
	}
	if got := g.Quantile(1); got != math.MaxInt {
		t.Errorf("geomfail quantile(1) = %d, want MaxInt", got)
	}
	// P == 1: all mass at 0.
	approx(t, GeometricFailures{P: 1}.PMF(0), 1, tol, "geomfail p=1 pmf0")
	if got := (GeometricFailures{P: 1}).Quantile(1); got != 0 {
		t.Errorf("geomfail p=1 quantile(1) = %d, want 0", got)
	}
}

func TestNegativeBinomialIntDist(t *testing.T) {
	nb := NegativeBinomialInt{R: 3, P: 0.5}

	// PMF: exp(gammaLn(k+R)-gammaLn(k+1)-gammaLn(R)+R ln P + k ln(1-P)).
	approx(t, nb.PMF(-1), 0, 0, "nbi pmf neg")
	approx(t, nb.PMF(0), 0.125, tol, "nbi pmf0")  // P^R = 0.5^3
	approx(t, nb.PMF(1), 0.1875, tol, "nbi pmf1") // C(3,1) P^3 (1-P)
	approx(t, nb.PMF(2), 0.1875, tol, "nbi pmf2") // C(4,2) P^3 (1-P)^2

	// CDF via regularized incomplete beta matches cumulative PMF.
	var sum float64
	for k := 0; k <= 5; k++ {
		sum += nb.PMF(k)
		approx(t, nb.CDF(k), sum, 1e-9, "nbi cdf==sum")
	}
	approx(t, nb.CDF(0), 0.125, tol, "nbi cdf0")
	approx(t, nb.CDF(1), 0.3125, tol, "nbi cdf1")

	// Quantile: smallest k with CDF(k) >= p.
	q := []struct {
		p    float64
		want int
	}{
		{0, 0},
		{0.125, 0},
		{0.2, 1},
		{0.3125, 1},
		{0.4, 2},
	}
	for _, c := range q {
		if got := nb.Quantile(c.p); got != c.want {
			t.Errorf("nbi quantile(%v) = %d, want %d", c.p, got, c.want)
		}
	}

	approx(t, nb.Mean(), 3, tol, "nbi mean")    // R(1-P)/P
	approx(t, nb.Variance(), 6, tol, "nbi var") // R(1-P)/P^2

	// R == 1 reduces to the failures-geometric.
	nb1 := NegativeBinomialInt{R: 1, P: 0.5}
	g := GeometricFailures{P: 0.5}
	for k := 0; k <= 4; k++ {
		approx(t, nb1.PMF(k), g.PMF(k), 1e-12, "nbi R=1 == geom")
	}

	// Invalid parameters and out-of-range p.
	approx(t, NegativeBinomialInt{R: 0, P: 0.5}.PMF(0), math.NaN(), 0, "nbi invalid R")
	approx(t, NegativeBinomialInt{R: 2, P: 0}.Mean(), math.NaN(), 0, "nbi invalid P")
	if got := nb.Quantile(2); got != 0 {
		t.Errorf("nbi quantile(2) = %d, want 0", got)
	}
	if got := nb.Quantile(1); got != math.MaxInt {
		t.Errorf("nbi quantile(1) = %d, want MaxInt", got)
	}
	// P == 1: certainly 0 failures.
	approx(t, NegativeBinomialInt{R: 3, P: 1}.PMF(0), 1, tol, "nbi p=1 pmf0")
	approx(t, NegativeBinomialInt{R: 3, P: 1}.PMF(1), 0, tol, "nbi p=1 pmf1")
}

func TestHypergeometricDist(t *testing.T) {
	// N=4, K=2, Draws=2: support 0..2.
	h := Hypergeometric{N: 4, K: 2, Draws: 2}
	approx(t, h.PMF(-1), 0, 0, "hyper pmf neg")
	approx(t, h.PMF(0), 1.0/6, tol, "hyper pmf0")
	approx(t, h.PMF(1), 4.0/6, tol, "hyper pmf1")
	approx(t, h.PMF(2), 1.0/6, tol, "hyper pmf2")
	approx(t, h.PMF(3), 0, 0, "hyper pmf out")

	approx(t, h.CDF(-1), 0, 0, "hyper cdf neg")
	approx(t, h.CDF(0), 1.0/6, tol, "hyper cdf0")
	approx(t, h.CDF(1), 5.0/6, tol, "hyper cdf1")
	approx(t, h.CDF(2), 1, tol, "hyper cdf2")
	approx(t, h.CDF(5), 1, tol, "hyper cdf above")

	approx(t, h.Mean(), 1, tol, "hyper mean")        // Draws*K/N = 2*2/4
	approx(t, h.Variance(), 1.0/3, tol, "hyper var") // 2*.5*.5*(2)/(3)

	q := []struct {
		p    float64
		want int
	}{
		{0, 0},
		{0.1, 0},
		{0.5, 1},
		{0.9, 2},
		{1, 2},
	}
	for _, c := range q {
		if got := h.Quantile(c.p); got != c.want {
			t.Errorf("hyper quantile(%v) = %d, want %d", c.p, got, c.want)
		}
	}

	// A distribution with a nonzero lower support bound: N=10, K=8, Draws=5
	// gives support 3..5. Quantile search must start at the lower bound.
	h2 := Hypergeometric{N: 10, K: 8, Draws: 5}
	approx(t, h2.PMF(2), 0, 0, "hyper2 below support")
	approx(t, h2.PMF(3), 56.0/252, tol, "hyper2 pmf3")
	approx(t, h2.PMF(4), 140.0/252, tol, "hyper2 pmf4")
	approx(t, h2.PMF(5), 56.0/252, tol, "hyper2 pmf5")
	approx(t, h2.CDF(2), 0, 0, "hyper2 cdf2")
	approx(t, h2.CDF(3), 56.0/252, tol, "hyper2 cdf3")
	if got := h2.Quantile(0.1); got != 3 {
		t.Errorf("hyper2 quantile(0.1) = %d, want 3", got)
	}
	if got := h2.Quantile(0); got != 3 {
		t.Errorf("hyper2 quantile(0) = %d, want 3", got)
	}

	// Invalid parameters and out-of-range p.
	approx(t, Hypergeometric{N: 5, K: 6, Draws: 2}.PMF(1), math.NaN(), 0, "hyper invalid K")
	approx(t, Hypergeometric{N: 5, K: 2, Draws: 9}.Mean(), math.NaN(), 0, "hyper invalid Draws")
	if got := h.Quantile(1.5); got != 0 {
		t.Errorf("hyper quantile(1.5) = %d, want 0", got)
	}
	// N == 1 has no variability.
	approx(t, Hypergeometric{N: 1, K: 1, Draws: 1}.Variance(), 0, tol, "hyper N=1 var")
}

func BenchmarkNegativeBinomialIntCDF(b *testing.B) {
	nb := NegativeBinomialInt{R: 10, P: 0.3}
	var sink float64
	for i := 0; i < b.N; i++ {
		sink += nb.CDF(50)
	}
	_ = sink
}

func BenchmarkNegativeBinomialIntQuantile(b *testing.B) {
	nb := NegativeBinomialInt{R: 10, P: 0.3}
	var sink int
	for i := 0; i < b.N; i++ {
		sink += nb.Quantile(0.95)
	}
	_ = sink
}

func BenchmarkHypergeometricQuantile(b *testing.B) {
	h := Hypergeometric{N: 500, K: 200, Draws: 100}
	var sink int
	for i := 0; i < b.N; i++ {
		sink += h.Quantile(0.9)
	}
	_ = sink
}
