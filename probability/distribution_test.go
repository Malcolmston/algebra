package probability

import (
	"math"
	"testing"
)

const testTol = 1e-9

func approx(a, b, tol float64) bool { return math.Abs(a-b) <= tol }

// mustDist unwraps a constructor result, panicking (and thus failing the test)
// if construction or validation fails. Taking exactly the (Distribution, error)
// pair lets it be called directly on a constructor, e.g. mustDist(Bernoulli(p)).
func mustDist(d Distribution, err error) Distribution {
	if err != nil {
		panic("probability test: unexpected error: " + err.Error())
	}
	if err := d.Validate(); err != nil {
		panic("probability test: distribution failed validation: " + err.Error())
	}
	return d
}

func TestNewDistributionErrors(t *testing.T) {
	if _, err := NewDistribution([]float64{1, 2}, []float64{0.5}); err == nil {
		t.Error("expected length-mismatch error")
	}
	if _, err := NewDistribution(nil, nil); err == nil {
		t.Error("expected empty-support error")
	}
	if _, err := NewDistribution([]float64{1, 2}, []float64{0.5, 0.4}); err == nil {
		t.Error("expected sum-not-one error")
	}
	if _, err := NewDistribution([]float64{1, 2}, []float64{-0.1, 1.1}); err == nil {
		t.Error("expected negative-probability error")
	}
}

func TestNewDistributionMergesDuplicates(t *testing.T) {
	d := mustDist(NewDistribution([]float64{2, 1, 2}, []float64{0.25, 0.5, 0.25}))
	if d.Len() != 2 {
		t.Fatalf("expected 2 distinct outcomes, got %d", d.Len())
	}
	if !approx(d.PMF(2), 0.5, testTol) {
		t.Errorf("merged PMF(2)=%g, want 0.5", d.PMF(2))
	}
	if d.Outcomes[0] != 1 || d.Outcomes[1] != 2 {
		t.Errorf("outcomes not sorted: %v", d.Outcomes)
	}
}

func TestBernoulliMoments(t *testing.T) {
	for _, p := range []float64{0.0, 0.2, 0.5, 0.8, 1.0} {
		d := mustDist(Bernoulli(p))
		if !approx(d.Mean(), p, testTol) {
			t.Errorf("Bernoulli(%g) mean=%g want %g", p, d.Mean(), p)
		}
		if !approx(d.Variance(), p*(1-p), testTol) {
			t.Errorf("Bernoulli(%g) var=%g want %g", p, d.Variance(), p*(1-p))
		}
	}
}

func TestBinomialClosedForm(t *testing.T) {
	n, p := 10, 0.3
	d := mustDist(Binomial(n, p))
	np := float64(n) * p
	wantMean := np
	wantVar := np * (1 - p)
	wantSkew := (1 - 2*p) / math.Sqrt(np*(1-p))
	wantKurt := (1 - 6*p*(1-p)) / (np * (1 - p))
	if !approx(d.Mean(), wantMean, testTol) {
		t.Errorf("mean=%g want %g", d.Mean(), wantMean)
	}
	if !approx(d.Variance(), wantVar, testTol) {
		t.Errorf("var=%g want %g", d.Variance(), wantVar)
	}
	if !approx(d.Skewness(), wantSkew, 1e-9) {
		t.Errorf("skew=%g want %g", d.Skewness(), wantSkew)
	}
	if !approx(d.Kurtosis(), wantKurt, 1e-9) {
		t.Errorf("excess kurtosis=%g want %g", d.Kurtosis(), wantKurt)
	}
	// PGF of a binomial is (1-p+pz)^n.
	for _, z := range []float64{0.5, 1.0, 1.5} {
		want := math.Pow(1-p+p*z, float64(n))
		if !approx(d.PGF(z), want, 1e-9) {
			t.Errorf("PGF(%g)=%g want %g", z, d.PGF(z), want)
		}
	}
}

func TestDiscreteUniformDie(t *testing.T) {
	d := mustDist(DiscreteUniform(1, 6))
	if !approx(d.Mean(), 3.5, testTol) {
		t.Errorf("die mean=%g want 3.5", d.Mean())
	}
	if !approx(d.Variance(), 35.0/12.0, testTol) {
		t.Errorf("die var=%g want %g", d.Variance(), 35.0/12.0)
	}
	if !approx(d.Median(), 3, testTol) && !approx(d.Median(), 4, testTol) {
		t.Errorf("die median=%g", d.Median())
	}
	if !approx(d.CDF(3), 0.5, testTol) {
		t.Errorf("die CDF(3)=%g want 0.5", d.CDF(3))
	}
	if !approx(d.Quantile(0.5), 3, testTol) {
		t.Errorf("die quantile(0.5)=%g want 3", d.Quantile(0.5))
	}
	if d.Min() != 1 || d.Max() != 6 {
		t.Errorf("die support [%g,%g] want [1,6]", d.Min(), d.Max())
	}
}

func TestPoissonMoments(t *testing.T) {
	lambda := 4.0
	d := mustDist(Poisson(lambda, 60))
	if !approx(d.Mean(), lambda, 1e-7) {
		t.Errorf("Poisson mean=%g want %g", d.Mean(), lambda)
	}
	if !approx(d.Variance(), lambda, 1e-7) {
		t.Errorf("Poisson var=%g want %g", d.Variance(), lambda)
	}
}

func TestGeometricMean(t *testing.T) {
	p := 0.25
	d := mustDist(Geometric(p, 400))
	if !approx(d.Mean(), 1/p, 1e-6) {
		t.Errorf("Geometric mean=%g want %g", d.Mean(), 1/p)
	}
	if !approx(d.Variance(), (1-p)/(p*p), 1e-5) {
		t.Errorf("Geometric var=%g want %g", d.Variance(), (1-p)/(p*p))
	}
}

func TestModeAndSupport(t *testing.T) {
	d := mustDist(NewDistribution([]float64{0, 1, 2}, []float64{0.2, 0.5, 0.3}))
	if d.Mode() != 1 {
		t.Errorf("mode=%g want 1", d.Mode())
	}
	sup := d.Support()
	sup[0] = 99 // ensure Support returns a copy
	if d.Outcomes[0] == 99 {
		t.Error("Support did not return a copy")
	}
}

func TestNormalize(t *testing.T) {
	d := Distribution{Outcomes: []float64{0, 1}, Probs: []float64{1, 3}}
	nd, err := d.Normalize()
	if err != nil {
		t.Fatal(err)
	}
	if !approx(nd.Probs[0], 0.25, testTol) || !approx(nd.Probs[1], 0.75, testTol) {
		t.Errorf("normalized probs=%v", nd.Probs)
	}
}
