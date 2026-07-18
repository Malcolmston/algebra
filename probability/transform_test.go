package probability

import "testing"

func TestAffineTransform(t *testing.T) {
	d := mustDist(DiscreteUniform(1, 6))
	// Y = 2X + 3: mean scales/shifts, variance scales by a^2.
	y := d.Affine(2, 3)
	if err := y.Validate(); err != nil {
		t.Fatal(err)
	}
	if !approx(y.Mean(), 2*d.Mean()+3, testTol) {
		t.Errorf("affine mean=%g want %g", y.Mean(), 2*d.Mean()+3)
	}
	if !approx(y.Variance(), 4*d.Variance(), testTol) {
		t.Errorf("affine var=%g want %g", y.Variance(), 4*d.Variance())
	}
}

func TestScaleAndShift(t *testing.T) {
	d := mustDist(Bernoulli(0.3))
	if !approx(d.Scale(10).Mean(), 3, testTol) {
		t.Errorf("scale mean=%g want 3", d.Scale(10).Mean())
	}
	if !approx(d.Shift(5).Mean(), 5.3, testTol) {
		t.Errorf("shift mean=%g want 5.3", d.Shift(5).Mean())
	}
}

func TestTransformMergesCollisions(t *testing.T) {
	// Y = X^2 maps -1 and 1 to the same outcome; they must merge.
	d := mustDist(NewDistribution([]float64{-1, 0, 1}, []float64{0.25, 0.5, 0.25}))
	y := d.Transform(func(x float64) float64 { return x * x })
	if y.Len() != 2 {
		t.Fatalf("expected 2 outcomes after squaring, got %d", y.Len())
	}
	if !approx(y.PMF(1), 0.5, testTol) {
		t.Errorf("PMF(1)=%g want 0.5", y.PMF(1))
	}
}

func TestStandardize(t *testing.T) {
	d := mustDist(DiscreteUniform(1, 6))
	z := mustDist(d.Standardize())
	if !approx(z.Mean(), 0, 1e-9) {
		t.Errorf("standardized mean=%g want 0", z.Mean())
	}
	if !approx(z.Variance(), 1, 1e-9) {
		t.Errorf("standardized var=%g want 1", z.Variance())
	}
}

func TestConvolveTwoDice(t *testing.T) {
	die := mustDist(DiscreteUniform(1, 6))
	sum := die.Convolve(die)
	if err := sum.Validate(); err != nil {
		t.Fatal(err)
	}
	if !approx(sum.PMF(2), 1.0/36.0, testTol) {
		t.Errorf("P(sum=2)=%g want 1/36", sum.PMF(2))
	}
	if !approx(sum.PMF(7), 6.0/36.0, testTol) {
		t.Errorf("P(sum=7)=%g want 6/36", sum.PMF(7))
	}
	if !approx(sum.Mean(), 7, testTol) {
		t.Errorf("sum mean=%g want 7", sum.Mean())
	}
	// Variance of a sum of independents is the sum of variances.
	if !approx(sum.Variance(), 2*die.Variance(), testTol) {
		t.Errorf("sum var=%g want %g", sum.Variance(), 2*die.Variance())
	}
}

func TestConvolvePowerEqualsBinomial(t *testing.T) {
	p := 0.3
	n := 10
	b := mustDist(Bernoulli(p))
	got, err := b.ConvolvePower(n)
	if err != nil {
		t.Fatal(err)
	}
	want := mustDist(Binomial(n, p))
	for k := 0; k <= n; k++ {
		if !approx(got.PMF(float64(k)), want.PMF(float64(k)), 1e-9) {
			t.Errorf("k=%d convolved=%g binomial=%g", k, got.PMF(float64(k)), want.PMF(float64(k)))
		}
	}
	// ConvolvePower(0) is a point mass at zero.
	z, _ := b.ConvolvePower(0)
	if z.Len() != 1 || z.Outcomes[0] != 0 || !approx(z.Probs[0], 1, testTol) {
		t.Errorf("ConvolvePower(0) is not a point mass at 0: %+v", z)
	}
}

func TestMixture(t *testing.T) {
	a := mustDist(Uniform([]float64{0}))
	b := mustDist(Uniform([]float64{10}))
	mix, err := Mixture([]float64{0.5, 0.5}, []Distribution{a, b})
	if err != nil {
		t.Fatal(err)
	}
	if !approx(mix.PMF(0), 0.5, testTol) || !approx(mix.PMF(10), 0.5, testTol) {
		t.Errorf("mixture probs wrong: %+v", mix)
	}
	if !approx(mix.Mean(), 5, testTol) {
		t.Errorf("mixture mean=%g want 5", mix.Mean())
	}
	if _, err := Mixture([]float64{0.5}, []Distribution{a, b}); err == nil {
		t.Error("expected length-mismatch error")
	}
}
