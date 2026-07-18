package probability

import (
	"math"
	"math/cmplx"
	"testing"
)

func TestMGFBernoulli(t *testing.T) {
	p := 0.4
	d := mustDist(Bernoulli(p))
	for _, tt := range []float64{-1, 0, 0.5, 1} {
		want := 1 - p + p*math.Exp(tt)
		if !approx(d.MGF(tt), want, testTol) {
			t.Errorf("MGF(%g)=%g want %g", tt, d.MGF(tt), want)
		}
	}
	if !approx(d.MGF(0), 1, testTol) {
		t.Errorf("MGF(0)=%g want 1", d.MGF(0))
	}
	// CGF is the log of the MGF.
	if !approx(d.CGF(0.5), math.Log(d.MGF(0.5)), testTol) {
		t.Errorf("CGF mismatch")
	}
}

func TestMomentsRawAndCentral(t *testing.T) {
	d := mustDist(DiscreteUniform(1, 6))
	if !approx(d.Moment(0), 1, testTol) {
		t.Errorf("Moment(0)=%g want 1", d.Moment(0))
	}
	if !approx(d.Moment(1), d.Mean(), testTol) {
		t.Errorf("Moment(1)=%g want mean %g", d.Moment(1), d.Mean())
	}
	if !approx(d.CentralMoment(2), d.Variance(), testTol) {
		t.Errorf("CentralMoment(2)=%g want var %g", d.CentralMoment(2), d.Variance())
	}
	// E[X^2] = Var + Mean^2.
	if !approx(d.Moment(2), d.Variance()+d.Mean()*d.Mean(), testTol) {
		t.Errorf("Moment(2) identity failed")
	}
}

func TestMGFDerivativeMomentMatchesRaw(t *testing.T) {
	d := mustDist(Binomial(6, 0.4))
	for k := 1; k <= 3; k++ {
		got := d.MGFDerivativeMoment(k, 1e-3)
		want := d.Moment(k)
		if !approx(got, want, 1e-3) {
			t.Errorf("MGFDerivativeMoment(%d)=%g want ~%g", k, got, want)
		}
	}
}

func TestCharacteristicFunctionAtZero(t *testing.T) {
	d := mustDist(Binomial(5, 0.5))
	if v := d.CharacteristicFunction(0); cmplx.Abs(v-1) > testTol {
		t.Errorf("phi(0)=%v want 1", v)
	}
}

func TestEntropyFairCoin(t *testing.T) {
	d := mustDist(Bernoulli(0.5))
	if !approx(d.Entropy(), math.Ln2, testTol) {
		t.Errorf("entropy=%g want ln2=%g", d.Entropy(), math.Ln2)
	}
	if !approx(d.EntropyBits(), 1, testTol) {
		t.Errorf("entropy bits=%g want 1", d.EntropyBits())
	}
}

func TestEntropyUniformDie(t *testing.T) {
	d := mustDist(DiscreteUniform(1, 6))
	if !approx(d.EntropyBits(), math.Log2(6), testTol) {
		t.Errorf("die entropy bits=%g want %g", d.EntropyBits(), math.Log2(6))
	}
}

func TestExpectationArbitraryFunction(t *testing.T) {
	d := mustDist(DiscreteUniform(1, 6))
	// E[X^2] via Expectation should equal the second raw moment.
	got := d.Expectation(func(x float64) float64 { return x * x })
	if !approx(got, d.Moment(2), testTol) {
		t.Errorf("E[X^2]=%g want %g", got, d.Moment(2))
	}
}
