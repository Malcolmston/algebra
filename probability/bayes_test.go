package probability

import (
	"math"
	"testing"
)

func TestBayesDiseaseTest(t *testing.T) {
	// P(D)=0.01, P(+|D)=0.99, P(+|~D)=0.05.
	priors := []float64{0.01, 0.99}
	likelihoods := []float64{0.99, 0.05}
	pB, err := TotalProbability(priors, likelihoods)
	if err != nil {
		t.Fatal(err)
	}
	wantPB := 0.99*0.01 + 0.05*0.99
	if !approx(pB, wantPB, testTol) {
		t.Errorf("total prob=%g want %g", pB, wantPB)
	}
	post, err := BayesPosterior(priors, likelihoods)
	if err != nil {
		t.Fatal(err)
	}
	want := 0.99 * 0.01 / wantPB // ≈ 0.16666...
	if !approx(post[0], want, testTol) {
		t.Errorf("posterior P(D|+)=%g want %g", post[0], want)
	}
	if !approx(post[0]+post[1], 1, testTol) {
		t.Error("posterior not normalized")
	}
	// Cross-check against the two-argument Bayes formula.
	single, err := Bayes(0.99, 0.01, wantPB)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(single, want, testTol) {
		t.Errorf("Bayes=%g want %g", single, want)
	}
}

func TestConditionalProbability(t *testing.T) {
	got, err := ConditionalProbability(0.12, 0.3)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(got, 0.4, testTol) {
		t.Errorf("P(A|B)=%g want 0.4", got)
	}
	if _, err := ConditionalProbability(0.5, 0); err == nil {
		t.Error("expected zero-conditioning error")
	}
	if _, err := ConditionalProbability(0.6, 0.5); err == nil {
		t.Error("expected joint>marginal error")
	}
}

func TestElementaryLaws(t *testing.T) {
	if !approx(Complement(0.3), 0.7, testTol) {
		t.Errorf("Complement(0.3)=%g", Complement(0.3))
	}
	if !approx(UnionProbability(0.5, 0.4, 0.2), 0.7, testTol) {
		t.Errorf("union=%g want 0.7", UnionProbability(0.5, 0.4, 0.2))
	}
	if !approx(JointProbabilityIndependent(0.5, 0.4), 0.2, testTol) {
		t.Errorf("joint indep=%g want 0.2", JointProbabilityIndependent(0.5, 0.4))
	}
	if !AreIndependent(0.5, 0.4, 0.2) {
		t.Error("0.5,0.4,0.2 should be independent")
	}
	if AreIndependent(0.5, 0.4, 0.25) {
		t.Error("0.5,0.4,0.25 should not be independent")
	}
}

func TestOddsRoundTrip(t *testing.T) {
	for _, p := range []float64{0.1, 0.25, 0.5, 0.9} {
		o := OddsFromProbability(p)
		back := ProbabilityFromOdds(o)
		if !approx(back, p, testTol) {
			t.Errorf("odds round-trip p=%g got %g", p, back)
		}
	}
	if !math.IsInf(OddsFromProbability(1), 1) {
		t.Error("odds at p=1 should be +Inf")
	}
	if !approx(ProbabilityFromOdds(math.Inf(1)), 1, testTol) {
		t.Error("prob at infinite odds should be 1")
	}
}
