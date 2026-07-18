package probability

import (
	"math"
	"testing"
)

func testJoint(t *testing.T) JointDistribution {
	t.Helper()
	j, err := NewJointDistribution(
		[]float64{0, 1},
		[]float64{0, 1},
		[][]float64{{0.1, 0.2}, {0.3, 0.4}},
	)
	if err != nil {
		t.Fatal(err)
	}
	return j
}

func TestJointMarginals(t *testing.T) {
	j := testJoint(t)
	mx := j.MarginalX()
	my := j.MarginalY()
	if !approx(mx.Probs[0], 0.3, testTol) || !approx(mx.Probs[1], 0.7, testTol) {
		t.Errorf("marginal X=%v want [0.3,0.7]", mx.Probs)
	}
	if !approx(my.Probs[0], 0.4, testTol) || !approx(my.Probs[1], 0.6, testTol) {
		t.Errorf("marginal Y=%v want [0.4,0.6]", my.Probs)
	}
}

func TestJointCovarianceCorrelation(t *testing.T) {
	j := testJoint(t)
	// E[XY]=0.4 (only the (1,1) cell contributes), E[X]=0.7, E[Y]=0.6.
	wantCov := 0.4 - 0.7*0.6
	if !approx(j.Covariance(), wantCov, testTol) {
		t.Errorf("cov=%g want %g", j.Covariance(), wantCov)
	}
	wantCorr := wantCov / math.Sqrt(0.21*0.24)
	if !approx(j.Correlation(), wantCorr, testTol) {
		t.Errorf("corr=%g want %g", j.Correlation(), wantCorr)
	}
	if j.Independent() {
		t.Error("joint should not be independent")
	}
}

func TestJointConditional(t *testing.T) {
	j := testJoint(t)
	// P(X|Y=0): cells [0.1,0.3] over P(Y=0)=0.4 → [0.25,0.75].
	cx, err := j.ConditionalXGivenY(0)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(cx.Probs[0], 0.25, testTol) || !approx(cx.Probs[1], 0.75, testTol) {
		t.Errorf("P(X|Y=0)=%v want [0.25,0.75]", cx.Probs)
	}
	// P(Y|X=1): cells [0.3,0.4] over P(X=1)=0.7.
	cy, err := j.ConditionalYGivenX(1)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(cy.Probs[0], 0.3/0.7, testTol) || !approx(cy.Probs[1], 0.4/0.7, testTol) {
		t.Errorf("P(Y|X=1)=%v", cy.Probs)
	}
}

func TestIndependentJoint(t *testing.T) {
	a := mustDist(Bernoulli(0.7))
	b := mustDist(Bernoulli(0.6))
	j := IndependentJoint(a, b)
	if !j.Independent() {
		t.Error("constructed joint should be independent")
	}
	if !approx(j.Covariance(), 0, testTol) {
		t.Errorf("independent covariance=%g want 0", j.Covariance())
	}
	// Marginals should match the inputs.
	if !approx(j.MarginalX().Mean(), 0.7, testTol) {
		t.Errorf("marginal X mean=%g want 0.7", j.MarginalX().Mean())
	}
}

func TestNewJointDistributionErrors(t *testing.T) {
	if _, err := NewJointDistribution([]float64{0}, []float64{0}, [][]float64{{0.5}}); err == nil {
		t.Error("expected sum-not-one error")
	}
	if _, err := NewJointDistribution([]float64{0, 1}, []float64{0}, [][]float64{{1.0}}); err == nil {
		t.Error("expected row-count error")
	}
}
