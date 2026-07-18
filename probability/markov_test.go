package probability

import "testing"

func twoState(t *testing.T) MarkovChain {
	t.Helper()
	m, err := NewMarkovChain([][]float64{
		{0.9, 0.1},
		{0.5, 0.5},
	})
	if err != nil {
		t.Fatal(err)
	}
	return m
}

// gamblersRuin builds the symmetric random walk on {0,1,2,3} with absorbing
// barriers at states 0 and 3.
func gamblersRuin(t *testing.T) MarkovChain {
	t.Helper()
	m, err := NewMarkovChain([][]float64{
		{1, 0, 0, 0},
		{0.5, 0, 0.5, 0},
		{0, 0.5, 0, 0.5},
		{0, 0, 0, 1},
	})
	if err != nil {
		t.Fatal(err)
	}
	return m
}

func TestMarkovNStep(t *testing.T) {
	m := twoState(t)
	p2, err := m.NStep(2)
	if err != nil {
		t.Fatal(err)
	}
	// Hand-computed P^2.
	want := [][]float64{{0.86, 0.14}, {0.70, 0.30}}
	for i := range want {
		for j := range want[i] {
			if !approx(p2[i][j], want[i][j], testTol) {
				t.Errorf("P^2[%d][%d]=%g want %g", i, j, p2[i][j], want[i][j])
			}
		}
	}
	// NStep(0) is the identity.
	p0, _ := m.NStep(0)
	if !approx(p0[0][0], 1, testTol) || !approx(p0[0][1], 0, testTol) {
		t.Errorf("P^0 not identity: %v", p0)
	}
}

func TestMarkovDistributionAfter(t *testing.T) {
	m := twoState(t)
	got, err := m.DistributionAfter([]float64{1, 0}, 2)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(got[0], 0.86, testTol) || !approx(got[1], 0.14, testTol) {
		t.Errorf("dist after 2 = %v want [0.86,0.14]", got)
	}
}

func TestMarkovStationary(t *testing.T) {
	m := twoState(t)
	pi, err := m.StationaryDistribution()
	if err != nil {
		t.Fatal(err)
	}
	if !approx(pi[0], 5.0/6.0, 1e-9) || !approx(pi[1], 1.0/6.0, 1e-9) {
		t.Errorf("stationary=%v want [5/6,1/6]", pi)
	}
	// Stationarity: pi·P == pi.
	next, _ := m.Step(pi)
	for i := range pi {
		if !approx(next[i], pi[i], 1e-9) {
			t.Errorf("pi not stationary at %d: %g vs %g", i, next[i], pi[i])
		}
	}
	// Mean recurrence times are the reciprocals of the stationary probabilities.
	mr, err := m.MeanRecurrenceTimes()
	if err != nil {
		t.Fatal(err)
	}
	if !approx(mr[0], 1.2, 1e-9) || !approx(mr[1], 6, 1e-9) {
		t.Errorf("mean recurrence=%v want [1.2,6]", mr)
	}
}

func TestMarkovStructuralProperties(t *testing.T) {
	m := twoState(t)
	if !m.IsIrreducible() {
		t.Error("two-state chain should be irreducible")
	}
	if !m.IsRegular() {
		t.Error("two-state chain should be regular")
	}
	if m.IsAbsorbing() {
		t.Error("two-state chain should not be absorbing")
	}
	g := gamblersRuin(t)
	if g.IsIrreducible() {
		t.Error("gambler's ruin is not irreducible")
	}
	if !g.IsAbsorbing() {
		t.Error("gambler's ruin should be absorbing")
	}
}

func TestMarkovAbsorbing(t *testing.T) {
	g := gamblersRuin(t)
	if got := g.AbsorbingStates(); len(got) != 2 || got[0] != 0 || got[1] != 3 {
		t.Errorf("absorbing states=%v want [0 3]", got)
	}
	if got := g.TransientStates(); len(got) != 2 || got[0] != 1 || got[1] != 2 {
		t.Errorf("transient states=%v want [1 2]", got)
	}
	steps, err := g.ExpectedStepsToAbsorption()
	if err != nil {
		t.Fatal(err)
	}
	// Absorbing states take 0 steps; from states 1 and 2 the expected number is 2.
	if !approx(steps[0], 0, testTol) || !approx(steps[3], 0, testTol) {
		t.Errorf("absorbing steps nonzero: %v", steps)
	}
	if !approx(steps[1], 2, 1e-9) || !approx(steps[2], 2, 1e-9) {
		t.Errorf("expected steps=%v want 2 at states 1,2", steps)
	}
	probs, err := g.AbsorptionProbabilities()
	if err != nil {
		t.Fatal(err)
	}
	// From state 1: absorbed at 0 with prob 2/3, at 3 with prob 1/3.
	if !approx(probs[0][0], 2.0/3.0, 1e-9) || !approx(probs[0][1], 1.0/3.0, 1e-9) {
		t.Errorf("absorption from state1=%v want [2/3,1/3]", probs[0])
	}
	if !approx(probs[1][0], 1.0/3.0, 1e-9) || !approx(probs[1][1], 2.0/3.0, 1e-9) {
		t.Errorf("absorption from state2=%v want [1/3,2/3]", probs[1])
	}
}

func TestMarkovReachable(t *testing.T) {
	g := gamblersRuin(t)
	if !g.Reachable(1, 0) {
		t.Error("state 0 should be reachable from state 1")
	}
	if g.Reachable(0, 1) {
		t.Error("state 1 should not be reachable from absorbing state 0")
	}
}

func TestNewMarkovChainErrors(t *testing.T) {
	if _, err := NewMarkovChain([][]float64{{0.5, 0.4}, {0.5, 0.5}}); err == nil {
		t.Error("expected non-stochastic-row error")
	}
	if _, err := NewMarkovChain([][]float64{{1, 0}}); err == nil {
		t.Error("expected non-square error")
	}
}

func BenchmarkConvolvePower(b *testing.B) {
	die, _ := DiscreteUniform(1, 6)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, _ = die.ConvolvePower(50)
	}
}
