package gametheory

import (
	"testing"
)

func TestPureNashPrisonersDilemma(t *testing.T) {
	pd := StandardPrisonersDilemma().Game()
	eq := pd.PureNashEquilibria(1e-9)
	if len(eq) != 1 || eq[0] != (PureProfile{Row: Defect, Col: Defect}) {
		t.Fatalf("PD pure Nash = %v, want [{1 1}]", eq)
	}
	if !pd.HasPureNash(1e-9) {
		t.Fatal("PD should have a pure Nash")
	}
}

func TestPureNashCoordination(t *testing.T) {
	// Stag Hunt has two pure equilibria: (Stag,Stag) and (Hare,Hare).
	eq := StagHunt().PureNashEquilibria(1e-9)
	want := map[PureProfile]bool{{0, 0}: true, {1, 1}: true}
	if len(eq) != 2 {
		t.Fatalf("StagHunt pure Nash = %v, want 2 equilibria", eq)
	}
	for _, e := range eq {
		if !want[e] {
			t.Fatalf("unexpected equilibrium %v", e)
		}
	}
}

func TestNoPureNashMatchingPennies(t *testing.T) {
	if MatchingPennies().HasPureNash(1e-9) {
		t.Fatal("matching pennies has no pure Nash")
	}
}

func TestMixedNashMatchingPennies(t *testing.T) {
	eqs := MatchingPennies().MixedNashEquilibria(1e-9)
	// Unique equilibrium: both players mix 50/50, value 0 to each.
	found := false
	for _, e := range eqs {
		if gametheoryVecApprox(e.Row, []float64{0.5, 0.5}, 1e-6) &&
			gametheoryVecApprox(e.Col, []float64{0.5, 0.5}, 1e-6) {
			found = true
			if !gametheoryApprox(e.RowValue, 0, 1e-6) || !gametheoryApprox(e.ColValue, 0, 1e-6) {
				t.Fatalf("matching-pennies values = (%v,%v), want (0,0)", e.RowValue, e.ColValue)
			}
		}
	}
	if !found {
		t.Fatalf("expected the 50/50 mixed equilibrium, got %v", eqs)
	}
}

func TestMixedNashBattleOfSexes(t *testing.T) {
	// Battle of the Sexes: payoffs (2,1) and (1,2). The mixed equilibrium has
	// the row player play strategy 0 with probability 2/3 and the column player
	// play strategy 0 with probability 1/3 (closed form for this payoff pattern).
	eqs := BattleOfTheSexes().MixedNashEquilibria(1e-9)
	var mixed *MixedEquilibrium
	for i := range eqs {
		if len(eqs[i].Row.Support(1e-9)) == 2 {
			mixed = &eqs[i]
			break
		}
	}
	if mixed == nil {
		t.Fatalf("expected a fully-mixed equilibrium, got %v", eqs)
	}
	if !gametheoryVecApprox(mixed.Row, []float64{2.0 / 3, 1.0 / 3}, 1e-6) {
		t.Fatalf("row mix = %v, want [2/3 1/3]", mixed.Row)
	}
	if !gametheoryVecApprox(mixed.Col, []float64{1.0 / 3, 2.0 / 3}, 1e-6) {
		t.Fatalf("col mix = %v, want [1/3 2/3]", mixed.Col)
	}
	// Two pure and one mixed equilibrium expected in total.
	if len(eqs) != 3 {
		t.Fatalf("BoS equilibria count = %d, want 3", len(eqs))
	}
}

func TestSolveLinearSanity(t *testing.T) {
	// 2x2 system: x + y = 3, x - y = 1  =>  x = 2, y = 1.
	sol, ok := gametheorySolveLinear([][]float64{{1, 1}, {1, -1}}, []float64{3, 1})
	if !ok || !gametheoryVecApprox(sol, []float64{2, 1}, 1e-12) {
		t.Fatalf("linear solve = %v (ok=%v), want [2 1]", sol, ok)
	}
	// Singular system reports failure.
	if _, ok := gametheorySolveLinear([][]float64{{1, 1}, {2, 2}}, []float64{1, 2}); ok {
		t.Fatal("expected singular system to fail")
	}
}
