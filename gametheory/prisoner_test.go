package gametheory

import "testing"

func TestPrisonersDilemmaValidation(t *testing.T) {
	if _, err := NewPrisonersDilemma(3, 5, 1, 0); err == nil {
		t.Fatal("expected error: R>T violates ordering")
	}
	pd, err := NewPrisonersDilemma(5, 3, 1, 0)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if !pd.IsValid() || !pd.IsIterationStable() {
		t.Fatal("standard PD should be valid and iteration-stable")
	}
}

func TestIteratedPDTitForTatVsAlwaysDefect(t *testing.T) {
	pd := StandardPrisonersDilemma()
	// TitForTat cooperates once then retaliates; AlwaysDefect defects always.
	// Round 1: (C,D) -> TfT gets S=0, AllD gets T=5.
	// Rounds 2..n: (D,D) -> both get P=1.
	rounds := 10
	sa, sb, _, _ := PlayIteratedPD(pd, TitForTat, AlwaysDefect, rounds)
	wantA := pd.S + float64(rounds-1)*pd.P // 0 + 9 = 9
	wantB := pd.T + float64(rounds-1)*pd.P // 5 + 9 = 14
	if !gametheoryApprox(sa, wantA, 1e-12) || !gametheoryApprox(sb, wantB, 1e-12) {
		t.Fatalf("scores = (%v,%v), want (%v,%v)", sa, sb, wantA, wantB)
	}
}

func TestIteratedPDMutualCooperation(t *testing.T) {
	pd := StandardPrisonersDilemma()
	rounds := 5
	sa, sb, _, _ := PlayIteratedPD(pd, TitForTat, TitForTat, rounds)
	// Two TfT players cooperate every round: each earns R per round.
	want := float64(rounds) * pd.R
	if !gametheoryApprox(sa, want, 1e-12) || !gametheoryApprox(sb, want, 1e-12) {
		t.Fatalf("TfT vs TfT scores = (%v,%v), want (%v,%v)", sa, sb, want, want)
	}
}

func TestGrimTrigger(t *testing.T) {
	// After a single defection GrimTrigger defects forever.
	if GrimTrigger([]int{Cooperate}, []int{Cooperate, Defect}) != Defect {
		t.Fatal("GrimTrigger should defect after a defection")
	}
	if GrimTrigger(nil, nil) != Cooperate {
		t.Fatal("GrimTrigger should open with cooperation")
	}
}

func TestWinStayLoseShift(t *testing.T) {
	// After mutual cooperation (a win), stay on Cooperate.
	if WinStayLoseShift([]int{Cooperate}, []int{Cooperate}) != Cooperate {
		t.Fatal("WSLS should stay after mutual cooperation")
	}
	// After being suckered (C vs D, a loss), shift to Defect.
	if WinStayLoseShift([]int{Cooperate}, []int{Defect}) != Defect {
		t.Fatal("WSLS should shift after a loss")
	}
	// After successfully defecting on a cooperator (a win), stay on Defect.
	if WinStayLoseShift([]int{Defect}, []int{Cooperate}) != Defect {
		t.Fatal("WSLS should stay on Defect after exploiting a cooperator")
	}
}

func TestCanonicalGamesEquilibria(t *testing.T) {
	// Chicken has exactly two pure Nash equilibria (the off-diagonal cells).
	eq := GameOfChicken().PureNashEquilibria(1e-9)
	if len(eq) != 2 {
		t.Fatalf("Chicken pure Nash count = %d, want 2", len(eq))
	}
}

// BenchmarkMixedNashBattleOfSexes exercises the heaviest routine: full support
// enumeration for mixed Nash equilibria (subset pairs, linear solves, and
// best-response verification).
func BenchmarkMixedNashBattleOfSexes(b *testing.B) {
	// A 4x4 game gives support enumeration a non-trivial amount of work.
	g, _ := NewGame(
		[][]float64{
			{4, 1, 2, 0},
			{0, 3, 1, 2},
			{2, 0, 3, 1},
			{1, 2, 0, 4},
		},
		[][]float64{
			{3, 0, 1, 2},
			{1, 4, 0, 1},
			{2, 1, 3, 0},
			{0, 2, 1, 3},
		},
	)
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		_ = g.MixedNashEquilibria(1e-9)
	}
}
