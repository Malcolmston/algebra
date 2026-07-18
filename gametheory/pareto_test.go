package gametheory

import "testing"

func TestParetoDominates(t *testing.T) {
	if !ParetoDominates([]float64{3, 3}, []float64{3, 1}) {
		t.Fatal("(3,3) should Pareto-dominate (3,1)")
	}
	if ParetoDominates([]float64{3, 1}, []float64{1, 3}) {
		t.Fatal("(3,1) should not Pareto-dominate (1,3)")
	}
	if ParetoDominates([]float64{2, 2}, []float64{2, 2}) {
		t.Fatal("equal vectors do not Pareto-dominate")
	}
}

func TestParetoOptimalIndices(t *testing.T) {
	pts := [][]float64{
		{1, 1}, // dominated by {2,2}
		{2, 2},
		{0, 3},
		{3, 0},
	}
	got := ParetoOptimalIndices(pts)
	want := map[int]bool{1: true, 2: true, 3: true}
	if len(got) != 3 {
		t.Fatalf("optimal indices = %v, want 3 of them", got)
	}
	for _, i := range got {
		if !want[i] {
			t.Fatalf("unexpected optimal index %d", i)
		}
	}
}

func TestParetoFrontierPrisonersDilemma(t *testing.T) {
	pd := StandardPrisonersDilemma().Game()
	// (Defect,Defect)=(1,1) is dominated by (Cooperate,Cooperate)=(3,3), so it is
	// not on the Pareto frontier even though it is the unique Nash equilibrium.
	front := pd.ParetoFrontier()
	for _, p := range front {
		if p == (PureProfile{Defect, Defect}) {
			t.Fatal("(Defect,Defect) must not be Pareto-optimal")
		}
	}
	if !pd.IsParetoOptimal(Cooperate, Cooperate) {
		t.Fatal("(Cooperate,Cooperate) should be Pareto-optimal")
	}
	if pd.IsParetoOptimal(Defect, Defect) {
		t.Fatal("(Defect,Defect) should not be Pareto-optimal")
	}
}

func TestSocialOptimum(t *testing.T) {
	pd := StandardPrisonersDilemma().Game()
	prof, sum := pd.SocialOptimum()
	if prof != (PureProfile{Cooperate, Cooperate}) || !gametheoryApprox(sum, 6, 1e-12) {
		t.Fatalf("social optimum = (%v, %v), want ({0 0}, 6)", prof, sum)
	}
}
