package gametheory

import "testing"

func TestShapleyTwoPlayerClosedForm(t *testing.T) {
	// v(1)=a, v(2)=b, v(12)=c, v(empty)=0.
	// Closed form: phi1 = (a + c - b)/2, phi2 = (b + c - a)/2.
	a, b, c := 1.0, 2.0, 4.0
	cg, err := NewCooperativeGame(2, func(mask uint) float64 {
		switch mask {
		case 0b01:
			return a
		case 0b10:
			return b
		case 0b11:
			return c
		default:
			return 0
		}
	})
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	phi := cg.ShapleyValue()
	want := []float64{(a + c - b) / 2, (b + c - a) / 2}
	if !gametheoryVecApprox(phi, want, 1e-12) {
		t.Fatalf("Shapley = %v, want %v", phi, want)
	}
	// Efficiency: values sum to v(N).
	if !gametheoryApprox(phi[0]+phi[1], c, 1e-12) {
		t.Fatalf("Shapley sum = %v, want %v", phi[0]+phi[1], c)
	}
}

func TestShapleyWeightedVoting(t *testing.T) {
	// Weighted majority game [q=3; w=2,1,1]. Known Shapley = (2/3, 1/6, 1/6).
	cg, err := WeightedVotingGame([]float64{2, 1, 1}, 3)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	phi := cg.ShapleyValue()
	want := []float64{2.0 / 3, 1.0 / 6, 1.0 / 6}
	if !gametheoryVecApprox(phi, want, 1e-12) {
		t.Fatalf("Shapley = %v, want %v", phi, want)
	}
}

func TestBanzhafWeightedVoting(t *testing.T) {
	// Same game: Banzhaf value = (3,1,1)/2^(n-1) = (0.75, 0.25, 0.25).
	cg, _ := WeightedVotingGame([]float64{2, 1, 1}, 3)
	beta := cg.BanzhafValue()
	want := []float64{0.75, 0.25, 0.25}
	if !gametheoryVecApprox(beta, want, 1e-12) {
		t.Fatalf("Banzhaf = %v, want %v", beta, want)
	}
}

func TestCooperativePredicates(t *testing.T) {
	// A convex (hence superadditive and monotone) game: v(S) = |S|^2.
	cg, _ := NewCooperativeGame(3, func(mask uint) float64 {
		n := 0
		for i := 0; i < 3; i++ {
			if mask&(1<<uint(i)) != 0 {
				n++
			}
		}
		return float64(n * n)
	})
	if !cg.IsConvex(1e-12) {
		t.Fatal("|S|^2 game should be convex")
	}
	if !cg.IsSuperadditive(1e-12) {
		t.Fatal("|S|^2 game should be superadditive")
	}
	if !cg.IsMonotone(1e-12) {
		t.Fatal("|S|^2 game should be monotone")
	}
	if !cg.IsEssential(1e-12) {
		t.Fatal("|S|^2 game should be essential")
	}
	if !gametheoryApprox(cg.GrandCoalitionValue(), 9, 1e-12) {
		t.Fatalf("v(N) = %v, want 9", cg.GrandCoalitionValue())
	}
	// The Shapley value of a convex game lies in the core.
	if !cg.IsInCore(cg.ShapleyValue(), 1e-9) {
		t.Fatal("Shapley value of a convex game must be in the core")
	}
}

func TestMarginalContribution(t *testing.T) {
	cg, _ := WeightedVotingGame([]float64{2, 1, 1}, 3)
	// Player 0 joining {1} (mask 0b010) turns a losing coalition into a winner.
	if !gametheoryApprox(cg.MarginalContribution(0, 0b010), 1, 1e-12) {
		t.Fatalf("marginal contribution = %v, want 1", cg.MarginalContribution(0, 0b010))
	}
}

func TestCoreRejectsInfeasible(t *testing.T) {
	cg, _ := WeightedVotingGame([]float64{2, 1, 1}, 3)
	// Allocation that under-pays a winning coalition is not in the core.
	if cg.IsInCore([]float64{0, 0.5, 0.5}, 1e-9) {
		t.Fatal("allocation should be rejected: coalition {1,2} is losing but {0,1} winning is starved")
	}
}
