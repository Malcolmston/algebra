package gametheory

import (
	"math"
	"testing"
)

// gametheoryApprox reports whether a and b agree within tol.
func gametheoryApprox(a, b, tol float64) bool { return math.Abs(a-b) <= tol }

// gametheoryVecApprox reports whether a and b agree componentwise within tol.
func gametheoryVecApprox(a, b []float64, tol float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if math.Abs(a[i]-b[i]) > tol {
			return false
		}
	}
	return true
}

func TestNewGameShape(t *testing.T) {
	if _, err := NewGame([][]float64{{1, 2}}, [][]float64{{1}}); err == nil {
		t.Fatal("expected shape error for mismatched matrices")
	}
	if _, err := NewGame(nil, nil); err == nil {
		t.Fatal("expected shape error for empty matrices")
	}
	g, err := NewGame([][]float64{{1, 2}}, [][]float64{{3, 4}})
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if g.NumRowStrategies() != 1 || g.NumColStrategies() != 2 {
		t.Fatalf("dims = %d x %d, want 1 x 2", g.NumRowStrategies(), g.NumColStrategies())
	}
}

func TestZeroSumAndConstantSum(t *testing.T) {
	mp := MatchingPennies()
	if !mp.IsZeroSum(1e-12) {
		t.Fatal("matching pennies should be zero-sum")
	}
	pd := StandardPrisonersDilemma().Game()
	if pd.IsZeroSum(1e-9) {
		t.Fatal("prisoner's dilemma is not zero-sum")
	}
	// Constant-sum game with constant 10.
	g, _ := NewGame([][]float64{{3, 7}, {5, 2}}, [][]float64{{7, 3}, {5, 8}})
	c, ok := g.IsConstantSum(1e-12)
	if !ok || !gametheoryApprox(c, 10, 1e-12) {
		t.Fatalf("constant-sum = (%v, %v), want (10, true)", c, ok)
	}
}

func TestSymmetryAndTranspose(t *testing.T) {
	if !StandardPrisonersDilemma().Game().IsSymmetric(1e-12) {
		t.Fatal("prisoner's dilemma should be symmetric")
	}
	// Matching pennies is zero-sum but not symmetric (Col != Row^T).
	if MatchingPennies().IsSymmetric(1e-12) {
		t.Fatal("matching pennies should not be symmetric")
	}
	g, _ := NewGame([][]float64{{1, 2, 3}}, [][]float64{{4, 5, 6}})
	tr := g.Transpose()
	if tr.NumRowStrategies() != 3 || tr.NumColStrategies() != 1 {
		t.Fatalf("transpose dims = %d x %d, want 3 x 1", tr.NumRowStrategies(), tr.NumColStrategies())
	}
	// Row player of the transpose earns the original column payoff.
	if !gametheoryApprox(tr.RowPayoff(0, 0), 4, 1e-12) {
		t.Fatalf("transposed row payoff = %v, want 4", tr.RowPayoff(0, 0))
	}
}

func TestExpectedPayoff(t *testing.T) {
	mp := MatchingPennies()
	p := UniformStrategy(2)
	q := UniformStrategy(2)
	if !gametheoryApprox(mp.ExpectedRowPayoff(p, q), 0, 1e-12) {
		t.Fatalf("expected payoff under uniform play = %v, want 0", mp.ExpectedRowPayoff(p, q))
	}
	// Pure play (Heads, Heads) gives the row player +1.
	if !gametheoryApprox(mp.ExpectedRowPayoff(PureStrategy(2, 0), PureStrategy(2, 0)), 1, 1e-12) {
		t.Fatal("pure (H,H) should pay row player +1")
	}
}

func TestMixedStrategyValidation(t *testing.T) {
	if _, err := NewMixedStrategy([]float64{0.5, 0.4}); err == nil {
		t.Fatal("expected error for probabilities not summing to 1")
	}
	m, err := NewMixedStrategy([]float64{0.25, 0.75})
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if !m.IsValid(1e-9) {
		t.Fatal("strategy should be valid")
	}
	// Entropy of a fair coin is 1 bit.
	if !gametheoryApprox(UniformStrategy(2).Entropy(), 1, 1e-12) {
		t.Fatalf("fair-coin entropy = %v, want 1", UniformStrategy(2).Entropy())
	}
	sup := MixedStrategy{0, 0.5, 0, 0.5}.Support(1e-9)
	if len(sup) != 2 || sup[0] != 1 || sup[1] != 3 {
		t.Fatalf("support = %v, want [1 3]", sup)
	}
}

func TestBestResponses(t *testing.T) {
	pd := StandardPrisonersDilemma().Game()
	// Against any column strategy the row player's best response is Defect.
	br := pd.RowBestResponses(UniformStrategy(2), 1e-9)
	if len(br) != 1 || br[0] != Defect {
		t.Fatalf("row best response = %v, want [Defect]", br)
	}
}
