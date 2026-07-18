package gametheory

import "testing"

func TestDominancePrisonersDilemma(t *testing.T) {
	pd := StandardPrisonersDilemma().Game()
	// Defect strictly dominates Cooperate for both players.
	if !pd.RowStrictlyDominates(Defect, Cooperate) {
		t.Fatal("Defect should strictly dominate Cooperate for the row player")
	}
	if !pd.ColStrictlyDominates(Defect, Cooperate) {
		t.Fatal("Defect should strictly dominate Cooperate for the column player")
	}
	if r, ok := pd.DominantRowStrategy(); !ok || r != Defect {
		t.Fatalf("dominant row = (%d,%v), want (Defect,true)", r, ok)
	}
	rows, cols := pd.IteratedStrictDominance()
	if len(rows) != 1 || rows[0] != Defect || len(cols) != 1 || cols[0] != Defect {
		t.Fatalf("iterated strict dominance survivors = (%v,%v), want ([1],[1])", rows, cols)
	}
}

func TestStrictlyDominatedLists(t *testing.T) {
	pd := StandardPrisonersDilemma().Game()
	if got := pd.StrictlyDominatedRowStrategies(); len(got) != 1 || got[0] != Cooperate {
		t.Fatalf("strictly dominated rows = %v, want [0]", got)
	}
	if got := pd.StrictlyDominatedColStrategies(); len(got) != 1 || got[0] != Cooperate {
		t.Fatalf("strictly dominated cols = %v, want [0]", got)
	}
}

func TestIteratedDominanceCascade(t *testing.T) {
	// A 3x3 game reducible by iterated strict dominance to a single cell.
	// Row payoffs; column payoffs mirror so both eliminate down to index 0.
	row := [][]float64{
		{3, 3, 3},
		{2, 2, 2},
		{1, 1, 1},
	}
	col := [][]float64{
		{3, 2, 1},
		{3, 2, 1},
		{3, 2, 1},
	}
	g, _ := NewGame(row, col)
	rows, cols := g.IteratedStrictDominance()
	if len(rows) != 1 || rows[0] != 0 || len(cols) != 1 || cols[0] != 0 {
		t.Fatalf("survivors = (%v,%v), want ([0],[0])", rows, cols)
	}
}

func TestWeakDominance(t *testing.T) {
	// Row 0 weakly (not strictly) dominates row 1: equal in col 0, better in col 1.
	row := [][]float64{
		{2, 3},
		{2, 1},
	}
	col := [][]float64{
		{0, 0},
		{0, 0},
	}
	g, _ := NewGame(row, col)
	if g.RowStrictlyDominates(0, 1) {
		t.Fatal("row 0 should not strictly dominate row 1")
	}
	if !g.RowWeaklyDominates(0, 1) {
		t.Fatal("row 0 should weakly dominate row 1")
	}
}
