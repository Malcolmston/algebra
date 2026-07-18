package gametheory

// MatchingPennies returns the classic zero-sum game of matching pennies. Each
// player chooses Heads (strategy 0) or Tails (strategy 1); the row player wins
// (+1) when the pennies match and loses (-1) otherwise, and the column player's
// payoffs are the negation. Its unique Nash equilibrium is both players mixing
// 50/50, with game value 0.
func MatchingPennies() Game {
	row := [][]float64{
		{1, -1},
		{-1, 1},
	}
	col := [][]float64{
		{-1, 1},
		{1, -1},
	}
	return Game{Row: row, Col: col}
}

// StagHunt returns the Stag Hunt coordination game with the common payoffs
// (Stag,Stag)=(4,4), (Stag,Hare)=(1,3), (Hare,Stag)=(3,1), (Hare,Hare)=(3,3),
// where strategy 0 is Stag and 1 is Hare. It has two pure Nash equilibria, the
// payoff-dominant (Stag,Stag) and the risk-dominant (Hare,Hare).
func StagHunt() Game {
	row := [][]float64{
		{4, 1},
		{3, 3},
	}
	col := [][]float64{
		{4, 3},
		{1, 3},
	}
	return Game{Row: row, Col: col}
}

// BattleOfTheSexes returns the Battle of the Sexes coordination game. Each
// player prefers a different one of two events (strategy 0 and strategy 1) but
// both prefer being together to being apart: the row player earns 2 at
// (0,0) and 1 at (1,1), the column player 1 at (0,0) and 2 at (1,1), and both
// earn 0 when they miscoordinate. It has two pure Nash equilibria and one mixed.
func BattleOfTheSexes() Game {
	row := [][]float64{
		{2, 0},
		{0, 1},
	}
	col := [][]float64{
		{1, 0},
		{0, 2},
	}
	return Game{Row: row, Col: col}
}

// GameOfChicken returns the game of Chicken (also called Hawk-Dove). Each player
// chooses Dare/Swerve (strategy 0 = Swerve, 1 = Dare): mutual swerving pays
// (3,3), swerving against a darer pays the swerver 1 and the darer 4, and mutual
// daring is the disaster (0,0). It has two pure Nash equilibria, one for each
// player daring while the other swerves, and one symmetric mixed equilibrium.
func GameOfChicken() Game {
	row := [][]float64{
		{3, 1},
		{4, 0},
	}
	col := [][]float64{
		{3, 4},
		{1, 0},
	}
	return Game{Row: row, Col: col}
}
