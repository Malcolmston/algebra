package gametheory

import (
	"errors"
	"math"
)

// Game is a finite two-player game in normal (strategic) form. The row player
// chooses a row i in [0, m) and the column player a column j in [0, n); the
// resulting payoffs are Row[i][j] to the row player and Col[i][j] to the
// column player. Both matrices are m-by-n and are shared by value on copy, so
// callers should treat a Game as immutable once constructed.
type Game struct {
	// Row holds the row player's payoffs, indexed Row[i][j].
	Row [][]float64
	// Col holds the column player's payoffs, indexed Col[i][j].
	Col [][]float64
}

// ErrShape is returned when payoff matrices are empty, ragged, or mismatched.
var ErrShape = errors.New("gametheory: payoff matrices must be non-empty, rectangular and of equal shape")

// NewGame constructs a Game from the row and column payoff matrices, which must
// be non-empty, rectangular, and of identical dimensions. The input slices are
// deep-copied so later mutation of the arguments does not affect the Game.
func NewGame(row, col [][]float64) (Game, error) {
	if err := gametheoryValidateShape(row, col); err != nil {
		return Game{}, err
	}
	return Game{Row: gametheoryCopyMatrix(row), Col: gametheoryCopyMatrix(col)}, nil
}

// NewZeroSumGame constructs a zero-sum Game from a single payoff matrix
// interpreted as the row player's payoff; the column player's payoff is its
// negation. The matrix must be non-empty and rectangular.
func NewZeroSumGame(payoff [][]float64) (Game, error) {
	if err := gametheoryValidateMatrix(payoff); err != nil {
		return Game{}, err
	}
	row := gametheoryCopyMatrix(payoff)
	col := make([][]float64, len(payoff))
	for i := range payoff {
		col[i] = make([]float64, len(payoff[i]))
		for j := range payoff[i] {
			col[i][j] = -payoff[i][j]
		}
	}
	return Game{Row: row, Col: col}, nil
}

// Validate reports whether the receiver's payoff matrices are well-formed
// (non-empty, rectangular and of equal shape), returning ErrShape otherwise.
func (g Game) Validate() error {
	return gametheoryValidateShape(g.Row, g.Col)
}

// NumRowStrategies returns m, the number of pure strategies available to the
// row player.
func (g Game) NumRowStrategies() int { return len(g.Row) }

// NumColStrategies returns n, the number of pure strategies available to the
// column player.
func (g Game) NumColStrategies() int {
	if len(g.Row) == 0 {
		return 0
	}
	return len(g.Row[0])
}

// RowPayoff returns the row player's payoff when the row player plays i and the
// column player plays j.
func (g Game) RowPayoff(i, j int) float64 { return g.Row[i][j] }

// ColPayoff returns the column player's payoff when the row player plays i and
// the column player plays j.
func (g Game) ColPayoff(i, j int) float64 { return g.Col[i][j] }

// IsZeroSum reports whether the game is zero-sum, i.e. Row[i][j]+Col[i][j] is
// zero for every cell within the tolerance tol.
func (g Game) IsZeroSum(tol float64) bool {
	for i := range g.Row {
		for j := range g.Row[i] {
			if math.Abs(g.Row[i][j]+g.Col[i][j]) > tol {
				return false
			}
		}
	}
	return true
}

// IsConstantSum reports whether the game is constant-sum, i.e. Row[i][j] +
// Col[i][j] equals the same constant in every cell within tolerance tol. The
// common sum and true are returned when the game is constant-sum; otherwise the
// sum of the first cell and false are returned.
func (g Game) IsConstantSum(tol float64) (float64, bool) {
	if len(g.Row) == 0 || len(g.Row[0]) == 0 {
		return 0, false
	}
	c := g.Row[0][0] + g.Col[0][0]
	for i := range g.Row {
		for j := range g.Row[i] {
			if math.Abs(g.Row[i][j]+g.Col[i][j]-c) > tol {
				return c, false
			}
		}
	}
	return c, true
}

// IsSymmetric reports whether the game is symmetric: it must be square and
// satisfy Col[i][j] == Row[j][i] for all i, j within tolerance tol, so that the
// two players face identical strategic situations.
func (g Game) IsSymmetric(tol float64) bool {
	m, n := g.NumRowStrategies(), g.NumColStrategies()
	if m != n {
		return false
	}
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			if math.Abs(g.Col[i][j]-g.Row[j][i]) > tol {
				return false
			}
		}
	}
	return true
}

// Transpose returns the game viewed from the column player's seat: the players
// swap roles, so the returned game has n rows and m columns with row payoffs
// Col[j][i] and column payoffs Row[j][i].
func (g Game) Transpose() Game {
	m, n := g.NumRowStrategies(), g.NumColStrategies()
	row := make([][]float64, n)
	col := make([][]float64, n)
	for j := 0; j < n; j++ {
		row[j] = make([]float64, m)
		col[j] = make([]float64, m)
		for i := 0; i < m; i++ {
			row[j][i] = g.Col[i][j]
			col[j][i] = g.Row[i][j]
		}
	}
	return Game{Row: row, Col: col}
}

// PureProfile is a pair of pure strategies, one for each player: Row is the row
// player's strategy index and Col the column player's.
type PureProfile struct {
	Row int
	Col int
}

// ExpectedRowPayoff returns the row player's expected payoff when the row player
// plays mixed strategy p and the column player plays mixed strategy q.
func (g Game) ExpectedRowPayoff(p, q MixedStrategy) float64 {
	return ExpectedPayoff(g.Row, p, q)
}

// ExpectedColPayoff returns the column player's expected payoff when the row
// player plays mixed strategy p and the column player plays mixed strategy q.
func (g Game) ExpectedColPayoff(p, q MixedStrategy) float64 {
	return ExpectedPayoff(g.Col, p, q)
}

// ExpectedPayoff returns the bilinear form p^T A q, the expected value of the
// payoff matrix A when the row is drawn from distribution p and the column from
// distribution q. It panics if the lengths of p, q do not match A's shape.
func ExpectedPayoff(a [][]float64, p, q MixedStrategy) float64 {
	var sum float64
	for i := range a {
		if p[i] == 0 {
			continue
		}
		var rowSum float64
		for j := range a[i] {
			rowSum += a[i][j] * q[j]
		}
		sum += p[i] * rowSum
	}
	return sum
}

// gametheoryValidateShape checks that row and col are non-empty, rectangular
// and of equal dimension.
func gametheoryValidateShape(row, col [][]float64) error {
	if err := gametheoryValidateMatrix(row); err != nil {
		return err
	}
	if err := gametheoryValidateMatrix(col); err != nil {
		return err
	}
	if len(row) != len(col) {
		return ErrShape
	}
	for i := range row {
		if len(row[i]) != len(col[i]) {
			return ErrShape
		}
	}
	return nil
}

// gametheoryValidateMatrix checks that a is non-empty and rectangular.
func gametheoryValidateMatrix(a [][]float64) error {
	if len(a) == 0 || len(a[0]) == 0 {
		return ErrShape
	}
	n := len(a[0])
	for i := range a {
		if len(a[i]) != n {
			return ErrShape
		}
	}
	return nil
}

// gametheoryCopyMatrix returns a deep copy of a.
func gametheoryCopyMatrix(a [][]float64) [][]float64 {
	out := make([][]float64, len(a))
	for i := range a {
		out[i] = append([]float64(nil), a[i]...)
	}
	return out
}
