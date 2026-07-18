package gametheory

import (
	"errors"
	"math"
)

// ZeroSumSolution is the solution of a two-player zero-sum game: Value is the
// value of the game to the row (maximizing) player, RowStrategy is an optimal
// maximin mixed strategy for the row player, and ColStrategy is an optimal
// minimax mixed strategy for the column player.
type ZeroSumSolution struct {
	Value       float64
	RowStrategy MixedStrategy
	ColStrategy MixedStrategy
}

// ErrNotZeroSum is returned by Game.SolveZeroSum when the game is not zero-sum.
var ErrNotZeroSum = errors.New("gametheory: game is not zero-sum")

// ErrSolve is returned when the internal linear program fails to solve.
var ErrSolve = errors.New("gametheory: failed to solve linear program")

// SolveZeroSum solves the matrix game whose payoff matrix (to the maximizing
// row player) is payoff, returning the value of the game and a pair of optimal
// minimax mixed strategies. It uses the linear-programming formulation of a
// matrix game, solved by a self-contained simplex method, and is exact for
// rational data up to floating-point round-off. The matrix must be non-empty
// and rectangular.
func SolveZeroSum(payoff [][]float64) (ZeroSumSolution, error) {
	if err := gametheoryValidateMatrix(payoff); err != nil {
		return ZeroSumSolution{}, err
	}
	m, n := len(payoff), len(payoff[0])

	// Shift so every entry is strictly positive: A' = A + shift, shift chosen so
	// min(A') = 1. This keeps the game value positive without changing optimal
	// strategies.
	minv := math.Inf(1)
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			if payoff[i][j] < minv {
				minv = payoff[i][j]
			}
		}
	}
	shift := 0.0
	if minv < 1 {
		shift = 1 - minv
	}
	ap := make([][]float64, m)
	for i := 0; i < m; i++ {
		ap[i] = make([]float64, n)
		for j := 0; j < n; j++ {
			ap[i][j] = payoff[i][j] + shift
		}
	}

	// Column player's LP: maximize sum_j y'_j subject to, for each row i,
	// sum_j A'[i][j] y'_j <= 1, y' >= 0. The optimal dual variables give the row
	// player's normalized strategy.
	c := make([]float64, n)
	for j := range c {
		c[j] = 1
	}
	b := make([]float64, m)
	for i := range b {
		b[i] = 1
	}
	yprime, dual, opt, ok := gametheorySimplexMax(ap, b, c)
	if !ok || opt <= 1e-12 {
		return ZeroSumSolution{}, ErrSolve
	}
	vShifted := 1 / opt // value of the shifted game

	colStrat := make(MixedStrategy, n)
	for j := 0; j < n; j++ {
		colStrat[j] = yprime[j] * vShifted
	}
	rowStrat := make(MixedStrategy, m)
	for i := 0; i < m; i++ {
		rowStrat[i] = dual[i] * vShifted
	}
	gametheoryRenormalize(rowStrat)
	gametheoryRenormalize(colStrat)

	return ZeroSumSolution{
		Value:       vShifted - shift,
		RowStrategy: rowStrat,
		ColStrategy: colStrat,
	}, nil
}

// SolveZeroSum solves the receiver as a zero-sum game, returning ErrNotZeroSum
// if Row[i][j]+Col[i][j] is not zero (within 1e-9) in every cell. The row
// player's payoff matrix is used as the matrix game.
func (g Game) SolveZeroSum() (ZeroSumSolution, error) {
	if !g.IsZeroSum(1e-9) {
		return ZeroSumSolution{}, ErrNotZeroSum
	}
	return SolveZeroSum(g.Row)
}

// GameValue returns the value (to the maximizing row player) of the matrix game
// with the given row-player payoff matrix.
func GameValue(payoff [][]float64) (float64, error) {
	sol, err := SolveZeroSum(payoff)
	if err != nil {
		return 0, err
	}
	return sol.Value, nil
}

// Maximin returns the row player's pure maximin strategy and its security
// value: the row that maximizes the row's worst-case (minimum over columns)
// payoff, and that maximized minimum.
func Maximin(payoff [][]float64) (int, float64) {
	bestRow, bestVal := -1, math.Inf(-1)
	for i := range payoff {
		rowMin := math.Inf(1)
		for j := range payoff[i] {
			if payoff[i][j] < rowMin {
				rowMin = payoff[i][j]
			}
		}
		if rowMin > bestVal {
			bestVal = rowMin
			bestRow = i
		}
	}
	return bestRow, bestVal
}

// Minimax returns the column player's pure minimax strategy and its value: the
// column that minimizes the column's worst-case (maximum over rows) payoff to
// the row player, and that minimized maximum.
func Minimax(payoff [][]float64) (int, float64) {
	bestCol, bestVal := -1, math.Inf(1)
	n := 0
	if len(payoff) > 0 {
		n = len(payoff[0])
	}
	for j := 0; j < n; j++ {
		colMax := math.Inf(-1)
		for i := range payoff {
			if payoff[i][j] > colMax {
				colMax = payoff[i][j]
			}
		}
		if colMax < bestVal {
			bestVal = colMax
			bestCol = j
		}
	}
	return bestCol, bestVal
}

// SaddlePoint reports whether the matrix game has a pure-strategy saddle point,
// i.e. a cell that is simultaneously the minimum of its row and the maximum of
// its column, equivalently where the pure maximin equals the pure minimax. When
// one exists it returns its row and column indices and true; otherwise it
// returns -1, -1, false. The saddle-point value equals payoff[i][j].
func SaddlePoint(payoff [][]float64) (int, int, bool) {
	for i := range payoff {
		rowMin := math.Inf(1)
		for j := range payoff[i] {
			if payoff[i][j] < rowMin {
				rowMin = payoff[i][j]
			}
		}
		for j := range payoff[i] {
			if payoff[i][j] != rowMin {
				continue
			}
			colMax := math.Inf(-1)
			for r := range payoff {
				if payoff[r][j] > colMax {
					colMax = payoff[r][j]
				}
			}
			if payoff[i][j] == colMax {
				return i, j, true
			}
		}
	}
	return -1, -1, false
}

// gametheorySimplexMax maximizes c·x subject to A x <= b and x >= 0, with all
// entries of b non-negative, using Bland's rule to guarantee termination. It
// returns the optimal primal solution x, the optimal dual solution (shadow
// prices, one per constraint), the optimal objective value, and true; it
// returns ok=false if the problem is unbounded or malformed.
func gametheorySimplexMax(a [][]float64, b, c []float64) (x, dual []float64, opt float64, ok bool) {
	m := len(a)
	if m == 0 {
		return nil, nil, 0, false
	}
	n := len(c)
	total := n + m
	t := make([][]float64, m+1)
	for i := range t {
		t[i] = make([]float64, total+1)
	}
	basis := make([]int, m)
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			t[i][j] = a[i][j]
		}
		t[i][n+i] = 1
		t[i][total] = b[i]
		basis[i] = n + i
	}
	for j := 0; j < n; j++ {
		t[m][j] = -c[j]
	}
	const eps = 1e-12
	for iter := 0; iter < 100000; iter++ {
		pivotCol := -1
		for j := 0; j < total; j++ {
			if t[m][j] < -eps {
				pivotCol = j
				break
			}
		}
		if pivotCol == -1 {
			break
		}
		pivotRow := -1
		best := math.Inf(1)
		for i := 0; i < m; i++ {
			if t[i][pivotCol] > eps {
				ratio := t[i][total] / t[i][pivotCol]
				if ratio < best-eps ||
					(math.Abs(ratio-best) <= eps && (pivotRow == -1 || basis[i] < basis[pivotRow])) {
					best = ratio
					pivotRow = i
				}
			}
		}
		if pivotRow == -1 {
			return nil, nil, 0, false // unbounded
		}
		gametheoryPivot(t, pivotRow, pivotCol)
		basis[pivotRow] = pivotCol
	}
	x = make([]float64, n)
	for i := 0; i < m; i++ {
		if basis[i] < n {
			x[basis[i]] = t[i][total]
		}
	}
	dual = make([]float64, m)
	for i := 0; i < m; i++ {
		dual[i] = t[m][n+i]
	}
	return x, dual, t[m][total], true
}

// gametheoryPivot performs a simplex pivot on tableau t about (row, col),
// normalizing the pivot row and eliminating the pivot column from all others.
func gametheoryPivot(t [][]float64, row, col int) {
	pv := t[row][col]
	for j := range t[row] {
		t[row][j] /= pv
	}
	for r := range t {
		if r == row {
			continue
		}
		factor := t[r][col]
		if factor == 0 {
			continue
		}
		for j := range t[r] {
			t[r][j] -= factor * t[row][j]
		}
	}
}
