package gametheory

import (
	"math"
	"math/bits"
	"sort"
)

// IsPureNash reports whether the pure-strategy profile (i, j) is a Nash
// equilibrium: neither player can strictly improve their payoff by a unilateral
// deviation to another pure strategy, judged within tolerance tol.
func (g Game) IsPureNash(i, j int, tol float64) bool {
	rp := g.Row[i][j]
	for k := range g.Row {
		if g.Row[k][j] > rp+tol {
			return false
		}
	}
	cp := g.Col[i][j]
	for l := range g.Col[i] {
		if g.Col[i][l] > cp+tol {
			return false
		}
	}
	return true
}

// PureNashEquilibria returns all pure-strategy Nash equilibria of the game as
// PureProfile values, in lexicographic (row, column) order. Strategies within
// tolerance tol are treated as tied when testing for profitable deviations.
func (g Game) PureNashEquilibria(tol float64) []PureProfile {
	var eq []PureProfile
	for i := range g.Row {
		for j := range g.Row[i] {
			if g.IsPureNash(i, j, tol) {
				eq = append(eq, PureProfile{Row: i, Col: j})
			}
		}
	}
	return eq
}

// HasPureNash reports whether the game has at least one pure-strategy Nash
// equilibrium, testing deviations within tolerance tol.
func (g Game) HasPureNash(tol float64) bool {
	for i := range g.Row {
		for j := range g.Row[i] {
			if g.IsPureNash(i, j, tol) {
				return true
			}
		}
	}
	return false
}

// MixedEquilibrium is a mixed-strategy Nash equilibrium: Row and Col are the two
// players' equilibrium mixed strategies, and RowValue and ColValue are the
// corresponding equilibrium expected payoffs.
type MixedEquilibrium struct {
	Row      MixedStrategy
	Col      MixedStrategy
	RowValue float64
	ColValue float64
}

// MixedNashEquilibriaWithSupport attempts to find a mixed Nash equilibrium whose
// row-player support is rowSupport and whose column-player support is
// colSupport. The two supports must have equal size. It solves the indifference
// conditions for each player, verifies non-negativity and the best-response
// (no-profitable-deviation) conditions within tolerance tol, and returns the
// equilibrium together with true on success, or the zero value and false when
// no equilibrium exists on the given supports.
func (g Game) MixedNashEquilibriaWithSupport(rowSupport, colSupport []int, tol float64) (MixedEquilibrium, bool) {
	k := len(rowSupport)
	if k == 0 || k != len(colSupport) {
		return MixedEquilibrium{}, false
	}
	m, n := g.NumRowStrategies(), g.NumColStrategies()

	// Column strategy y over colSupport makes the row player indifferent across
	// rowSupport at value w. Unknowns: y_0..y_{k-1}, w.
	y, w, ok := gametheoryIndifference(g.Row, rowSupport, colSupport, tol)
	if !ok {
		return MixedEquilibrium{}, false
	}
	// Row strategy x over rowSupport makes the column player indifferent across
	// colSupport at value v. Column player's payoff Col is indexed [i][j]; from
	// the column player's perspective the "rows" are its strategies, so we work
	// with the transpose so that support/opponent roles line up.
	x, v, ok := gametheoryIndifferenceCol(g.Col, rowSupport, colSupport, tol)
	if !ok {
		return MixedEquilibrium{}, false
	}

	// Assemble full-length strategies, clamping tiny negatives to zero.
	rowStrat := make(MixedStrategy, m)
	for idx, i := range rowSupport {
		if x[idx] < -tol {
			return MixedEquilibrium{}, false
		}
		rowStrat[i] = math.Max(0, x[idx])
	}
	colStrat := make(MixedStrategy, n)
	for idx, j := range colSupport {
		if y[idx] < -tol {
			return MixedEquilibrium{}, false
		}
		colStrat[j] = math.Max(0, y[idx])
	}
	gametheoryRenormalize(rowStrat)
	gametheoryRenormalize(colStrat)

	// Best-response check: no pure strategy off the support does strictly better.
	for i := 0; i < m; i++ {
		var payoff float64
		for j := 0; j < n; j++ {
			payoff += g.Row[i][j] * colStrat[j]
		}
		if payoff > w+1e-7 {
			return MixedEquilibrium{}, false
		}
	}
	for j := 0; j < n; j++ {
		var payoff float64
		for i := 0; i < m; i++ {
			payoff += g.Col[i][j] * rowStrat[i]
		}
		if payoff > v+1e-7 {
			return MixedEquilibrium{}, false
		}
	}
	return MixedEquilibrium{Row: rowStrat, Col: colStrat, RowValue: w, ColValue: v}, true
}

// MixedNashEquilibria returns mixed-strategy Nash equilibria of the two-player
// game found by support enumeration over all pairs of equal-sized supports.
// Pure equilibria appear as singleton-support equilibria, so the result is a
// superset of PureNashEquilibria expressed as mixed strategies. Duplicate
// equilibria (which can arise in degenerate games) are removed within tolerance
// tol. For non-degenerate games this returns the complete finite set of
// equilibria; degenerate games may have a continuum, of which only the
// enumerated vertices are reported.
func (g Game) MixedNashEquilibria(tol float64) []MixedEquilibrium {
	m, n := g.NumRowStrategies(), g.NumColStrategies()
	var out []MixedEquilibrium
	for rowMask := 1; rowMask < (1 << m); rowMask++ {
		rk := bits.OnesCount(uint(rowMask))
		rowSupport := gametheoryMaskToSlice(rowMask, m)
		for colMask := 1; colMask < (1 << n); colMask++ {
			if bits.OnesCount(uint(colMask)) != rk {
				continue
			}
			colSupport := gametheoryMaskToSlice(colMask, n)
			eq, ok := g.MixedNashEquilibriaWithSupport(rowSupport, colSupport, tol)
			if !ok {
				continue
			}
			if !gametheoryContainsEquilibrium(out, eq, 1e-6) {
				out = append(out, eq)
			}
		}
	}
	sort.Slice(out, func(a, b int) bool {
		la, lb := len(out[a].Row.Support(1e-9)), len(out[b].Row.Support(1e-9))
		if la != lb {
			return la < lb
		}
		return out[a].RowValue < out[b].RowValue
	})
	return out
}

// gametheoryIndifference solves for the opponent (column) distribution y over
// colSupport that makes the row player indifferent across rowSupport, returning
// y (aligned with colSupport), the common row payoff w, and ok.
func gametheoryIndifference(payoff [][]float64, rowSupport, colSupport []int, tol float64) ([]float64, float64, bool) {
	k := len(rowSupport)
	// Unknowns: y_0..y_{k-1}, w. Equations: for each r in rowSupport,
	// sum_c payoff[rowSupport[r]][colSupport[c]] * y_c - w = 0; plus sum y_c = 1.
	a := make([][]float64, k+1)
	b := make([]float64, k+1)
	for r := 0; r < k; r++ {
		a[r] = make([]float64, k+1)
		for c := 0; c < k; c++ {
			a[r][c] = payoff[rowSupport[r]][colSupport[c]]
		}
		a[r][k] = -1
		b[r] = 0
	}
	a[k] = make([]float64, k+1)
	for c := 0; c < k; c++ {
		a[k][c] = 1
	}
	b[k] = 1
	sol, ok := gametheorySolveLinear(a, b)
	if !ok {
		return nil, 0, false
	}
	return sol[:k], sol[k], true
}

// gametheoryIndifferenceCol solves for the row distribution x over rowSupport
// that makes the column player indifferent across colSupport, given the column
// payoff matrix (indexed [row][col]).
func gametheoryIndifferenceCol(colPayoff [][]float64, rowSupport, colSupport []int, tol float64) ([]float64, float64, bool) {
	k := len(rowSupport)
	// Unknowns: x_0..x_{k-1}, v. Equations: for each c in colSupport,
	// sum_r colPayoff[rowSupport[r]][colSupport[c]] * x_r - v = 0; plus sum x = 1.
	a := make([][]float64, k+1)
	b := make([]float64, k+1)
	for c := 0; c < k; c++ {
		a[c] = make([]float64, k+1)
		for r := 0; r < k; r++ {
			a[c][r] = colPayoff[rowSupport[r]][colSupport[c]]
		}
		a[c][k] = -1
		b[c] = 0
	}
	a[k] = make([]float64, k+1)
	for r := 0; r < k; r++ {
		a[k][r] = 1
	}
	b[k] = 1
	sol, ok := gametheorySolveLinear(a, b)
	if !ok {
		return nil, 0, false
	}
	return sol[:k], sol[k], true
}

// gametheorySolveLinear solves the square system a x = b by Gaussian
// elimination with partial pivoting, returning the solution and true, or nil
// and false when the system is singular.
func gametheorySolveLinear(a [][]float64, b []float64) ([]float64, bool) {
	n := len(a)
	m := make([][]float64, n)
	for i := range a {
		m[i] = append([]float64(nil), a[i]...)
		m[i] = append(m[i], b[i])
	}
	for col := 0; col < n; col++ {
		pivot := col
		best := math.Abs(m[col][col])
		for r := col + 1; r < n; r++ {
			if v := math.Abs(m[r][col]); v > best {
				best = v
				pivot = r
			}
		}
		if best < 1e-12 {
			return nil, false
		}
		m[col], m[pivot] = m[pivot], m[col]
		pv := m[col][col]
		for r := 0; r < n; r++ {
			if r == col {
				continue
			}
			factor := m[r][col] / pv
			if factor == 0 {
				continue
			}
			for c := col; c <= n; c++ {
				m[r][c] -= factor * m[col][c]
			}
		}
	}
	x := make([]float64, n)
	for i := 0; i < n; i++ {
		x[i] = m[i][n] / m[i][i]
	}
	return x, true
}

// gametheoryRenormalize scales m so its entries sum to one, if the sum is
// positive.
func gametheoryRenormalize(m MixedStrategy) {
	var sum float64
	for _, p := range m {
		sum += p
	}
	if sum > 0 {
		for i := range m {
			m[i] /= sum
		}
	}
}

// gametheoryMaskToSlice returns the sorted indices in [0, n) whose bit is set in
// mask.
func gametheoryMaskToSlice(mask, n int) []int {
	var s []int
	for i := 0; i < n; i++ {
		if mask&(1<<i) != 0 {
			s = append(s, i)
		}
	}
	return s
}

// gametheoryContainsEquilibrium reports whether list already holds an
// equilibrium whose row and column strategies match eq within tol.
func gametheoryContainsEquilibrium(list []MixedEquilibrium, eq MixedEquilibrium, tol float64) bool {
	for _, e := range list {
		if gametheoryVecClose(e.Row, eq.Row, tol) && gametheoryVecClose(e.Col, eq.Col, tol) {
			return true
		}
	}
	return false
}

// gametheoryVecClose reports whether a and b agree componentwise within tol.
func gametheoryVecClose(a, b MixedStrategy, tol float64) bool {
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
