package gametheory

// RowStrictlyDominates reports whether row strategy i strictly dominates row
// strategy k for the row player: Row[i][j] > Row[k][j] for every column j.
func (g Game) RowStrictlyDominates(i, k int) bool {
	if i == k {
		return false
	}
	for j := range g.Row[i] {
		if g.Row[i][j] <= g.Row[k][j] {
			return false
		}
	}
	return true
}

// RowWeaklyDominates reports whether row strategy i weakly dominates row
// strategy k: Row[i][j] >= Row[k][j] for every column j, with strict inequality
// in at least one column.
func (g Game) RowWeaklyDominates(i, k int) bool {
	if i == k {
		return false
	}
	strict := false
	for j := range g.Row[i] {
		if g.Row[i][j] < g.Row[k][j] {
			return false
		}
		if g.Row[i][j] > g.Row[k][j] {
			strict = true
		}
	}
	return strict
}

// ColStrictlyDominates reports whether column strategy j strictly dominates
// column strategy l for the column player: Col[i][j] > Col[i][l] for every row i.
func (g Game) ColStrictlyDominates(j, l int) bool {
	if j == l {
		return false
	}
	for i := range g.Col {
		if g.Col[i][j] <= g.Col[i][l] {
			return false
		}
	}
	return true
}

// ColWeaklyDominates reports whether column strategy j weakly dominates column
// strategy l: Col[i][j] >= Col[i][l] for every row i, with strict inequality in
// at least one row.
func (g Game) ColWeaklyDominates(j, l int) bool {
	if j == l {
		return false
	}
	strict := false
	for i := range g.Col {
		if g.Col[i][j] < g.Col[i][l] {
			return false
		}
		if g.Col[i][j] > g.Col[i][l] {
			strict = true
		}
	}
	return strict
}

// DominantRowStrategy returns a row strategy that strictly dominates every other
// row strategy, together with true, if one exists; otherwise it returns -1 and
// false.
func (g Game) DominantRowStrategy() (int, bool) {
	m := g.NumRowStrategies()
	for i := 0; i < m; i++ {
		dominant := true
		for k := 0; k < m; k++ {
			if k != i && !g.RowStrictlyDominates(i, k) {
				dominant = false
				break
			}
		}
		if dominant {
			return i, true
		}
	}
	return -1, false
}

// DominantColStrategy returns a column strategy that strictly dominates every
// other column strategy, together with true, if one exists; otherwise it
// returns -1 and false.
func (g Game) DominantColStrategy() (int, bool) {
	n := g.NumColStrategies()
	for j := 0; j < n; j++ {
		dominant := true
		for l := 0; l < n; l++ {
			if l != j && !g.ColStrictlyDominates(j, l) {
				dominant = false
				break
			}
		}
		if dominant {
			return j, true
		}
	}
	return -1, false
}

// StrictlyDominatedRowStrategies returns the sorted indices of row strategies
// that are strictly dominated by some other single (pure) row strategy.
func (g Game) StrictlyDominatedRowStrategies() []int {
	m := g.NumRowStrategies()
	var out []int
	for k := 0; k < m; k++ {
		for i := 0; i < m; i++ {
			if i != k && g.RowStrictlyDominates(i, k) {
				out = append(out, k)
				break
			}
		}
	}
	return out
}

// StrictlyDominatedColStrategies returns the sorted indices of column strategies
// that are strictly dominated by some other single (pure) column strategy.
func (g Game) StrictlyDominatedColStrategies() []int {
	n := g.NumColStrategies()
	var out []int
	for l := 0; l < n; l++ {
		for j := 0; j < n; j++ {
			if j != l && g.ColStrictlyDominates(j, l) {
				out = append(out, l)
				break
			}
		}
	}
	return out
}

// IteratedStrictDominance performs iterated elimination of strictly dominated
// pure strategies (by other pure strategies) until no further elimination is
// possible, and returns the sorted indices of the surviving row and column
// strategies. Because only strictly dominated strategies are removed, the
// result is independent of the order of elimination.
func (g Game) IteratedStrictDominance() (rows, cols []int) {
	return g.iteratedDominance(true)
}

// IteratedWeakDominance performs iterated elimination of weakly dominated pure
// strategies (by other pure strategies) until no further elimination is
// possible, and returns the sorted indices of the surviving row and column
// strategies. Note that, unlike strict dominance, the surviving set for weak
// dominance can depend on the order of elimination; this method uses a fixed,
// deterministic lowest-index-first order.
func (g Game) IteratedWeakDominance() (rows, cols []int) {
	return g.iteratedDominance(false)
}

// iteratedDominance is the shared engine for iterated elimination; strict
// selects strict versus weak domination.
func (g Game) iteratedDominance(strict bool) (rows, cols []int) {
	activeRows := gametheoryRange(g.NumRowStrategies())
	activeCols := gametheoryRange(g.NumColStrategies())
	for {
		changed := false
		// Eliminate a dominated row.
		for a := 0; a < len(activeRows); a++ {
			k := activeRows[a]
			dominated := false
			for b := 0; b < len(activeRows); b++ {
				i := activeRows[b]
				if i == k {
					continue
				}
				if gametheoryRowDom(g.Row, i, k, activeCols, strict) {
					dominated = true
					break
				}
			}
			if dominated {
				activeRows = append(activeRows[:a], activeRows[a+1:]...)
				changed = true
				break
			}
		}
		if changed {
			continue
		}
		// Eliminate a dominated column.
		for a := 0; a < len(activeCols); a++ {
			l := activeCols[a]
			dominated := false
			for b := 0; b < len(activeCols); b++ {
				j := activeCols[b]
				if j == l {
					continue
				}
				if gametheoryColDom(g.Col, j, l, activeRows, strict) {
					dominated = true
					break
				}
			}
			if dominated {
				activeCols = append(activeCols[:a], activeCols[a+1:]...)
				changed = true
				break
			}
		}
		if !changed {
			break
		}
	}
	return activeRows, activeCols
}

// gametheoryRowDom reports whether row i dominates row k over the given active
// columns, strictly or weakly.
func gametheoryRowDom(row [][]float64, i, k int, cols []int, strict bool) bool {
	someStrict := false
	for _, j := range cols {
		if strict {
			if row[i][j] <= row[k][j] {
				return false
			}
		} else {
			if row[i][j] < row[k][j] {
				return false
			}
			if row[i][j] > row[k][j] {
				someStrict = true
			}
		}
	}
	if strict {
		return true
	}
	return someStrict
}

// gametheoryColDom reports whether column j dominates column l over the given
// active rows, strictly or weakly.
func gametheoryColDom(col [][]float64, j, l int, rows []int, strict bool) bool {
	someStrict := false
	for _, i := range rows {
		if strict {
			if col[i][j] <= col[i][l] {
				return false
			}
		} else {
			if col[i][j] < col[i][l] {
				return false
			}
			if col[i][j] > col[i][l] {
				someStrict = true
			}
		}
	}
	if strict {
		return true
	}
	return someStrict
}

// gametheoryRange returns the slice [0, 1, ..., n-1].
func gametheoryRange(n int) []int {
	s := make([]int, n)
	for i := range s {
		s[i] = i
	}
	return s
}
