package tropical

// StarSeries returns the bounded tropical closure I (+) A (+) A^2 (+) ... (+)
// A^k, the tropical sum of the first k+1 powers of a square matrix. When the
// full Kleene star converges this equals it for every k >= n-1. It returns
// ErrNotSquare for a non-square matrix.
func (m Matrix) StarSeries(k int) (Matrix, error) {
	if !m.IsSquare() {
		return Matrix{}, ErrNotSquare
	}
	acc := Identity(m.sr, m.rows)
	term := Identity(m.sr, m.rows)
	for i := 1; i <= k; i++ {
		term, _ = term.Mul(m)
		acc, _ = acc.Add(term)
	}
	return acc, nil
}

// PlusClosure returns the tropical Kleene plus A^+ = A (+) A^2 (+) A^3 (+) ...,
// the tropical sum of all positive powers of a square matrix, computed with a
// Floyd-Warshall style relaxation. Over min-plus the (i,j) entry is the weight
// of the shortest walk of length at least one from i to j. It returns
// ErrDivergent if a bad cycle makes the closure diverge and ErrNotSquare for a
// non-square matrix.
func (m Matrix) PlusClosure() (Matrix, error) {
	if !m.IsSquare() {
		return Matrix{}, ErrNotSquare
	}
	n := m.rows
	d := m.Clone()
	for k := 0; k < n; k++ {
		for i := 0; i < n; i++ {
			dik := d.data[i][k]
			if m.sr.IsZero(dik) {
				continue
			}
			for j := 0; j < n; j++ {
				d.data[i][j] = m.sr.Add(d.data[i][j], m.sr.Mul(dik, d.data[k][j]))
			}
		}
	}
	// A bad cycle shows up as a diagonal entry strictly better than the
	// tropical one (0): negative for min-plus, positive for max-plus.
	for i := 0; i < n; i++ {
		if d.data[i][i] != 0 && m.sr.AtLeastAsGood(d.data[i][i], 0) {
			return d, ErrDivergent
		}
	}
	return d, nil
}

// Star returns the tropical Kleene star A* = I (+) A (+) A^2 (+) ..., the
// reflexive-transitive closure of a square matrix. Over min-plus it is the
// all-pairs shortest-path matrix (empty walks allowed, so the diagonal is at
// most 0). It returns ErrDivergent if a bad cycle makes the closure diverge and
// ErrNotSquare for a non-square matrix.
func (m Matrix) Star() (Matrix, error) {
	plus, err := m.PlusClosure()
	if err != nil {
		return plus, err
	}
	for i := 0; i < m.rows; i++ {
		plus.data[i][i] = m.sr.Add(plus.data[i][i], 0)
	}
	return plus, nil
}

// HasBadCycle reports whether the matrix has a cycle that makes its tropical
// closure diverge: a negative-weight cycle for min-plus or a positive-weight
// cycle for max-plus. It returns ErrNotSquare for a non-square matrix.
func (m Matrix) HasBadCycle() (bool, error) {
	_, err := m.PlusClosure()
	if err == ErrDivergent {
		return true, nil
	}
	if err != nil {
		return false, err
	}
	return false, nil
}

// ShortestPaths returns the all-pairs shortest-path matrix of a min-plus
// weighted adjacency matrix, treating the tropical zero (+Inf) as the absence
// of an edge. Diagonal entries are 0 unless a shorter cyclic walk exists. It
// returns ErrDivergent when a negative cycle is present. It panics if the
// matrix is not min-plus.
func (m Matrix) ShortestPaths() (Matrix, error) {
	if !m.sr.IsMinPlus() {
		panic("tropical: ShortestPaths requires a min-plus matrix")
	}
	return m.Star()
}

// LongestPaths returns the all-pairs longest-path matrix of a max-plus weighted
// adjacency matrix, treating the tropical zero (-Inf) as the absence of an
// edge. It returns ErrDivergent when a positive cycle is present. It panics if
// the matrix is not max-plus.
func (m Matrix) LongestPaths() (Matrix, error) {
	if !m.sr.IsMaxPlus() {
		panic("tropical: LongestPaths requires a max-plus matrix")
	}
	return m.Star()
}

// Reachability returns a boolean matrix whose (i,j) entry reports whether j is
// reachable from i in at least one step, i.e. whether the plus-closure entry is
// not the tropical zero. Bad cycles do not affect reachability, so no error is
// returned for divergence.
func (m Matrix) Reachability() ([][]bool, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	n := m.rows
	reach := make([][]bool, n)
	for i := range reach {
		reach[i] = make([]bool, n)
		for j := 0; j < n; j++ {
			reach[i][j] = !m.sr.IsZero(m.data[i][j])
		}
	}
	for k := 0; k < n; k++ {
		for i := 0; i < n; i++ {
			if !reach[i][k] {
				continue
			}
			for j := 0; j < n; j++ {
				if reach[k][j] {
					reach[i][j] = true
				}
			}
		}
	}
	return reach, nil
}
