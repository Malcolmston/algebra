package socialchoice

// PairwiseMatrix holds pairwise majority tallies: entry [i][j] is the number of
// voters who strictly prefer candidate i to candidate j.
type PairwiseMatrix [][]int

// Pairwise builds the PairwiseMatrix of the profile, summing ballot weights.
func (p *Profile) Pairwise() PairwiseMatrix {
	n := p.Candidates
	m := make(PairwiseMatrix, n)
	for i := range m {
		m[i] = make([]int, n)
	}
	for bi, b := range p.Ballots {
		w := p.Weight(bi)
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				if i == j {
					continue
				}
				if b.Prefers(i, j) {
					m[i][j] += w
				}
			}
		}
	}
	return m
}

// N returns the number of candidates covered by the matrix.
func (m PairwiseMatrix) N() int { return len(m) }

// Prefer returns the number of voters preferring i to j.
func (m PairwiseMatrix) Prefer(i, j int) int { return m[i][j] }

// Margin returns the signed majority margin of i over j, m[i][j]-m[j][i].
func (m PairwiseMatrix) Margin(i, j int) int { return m[i][j] - m[j][i] }

// Beats reports whether a strict majority prefers i to j.
func (m PairwiseMatrix) Beats(i, j int) bool { return m[i][j] > m[j][i] }

// Ties reports whether i and j are pairwise tied.
func (m PairwiseMatrix) Ties(i, j int) bool { return m[i][j] == m[j][i] }

// Wins returns the number of candidates that i beats pairwise.
func (m PairwiseMatrix) Wins(i int) int {
	c := 0
	for j := range m {
		if j != i && m.Beats(i, j) {
			c++
		}
	}
	return c
}

// Losses returns the number of candidates that beat i pairwise.
func (m PairwiseMatrix) Losses(i int) int {
	c := 0
	for j := range m {
		if j != i && m.Beats(j, i) {
			c++
		}
	}
	return c
}

// TieCount returns the number of candidates pairwise tied with i (excluding i).
func (m PairwiseMatrix) TieCount(i int) int {
	c := 0
	for j := range m {
		if j != i && m.Ties(i, j) {
			c++
		}
	}
	return c
}

// CondorcetWinner returns a candidate that beats every other candidate pairwise
// and true, or 0 and false when none exists.
func (m PairwiseMatrix) CondorcetWinner() (int, bool) {
	n := len(m)
	for i := 0; i < n; i++ {
		if m.Wins(i) == n-1 {
			return i, true
		}
	}
	return 0, false
}

// CondorcetLoser returns a candidate that loses to every other candidate
// pairwise and true, or 0 and false when none exists.
func (m PairwiseMatrix) CondorcetLoser() (int, bool) {
	n := len(m)
	for i := 0; i < n; i++ {
		if m.Losses(i) == n-1 {
			return i, true
		}
	}
	return 0, false
}

// CondorcetWinner returns the profile's Condorcet winner, if any.
func (p *Profile) CondorcetWinner() (int, bool) { return p.Pairwise().CondorcetWinner() }

// CondorcetLoser returns the profile's Condorcet loser, if any.
func (p *Profile) CondorcetLoser() (int, bool) { return p.Pairwise().CondorcetLoser() }

// HasCondorcetWinner reports whether the profile has a Condorcet winner.
func (p *Profile) HasCondorcetWinner() bool {
	_, ok := p.CondorcetWinner()
	return ok
}

// IsCondorcetWinner reports whether candidate c beats every other candidate
// pairwise in the profile.
func (p *Profile) IsCondorcetWinner(c int) bool {
	m := p.Pairwise()
	return c >= 0 && c < len(m) && m.Wins(c) == len(m)-1
}

// CopelandScores returns each candidate's Copeland score, counting 1 for every
// pairwise win and 0.5 for every pairwise tie.
func (m PairwiseMatrix) CopelandScores() []float64 {
	n := len(m)
	s := make([]float64, n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if i == j {
				continue
			}
			switch {
			case m.Beats(i, j):
				s[i]++
			case m.Ties(i, j):
				s[i] += 0.5
			}
		}
	}
	return s
}

// CopelandNet returns each candidate's net Copeland score, wins minus losses.
func (m PairwiseMatrix) CopelandNet() []int {
	n := len(m)
	s := make([]int, n)
	for i := 0; i < n; i++ {
		s[i] = m.Wins(i) - m.Losses(i)
	}
	return s
}

// CopelandScores returns the profile's Copeland scores (win=1, tie=0.5).
func (p *Profile) CopelandScores() []float64 { return p.Pairwise().CopelandScores() }

// CopelandWinner returns the highest Copeland-scoring candidate.
func (p *Profile) CopelandWinner() int { return ArgMaxFloat(p.CopelandScores()) }

// CopelandRanking returns candidates ordered by descending Copeland score.
func (p *Profile) CopelandRanking() []int { return RankingFromScores(p.CopelandScores()) }

// MiniMaxScores returns, for each candidate, the largest margin by which it is
// defeated in any pairwise contest (its worst pairwise loss); a value of zero or
// below means the candidate is never beaten. This is the margins variant of the
// Simpson–Kramer minimax rule.
func (m PairwiseMatrix) MiniMaxScores() []int {
	n := len(m)
	worst := make([]int, n)
	for i := 0; i < n; i++ {
		w := 0
		for j := 0; j < n; j++ {
			if i == j {
				continue
			}
			if d := m.Margin(j, i); d > w {
				w = d
			}
		}
		worst[i] = w
	}
	return worst
}

// MiniMaxWinner returns the candidate whose worst pairwise defeat margin is
// smallest, breaking ties toward the lower index.
func (p *Profile) MiniMaxWinner() int {
	worst := p.Pairwise().MiniMaxScores()
	best := 0
	for i := 1; i < len(worst); i++ {
		if worst[i] < worst[best] {
			best = i
		}
	}
	if len(worst) == 0 {
		return -1
	}
	return best
}

// transitiveClosure returns the reflexive-free transitive closure of the boolean
// adjacency matrix adj via the Floyd–Warshall algorithm. reach[i][j] is true
// when j is reachable from i through one or more edges.
func transitiveClosure(adj [][]bool) [][]bool {
	n := len(adj)
	reach := make([][]bool, n)
	for i := range reach {
		reach[i] = make([]bool, n)
		copy(reach[i], adj[i])
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
	return reach
}

// SmithSet returns the Smith set: the smallest non-empty set of candidates such
// that every member beats every non-member pairwise. It is computed as the
// maximal elements of the transitive closure of the "at least as good as"
// relation (i does not lose to j). The result is sorted ascending.
func (m PairwiseMatrix) SmithSet() []int {
	n := len(m)
	adj := make([][]bool, n)
	for i := 0; i < n; i++ {
		adj[i] = make([]bool, n)
		for j := 0; j < n; j++ {
			if i == j {
				continue
			}
			// i is at least as good as j when it does not lose to j.
			adj[i][j] = m[i][j] >= m[j][i]
		}
	}
	reach := transitiveClosure(adj)
	var set []int
	for i := 0; i < n; i++ {
		all := true
		for j := 0; j < n; j++ {
			if i == j {
				continue
			}
			if !reach[i][j] {
				all = false
				break
			}
		}
		if all {
			set = append(set, i)
		}
	}
	return set
}

// SchwartzSet returns the Schwartz set: the union of the minimal undominated
// sets under the strict beat relation (ties do not dominate). The result is
// sorted ascending.
func (m PairwiseMatrix) SchwartzSet() []int {
	n := len(m)
	adj := make([][]bool, n)
	for i := 0; i < n; i++ {
		adj[i] = make([]bool, n)
		for j := 0; j < n; j++ {
			if i != j {
				adj[i][j] = m.Beats(i, j)
			}
		}
	}
	reach := transitiveClosure(adj)
	var set []int
	for x := 0; x < n; x++ {
		inSchwartz := true
		for y := 0; y < n; y++ {
			if y == x {
				continue
			}
			// If y can reach x but x cannot reach y, x's component is dominated.
			if reach[y][x] && !reach[x][y] {
				inSchwartz = false
				break
			}
		}
		if inSchwartz {
			set = append(set, x)
		}
	}
	return set
}

// SmithSet returns the profile's Smith set.
func (p *Profile) SmithSet() []int { return p.Pairwise().SmithSet() }

// SchwartzSet returns the profile's Schwartz set.
func (p *Profile) SchwartzSet() []int { return p.Pairwise().SchwartzSet() }

// InSet reports whether candidate c appears in the (sorted or unsorted) set.
func InSet(set []int, c int) bool {
	for _, x := range set {
		if x == c {
			return true
		}
	}
	return false
}
