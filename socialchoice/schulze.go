package socialchoice

import "sort"

// StrongestPaths returns the Schulze strongest-beatpath matrix p, where p[i][j]
// is the strength of the strongest path from i to j. A direct link i→j has
// strength m[i][j] when i beats j pairwise and 0 otherwise; a path's strength is
// its weakest link, and the matrix maximises over all paths via Floyd–Warshall.
func (m PairwiseMatrix) StrongestPaths() [][]int {
	n := len(m)
	p := make([][]int, n)
	for i := 0; i < n; i++ {
		p[i] = make([]int, n)
		for j := 0; j < n; j++ {
			if i != j && m[i][j] > m[j][i] {
				p[i][j] = m[i][j]
			}
		}
	}
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if i == j {
				continue
			}
			for k := 0; k < n; k++ {
				if i == k || j == k {
					continue
				}
				s := p[j][i]
				if p[i][k] < s {
					s = p[i][k]
				}
				if s > p[j][k] {
					p[j][k] = s
				}
			}
		}
	}
	return p
}

// SchulzeWins returns each candidate's number of Schulze pairwise victories: the
// count of opponents j for which the strongest path from i to j is stronger than
// the strongest path from j to i.
func (m PairwiseMatrix) SchulzeWins() []int {
	p := m.StrongestPaths()
	n := len(m)
	wins := make([]int, n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if i != j && p[i][j] > p[j][i] {
				wins[i]++
			}
		}
	}
	return wins
}

// SchulzeWinner returns a Schulze-method winner: a candidate whose strongest path
// to every other candidate is at least as strong as the reverse. Ties are
// resolved toward the lower index.
func (m PairwiseMatrix) SchulzeWinner() int {
	p := m.StrongestPaths()
	n := len(m)
	for i := 0; i < n; i++ {
		ok := true
		for j := 0; j < n; j++ {
			if i != j && p[j][i] > p[i][j] {
				ok = false
				break
			}
		}
		if ok {
			return i
		}
	}
	return -1
}

// SchulzeRanking returns candidates ordered by descending number of Schulze
// pairwise victories, breaking ties toward the lower index.
func (m PairwiseMatrix) SchulzeRanking() []int {
	wins := m.SchulzeWins()
	idx := make([]int, len(wins))
	for i := range idx {
		idx[i] = i
	}
	sort.SliceStable(idx, func(a, b int) bool {
		if wins[idx[a]] != wins[idx[b]] {
			return wins[idx[a]] > wins[idx[b]]
		}
		return idx[a] < idx[b]
	})
	return idx
}

// SchulzeWinner returns the profile's Schulze-method winner.
func (p *Profile) SchulzeWinner() int { return p.Pairwise().SchulzeWinner() }

// SchulzeRanking returns the profile's Schulze ranking.
func (p *Profile) SchulzeRanking() []int { return p.Pairwise().SchulzeRanking() }
