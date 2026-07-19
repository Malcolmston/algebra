package socialchoice

import "sort"

// MajorityPair is a directed pairwise victory of Winner over Loser with the
// given majority margin and winning-vote count, used by ranked pairs.
type MajorityPair struct {
	Winner int
	Loser  int
	Margin int
	Votes  int
}

// SortedMajorities returns every pairwise victory in the matrix ordered by
// descending margin; ties in margin are broken by descending winning votes and
// then by the (Winner, Loser) index pair to make the order deterministic.
func (m PairwiseMatrix) SortedMajorities() []MajorityPair {
	n := len(m)
	var pairs []MajorityPair
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if i != j && m.Beats(i, j) {
				pairs = append(pairs, MajorityPair{
					Winner: i,
					Loser:  j,
					Margin: m.Margin(i, j),
					Votes:  m[i][j],
				})
			}
		}
	}
	sort.SliceStable(pairs, func(a, b int) bool {
		pa, pb := pairs[a], pairs[b]
		if pa.Margin != pb.Margin {
			return pa.Margin > pb.Margin
		}
		if pa.Votes != pb.Votes {
			return pa.Votes > pb.Votes
		}
		if pa.Winner != pb.Winner {
			return pa.Winner < pb.Winner
		}
		return pa.Loser < pb.Loser
	})
	return pairs
}

// pathExists reports whether j is reachable from i in the locked-edge adjacency.
func pathExists(locked [][]bool, i, j int) bool {
	n := len(locked)
	visited := make([]bool, n)
	stack := []int{i}
	for len(stack) > 0 {
		x := stack[len(stack)-1]
		stack = stack[:len(stack)-1]
		if x == j {
			return true
		}
		if visited[x] {
			continue
		}
		visited[x] = true
		for k := 0; k < n; k++ {
			if locked[x][k] && !visited[k] {
				stack = append(stack, k)
			}
		}
	}
	return false
}

// RankedPairsLocked returns the acyclic set of locked pairwise victories under
// Tideman's ranked-pairs method: majorities are considered strongest first and a
// pair is locked in unless doing so would create a cycle with the pairs already
// locked. locked[i][j] is true when the edge i→j is locked.
func (m PairwiseMatrix) RankedPairsLocked() [][]bool {
	n := len(m)
	locked := make([][]bool, n)
	for i := range locked {
		locked[i] = make([]bool, n)
	}
	for _, pr := range m.SortedMajorities() {
		// Lock pr.Winner→pr.Loser unless pr.Loser already reaches pr.Winner.
		if !pathExists(locked, pr.Loser, pr.Winner) {
			locked[pr.Winner][pr.Loser] = true
		}
	}
	return locked
}

// RankedPairsRanking returns the ranked-pairs ordering: a topological order of
// the locked-victory graph, most preferred first. Ties are broken toward the
// lower index.
func (m PairwiseMatrix) RankedPairsRanking() []int {
	n := len(m)
	locked := m.RankedPairsLocked()
	indeg := make([]int, n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if locked[i][j] {
				indeg[j]++
			}
		}
	}
	placed := make([]bool, n)
	order := make([]int, 0, n)
	for len(order) < n {
		next := -1
		for i := 0; i < n; i++ {
			if !placed[i] && indeg[i] == 0 {
				next = i
				break
			}
		}
		if next == -1 {
			// Should not happen for an acyclic graph; place remaining by index.
			for i := 0; i < n; i++ {
				if !placed[i] {
					next = i
					break
				}
			}
		}
		placed[next] = true
		order = append(order, next)
		for j := 0; j < n; j++ {
			if locked[next][j] {
				indeg[j]--
			}
		}
	}
	return order
}

// RankedPairsWinner returns the source of the locked-victory graph, the
// ranked-pairs winner.
func (m PairwiseMatrix) RankedPairsWinner() int {
	order := m.RankedPairsRanking()
	if len(order) == 0 {
		return -1
	}
	return order[0]
}

// RankedPairsWinner returns the profile's ranked-pairs (Tideman) winner.
func (p *Profile) RankedPairsWinner() int { return p.Pairwise().RankedPairsWinner() }

// RankedPairsRanking returns the profile's ranked-pairs ordering.
func (p *Profile) RankedPairsRanking() []int { return p.Pairwise().RankedPairsRanking() }
