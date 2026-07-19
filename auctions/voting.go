package auctions

import "errors"

// WeightedVotingGame is a simple game [q; w_0,...,w_{n-1}]: a coalition wins
// when the total weight of its members reaches the quota Quota.
type WeightedVotingGame struct {
	// Quota is the threshold weight required for a coalition to win.
	Quota float64
	// Weights holds each player's voting weight.
	Weights []float64
}

// NewWeightedVotingGame builds a weighted voting game from a quota and weights.
// It returns an error if no weights are given or the quota is non-positive.
func NewWeightedVotingGame(quota float64, weights []float64) (WeightedVotingGame, error) {
	if len(weights) == 0 {
		return WeightedVotingGame{}, errors.New("auctions: weighted voting game needs at least one player")
	}
	if quota <= 0 {
		return WeightedVotingGame{}, errors.New("auctions: quota must be positive")
	}
	cp := make([]float64, len(weights))
	copy(cp, weights)
	return WeightedVotingGame{Quota: quota, Weights: cp}, nil
}

// Players returns the number of players in the game.
func (w WeightedVotingGame) Players() int { return len(w.Weights) }

// CoalitionWeight returns the total weight of the members of coalition s.
func (w WeightedVotingGame) CoalitionWeight(s Coalition) float64 {
	var total float64
	for i := range w.Weights {
		if s.Contains(i) {
			total += w.Weights[i]
		}
	}
	return total
}

// IsWinning reports whether coalition s meets or exceeds the quota.
func (w WeightedVotingGame) IsWinning(s Coalition) bool {
	return w.CoalitionWeight(s) >= w.Quota-1e-12
}

// CoopGame converts the weighted voting game into the corresponding simple
// cooperative game with worth 1 for winning coalitions and 0 otherwise.
func (w WeightedVotingGame) CoopGame() CoopGame {
	n := w.Players()
	size := 1 << uint(n)
	vals := make([]float64, size)
	for m := 0; m < size; m++ {
		if w.IsWinning(Coalition(m)) {
			vals[m] = 1
		}
	}
	return CoopGame{Players: n, values: vals}
}

// IsCritical reports whether player i is critical (a swing) in coalition s:
// s is winning but s without i is losing.
func (w WeightedVotingGame) IsCritical(i int, s Coalition) bool {
	if !s.Contains(i) {
		return false
	}
	return w.IsWinning(s) && !w.IsWinning(s.Remove(i))
}

// SwingCount returns the number of coalitions in which player i is critical.
func (w WeightedVotingGame) SwingCount(i int) int {
	n := w.Players()
	size := 1 << uint(n)
	count := 0
	for m := 0; m < size; m++ {
		if w.IsCritical(i, Coalition(m)) {
			count++
		}
	}
	return count
}

// BanzhafIndex returns the normalized Banzhaf power index of each player: their
// swing count divided by the total number of swings across all players.
func (w WeightedVotingGame) BanzhafIndex() []float64 {
	n := w.Players()
	swings := make([]int, n)
	total := 0
	for i := 0; i < n; i++ {
		swings[i] = w.SwingCount(i)
		total += swings[i]
	}
	out := make([]float64, n)
	if total == 0 {
		return out
	}
	for i := 0; i < n; i++ {
		out[i] = float64(swings[i]) / float64(total)
	}
	return out
}

// ShapleyShubikIndex returns the Shapley-Shubik power index of each player: the
// fraction of the n! orderings in which the player is pivotal (turns the
// preceding coalition from losing to winning).
func (w WeightedVotingGame) ShapleyShubikIndex() []float64 {
	return w.CoopGame().ShapleyValue()
}

// MinimalWinningCoalitions returns the winning coalitions that become losing if
// any single member leaves.
func (w WeightedVotingGame) MinimalWinningCoalitions() []Coalition {
	n := w.Players()
	size := 1 << uint(n)
	var out []Coalition
	for m := 1; m < size; m++ {
		s := Coalition(m)
		if !w.IsWinning(s) {
			continue
		}
		minimal := true
		for i := 0; i < n; i++ {
			if s.Contains(i) && w.IsWinning(s.Remove(i)) {
				minimal = false
				break
			}
		}
		if minimal {
			out = append(out, s)
		}
	}
	return out
}

// VetoPlayers returns the players present in every winning coalition; such a
// player can unilaterally block any decision.
func (w WeightedVotingGame) VetoPlayers() []int {
	n := w.Players()
	size := 1 << uint(n)
	var out []int
	for i := 0; i < n; i++ {
		veto := true
		for m := 0; m < size; m++ {
			s := Coalition(m)
			if w.IsWinning(s) && !s.Contains(i) {
				veto = false
				break
			}
		}
		if veto {
			out = append(out, i)
		}
	}
	return out
}

// Dummies returns the players who are never critical in any coalition.
func (w WeightedVotingGame) Dummies() []int {
	n := w.Players()
	var out []int
	for i := 0; i < n; i++ {
		if w.SwingCount(i) == 0 {
			out = append(out, i)
		}
	}
	return out
}

// Dictator returns the index of a dictator (a player who wins alone and whose
// absence makes every coalition lose) and true, or -1 and false if there is
// none.
func (w WeightedVotingGame) Dictator() (int, bool) {
	n := w.Players()
	for i := 0; i < n; i++ {
		if !w.IsWinning(SingletonCoalition(i)) {
			continue
		}
		isDictator := true
		size := 1 << uint(n)
		for m := 0; m < size; m++ {
			s := Coalition(m)
			if w.IsWinning(s) && !s.Contains(i) {
				isDictator = false
				break
			}
		}
		if isDictator {
			return i, true
		}
	}
	return -1, false
}
