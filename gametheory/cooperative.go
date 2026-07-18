package gametheory

import (
	"errors"
	"math/bits"
)

// CooperativeGame is a transferable-utility cooperative game on a fixed set of
// players {0, 1, ..., Players-1}. A coalition is encoded as a bitmask whose bit
// i is set when player i belongs to the coalition. The worth of every coalition
// is precomputed and stored in a slice of length 2^Players indexed by that
// bitmask, so lookups and the solution concepts run without re-invoking the
// characteristic function.
type CooperativeGame struct {
	// Players is the number of players n.
	Players int
	// values[mask] is the worth v(S) of the coalition encoded by mask.
	values []float64
}

// ErrPlayers is returned when a player count is out of the supported range.
var ErrPlayers = errors.New("gametheory: player count must be between 1 and 20")

// NewCooperativeGame builds a CooperativeGame with the given number of players
// by evaluating the characteristic function on every one of the 2^players
// coalition bitmasks. The characteristic function should return the worth
// v(coalition); the empty coalition's value is stored as evaluated (callers
// conventionally return 0 for it). players must be between 1 and 20 inclusive.
func NewCooperativeGame(players int, characteristic func(coalition uint) float64) (CooperativeGame, error) {
	if players < 1 || players > 20 {
		return CooperativeGame{}, ErrPlayers
	}
	size := 1 << players
	values := make([]float64, size)
	for mask := 0; mask < size; mask++ {
		values[mask] = characteristic(uint(mask))
	}
	return CooperativeGame{Players: players, values: values}, nil
}

// NewCooperativeGameFromValues builds a CooperativeGame from a fully specified
// worth table, where values[mask] is v(mask) for each coalition bitmask. The
// slice length must equal 2^players for some players in [1, 20]. The slice is
// copied.
func NewCooperativeGameFromValues(values []float64) (CooperativeGame, error) {
	n := bits.TrailingZeros(uint(len(values)))
	if len(values) == 0 || (1<<n) != len(values) || n < 1 || n > 20 {
		return CooperativeGame{}, ErrPlayers
	}
	return CooperativeGame{Players: n, values: append([]float64(nil), values...)}, nil
}

// Worth returns v(coalition), the worth of the coalition encoded by the given
// bitmask.
func (cg CooperativeGame) Worth(coalition uint) float64 {
	return cg.values[coalition]
}

// GrandCoalitionValue returns v(N), the worth of the coalition of all players.
func (cg CooperativeGame) GrandCoalitionValue() float64 {
	return cg.values[(1<<cg.Players)-1]
}

// MarginalContribution returns v(coalition ∪ {player}) - v(coalition), the value
// player adds when joining the coalition. The coalition must not already
// contain player.
func (cg CooperativeGame) MarginalContribution(player int, coalition uint) float64 {
	with := coalition | (1 << uint(player))
	return cg.values[with] - cg.values[coalition]
}

// ShapleyValue returns the Shapley value of the game: a fair allocation of the
// grand coalition's worth in which each player receives the probability-weighted
// average of their marginal contributions over all coalitions they may join,
//
//	φ_i = Σ_{S ⊆ N\{i}} |S|! (n-|S|-1)! / n! · (v(S ∪ {i}) - v(S)).
//
// The returned slice has one entry per player and, for any game, sums to v(N).
func (cg CooperativeGame) ShapleyValue() []float64 {
	n := cg.Players
	phi := make([]float64, n)
	// weight[s] = s! (n-s-1)! / n!, the probability weight of a size-s coalition.
	weight := make([]float64, n)
	fact := gametheoryFactorials(n)
	for s := 0; s < n; s++ {
		weight[s] = float64(fact[s]) * float64(fact[n-s-1]) / float64(fact[n])
	}
	full := 1 << n
	for i := 0; i < n; i++ {
		bit := uint(1) << uint(i)
		var sum float64
		for mask := 0; mask < full; mask++ {
			if uint(mask)&bit != 0 {
				continue // S must exclude i
			}
			s := bits.OnesCount(uint(mask))
			marginal := cg.values[uint(mask)|bit] - cg.values[mask]
			sum += weight[s] * marginal
		}
		phi[i] = sum
	}
	return phi
}

// BanzhafValue returns the (non-normalized) Banzhaf value of the game: for each
// player, the average marginal contribution taken with equal weight over all
// 2^(n-1) coalitions not containing that player,
//
//	β_i = 1 / 2^(n-1) · Σ_{S ⊆ N\{i}} (v(S ∪ {i}) - v(S)).
//
// Unlike the Shapley value, the Banzhaf value need not sum to v(N).
func (cg CooperativeGame) BanzhafValue() []float64 {
	n := cg.Players
	beta := make([]float64, n)
	denom := float64(int(1) << (n - 1))
	full := 1 << n
	for i := 0; i < n; i++ {
		bit := uint(1) << uint(i)
		var sum float64
		for mask := 0; mask < full; mask++ {
			if uint(mask)&bit != 0 {
				continue
			}
			sum += cg.values[uint(mask)|bit] - cg.values[mask]
		}
		beta[i] = sum / denom
	}
	return beta
}

// IsSuperadditive reports whether the game is superadditive: for every pair of
// disjoint coalitions S and T, v(S ∪ T) >= v(S) + v(T) within tolerance tol.
func (cg CooperativeGame) IsSuperadditive(tol float64) bool {
	full := 1 << cg.Players
	for s := 0; s < full; s++ {
		// Iterate over subsets t of the complement of s.
		comp := (full - 1) &^ s
		for t := comp; ; t = (t - 1) & comp {
			if t != 0 {
				if cg.values[s|t] < cg.values[s]+cg.values[t]-tol {
					return false
				}
			}
			if t == 0 {
				break
			}
		}
	}
	return true
}

// IsMonotone reports whether the game is monotone: enlarging a coalition never
// decreases its worth, i.e. v(S) <= v(T) whenever S ⊆ T, within tolerance tol.
func (cg CooperativeGame) IsMonotone(tol float64) bool {
	full := 1 << cg.Players
	for s := 0; s < full; s++ {
		comp := (full - 1) &^ s
		for add := comp; ; add = (add - 1) & comp {
			if add != 0 {
				if cg.values[s|add] < cg.values[s]-tol {
					return false
				}
			}
			if add == 0 {
				break
			}
		}
	}
	return true
}

// IsConvex reports whether the game is convex (supermodular): for all coalitions
// S and T, v(S ∪ T) + v(S ∩ T) >= v(S) + v(T) within tolerance tol. Convex games
// always have a non-empty core containing the Shapley value.
func (cg CooperativeGame) IsConvex(tol float64) bool {
	full := 1 << cg.Players
	for s := 0; s < full; s++ {
		for t := 0; t < full; t++ {
			union := s | t
			inter := s & t
			if cg.values[union]+cg.values[inter] < cg.values[s]+cg.values[t]-tol {
				return false
			}
		}
	}
	return true
}

// IsInCore reports whether allocation x is in the core of the game: it must be
// efficient (Σ x_i = v(N)) and coalitionally rational (Σ_{i∈S} x_i >= v(S) for
// every coalition S), all within tolerance tol. len(x) must equal Players.
func (cg CooperativeGame) IsInCore(x []float64, tol float64) bool {
	if len(x) != cg.Players {
		return false
	}
	full := 1 << cg.Players
	var total float64
	for _, v := range x {
		total += v
	}
	if total < cg.GrandCoalitionValue()-tol || total > cg.GrandCoalitionValue()+tol {
		return false
	}
	for mask := 1; mask < full; mask++ {
		var sum float64
		for i := 0; i < cg.Players; i++ {
			if mask&(1<<i) != 0 {
				sum += x[i]
			}
		}
		if sum < cg.values[mask]-tol {
			return false
		}
	}
	return true
}

// IsEssential reports whether the game is essential: the grand coalition is
// worth strictly more (by more than tol) than the players acting alone, i.e.
// v(N) > Σ_i v({i}). Only essential games have a non-trivial solution.
func (cg CooperativeGame) IsEssential(tol float64) bool {
	var soloSum float64
	for i := 0; i < cg.Players; i++ {
		soloSum += cg.values[1<<uint(i)]
	}
	return cg.GrandCoalitionValue() > soloSum+tol
}

// WeightedVotingGame constructs the simple cooperative game in which a coalition
// wins (worth 1) when the total weight of its members is at least quota, and
// loses (worth 0) otherwise. weights must have between 1 and 20 entries.
func WeightedVotingGame(weights []float64, quota float64) (CooperativeGame, error) {
	n := len(weights)
	if n < 1 || n > 20 {
		return CooperativeGame{}, ErrPlayers
	}
	w := append([]float64(nil), weights...)
	return NewCooperativeGame(n, func(coalition uint) float64 {
		var sum float64
		for i := 0; i < n; i++ {
			if coalition&(1<<uint(i)) != 0 {
				sum += w[i]
			}
		}
		if sum >= quota {
			return 1
		}
		return 0
	})
}

// gametheoryFactorials returns the slice of factorials 0!..n! as int64.
func gametheoryFactorials(n int) []int64 {
	f := make([]int64, n+1)
	f[0] = 1
	for i := 1; i <= n; i++ {
		f[i] = f[i-1] * int64(i)
	}
	return f
}
