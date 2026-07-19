package auctions

import (
	"errors"
	"math/rand"
)

// ErrPlayers is returned when a player count is outside the supported range.
var ErrPlayers = errors.New("auctions: player count must be between 1 and 20")

// CoopGame is a transferable-utility cooperative game on the fixed player set
// {0,...,Players-1}. The worth v(S) of every coalition S is precomputed and
// stored in a slice of length 2^Players indexed by the coalition bitmask.
type CoopGame struct {
	// Players is the number of players n.
	Players int
	values  []float64
}

// NewCoopGame builds a CoopGame by evaluating the characteristic function on
// every coalition bitmask. worth should return v(S); callers conventionally
// return 0 for the empty coalition. players must be in the range [1, 20].
func NewCoopGame(players int, worth func(Coalition) float64) (CoopGame, error) {
	if players < 1 || players > 20 {
		return CoopGame{}, ErrPlayers
	}
	size := 1 << uint(players)
	vals := make([]float64, size)
	for m := 0; m < size; m++ {
		vals[m] = worth(Coalition(m))
	}
	return CoopGame{Players: players, values: vals}, nil
}

// NewCoopGameValues builds a CoopGame from a precomputed table of coalition
// worths. values must have length exactly 2^players and be indexed by coalition
// bitmask.
func NewCoopGameValues(players int, values []float64) (CoopGame, error) {
	if players < 1 || players > 20 {
		return CoopGame{}, ErrPlayers
	}
	if len(values) != 1<<uint(players) {
		return CoopGame{}, errors.New("auctions: values length must equal 2^players")
	}
	cp := make([]float64, len(values))
	copy(cp, values)
	return CoopGame{Players: players, values: cp}, nil
}

// Value returns the worth v(S) of the coalition S.
func (g CoopGame) Value(s Coalition) float64 { return g.values[int(s)] }

// GrandCoalition returns the coalition of all players.
func (g CoopGame) GrandCoalition() Coalition { return FullCoalition(g.Players) }

// GrandValue returns v(N), the worth of the grand coalition.
func (g CoopGame) GrandValue() float64 { return g.values[len(g.values)-1] }

// MarginalContribution returns v(S ∪ {i}) - v(S), the value player i adds to
// coalition S (which need not exclude i; if i ∈ S the result is 0).
func (g CoopGame) MarginalContribution(i int, s Coalition) float64 {
	return g.Value(s.Add(i)) - g.Value(s.Remove(i))
}

// factorial returns k! as a float64.
func factorial(k int) float64 {
	r := 1.0
	for i := 2; i <= k; i++ {
		r *= float64(i)
	}
	return r
}

// ShapleyValue returns the Shapley value φ, the unique efficient allocation
// that averages each player's marginal contribution over all coalition
// orderings. φ_i = Σ_{S⊆N\{i}} |S|!(n-|S|-1)!/n! · (v(S∪{i}) - v(S)).
func (g CoopGame) ShapleyValue() []float64 {
	n := g.Players
	phi := make([]float64, n)
	nfact := factorial(n)
	size := 1 << uint(n)
	for i := 0; i < n; i++ {
		var sum float64
		for m := 0; m < size; m++ {
			s := Coalition(m)
			if s.Contains(i) {
				continue
			}
			ss := s.Size()
			w := factorial(ss) * factorial(n-ss-1) / nfact
			sum += w * (g.Value(s.Add(i)) - g.Value(s))
		}
		phi[i] = sum
	}
	return phi
}

// ShapleyValueMonteCarlo estimates the Shapley value by sampling random player
// orderings, using a math/rand source seeded by seed so that the estimate is
// reproducible. samples is the number of random permutations drawn.
func (g CoopGame) ShapleyValueMonteCarlo(samples int, seed int64) []float64 {
	n := g.Players
	phi := make([]float64, n)
	if samples <= 0 {
		return phi
	}
	rng := rand.New(rand.NewSource(seed))
	perm := make([]int, n)
	for i := range perm {
		perm[i] = i
	}
	for s := 0; s < samples; s++ {
		rng.Shuffle(n, func(a, b int) { perm[a], perm[b] = perm[b], perm[a] })
		var built Coalition
		for _, p := range perm {
			phi[p] += g.Value(built.Add(p)) - g.Value(built)
			built = built.Add(p)
		}
	}
	for i := range phi {
		phi[i] /= float64(samples)
	}
	return phi
}

// BanzhafValue returns the (raw, probabilistic) Banzhaf value: each player's
// average marginal contribution over the 2^(n-1) coalitions not containing it,
// β_i = (1/2^(n-1)) Σ_{S⊆N\{i}} (v(S∪{i}) - v(S)).
func (g CoopGame) BanzhafValue() []float64 {
	n := g.Players
	beta := make([]float64, n)
	denom := float64(int(1) << uint(n-1))
	size := 1 << uint(n)
	for i := 0; i < n; i++ {
		var sum float64
		for m := 0; m < size; m++ {
			s := Coalition(m)
			if s.Contains(i) {
				continue
			}
			sum += g.Value(s.Add(i)) - g.Value(s)
		}
		beta[i] = sum / denom
	}
	return beta
}

// BanzhafPowerIndex returns the normalized Banzhaf power index: the raw Banzhaf
// values rescaled to sum to one. If every raw value is zero the result is the
// zero vector.
func (g CoopGame) BanzhafPowerIndex() []float64 {
	raw := g.BanzhafValue()
	var total float64
	for _, r := range raw {
		total += r
	}
	out := make([]float64, len(raw))
	if total == 0 {
		return out
	}
	for i, r := range raw {
		out[i] = r / total
	}
	return out
}

// Excess returns the excess e(S, x) = v(S) - Σ_{i∈S} x_i of coalition S at the
// allocation x. A positive excess means S is dissatisfied with x.
func (g CoopGame) Excess(s Coalition, x []float64) float64 {
	sum := 0.0
	for i := 0; i < g.Players; i++ {
		if s.Contains(i) {
			sum += x[i]
		}
	}
	return g.Value(s) - sum
}

// MaxExcess returns the largest excess over all proper non-empty coalitions at
// the allocation x.
func (g CoopGame) MaxExcess(x []float64) float64 {
	n := g.Players
	size := 1 << uint(n)
	best := 0.0
	first := true
	for m := 1; m < size-1; m++ {
		e := g.Excess(Coalition(m), x)
		if first || e > best {
			best = e
			first = false
		}
	}
	return best
}

// IsEfficient reports whether the allocation x distributes exactly v(N) among
// the players, within the tolerance tol.
func (g CoopGame) IsEfficient(x []float64, tol float64) bool {
	sum := 0.0
	for _, v := range x {
		sum += v
	}
	d := sum - g.GrandValue()
	return d <= tol && d >= -tol
}

// IsImputation reports whether x is an imputation: efficient and individually
// rational (x_i >= v({i}) for every player i), within tolerance tol.
func (g CoopGame) IsImputation(x []float64, tol float64) bool {
	if !g.IsEfficient(x, tol) {
		return false
	}
	for i := 0; i < g.Players; i++ {
		if x[i] < g.Value(SingletonCoalition(i))-tol {
			return false
		}
	}
	return true
}

// InCore reports whether the allocation x lies in the core: efficient and
// coalitionally rational (Σ_{i∈S} x_i >= v(S) for every coalition S), within
// tolerance tol.
func (g CoopGame) InCore(x []float64, tol float64) bool {
	if !g.IsEfficient(x, tol) {
		return false
	}
	n := g.Players
	size := 1 << uint(n)
	for m := 1; m < size; m++ {
		if g.Excess(Coalition(m), x) > tol {
			return false
		}
	}
	return true
}
