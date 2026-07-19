package auctions

import (
	"errors"
	"math"
	"math/rand"
	"sort"
)

// EqualSplit returns the egalitarian allocation that divides v(N) equally among
// the players.
func (g CoopGame) EqualSplit() []float64 {
	n := g.Players
	out := make([]float64, n)
	share := g.GrandValue() / float64(n)
	for i := range out {
		out[i] = share
	}
	return out
}

// IsAdditive reports whether the game is additive (inessential):
// v(S) = Σ_{i∈S} v({i}) for every coalition S.
func (g CoopGame) IsAdditive() bool {
	n := g.Players
	size := 1 << uint(n)
	single := make([]float64, n)
	for i := 0; i < n; i++ {
		single[i] = g.Value(SingletonCoalition(i))
	}
	for m := 0; m < size; m++ {
		s := Coalition(m)
		var sum float64
		for i := 0; i < n; i++ {
			if s.Contains(i) {
				sum += single[i]
			}
		}
		if !approxEqual(g.Value(s), sum, 1e-9) {
			return false
		}
	}
	return true
}

// IndividualRational reports whether every player receives at least their
// singleton worth under x (within tolerance tol).
func (g CoopGame) IndividualRational(x []float64, tol float64) bool {
	for i := 0; i < g.Players; i++ {
		if x[i] < g.Value(SingletonCoalition(i))-tol {
			return false
		}
	}
	return true
}

// UtopiaPayoff returns the utopia (marginal) payoff vector M with
// M_i = v(N) - v(N\{i}), the most player i could hope to receive if every other
// player were paid their marginal contribution to the grand coalition.
func (g CoopGame) UtopiaPayoff() []float64 {
	n := g.Players
	full := FullCoalition(n)
	vN := g.GrandValue()
	out := make([]float64, n)
	for i := 0; i < n; i++ {
		out[i] = vN - g.Value(full.Remove(i))
	}
	return out
}

// MinimalRights returns the minimal-rights (concession) vector m used by the
// tau-value: m_i = max over coalitions S containing i of
// v(S) - Σ_{j∈S, j≠i} M_j, where M is the utopia payoff.
func (g CoopGame) MinimalRights() []float64 {
	n := g.Players
	m := g.UtopiaPayoff()
	size := 1 << uint(n)
	out := make([]float64, n)
	for i := 0; i < n; i++ {
		best := math.Inf(-1)
		for mask := 0; mask < size; mask++ {
			s := Coalition(mask)
			if !s.Contains(i) {
				continue
			}
			val := g.Value(s)
			for j := 0; j < n; j++ {
				if j != i && s.Contains(j) {
					val -= m[j]
				}
			}
			if val > best {
				best = val
			}
		}
		out[i] = best
	}
	return out
}

// TauValue returns the tau-value (compromise value) of the game: the efficient
// allocation on the segment between the minimal-rights vector m and the utopia
// payoff M, tau = m + t(M - m) with t chosen so that Σ tau = v(N). It requires
// the game to be quasi-balanced (m_i <= M_i for all i and Σ m <= v(N) <= Σ M).
func (g CoopGame) TauValue() ([]float64, error) {
	n := g.Players
	M := g.UtopiaPayoff()
	m := g.MinimalRights()
	var sumM, summ float64
	for i := 0; i < n; i++ {
		if m[i] > M[i]+1e-9 {
			return nil, errors.New("auctions: game is not quasi-balanced (m_i > M_i)")
		}
		sumM += M[i]
		summ += m[i]
	}
	vN := g.GrandValue()
	if vN < summ-1e-9 || vN > sumM+1e-9 {
		return nil, errors.New("auctions: game is not quasi-balanced (v(N) outside [Σm, ΣM])")
	}
	denom := sumM - summ
	t := 0.0
	if math.Abs(denom) > 1e-12 {
		t = (vN - summ) / denom
	}
	out := make([]float64, n)
	for i := 0; i < n; i++ {
		out[i] = m[i] + t*(M[i]-m[i])
	}
	return out, nil
}

// SortedExcesses returns the excesses of every proper non-empty coalition at
// the allocation x, sorted in decreasing order. The nucleolus lexicographically
// minimizes this vector.
func (g CoopGame) SortedExcesses(x []float64) []float64 {
	n := g.Players
	size := 1 << uint(n)
	out := make([]float64, 0, size-2)
	for m := 1; m < size-1; m++ {
		out = append(out, g.Excess(Coalition(m), x))
	}
	sort.Sort(sort.Reverse(sort.Float64Slice(out)))
	return out
}

// ShapleyShubikIndex returns the Shapley-Shubik power index of a simple game,
// which equals its Shapley value.
func (g CoopGame) ShapleyShubikIndex() []float64 { return g.ShapleyValue() }

// Add returns the game v + w, the coalition-wise sum of two games on the same
// player set.
func (g CoopGame) Add(other CoopGame) (CoopGame, error) {
	if g.Players != other.Players {
		return CoopGame{}, errors.New("auctions: games have different player counts")
	}
	vals := make([]float64, len(g.values))
	for i := range vals {
		vals[i] = g.values[i] + other.values[i]
	}
	return CoopGame{Players: g.Players, values: vals}, nil
}

// Scale returns the game c·v obtained by multiplying every coalition worth by c.
func (g CoopGame) Scale(c float64) CoopGame {
	vals := make([]float64, len(g.values))
	for i := range vals {
		vals[i] = c * g.values[i]
	}
	return CoopGame{Players: g.Players, values: vals}
}

// SortBidsDescending returns a copy of the bids sorted by decreasing value, ties
// broken by increasing bidder index.
func SortBidsDescending(bids []Bid) []Bid { return sortedByValue(bids) }

// TotalBidValue returns the sum of all bid values.
func TotalBidValue(bids []Bid) float64 {
	var total float64
	for _, b := range bids {
		total += b.Value
	}
	return total
}

// AllocativeEfficiency returns the ratio of the realised valuation to the
// maximum possible valuation for a single-item auction: 1 when the item is
// awarded to a highest-value bidder, less otherwise. valuations is indexed by
// bidder. It returns 0 when no valuations are supplied.
func AllocativeEfficiency(valuations map[int]float64, out AuctionOutcome) float64 {
	var best float64
	first := true
	for _, v := range valuations {
		if first || v > best {
			best = v
			first = false
		}
	}
	if first || best == 0 {
		return 0
	}
	if out.Winner < 0 {
		return 0
	}
	return valuations[out.Winner] / best
}

// FeasibleContains reports whether the point p lies inside (or on the boundary
// of) the convex hull of the feasible set pts.
func FeasibleContains(pts []Point, p Point) bool {
	hull := ConvexHull2D(pts)
	m := len(hull)
	if m == 1 {
		return hull[0] == p
	}
	if m == 2 {
		// collinear: check p on the segment's bounding box and line
		a, b := hull[0], hull[1]
		if math.Abs(cross(a, b, p)) > 1e-9 {
			return false
		}
		return p.X >= math.Min(a.X, b.X)-1e-9 && p.X <= math.Max(a.X, b.X)+1e-9 &&
			p.Y >= math.Min(a.Y, b.Y)-1e-9 && p.Y <= math.Max(a.Y, b.Y)+1e-9
	}
	for i := 0; i < m; i++ {
		a := hull[i]
		b := hull[(i+1)%m]
		// CCW polygon: interior is to the left of each directed edge.
		if cross(a, b, p) < -1e-9 {
			return false
		}
	}
	return true
}

// WeightedNashBargainingSolution returns the asymmetric Nash bargaining
// solution with bargaining power alpha in (0,1) for player 0 and 1-alpha for
// player 1: the feasible point maximizing (x-d.X)^alpha (y-d.Y)^(1-alpha).
func WeightedNashBargainingSolution(pts []Point, d Point, alpha float64) (Point, error) {
	if alpha <= 0 || alpha >= 1 {
		return Point{}, errors.New("auctions: alpha must lie strictly between 0 and 1")
	}
	hull := ConvexHull2D(pts)
	if len(hull) == 0 {
		return Point{}, errors.New("auctions: empty feasible set")
	}
	best := Point{}
	bestVal := math.Inf(-1)
	found := false
	obj := func(p Point) (float64, bool) {
		u := p.X - d.X
		w := p.Y - d.Y
		if u <= 0 || w <= 0 {
			return 0, false
		}
		return alpha*math.Log(u) + (1-alpha)*math.Log(w), true
	}
	consider := func(p Point) {
		if v, ok := obj(p); ok && v > bestVal {
			bestVal = v
			best = p
			found = true
		}
	}
	m := len(hull)
	for i := 0; i < m; i++ {
		a := hull[i]
		b := hull[(i+1)%m]
		dx := b.X - a.X
		dy := b.Y - a.Y
		ax := a.X - d.X
		ay := a.Y - d.Y
		if math.Abs(dx*dy) > 1e-15 {
			t := -(alpha*dx*ay + (1-alpha)*dy*ax) / (dx * dy)
			if t >= 0 && t <= 1 {
				consider(Point{X: a.X + t*dx, Y: a.Y + t*dy})
			}
		}
		consider(a)
		consider(b)
	}
	if !found {
		return Point{}, errors.New("auctions: no feasible point strictly dominates the disagreement point")
	}
	return best, nil
}

// AgentRank returns the position (0 = most preferred) of item in agent's
// preference list, or -1 if the item is not listed.
func AgentRank(prefs [][]int, agent, item int) int {
	if agent < 0 || agent >= len(prefs) {
		return -1
	}
	for pos, it := range prefs[agent] {
		if it == item {
			return pos
		}
	}
	return -1
}

// RandomSerialDictatorship runs serial dictatorship with a uniformly random
// priority order drawn from a math/rand source seeded by seed, so the outcome
// is reproducible. It returns the assignment and the order used.
func RandomSerialDictatorship(prefs [][]int, seed int64) ([]int, []int, error) {
	n := len(prefs)
	order := make([]int, n)
	for i := range order {
		order[i] = i
	}
	rng := rand.New(rand.NewSource(seed))
	rng.Shuffle(n, func(i, j int) { order[i], order[j] = order[j], order[i] })
	assignment, err := SerialDictatorship(prefs, order)
	if err != nil {
		return nil, nil, err
	}
	return assignment, order, nil
}
