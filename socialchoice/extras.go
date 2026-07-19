package socialchoice

import (
	"math"
	"sort"
)

// Clone returns a copy of the ballot.
func (b Ballot) Clone() Ballot {
	c := make(Ballot, len(b))
	copy(c, b)
	return c
}

// Copy returns a deep copy of the profile.
func (p *Profile) Copy() *Profile {
	ballots := make([]Ballot, len(p.Ballots))
	for i, b := range p.Ballots {
		ballots[i] = b.Clone()
	}
	counts := make([]int, len(p.Counts))
	copy(counts, p.Counts)
	return &Profile{Candidates: p.Candidates, Ballots: ballots, Counts: counts}
}

// RemoveCandidate returns a new profile identical to p but with candidate c
// removed from every ballot; the candidate count is unchanged so remaining
// indices stay stable (c simply never appears).
func (p *Profile) RemoveCandidate(c int) *Profile {
	q := &Profile{Candidates: p.Candidates, Counts: append([]int(nil), p.Counts...)}
	for _, b := range p.Ballots {
		nb := make(Ballot, 0, len(b))
		for _, x := range b {
			if x != c {
				nb = append(nb, x)
			}
		}
		q.Ballots = append(q.Ballots, nb)
	}
	return q
}

// PluralityScore returns the first-preference tally of a single candidate.
func (p *Profile) PluralityScore(c int) float64 {
	s := p.PluralityScores()
	if c < 0 || c >= len(s) {
		return 0
	}
	return s[c]
}

// BordaScore returns the Borda score of a single candidate.
func (p *Profile) BordaScore(c int) float64 {
	s := p.BordaScores()
	if c < 0 || c >= len(s) {
		return 0
	}
	return s[c]
}

// PairwiseMargin returns the signed majority margin of i over j in the profile.
func (p *Profile) PairwiseMargin(i, j int) int { return p.Pairwise().Margin(i, j) }

// AntiPluralityRanking returns candidates ordered by descending anti-plurality
// score.
func (p *Profile) AntiPluralityRanking() []int { return RankingFromScores(p.AntiPluralityScores()) }

// DowdallRanking returns candidates ordered by descending Dowdall score.
func (p *Profile) DowdallRanking() []int { return RankingFromScores(p.DowdallScores()) }

// PositionalRanking returns candidates ordered by descending positional score
// under the given weights.
func (p *Profile) PositionalRanking(weights []float64) []int {
	return RankingFromScores(p.PositionalScores(weights))
}

// MiniMaxScores returns the profile's minimax worst-defeat margins.
func (p *Profile) MiniMaxScores() []int { return p.Pairwise().MiniMaxScores() }

// CopelandNet returns the profile's net Copeland scores (wins minus losses).
func (p *Profile) CopelandNet() []int { return p.Pairwise().CopelandNet() }

// SchulzeWins returns the profile's per-candidate Schulze victory counts.
func (p *Profile) SchulzeWins() []int { return p.Pairwise().SchulzeWins() }

// Transpose returns the transpose of the pairwise matrix.
func (m PairwiseMatrix) Transpose() PairwiseMatrix {
	n := len(m)
	t := make(PairwiseMatrix, n)
	for i := 0; i < n; i++ {
		t[i] = make([]int, n)
		for j := 0; j < n; j++ {
			t[i][j] = m[j][i]
		}
	}
	return t
}

// NetMargins returns the antisymmetric margin matrix, entry [i][j] = margin of i
// over j.
func (m PairwiseMatrix) NetMargins() [][]int {
	n := len(m)
	d := make([][]int, n)
	for i := 0; i < n; i++ {
		d[i] = make([]int, n)
		for j := 0; j < n; j++ {
			d[i][j] = m[i][j] - m[j][i]
		}
	}
	return d
}

// HasTie reports whether any pair of distinct candidates is pairwise tied.
func (m PairwiseMatrix) HasTie() bool {
	n := len(m)
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if m.Ties(i, j) {
				return true
			}
		}
	}
	return false
}

// IsTournament reports whether the majority relation is a strict tournament: no
// pair of distinct candidates is tied.
func (m PairwiseMatrix) IsTournament() bool { return !m.HasTie() }

// HasCondorcetCycle reports whether the majority relation contains a cycle,
// equivalently whether the Smith set has more than one candidate.
func (m PairwiseMatrix) HasCondorcetCycle() bool { return len(m.SmithSet()) > 1 }

// HasCondorcetCycle reports whether the profile's majority relation is cyclic.
func (p *Profile) HasCondorcetCycle() bool { return p.Pairwise().HasCondorcetCycle() }

// CondorcetParadox reports whether the profile exhibits the Condorcet paradox:
// no Condorcet winner exists.
func (p *Profile) CondorcetParadox() bool { return !p.HasCondorcetWinner() }

// TopCycle returns the Smith set, also known as the top cycle of the majority
// tournament.
func (m PairwiseMatrix) TopCycle() []int { return m.SmithSet() }

// Covers reports whether candidate i covers candidate j: i beats j and beats
// every candidate that j beats.
func (m PairwiseMatrix) Covers(i, j int) bool {
	if i == j || !m.Beats(i, j) {
		return false
	}
	n := len(m)
	for k := 0; k < n; k++ {
		if k == i || k == j {
			continue
		}
		if m.Beats(j, k) && !m.Beats(i, k) {
			return false
		}
	}
	return true
}

// UncoveredSet returns the uncovered set: the candidates that no other candidate
// covers. The result is sorted ascending.
func (m PairwiseMatrix) UncoveredSet() []int {
	n := len(m)
	var set []int
	for x := 0; x < n; x++ {
		covered := false
		for y := 0; y < n; y++ {
			if y != x && m.Covers(y, x) {
				covered = true
				break
			}
		}
		if !covered {
			set = append(set, x)
		}
	}
	return set
}

// UncoveredSet returns the profile's uncovered set.
func (p *Profile) UncoveredSet() []int { return p.Pairwise().UncoveredSet() }

// BlackWinner returns the Black-method winner: the Condorcet winner when one
// exists, otherwise the Borda winner.
func (p *Profile) BlackWinner() int {
	if c, ok := p.CondorcetWinner(); ok {
		return c
	}
	return p.BordaWinner()
}

// permutations calls yield with every permutation of 0..n-1.
func permutations(n int, yield func(perm []int)) {
	perm := make([]int, n)
	for i := range perm {
		perm[i] = i
	}
	var rec func(k int)
	rec = func(k int) {
		if k == n {
			yield(perm)
			return
		}
		for i := k; i < n; i++ {
			perm[k], perm[i] = perm[i], perm[k]
			rec(k + 1)
			perm[k], perm[i] = perm[i], perm[k]
		}
	}
	rec(0)
}

// KemenyScore returns the Kemeny agreement score of a full ranking: the summed
// pairwise support for every ordered pair placed in agreement with the ranking.
// The Kemeny–Young optimal ranking maximises this score.
func (m PairwiseMatrix) KemenyScore(ranking []int) int {
	score := 0
	for a := 0; a < len(ranking); a++ {
		for b := a + 1; b < len(ranking); b++ {
			score += m[ranking[a]][ranking[b]]
		}
	}
	return score
}

// KemenyRanking returns a Kemeny–Young optimal ranking, found by exhaustive
// search over all orderings (suitable for a modest number of candidates). Ties
// are broken toward the lexicographically smallest ranking.
func (m PairwiseMatrix) KemenyRanking() []int {
	n := len(m)
	best := make([]int, n)
	for i := range best {
		best[i] = i
	}
	bestScore := m.KemenyScore(best)
	permutations(n, func(perm []int) {
		s := m.KemenyScore(perm)
		if s > bestScore {
			bestScore = s
			copy(best, perm)
		}
	})
	return best
}

// KemenyRanking returns the profile's Kemeny–Young optimal ranking.
func (p *Profile) KemenyRanking() []int { return p.Pairwise().KemenyRanking() }

// KemenyWinner returns the top candidate of the Kemeny–Young optimal ranking.
func (p *Profile) KemenyWinner() int {
	r := p.KemenyRanking()
	if len(r) == 0 {
		return -1
	}
	return r[0]
}

// NumBallots returns the number of approval ballots.
func (a ApprovalProfile) NumBallots() int { return len(a.Ballots) }

// TotalApprovals returns the total number of approvals cast across all ballots.
func (a ApprovalProfile) TotalApprovals() int {
	t := 0
	for _, b := range a.Ballots {
		for c := 0; c < a.Candidates && c < len(b); c++ {
			if b[c] {
				t++
			}
		}
	}
	return t
}

// NumBallots returns the number of score ballots.
func (s ScoreProfile) NumBallots() int { return len(s.Ballots) }

// ApprovalFromRanked derives an approval profile from a ranked profile by
// approving each ballot's top k candidates. Ballot weights are expanded into
// repeated approval ballots.
func ApprovalFromRanked(p *Profile, k int) ApprovalProfile {
	var ballots [][]bool
	for bi, b := range p.Ballots {
		row := make([]bool, p.Candidates)
		limit := k
		if limit > len(b) {
			limit = len(b)
		}
		for i := 0; i < limit; i++ {
			row[b[i]] = true
		}
		for w := 0; w < p.Weight(bi); w++ {
			cp := make([]bool, p.Candidates)
			copy(cp, row)
			ballots = append(ballots, cp)
		}
	}
	return ApprovalProfile{Candidates: p.Candidates, Ballots: ballots}
}

// RangeFromRanked derives a score profile from a ranked profile using Borda-style
// position scores: a candidate at zero-based rank r on an n-candidate ballot
// receives n-1-r points, and unranked candidates receive zero.
func RangeFromRanked(p *Profile) ScoreProfile {
	var ballots [][]float64
	maxScore := float64(p.Candidates - 1)
	for bi, b := range p.Ballots {
		row := make([]float64, p.Candidates)
		for r, c := range b {
			row[c] = float64(p.Candidates - 1 - r)
		}
		for w := 0; w < p.Weight(bi); w++ {
			cp := make([]float64, p.Candidates)
			copy(cp, row)
			ballots = append(ballots, cp)
		}
	}
	return ScoreProfile{Candidates: p.Candidates, Ballots: ballots, Max: maxScore}
}

// TotalSeats returns the sum of an allocation.
func TotalSeats(alloc []int) int { return sumInts(alloc) }

// GallagherIndex returns the Gallagher (least-squares) disproportionality index
// between vote shares and seat shares, expressed as a percentage in [0,100].
func GallagherIndex(votes, seats []int) float64 {
	totalV := float64(sumInts(votes))
	totalS := float64(sumInts(seats))
	if totalV == 0 || totalS == 0 {
		return 0
	}
	var sum float64
	n := len(votes)
	if len(seats) < n {
		n = len(seats)
	}
	for i := 0; i < n; i++ {
		vs := float64(votes[i]) / totalV * 100
		ss := float64(seats[i]) / totalS * 100
		d := vs - ss
		sum += d * d
	}
	return math.Sqrt(sum / 2)
}

// LoosemoreHanbyIndex returns the Loosemore–Hanby disproportionality index, half
// the summed absolute differences between vote and seat shares, as a percentage.
func LoosemoreHanbyIndex(votes, seats []int) float64 {
	totalV := float64(sumInts(votes))
	totalS := float64(sumInts(seats))
	if totalV == 0 || totalS == 0 {
		return 0
	}
	var sum float64
	n := len(votes)
	if len(seats) < n {
		n = len(seats)
	}
	for i := 0; i < n; i++ {
		vs := float64(votes[i]) / totalV * 100
		ss := float64(seats[i]) / totalS * 100
		sum += math.Abs(vs - ss)
	}
	return sum / 2
}

// EffectiveNumberOfParties returns the Laakso–Taagepera effective number of
// parties for the given vote (or seat) shares: the inverse of the sum of squared
// shares.
func EffectiveNumberOfParties(counts []int) float64 {
	total := float64(sumInts(counts))
	if total == 0 {
		return 0
	}
	var sumSq float64
	for _, c := range counts {
		s := float64(c) / total
		sumSq += s * s
	}
	if sumSq == 0 {
		return 0
	}
	return 1 / sumSq
}

// SortedInts returns a sorted copy of the slice, a small convenience used when
// comparing candidate sets.
func SortedInts(xs []int) []int {
	c := append([]int(nil), xs...)
	sort.Ints(c)
	return c
}
