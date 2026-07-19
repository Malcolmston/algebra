package socialchoice

import "math"

// MajorityThreshold returns the smallest tally that constitutes a strict
// majority of total voters, floor(total/2)+1.
func MajorityThreshold(total int) int { return total/2 + 1 }

// HasMajorityWinner reports whether some candidate holds a strict majority of
// first preferences, returning that candidate and true when so.
func (p *Profile) HasMajorityWinner() (int, bool) {
	scores := p.PluralityScores()
	thr := float64(MajorityThreshold(p.TotalVoters()))
	c := ArgMaxFloat(scores)
	if c >= 0 && scores[c] >= thr {
		return c, true
	}
	return 0, false
}

// MajorityWinner returns the candidate holding a first-preference majority, or
// -1 when none does.
func (p *Profile) MajorityWinner() int {
	if c, ok := p.HasMajorityWinner(); ok {
		return c
	}
	return -1
}

// RunoffRound records the continuing candidates and their tallies at one round
// of an elimination rule, together with the candidate eliminated that round
// (-1 when the round is decisive).
type RunoffRound struct {
	Tallies    []float64
	Continuing []bool
	Eliminated int
}

// tallyTop counts, for each continuing candidate, the weight of ballots whose
// most preferred continuing candidate is that candidate. Exhausted ballots
// (listing no continuing candidate) are ignored.
func (p *Profile) tallyTop(continuing []bool) []float64 {
	scores := make([]float64, p.Candidates)
	for bi, b := range p.Ballots {
		if c := b.TopAmong(continuing); c >= 0 {
			scores[c] += float64(p.Weight(bi))
		}
	}
	return scores
}

// InstantRunoffRounds runs instant-runoff voting (IRV) and returns the ordered
// rounds together with the winner. Each round eliminates the continuing
// candidate with the fewest active first preferences (ties toward the lower
// index) until a candidate holds a majority of the active vote or one remains.
func (p *Profile) InstantRunoffRounds() ([]RunoffRound, int) {
	continuing := make([]bool, p.Candidates)
	remaining := p.Candidates
	for i := range continuing {
		continuing[i] = true
	}
	var rounds []RunoffRound
	for {
		tallies := p.tallyTop(continuing)
		active := 0.0
		for c := 0; c < p.Candidates; c++ {
			if continuing[c] {
				active += tallies[c]
			}
		}
		// Check for a majority of the active vote or a last survivor.
		leader, leadScore := -1, -1.0
		for c := 0; c < p.Candidates; c++ {
			if continuing[c] && tallies[c] > leadScore {
				leader, leadScore = c, tallies[c]
			}
		}
		if remaining == 1 || leadScore*2 > active {
			snap := make([]bool, p.Candidates)
			copy(snap, continuing)
			rounds = append(rounds, RunoffRound{Tallies: tallies, Continuing: snap, Eliminated: -1})
			return rounds, leader
		}
		// Eliminate the weakest continuing candidate.
		worst := -1
		for c := 0; c < p.Candidates; c++ {
			if !continuing[c] {
				continue
			}
			if worst == -1 || tallies[c] < tallies[worst] {
				worst = c
			}
		}
		snap := make([]bool, p.Candidates)
		copy(snap, continuing)
		rounds = append(rounds, RunoffRound{Tallies: tallies, Continuing: snap, Eliminated: worst})
		continuing[worst] = false
		remaining--
	}
}

// InstantRunoffWinner returns the IRV winner of the profile.
func (p *Profile) InstantRunoffWinner() int {
	_, w := p.InstantRunoffRounds()
	return w
}

// CoombsWinner returns the winner of the Coombs rule: repeatedly eliminate the
// continuing candidate with the most active last-place votes until a candidate
// holds a majority of the active first-preference vote or one remains.
func (p *Profile) CoombsWinner() int {
	continuing := make([]bool, p.Candidates)
	remaining := p.Candidates
	for i := range continuing {
		continuing[i] = true
	}
	for {
		tallies := p.tallyTop(continuing)
		active := 0.0
		leader, leadScore := -1, -1.0
		for c := 0; c < p.Candidates; c++ {
			if continuing[c] {
				active += tallies[c]
				if tallies[c] > leadScore {
					leader, leadScore = c, tallies[c]
				}
			}
		}
		if remaining == 1 || leadScore*2 > active {
			return leader
		}
		// Count active last-place votes.
		last := make([]float64, p.Candidates)
		for bi, b := range p.Ballots {
			if c := b.BottomAmong(continuing); c >= 0 {
				last[c] += float64(p.Weight(bi))
			}
		}
		worst := -1
		for c := 0; c < p.Candidates; c++ {
			if !continuing[c] {
				continue
			}
			if worst == -1 || last[c] > last[worst] {
				worst = c
			}
		}
		continuing[worst] = false
		remaining--
	}
}

// BucklinScores returns each candidate's Bucklin tally at the given depth: the
// weight of ballots ranking the candidate within their top (depth+1) positions.
func (p *Profile) BucklinScores(depth int) []float64 {
	s := make([]float64, p.Candidates)
	for bi, b := range p.Ballots {
		w := float64(p.Weight(bi))
		limit := depth + 1
		if limit > len(b) {
			limit = len(b)
		}
		for k := 0; k < limit; k++ {
			s[b[k]] += w
		}
	}
	return s
}

// BucklinWinner returns the winner of the Bucklin rule: deepen the count one
// preference level at a time until at least one candidate's tally reaches a
// majority, then return the highest-tallying such candidate.
func (p *Profile) BucklinWinner() int {
	thr := float64(MajorityThreshold(p.TotalVoters()))
	for depth := 0; depth < p.Candidates; depth++ {
		s := p.BucklinScores(depth)
		best, bestScore := -1, -1.0
		anyMajority := false
		for c := 0; c < p.Candidates; c++ {
			if s[c] >= thr {
				anyMajority = true
			}
			if s[c] > bestScore {
				best, bestScore = c, s[c]
			}
		}
		if anyMajority {
			return best
		}
	}
	return ArgMaxFloat(p.BucklinScores(p.Candidates - 1))
}

// ContingentVoteWinner returns the winner of the contingent vote: if a candidate
// holds a first-preference majority they win, otherwise all but the top two
// first-preference candidates are eliminated and the ballots are recounted among
// those two.
func (p *Profile) ContingentVoteWinner() int {
	if c, ok := p.HasMajorityWinner(); ok {
		return c
	}
	return p.topTwoRunoff(2)
}

// SupplementaryVoteWinner returns the winner of the supplementary vote, the
// contingent vote restricted to each ballot's first two preferences.
func (p *Profile) SupplementaryVoteWinner() int {
	if c, ok := p.HasMajorityWinner(); ok {
		return c
	}
	// Keep only the top two candidates and count using first two preferences.
	first := p.PluralityScores()
	keep := topKIndices(first, 2)
	continuing := make([]bool, p.Candidates)
	for _, c := range keep {
		continuing[c] = true
	}
	scores := make([]float64, p.Candidates)
	for bi, b := range p.Ballots {
		w := float64(p.Weight(bi))
		limit := len(b)
		if limit > 2 {
			limit = 2
		}
		for k := 0; k < limit; k++ {
			if continuing[b[k]] {
				scores[b[k]] += w
				break
			}
		}
	}
	best := -1
	for _, c := range keep {
		if best == -1 || scores[c] > scores[best] {
			best = c
		}
	}
	return best
}

// TwoRoundWinner returns the winner of the two-round (runoff) system: a
// first-preference majority wins outright, otherwise the top two first-round
// candidates meet in a pairwise runoff decided by the full ballots.
func (p *Profile) TwoRoundWinner() int {
	if c, ok := p.HasMajorityWinner(); ok {
		return c
	}
	return p.topTwoRunoff(2)
}

// topTwoRunoff keeps the k highest first-preference candidates and returns the
// one preferred by the most ballots when the field is reduced to them.
func (p *Profile) topTwoRunoff(k int) int {
	first := p.PluralityScores()
	keep := topKIndices(first, k)
	continuing := make([]bool, p.Candidates)
	for _, c := range keep {
		continuing[c] = true
	}
	scores := p.tallyTop(continuing)
	best := -1
	for _, c := range keep {
		if best == -1 || scores[c] > scores[best] {
			best = c
		}
	}
	return best
}

// topKIndices returns the indices of the k largest scores, breaking ties toward
// the lower index.
func topKIndices(scores []float64, k int) []int {
	order := RankingFromScores(scores)
	if k > len(order) {
		k = len(order)
	}
	return order[:k]
}

// HareQuota returns the Hare quota total/seats as a float.
func HareQuota(total, seats int) float64 {
	if seats <= 0 {
		return math.Inf(1)
	}
	return float64(total) / float64(seats)
}

// DroopQuota returns the integer Droop quota floor(total/(seats+1))+1 used by the
// single transferable vote.
func DroopQuota(total, seats int) int {
	if seats < 0 {
		return 0
	}
	return int(math.Floor(float64(total)/float64(seats+1))) + 1
}

// HagenbachBischoffQuota returns the exact Hagenbach–Bischoff quota
// total/(seats+1).
func HagenbachBischoffQuota(total, seats int) float64 {
	return float64(total) / float64(seats+1)
}

// STV runs the single transferable vote with the Droop quota and fractional
// (Gregory) surplus transfers, filling the requested number of seats. It returns
// the elected candidates in the order they reached quota or survived. Ballot
// weights start at 1 and are scaled by surplus/votes when a candidate is
// elected; the lowest candidate is eliminated when no one reaches quota.
func (p *Profile) STV(seats int) []int {
	if seats <= 0 {
		return nil
	}
	if seats > p.Candidates {
		seats = p.Candidates
	}
	quota := float64(DroopQuota(p.TotalVoters(), seats))

	// Per-ballot fractional weight.
	weight := make([]float64, len(p.Ballots))
	for i := range weight {
		weight[i] = float64(p.Weight(i))
	}
	continuing := make([]bool, p.Candidates)
	elected := make([]bool, p.Candidates)
	for i := range continuing {
		continuing[i] = true
	}
	remaining := p.Candidates
	var winners []int

	tally := func() []float64 {
		s := make([]float64, p.Candidates)
		for bi, b := range p.Ballots {
			if c := b.TopAmong(continuing); c >= 0 {
				s[c] += weight[bi]
			}
		}
		return s
	}

	for len(winners) < seats {
		if remaining == seats-len(winners) {
			// Every remaining candidate is elected to fill the seats.
			for c := 0; c < p.Candidates; c++ {
				if continuing[c] {
					continuing[c] = false
					elected[c] = true
					winners = append(winners, c)
				}
			}
			break
		}
		scores := tally()
		// Find the strongest continuing candidate.
		leader := -1
		for c := 0; c < p.Candidates; c++ {
			if continuing[c] && (leader == -1 || scores[c] > scores[leader]) {
				leader = c
			}
		}
		if leader >= 0 && scores[leader] >= quota-1e-9 {
			// Elect the leader and transfer the surplus.
			surplus := scores[leader] - quota
			factor := 0.0
			if scores[leader] > 0 {
				factor = surplus / scores[leader]
			}
			for bi, b := range p.Ballots {
				if b.TopAmong(continuing) == leader {
					weight[bi] *= factor
				}
			}
			continuing[leader] = false
			elected[leader] = true
			winners = append(winners, leader)
			remaining--
			continue
		}
		// No one reached quota: eliminate the weakest candidate.
		worst := -1
		for c := 0; c < p.Candidates; c++ {
			if continuing[c] && (worst == -1 || scores[c] < scores[worst]) {
				worst = c
			}
		}
		if worst == -1 {
			break
		}
		continuing[worst] = false
		remaining--
	}
	return winners
}
