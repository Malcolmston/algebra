package socialchoice

// PluralityScores returns each candidate's first-preference tally: the total
// weight of ballots ranking that candidate top.
func (p *Profile) PluralityScores() []float64 {
	s := make([]float64, p.Candidates)
	for bi, b := range p.Ballots {
		if c, ok := b.Top(); ok {
			s[c] += float64(p.Weight(bi))
		}
	}
	return s
}

// PluralityWinner returns the candidate with the most first preferences.
func (p *Profile) PluralityWinner() int { return ArgMaxFloat(p.PluralityScores()) }

// PluralityRanking returns candidates ordered by descending first-preference
// tally.
func (p *Profile) PluralityRanking() []int { return RankingFromScores(p.PluralityScores()) }

// FirstPreferenceCounts is an alias for PluralityScores, returning the top-choice
// tally of every candidate.
func (p *Profile) FirstPreferenceCounts() []float64 { return p.PluralityScores() }

// LastPreferenceCounts returns each candidate's tally of last-place (bottom
// listed) rankings, the basis of the Coombs rule.
func (p *Profile) LastPreferenceCounts() []float64 {
	s := make([]float64, p.Candidates)
	for bi, b := range p.Ballots {
		if c, ok := b.Bottom(); ok {
			s[c] += float64(p.Weight(bi))
		}
	}
	return s
}

// AntiPluralityScores returns each candidate's anti-plurality (veto) score: the
// total weight of ballots that do not rank the candidate last. Each ballot
// vetoes only its bottom listed candidate.
func (p *Profile) AntiPluralityScores() []float64 {
	total := float64(p.TotalVoters())
	veto := p.LastPreferenceCounts()
	s := make([]float64, p.Candidates)
	for c := range s {
		s[c] = total - veto[c]
	}
	return s
}

// AntiPluralityWinner returns the candidate vetoed by the fewest ballots.
func (p *Profile) AntiPluralityWinner() int { return ArgMaxFloat(p.AntiPluralityScores()) }

// PositionalScores returns each candidate's score under an arbitrary positional
// rule: a candidate ranked at position k on a ballot earns weights[k] points,
// scaled by the ballot weight. Candidates not listed on a ballot earn the score
// of the final position, weights[len(weights)-1]. weights must have length equal
// to the candidate count.
func (p *Profile) PositionalScores(weights []float64) []float64 {
	s := make([]float64, p.Candidates)
	if len(weights) == 0 {
		return s
	}
	last := weights[len(weights)-1]
	for bi, b := range p.Ballots {
		w := float64(p.Weight(bi))
		listed := make([]bool, p.Candidates)
		for pos, c := range b {
			pw := last
			if pos < len(weights) {
				pw = weights[pos]
			}
			s[c] += w * pw
			listed[c] = true
		}
		for c := 0; c < p.Candidates; c++ {
			if !listed[c] {
				s[c] += w * last
			}
		}
	}
	return s
}

// PositionalWinner returns the winner under the positional rule with the given
// position weights.
func (p *Profile) PositionalWinner(weights []float64) int {
	return ArgMaxFloat(p.PositionalScores(weights))
}

// BordaScores returns each candidate's Borda score, the total weight of voters
// preferring the candidate to each other candidate summed over all opponents.
// For complete strict ballots this equals the usual (n-1, n-2, …, 0) tally, and
// the pairwise formulation extends it consistently to truncated ballots.
func (p *Profile) BordaScores() []float64 {
	m := p.Pairwise()
	n := p.Candidates
	s := make([]float64, n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if i != j {
				s[i] += float64(m[i][j])
			}
		}
	}
	return s
}

// BordaWinner returns the highest Borda-scoring candidate.
func (p *Profile) BordaWinner() int { return ArgMaxFloat(p.BordaScores()) }

// BordaRanking returns candidates ordered by descending Borda score.
func (p *Profile) BordaRanking() []int { return RankingFromScores(p.BordaScores()) }

// DowdallScores returns each candidate's Dowdall (Nauru) score, awarding
// 1/(k+1) points for a rank at zero-based position k.
func (p *Profile) DowdallScores() []float64 {
	s := make([]float64, p.Candidates)
	for bi, b := range p.Ballots {
		w := float64(p.Weight(bi))
		for pos, c := range b {
			s[c] += w / float64(pos+1)
		}
	}
	return s
}

// DowdallWinner returns the highest Dowdall-scoring candidate.
func (p *Profile) DowdallWinner() int { return ArgMaxFloat(p.DowdallScores()) }

// bordaAmong returns the Borda score of each continuing candidate within the
// sub-election restricted to the continuing candidates, using the fixed pairwise
// matrix (eliminating a candidate only removes its pairwise comparisons).
func bordaAmong(m PairwiseMatrix, continuing []bool) []float64 {
	n := len(continuing)
	s := make([]float64, n)
	for i := 0; i < n; i++ {
		if !continuing[i] {
			continue
		}
		for j := 0; j < n; j++ {
			if i != j && continuing[j] {
				s[i] += float64(m[i][j])
			}
		}
	}
	return s
}

// NansonWinner returns the winner of the Nanson rule: repeatedly eliminate every
// continuing candidate whose Borda score is strictly below the average Borda
// score until one candidate remains (ties resolved toward the lower index).
func (p *Profile) NansonWinner() int {
	m := p.Pairwise()
	continuing := make([]bool, p.Candidates)
	remaining := p.Candidates
	for i := range continuing {
		continuing[i] = true
	}
	for remaining > 1 {
		scores := bordaAmong(m, continuing)
		var sum float64
		for i := range scores {
			if continuing[i] {
				sum += scores[i]
			}
		}
		avg := sum / float64(remaining)
		eliminated := 0
		toDrop := make([]bool, p.Candidates)
		for i := range scores {
			if continuing[i] && scores[i] < avg {
				toDrop[i] = true
				eliminated++
			}
		}
		if eliminated == 0 || eliminated == remaining {
			break
		}
		for i := range toDrop {
			if toDrop[i] {
				continuing[i] = false
				remaining--
			}
		}
	}
	for i, ok := range continuing {
		if ok {
			return i
		}
	}
	return -1
}

// BaldwinWinner returns the winner of the Baldwin rule: repeatedly eliminate the
// single continuing candidate with the lowest Borda score until one remains
// (ties resolved toward the lower index).
func (p *Profile) BaldwinWinner() int {
	m := p.Pairwise()
	continuing := make([]bool, p.Candidates)
	remaining := p.Candidates
	for i := range continuing {
		continuing[i] = true
	}
	for remaining > 1 {
		scores := bordaAmong(m, continuing)
		worst := -1
		for i := range scores {
			if !continuing[i] {
				continue
			}
			if worst == -1 || scores[i] < scores[worst] {
				worst = i
			}
		}
		continuing[worst] = false
		remaining--
	}
	for i, ok := range continuing {
		if ok {
			return i
		}
	}
	return -1
}
