package socialchoice

import "sort"

// ApprovalProfile is a set of approval ballots over Candidates candidates; each
// ballot approves the candidates whose entry is true.
type ApprovalProfile struct {
	Candidates int
	Ballots    [][]bool
}

// Scores returns each candidate's approval tally, the number of ballots
// approving it.
func (a ApprovalProfile) Scores() []float64 {
	s := make([]float64, a.Candidates)
	for _, b := range a.Ballots {
		for c := 0; c < a.Candidates && c < len(b); c++ {
			if b[c] {
				s[c]++
			}
		}
	}
	return s
}

// Winner returns the most approved candidate.
func (a ApprovalProfile) Winner() int { return ArgMaxFloat(a.Scores()) }

// Ranking returns candidates ordered by descending approval.
func (a ApprovalProfile) Ranking() []int { return RankingFromScores(a.Scores()) }

// ScoreProfile is a set of range/score ballots; Ballots[v][c] is the score voter
// v assigns candidate c, bounded above by Max.
type ScoreProfile struct {
	Candidates int
	Ballots    [][]float64
	Max        float64
}

// Totals returns the sum of scores each candidate receives.
func (s ScoreProfile) Totals() []float64 {
	t := make([]float64, s.Candidates)
	for _, b := range s.Ballots {
		for c := 0; c < s.Candidates && c < len(b); c++ {
			t[c] += b[c]
		}
	}
	return t
}

// Averages returns each candidate's mean score over all ballots.
func (s ScoreProfile) Averages() []float64 {
	t := s.Totals()
	n := float64(len(s.Ballots))
	if n == 0 {
		return t
	}
	for c := range t {
		t[c] /= n
	}
	return t
}

// Winner returns the highest total-score candidate (range/score voting).
func (s ScoreProfile) Winner() int { return ArgMaxFloat(s.Totals()) }

// Ranking returns candidates ordered by descending total score.
func (s ScoreProfile) Ranking() []int { return RankingFromScores(s.Totals()) }

// STARWinner returns the STAR (Score Then Automatic Runoff) winner: the two
// highest-scoring candidates advance to an automatic runoff won by whichever is
// preferred (scored higher) on more ballots, ties broken by total score and then
// the lower index.
func (s ScoreProfile) STARWinner() int {
	if s.Candidates == 0 {
		return -1
	}
	totals := s.Totals()
	if s.Candidates == 1 {
		return 0
	}
	top := topKIndices(totals, 2)
	a, b := top[0], top[1]
	var prefA, prefB float64
	for _, bal := range s.Ballots {
		sa, sb := 0.0, 0.0
		if a < len(bal) {
			sa = bal[a]
		}
		if b < len(bal) {
			sb = bal[b]
		}
		switch {
		case sa > sb:
			prefA++
		case sb > sa:
			prefB++
		}
	}
	switch {
	case prefA > prefB:
		return a
	case prefB > prefA:
		return b
	case totals[a] >= totals[b]:
		return a
	default:
		return b
	}
}

// CumulativeProfile is a set of cumulative-voting ballots; Ballots[v][c] is the
// number of points voter v allocates to candidate c.
type CumulativeProfile struct {
	Candidates int
	Ballots    [][]float64
}

// Scores returns the total points each candidate accumulates.
func (c CumulativeProfile) Scores() []float64 {
	t := make([]float64, c.Candidates)
	for _, b := range c.Ballots {
		for i := 0; i < c.Candidates && i < len(b); i++ {
			t[i] += b[i]
		}
	}
	return t
}

// Winner returns the candidate with the most cumulative points.
func (c CumulativeProfile) Winner() int { return ArgMaxFloat(c.Scores()) }

// JudgmentProfile is a set of graded ballots for majority-judgment voting;
// Grades[v][c] is the grade in 0..NumGrades-1 that voter v assigns candidate c,
// with higher grades better.
type JudgmentProfile struct {
	Candidates int
	NumGrades  int
	Grades     [][]int
}

// gradeCounts returns, for each candidate, a histogram of grades.
func (j JudgmentProfile) gradeCounts() [][]int {
	counts := make([][]int, j.Candidates)
	for c := range counts {
		counts[c] = make([]int, j.NumGrades)
	}
	for _, ballot := range j.Grades {
		for c := 0; c < j.Candidates && c < len(ballot); c++ {
			g := ballot[c]
			if g >= 0 && g < j.NumGrades {
				counts[c][g]++
			}
		}
	}
	return counts
}

// lowerMedian returns the lower-median grade of a grade histogram with the given
// total, or -1 when empty.
func lowerMedian(counts []int, total int) int {
	if total == 0 {
		return -1
	}
	r := (total - 1) / 2
	cum := 0
	for g := 0; g < len(counts); g++ {
		cum += counts[g]
		if cum > r {
			return g
		}
	}
	return len(counts) - 1
}

// MedianGrades returns each candidate's lower-median grade.
func (j JudgmentProfile) MedianGrades() []int {
	counts := j.gradeCounts()
	med := make([]int, j.Candidates)
	for c := range counts {
		total := 0
		for _, x := range counts[c] {
			total += x
		}
		med[c] = lowerMedian(counts[c], total)
	}
	return med
}

// Winner returns the majority-judgment winner: the highest median grade wins,
// and ties are broken by repeatedly removing one median grade from each tied
// candidate and recomparing medians. A residual tie resolves toward the lower
// index.
func (j JudgmentProfile) Winner() int {
	if j.Candidates == 0 {
		return -1
	}
	counts := j.gradeCounts()
	totals := make([]int, j.Candidates)
	for c := range counts {
		for _, x := range counts[c] {
			totals[c] += x
		}
	}
	// Candidates still in contention.
	live := make([]bool, j.Candidates)
	for c := range live {
		live[c] = true
	}
	for {
		// Compute current medians for live candidates.
		best, bestMed := -1, -1
		var tied []int
		for c := 0; c < j.Candidates; c++ {
			if !live[c] {
				continue
			}
			med := lowerMedian(counts[c], totals[c])
			if med > bestMed {
				bestMed, best = med, c
				tied = []int{c}
			} else if med == bestMed {
				tied = append(tied, c)
			}
		}
		if len(tied) <= 1 {
			return best
		}
		// Drop non-tied candidates from contention.
		for c := range live {
			live[c] = false
		}
		progressed := false
		for _, c := range tied {
			live[c] = true
			med := lowerMedian(counts[c], totals[c])
			if med >= 0 && counts[c][med] > 0 && totals[c] > 0 {
				counts[c][med]--
				totals[c]--
				progressed = true
			}
		}
		if !progressed {
			// All tied candidates exhausted; resolve toward the lower index.
			sort.Ints(tied)
			return tied[0]
		}
	}
}
