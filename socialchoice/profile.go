package socialchoice

import (
	"errors"
	"fmt"
	"sort"
)

// Ballot is a strict preference order over candidate indices, most preferred
// first. A ballot may be truncated: any candidate not listed is considered
// ranked below every listed candidate and tied with the other omitted ones.
type Ballot []int

// Profile is a weighted collection of ranked ballots over a fixed number of
// candidates identified by the indices 0..Candidates-1. Counts[i] is the number
// of voters casting Ballots[i]; a nil or short Counts slice treats the missing
// weights as 1.
type Profile struct {
	Candidates int
	Ballots    []Ballot
	Counts     []int
}

// NewProfile builds a Profile over the given number of candidates and validates
// it. counts may be nil, in which case every ballot has weight 1.
func NewProfile(candidates int, ballots []Ballot, counts []int) (*Profile, error) {
	p := &Profile{Candidates: candidates, Ballots: ballots, Counts: counts}
	if err := p.Validate(); err != nil {
		return nil, err
	}
	return p, nil
}

// Validate reports whether the profile is well formed: a positive candidate
// count, in-range and duplicate-free ballot entries, and non-negative counts of
// a length not exceeding the number of ballots.
func (p *Profile) Validate() error {
	if p.Candidates <= 0 {
		return errors.New("socialchoice: Candidates must be positive")
	}
	if len(p.Counts) > len(p.Ballots) {
		return errors.New("socialchoice: more counts than ballots")
	}
	for bi, b := range p.Ballots {
		seen := make([]bool, p.Candidates)
		for _, c := range b {
			if c < 0 || c >= p.Candidates {
				return fmt.Errorf("socialchoice: ballot %d references candidate %d out of range", bi, c)
			}
			if seen[c] {
				return fmt.Errorf("socialchoice: ballot %d lists candidate %d twice", bi, c)
			}
			seen[c] = true
		}
	}
	for i, w := range p.Counts {
		if w < 0 {
			return fmt.Errorf("socialchoice: negative count %d at ballot %d", w, i)
		}
	}
	return nil
}

// Weight returns the number of voters casting ballot i, defaulting to 1 when
// Counts does not cover that ballot.
func (p *Profile) Weight(i int) int {
	if i < len(p.Counts) {
		return p.Counts[i]
	}
	return 1
}

// NumBallots returns the number of distinct ballots recorded in the profile.
func (p *Profile) NumBallots() int { return len(p.Ballots) }

// TotalVoters returns the sum of the ballot weights, i.e. the electorate size.
func (p *Profile) TotalVoters() int {
	total := 0
	for i := range p.Ballots {
		total += p.Weight(i)
	}
	return total
}

// AddBallot appends a ballot with the given weight and returns the updated
// profile pointer for chaining. It does not re-validate; call Validate when
// building profiles from untrusted input.
func (p *Profile) AddBallot(b Ballot, weight int) *Profile {
	p.Ballots = append(p.Ballots, b)
	// Grow Counts to keep indices aligned.
	for len(p.Counts) < len(p.Ballots)-1 {
		p.Counts = append(p.Counts, 1)
	}
	p.Counts = append(p.Counts, weight)
	return p
}

// Position returns the zero-based rank of candidate c on the ballot and whether
// c appears at all. Rank 0 is the most preferred candidate.
func (b Ballot) Position(c int) (int, bool) {
	for i, x := range b {
		if x == c {
			return i, true
		}
	}
	return 0, false
}

// Ranks reports whether candidate c appears anywhere on the ballot.
func (b Ballot) Ranks(c int) bool {
	_, ok := b.Position(c)
	return ok
}

// Prefers reports whether the ballot strictly prefers candidate a to candidate
// b. A listed candidate is preferred to any unlisted one; two unlisted
// candidates yield no preference (false in both directions).
func (b Ballot) Prefers(a, c int) bool {
	pa, oka := b.Position(a)
	pc, okc := b.Position(c)
	switch {
	case oka && okc:
		return pa < pc
	case oka:
		return true
	default:
		return false
	}
}

// Top returns the most preferred candidate on the ballot and whether the ballot
// is non-empty.
func (b Ballot) Top() (int, bool) {
	if len(b) == 0 {
		return 0, false
	}
	return b[0], true
}

// Bottom returns the least preferred listed candidate and whether the ballot is
// non-empty. Unlisted candidates, though notionally lower, are not returned.
func (b Ballot) Bottom() (int, bool) {
	if len(b) == 0 {
		return 0, false
	}
	return b[len(b)-1], true
}

// Reversed returns a new ballot with the preference order inverted.
func (b Ballot) Reversed() Ballot {
	r := make(Ballot, len(b))
	for i, x := range b {
		r[len(b)-1-i] = x
	}
	return r
}

// TopAmong returns the most preferred candidate on the ballot for which
// continuing[c] is true, or -1 if the ballot lists no continuing candidate.
func (b Ballot) TopAmong(continuing []bool) int {
	for _, c := range b {
		if c >= 0 && c < len(continuing) && continuing[c] {
			return c
		}
	}
	return -1
}

// BottomAmong returns the least preferred listed candidate for which
// continuing[c] is true, or -1 if none is listed.
func (b Ballot) BottomAmong(continuing []bool) int {
	last := -1
	for _, c := range b {
		if c >= 0 && c < len(continuing) && continuing[c] {
			last = c
		}
	}
	return last
}

// ArgMaxFloat returns the index of the maximum value, breaking ties toward the
// lowest index, and -1 for an empty slice.
func ArgMaxFloat(xs []float64) int {
	if len(xs) == 0 {
		return -1
	}
	best := 0
	for i := 1; i < len(xs); i++ {
		if xs[i] > xs[best] {
			best = i
		}
	}
	return best
}

// ArgMinFloat returns the index of the minimum value, breaking ties toward the
// lowest index, and -1 for an empty slice.
func ArgMinFloat(xs []float64) int {
	if len(xs) == 0 {
		return -1
	}
	best := 0
	for i := 1; i < len(xs); i++ {
		if xs[i] < xs[best] {
			best = i
		}
	}
	return best
}

// ArgMaxInt returns the index of the maximum value, breaking ties toward the
// lowest index, and -1 for an empty slice.
func ArgMaxInt(xs []int) int {
	if len(xs) == 0 {
		return -1
	}
	best := 0
	for i := 1; i < len(xs); i++ {
		if xs[i] > xs[best] {
			best = i
		}
	}
	return best
}

// RankingFromScores returns candidate indices ordered by descending score,
// breaking ties toward the lower index.
func RankingFromScores(scores []float64) []int {
	idx := make([]int, len(scores))
	for i := range idx {
		idx[i] = i
	}
	sort.SliceStable(idx, func(a, b int) bool {
		if scores[idx[a]] != scores[idx[b]] {
			return scores[idx[a]] > scores[idx[b]]
		}
		return idx[a] < idx[b]
	})
	return idx
}

// WinnerFromScores returns the index of the highest-scoring candidate, breaking
// ties toward the lower index.
func WinnerFromScores(scores []float64) int { return ArgMaxFloat(scores) }
