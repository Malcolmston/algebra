package socialchoice

import (
	"fmt"
	"math"
	"reflect"
	"testing"
)

// tennessee returns the standard Tennessee-capital ranked-ballot example with
// candidates 0=Memphis, 1=Nashville, 2=Chattanooga, 3=Knoxville.
func tennessee(t *testing.T) *Profile {
	t.Helper()
	p, err := NewProfile(4,
		[]Ballot{{0, 1, 2, 3}, {1, 2, 3, 0}, {2, 3, 1, 0}, {3, 2, 1, 0}},
		[]int{42, 26, 15, 17},
	)
	if err != nil {
		t.Fatalf("NewProfile: %v", err)
	}
	return p
}

func approx(a, b, tol float64) bool { return math.Abs(a-b) <= tol }

func TestProfileBasics(t *testing.T) {
	p := tennessee(t)
	if got := p.TotalVoters(); got != 100 {
		t.Fatalf("TotalVoters = %d, want 100", got)
	}
	if got := p.NumBallots(); got != 4 {
		t.Fatalf("NumBallots = %d, want 4", got)
	}
	if !p.Ballots[0].Prefers(0, 3) || p.Ballots[0].Prefers(3, 0) {
		t.Fatal("ballot 0 should prefer 0 to 3")
	}
	if got, ok := p.Ballots[1].Top(); !ok || got != 1 {
		t.Fatalf("ballot 1 top = %d,%v", got, ok)
	}
	if got, ok := p.Ballots[1].Bottom(); !ok || got != 0 {
		t.Fatalf("ballot 1 bottom = %d,%v", got, ok)
	}
}

func TestProfileValidate(t *testing.T) {
	if _, err := NewProfile(0, nil, nil); err == nil {
		t.Fatal("want error for zero candidates")
	}
	if _, err := NewProfile(3, []Ballot{{0, 0}}, nil); err == nil {
		t.Fatal("want error for duplicate candidate on a ballot")
	}
	if _, err := NewProfile(3, []Ballot{{3}}, nil); err == nil {
		t.Fatal("want error for out-of-range candidate")
	}
	if _, err := NewProfile(3, []Ballot{{0}}, []int{-1}); err == nil {
		t.Fatal("want error for negative count")
	}
}

func TestSingleWinnerRules(t *testing.T) {
	p := tennessee(t)
	cases := []struct {
		name string
		got  int
		want int
	}{
		{"plurality", p.PluralityWinner(), 0},
		{"antiplurality", p.AntiPluralityWinner(), 1},
		{"borda", p.BordaWinner(), 1},
		{"copeland", p.CopelandWinner(), 1},
		{"minimax", p.MiniMaxWinner(), 1},
		{"schulze", p.SchulzeWinner(), 1},
		{"rankedpairs", p.RankedPairsWinner(), 1},
		{"kemeny", p.KemenyWinner(), 1},
		{"black", p.BlackWinner(), 1},
		{"nanson", p.NansonWinner(), 1},
		{"baldwin", p.BaldwinWinner(), 1},
		{"irv", p.InstantRunoffWinner(), 3},
		{"coombs", p.CoombsWinner(), 1},
		{"bucklin", p.BucklinWinner(), 1},
		{"tworound", p.TwoRoundWinner(), 1},
		{"contingent", p.ContingentVoteWinner(), 1},
		{"supplementary", p.SupplementaryVoteWinner(), 0},
	}
	for _, c := range cases {
		if c.got != c.want {
			t.Errorf("%s winner = %d, want %d", c.name, c.got, c.want)
		}
	}
}

func TestPairwiseAndScores(t *testing.T) {
	p := tennessee(t)
	m := p.Pairwise()
	if m[1][0] != 58 || m[0][1] != 42 {
		t.Fatalf("pairwise 0/1 = %d/%d, want 58/42", m[1][0], m[0][1])
	}
	if !m.Beats(1, 0) || !m.Beats(1, 2) || !m.Beats(1, 3) {
		t.Fatal("Nashville should beat all others")
	}
	if cw, ok := m.CondorcetWinner(); !ok || cw != 1 {
		t.Fatalf("CondorcetWinner = %d,%v want 1,true", cw, ok)
	}
	if cl, ok := m.CondorcetLoser(); !ok || cl != 0 {
		t.Fatalf("CondorcetLoser = %d,%v want 0,true", cl, ok)
	}
	wantBorda := []float64{126, 194, 173, 107}
	if got := p.BordaScores(); !reflect.DeepEqual(got, wantBorda) {
		t.Fatalf("BordaScores = %v, want %v", got, wantBorda)
	}
	wantCope := []float64{0, 3, 2, 1}
	if got := p.CopelandScores(); !reflect.DeepEqual(got, wantCope) {
		t.Fatalf("CopelandScores = %v, want %v", got, wantCope)
	}
	wantMM := []int{16, 0, 36, 66}
	if got := m.MiniMaxScores(); !reflect.DeepEqual(got, wantMM) {
		t.Fatalf("MiniMaxScores = %v, want %v", got, wantMM)
	}
}

func TestSetsTennessee(t *testing.T) {
	p := tennessee(t)
	if got := p.SmithSet(); !reflect.DeepEqual(got, []int{1}) {
		t.Fatalf("SmithSet = %v, want [1]", got)
	}
	if got := p.SchwartzSet(); !reflect.DeepEqual(got, []int{1}) {
		t.Fatalf("SchwartzSet = %v, want [1]", got)
	}
	if got := p.UncoveredSet(); !reflect.DeepEqual(got, []int{1}) {
		t.Fatalf("UncoveredSet = %v, want [1]", got)
	}
	if p.HasCondorcetCycle() {
		t.Fatal("Tennessee has a Condorcet winner, so no cycle")
	}
	if !p.HasCondorcetWinner() {
		t.Fatal("expected a Condorcet winner")
	}
}

func TestCondorcetParadox(t *testing.T) {
	p, err := NewProfile(3, []Ballot{{0, 1, 2}, {1, 2, 0}, {2, 0, 1}}, nil)
	if err != nil {
		t.Fatal(err)
	}
	if p.HasCondorcetWinner() {
		t.Fatal("cyclic profile should have no Condorcet winner")
	}
	if !p.CondorcetParadox() {
		t.Fatal("expected Condorcet paradox")
	}
	if !p.HasCondorcetCycle() {
		t.Fatal("expected a Condorcet cycle")
	}
	if got := p.SmithSet(); !reflect.DeepEqual(got, []int{0, 1, 2}) {
		t.Fatalf("SmithSet = %v, want all three", got)
	}
}

func TestIRVRounds(t *testing.T) {
	p := tennessee(t)
	rounds, w := p.InstantRunoffRounds()
	if w != 3 {
		t.Fatalf("IRV winner = %d, want 3", w)
	}
	// First eliminated is Chattanooga (fewest first preferences).
	if rounds[0].Eliminated != 2 {
		t.Fatalf("first eliminated = %d, want 2", rounds[0].Eliminated)
	}
}

func TestSTV(t *testing.T) {
	p, err := NewProfile(3,
		[]Ballot{{0, 1, 2}, {0, 2, 1}, {1, 2, 0}, {2, 1, 0}},
		[]int{40, 10, 30, 20},
	)
	if err != nil {
		t.Fatal(err)
	}
	if q := DroopQuota(100, 2); q != 34 {
		t.Fatalf("DroopQuota = %d, want 34", q)
	}
	got := p.STV(2)
	if !reflect.DeepEqual(got, []int{0, 1}) {
		t.Fatalf("STV winners = %v, want [0 1]", got)
	}
}

func TestApproval(t *testing.T) {
	a := ApprovalProfile{Candidates: 3, Ballots: [][]bool{
		{true, true, false},
		{false, true, true},
		{false, true, false},
	}}
	if got := a.Scores(); !reflect.DeepEqual(got, []float64{1, 3, 1}) {
		t.Fatalf("approval scores = %v", got)
	}
	if a.Winner() != 1 {
		t.Fatalf("approval winner = %d, want 1", a.Winner())
	}
	if a.TotalApprovals() != 5 {
		t.Fatalf("total approvals = %d, want 5", a.TotalApprovals())
	}
}

func TestScoreAndSTAR(t *testing.T) {
	s := ScoreProfile{Candidates: 3, Max: 5, Ballots: [][]float64{
		{5, 3, 0},
		{4, 5, 1},
		{0, 4, 2},
	}}
	if got := s.Totals(); !reflect.DeepEqual(got, []float64{9, 12, 3}) {
		t.Fatalf("score totals = %v", got)
	}
	if s.Winner() != 1 {
		t.Fatalf("score winner = %d, want 1", s.Winner())
	}
	if got := s.STARWinner(); got != 1 {
		t.Fatalf("STAR winner = %d, want 1", got)
	}
}

func TestMajorityJudgment(t *testing.T) {
	j := JudgmentProfile{Candidates: 2, NumGrades: 3, Grades: [][]int{
		{2, 1}, {2, 0}, {1, 2},
	}}
	if got := j.MedianGrades(); !reflect.DeepEqual(got, []int{2, 1}) {
		t.Fatalf("median grades = %v, want [2 1]", got)
	}
	if j.Winner() != 0 {
		t.Fatalf("MJ winner = %d, want 0", j.Winner())
	}
}

func TestCumulative(t *testing.T) {
	c := CumulativeProfile{Candidates: 3, Ballots: [][]float64{{3, 0, 0}, {1, 2, 0}}}
	if got := c.Scores(); !reflect.DeepEqual(got, []float64{4, 2, 0}) {
		t.Fatalf("cumulative scores = %v", got)
	}
	if c.Winner() != 0 {
		t.Fatalf("cumulative winner = %d, want 0", c.Winner())
	}
}

func TestApportionment(t *testing.T) {
	votes := []int{100000, 80000, 30000, 20000}
	seats := 8
	cases := []struct {
		name string
		fn   func([]int, int) ([]int, error)
		want []int
	}{
		{"dhondt", DHondt, []int{4, 3, 1, 0}},
		{"saintelague", SainteLague, []int{3, 3, 1, 1}},
		{"hamilton", Hamilton, []int{3, 3, 1, 1}},
		{"huntinghill", HuntingtonHill, []int{3, 3, 1, 1}},
		{"adams", Adams, []int{3, 3, 1, 1}},
	}
	for _, c := range cases {
		got, err := c.fn(votes, seats)
		if err != nil {
			t.Fatalf("%s: %v", c.name, err)
		}
		if !reflect.DeepEqual(got, c.want) {
			t.Errorf("%s = %v, want %v", c.name, got, c.want)
		}
		if TotalSeats(got) != seats {
			t.Errorf("%s allocated %d seats, want %d", c.name, TotalSeats(got), seats)
		}
	}
}

func TestQuotaRule(t *testing.T) {
	votes := []int{100000, 80000, 30000, 20000}
	seats := 8
	dh, _ := DHondt(votes, seats)
	if !SatisfiesQuota(dh, votes, seats) {
		t.Fatal("D'Hondt allocation should satisfy the quota rule here")
	}
	ham, _ := Hamilton(votes, seats)
	if !SatisfiesQuota(ham, votes, seats) {
		t.Fatal("Hamilton always satisfies the quota rule")
	}
	if lq := LowerQuotas(votes, seats); !reflect.DeepEqual(lq, []int{3, 2, 1, 0}) {
		t.Fatalf("LowerQuotas = %v", lq)
	}
}

func TestAlabamaParadox(t *testing.T) {
	if !HasAlabamaParadox([]int{6, 6, 2}, 10) {
		t.Fatal("expected Alabama paradox for {6,6,2} at 10 -> 11 seats")
	}
	if HasAlabamaParadox([]int{100, 100, 100}, 3) {
		t.Fatal("did not expect Alabama paradox for equal parties")
	}
	if !IsHouseMonotone([]int{6, 6, 2}, 12, func(s int) float64 { return float64(s + 1) }) {
		t.Fatal("D'Hondt should be house-monotone")
	}
}

func TestDisproportionality(t *testing.T) {
	votes := []int{100000, 80000, 30000, 20000}
	seats := []int{4, 3, 1, 0}
	g := GallagherIndex(votes, seats)
	if !approx(g, 7.93, 0.1) {
		t.Fatalf("Gallagher index = %.3f, want ~7.93", g)
	}
	enp := EffectiveNumberOfParties(votes)
	if !approx(enp, 2.99, 0.05) {
		t.Fatalf("ENP = %.3f, want ~2.99", enp)
	}
	if lh := LoosemoreHanbyIndex(votes, seats); lh < 0 {
		t.Fatalf("Loosemore-Hanby negative: %v", lh)
	}
}

func TestPositionalGeneric(t *testing.T) {
	p := tennessee(t)
	// Borda weights (3,2,1,0) must reproduce the Borda scores.
	if got := p.PositionalScores([]float64{3, 2, 1, 0}); !reflect.DeepEqual(got, p.BordaScores()) {
		t.Fatalf("positional Borda mismatch: %v vs %v", got, p.BordaScores())
	}
	if w := p.DowdallWinner(); w != 0 && w != 1 {
		t.Fatalf("unexpected Dowdall winner %d", w)
	}
}

func TestDerivedProfiles(t *testing.T) {
	p := tennessee(t)
	a := ApprovalFromRanked(p, 1)
	if got := a.Scores(); !reflect.DeepEqual(got, p.PluralityScores()) {
		t.Fatalf("top-1 approval should equal plurality: %v vs %v", got, p.PluralityScores())
	}
	r := RangeFromRanked(p)
	if r.Winner() != p.BordaWinner() {
		t.Fatalf("range-from-ranked winner %d, Borda %d", r.Winner(), p.BordaWinner())
	}
}

func ExampleProfile() {
	p, _ := NewProfile(4,
		[]Ballot{{0, 1, 2, 3}, {1, 2, 3, 0}, {2, 3, 1, 0}, {3, 2, 1, 0}},
		[]int{42, 26, 15, 17},
	)
	cw, _ := p.CondorcetWinner()
	fmt.Println("plurality:", p.PluralityWinner())
	fmt.Println("condorcet:", cw)
	fmt.Println("irv:", p.InstantRunoffWinner())
	// Output:
	// plurality: 0
	// condorcet: 1
	// irv: 3
}
