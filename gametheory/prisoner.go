package gametheory

import "errors"

// Pure-strategy indices shared by the Prisoner's-dilemma helpers and iterated
// play: Cooperate is strategy 0 and Defect is strategy 1.
const (
	// Cooperate is the cooperative move (strategy index 0).
	Cooperate = 0
	// Defect is the defecting move (strategy index 1).
	Defect = 1
)

// PrisonersDilemma holds the four payoffs of a symmetric Prisoner's dilemma from
// one player's perspective: T is the temptation payoff (defect against a
// cooperator), R the reward for mutual cooperation, P the punishment for mutual
// defection, and S the sucker's payoff (cooperate against a defector). A proper
// dilemma requires T > R > P > S.
type PrisonersDilemma struct {
	T float64
	R float64
	P float64
	S float64
}

// ErrPayoffOrder is returned when Prisoner's-dilemma payoffs violate T>R>P>S.
var ErrPayoffOrder = errors.New("gametheory: Prisoner's dilemma requires T > R > P > S")

// NewPrisonersDilemma constructs a PrisonersDilemma, returning ErrPayoffOrder
// unless the payoffs satisfy the defining ordering t > r > p > s.
func NewPrisonersDilemma(t, r, p, s float64) (PrisonersDilemma, error) {
	pd := PrisonersDilemma{T: t, R: r, P: p, S: s}
	if !pd.IsValid() {
		return PrisonersDilemma{}, ErrPayoffOrder
	}
	return pd, nil
}

// StandardPrisonersDilemma returns the canonical Prisoner's dilemma with the
// textbook payoffs T=5, R=3, P=1, S=0.
func StandardPrisonersDilemma() PrisonersDilemma {
	return PrisonersDilemma{T: 5, R: 3, P: 1, S: 0}
}

// IsValid reports whether the payoffs satisfy the Prisoner's-dilemma ordering
// T > R > P > S.
func (pd PrisonersDilemma) IsValid() bool {
	return pd.T > pd.R && pd.R > pd.P && pd.P > pd.S
}

// IsIterationStable reports whether the payoffs satisfy 2R > T + S, the extra
// condition ensuring that in the iterated game sustained mutual cooperation
// out-scores alternating exploitation.
func (pd PrisonersDilemma) IsIterationStable() bool {
	return 2*pd.R > pd.T+pd.S
}

// Game returns the two-player normal-form game for the Prisoner's dilemma, with
// row and column strategy 0 = Cooperate and 1 = Defect. The row player's payoff
// matrix is [[R, S], [T, P]] and the column player's is its transpose.
func (pd PrisonersDilemma) Game() Game {
	row := [][]float64{
		{pd.R, pd.S},
		{pd.T, pd.P},
	}
	col := [][]float64{
		{pd.R, pd.T},
		{pd.S, pd.P},
	}
	return Game{Row: row, Col: col}
}

// IPDStrategy is a deterministic strategy for the iterated Prisoner's dilemma.
// Given the caller's own past moves (myHistory) and the opponent's past moves
// (oppHistory), both equal-length slices of Cooperate/Defect ordered oldest to
// newest, it returns the move to play this round.
type IPDStrategy func(myHistory, oppHistory []int) int

// AlwaysCooperate is the IPD strategy that cooperates every round.
func AlwaysCooperate(myHistory, oppHistory []int) int { return Cooperate }

// AlwaysDefect is the IPD strategy that defects every round.
func AlwaysDefect(myHistory, oppHistory []int) int { return Defect }

// TitForTat is the IPD strategy that cooperates on the first round and
// thereafter copies the opponent's most recent move.
func TitForTat(myHistory, oppHistory []int) int {
	if len(oppHistory) == 0 {
		return Cooperate
	}
	return oppHistory[len(oppHistory)-1]
}

// GrimTrigger is the IPD strategy that cooperates until the opponent defects
// even once, after which it defects forever.
func GrimTrigger(myHistory, oppHistory []int) int {
	for _, m := range oppHistory {
		if m == Defect {
			return Defect
		}
	}
	return Cooperate
}

// WinStayLoseShift is the IPD strategy (also called Pavlov) that cooperates on
// the first round, then repeats its previous move after a "good" outcome (mutual
// cooperation or successfully defecting against a cooperator) and switches its
// move after a "bad" outcome.
func WinStayLoseShift(myHistory, oppHistory []int) int {
	if len(myHistory) == 0 {
		return Cooperate
	}
	last := len(myHistory) - 1
	myLast, oppLast := myHistory[last], oppHistory[last]
	good := (myLast == Cooperate && oppLast == Cooperate) ||
		(myLast == Defect && oppLast == Cooperate)
	if good {
		return myLast
	}
	if myLast == Cooperate {
		return Defect
	}
	return Cooperate
}

// PlayIteratedPD plays the iterated Prisoner's dilemma between strategies a and
// b for the given number of rounds and returns each player's total accumulated
// payoff along with the full move histories (oldest to newest). Each round both
// strategies choose simultaneously based on the history so far; payoffs use the
// PrisonersDilemma's T, R, P, S. rounds must be non-negative.
func PlayIteratedPD(pd PrisonersDilemma, a, b IPDStrategy, rounds int) (scoreA, scoreB float64, historyA, historyB []int) {
	historyA = make([]int, 0, rounds)
	historyB = make([]int, 0, rounds)
	for r := 0; r < rounds; r++ {
		moveA := a(historyA, historyB)
		moveB := b(historyB, historyA)
		pa, pb := gametheoryPDPayoff(pd, moveA, moveB)
		scoreA += pa
		scoreB += pb
		historyA = append(historyA, moveA)
		historyB = append(historyB, moveB)
	}
	return scoreA, scoreB, historyA, historyB
}

// gametheoryPDPayoff returns the payoffs to the two players for a single round
// given their moves.
func gametheoryPDPayoff(pd PrisonersDilemma, moveA, moveB int) (float64, float64) {
	switch {
	case moveA == Cooperate && moveB == Cooperate:
		return pd.R, pd.R
	case moveA == Cooperate && moveB == Defect:
		return pd.S, pd.T
	case moveA == Defect && moveB == Cooperate:
		return pd.T, pd.S
	default:
		return pd.P, pd.P
	}
}
