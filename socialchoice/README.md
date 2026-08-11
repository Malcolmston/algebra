# socialchoice — implements classical social-choice theory: ranked and cardinal voting rules

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/socialchoice.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/socialchoice)

Package socialchoice implements classical social-choice theory: ranked and
cardinal voting rules, pairwise/tournament analysis, and seat apportionment.

- positional rules: Plurality, AntiPlurality, Borda, Dowdall and the generic
PositionalScores, plus the iterative Borda methods Nanson and Baldwin; -
pairwise analysis: the PairwiseMatrix of majority preferences, Condorcet
winners and losers, Copeland and MiniMax scores, and the Smith and Schwartz
sets computed from the majority tournament; - Condorcet-consistent
completions: the Schulze beatpath method and Tideman's ranked pairs; -
elimination rules: instant-runoff voting, the single transferable vote with
Droop-quota Gregory surplus transfers, Coombs, Bucklin, the contingent and
supplementary votes, and the two-round runoff; - cardinal rules: approval,
range (score), STAR, cumulative and majority-judgment voting.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/socialchoice"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func Adams(votes []int, seats int) ([]int, error)` | Adams allocates seats by the Adams method, divisor s, so every party with positive votes receives at least one seat when seats permit. |
| `func ArgMaxFloat(xs []float64) int` | ArgMaxFloat returns the index of the maximum value, breaking ties toward the lowest index, and -1 for an empty slice. |
| `func ArgMaxInt(xs []int) int` | ArgMaxInt returns the index of the maximum value, breaking ties toward the lowest index, and -1 for an empty slice. |
| `func ArgMinFloat(xs []float64) int` | ArgMinFloat returns the index of the minimum value, breaking ties toward the lowest index, and -1 for an empty slice. |
| `func DHondt(votes []int, seats int) ([]int, error)` | DHondt allocates seats by the D'Hondt (Jefferson) method, divisor s+1. |
| `func Dean(votes []int, seats int) ([]int, error)` | Dean allocates seats by the Dean method, whose divisor is the harmonic mean of s and s+1. |
| `func DroopLargestRemainder(votes []int, seats int) ([]int, error)` | DroopLargestRemainder allocates seats by the largest-remainder method with the Hagenbach–Bischoff (exact Droop) quota total/(seats+1). |
| `func DroopQuota(total, seats int) int` | DroopQuota returns the integer Droop quota floor(total/(seats+1))+1 used by the single transferable vote. |
| `func EffectiveNumberOfParties(counts []int) float64` | EffectiveNumberOfParties returns the Laakso–Taagepera effective number of parties for the given vote (or seat) shares: the inverse of the sum of… |
| `func ExactQuotas(votes []int, seats int) []float64` | ExactQuotas returns each party's exact seat entitlement votes*seats/total. |
| `func GallagherIndex(votes, seats []int) float64` | GallagherIndex returns the Gallagher (least-squares) disproportionality index between vote shares and seat shares, expressed as a percentage in… |
| `func HagenbachBischoffQuota(total, seats int) float64` | HagenbachBischoffQuota returns the exact Hagenbach–Bischoff quota total/(seats+1). |
| `func Hamilton(votes []int, seats int) ([]int, error)` | Hamilton allocates seats by Hamilton's method, the largest-remainder method with the Hare quota total/seats. |
| `func HareQuota(total, seats int) float64` | HareQuota returns the Hare quota total/seats as a float. |
| `func HasAlabamaParadox(votes []int, seats int) bool` | HasAlabamaParadox reports whether increasing the house size from seats to seats+1 causes any party to lose a seat under Hamilton's method — the… |
| `func HighestAverages(votes []int, seats int, divisor DivisorFunc) ([]int, error)` | HighestAverages allocates seats among parties by the generic highest-averages (divisor) method: it repeatedly awards the next seat to the party with… |
| `func HuntingtonHill(votes []int, seats int) ([]int, error)` | HuntingtonHill allocates seats by the Huntington–Hill (equal proportions) method, whose divisor is the geometric mean sqrt(s(s+1)); every party… |
| `func InSet(set []int, c int) bool` | InSet reports whether candidate c appears in the (sorted or unsorted) set. |
| `func IsHouseMonotone(votes []int, maxSeats int, divisor DivisorFunc) bool` | IsHouseMonotone reports whether the given divisor method is free of the Alabama paradox across house sizes 1..maxSeats for these votes; every divisor… |
| `func Jefferson(votes []int, seats int) ([]int, error)` | Jefferson is an alias for DHondt. |
| `func LargestRemainder(votes []int, seats int, quota float64) ([]int, error)` | LargestRemainder allocates seats by the largest-remainder (quota) method for a given quota: each party first receives the floor of its votes/quota,… |
| `func LoosemoreHanbyIndex(votes, seats []int) float64` | LoosemoreHanbyIndex returns the Loosemore–Hanby disproportionality index, half the summed absolute differences between vote and seat shares, as a… |
| `func LowerQuotas(votes []int, seats int) []int` | LowerQuotas returns the floor of each party's exact seat entitlement. |
| `func MajorityThreshold(total int) int` | MajorityThreshold returns the smallest tally that constitutes a strict majority of total voters, floor(total/2)+1. |
| `func ModifiedSainteLague(votes []int, seats int) ([]int, error)` | ModifiedSainteLague allocates seats by the modified Sainte-Laguë method, whose first divisor is 1.4 and whose later divisors are 3, 5, 7, …. |
| `func RankingFromScores(scores []float64) []int` | RankingFromScores returns candidate indices ordered by descending score, breaking ties toward the lower index. |
| `func SainteLague(votes []int, seats int) ([]int, error)` | SainteLague allocates seats by the Sainte-Laguë (Webster) method, divisor 2s+1 (equivalently the odd-number divisors 1, 3, 5, …). |
| `func SatisfiesQuota(alloc, votes []int, seats int) bool` | SatisfiesQuota reports whether every party's allocation lies within its lower and upper quota, i.e. |
| `func SortedInts(xs []int) []int` | SortedInts returns a sorted copy of the slice, a small convenience used when comparing candidate sets. |
| `func TotalSeats(alloc []int) int` | TotalSeats returns the sum of an allocation. |
| `func UpperQuotas(votes []int, seats int) []int` | UpperQuotas returns the ceiling of each party's exact seat entitlement. |
| `func ViolatesLowerQuota(alloc, votes []int, seats int) bool` | ViolatesLowerQuota reports whether any party receives fewer seats than the floor of its exact entitlement. |
| `func ViolatesUpperQuota(alloc, votes []int, seats int) bool` | ViolatesUpperQuota reports whether any party receives more seats than the ceiling of its exact entitlement. |
| `func Webster(votes []int, seats int) ([]int, error)` | Webster is an alias for SainteLague. |
| `func WinnerFromScores(scores []float64) int` | WinnerFromScores returns the index of the highest-scoring candidate, breaking ties toward the lower index. |

### Types

| Type | What it is |
| --- | --- |
| `ApprovalProfile` | ApprovalProfile is a set of approval ballots over Candidates candidates; each ballot approves the candidates whose entry is true. |
| `Ballot` | Ballot is a strict preference order over candidate indices, most preferred first. |
| `CumulativeProfile` | CumulativeProfile is a set of cumulative-voting ballots; Ballots`v``c` is the number of points voter v allocates to candidate c. |
| `DivisorFunc` | DivisorFunc maps the number of seats a party currently holds to the divisor applied to its votes when computing the next-seat priority in a highest-… |
| `JudgmentProfile` | JudgmentProfile is a set of graded ballots for majority-judgment voting; Grades`v``c` is the grade in 0..NumGrades-1 that voter v assigns candidate… |
| `MajorityPair` | MajorityPair is a directed pairwise victory of Winner over Loser with the given majority margin and winning-vote count, used by ranked pairs. |
| `PairwiseMatrix` | PairwiseMatrix holds pairwise majority tallies: entry `i``j` is the number of voters who strictly prefer candidate i to candidate j. |
| `Profile` | Profile is a weighted collection of ranked ballots over a fixed number of candidates identified by the indices 0..Candidates-1. |
| `RunoffRound` | RunoffRound records the continuing candidates and their tallies at one round of an elimination rule, together with the candidate eliminated that… |
| `ScoreProfile` | ScoreProfile is a set of range/score ballots; Ballots`v``c` is the score voter v assigns candidate c, bounded above by Max. |

<details>
<summary><code>ApprovalProfile</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ApprovalFromRanked(p *Profile, k int) ApprovalProfile` | ApprovalFromRanked derives an approval profile from a ranked profile by approving each ballot's top k candidates. |
| `func (a ApprovalProfile) NumBallots() int` | NumBallots returns the number of approval ballots. |
| `func (a ApprovalProfile) Ranking() []int` | Ranking returns candidates ordered by descending approval. |
| `func (a ApprovalProfile) Scores() []float64` | Scores returns each candidate's approval tally, the number of ballots approving it. |
| `func (a ApprovalProfile) TotalApprovals() int` | TotalApprovals returns the total number of approvals cast across all ballots. |
| `func (a ApprovalProfile) Winner() int` | Winner returns the most approved candidate. |

</details>

<details>
<summary><code>Ballot</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (b Ballot) Bottom() (int, bool)` | Bottom returns the least preferred listed candidate and whether the ballot is non-empty. |
| `func (b Ballot) BottomAmong(continuing []bool) int` | BottomAmong returns the least preferred listed candidate for which continuing`c` is true, or -1 if none is listed. |
| `func (b Ballot) Clone() Ballot` | Clone returns a copy of the ballot. |
| `func (b Ballot) Position(c int) (int, bool)` | Position returns the zero-based rank of candidate c on the ballot and whether c appears at all. |
| `func (b Ballot) Prefers(a, c int) bool` | Prefers reports whether the ballot strictly prefers candidate a to candidate b. |
| `func (b Ballot) Ranks(c int) bool` | Ranks reports whether candidate c appears anywhere on the ballot. |
| `func (b Ballot) Reversed() Ballot` | Reversed returns a new ballot with the preference order inverted. |
| `func (b Ballot) Top() (int, bool)` | Top returns the most preferred candidate on the ballot and whether the ballot is non-empty. |
| `func (b Ballot) TopAmong(continuing []bool) int` | TopAmong returns the most preferred candidate on the ballot for which continuing`c` is true, or -1 if the ballot lists no continuing candidate. |

</details>

<details>
<summary><code>CumulativeProfile</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (c CumulativeProfile) Scores() []float64` | Scores returns the total points each candidate accumulates. |
| `func (c CumulativeProfile) Winner() int` | Winner returns the candidate with the most cumulative points. |

</details>

<details>
<summary><code>JudgmentProfile</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (j JudgmentProfile) MedianGrades() []int` | MedianGrades returns each candidate's lower-median grade. |
| `func (j JudgmentProfile) Winner() int` | Winner returns the majority-judgment winner: the highest median grade wins, and ties are broken by repeatedly removing one median grade from each… |

</details>

<details>
<summary><code>PairwiseMatrix</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (m PairwiseMatrix) Beats(i, j int) bool` | Beats reports whether a strict majority prefers i to j. |
| `func (m PairwiseMatrix) CondorcetLoser() (int, bool)` | CondorcetLoser returns a candidate that loses to every other candidate pairwise and true, or 0 and false when none exists. |
| `func (m PairwiseMatrix) CondorcetWinner() (int, bool)` | CondorcetWinner returns a candidate that beats every other candidate pairwise and true, or 0 and false when none exists. |
| `func (m PairwiseMatrix) CopelandNet() []int` | CopelandNet returns each candidate's net Copeland score, wins minus losses. |
| `func (m PairwiseMatrix) CopelandScores() []float64` | CopelandScores returns each candidate's Copeland score, counting 1 for every pairwise win and 0.5 for every pairwise tie. |
| `func (m PairwiseMatrix) Covers(i, j int) bool` | Covers reports whether candidate i covers candidate j: i beats j and beats every candidate that j beats. |
| `func (m PairwiseMatrix) HasCondorcetCycle() bool` | HasCondorcetCycle reports whether the majority relation contains a cycle, equivalently whether the Smith set has more than one candidate. |
| `func (m PairwiseMatrix) HasTie() bool` | HasTie reports whether any pair of distinct candidates is pairwise tied. |
| `func (m PairwiseMatrix) IsTournament() bool` | IsTournament reports whether the majority relation is a strict tournament: no pair of distinct candidates is tied. |
| `func (m PairwiseMatrix) KemenyRanking() []int` | KemenyRanking returns a Kemeny–Young optimal ranking, found by exhaustive search over all orderings (suitable for a modest number of candidates). |
| `func (m PairwiseMatrix) KemenyScore(ranking []int) int` | KemenyScore returns the Kemeny agreement score of a full ranking: the summed pairwise support for every ordered pair placed in agreement with the… |
| `func (m PairwiseMatrix) Losses(i int) int` | Losses returns the number of candidates that beat i pairwise. |
| `func (m PairwiseMatrix) Margin(i, j int) int` | Margin returns the signed majority margin of i over j, m`i``j`-m`j``i`. |
| `func (m PairwiseMatrix) MiniMaxScores() []int` | MiniMaxScores returns, for each candidate, the largest margin by which it is defeated in any pairwise contest (its worst pairwise loss); a value of… |
| `func (m PairwiseMatrix) N() int` | N returns the number of candidates covered by the matrix. |
| `func (m PairwiseMatrix) NetMargins() [][]int` | NetMargins returns the antisymmetric margin matrix, entry `i``j` = margin of i over j. |
| `func (m PairwiseMatrix) Prefer(i, j int) int` | Prefer returns the number of voters preferring i to j. |
| `func (m PairwiseMatrix) RankedPairsLocked() [][]bool` | RankedPairsLocked returns the acyclic set of locked pairwise victories under Tideman's ranked-pairs method: majorities are considered strongest first… |
| `func (m PairwiseMatrix) RankedPairsRanking() []int` | RankedPairsRanking returns the ranked-pairs ordering: a topological order of the locked-victory graph, most preferred first. |
| `func (m PairwiseMatrix) RankedPairsWinner() int` | RankedPairsWinner returns the source of the locked-victory graph, the ranked-pairs winner. |
| `func (m PairwiseMatrix) SchulzeRanking() []int` | SchulzeRanking returns candidates ordered by descending number of Schulze pairwise victories, breaking ties toward the lower index. |
| `func (m PairwiseMatrix) SchulzeWinner() int` | SchulzeWinner returns a Schulze-method winner: a candidate whose strongest path to every other candidate is at least as strong as the reverse. |
| `func (m PairwiseMatrix) SchulzeWins() []int` | SchulzeWins returns each candidate's number of Schulze pairwise victories: the count of opponents j for which the strongest path from i to j is… |
| `func (m PairwiseMatrix) SchwartzSet() []int` | SchwartzSet returns the Schwartz set: the union of the minimal undominated sets under the strict beat relation (ties do not dominate). |
| `func (m PairwiseMatrix) SmithSet() []int` | SmithSet returns the Smith set: the smallest non-empty set of candidates such that every member beats every non-member pairwise. |
| `func (m PairwiseMatrix) SortedMajorities() []MajorityPair` | SortedMajorities returns every pairwise victory in the matrix ordered by descending margin; ties in margin are broken by descending winning votes and… |
| `func (m PairwiseMatrix) StrongestPaths() [][]int` | StrongestPaths returns the Schulze strongest-beatpath matrix p, where p`i``j` is the strength of the strongest path from i to j. |
| `func (m PairwiseMatrix) TieCount(i int) int` | TieCount returns the number of candidates pairwise tied with i (excluding i). |
| `func (m PairwiseMatrix) Ties(i, j int) bool` | Ties reports whether i and j are pairwise tied. |
| `func (m PairwiseMatrix) TopCycle() []int` | TopCycle returns the Smith set, also known as the top cycle of the majority tournament. |
| `func (m PairwiseMatrix) Transpose() PairwiseMatrix` | Transpose returns the transpose of the pairwise matrix. |
| `func (m PairwiseMatrix) UncoveredSet() []int` | UncoveredSet returns the uncovered set: the candidates that no other candidate covers. |
| `func (m PairwiseMatrix) Wins(i int) int` | Wins returns the number of candidates that i beats pairwise. |

</details>

<details>
<summary><code>Profile</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewProfile(candidates int, ballots []Ballot, counts []int) (*Profile, error)` | NewProfile builds a Profile over the given number of candidates and validates it. |
| `func (p *Profile) AddBallot(b Ballot, weight int) *Profile` | AddBallot appends a ballot with the given weight and returns the updated profile pointer for chaining. |
| `func (p *Profile) AntiPluralityRanking() []int` | AntiPluralityRanking returns candidates ordered by descending anti-plurality score. |
| `func (p *Profile) AntiPluralityScores() []float64` | AntiPluralityScores returns each candidate's anti-plurality (veto) score: the total weight of ballots that do not rank the candidate last. |
| `func (p *Profile) AntiPluralityWinner() int` | AntiPluralityWinner returns the candidate vetoed by the fewest ballots. |
| `func (p *Profile) BaldwinWinner() int` | BaldwinWinner returns the winner of the Baldwin rule: repeatedly eliminate the single continuing candidate with the lowest Borda score until one… |
| `func (p *Profile) BlackWinner() int` | BlackWinner returns the Black-method winner: the Condorcet winner when one exists, otherwise the Borda winner. |
| `func (p *Profile) BordaRanking() []int` | BordaRanking returns candidates ordered by descending Borda score. |
| `func (p *Profile) BordaScore(c int) float64` | BordaScore returns the Borda score of a single candidate. |
| `func (p *Profile) BordaScores() []float64` | BordaScores returns each candidate's Borda score, the total weight of voters preferring the candidate to each other candidate summed over all… |
| `func (p *Profile) BordaWinner() int` | BordaWinner returns the highest Borda-scoring candidate. |
| `func (p *Profile) BucklinScores(depth int) []float64` | BucklinScores returns each candidate's Bucklin tally at the given depth: the weight of ballots ranking the candidate within their top (depth+1)… |
| `func (p *Profile) BucklinWinner() int` | BucklinWinner returns the winner of the Bucklin rule: deepen the count one preference level at a time until at least one candidate's tally reaches a… |
| `func (p *Profile) CondorcetLoser() (int, bool)` | CondorcetLoser returns the profile's Condorcet loser, if any. |
| `func (p *Profile) CondorcetParadox() bool` | CondorcetParadox reports whether the profile exhibits the Condorcet paradox: no Condorcet winner exists. |
| `func (p *Profile) CondorcetWinner() (int, bool)` | CondorcetWinner returns the profile's Condorcet winner, if any. |
| `func (p *Profile) ContingentVoteWinner() int` | ContingentVoteWinner returns the winner of the contingent vote: if a candidate holds a first-preference majority they win, otherwise all but the top… |
| `func (p *Profile) CoombsWinner() int` | CoombsWinner returns the winner of the Coombs rule: repeatedly eliminate the continuing candidate with the most active last-place votes until a… |
| `func (p *Profile) CopelandNet() []int` | CopelandNet returns the profile's net Copeland scores (wins minus losses). |
| `func (p *Profile) CopelandRanking() []int` | CopelandRanking returns candidates ordered by descending Copeland score. |
| `func (p *Profile) CopelandScores() []float64` | CopelandScores returns the profile's Copeland scores (win=1, tie=0.5). |
| `func (p *Profile) CopelandWinner() int` | CopelandWinner returns the highest Copeland-scoring candidate. |
| `func (p *Profile) Copy() *Profile` | Copy returns a deep copy of the profile. |
| `func (p *Profile) DowdallRanking() []int` | DowdallRanking returns candidates ordered by descending Dowdall score. |
| `func (p *Profile) DowdallScores() []float64` | DowdallScores returns each candidate's Dowdall (Nauru) score, awarding 1/(k+1) points for a rank at zero-based position k. |
| `func (p *Profile) DowdallWinner() int` | DowdallWinner returns the highest Dowdall-scoring candidate. |
| `func (p *Profile) FirstPreferenceCounts() []float64` | FirstPreferenceCounts is an alias for PluralityScores, returning the top-choice tally of every candidate. |
| `func (p *Profile) HasCondorcetCycle() bool` | HasCondorcetCycle reports whether the profile's majority relation is cyclic. |
| `func (p *Profile) HasCondorcetWinner() bool` | HasCondorcetWinner reports whether the profile has a Condorcet winner. |
| `func (p *Profile) HasMajorityWinner() (int, bool)` | HasMajorityWinner reports whether some candidate holds a strict majority of first preferences, returning that candidate and true when so. |
| `func (p *Profile) InstantRunoffRounds() ([]RunoffRound, int)` | InstantRunoffRounds runs instant-runoff voting (IRV) and returns the ordered rounds together with the winner. |
| `func (p *Profile) InstantRunoffWinner() int` | InstantRunoffWinner returns the IRV winner of the profile. |
| `func (p *Profile) IsCondorcetWinner(c int) bool` | IsCondorcetWinner reports whether candidate c beats every other candidate pairwise in the profile. |
| `func (p *Profile) KemenyRanking() []int` | KemenyRanking returns the profile's Kemeny–Young optimal ranking. |
| `func (p *Profile) KemenyWinner() int` | KemenyWinner returns the top candidate of the Kemeny–Young optimal ranking. |
| `func (p *Profile) LastPreferenceCounts() []float64` | LastPreferenceCounts returns each candidate's tally of last-place (bottom listed) rankings, the basis of the Coombs rule. |
| `func (p *Profile) MajorityWinner() int` | MajorityWinner returns the candidate holding a first-preference majority, or -1 when none does. |
| `func (p *Profile) MiniMaxScores() []int` | MiniMaxScores returns the profile's minimax worst-defeat margins. |
| `func (p *Profile) MiniMaxWinner() int` | MiniMaxWinner returns the candidate whose worst pairwise defeat margin is smallest, breaking ties toward the lower index. |
| `func (p *Profile) NansonWinner() int` | NansonWinner returns the winner of the Nanson rule: repeatedly eliminate every continuing candidate whose Borda score is strictly below the average… |
| `func (p *Profile) NumBallots() int` | NumBallots returns the number of distinct ballots recorded in the profile. |
| `func (p *Profile) Pairwise() PairwiseMatrix` | Pairwise builds the PairwiseMatrix of the profile, summing ballot weights. |
| `func (p *Profile) PairwiseMargin(i, j int) int` | PairwiseMargin returns the signed majority margin of i over j in the profile. |
| `func (p *Profile) PluralityRanking() []int` | PluralityRanking returns candidates ordered by descending first-preference tally. |
| `func (p *Profile) PluralityScore(c int) float64` | PluralityScore returns the first-preference tally of a single candidate. |
| `func (p *Profile) PluralityScores() []float64` | PluralityScores returns each candidate's first-preference tally: the total weight of ballots ranking that candidate top. |
| `func (p *Profile) PluralityWinner() int` | PluralityWinner returns the candidate with the most first preferences. |
| `func (p *Profile) PositionalRanking(weights []float64) []int` | PositionalRanking returns candidates ordered by descending positional score under the given weights. |
| `func (p *Profile) PositionalScores(weights []float64) []float64` | PositionalScores returns each candidate's score under an arbitrary positional rule: a candidate ranked at position k on a ballot earns weights`k`… |
| `func (p *Profile) PositionalWinner(weights []float64) int` | PositionalWinner returns the winner under the positional rule with the given position weights. |
| `func (p *Profile) RankedPairsRanking() []int` | RankedPairsRanking returns the profile's ranked-pairs ordering. |
| `func (p *Profile) RankedPairsWinner() int` | RankedPairsWinner returns the profile's ranked-pairs (Tideman) winner. |
| `func (p *Profile) RemoveCandidate(c int) *Profile` | RemoveCandidate returns a new profile identical to p but with candidate c removed from every ballot; the candidate count is unchanged so remaining… |
| `func (p *Profile) STV(seats int) []int` | STV runs the single transferable vote with the Droop quota and fractional (Gregory) surplus transfers, filling the requested number of seats. |
| `func (p *Profile) SchulzeRanking() []int` | SchulzeRanking returns the profile's Schulze ranking. |
| `func (p *Profile) SchulzeWinner() int` | SchulzeWinner returns the profile's Schulze-method winner. |
| `func (p *Profile) SchulzeWins() []int` | SchulzeWins returns the profile's per-candidate Schulze victory counts. |
| `func (p *Profile) SchwartzSet() []int` | SchwartzSet returns the profile's Schwartz set. |
| `func (p *Profile) SmithSet() []int` | SmithSet returns the profile's Smith set. |
| `func (p *Profile) SupplementaryVoteWinner() int` | SupplementaryVoteWinner returns the winner of the supplementary vote, the contingent vote restricted to each ballot's first two preferences. |
| `func (p *Profile) TotalVoters() int` | TotalVoters returns the sum of the ballot weights, i.e. |
| `func (p *Profile) TwoRoundWinner() int` | TwoRoundWinner returns the winner of the two-round (runoff) system: a first-preference majority wins outright, otherwise the top two first-round… |
| `func (p *Profile) UncoveredSet() []int` | UncoveredSet returns the profile's uncovered set. |
| `func (p *Profile) Validate() error` | Validate reports whether the profile is well formed: a positive candidate count, in-range and duplicate-free ballot entries, and non-negative counts… |
| `func (p *Profile) Weight(i int) int` | Weight returns the number of voters casting ballot i, defaulting to 1 when Counts does not cover that ballot. |

</details>

<details>
<summary><code>ScoreProfile</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func RangeFromRanked(p *Profile) ScoreProfile` | RangeFromRanked derives a score profile from a ranked profile using Borda-style position scores: a candidate at zero-based rank r on an n-candidate… |
| `func (s ScoreProfile) Averages() []float64` | Averages returns each candidate's mean score over all ballots. |
| `func (s ScoreProfile) NumBallots() int` | NumBallots returns the number of score ballots. |
| `func (s ScoreProfile) Ranking() []int` | Ranking returns candidates ordered by descending total score. |
| `func (s ScoreProfile) STARWinner() int` | STARWinner returns the STAR (Score Then Automatic Runoff) winner: the two highest-scoring candidates advance to an automatic runoff won by whichever… |
| `func (s ScoreProfile) Totals() []float64` | Totals returns the sum of scores each candidate receives. |
| `func (s ScoreProfile) Winner() int` | Winner returns the highest total-score candidate (range/score voting). |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/socialchoice).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
