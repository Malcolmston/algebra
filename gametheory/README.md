# gametheory — implements classical game theory: normal-form games, Nash equilibria (pure and mixed)

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/gametheory.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/gametheory)

Package gametheory implements classical game theory: normal-form games, Nash
equilibria (pure and mixed), zero-sum solutions via linear programming,
dominated-strategy elimination, Pareto frontiers, cooperative solution
concepts (the Shapley and Banzhaf values, the core), and a small toolkit of
Prisoner's-dilemma and canonical-game helpers.

- equilibrium analysis: PureNashEquilibria and, via support enumeration, the
full set of mixed equilibria (MixedNashEquilibria); - zero-sum solving:
SolveZeroSum computes the value of the game and a pair of optimal minimax
strategies by linear programming (a self-contained Bland's-rule simplex), with
Maximin, Minimax and SaddlePoint helpers for the pure case; - dominance:
strict and weak domination tests and iterated elimination; - efficiency:
Pareto frontiers and the social optimum; - cooperation: CooperativeGame with
the Shapley value, the Banzhaf value, core membership tests and structural
predicates (superadditivity, monotonicity, convexity), plus
WeightedVotingGame; - the Prisoner's dilemma: PrisonersDilemma with
iterated-play simulation and stock strategies (TitForTat, GrimTrigger,
WinStayLoseShift, ...), and constructors for the canonical games
MatchingPennies, StagHunt, BattleOfTheSexes and GameOfChicken.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/gametheory"
```

## Usage

This is the package's own `ExampleCooperativeGame_IsInCore`, so it compiles and its output is
asserted on every `go test ./gametheory/`.

```go
cg, err := gametheory.NewCooperativeGameFromValues([]float64{0, 1, 1, 3})
	if err != nil {
		panic(err)
	}
	fmt.Println(cg.IsInCore([]float64{1.5, 1.5}, 1e-9))
	fmt.Println(cg.IsInCore([]float64{3, 0}, 1e-9))
```

```
true
false
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AlwaysCooperate(myHistory, oppHistory []int) int` | AlwaysCooperate is the IPD strategy that cooperates every round. |
| `func AlwaysDefect(myHistory, oppHistory []int) int` | AlwaysDefect is the IPD strategy that defects every round. |
| `func ExpectedPayoff(a [][]float64, p, q MixedStrategy) float64` | ExpectedPayoff returns the bilinear form p^T A q, the expected value of the payoff matrix A when the row is drawn from distribution p and the column… |
| `func GameValue(payoff [][]float64) (float64, error)` | GameValue returns the value (to the maximizing row player) of the matrix game with the given row-player payoff matrix. |
| `func GrimTrigger(myHistory, oppHistory []int) int` | GrimTrigger is the IPD strategy that cooperates until the opponent defects even once, after which it defects forever. |
| `func Maximin(payoff [][]float64) (int, float64)` | Maximin returns the row player's pure maximin strategy and its security value: the row that maximizes the row's worst-case (minimum over columns)… |
| `func Minimax(payoff [][]float64) (int, float64)` | Minimax returns the column player's pure minimax strategy and its value: the column that minimizes the column's worst-case (maximum over rows) payoff… |
| `func ParetoDominates(a, b []float64) bool` | ParetoDominates reports whether payoff vector a Pareto-dominates payoff vector b: a is at least as large as b in every coordinate and strictly larger… |
| `func ParetoOptimalIndices(points [][]float64) []int` | ParetoOptimalIndices returns the sorted indices of the Pareto-optimal points among the given points (each a vector of objectives, all maximized): a… |
| `func PlayIteratedPD(pd PrisonersDilemma, a, b IPDStrategy, rounds int) (scoreA, scoreB float64, historyA, historyB []int)` | PlayIteratedPD plays the iterated Prisoner's dilemma between strategies a and b for the given number of rounds and returns each player's total… |
| `func SaddlePoint(payoff [][]float64) (int, int, bool)` | SaddlePoint reports whether the matrix game has a pure-strategy saddle point, i.e. |
| `func TitForTat(myHistory, oppHistory []int) int` | TitForTat is the IPD strategy that cooperates on the first round and thereafter copies the opponent's most recent move. |
| `func WinStayLoseShift(myHistory, oppHistory []int) int` | WinStayLoseShift is the IPD strategy (also called Pavlov) that cooperates on the first round, then repeats its previous move after a "good" outcome… |

### Types

| Type | What it is |
| --- | --- |
| `CooperativeGame` | CooperativeGame is a transferable-utility cooperative game on a fixed set of players {0, 1, ..., Players-1}. |
| `Game` | Game is a finite two-player game in normal (strategic) form. |
| `IPDStrategy` | IPDStrategy is a deterministic strategy for the iterated Prisoner's dilemma. |
| `MixedEquilibrium` | MixedEquilibrium is a mixed-strategy Nash equilibrium: Row and Col are the two players' equilibrium mixed strategies, and RowValue and ColValue are… |
| `MixedStrategy` | MixedStrategy is a probability distribution over a player's pure strategies: entry i is the probability of playing pure strategy i. |
| `PrisonersDilemma` | PrisonersDilemma holds the four payoffs of a symmetric Prisoner's dilemma from one player's perspective: T is the temptation payoff (defect against a… |
| `PureProfile` | PureProfile is a pair of pure strategies, one for each player: Row is the row player's strategy index and Col the column player's. |
| `ZeroSumSolution` | ZeroSumSolution is the solution of a two-player zero-sum game: Value is the value of the game to the row (maximizing) player, RowStrategy is an… |

<details>
<summary><code>CooperativeGame</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewCooperativeGame(players int, characteristic func(coalition uint) float64) (CooperativeGame, error)` | NewCooperativeGame builds a CooperativeGame with the given number of players by evaluating the characteristic function on every one of the 2^players… |
| `func NewCooperativeGameFromValues(values []float64) (CooperativeGame, error)` | NewCooperativeGameFromValues builds a CooperativeGame from a fully specified worth table, where values`mask` is v(mask) for each coalition bitmask. |
| `func WeightedVotingGame(weights []float64, quota float64) (CooperativeGame, error)` | WeightedVotingGame constructs the simple cooperative game in which a coalition wins (worth 1) when the total weight of its members is at least quota,… |
| `func (cg CooperativeGame) BanzhafValue() []float64` | BanzhafValue returns the (non-normalized) Banzhaf value of the game: for each player, the average marginal contribution taken with equal weight over… |
| `func (cg CooperativeGame) GrandCoalitionValue() float64` | GrandCoalitionValue returns v(N), the worth of the coalition of all players. |
| `func (cg CooperativeGame) IsConvex(tol float64) bool` | IsConvex reports whether the game is convex (supermodular): for all coalitions S and T, v(S ∪ T) + v(S ∩ T) >= v(S) + v(T) within tolerance tol. |
| `func (cg CooperativeGame) IsEssential(tol float64) bool` | IsEssential reports whether the game is essential: the grand coalition is worth strictly more (by more than tol) than the players acting alone, i.e. |
| `func (cg CooperativeGame) IsInCore(x []float64, tol float64) bool` | IsInCore reports whether allocation x is in the core of the game: it must be efficient (Σ x_i = v(N)) and coalitionally rational (Σ_{i∈S} x_i >=… |
| `func (cg CooperativeGame) IsMonotone(tol float64) bool` | IsMonotone reports whether the game is monotone: enlarging a coalition never decreases its worth, i.e. |
| `func (cg CooperativeGame) IsSuperadditive(tol float64) bool` | IsSuperadditive reports whether the game is superadditive: for every pair of disjoint coalitions S and T, v(S ∪ T) >= v(S) + v(T) within tolerance… |
| `func (cg CooperativeGame) MarginalContribution(player int, coalition uint) float64` | MarginalContribution returns v(coalition ∪ {player}) - v(coalition), the value player adds when joining the coalition. |
| `func (cg CooperativeGame) ShapleyValue() []float64` | ShapleyValue returns the Shapley value of the game: a fair allocation of the grand coalition's worth in which each player receives the… |
| `func (cg CooperativeGame) Worth(coalition uint) float64` | Worth returns v(coalition), the worth of the coalition encoded by the given bitmask. |

</details>

<details>
<summary><code>Game</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BattleOfTheSexes() Game` | BattleOfTheSexes returns the Battle of the Sexes coordination game. |
| `func GameOfChicken() Game` | GameOfChicken returns the game of Chicken (also called Hawk-Dove). |
| `func MatchingPennies() Game` | MatchingPennies returns the classic zero-sum game of matching pennies. |
| `func NewGame(row, col [][]float64) (Game, error)` | NewGame constructs a Game from the row and column payoff matrices, which must be non-empty, rectangular, and of identical dimensions. |
| `func NewZeroSumGame(payoff [][]float64) (Game, error)` | NewZeroSumGame constructs a zero-sum Game from a single payoff matrix interpreted as the row player's payoff; the column player's payoff is its… |
| `func StagHunt() Game` | StagHunt returns the Stag Hunt coordination game with the common payoffs (Stag,Stag)=(4,4), (Stag,Hare)=(1,3), (Hare,Stag)=(3,1), (Hare,Hare)=(3,3),… |
| `func (g Game) BestResponseColValue(p MixedStrategy) float64` | BestResponseColValue returns the maximum expected payoff the column player can obtain against the row player's mixed strategy p. |
| `func (g Game) BestResponseRowValue(q MixedStrategy) float64` | BestResponseRowValue returns the maximum expected payoff the row player can obtain against the column player's mixed strategy q. |
| `func (g Game) ColBestResponses(p MixedStrategy, tol float64) []int` | ColBestResponses returns the sorted column-player pure strategies that maximize the column player's expected payoff against the row player's mixed… |
| `func (g Game) ColPayoff(i, j int) float64` | ColPayoff returns the column player's payoff when the row player plays i and the column player plays j. |
| `func (g Game) ColStrictlyDominates(j, l int) bool` | ColStrictlyDominates reports whether column strategy j strictly dominates column strategy l for the column player: Col`i``j` > Col`i``l` for every… |
| `func (g Game) ColWeaklyDominates(j, l int) bool` | ColWeaklyDominates reports whether column strategy j weakly dominates column strategy l: Col`i``j` >= Col`i``l` for every row i, with strict… |
| `func (g Game) DominantColStrategy() (int, bool)` | DominantColStrategy returns a column strategy that strictly dominates every other column strategy, together with true, if one exists; otherwise it… |
| `func (g Game) DominantRowStrategy() (int, bool)` | DominantRowStrategy returns a row strategy that strictly dominates every other row strategy, together with true, if one exists; otherwise it returns… |
| `func (g Game) ExpectedColPayoff(p, q MixedStrategy) float64` | ExpectedColPayoff returns the column player's expected payoff when the row player plays mixed strategy p and the column player plays mixed strategy q. |
| `func (g Game) ExpectedRowPayoff(p, q MixedStrategy) float64` | ExpectedRowPayoff returns the row player's expected payoff when the row player plays mixed strategy p and the column player plays mixed strategy q. |
| `func (g Game) HasPureNash(tol float64) bool` | HasPureNash reports whether the game has at least one pure-strategy Nash equilibrium, testing deviations within tolerance tol. |
| `func (g Game) IsConstantSum(tol float64) (float64, bool)` | IsConstantSum reports whether the game is constant-sum, i.e. |
| `func (g Game) IsParetoOptimal(i, j int) bool` | IsParetoOptimal reports whether the pure-strategy profile (i, j) is Pareto-optimal among all pure profiles: no other profile yields both players at… |
| `func (g Game) IsPureNash(i, j int, tol float64) bool` | IsPureNash reports whether the pure-strategy profile (i, j) is a Nash equilibrium: neither player can strictly improve their payoff by a unilateral… |
| `func (g Game) IsSymmetric(tol float64) bool` | IsSymmetric reports whether the game is symmetric: it must be square and satisfy Col`i``j` == Row`j``i` for all i, j within tolerance tol, so that… |
| `func (g Game) IsZeroSum(tol float64) bool` | IsZeroSum reports whether the game is zero-sum, i.e. |
| `func (g Game) IteratedStrictDominance() (rows, cols []int)` | IteratedStrictDominance performs iterated elimination of strictly dominated pure strategies (by other pure strategies) until no further elimination… |
| `func (g Game) IteratedWeakDominance() (rows, cols []int)` | IteratedWeakDominance performs iterated elimination of weakly dominated pure strategies (by other pure strategies) until no further elimination is… |
| `func (g Game) MixedNashEquilibria(tol float64) []MixedEquilibrium` | MixedNashEquilibria returns mixed-strategy Nash equilibria of the two-player game found by support enumeration over all pairs of equal-sized supports. |
| `func (g Game) MixedNashEquilibriaWithSupport(rowSupport, colSupport []int, tol float64) (MixedEquilibrium, bool)` | MixedNashEquilibriaWithSupport attempts to find a mixed Nash equilibrium whose row-player support is rowSupport and whose column-player support is… |
| `func (g Game) NumColStrategies() int` | NumColStrategies returns n, the number of pure strategies available to the column player. |
| `func (g Game) NumRowStrategies() int` | NumRowStrategies returns m, the number of pure strategies available to the row player. |
| `func (g Game) ParetoFrontier() []PureProfile` | ParetoFrontier returns the pure-strategy profiles whose payoff pair (RowPayoff, ColPayoff) is Pareto-optimal: no other pure profile gives both… |
| `func (g Game) PureNashEquilibria(tol float64) []PureProfile` | PureNashEquilibria returns all pure-strategy Nash equilibria of the game as PureProfile values, in lexicographic (row, column) order. |
| `func (g Game) RowBestResponses(q MixedStrategy, tol float64) []int` | RowBestResponses returns the sorted row-player pure strategies that maximize the row player's expected payoff against the column player's mixed… |
| `func (g Game) RowPayoff(i, j int) float64` | RowPayoff returns the row player's payoff when the row player plays i and the column player plays j. |
| `func (g Game) RowStrictlyDominates(i, k int) bool` | RowStrictlyDominates reports whether row strategy i strictly dominates row strategy k for the row player: Row`i``j` > Row`k``j` for every column j. |
| `func (g Game) RowWeaklyDominates(i, k int) bool` | RowWeaklyDominates reports whether row strategy i weakly dominates row strategy k: Row`i``j` >= Row`k``j` for every column j, with strict inequality… |
| `func (g Game) SocialOptimum() (PureProfile, float64)` | SocialOptimum returns the pure-strategy profile maximizing the sum of the two players' payoffs (the utilitarian social welfare) and that maximal sum. |
| `func (g Game) SolveZeroSum() (ZeroSumSolution, error)` | SolveZeroSum solves the receiver as a zero-sum game, returning ErrNotZeroSum if Row`i``j`+Col`i``j` is not zero (within 1e-9) in every cell. |
| `func (g Game) StrictlyDominatedColStrategies() []int` | StrictlyDominatedColStrategies returns the sorted indices of column strategies that are strictly dominated by some other single (pure) column… |
| `func (g Game) StrictlyDominatedRowStrategies() []int` | StrictlyDominatedRowStrategies returns the sorted indices of row strategies that are strictly dominated by some other single (pure) row strategy. |
| `func (g Game) Transpose() Game` | Transpose returns the game viewed from the column player's seat: the players swap roles, so the returned game has n rows and m columns with row… |
| `func (g Game) Validate() error` | Validate reports whether the receiver's payoff matrices are well-formed (non-empty, rectangular and of equal shape), returning ErrShape otherwise. |

</details>

<details>
<summary><code>MixedStrategy</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewMixedStrategy(probs []float64) (MixedStrategy, error)` | NewMixedStrategy validates probs as a probability distribution (non-negative entries summing to one within a small tolerance) and returns a copy of… |
| `func PureStrategy(n, i int) MixedStrategy` | PureStrategy returns the degenerate mixed strategy over n pure strategies that plays pure strategy i with probability one. |
| `func UniformStrategy(n int) MixedStrategy` | UniformStrategy returns the uniform mixed strategy over n pure strategies, assigning probability 1/n to each. |
| `func (m MixedStrategy) Entropy() float64` | Entropy returns the Shannon entropy of the distribution in bits (base-2 log). |
| `func (m MixedStrategy) IsValid(tol float64) bool` | IsValid reports whether the receiver is a valid probability distribution: every entry is at least -tol and the entries sum to one within tol. |
| `func (m MixedStrategy) Support(tol float64) []int` | Support returns the sorted indices of pure strategies played with probability greater than tol. |

</details>

<details>
<summary><code>PrisonersDilemma</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewPrisonersDilemma(t, r, p, s float64) (PrisonersDilemma, error)` | NewPrisonersDilemma constructs a PrisonersDilemma, returning ErrPayoffOrder unless the payoffs satisfy the defining ordering t > r > p > s. |
| `func StandardPrisonersDilemma() PrisonersDilemma` | StandardPrisonersDilemma returns the canonical Prisoner's dilemma with the textbook payoffs T=5, R=3, P=1, S=0. |
| `func (pd PrisonersDilemma) Game() Game` | Game returns the two-player normal-form game for the Prisoner's dilemma, with row and column strategy 0 = Cooperate and 1 = Defect. |
| `func (pd PrisonersDilemma) IsIterationStable() bool` | IsIterationStable reports whether the payoffs satisfy 2R > T + S, the extra condition ensuring that in the iterated game sustained mutual cooperation… |
| `func (pd PrisonersDilemma) IsValid() bool` | IsValid reports whether the payoffs satisfy the Prisoner's-dilemma ordering T > R > P > S. |

</details>

<details>
<summary><code>ZeroSumSolution</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func SolveZeroSum(payoff [][]float64) (ZeroSumSolution, error)` | SolveZeroSum solves the matrix game whose payoff matrix (to the maximizing row player) is payoff, returning the value of the game and a pair of… |

</details>

### Constants

`Cooperate`, `Defect`

### Variables

`ErrDistribution`, `ErrNilCharacteristic`, `ErrNotZeroSum`, `ErrPayoffOrder`, `ErrPlayers`, `ErrShape`, `ErrSolve`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/gametheory).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
