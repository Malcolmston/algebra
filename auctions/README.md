# auctions — implements mechanism design and cooperative game theory in pure Go

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/auctions.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/auctions)

Package auctions implements mechanism design and cooperative game theory in
pure Go: sealed-bid auctions, the Vickrey-Clarke-Groves (VCG) mechanism with
combinatorial winner determination, and the classical solution concepts of
transferable-utility cooperative games together with two-sided matching.

- sealed-bid auctions: first-price, second-price (Vickrey), third-price and
all-pay single-item formats with optional reserve prices, plus multi-unit
uniform-price, pay-as-bid and multi-unit Vickrey clearings, and the
open-outcry English and Dutch formats through their strategic equivalences.
See Bid, FirstPriceAuction, VickreyAuction, EnglishAuction, DutchAuction and
MultiUnitVickrey.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/auctions"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AgentRank(prefs [][]int, agent, item int) int` | AgentRank returns the position (0 = most preferred) of item in agent's preference list, or -1 if the item is not listed. |
| `func AllPayUtility(valuation, ownBid float64, out AuctionOutcome, bidder int) float64` | AllPayUtility returns a bidder's payoff in an all-pay auction: their valuation minus their own bid if they won, otherwise minus their bid. |
| `func AllocativeEfficiency(valuations map[int]float64, out AuctionOutcome) float64` | AllocativeEfficiency returns the ratio of the realised valuation to the maximum possible valuation for a single-item auction: 1 when the item is… |
| `func BidderUtility(valuation float64, out AuctionOutcome, bidder int) float64` | BidderUtility returns the payoff of the winner of a single-item auction given their private valuation and the outcome: valuation minus price if they… |
| `func BlockingPairs(menPrefs, womenPrefs [][]int, m Matching) [][2]int` | BlockingPairs returns every (man, woman) pair that blocks the matching: both strictly prefer each other to their assigned partner. |
| `func CountBlockingPairs(menPrefs, womenPrefs [][]int, m Matching) int` | CountBlockingPairs returns the number of blocking pairs of the matching. |
| `func Dominates(a, b Point) bool` | Dominates reports whether a weakly dominates b (a.X >= b.X and a.Y >= b.Y) with at least one strict inequality. |
| `func FeasibleContains(pts []Point, p Point) bool` | FeasibleContains reports whether the point p lies inside (or on the boundary of) the convex hull of the feasible set pts. |
| `func IsEfficientAllocation(valuations map[int]float64, out AuctionOutcome) bool` | IsEfficientAllocation reports whether the auction outcome awards the item to a bidder holding the maximum valuation among the supplied valuations,… |
| `func IsParetoEfficient(pts []Point, p Point) bool` | IsParetoEfficient reports whether p is not strictly dominated by any point in the feasible set. |
| `func IsPerfectMatching(m Matching) bool` | IsPerfectMatching reports whether every man and woman is matched. |
| `func IsStableMatching(menPrefs, womenPrefs [][]int, m Matching) bool` | IsStableMatching reports whether the matching has no blocking pair. |
| `func NashProduct(p, d Point) float64` | NashProduct returns the Nash product (p.X - d.X)(p.Y - d.Y) relative to the disagreement point d. |
| `func RandomSerialDictatorship(prefs [][]int, seed int64) ([]int, []int, error)` | RandomSerialDictatorship runs serial dictatorship with a uniformly random priority order drawn from a math/rand source seeded by seed, so the outcome… |
| `func SerialDictatorship(prefs [][]int, order []int) ([]int, error)` | SerialDictatorship allocates houses by letting agents choose in the given priority order, each taking their most-preferred still-available house. |
| `func SocialWelfare(bids []CombinatorialBid, alloc Allocation) float64` | SocialWelfare returns the total value of an allocation over the given bids. |
| `func TopTradingCycles(prefs [][]int) ([]int, error)` | TopTradingCycles runs Gale's top-trading-cycles algorithm for a housing market: agent i initially owns house i and prefs`i` ranks all houses. |
| `func TotalBidValue(bids []Bid) float64` | TotalBidValue returns the sum of all bid values. |
| `func TotalSurplus(winners []int, valuations map[int]float64) float64` | TotalSurplus returns the total value realised by an allocation given the winners' valuations indexed by bidder. |

### Types

| Type | What it is |
| --- | --- |
| `Allocation` | Allocation is the outcome of a combinatorial winner-determination problem: the indices (into the original bid slice) of the accepted bids and the… |
| `AuctionOutcome` | AuctionOutcome is the result of a single-item sealed-bid auction: the winning bidder, the price they pay and the seller's revenue. |
| `Bid` | Bid is a single sealed bid: the identity of the bidder and the value they bid. |
| `Coalition` | Coalition is a subset of players encoded as a bitmask: bit i is set when player i belongs to the coalition. |
| `CombinatorialBid` | CombinatorialBid is a bid in a combinatorial auction: bidder Bidder offers Value for the exact bundle of items Items (item indices). |
| `CoopGame` | CoopGame is a transferable-utility cooperative game on the fixed player set {0,...,Players-1}. |
| `Matching` | Matching is a two-sided assignment between men and women (or, generally, two equal-sized sides). |
| `MultiUnitOutcome` | MultiUnitOutcome is the result of a multi-unit auction in which every winning bidder demands a single identical unit. |
| `Point` | Point is a two-dimensional payoff pair (utilities of players 0 and 1). |
| `VCGResult` | VCGResult is the outcome of running the VCG mechanism on a combinatorial auction: the welfare-maximizing allocation, and for every bidder that wins… |
| `WeightedVotingGame` | WeightedVotingGame is a simple game [q; w_0,...,w_{n-1}]: a coalition wins when the total weight of its members reaches the quota Quota. |

<details>
<summary><code>Allocation</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func WinnerDetermination(bids []CombinatorialBid, numItems int) (Allocation, error)` | WinnerDetermination solves the combinatorial-auction winner-determination problem under XOR bidding (each bidder wins at most one bid): it selects a… |
| `func WinnerDeterminationOR(bids []CombinatorialBid, numItems int) (Allocation, error)` | WinnerDeterminationOR solves the winner-determination problem under OR bidding, where a bidder may be awarded several of their disjoint bids. |

</details>

<details>
<summary><code>AuctionOutcome</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AllPayAuction(bids []Bid, reserve float64) (AuctionOutcome, error)` | AllPayAuction runs an all-pay auction: the highest bidder wins the item but every bidder pays their own bid, so the seller's revenue is the sum of… |
| `func DutchAuction(bids []Bid, reserve float64) (AuctionOutcome, error)` | DutchAuction runs a descending-clock (Dutch) auction with an optional reserve. |
| `func EnglishAuction(bids []Bid, reserve float64) (AuctionOutcome, error)` | EnglishAuction runs an ascending open-outcry (English) auction with an optional reserve. |
| `func FirstPriceAuction(bids []Bid, reserve float64) (AuctionOutcome, error)` | FirstPriceAuction runs a first-price sealed-bid auction with an optional reserve price: the highest bidder wins and pays their own bid, provided it… |
| `func SecondPriceAuction(bids []Bid, reserve float64) (AuctionOutcome, error)` | SecondPriceAuction runs a second-price (Vickrey) sealed-bid auction with an optional reserve: the highest bidder wins and pays the larger of the… |
| `func ThirdPriceAuction(bids []Bid, reserve float64) (AuctionOutcome, error)` | ThirdPriceAuction runs a third-price sealed-bid auction: the highest bidder wins and pays the third-highest bid (or the reserve, whichever is larger). |
| `func VickreyAuction(bids []Bid, reserve float64) (AuctionOutcome, error)` | VickreyAuction is an alias for SecondPriceAuction, the single-item Vickrey mechanism. |

</details>

<details>
<summary><code>Bid</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func HighestBid(bids []Bid) (Bid, error)` | HighestBid returns the largest bid value and its bidder. |
| `func SecondHighestBid(bids []Bid) (Bid, error)` | SecondHighestBid returns the second-largest bid value and its bidder. |
| `func SortBidsDescending(bids []Bid) []Bid` | SortBidsDescending returns a copy of the bids sorted by decreasing value, ties broken by increasing bidder index. |

</details>

<details>
<summary><code>Coalition</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AllCoalitions(n int) []Coalition` | AllCoalitions returns every one of the 2^n coalitions over n players, in increasing bitmask order (starting with the empty coalition). |
| `func CoalitionFromMembers(members []int) Coalition` | CoalitionFromMembers builds a coalition from an explicit list of player indices. |
| `func EmptyCoalition() Coalition` | EmptyCoalition returns the empty coalition. |
| `func FullCoalition(n int) Coalition` | FullCoalition returns the grand coalition {0,...,n-1}. |
| `func ProperSubsets(c Coalition, n int) []Coalition` | ProperSubsets returns every non-empty proper subset of the coalition c over a universe of n players (excluding the empty set and c itself), sorted by… |
| `func SingletonCoalition(i int) Coalition` | SingletonCoalition returns the coalition {i}. |
| `func (c Coalition) Add(i int) Coalition` | Add returns the coalition with player i added. |
| `func (c Coalition) Complement(n int) Coalition` | Complement returns the coalition of all players in {0,...,n-1} not in c. |
| `func (c Coalition) Contains(i int) bool` | Contains reports whether player i belongs to the coalition. |
| `func (c Coalition) Difference(other Coalition) Coalition` | Difference returns the players in c that are not in other. |
| `func (c Coalition) Intersect(other Coalition) Coalition` | Intersect returns the coalition of players in both c and other. |
| `func (c Coalition) IsEmpty() bool` | IsEmpty reports whether the coalition has no members. |
| `func (c Coalition) IsSubsetOf(other Coalition) bool` | IsSubsetOf reports whether every member of c is also a member of other. |
| `func (c Coalition) Members(n int) []int` | Members returns the sorted player indices contained in the coalition, given a universe of n players. |
| `func (c Coalition) Remove(i int) Coalition` | Remove returns the coalition with player i removed. |
| `func (c Coalition) Size() int` | Size returns the number of players in the coalition. |
| `func (c Coalition) Union(other Coalition) Coalition` | Union returns the coalition containing every player in c or in other. |

</details>

<details>
<summary><code>CombinatorialBid</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func SingleItemToCombinatorial(bids []Bid) []CombinatorialBid` | SingleItemToCombinatorial converts single-item sealed bids into combinatorial bids for the lone item 0, so that VCGMechanism reproduces the… |

</details>

<details>
<summary><code>CoopGame</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewCoopGame(players int, worth func(Coalition) float64) (CoopGame, error)` | NewCoopGame builds a CoopGame by evaluating the characteristic function on every coalition bitmask. |
| `func NewCoopGameValues(players int, values []float64) (CoopGame, error)` | NewCoopGameValues builds a CoopGame from a precomputed table of coalition worths. |
| `func (g CoopGame) Add(other CoopGame) (CoopGame, error)` | Add returns the game v + w, the coalition-wise sum of two games on the same player set. |
| `func (g CoopGame) AreSymmetric(i, j int) bool` | AreSymmetric reports whether players i and j are interchangeable: swapping them leaves the worth of every coalition unchanged. |
| `func (g CoopGame) BanzhafPowerIndex() []float64` | BanzhafPowerIndex returns the normalized Banzhaf power index: the raw Banzhaf values rescaled to sum to one. |
| `func (g CoopGame) BanzhafValue() []float64` | BanzhafValue returns the (raw, probabilistic) Banzhaf value: each player's average marginal contribution over the 2^(n-1) coalitions not containing… |
| `func (g CoopGame) CoreConstraintViolations(x []float64, tol float64) []Coalition` | CoreConstraintViolations returns the coalitions whose core inequality x(S) >= v(S) is violated by more than tol at the allocation x. |
| `func (g CoopGame) CoreIsNonEmpty() bool` | CoreIsNonEmpty reports whether the core of the game is non-empty. |
| `func (g CoopGame) DualGame() CoopGame` | DualGame returns the dual game v*(S) = v(N) - v(N\S). |
| `func (g CoopGame) EqualSplit() []float64` | EqualSplit returns the egalitarian allocation that divides v(N) equally among the players. |
| `func (g CoopGame) Excess(s Coalition, x []float64) float64` | Excess returns the excess e(S, x) = v(S) - Σ_{i∈S} x_i of coalition S at the allocation x. |
| `func (g CoopGame) GrandCoalition() Coalition` | GrandCoalition returns the coalition of all players. |
| `func (g CoopGame) GrandValue() float64` | GrandValue returns v(N), the worth of the grand coalition. |
| `func (g CoopGame) HarsanyiDividends() []float64` | HarsanyiDividends returns the Harsanyi dividends d(S) of the game, the Möbius transform of the characteristic function: d(S) = Σ_{T⊆S}… |
| `func (g CoopGame) InCore(x []float64, tol float64) bool` | InCore reports whether the allocation x lies in the core: efficient and coalitionally rational (Σ_{i∈S} x_i >= v(S) for every coalition S), within… |
| `func (g CoopGame) InEpsilonCore(x []float64, eps, tol float64) bool` | InEpsilonCore reports whether x lies in the (strong) ε-core: efficient with every coalition's excess at most eps (within tolerance tol). |
| `func (g CoopGame) IndividualRational(x []float64, tol float64) bool` | IndividualRational reports whether every player receives at least their singleton worth under x (within tolerance tol). |
| `func (g CoopGame) IsAdditive() bool` | IsAdditive reports whether the game is additive (inessential): v(S) = Σ_{i∈S} v({i}) for every coalition S. |
| `func (g CoopGame) IsBalanced() bool` | IsBalanced reports whether the game is balanced, which by the Bondareva-Shapley theorem is equivalent to having a non-empty core. |
| `func (g CoopGame) IsConstantSum() bool` | IsConstantSum reports whether v(S) + v(N\S) = v(N) for every coalition S. |
| `func (g CoopGame) IsConvex() bool` | IsConvex reports whether the game is convex (supermodular): v(S∪T) + v(S∩T) >= v(S) + v(T) for all coalitions S and T. Convex games have a… |
| `func (g CoopGame) IsDummy(i int) bool` | IsDummy reports whether player i is a dummy: its marginal contribution equals its singleton worth for every coalition, so it adds nothing beyond… |
| `func (g CoopGame) IsEfficient(x []float64, tol float64) bool` | IsEfficient reports whether the allocation x distributes exactly v(N) among the players, within the tolerance tol. |
| `func (g CoopGame) IsEssential() bool` | IsEssential reports whether v(N) > Σ_i v({i}); an essential game admits a non-degenerate set of imputations. |
| `func (g CoopGame) IsImputation(x []float64, tol float64) bool` | IsImputation reports whether x is an imputation: efficient and individually rational (x_i >= v({i}) for every player i), within tolerance tol. |
| `func (g CoopGame) IsMonotone() bool` | IsMonotone reports whether v(S) <= v(T) whenever S ⊆ T. |
| `func (g CoopGame) IsNullPlayer(i int) bool` | IsNullPlayer reports whether player i is a null player: its marginal contribution to every coalition is zero. |
| `func (g CoopGame) IsSimple() bool` | IsSimple reports whether the game is a simple game: monotone with worths in {0,1} and v(N)=1. |
| `func (g CoopGame) IsSuperadditive() bool` | IsSuperadditive reports whether v(S∪T) >= v(S) + v(T) for all disjoint coalitions S and T. |
| `func (g CoopGame) IsZeroNormalized() bool` | IsZeroNormalized reports whether every singleton coalition is worth zero. |
| `func (g CoopGame) LeastCore() (float64, []float64, bool)` | LeastCore returns the least-core value ε* together with an allocation that attains it: the allocation minimizing the maximum excess, subject to… |
| `func (g CoopGame) LeastCoreValue() (float64, bool)` | LeastCoreValue returns ε*, the least-core value: the minimum, over all efficient allocations, of the maximum coalition excess. |
| `func (g CoopGame) MarginalContribution(i int, s Coalition) float64` | MarginalContribution returns v(S ∪ {i}) - v(S \ {i}), the value player i adds to the coalition S. S need not exclude i: player i is removed from… |
| `func (g CoopGame) MaxExcess(x []float64) float64` | MaxExcess returns the largest excess over all proper non-empty coalitions at the allocation x. |
| `func (g CoopGame) MinimalRights() []float64` | MinimalRights returns the minimal-rights (concession) vector m used by the tau-value: m_i = max over coalitions S containing i of v(S) - Σ_{j∈S,… |
| `func (g CoopGame) Nucleolus() ([]float64, error)` | Nucleolus returns the nucleolus of the game: the unique imputation that lexicographically minimizes the sorted (descending) vector of coalition… |
| `func (g CoopGame) Prenucleolus() ([]float64, error)` | Prenucleolus returns the prenucleolus: the lexicographic minimizer of the sorted excess vector over all efficient allocations (without the… |
| `func (g CoopGame) Scale(c float64) CoopGame` | Scale returns the game c·v obtained by multiplying every coalition worth by c. |
| `func (g CoopGame) ShapleyFromDividends() []float64` | ShapleyFromDividends returns the Shapley value computed from the Harsanyi dividends, φ_i = Σ_{S ∋ i} d(S)/\|S\|. |
| `func (g CoopGame) ShapleyShubikIndex() []float64` | ShapleyShubikIndex returns the Shapley-Shubik power index of a simple game, which equals its Shapley value. |
| `func (g CoopGame) ShapleyValue() []float64` | ShapleyValue returns the Shapley value φ, the unique efficient allocation that averages each player's marginal contribution over all coalition… |
| `func (g CoopGame) ShapleyValueMonteCarlo(samples int, seed int64) []float64` | ShapleyValueMonteCarlo estimates the Shapley value by sampling random player orderings, using a math/rand source seeded by seed so that the estimate… |
| `func (g CoopGame) SortedExcesses(x []float64) []float64` | SortedExcesses returns the excesses of every proper non-empty coalition at the allocation x, sorted in decreasing order. |
| `func (g CoopGame) TauValue() ([]float64, error)` | TauValue returns the tau-value (compromise value) of the game: the efficient allocation on the segment between the minimal-rights vector m and the… |
| `func (g CoopGame) UtopiaPayoff() []float64` | UtopiaPayoff returns the utopia (marginal) payoff vector M with M_i = v(N) - v(N\{i}), the most player i could hope to receive if every other player… |
| `func (g CoopGame) Value(s Coalition) float64` | Value returns the worth v(S) of the coalition S. |
| `func (g CoopGame) ZeroNormalized() CoopGame` | ZeroNormalized returns the strategically equivalent zero-normalized game obtained by subtracting each player's singleton worth: w(S) = v(S) -… |

</details>

<details>
<summary><code>Matching</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func GaleShapley(menPrefs, womenPrefs [][]int) (Matching, error)` | GaleShapley computes the man-optimal stable matching by the Gale-Shapley deferred-acceptance algorithm with men proposing. |
| `func ManOptimalMatching(menPrefs, womenPrefs [][]int) (Matching, error)` | ManOptimalMatching is Gale-Shapley with men proposing; it returns the stable matching most preferred by every man. |
| `func WomanOptimalMatching(menPrefs, womenPrefs [][]int) (Matching, error)` | WomanOptimalMatching is Gale-Shapley with women proposing; it returns the stable matching most preferred by every woman. |

</details>

<details>
<summary><code>MultiUnitOutcome</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func MultiUnitVickrey(bids []Bid, units int) (MultiUnitOutcome, error)` | MultiUnitVickrey runs the multi-unit Vickrey (VCG) auction under single-unit demand: the units highest bidders each win one unit, and each winner… |
| `func PayAsBidAuction(bids []Bid, units int) (MultiUnitOutcome, error)` | PayAsBidAuction (discriminatory auction) sells units identical goods under single-unit demand: the units highest bidders each win one unit and pay… |
| `func UniformPriceAuction(bids []Bid, units int) (MultiUnitOutcome, error)` | UniformPriceAuction sells units identical goods under single-unit demand: the units highest bidders each win one unit and all pay the same clearing… |

</details>

<details>
<summary><code>Point</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ConvexHull2D(pts []Point) []Point` | ConvexHull2D returns the vertices of the convex hull of the input points in counter-clockwise order using Andrew's monotone-chain algorithm. |
| `func EgalitarianSolution(pts []Point, d Point) (Point, error)` | EgalitarianSolution returns the egalitarian bargaining solution: the Pareto-efficient point giving both players an equal gain over the disagreement… |
| `func IdealPoint(pts []Point) (Point, error)` | IdealPoint returns the utopia point of the feasible set: the maximum achievable X paired with the maximum achievable Y (generally infeasible… |
| `func KalaiSmorodinskySolution(pts []Point, d Point) (Point, error)` | KalaiSmorodinskySolution returns the Kalai-Smorodinsky bargaining solution: the Pareto-efficient point on the segment from the disagreement point d… |
| `func NashBargainingDiscrete(pts []Point, d Point) (Point, error)` | NashBargainingDiscrete returns the Nash bargaining solution restricted to the finite set of feasible points (no convexification): the point… |
| `func NashBargainingSolution(pts []Point, d Point) (Point, error)` | NashBargainingSolution returns the Nash bargaining solution over the convex hull of the feasible payoff set: the point maximizing the product of… |
| `func ParetoFrontier2D(pts []Point) []Point` | ParetoFrontier2D returns the Pareto-efficient points of the input set — those not dominated by any other point — sorted by increasing X. |
| `func UtilitarianSolution(pts []Point) (Point, error)` | UtilitarianSolution returns a feasible point maximizing the sum of utilities X + Y. Because the objective is linear, an optimum is attained at a hull… |
| `func WeightedNashBargainingSolution(pts []Point, d Point, alpha float64) (Point, error)` | WeightedNashBargainingSolution returns the asymmetric Nash bargaining solution with bargaining power alpha in (0,1) for player 0 and 1-alpha for… |

</details>

<details>
<summary><code>VCGResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func VCGMechanism(bids []CombinatorialBid, numItems int) (VCGResult, error)` | VCGMechanism runs the Vickrey-Clarke-Groves mechanism under XOR bidding: it computes the efficient allocation and charges each winner the externality… |

</details>

<details>
<summary><code>WeightedVotingGame</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewWeightedVotingGame(quota float64, weights []float64) (WeightedVotingGame, error)` | NewWeightedVotingGame builds a weighted voting game from a quota and weights. |
| `func (w WeightedVotingGame) BanzhafIndex() []float64` | BanzhafIndex returns the normalized Banzhaf power index of each player: their swing count divided by the total number of swings across all players. |
| `func (w WeightedVotingGame) CoalitionWeight(s Coalition) float64` | CoalitionWeight returns the total weight of the members of coalition s. |
| `func (w WeightedVotingGame) CoopGame() CoopGame` | CoopGame converts the weighted voting game into the corresponding simple cooperative game with worth 1 for winning coalitions and 0 otherwise. |
| `func (w WeightedVotingGame) Dictator() (int, bool)` | Dictator returns the index of a dictator (a player who wins alone and whose absence makes every coalition lose) and true, or -1 and false if there is… |
| `func (w WeightedVotingGame) Dummies() []int` | Dummies returns the players who are never critical in any coalition. |
| `func (w WeightedVotingGame) IsCritical(i int, s Coalition) bool` | IsCritical reports whether player i is critical (a swing) in coalition s: s is winning but s without i is losing. |
| `func (w WeightedVotingGame) IsWinning(s Coalition) bool` | IsWinning reports whether coalition s meets or exceeds the quota. |
| `func (w WeightedVotingGame) MinimalWinningCoalitions() []Coalition` | MinimalWinningCoalitions returns the winning coalitions that become losing if any single member leaves. |
| `func (w WeightedVotingGame) Players() int` | Players returns the number of players in the game. |
| `func (w WeightedVotingGame) ShapleyShubikIndex() []float64` | ShapleyShubikIndex returns the Shapley-Shubik power index of each player: the fraction of the n! orderings in which the player is pivotal (turns the… |
| `func (w WeightedVotingGame) SwingCount(i int) int` | SwingCount returns the number of coalitions in which player i is critical. |
| `func (w WeightedVotingGame) VetoPlayers() []int` | VetoPlayers returns the players present in every winning coalition; such a player can unilaterally block any decision. |

</details>

### Variables

`ErrNoBids`, `ErrPlayers`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/auctions).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
