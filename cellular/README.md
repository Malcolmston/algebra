# cellular — implements cellular automata in pure Go using only the standard library

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/cellular.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/cellular)

Package cellular implements cellular automata in pure Go using only the
standard library.

- Elementary one-dimensional automata (`ElementaryRule`) with the standard
Wolfram numbering 0–255, single-step and spacetime evolution, the four
symmetry conjugates (mirror, complement and mirror-complement), the 88
inequivalent equivalence classes, additivity and totalistic tests, the Langton
lambda parameter and a Wolfram-class heuristic (`ElementaryRule.Class`).

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/cellular"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func BinaryEntropy(p float64) float64` | BinaryEntropy returns the Shannon entropy in bits of a Bernoulli(p) distribution, -p*log2(p) - (1-p)*log2(1-p). |
| `func BlockEntropy(state []int, blockLen, k int) float64` | BlockEntropy returns the Shannon entropy in bits of the distribution of length blockLen blocks read with wraparound across state over the alphabet… |
| `func BoundingBox(g *Grid) (minRow, minCol, maxRow, maxCol int, ok bool)` | BoundingBox returns the smallest rectangle containing every live cell as (minRow, minCol, maxRow, maxCol). |
| `func ClassName(class int) string` | ClassName returns the textual name of a Wolfram class number. |
| `func CloneState(state []int) []int` | CloneState returns an independent copy of state. |
| `func CountNeighbours(g *Grid, r, c int, nb Neighbourhood, bc Boundary) int` | CountNeighbours returns the number of live cells adjacent to (r, c) under the given neighbourhood and boundary condition. |
| `func CountOnes(rows [][]int) int` | CountOnes returns the total number of non-zero cells across a spacetime diagram. |
| `func DamageSpread(rule Rule1D, initial []int, flipAt, steps int, bc Boundary) []int` | DamageSpread evolves two configurations that differ by a single flipped cell and returns, for each of the steps time steps, the Hamming distance… |
| `func DecodeBaseK(code, digits, k int) ([]int, error)` | DecodeBaseK returns the digits digits-long base-k representation of code, most significant digit first. |
| `func Density(state []int) float64` | Density returns the fraction of non-zero cells in state, or 0 for an empty slice. |
| `func DiagramDimensions(rows [][]int) (height, width int)` | DiagramDimensions returns the number of rows and the (first-row) width of a diagram. |
| `func DiagramString(rows [][]int) string` | DiagramString is a convenience that renders a binary diagram with '.' for 0 and '#' for 1. |
| `func ElementaryEvolve(r ElementaryRule, initial []int, steps int, bc Boundary) [][]int` | ElementaryEvolve returns the steps+1 row spacetime diagram of the elementary rule starting from initial. |
| `func ElementaryEvolveSeed(r ElementaryRule, width, steps int) [][]int` | ElementaryEvolveSeed evolves the rule for steps steps from a single central seed on a background of width cells, using a fixed-zero boundary so the… |
| `func ElementaryStep(r ElementaryRule, state []int, bc Boundary) []int` | ElementaryStep advances a binary configuration by one time step under the elementary rule and boundary condition, a thin wrapper over Step1D. |
| `func EncodeBaseK(vals []int, k int) (int, error)` | EncodeBaseK interprets vals as the digits of a base-k number, most-significant digit first, and returns its integer value. |
| `func EqualState(a, b []int) bool` | EqualState reports whether a and b have identical length and contents. |
| `func Evolve1D(rule Rule1D, initial []int, steps int, bc Boundary) [][]int` | Evolve1D returns the spacetime diagram produced by iterating rule for steps time steps starting from initial. |
| `func FinalRow(rows [][]int) []int` | FinalRow returns the last row of a spacetime diagram, or nil if it is empty. |
| `func FindPeriod(l *LifeRule, g *Grid, nb Neighbourhood, bc Boundary, maxSteps int) (period, preperiod int)` | FindPeriod runs a life-like rule from g under a fixed boundary and reports the smallest period p in 1..maxPeriod for which the configuration recurs,… |
| `func FormatRuleSummary(r ElementaryRule) string` | FormatRuleSummary returns a one-line summary of an elementary rule combining its number, symmetry conjugates, lambda parameter and heuristic class. |
| `func GridPopulations(frames []*Grid) []int` | GridPopulations returns the population of each grid in a sequence of frames. |
| `func HammingDistance(a, b []int) (int, error)` | HammingDistance returns the number of positions at which a and b differ. |
| `func IterateState(rule Rule1D, initial []int, steps int, bc Boundary) []int` | IterateState applies rule steps times and returns only the final configuration, allocating one working buffer instead of a full diagram. |
| `func LangtonLambda(rule *GeneralRule) float64` | LangtonLambda returns Langton's lambda parameter for a general rule: the fraction of neighbourhoods that map to a non-quiescent state, where state 0… |
| `func LightConeWidth(rows [][]int) []int` | LightConeWidth returns, for each row of a diagram grown from a central seed, the distance from the centre to the outermost non-zero cell, i.e. |
| `func LiveCells(g *Grid) [][2]int` | LiveCells returns the sorted coordinates of the live cells of a grid as [row, col] pairs, ordered by row then column. |
| `func MeanSpatialEntropy(rows [][]int) float64` | MeanSpatialEntropy returns the average SpatialEntropy over the rows of a spacetime diagram. |
| `func NeighborhoodCount(k, r int) *big.Int` | NeighborhoodCount returns the number of distinct neighbourhoods for a k-state, radius-r automaton, namely k^(2r+1). |
| `func NumOuterTotalisticRules(k, r int) *big.Int` | NumOuterTotalisticRules returns the number of distinct k-state, radius-r outer-totalistic rules, k^(k*(2r(k-1)+1)). |
| `func NumRules(k, r int) *big.Int` | NumRules returns the number of distinct k-state, radius-r rules as an exact big integer, k^(k^(2r+1)). |
| `func NumTotalisticRules(k, r int) *big.Int` | NumTotalisticRules returns the number of distinct k-state, radius-r totalistic rules, k^((2r+1)(k-1)+1). |
| `func ParseDiagram(s string, on rune) [][]int` | ParseDiagram parses a multi-line ASCII binary diagram (as produced by DiagramString) back into a spacetime diagram, mapping on to 1 and every other… |
| `func PopCountK(vals []int) int` | PopCountK returns the number of non-zero digits among vals. |
| `func Population(state []int) int` | Population returns the number of non-zero cells in state. |
| `func RandomBinaryState(n int, seed uint64) []int` | RandomBinaryState returns a reproducible length-n binary configuration, a convenience wrapper around RandomState with k = 2. |
| `func RandomState(n, k int, seed uint64) ([]int, error)` | RandomState returns a reproducible length-n configuration whose cells are uniformly distributed over 0..k-1, generated by a deterministic splitmix64… |
| `func RenderASCII(rows [][]int, off, on rune) string` | RenderASCII renders a binary spacetime diagram as a multi-line string, using on for non-zero cells and off for zero cells, one row per line. |
| `func RenderStates(rows [][]int, glyphs []rune) string` | RenderStates renders a spacetime diagram over an arbitrary alphabet using the runes in glyphs indexed by cell value. |
| `func RowDensities(rows [][]int) []float64` | RowDensities returns the density of each row of a spacetime diagram. |
| `func RowPopulations(rows [][]int) []int` | RowPopulations returns the population (non-zero count) of each row of a spacetime diagram. |
| `func RuleTableString(r ElementaryRule) string` | RuleTableString returns a multi-line human-readable transition table for an elementary rule, listing each neighbourhood pattern and its output. |
| `func RuleToBits(r ElementaryRule) []int` | RuleToBits returns the 8 output bits of an elementary rule, indexed by neighbourhood value 0..7 (identical to Table but returned as a slice). |
| `func SecondOrderEvolve(r ElementaryRule, prev, cur []int, steps int, bc Boundary) [][]int` | SecondOrderEvolve runs the second-order construction of a rule for steps generations from the two initial configurations and returns the spacetime… |
| `func SecondOrderReverse(r ElementaryRule, next, cur []int, bc Boundary) []int` | SecondOrderReverse computes one backward step of the second-order construction: given the current and next configurations it recovers the previous… |
| `func SecondOrderStep(r ElementaryRule, prev, cur []int, bc Boundary) []int` | SecondOrderStep computes one forward step of the second-order (Fredkin) construction of an elementary rule: the next configuration is the rule… |
| `func ShannonEntropyBits(probs []float64) float64` | ShannonEntropyBits returns the Shannon entropy in bits of a discrete probability distribution. |
| `func SingleSeedState(n int) []int` | SingleSeedState returns a length-n configuration of zeros with a single 1 in the middle cell (index n/2). |
| `func SingleSeedStateAt(n, pos int) ([]int, error)` | SingleSeedStateAt returns a length-n configuration of zeros with a single 1 at index pos. |
| `func SpacetimeDensity(rows [][]int) float64` | SpacetimeDensity returns the mean fraction of non-zero cells across every row of a spacetime diagram. |
| `func SpatialEntropy(state []int) float64` | SpatialEntropy returns the binary symbol entropy of state (equivalent to SymbolEntropy with k = 2 for a binary configuration). |
| `func Stamp(g, pattern *Grid, top, left int)` | Stamp copies the live cells of pattern into g with its top-left corner at (top, left), leaving the rest of g unchanged. |
| `func StateFromString(s string, on rune) []int` | StateFromString parses a string into a binary configuration, mapping each rune equal to on to 1 and every other rune to 0. |
| `func StateToString(state []int, off, on rune) string` | StateToString renders a configuration as a string using off for 0 and on for every non-zero value. |
| `func Step1D(rule Rule1D, state []int, bc Boundary) []int` | Step1D advances a one-dimensional configuration by one time step under rule and the given boundary condition, returning a freshly allocated slice of… |
| `func Sum(state []int) int` | Sum returns the arithmetic sum of the cell values in state. |
| `func SymbolEntropy(state []int, k int) float64` | SymbolEntropy returns the empirical Shannon entropy in bits of the symbol distribution of state over the alphabet 0..k-1. |
| `func TemporalPeriod(rows [][]int, maxPeriod int) int` | TemporalPeriod scans a spacetime diagram from the last row backwards and returns the smallest positive p (up to maxPeriod) such that the final row… |
| `func TotalisticRuleString(t *TotalisticRule) string` | TotalisticRuleString renders the sum-to-output table of a totalistic rule as a single line such as "0->0 1->1 2->0 3->1". |
| `func TransposeRows(rows [][]int) [][]int` | TransposeRows returns the transpose of a rectangular spacetime diagram, so that column j of the input becomes row j of the output. |
| `func Uniform(state []int) bool` | Uniform reports whether every cell of state holds the same value. |
| `func ValidateState(state []int, k int) error` | ValidateState reports an error if any cell of state lies outside 0..k-1. |

### Types

| Type | What it is |
| --- | --- |
| `Boundary` | Boundary selects how cells beyond the ends of a one-dimensional lattice (or the borders of a two-dimensional grid) are read when applying a rule. |
| `ElementaryRule` | ElementaryRule is an elementary (k=2, r=1) one-dimensional cellular-automaton rule identified by its Wolfram number 0..255. |
| `GeneralRule` | GeneralRule is a fully general k-state, radius-r rule specified by a complete lookup table over all k^(2r+1) neighbourhoods. |
| `Grid` | Grid is a rectangular two-dimensional binary configuration stored in row-major order. |
| `LifeRule` | LifeRule is a life-like (outer-totalistic binary two-dimensional) rule given by its birth and survival neighbour counts. |
| `MargolusRule` | MargolusRule is a reversible two-dimensional block rule acting on 2x2 Margolus-neighbourhood blocks of a binary grid. |
| `Neighbourhood` | Neighbourhood selects the cell adjacency used by a life-like rule. |
| `OuterTotalisticRule` | OuterTotalisticRule is a k-state, radius-r outer-totalistic rule: the new value of a cell depends on the centre cell together with the sum of the… |
| `Rule1D` | Rule1D is the common interface implemented by every one-dimensional cellular automaton rule in the package. |
| `SecondOrderCA` | SecondOrderCA holds the state of a running second-order elementary automaton: the underlying rule, the previous and current generations and the… |
| `TotalisticRule` | TotalisticRule is a k-state, radius-r totalistic one-dimensional rule: the new value of a cell depends only on the sum of the values in its… |

<details>
<summary><code>Boundary</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (b Boundary) String() string` | String returns a human-readable name for the boundary condition. |

</details>

<details>
<summary><code>ElementaryRule</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AdditiveElementaryRules() []ElementaryRule` | AdditiveElementaryRules returns the eight additive (GF(2)-linear) elementary rules in increasing order. |
| `func AllElementaryRules() []ElementaryRule` | AllElementaryRules returns all 256 elementary rules in numeric order. |
| `func InequivalentElementaryRules() []ElementaryRule` | InequivalentElementaryRules returns the 88 canonical elementary rules, one per equivalence class under reflection and complementation, in numeric… |
| `func MustElementaryRule(n int) ElementaryRule` | MustElementaryRule is like NewElementaryRule but panics on an out-of-range number. |
| `func NewElementaryRule(n int) (ElementaryRule, error)` | NewElementaryRule returns the ElementaryRule with Wolfram number n. |
| `func Rule102() ElementaryRule` | Rule102 returns elementary rule 102, an additive rule mirror-equivalent to rule 60. |
| `func Rule110() ElementaryRule` | Rule110 returns elementary rule 110, proved capable of universal computation. |
| `func Rule126() ElementaryRule` | Rule126 returns elementary rule 126, a chaotic rule. |
| `func Rule150() ElementaryRule` | Rule150 returns elementary rule 150, the additive sum-mod-2 rule of all three neighbours. |
| `func Rule18() ElementaryRule` | Rule18 returns elementary rule 18, a chaotic rule generating a Sierpinski-like fractal. |
| `func Rule184() ElementaryRule` | Rule184 returns elementary rule 184, the traffic/ballistic-annihilation rule. |
| `func Rule22() ElementaryRule` | Rule22 returns elementary rule 22, a chaotic outer-totalistic rule. |
| `func Rule250() ElementaryRule` | Rule250 returns elementary rule 250, which fills space from a single seed. |
| `func Rule30() ElementaryRule` | Rule30 returns elementary rule 30, a chaotic, additive-looking rule used as a random-number source. |
| `func Rule4() ElementaryRule` | Rule4 returns elementary rule 4, a class-2 rule whose isolated cells persist. |
| `func Rule45() ElementaryRule` | Rule45 returns elementary rule 45, a chaotic rule often used for random generation. |
| `func Rule54() ElementaryRule` | Rule54 returns elementary rule 54, a class-4 candidate with particle-like structures. |
| `func Rule60() ElementaryRule` | Rule60 returns elementary rule 60, an additive (Pascal-triangle mod 2) rule. |
| `func Rule90() ElementaryRule` | Rule90 returns elementary rule 90, the additive XOR rule that draws the Sierpinski triangle from a single seed. |
| `func RuleFromBits(bitsTable []int) (ElementaryRule, error)` | RuleFromBits builds an elementary rule from an 8-entry output table indexed by neighbourhood value 0..7. |
| `func SierpinskiRule() ElementaryRule` | SierpinskiRule is an alias for Rule90, whose single-seed evolution is the Sierpinski triangle. |
| `func TotalisticElementaryRules() []ElementaryRule` | TotalisticElementaryRules returns the elementary rules that are totalistic in the strict sense (output depends only on the neighbourhood sum). |
| `func TrafficRule() ElementaryRule` | TrafficRule is an alias for Rule184, the elementary traffic model. |
| `func (r ElementaryRule) ActiveCount() int` | ActiveCount returns the number of neighbourhoods that map to 1 (the population count of the rule number). |
| `func (r ElementaryRule) Apply(neighbourhood []int) int` | Apply implements the Rule1D interface for a length-3 neighbourhood slice (left, centre, right). |
| `func (r ElementaryRule) ApplyLCR(left, centre, right int) int` | ApplyLCR returns the next centre value for the ordered neighbourhood (left, centre, right), each 0 or 1. |
| `func (r ElementaryRule) Bit(p int) int` | Bit returns the output bit (0 or 1) for the neighbourhood value p in 0..7, where p = 4*left + 2*centre + right. |
| `func (r ElementaryRule) Class() int` | Class returns a heuristic Wolfram class (1, 2 or 3) for an elementary rule, obtained by evolving a fixed reproducible pseudo-random configuration on… |
| `func (r ElementaryRule) ComplementRule() ElementaryRule` | ComplementRule returns the rule obtained by exchanging the roles of 0 and 1 in both inputs and outputs. |
| `func (r ElementaryRule) EquivalenceClassSize() int` | EquivalenceClassSize returns the number of distinct rules in r's equivalence class (1, 2 or 4). |
| `func (r ElementaryRule) Equivalents() []ElementaryRule` | Equivalents returns the sorted set of rules in the same equivalence class as r under left-right reflection and 0/1 complementation (1, 2 or 4… |
| `func (r ElementaryRule) IsAdditive() bool` | IsAdditive reports whether the rule is additive (linear over GF(2)): its output equals a fixed XOR combination of the three neighbourhood cells. |
| `func (r ElementaryRule) IsCanonical() bool` | IsCanonical reports whether r is the minimal representative of its equivalence class. |
| `func (r ElementaryRule) IsLegal() bool` | IsLegal reports whether the rule is "legal" in Wolfram's original sense: quiescent (blank stays blank) and left-right symmetric. |
| `func (r ElementaryRule) IsOuterTotalistic() bool` | IsOuterTotalistic reports whether the rule depends only on the centre cell and the sum of its two neighbours (left+right), rather than on their… |
| `func (r ElementaryRule) IsQuiescent() bool` | IsQuiescent reports whether the all-zero neighbourhood maps to 0, so that a blank background stays blank. |
| `func (r ElementaryRule) IsSymmetric() bool` | IsSymmetric reports whether the rule is invariant under left-right reflection, i.e. |
| `func (r ElementaryRule) IsTotalistic() bool` | IsTotalistic reports whether the rule depends only on the total number of 1s in the whole three-cell neighbourhood. |
| `func (r ElementaryRule) LambdaParameter() float64` | LambdaParameter returns Langton's lambda for the rule: the fraction of neighbourhoods mapping to a non-quiescent state, i.e. |
| `func (r ElementaryRule) MinimalEquivalent() ElementaryRule` | MinimalEquivalent returns the smallest rule number in r's equivalence class, the canonical representative used to enumerate the 88 inequivalent… |
| `func (r ElementaryRule) MirrorComplementRule() ElementaryRule` | MirrorComplementRule returns the composition of MirrorRule and ComplementRule. |
| `func (r ElementaryRule) MirrorRule() ElementaryRule` | MirrorRule returns the rule obtained by reflecting every neighbourhood left-to-right. |
| `func (r ElementaryRule) Number() int` | Number returns the Wolfram rule number as an int. |
| `func (r ElementaryRule) Radius() int` | Radius reports the neighbourhood radius, which is always 1 for an elementary rule. |
| `func (r ElementaryRule) States() int` | States reports the number of cell states, which is always 2 for an elementary rule. |
| `func (r ElementaryRule) String() string` | String returns the rule in the form "rule 110". |
| `func (r ElementaryRule) Table() [8]int` | Table returns the 8-entry output table indexed by neighbourhood value 0..7. |

</details>

<details>
<summary><code>GeneralRule</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ElementaryAsGeneral(r ElementaryRule) *GeneralRule` | ElementaryAsGeneral converts an elementary rule to the equivalent GeneralRule with k = 2 and r = 1, useful for exercising the general machinery. |
| `func GeneralRuleFromCode(k, r, code int) (*GeneralRule, error)` | GeneralRuleFromCode builds a general rule whose table is the base-k representation of code, least significant digit indexing neighbourhood 0. |
| `func NewGeneralRule(k, r int, table []int) (*GeneralRule, error)` | NewGeneralRule builds a general rule from an explicit table whose length must equal k^(2r+1). |
| `func (g *GeneralRule) Apply(neighbourhood []int) int` | Apply implements Rule1D by indexing the full lookup table. |
| `func (g *GeneralRule) Index(neighbourhood []int) int` | Index returns the base-k index of a neighbourhood slice (leftmost cell most significant), the position used to look it up in the table. |
| `func (g *GeneralRule) Radius() int` | Radius reports the neighbourhood radius r. |
| `func (g *GeneralRule) States() int` | States reports the number of cell states k. |
| `func (g *GeneralRule) String() string` | String returns a description such as "general k=2 r=1". |
| `func (g *GeneralRule) Table() []int` | Table returns a copy of the full lookup table. |

</details>

<details>
<summary><code>Grid</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Acorn() *Grid` | Acorn returns the acorn methuselah, which stabilises only after thousands of Conway generations. |
| `func Beacon() *Grid` | Beacon returns the beacon, a Conway period-2 oscillator. |
| `func Beehive() *Grid` | Beehive returns the beehive still life. |
| `func Blinker() *Grid` | Blinker returns the horizontal blinker, a Conway period-2 oscillator. |
| `func Block() *Grid` | Block returns the 2x2 still life (Conway period 1). |
| `func CenteredPattern(pattern *Grid, rows, cols int) *Grid` | CenteredPattern returns a rows x cols grid with pattern stamped so its bounding box is centred. |
| `func Diehard() *Grid` | Diehard returns the diehard methuselah, which vanishes after 130 Conway generations. |
| `func Glider() *Grid` | Glider returns the glider, the smallest Conway spaceship, travelling diagonally with period 4. |
| `func GliderOn(rows, cols int) *Grid` | GliderOn returns a rows x cols grid containing a single glider stamped near the top-left corner. |
| `func GridFromStrings(lines []string, on rune) (*Grid, error)` | GridFromStrings builds a grid from a slice of equal-length strings, mapping each rune equal to on to 1 and all others to 0. |
| `func LWSS() *Grid` | LWSS returns the lightweight spaceship, a Conway period-4 orthogonal spaceship. |
| `func LifeEvolve(l *LifeRule, g *Grid, steps int, nb Neighbourhood, bc Boundary) []*Grid` | LifeEvolve runs a life-like rule for steps generations and returns steps+1 grids beginning with a copy of g. |
| `func LifeStep(l *LifeRule, g *Grid, nb Neighbourhood, bc Boundary) *Grid` | LifeStep advances a grid one generation under a life-like rule, neighbourhood and boundary condition, returning a new grid. |
| `func MargolusEvolve(m MargolusRule, g *Grid, steps int) ([]*Grid, error)` | MargolusEvolve runs a Margolus block rule for steps time steps, alternating the block partition each step (starting with the even partition), and… |
| `func MargolusStep(m MargolusRule, g *Grid, even bool) (*Grid, error)` | MargolusStep applies a Margolus block rule to a binary grid for one time step. |
| `func NewGrid(rows, cols int) *Grid` | NewGrid returns a rows x cols grid of zeros. |
| `func NewGridFrom(cells [][]int) (*Grid, error)` | NewGridFrom builds a grid from a rectangular slice of rows, mapping non-zero entries to 1. |
| `func RPentomino() *Grid` | RPentomino returns the R-pentomino, a famous long-lived Conway methuselah. |
| `func Toad() *Grid` | Toad returns the toad, a Conway period-2 oscillator. |
| `func (g *Grid) At(r, c int) int` | At returns the value of cell (r, c). |
| `func (g *Grid) Cells() [][]int` | Cells returns a fresh row-major two-dimensional copy of the grid contents. |
| `func (g *Grid) Clone() *Grid` | Clone returns an independent copy of the grid. |
| `func (g *Grid) Cols() int` | Cols returns the number of columns. |
| `func (g *Grid) Equal(other *Grid) bool` | Equal reports whether two grids have identical dimensions and contents. |
| `func (g *Grid) Population() int` | Population returns the number of live cells. |
| `func (g *Grid) Rows() int` | Rows returns the number of rows. |
| `func (g *Grid) Set(r, c, v int)` | Set writes v (any non-zero value becomes 1) into cell (r, c). |
| `func (g *Grid) String() string` | String renders the grid with '.' for dead and '#' for live cells. |

</details>

<details>
<summary><code>LifeRule</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Amoeba() *LifeRule` | Amoeba returns the Amoeba rule B357/S1358. |
| `func Anneal() *LifeRule` | Anneal returns the Anneal (twisted majority) rule B4678/S35678. |
| `func Conway() *LifeRule` | Conway returns Conway's Game of Life, B3/S23. |
| `func Coral() *LifeRule` | Coral returns the Coral rule B3/S45678. |
| `func DayAndNight() *LifeRule` | DayAndNight returns the Day & Night rule B3678/S34678, which is symmetric under exchange of live and dead cells. |
| `func Diamoeba() *LifeRule` | Diamoeba returns the Diamoeba rule B35678/S5678. |
| `func HighLife() *LifeRule` | HighLife returns the HighLife rule B36/S23, notable for its replicator. |
| `func LifeWithoutDeath() *LifeRule` | LifeWithoutDeath returns B3/S012345678, in which cells never die (Life without Death). |
| `func Maze() *LifeRule` | Maze returns the Maze rule B3/S12345. |
| `func Move() *LifeRule` | Move returns the Move rule B368/S245, which supports slow-moving spaceships. |
| `func NewLifeRule(born, survive []int) (*LifeRule, error)` | NewLifeRule builds a life-like rule from explicit birth and survival counts. |
| `func ParseRuleString(s string) (*LifeRule, error)` | ParseRuleString parses a life-like rulestring in B/S form (for example "B3/S23") or in the survival/birth Golly form ("23/3"). |
| `func Replicator() *LifeRule` | Replicator returns the Replicator rule B1357/S1357, which copies any pattern. |
| `func Seeds() *LifeRule` | Seeds returns the Seeds rule B2/S, in which every live cell dies each generation. |
| `func (l *LifeRule) Next(alive, neighbours int) int` | Next returns the next value of a cell given its current value and its live neighbour count. |
| `func (l *LifeRule) String() string` | String returns the rule in canonical B/S form, for example "B3/S23". |

</details>

<details>
<summary><code>MargolusRule</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func CrittersRule() MargolusRule` | CrittersRule returns an involutive, Critters-style reversible block rule: blocks with exactly two live cells are left unchanged, and every other… |
| `func Rotate180Rule() MargolusRule` | Rotate180Rule returns the Margolus rule that rotates every 2x2 block by 180 degrees (nw<->se, ne<->sw). |
| `func (m MargolusRule) Inverse() (MargolusRule, error)` | Inverse returns the inverse Margolus rule. |
| `func (m MargolusRule) IsPermutation() bool` | IsPermutation reports whether the rule is a bijection of the 16 block states and therefore reversible. |

</details>

<details>
<summary><code>Neighbourhood</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (n Neighbourhood) String() string` | String names the neighbourhood. |

</details>

<details>
<summary><code>OuterTotalisticRule</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewOuterTotalisticRule(k, r, code int) (*OuterTotalisticRule, error)` | NewOuterTotalisticRule builds a k-state, radius-r outer-totalistic rule from its code. |
| `func (o *OuterTotalisticRule) Apply(neighbourhood []int) int` | Apply implements Rule1D. neighbourhood must have length 2r+1; the centre cell is the middle element. |
| `func (o *OuterTotalisticRule) Code() int` | Code returns the base-k code identifying the rule. |
| `func (o *OuterTotalisticRule) OuterMax() int` | OuterMax returns the largest possible outer-neighbour sum, 2r(k-1). |
| `func (o *OuterTotalisticRule) Radius() int` | Radius reports the neighbourhood radius r. |
| `func (o *OuterTotalisticRule) States() int` | States reports the number of cell states k. |
| `func (o *OuterTotalisticRule) String() string` | String returns a description such as "outer-totalistic k=2 r=1 code=6". |
| `func (o *OuterTotalisticRule) Table() []int` | Table returns a copy of the output table. |

</details>

<details>
<summary><code>SecondOrderCA</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewSecondOrderCA(r ElementaryRule, prev, cur []int, bc Boundary) (*SecondOrderCA, error)` | NewSecondOrderCA creates a second-order automaton from a rule and two initial generations. |
| `func (s *SecondOrderCA) Clone() *SecondOrderCA` | Clone returns an independent copy of the automaton. |
| `func (s *SecondOrderCA) State() []int` | State returns a copy of the current configuration. |
| `func (s *SecondOrderCA) Step() []int` | Step advances the automaton one generation forward and returns the new current configuration. |
| `func (s *SecondOrderCA) StepBack() []int` | StepBack reverses the automaton one generation and returns the recovered configuration, exactly undoing a call to Step. |

</details>

<details>
<summary><code>TotalisticRule</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewTotalisticRule(k, r, code int) (*TotalisticRule, error)` | NewTotalisticRule builds a k-state, radius-r totalistic rule from its code. |
| `func (t *TotalisticRule) Apply(neighbourhood []int) int` | Apply implements Rule1D by summing the neighbourhood and looking up the sum in the table. |
| `func (t *TotalisticRule) Code() int` | Code returns the base-k code identifying the rule. |
| `func (t *TotalisticRule) MaxSum() int` | MaxSum returns the largest possible neighbourhood sum, (2r+1)(k-1). |
| `func (t *TotalisticRule) Radius() int` | Radius reports the neighbourhood radius r. |
| `func (t *TotalisticRule) States() int` | States reports the number of cell states k. |
| `func (t *TotalisticRule) String() string` | String returns a description such as "totalistic k=2 r=1 code=90". |
| `func (t *TotalisticRule) Table() []int` | Table returns a copy of the output table indexed by neighbourhood sum. |

</details>

### Variables

`ErrBadState`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/cellular).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
