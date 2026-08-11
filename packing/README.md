# packing — implements classical packing and covering algorithms

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/packing.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/packing)

Package packing implements classical packing and covering algorithms.

- One dimensional bin packing. The online heuristics next-fit, first-fit,
best-fit and worst-fit together with their offline "decreasing" variants
(items sorted by non-increasing size before packing), plus the classical
Martello-Toth lower bounds L1 and L2 and the known worst-case approximation
ratios.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/packing"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AnCenterDensity(n int) float64` | AnCenterDensity returns the center density of A_n, 2^(-n/2) / sqrt(n+1). |
| `func AnCoveringRadius(n int) float64` | AnCoveringRadius returns the covering radius of A_n. |
| `func AnCovolume(n int) float64` | AnCovolume returns the covolume of the root lattice A_n, sqrt(n+1). |
| `func AnDensity(n int) float64` | AnDensity returns the sphere-packing density of A_n. |
| `func AnKissingNumber(n int) int` | AnKissingNumber returns the kissing number of A_n, n(n+1). |
| `func AnMinimalNorm(n int) float64` | AnMinimalNorm returns the minimal squared vector length of A_n, which is 2. |
| `func AnPackingRadius(n int) float64` | AnPackingRadius returns the packing radius of A_n, 1/sqrt(2). |
| `func AnStarCoveringRadius(n int) float64` | AnStarCoveringRadius returns the covering radius of A_n^*, whose square is n(n+2)/(12(n+1)). |
| `func AnStarCovolume(n int) float64` | AnStarCovolume returns the covolume of the dual lattice A_n^*, 1/sqrt(n+1). |
| `func AnStarThickness(n int) float64` | AnStarThickness returns the covering thickness of A_n^*. |
| `func AnThickness(n int) float64` | AnThickness returns the covering thickness of A_n. |
| `func BallSurface(n int, r float64) float64` | BallSurface returns the surface area of the (n-1)-sphere bounding an n-dimensional ball of radius r, n * V_n * r^(n-1). |
| `func BallVolume(n int, r float64) float64` | BallVolume returns the volume of an n-dimensional Euclidean ball of radius r, pi^(n/2) r^n / Gamma(n/2 + 1). |
| `func BestFitApproxRatio() float64` | BestFitApproxRatio returns the asymptotic worst-case approximation ratio of best-fit, which is 1.7 (identical to first-fit). |
| `func BestLatticeCenterDensity(n int) float64` | BestLatticeCenterDensity returns the center density corresponding to `BestLatticeDensity` in dimension n. |
| `func BestLatticeDensity(n int) float64` | BestLatticeDensity returns the packing density of the densest lattice among the standard family (Z^n, A_n, D_n, E6, E7, E8, Leech) implemented here… |
| `func BinPackingLowerBound(sizes []float64, capacity float64) int` | BinPackingLowerBound returns the strongest lower bound implemented here on the number of bins, the maximum of `LowerBoundL1` and `LowerBoundL2`. |
| `func CenterDensity(rho, covolume float64, n int) float64` | CenterDensity returns the center density delta = rho^n / V of a lattice with packing radius rho and covolume V in dimension n. |
| `func CircleInCircleDensity(n int) (density float64, ok bool)` | CircleInCircleDensity returns the fraction of the enclosing circle covered by n equal circles packed with the best known ratio, n / (R/r)^2. |
| `func CircleInCircleKnown(n int) bool` | CircleInCircleKnown reports whether a best known packing of n equal circles in a circle is tabulated. |
| `func CircleInCircleRatio(n int) (ratio float64, ok bool)` | CircleInCircleRatio returns the best known ratio R/r for packing n equal circles of radius r inside a circle of radius R, together with ok reporting… |
| `func CircleInSquareDensity(n int) (density float64, ok bool)` | CircleInSquareDensity returns the fraction of a unit square covered by n equal circles packed with the best known radius, n*pi*r^2. |
| `func CircleInSquareKnown(n int) bool` | CircleInSquareKnown reports whether a best known packing of n equal circles in a square is tabulated (including the analytic grid and single-circle… |
| `func CircleInSquareRadius(n int) (r float64, ok bool)` | CircleInSquareRadius returns the radius of n equal circles packed in a unit square using the best known spread from `CircleInSquareSpread`. |
| `func CircleInSquareSpread(n int) (m float64, ok bool)` | CircleInSquareSpread returns the best known optimal point spread m(n) for n points in the unit square (the largest achievable minimum pairwise… |
| `func DensityToCenterDensity(density float64, n int) float64` | DensityToCenterDensity converts a packing density Delta back to a center density delta = Delta / V_n in dimension n. |
| `func DnCenterDensity(n int) float64` | DnCenterDensity returns the center density of D_n, 2^(-(n+2)/2). |
| `func DnCoveringRadius(n int) float64` | DnCoveringRadius returns the covering radius of D_n. |
| `func DnCovolume(n int) float64` | DnCovolume returns the covolume of the checkerboard lattice D_n, which is 2. |
| `func DnDensity(n int) float64` | DnDensity returns the sphere-packing density of D_n. |
| `func DnKissingNumber(n int) int` | DnKissingNumber returns the kissing number of D_n, 2n(n-1). |
| `func DnMinimalNorm(n int) float64` | DnMinimalNorm returns the minimal squared vector length of D_n, which is 2. |
| `func DnPackingRadius(n int) float64` | DnPackingRadius returns the packing radius of D_n, 1/sqrt(2). |
| `func DnThickness(n int) float64` | DnThickness returns the covering thickness of D_n. |
| `func E6CenterDensity() float64` | E6CenterDensity returns the center density of E6, 1/(8 sqrt(3)). |
| `func E6CoveringRadius() float64` | E6CoveringRadius returns the covering radius of E6, whose square is 4/3. |
| `func E6Covolume() float64` | E6Covolume returns the covolume of E6, sqrt(3). |
| `func E6Density() float64` | E6Density returns the sphere-packing density of E6. |
| `func E6Dimension() int` | E6Dimension returns the dimension of the lattice E6, which is 6. |
| `func E6KissingNumber() int` | E6KissingNumber returns the kissing number of E6, which is 72. |
| `func E6MinimalNorm() float64` | E6MinimalNorm returns the minimal squared vector length of E6, which is 2. |
| `func E6PackingRadius() float64` | E6PackingRadius returns the packing radius of E6, 1/sqrt(2). |
| `func E6Thickness() float64` | E6Thickness returns the covering thickness of E6. |
| `func E7CenterDensity() float64` | E7CenterDensity returns the center density of E7, 1/16. |
| `func E7CoveringRadius() float64` | E7CoveringRadius returns the covering radius of E7, whose square is 3/2. |
| `func E7Covolume() float64` | E7Covolume returns the covolume of E7, sqrt(2). |
| `func E7Density() float64` | E7Density returns the sphere-packing density of E7. |
| `func E7Dimension() int` | E7Dimension returns the dimension of the lattice E7, which is 7. |
| `func E7KissingNumber() int` | E7KissingNumber returns the kissing number of E7, which is 126. |
| `func E7MinimalNorm() float64` | E7MinimalNorm returns the minimal squared vector length of E7, which is 2. |
| `func E7PackingRadius() float64` | E7PackingRadius returns the packing radius of E7, 1/sqrt(2). |
| `func E7Thickness() float64` | E7Thickness returns the covering thickness of E7. |
| `func E8CenterDensity() float64` | E8CenterDensity returns the center density of E8, which is 1/16. |
| `func E8CoveringRadius() float64` | E8CoveringRadius returns the covering radius of E8, which is 1. |
| `func E8Covolume() float64` | E8Covolume returns the covolume of E8, which is 1 (E8 is unimodular). |
| `func E8Density() float64` | E8Density returns the sphere-packing density of E8, pi^4/384. |
| `func E8Dimension() int` | E8Dimension returns the dimension of the lattice E8, which is 8. |
| `func E8KissingNumber() int` | E8KissingNumber returns the kissing number of E8, which is 240. |
| `func E8MinimalNorm() float64` | E8MinimalNorm returns the minimal squared vector length of E8, which is 2. |
| `func E8PackingRadius() float64` | E8PackingRadius returns the packing radius of E8, 1/sqrt(2). |
| `func E8Thickness() float64` | E8Thickness returns the covering thickness of E8, V_8 = pi^4/24. |
| `func FirstFitApproxRatio() float64` | FirstFitApproxRatio returns the asymptotic worst-case approximation ratio of first-fit and best-fit, which is 1.7. |
| `func FirstFitDecreasingAdditiveBound() float64` | FirstFitDecreasingAdditiveBound returns the additive constant in the tight bound FFD(I) <= (11/9)OPT(I) + 6/9 due to Doesa et al. |
| `func FirstFitDecreasingApproxRatio() float64` | FirstFitDecreasingApproxRatio returns the asymptotic worst-case approximation ratio of first-fit-decreasing, which is 11/9. |
| `func FractionalKnapsack(values, weights []float64, capacity float64) (value float64, fraction []float64, err error)` | FractionalKnapsack solves the fractional (continuous) knapsack problem exactly by the greedy value-to-weight rule: items are taken whole in order of… |
| `func GreedyKnapsack(values, weights []float64, capacity float64) (value float64, chosen []int, weight float64, err error)` | GreedyKnapsack applies the greedy value-to-weight heuristic to the 0/1 knapsack problem: it considers items in order of decreasing density and takes… |
| `func GreedyKnapsackBest(values, weights []float64, capacity float64) (value float64, chosen []int, weight float64, err error)` | GreedyKnapsackBest runs the greedy heuristic but also compares it against the single most valuable item that fits, returning whichever gives more… |
| `func GridSpread(k int) float64` | GridSpread returns the point spread of the k-by-k square grid of k^2 points in the unit square, 1/(k-1). |
| `func HarmonicNumber(n int) float64` | HarmonicNumber returns the n-th harmonic number H(n) = 1 + 1/2 + ... |
| `func HermiteConstant(n int) float64` | HermiteConstant returns the Hermite constant gamma_n, the maximum over all n-dimensional lattices of the ratio (minimal norm) / covolume^(2/n), for… |
| `func HermiteConstantKnown(n int) bool` | HermiteConstantKnown reports whether the Hermite constant is proven in dimension n. |
| `func HexagonalCoveringThickness() float64` | HexagonalCoveringThickness returns the thickness of the thinnest covering of the plane by equal circles, 2*pi/sqrt(27) ~ 1.2092, achieved by the… |
| `func HexagonalKissingNumber() int` | HexagonalKissingNumber returns the kissing number of the hexagonal circle packing, which is 6. |
| `func HexagonalPackingDensity() float64` | HexagonalPackingDensity returns the density of the densest packing of equal circles in the plane, pi/sqrt(12) = pi/(2*sqrt(3)) ~ 0.9069, achieved by… |
| `func IsSetCover(universe int, sets [][]int, chosen []int) bool` | IsSetCover reports whether the chosen sets together cover the whole universe {0,...,universe-1}. |
| `func KissingNumberKnown(n int) bool` | KissingNumberKnown reports whether the optimal kissing number is known in dimension n. |
| `func KissingNumberLowerBound(n int) int` | KissingNumberLowerBound returns the best kissing number attainable by the standard lattice family used in this package for dimension n: Z^n gives 2n,… |
| `func KissingNumberOptimal(n int) int` | KissingNumberOptimal returns the exact optimal kissing number tau(n) in dimension n for the dimensions in which it is known (1, 2, 3, 4, 8, 24), and… |
| `func KnapsackFeasible(weights []float64, chosen []int, capacity float64) bool` | KnapsackFeasible reports whether the chosen items fit within the capacity. |
| `func KnapsackUpperBound(values, weights []float64, capacity float64) (float64, error)` | KnapsackUpperBound returns the linear-programming (fractional) relaxation value of the 0/1 knapsack instance, an upper bound on the integer optimum. |
| `func KnapsackValue(values []float64, chosen []int) float64` | KnapsackValue returns the total value of the items whose indices are given. |
| `func KnapsackWeight(weights []float64, chosen []int) float64` | KnapsackWeight returns the total weight of the items whose indices are given. |
| `func LeechCenterDensity() float64` | LeechCenterDensity returns the center density of the Leech lattice, which is 1. |
| `func LeechCoveringRadius() float64` | LeechCoveringRadius returns the covering radius of the Leech lattice, sqrt(2) (its deep holes lie at distance sqrt(2) and come in 23 types). |
| `func LeechCovolume() float64` | LeechCovolume returns the covolume of the Leech lattice, which is 1 (it is unimodular). |
| `func LeechDensity() float64` | LeechDensity returns the sphere-packing density of the Leech lattice, the volume of the unit 24-ball pi^12/12!. |
| `func LeechDimension() int` | LeechDimension returns the dimension of the Leech lattice, which is 24. |
| `func LeechKissingNumber() int` | LeechKissingNumber returns the kissing number of the Leech lattice, 196560. |
| `func LeechMinimalNorm() float64` | LeechMinimalNorm returns the minimal squared vector length of the Leech lattice, which is 4. |
| `func LeechPackingRadius() float64` | LeechPackingRadius returns the packing radius of the Leech lattice, which is 1 (minimal distance 2). |
| `func LeechThickness() float64` | LeechThickness returns the covering thickness of the Leech lattice. |
| `func LowerBoundL1(sizes []float64, capacity float64) int` | LowerBoundL1 returns the trivial "continuous" lower bound on the number of bins: the total size divided by the capacity, rounded up. |
| `func LowerBoundL2(sizes []float64, capacity float64) int` | LowerBoundL2 returns the Martello-Toth L2 lower bound on the number of bins. |
| `func MaxSetSize(sets [][]int) int` | MaxSetSize returns the size of the largest set, used as the parameter of `SetCoverGreedyBound`. |
| `func MaxSize(sizes []float64) float64` | MaxSize returns the largest item size, or 0 for an empty slice. |
| `func MinimalDistance(minNorm float64) float64` | MinimalDistance returns the shortest distance between distinct lattice points given the minimal norm (minimal squared vector length), sqrt(minNorm). |
| `func MinkowskiHlawkaBound(n int) float64` | MinkowskiHlawkaBound returns the Minkowski-Hlawka lower bound on the optimal lattice packing density in dimension n, zeta(n) / 2^(n-1). |
| `func ModifiedGreedyKnapsackApproxRatio() float64` | ModifiedGreedyKnapsackApproxRatio returns the worst-case approximation ratio of the modified greedy 0/1 knapsack heuristic `GreedyKnapsackBest`,… |
| `func NextFitApproxRatio() float64` | NextFitApproxRatio returns the asymptotic worst-case approximation ratio of next-fit relative to the optimum, which is 2. |
| `func NextFitDecreasingApproxRatio() float64` | NextFitDecreasingApproxRatio returns the asymptotic worst-case approximation ratio of next-fit-decreasing, the constant h_infinity ~ 1.6910. |
| `func PackingDensity(centerDensity float64, n int) float64` | PackingDensity returns the sphere-packing density Delta = delta * V_n given a center density delta in dimension n. |
| `func PackingRadius(minNorm float64) float64` | PackingRadius returns the packing radius of a lattice whose minimal squared vector length (minimal norm) is minNorm, namely sqrt(minNorm)/2. |
| `func RingDensity(k int) float64` | RingDensity returns the packing density k / RingRatio(k)^2 of a single ring of k circles inside the enclosing circle. |
| `func RingRatio(k int) float64` | RingRatio returns the enclosing ratio R/r for k equal circles arranged in a single symmetric ring touching the boundary of the outer circle, 1 +… |
| `func SetCoverGreedy(universe int, sets [][]int) ([]int, error)` | SetCoverGreedy applies the greedy set-cover heuristic to cover the universe {0,...,universe-1}: repeatedly choose the set covering the most still… |
| `func SetCoverGreedyBound(largestSetSize int) float64` | SetCoverGreedyBound returns the greedy set-cover approximation guarantee H(k), the k-th harmonic number, where k is the size of the largest set. |
| `func SquareCoveringThickness() float64` | SquareCoveringThickness returns the thickness of the covering of the plane by equal circles centered on the square lattice, pi/2 ~ 1.5708. |
| `func SquareKissingNumber() int` | SquareKissingNumber returns the kissing number of the square circle packing, which is 4. |
| `func SquarePackingDensity() float64` | SquarePackingDensity returns the density of equal circles packed on the square lattice, pi/4 ~ 0.7854. |
| `func SumSizes(sizes []float64) float64` | SumSizes returns the total of the item sizes. |
| `func Thickness(coveringRadius, covolume float64, n int) float64` | Thickness returns the covering thickness Theta = V_n * R^n / V of a lattice with covering radius R and covolume V in dimension n. |
| `func UnitBallVolume(n int) float64` | UnitBallVolume returns the volume V_n of the unit ball in R^n, pi^(n/2) / Gamma(n/2 + 1). |
| `func UnitSphereSurface(n int) float64` | UnitSphereSurface returns the surface area of the unit (n-1)-sphere in R^n, 2 pi^(n/2) / Gamma(n/2). |
| `func WeightedSetCoverGreedy(universe int, sets [][]int, costs []float64) (chosen []int, totalCost float64, err error)` | WeightedSetCoverGreedy applies the greedy heuristic to the weighted set-cover problem: repeatedly choose the set minimizing cost divided by the… |
| `func WorstFitApproxRatio() float64` | WorstFitApproxRatio returns the asymptotic worst-case approximation ratio of worst-fit, which is 2 (identical to next-fit). |
| `func ZnCenterDensity(n int) float64` | ZnCenterDensity returns the center density of Z^n, 2^(-n). |
| `func ZnCoveringRadius(n int) float64` | ZnCoveringRadius returns the covering radius of Z^n, sqrt(n)/2 (the deep hole is the all-halves center of a unit cube). |
| `func ZnCovolume(n int) float64` | ZnCovolume returns the covolume of the integer lattice Z^n, which is 1. |
| `func ZnDensity(n int) float64` | ZnDensity returns the sphere-packing density of Z^n, V_n / 2^n. |
| `func ZnKissingNumber(n int) int` | ZnKissingNumber returns the kissing number of Z^n, 2n (the +-1 unit vectors). |
| `func ZnMinimalNorm(n int) float64` | ZnMinimalNorm returns the minimal squared vector length of Z^n, which is 1. |
| `func ZnPackingRadius(n int) float64` | ZnPackingRadius returns the packing radius of Z^n, which is 1/2. |
| `func ZnThickness(n int) float64` | ZnThickness returns the covering thickness of Z^n, V_n * (sqrt(n)/2)^n. |

### Types

| Type | What it is |
| --- | --- |
| `Bin` | Bin is a single bin of a `Packing`. |
| `LatticeInfo` | LatticeInfo bundles the standard packing and covering invariants of a lattice: its name, dimension, covolume, minimal norm (squared minimal vector… |
| `Packing` | Packing is the result of a bin packing routine: an assignment of every item to exactly one `Bin` together with the capacity and the original item… |

<details>
<summary><code>LatticeInfo</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AnInfo(n int) LatticeInfo` | AnInfo returns the `LatticeInfo` for the root lattice A_n. |
| `func DnInfo(n int) LatticeInfo` | DnInfo returns the `LatticeInfo` for the checkerboard lattice D_n. |
| `func E6Info() LatticeInfo` | E6Info returns the `LatticeInfo` for the exceptional lattice E6. |
| `func E7Info() LatticeInfo` | E7Info returns the `LatticeInfo` for the exceptional lattice E7. |
| `func E8Info() LatticeInfo` | E8Info returns the `LatticeInfo` for the exceptional lattice E8. |
| `func LeechInfo() LatticeInfo` | LeechInfo returns the `LatticeInfo` for the Leech lattice Lambda_24. |
| `func ZnInfo(n int) LatticeInfo` | ZnInfo returns the `LatticeInfo` for the integer lattice Z^n. |

</details>

<details>
<summary><code>Packing</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AlmostWorstFit(sizes []float64, capacity float64) (Packing, error)` | AlmostWorstFit packs the items using the almost-worst-fit heuristic. |
| `func AlmostWorstFitDecreasing(sizes []float64, capacity float64) (Packing, error)` | AlmostWorstFitDecreasing sorts the items by non-increasing size and packs them with almost-worst-fit. |
| `func BestFit(sizes []float64, capacity float64) (Packing, error)` | BestFit packs the items using the best-fit heuristic. |
| `func BestFitDecreasing(sizes []float64, capacity float64) (Packing, error)` | BestFitDecreasing sorts the items by non-increasing size and packs them with best-fit. |
| `func FirstFit(sizes []float64, capacity float64) (Packing, error)` | FirstFit packs the items using the first-fit heuristic. |
| `func FirstFitDecreasing(sizes []float64, capacity float64) (Packing, error)` | FirstFitDecreasing sorts the items by non-increasing size and packs them with first-fit. |
| `func NextFit(sizes []float64, capacity float64) (Packing, error)` | NextFit packs the items using the next-fit heuristic. |
| `func NextFitDecreasing(sizes []float64, capacity float64) (Packing, error)` | NextFitDecreasing sorts the items by non-increasing size and packs them with next-fit. |
| `func WorstFit(sizes []float64, capacity float64) (Packing, error)` | WorstFit packs the items using the worst-fit heuristic. |
| `func WorstFitDecreasing(sizes []float64, capacity float64) (Packing, error)` | WorstFitDecreasing sorts the items by non-increasing size and packs them with worst-fit. |
| `func (p Packing) Fullness() float64` | Fullness returns the average fraction of bin capacity that is used, in the range [0,1]. |
| `func (p Packing) Loads() []float64` | Loads returns the load (sum of item sizes) of each bin in bin order. |
| `func (p Packing) MaxLoad() float64` | MaxLoad returns the largest bin load, or 0 for an empty packing. |
| `func (p Packing) MinLoad() float64` | MinLoad returns the smallest bin load, or 0 for an empty packing. |
| `func (p Packing) NumBins() int` | NumBins reports the number of bins used by the packing. |
| `func (p Packing) NumItems() int` | NumItems reports the total number of packed items. |
| `func (p Packing) TotalSize() float64` | TotalSize returns the sum of all item sizes. |
| `func (p Packing) Valid() bool` | Valid reports whether the packing is internally consistent: every item index appears exactly once, every recorded load matches the sum of its items,… |
| `func (p Packing) Waste() float64` | Waste returns the total unused capacity, NumBins*Capacity - TotalSize. |

</details>

### Variables

`ErrCapacity`, `ErrDimensionMismatch`, `ErrItemTooLarge`, `ErrNegativeSize`, `ErrNegativeWeight`, `ErrNotCoverable`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/packing).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
