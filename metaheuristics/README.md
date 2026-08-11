# metaheuristics — self-contained collection of metaheuristic global optimizers over continuous (R^n)

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/metaheuristics.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/metaheuristics)

Package metaheuristics implements a self-contained collection of metaheuristic
global optimizers over continuous (R^n) and combinatorial search spaces,
together with a suite of standard benchmark objective functions used to
validate them.

- Genetic algorithms with configurable selection, crossover, mutation and
elitism operators. - Simulated annealing with a family of cooling schedules
(exponential, geometric, linear, logarithmic, Boltzmann and Cauchy). -
Particle swarm optimization (PSO). - Differential evolution (DE) with the
classic mutation strategies. - A compact CMA-ES (covariance matrix adaptation
evolution strategy). - Tabu search over neighbourhoods. - Ant colony
optimization (ACO) for the travelling salesman problem. - Harmony search. -
Hill climbing with random restarts.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/metaheuristics"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AcceptanceProbability(curF, candF, T float64) float64` | AcceptanceProbability returns the Metropolis acceptance probability of moving from an energy of curF to a candidate energy candF at temperature T.… |
| `func Ackley(x []float64) float64` | Ackley is the multimodal Ackley benchmark with a nearly flat outer region and a deep central basin. |
| `func AnnealDiscrete( energy func(state []int) float64, neighbor func(state []int, rng *RNG) []int, initial []int, T0, alpha float64, iterations int, rng *RNG, ) (best []int, bestEnergy float64, err error)` | AnnealDiscrete minimizes an objective over an abstract discrete state space using simulated annealing. |
| `func ArithmeticCrossover(a, b []float64, rng *RNG) ([]float64, []float64)` | ArithmeticCrossover returns the two convex combinations t*a+(1-t)*b and (1-t)*a+t*b for a random blending factor t in [0,1]. |
| `func Beale(x []float64) float64` | Beale is the 2-D Beale benchmark with global minimum 0 at (3, 0.5). |
| `func BitFlipMutation(g []int, rate float64, rng *RNG) []int` | BitFlipMutation flips each bit of a binary genome (values 0 or 1) with probability rate, returning a new genome. |
| `func BlendCrossover(a, b []float64, alpha float64, rng *RNG) ([]float64, []float64)` | BlendCrossover implements BLX-alpha: each child gene is drawn uniformly from an interval extended by alpha times the parents' gene spread on both… |
| `func Booth(x []float64) float64` | Booth is the 2-D Booth benchmark f(x,y) = (x+2y-7)^2 + (2x+y-5)^2 with global minimum 0 at (1, 3). |
| `func Bukin6(x []float64) float64` | Bukin6 is the 2-D Bukin function N.6 with global minimum 0 at (-10, 1). |
| `func Clamp(x, lo, hi float64) float64` | Clamp returns x clamped to the closed interval [lo, hi]. |
| `func ConstrictionFactor(c1, c2 float64) float64` | ConstrictionFactor returns Clerc's constriction coefficient chi for the sum of acceleration coefficients phi = c1 + c2 (which must exceed 4). |
| `func CrossInTray(x []float64) float64` | CrossInTray is the 2-D cross-in-tray benchmark with four global minima of value approximately -2.06261. |
| `func DistanceMatrix(cities []City) [][]float64` | DistanceMatrix builds the symmetric matrix of pairwise Euclidean distances for the given cities. |
| `func DixonPrice(x []float64) float64` | DixonPrice is a unimodal benchmark f(x) = (x_1-1)^2 + sum_{i>=2} i (2 x_i^2 - x_{i-1})^2. |
| `func DropWave(x []float64) float64` | DropWave is the 2-D drop-wave benchmark with global minimum -1 at the origin. |
| `func Easom(x []float64) float64` | Easom is the 2-D Easom benchmark with a single narrow global minimum of -1 at (pi, pi). |
| `func Eggholder(x []float64) float64` | Eggholder is the 2-D Eggholder benchmark with global minimum approximately -959.6407 at (512, 404.2319). |
| `func EuclideanDistance(a, b City) float64` | EuclideanDistance returns the straight-line distance between two cities. |
| `func GaussianMutation(g []float64, rate, sigma float64, b Bounds, rng *RNG) []float64` | GaussianMutation adds independent zero-mean Gaussian noise to each gene with probability rate; the standard deviation is sigma times the box width. |
| `func GoldsteinPrice(x []float64) float64` | GoldsteinPrice is the 2-D Goldstein-Price benchmark with global minimum 3 at (0, -1). |
| `func Griewank(x []float64) float64` | Griewank is a multimodal benchmark f(x) = 1 + sum x_i^2/4000 - prod cos(x_i/sqrt(i+1)). |
| `func Himmelblau(x []float64) float64` | Himmelblau is the 2-D Himmelblau benchmark f(x,y) = (x^2+y-11)^2 + (x+y^2-7)^2. |
| `func HolderTable(x []float64) float64` | HolderTable is the 2-D Holder table benchmark with four global minima of value approximately -19.2085. |
| `func JacobiEigen(a [][]float64) (values []float64, vectors [][]float64)` | JacobiEigen computes the eigenvalues and orthonormal eigenvectors of the symmetric matrix a using the cyclic Jacobi rotation algorithm. |
| `func Lerp(a, b, t float64) float64` | Lerp returns the linear interpolation (1-t)*a + t*b. |
| `func LeviN13(x []float64) float64` | LeviN13 is the 2-D Levi function N.13 with global minimum 0 at (1, 1). |
| `func Matyas(x []float64) float64` | Matyas is the 2-D Matyas benchmark f(x,y) = 0.26(x^2+y^2) - 0.48 x y with global minimum 0 at the origin. |
| `func McCormick(x []float64) float64` | McCormick is the 2-D McCormick benchmark with global minimum approximately -1.9133 at (-0.54719, -1.54719). |
| `func Michalewicz(x []float64) float64` | Michalewicz is the multimodal Michalewicz benchmark with steepness parameter m = 10: f(x) = -sum sin(x_i) sin^{2m}((i+1) x_i^2 / pi). |
| `func NearestNeighborTour(dist [][]float64, start int) []int` | NearestNeighborTour builds a tour greedily starting from the given city, always moving to the nearest unvisited city. |
| `func OnePointCrossover(a, b []float64, rng *RNG) ([]float64, []float64)` | OnePointCrossover recombines two parent genomes at a single random cut point, returning two children. |
| `func PolynomialMutation(g []float64, rate, eta float64, b Bounds, rng *RNG) []float64` | PolynomialMutation applies the polynomial mutation operator with distribution index eta to each gene with probability rate, respecting the box. |
| `func Powell(x []float64) float64` | Powell is the 4k-dimensional Powell benchmark with global minimum 0 at the origin. |
| `func RankSelect(p Population, sp float64, rng *RNG) int` | RankSelect returns an index using linear rank selection with selection pressure sp in [1, 2]: the best individual receives expected count sp and the… |
| `func Rastrigin(x []float64) float64` | Rastrigin is a highly multimodal benchmark f(x) = 10n + sum(x_i^2 - 10 cos(2*pi*x_i)). |
| `func Rosenbrock(x []float64) float64` | Rosenbrock is the classic banana-valley benchmark f(x) = sum_{i<n-1} [100 (x_{i+1}-x_i^2)^2 + (1-x_i)^2]. |
| `func RouletteSelect(p Population, rng *RNG) int` | RouletteSelect returns an index chosen with probability proportional to fitness after transforming the minimization objective into non-negative… |
| `func SBXCrossover(a, b []float64, eta float64, rng *RNG) ([]float64, []float64)` | SBXCrossover implements simulated binary crossover with distribution index eta (larger eta yields children closer to their parents). |
| `func SchafferN2(x []float64) float64` | SchafferN2 is the 2-D Schaffer function N.2 with global minimum 0 at the origin. |
| `func SchafferN4(x []float64) float64` | SchafferN4 is the 2-D Schaffer function N.4 with global minimum approximately 0.292579 at (0, 1.25313). |
| `func Schwefel(x []float64) float64` | Schwefel is a deceptive multimodal benchmark f(x) = 418.9829 n - sum x_i sin(sqrt(\|x_i\|)). |
| `func SixHumpCamel(x []float64) float64` | SixHumpCamel is the 2-D six-hump camel benchmark with two global minima of value approximately -1.0316. |
| `func Sphere(x []float64) float64` | Sphere is the sum-of-squares benchmark f(x) = sum x_i^2. |
| `func StochasticUniversalSampling(p Population, m int, rng *RNG) []int` | StochasticUniversalSampling returns m indices sampled from the population in a single sweep with minimum spread, using non-negative weights derived… |
| `func StyblinskiTang(x []float64) float64` | StyblinskiTang is a multimodal benchmark f(x) = 0.5 sum(x_i^4 - 16 x_i^2 + 5 x_i). |
| `func SumOfDifferentPowers(x []float64) float64` | SumOfDifferentPowers is the benchmark f(x) = sum \|x_i\|^(i+2). |
| `func SumSquares(x []float64) float64` | SumSquares is the weighted sum-of-squares benchmark f(x) = sum (i+1) x_i^2. |
| `func TabuSearchCombinatorial( initial []int, neighbors func(state []int) []TabuMove, cfg TabuSearchConfig, ) (best []int, bestValue float64, history []float64, err error)` | TabuSearchCombinatorial minimizes an objective over a discrete state space (states represented as int slices). |
| `func ThreeHumpCamel(x []float64) float64` | ThreeHumpCamel is the 2-D three-hump camel benchmark with global minimum 0 at the origin. |
| `func TourLength(tour []int, dist [][]float64) float64` | TourLength returns the total length of a closed tour visiting the cities in the order given by the permutation tour, using the distance matrix dist… |
| `func TournamentSelect(p Population, k int, rng *RNG) int` | TournamentSelect returns the index of the winner (lowest fitness) of a tournament among k uniformly-random members of the population. |
| `func Trid(x []float64) float64` | Trid is the Trid benchmark f(x) = sum (x_i-1)^2 - sum_{i>=2} x_i x_{i-1}. |
| `func TwoOpt(tour []int, dist [][]float64) ([]int, float64)` | TwoOpt improves a tour by repeatedly applying the best 2-opt move until no move shortens the tour, returning the locally optimal tour and its length. |
| `func TwoOptSwap(tour []int, i, j int) []int` | TwoOptSwap returns a new tour with the segment between positions i and j (inclusive) reversed. |
| `func TwoPointCrossover(a, b []float64, rng *RNG) ([]float64, []float64)` | TwoPointCrossover recombines two parents by swapping the segment between two random cut points. |
| `func UniformCrossover(a, b []float64, p float64, rng *RNG) ([]float64, []float64)` | UniformCrossover recombines two parents by independently swapping each gene with probability p. |
| `func UniformMutation(g []float64, rate float64, b Bounds, rng *RNG) []float64` | UniformMutation resets each gene, with probability rate, to a uniform value within its box coordinate. |
| `func VecAXPY(a []float64, s float64, b []float64) []float64` | VecAXPY returns a + s*b (a scaled vector accumulation). |
| `func VecAdd(a, b []float64) []float64` | VecAdd returns the element-wise sum a+b. |
| `func VecCopy(x []float64) []float64` | VecCopy returns a fresh copy of x. |
| `func VecDist(a, b []float64) float64` | VecDist returns the Euclidean distance between a and b. |
| `func VecDot(a, b []float64) float64` | VecDot returns the Euclidean inner product of a and b. |
| `func VecMean(xs [][]float64) []float64` | VecMean returns the coordinate-wise mean of the rows of xs. |
| `func VecNorm(x []float64) float64` | VecNorm returns the Euclidean (L2) norm of x. |
| `func VecNorm1(x []float64) float64` | VecNorm1 returns the L1 (taxicab) norm of x. |
| `func VecNormInf(x []float64) float64` | VecNormInf returns the L-infinity (maximum absolute value) norm of x. |
| `func VecScale(a []float64, s float64) []float64` | VecScale returns the vector s*a. |
| `func VecSub(a, b []float64) []float64` | VecSub returns the element-wise difference a-b. |
| `func Zakharov(x []float64) float64` | Zakharov is a unimodal benchmark f(x) = sum x_i^2 + (sum 0.5 i x_i)^2 + (sum 0.5 i x_i)^4. |

### Types

| Type | What it is |
| --- | --- |
| `ACOConfig` | ACOConfig configures `RunACO`, ant colony optimization for the TSP. |
| `ACOResult` | ACOResult holds the outcome of an ant colony run. |
| `AnnealConfig` | AnnealConfig configures `SimulatedAnnealing`. |
| `Benchmark` | Benchmark bundles a benchmark objective with metadata: a human-readable name, a recommended search box for a given dimension, and its known global… |
| `Bounds` | Bounds describes an axis-aligned box in R^n by its per-coordinate lower and upper limits. |
| `CMAESConfig` | CMAESConfig configures `RunCMAES`, the covariance matrix adaptation evolution strategy. |
| `City` | City is a point in the plane used to build travelling-salesman instances. |
| `ContinuousTabuConfig` | ContinuousTabuConfig configures `ContinuousTabuSearch`. |
| `CoolingSchedule` | CoolingSchedule maps an iteration index k (0-based) to a temperature, given the initial temperature T0 and the total number of iterations kMax. |
| `DEConfig` | DEConfig configures `RunDE`. |
| `DEStrategy` | DEStrategy names a differential-evolution mutation strategy. |
| `GAConfig` | GAConfig configures the real-coded genetic algorithm `RunGA`. |
| `HarmonyConfig` | HarmonyConfig configures `RunHarmonySearch`. |
| `HillClimbConfig` | HillClimbConfig configures the hill-climbing optimizers. |
| `Individual` | Individual is one member of a genetic-algorithm population: a real-valued genome together with its cached objective (fitness) value. |
| `ObjectiveFunc` | ObjectiveFunc is a real-valued function of a vector argument. |
| `PSOConfig` | PSOConfig configures `RunPSO`. |
| `Particle` | Particle is one member of a particle swarm: its current position and velocity, together with the best position it has personally visited. |
| `Population` | Population is a slice of individuals. |
| `RNG` | RNG is a deterministic pseudo-random number generator used by every stochastic routine in the package. |
| `Result` | Result records the outcome of an optimization run. |
| `TabuList` | TabuList is a bounded first-in-first-out set of recently visited moves, identified by string keys. |
| `TabuMove` | TabuMove describes a candidate neighbour in a combinatorial tabu search: the resulting state, its objective value, and a string key identifying the… |
| `TabuSearchConfig` | TabuSearchConfig configures `TabuSearchCombinatorial`. |

<details>
<summary><code>ACOConfig</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DefaultACOConfig() ACOConfig` | DefaultACOConfig returns a reasonable configuration. |

</details>

<details>
<summary><code>ACOResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func RunACO(dist [][]float64, cfg ACOConfig, rng *RNG) (ACOResult, error)` | RunACO solves a travelling salesman instance given by its distance matrix using ant colony optimization. |

</details>

<details>
<summary><code>AnnealConfig</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DefaultAnnealConfig(b Bounds) AnnealConfig` | DefaultAnnealConfig returns a reasonable configuration for the given box. |

</details>

<details>
<summary><code>Benchmark</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func StandardBenchmarks() []Benchmark` | StandardBenchmarks returns a slice of the scalable n-dimensional benchmark functions (Sphere, Rastrigin, Ackley, Rosenbrock, Griewank, Schwefel,… |

</details>

<details>
<summary><code>Bounds</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewBounds(lower, upper []float64) (Bounds, error)` | NewBounds constructs a `Bounds` from separate lower and upper slices. |
| `func UniformBounds(dim int, lo, hi float64) Bounds` | UniformBounds returns a `Bounds` of the given dimension in which every coordinate ranges over [lo, hi]. |
| `func (b Bounds) Center() []float64` | Center returns the geometric center of the box. |
| `func (b Bounds) Clip(x []float64) []float64` | Clip returns a copy of x with each coordinate clamped into the box. |
| `func (b Bounds) ClipInPlace(x []float64) []float64` | ClipInPlace clamps each coordinate of x into the box, modifying x, and returns x for convenience. |
| `func (b Bounds) Contains(x []float64) bool` | Contains reports whether x lies within the closed box. |
| `func (b Bounds) Dim() int` | Dim returns the dimension of the box. |
| `func (b Bounds) Reflect(x []float64) []float64` | Reflect returns a copy of x in which coordinates outside the box are reflected back inside across the violated face. |
| `func (b Bounds) Valid() bool` | Valid reports whether the bounds are well formed: non-empty, of equal length and with every lower limit not exceeding its upper limit. |
| `func (b Bounds) Width() []float64` | Width returns the per-coordinate widths (Upper-Lower) of the box. |

</details>

<details>
<summary><code>CMAESConfig</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DefaultCMAESConfig(b Bounds) CMAESConfig` | DefaultCMAESConfig returns a reasonable configuration for the given box. |

</details>

<details>
<summary><code>ContinuousTabuConfig</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DefaultContinuousTabuConfig(b Bounds) ContinuousTabuConfig` | DefaultContinuousTabuConfig returns a reasonable configuration for the box. |

</details>

<details>
<summary><code>CoolingSchedule</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BoltzmannCooling() CoolingSchedule` | BoltzmannCooling returns the schedule T_k = T0 / (1 + ln(1+k)). |
| `func CauchyCooling() CoolingSchedule` | CauchyCooling returns the fast Cauchy schedule T_k = T0 / (1 + k). |
| `func ExponentialCooling(alpha float64) CoolingSchedule` | ExponentialCooling returns a schedule T_k = T0 * alpha^k with 0 < alpha < 1. |
| `func GeometricCooling(alpha float64) CoolingSchedule` | GeometricCooling is a synonym for `ExponentialCooling`; each step multiplies the temperature by alpha. |
| `func LinearCooling() CoolingSchedule` | LinearCooling returns a schedule that decreases the temperature linearly from T0 at k=0 to (approximately) 0 at k=kMax. |
| `func LogarithmicCooling() CoolingSchedule` | LogarithmicCooling returns the classical Boltzmann schedule T_k = T0 / ln(k + e), which cools slowly enough to guarantee (in the limit) convergence… |
| `func QuadraticCooling() CoolingSchedule` | QuadraticCooling returns a schedule that decreases quadratically from T0 to 0 over kMax iterations. |

</details>

<details>
<summary><code>DEConfig</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DefaultDEConfig(b Bounds) DEConfig` | DefaultDEConfig returns a reasonable configuration for the given box. |

</details>

<details>
<summary><code>DEStrategy</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (s DEStrategy) String() string` | String returns the conventional name of the strategy, such as "DE/rand/1". |

</details>

<details>
<summary><code>GAConfig</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DefaultGAConfig(b Bounds) GAConfig` | DefaultGAConfig returns a reasonable configuration for the given box. |

</details>

<details>
<summary><code>HarmonyConfig</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DefaultHarmonyConfig(b Bounds) HarmonyConfig` | DefaultHarmonyConfig returns a reasonable configuration for the given box. |

</details>

<details>
<summary><code>HillClimbConfig</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DefaultHillClimbConfig(b Bounds) HillClimbConfig` | DefaultHillClimbConfig returns a reasonable configuration for the given search box. |

</details>

<details>
<summary><code>Individual</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (ind Individual) Clone() Individual` | Clone returns a deep copy of the individual. |

</details>

<details>
<summary><code>ObjectiveFunc</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Negate(f ObjectiveFunc) ObjectiveFunc` | Negate returns a new objective equal to -f(x). |
| `func Penalized(f ObjectiveFunc, b Bounds, weight float64) ObjectiveFunc` | Penalized returns a new objective that adds weight times the squared amount by which x violates the box b to the value of f. |
| `func Scale(f ObjectiveFunc, a, b float64) ObjectiveFunc` | Scale returns a new objective a*f(x)+b. |
| `func Shift(f ObjectiveFunc, offset []float64) ObjectiveFunc` | Shift returns a new objective f(x - offset), translating the landscape so that the original optimum at o moves to o+offset. |

</details>

<details>
<summary><code>PSOConfig</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DefaultPSOConfig(b Bounds) PSOConfig` | DefaultPSOConfig returns a reasonable configuration for the given box using the well-known constriction-inspired coefficients. |

</details>

<details>
<summary><code>Particle</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func InitSwarm(f ObjectiveFunc, b Bounds, n int, rng *RNG) []Particle` | InitSwarm creates a swarm of n particles with random positions and small random velocities, evaluating each particle's initial fitness. |

</details>

<details>
<summary><code>Population</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func InitPopulation(f ObjectiveFunc, b Bounds, n int, rng *RNG) Population` | InitPopulation creates a population of n individuals with genomes drawn uniformly from the box and their fitnesses evaluated. |
| `func (p Population) Best() Individual` | Best returns a copy of the fittest (minimum-fitness) individual. |
| `func (p Population) MeanFitness() float64` | MeanFitness returns the arithmetic mean of the population's fitness values. |
| `func (p Population) SortByFitness()` | SortByFitness sorts the population in ascending order of fitness (best, i.e. |

</details>

<details>
<summary><code>RNG</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewRNG(seed int64) *RNG` | NewRNG returns a deterministic `RNG` seeded with the supplied value. |
| `func (g *RNG) Choice(weights []float64) int` | Choice returns a pseudo-random element index in [0, n) selected with probability proportional to the supplied non-negative weights. |
| `func (g *RNG) Float64() float64` | Float64 returns a pseudo-random float64 in the half-open interval [0, 1). |
| `func (g *RNG) Float64Range(lo, hi float64) float64` | Float64Range returns a pseudo-random float64 uniformly distributed in the half-open interval [lo, hi). |
| `func (g *RNG) Gaussian(mean, std float64) float64` | Gaussian returns a normally distributed float64 with the given mean and standard deviation. |
| `func (g *RNG) GaussianVec(n int, mean, std float64) []float64` | GaussianVec returns a new vector of length n with each coordinate drawn from a normal distribution of the given mean and standard deviation. |
| `func (g *RNG) IntRange(lo, hi int) int` | IntRange returns a pseudo-random int in [lo, hi). |
| `func (g *RNG) Intn(n int) int` | Intn returns a pseudo-random int in [0, n). |
| `func (g *RNG) NormFloat64() float64` | NormFloat64 returns a normally distributed float64 with mean 0 and standard deviation 1. |
| `func (g *RNG) Perm(n int) []int` | Perm returns a pseudo-random permutation of the integers [0, n). |
| `func (g *RNG) Sample(n, k int) []int` | Sample returns k distinct indices drawn without replacement from [0, n). |
| `func (g *RNG) Shuffle(s []int)` | Shuffle pseudo-randomly permutes the elements of s in place using the Fisher-Yates algorithm. |
| `func (g *RNG) ShuffleFloat(s []float64)` | ShuffleFloat pseudo-randomly permutes the elements of s in place. |
| `func (g *RNG) UniformVec(b Bounds) []float64` | UniformVec returns a new vector of length b.Dim() with each coordinate drawn uniformly from the corresponding [lo, hi) interval of b. |

</details>

<details>
<summary><code>Result</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ContinuousTabuSearch(f ObjectiveFunc, cfg ContinuousTabuConfig, start []float64, rng *RNG) (Result, error)` | ContinuousTabuSearch minimizes f over the box by, at each step, sampling several Gaussian neighbours and moving to the best whose quantized grid cell… |
| `func CoordinateDescent(f ObjectiveFunc, cfg HillClimbConfig, start []float64) (Result, error)` | CoordinateDescent performs axis-parallel line minimization: each iteration it probes both directions along every coordinate with the current step and… |
| `func HillClimb(f ObjectiveFunc, cfg HillClimbConfig, start []float64, rng *RNG) (Result, error)` | HillClimb runs stochastic hill climbing from the supplied start point, returning the best point found. |
| `func HillClimbRestarts(f ObjectiveFunc, cfg HillClimbConfig, restarts int, rng *RNG) (Result, error)` | HillClimbRestarts runs HillClimb from restarts independent random starting points drawn uniformly from the box, returning the best result across all… |
| `func RandomSearch(f ObjectiveFunc, b Bounds, iterations int, rng *RNG) (Result, error)` | RandomSearch samples MaxIterations points uniformly from the box and returns the best. |
| `func RunCMAES(f ObjectiveFunc, cfg CMAESConfig, rng *RNG) (Result, error)` | RunCMAES minimizes f over R^n using CMA-ES, clipping candidate samples to the box for evaluation. |
| `func RunDE(f ObjectiveFunc, cfg DEConfig, rng *RNG) (Result, error)` | RunDE minimizes f over the box using differential evolution with binomial crossover. |
| `func RunGA(f ObjectiveFunc, cfg GAConfig, rng *RNG) (Result, error)` | RunGA evolves a real-coded genetic algorithm and returns the best individual found. |
| `func RunHarmonySearch(f ObjectiveFunc, cfg HarmonyConfig, rng *RNG) (Result, error)` | RunHarmonySearch minimizes f over the box using the harmony search metaheuristic. |
| `func RunPSO(f ObjectiveFunc, cfg PSOConfig, rng *RNG) (Result, error)` | RunPSO minimizes f over the box using particle swarm optimization with a global (gbest) topology. |
| `func SimulatedAnnealing(f ObjectiveFunc, cfg AnnealConfig, start []float64, rng *RNG) (Result, error)` | SimulatedAnnealing minimizes f over the continuous box using the Metropolis criterion and the configured cooling schedule, starting from the given… |

</details>

<details>
<summary><code>TabuList</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewTabuList(tenure int) *TabuList` | NewTabuList creates an empty tabu list that remembers up to tenure moves. |
| `func (t *TabuList) Add(key string)` | Add records a move key, evicting the oldest entry when the tenure is exceeded. |
| `func (t *TabuList) Clear()` | Clear empties the tabu list. |
| `func (t *TabuList) Contains(key string) bool` | Contains reports whether the given move key is currently tabu. |
| `func (t *TabuList) Len() int` | Len returns the number of moves currently in the list. |
| `func (t *TabuList) Tenure() int` | Tenure returns the maximum number of moves the list retains. |

</details>

### Variables

`ErrDimMismatch`, `ErrEmptyBounds`, `ErrInvalidConfig`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/metaheuristics).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
