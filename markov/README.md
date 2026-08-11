# markov — provides tools for finite-state Markov chains, hidden Markov models (HMMs)

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/markov.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/markov)

Package markov provides tools for finite-state Markov chains, hidden Markov
models (HMMs), and Markov-chain Monte Carlo (MCMC) using only the Go standard
library.

The central type for discrete-time chains is `MarkovChain`, a wrapper around a
row-stochastic transition matrix. It computes n-step transition matrices,
evolves probability distributions, finds stationary distributions (both by
power iteration and by exact linear solve), classifies states into
communicating classes, determines periodicity, irreducibility and ergodicity,
and handles absorbing chains via the fundamental matrix, giving expected steps
to absorption, absorption probabilities and their variances. It also computes
mean first-passage and mean recurrence times, hitting probabilities to
arbitrary target sets, the entropy rate, Kemeny's constant and the
time-reversed chain.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/markov"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AcceptanceRate(chain []float64) float64` | AcceptanceRate returns the fraction of consecutive samples in a scalar chain that changed value, an empirical proxy for the Metropolis acceptance… |
| `func ArgMax(v []float64) int` | ArgMax returns the index of the (first) maximal entry of v, or -1 if v is empty. |
| `func ArgMin(v []float64) int` | ArgMin returns the index of the (first) minimal entry of v, or -1 if v is empty. |
| `func Autocorrelation(chain []float64, lag int) float64` | Autocorrelation returns the sample autocorrelation of a scalar chain at the given lag (using the biased 1/N normalization). |
| `func AutocorrelationFunction(chain []float64, maxLag int) []float64` | AutocorrelationFunction returns the autocorrelation of the chain for lags 0 through maxLag inclusive. |
| `func BatchMeansVariance(chain []float64, numBatches int) float64` | BatchMeansVariance estimates the variance of the sample mean of a chain by the non-overlapping batch-means method with the given number of batches. |
| `func BhattacharyyaCoefficient(p, q []float64) float64` | BhattacharyyaCoefficient returns Σ_i sqrt(p_i q_i), a similarity measure in [0,1] for probability vectors. |
| `func BhattacharyyaDistance(p, q []float64) float64` | BhattacharyyaDistance returns -log of the Bhattacharyya coefficient. |
| `func Burnin(chain [][]float64, burn int) [][]float64` | Burnin returns the chain with the first burn samples discarded. |
| `func ChiSquareDistance(p, q []float64) float64` | ChiSquareDistance returns the Pearson chi-square distance Σ_i (p_i-q_i)²/q_i between probability vectors. |
| `func ColSums(a [][]float64) []float64` | ColSums returns the sum of each column of a. |
| `func Cols(a [][]float64) int` | Cols returns the number of columns of a (the length of its first row, or 0). |
| `func Column(chain [][]float64, k int) []float64` | Column extracts column k from a multidimensional chain (a slice of state vectors) as a scalar series, for use with the scalar diagnostics. |
| `func CopyMatrix(a [][]float64) [][]float64` | CopyMatrix returns a deep copy of a. |
| `func CopyVector(v []float64) []float64` | CopyVector returns a copy of v. |
| `func CovarianceMatrix(chain [][]float64) [][]float64` | CovarianceMatrix returns the unbiased sample covariance matrix of a multidimensional chain (rows are observations, columns are variables). |
| `func CrossEntropy(p, q []float64) float64` | CrossEntropy returns the cross entropy H(p, q) = -Σ_i p_i log q_i in nats. |
| `func CumulativeDistribution(p []float64) []float64` | CumulativeDistribution returns the cumulative sums of the (non-negative) probability vector p, normalized so the final entry is 1. |
| `func Determinant(a [][]float64) (float64, error)` | Determinant returns the determinant of the square matrix a. |
| `func Dot(a, b []float64) float64` | Dot returns the inner product of vectors a and b. |
| `func EffectiveSampleSize(chain []float64) float64` | EffectiveSampleSize estimates the effective sample size of a scalar chain as N / tau, where tau is the integrated autocorrelation time. |
| `func EuclideanDistance(p, q []float64) float64` | EuclideanDistance returns the L2 distance between p and q. |
| `func GelmanRubin(chains [][]float64) float64` | GelmanRubin returns the Gelman-Rubin potential scale reduction factor R-hat for several scalar chains (all of the same length m). |
| `func GibbsSampler(conditionals []func(x []float64, rng *rand.Rand) float64, x0 []float64, n int, rng *rand.Rand) [][]float64` | GibbsSampler runs a systematic-scan Gibbs sampler. |
| `func HadamardProduct(a, b [][]float64) [][]float64` | HadamardProduct returns the elementwise product of a and b. |
| `func HellingerDistance(p, q []float64) float64` | HellingerDistance returns the Hellinger distance between probability vectors p and q, in [0,1]. |
| `func Identity(n int) [][]float64` | Identity returns the n×n identity matrix. |
| `func IntegratedAutocorrelationTime(chain []float64, maxLag int) float64` | IntegratedAutocorrelationTime estimates the integrated autocorrelation time tau = 1 + 2 Σ_{k>=1} ρ(k), truncating the sum at the first lag where… |
| `func IsColumnStochastic(a [][]float64, tol float64) bool` | IsColumnStochastic reports whether a is square, non-negative, and every column sums to 1 within tol. |
| `func IsDoublyStochastic(a [][]float64, tol float64) bool` | IsDoublyStochastic reports whether a is both row- and column-stochastic within tol. |
| `func IsNonNegative(a [][]float64, tol float64) bool` | IsNonNegative reports whether every entry of a is >= -tol. |
| `func IsProbabilityVector(v []float64, tol float64) bool` | IsProbabilityVector reports whether v is non-negative and sums to 1 within tol. |
| `func IsSquare(a [][]float64) bool` | IsSquare reports whether a is a non-empty square matrix. |
| `func IsStochastic(a [][]float64, tol float64) bool` | IsStochastic reports whether a is (row) stochastic: square, non-negative, and every row sums to 1 within tol. |
| `func IsSubStochastic(a [][]float64, tol float64) bool` | IsSubStochastic reports whether a is non-negative and every row sums to at most 1+tol (rows may sum to less than 1). |
| `func IsSymmetricMatrix(a [][]float64, tol float64) bool` | IsSymmetricMatrix reports whether the square matrix a equals its transpose within tol. |
| `func JensenShannonDivergence(p, q []float64) float64` | JensenShannonDivergence returns the Jensen-Shannon divergence between p and q in nats: the symmetrized, smoothed average of the two KL divergences to… |
| `func KLDivergence(p, q []float64) float64` | KLDivergence returns the Kullback-Leibler divergence D(p‖q) in nats. |
| `func KroneckerProduct(a, b [][]float64) [][]float64` | KroneckerProduct returns the Kronecker product a⊗b. |
| `func L1Distance(p, q []float64) float64` | L1Distance returns the L1 (taxicab) distance between p and q. |
| `func LUDecompose(a [][]float64) (lu [][]float64, piv []int, sign float64, err error)` | LUDecompose returns the LU decomposition of the square matrix a with partial pivoting. |
| `func LogSumExp(v []float64) float64` | LogSumExp returns log(Σ exp(v_i)) computed stably by factoring out the maximum. |
| `func MatAdd(a, b [][]float64) [][]float64` | MatAdd returns a+b. |
| `func MatEqual(a, b [][]float64, tol float64) bool` | MatEqual reports whether a and b have the same shape and all corresponding entries differ by at most tol in absolute value. |
| `func MatInverse(a [][]float64) ([][]float64, error)` | MatInverse returns the inverse of the square matrix a, or ErrSingular if a is singular. |
| `func MatMul(a, b [][]float64) [][]float64` | MatMul returns the matrix product a·b. |
| `func MatNorm1(a [][]float64) float64` | MatNorm1 returns the induced 1-norm (maximum absolute column sum) of a. |
| `func MatNormFro(a [][]float64) float64` | MatNormFro returns the Frobenius norm of a. |
| `func MatNormInf(a [][]float64) float64` | MatNormInf returns the induced infinity-norm (maximum absolute row sum) of a. |
| `func MatPow(a [][]float64, k int) [][]float64` | MatPow returns a raised to the non-negative integer power k, computed by binary exponentiation. |
| `func MatScale(a [][]float64, s float64) [][]float64` | MatScale returns s·a. |
| `func MatSub(a, b [][]float64) [][]float64` | MatSub returns a-b. |
| `func MatVec(a [][]float64, v []float64) []float64` | MatVec returns the matrix-vector product a·v (v as a column). |
| `func MaxValue(v []float64) float64` | MaxValue returns the maximum entry of v, or NaN if v is empty. |
| `func MeanVector(chain [][]float64) []float64` | MeanVector returns the componentwise mean of a multidimensional chain (a slice of equal-length state vectors). |
| `func MetropolisHastings1D( logTarget func(float64) float64, propose func(x float64, rng *rand.Rand) float64, logProposal func(to, from float64) float64, x0 float64, n int, rng *rand.Rand, ) ([]float64, float64)` | MetropolisHastings1D runs a general (possibly asymmetric) Metropolis-Hastings sampler in one dimension. |
| `func MinValue(v []float64) float64` | MinValue returns the minimum entry of v, or NaN if v is empty. |
| `func Normalize(v []float64) []float64` | Normalize returns a copy of v scaled to sum to 1. |
| `func NormalizeRows(a [][]float64) [][]float64` | NormalizeRows returns a copy of a with each non-zero row scaled to sum to 1. |
| `func OuterProduct(u, v []float64) [][]float64` | OuterProduct returns the outer product u·vᵀ as a len(u)×len(v) matrix. |
| `func PointMass(n, i int) []float64` | PointMass returns the length-n probability vector placing all mass on state i. |
| `func RandomWalkMetropolis(logTarget func(float64) float64, x0, stepSize float64, n int, rng *rand.Rand) ([]float64, float64)` | RandomWalkMetropolis runs a symmetric random-walk Metropolis sampler in one dimension targeting the distribution whose log-density is logTarget… |
| `func RowSums(a [][]float64) []float64` | RowSums returns the sum of each row of a. |
| `func Rows(a [][]float64) int` | Rows returns the number of rows of a. |
| `func SampleCategorical(p []float64, rng *rand.Rand) int` | SampleCategorical draws an index from the categorical distribution p using rng. |
| `func SampleFromCDF(cdf []float64, u float64) int` | SampleFromCDF returns the smallest index i such that cdf`i` >= u, given a non-decreasing cumulative distribution cdf and a uniform draw u in [0,1). |
| `func SampleMean(chain []float64) float64` | SampleMean returns the mean of a scalar chain. |
| `func SampleStdDev(chain []float64) float64` | SampleStdDev returns the square root of the unbiased sample variance. |
| `func SampleVariance(chain []float64) float64` | SampleVariance returns the unbiased (n-1) sample variance of a scalar chain. |
| `func ShannonEntropy(p []float64) float64` | ShannonEntropy returns the Shannon entropy of the probability vector p in nats. |
| `func Softmax(v []float64) []float64` | Softmax returns the softmax (normalized exponentials) of v, computed stably. |
| `func SolveLinear(a [][]float64, b []float64) ([]float64, error)` | SolveLinear solves the linear system a·x = b for a square matrix a. |
| `func SolveMatrix(a, b [][]float64) ([][]float64, error)` | SolveMatrix solves a·X = B for a square matrix a and a matrix right-hand side B, returning X with the same number of columns as B. |
| `func Thin(chain [][]float64, step int) [][]float64` | Thin returns every step-th sample of the chain (step>=1). |
| `func TotalVariationDistance(p, q []float64) float64` | TotalVariationDistance returns the total-variation distance between two probability vectors p and q, defined as half their L1 distance. |
| `func Trace(a [][]float64) float64` | Trace returns the sum of the diagonal entries of the square matrix a. |
| `func Transpose(a [][]float64) [][]float64` | Transpose returns the transpose of the (possibly rectangular) matrix a. |
| `func UniformDistribution(n int) []float64` | UniformDistribution returns the uniform probability vector of length n. |
| `func VecAdd(a, b []float64) []float64` | VecAdd returns a+b, or nil if the lengths differ. |
| `func VecEqual(a, b []float64, tol float64) bool` | VecEqual reports whether a and b have equal length and every corresponding entry differs by at most tol. |
| `func VecMat(v []float64, a [][]float64) []float64` | VecMat returns the row-vector–matrix product v·a. |
| `func VecNorm1(v []float64) float64` | VecNorm1 returns the L1 norm (sum of absolute values) of v. |
| `func VecNorm2(v []float64) float64` | VecNorm2 returns the Euclidean (L2) norm of v. |
| `func VecNormInf(v []float64) float64` | VecNormInf returns the maximum absolute value among the entries of v. |
| `func VecScale(v []float64, s float64) []float64` | VecScale returns s·v. |
| `func VecSub(a, b []float64) []float64` | VecSub returns a-b, or nil if the lengths differ. |
| `func VecSum(v []float64) float64` | VecSum returns the sum of the entries of v. |
| `func Zeros(rows, cols int) [][]float64` | Zeros returns a rows×cols matrix filled with zeros. |

### Types

| Type | What it is |
| --- | --- |
| `AliasTable` | AliasTable implements Walker's alias method for O(1) sampling from a fixed discrete distribution. |
| `HMM` | HMM is a discrete hidden Markov model. |
| `MHResult` | MHResult holds the output of a Metropolis-Hastings run: the chain of samples and the fraction of proposed moves that were accepted. |
| `MarkovChain` | MarkovChain is a finite-state discrete-time Markov chain represented by a row-stochastic transition matrix P, where P`i``j` is the probability of… |

<details>
<summary><code>AliasTable</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewAliasTable(weights []float64) *AliasTable` | NewAliasTable builds an alias table for the discrete distribution weights (which need not be normalized). |
| `func (t *AliasTable) Len() int` | Len returns the number of outcomes in the alias table. |
| `func (t *AliasTable) Sample(rng *rand.Rand) int` | Sample draws an index from the alias table using rng. |

</details>

<details>
<summary><code>HMM</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewHMM(a, b [][]float64, pi []float64) (*HMM, error)` | NewHMM constructs a hidden Markov model from the transition matrix a, the emission matrix b, and the initial distribution pi. |
| `func (h *HMM) Backward(obs []int) (beta [][]float64, err error)` | Backward runs the unscaled backward algorithm and returns the beta trellis (T×N), where beta`t``i` = P(o_{t+1}..o_{T-1} \| x_t=i). |
| `func (h *HMM) BackwardScaled(obs []int, scale []float64) (beta [][]float64, err error)` | BackwardScaled runs the scaled backward algorithm using the scaling factors produced by ForwardScaled. |
| `func (h *HMM) BaumWelch(obs []int, maxIter int, tol float64) (*HMM, []float64, error)` | BaumWelch trains the model on a single observation sequence using the scaled Baum-Welch (EM) algorithm. |
| `func (h *HMM) BaumWelchMultiple(seqs [][]int, maxIter int, tol float64) (*HMM, []float64, error)` | BaumWelchMultiple trains the model on several observation sequences at once, accumulating the expected sufficient statistics across all sequences… |
| `func (h *HMM) Decode(obs []int) ([]int, float64, error)` | Decode is an alias for Viterbi: it returns the most likely hidden-state path and its log probability. |
| `func (h *HMM) EmissionMatrix() [][]float64` | EmissionMatrix returns a copy of the emission matrix B. |
| `func (h *HMM) EmissionProb(state, symbol int) float64` | EmissionProb returns the probability B`state``symbol` of emitting symbol in the given hidden state. |
| `func (h *HMM) ExpectedTransitionCounts(obs []int) ([][]float64, error)` | ExpectedTransitionCounts returns the matrix of expected numbers of i→j transitions over the sequence obs, Σ_t P(x_t=i, x_{t+1}=j \| obs), from the… |
| `func (h *HMM) Filter(obs []int) ([]float64, error)` | Filter returns the filtering distribution P(x_{T-1} \| o_0..o_{T-1}), i.e. |
| `func (h *HMM) Forward(obs []int) (alpha [][]float64, likelihood float64, err error)` | Forward runs the unscaled forward algorithm and returns the alpha trellis (T×N), where alpha`t``i` = P(o_0..o_t, x_t=i). |
| `func (h *HMM) ForwardScaled(obs []int) (alpha [][]float64, scale []float64, logLik float64, err error)` | ForwardScaled runs the numerically stable scaled forward algorithm. |
| `func (h *HMM) Generate(length int, rng *rand.Rand) (states, obs []int)` | Generate samples a state path and observation sequence of the given length from the model using rng. |
| `func (h *HMM) InitialDistribution() []float64` | InitialDistribution returns a copy of the initial state distribution Pi. |
| `func (h *HMM) Likelihood(obs []int) (float64, error)` | Likelihood returns P(obs), the total probability of the observation sequence under the model. |
| `func (h *HMM) LogLikelihood(obs []int) (float64, error)` | LogLikelihood returns log P(obs), computed stably via the scaled forward algorithm. |
| `func (h *HMM) MostLikelyStates(obs []int) ([]int, error)` | MostLikelyStates returns the sequence of individually most probable states, argmax_i P(x_t=i \| obs) for each t. |
| `func (h *HMM) NumStates() int` | NumStates returns the number of hidden states N. |
| `func (h *HMM) NumSymbols() int` | NumSymbols returns the number of observable symbols M. |
| `func (h *HMM) PosteriorMarginals(obs []int) ([][]float64, error)` | PosteriorMarginals returns the gamma matrix (T×N) of smoothed state posteriors, gamma`t``i` = P(x_t=i \| obs), using the scaled forward-backward… |
| `func (h *HMM) PredictNextObservationDistribution(obs []int) ([]float64, error)` | PredictNextObservationDistribution returns the predictive distribution over the next observed symbol given obs. |
| `func (h *HMM) PredictNextStateDistribution(obs []int) ([]float64, error)` | PredictNextStateDistribution returns the distribution over the hidden state at time T (one step past the last observation) given obs. |
| `func (h *HMM) StationaryDistribution() ([]float64, error)` | StationaryDistribution returns the stationary distribution of the HMM's hidden-state transition matrix A. |
| `func (h *HMM) TransitionMatrix() [][]float64` | TransitionMatrix returns a copy of the hidden-state transition matrix A. |
| `func (h *HMM) TransitionProb(i, j int) float64` | TransitionProb returns the hidden-state transition probability A`i``j`. |
| `func (h *HMM) Viterbi(obs []int) (path []int, logProb float64, err error)` | Viterbi returns the most likely hidden-state sequence for obs and the log probability of that path, using the log-domain Viterbi algorithm. |

</details>

<details>
<summary><code>MHResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func RandomWalkMetropolisND(logTarget func([]float64) float64, x0 []float64, stepSize float64, n int, rng *rand.Rand) *MHResult` | RandomWalkMetropolisND runs a random-walk Metropolis sampler in d dimensions with an isotropic Gaussian proposal of the given step size. |

</details>

<details>
<summary><code>MarkovChain</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewMarkovChain(p [][]float64) (*MarkovChain, error)` | NewMarkovChain builds a MarkovChain from the row-stochastic matrix p. |
| `func NewMarkovChainFromCounts(counts [][]float64) (*MarkovChain, error)` | NewMarkovChainFromCounts builds a MarkovChain by row-normalizing a matrix of transition counts. |
| `func NewMarkovChainUnchecked(p [][]float64) (*MarkovChain, error)` | NewMarkovChainUnchecked builds a MarkovChain from p without validating that it is stochastic. |
| `func (c *MarkovChain) AbsorbingStates() []int` | AbsorbingStates returns the sorted list of absorbing states (states i with P`i``i` = 1). |
| `func (c *MarkovChain) AbsorptionProbabilities() (b [][]float64, transient, absorbing []int, err error)` | AbsorptionProbabilities returns the matrix B = N·R of absorption probabilities: B`i``k` is the probability that a chain started in transient state i… |
| `func (c *MarkovChain) AbsorptionProbability(from, target int) (float64, error)` | AbsorptionProbability returns the probability that the chain started in state from is eventually absorbed in the absorbing state target. |
| `func (c *MarkovChain) Accessible(i, j int) bool` | Accessible reports whether state j is accessible from state i, i.e. |
| `func (c *MarkovChain) AdjacencyMatrix(tol float64) [][]bool` | AdjacencyMatrix returns the boolean adjacency matrix of the chain's transition graph: entry `i``j` is true iff P`i``j` > tol (state j is reachable… |
| `func (c *MarkovChain) CanonicalForm() (reordered [][]float64, transient, absorbing []int)` | CanonicalForm reorders the states of an absorbing chain into transient states followed by absorbing states and returns the reordered transition… |
| `func (c *MarkovChain) ClassifyState(i int) string` | ClassifyState returns a short label for state i: "absorbing", "recurrent", or "transient". |
| `func (c *MarkovChain) ClassifyStates() map[int]string` | ClassifyStates returns a map from every state to its classification label, one of "absorbing", "recurrent", or "transient". |
| `func (c *MarkovChain) Communicates(i, j int) bool` | Communicates reports whether states i and j communicate, i.e. |
| `func (c *MarkovChain) CommunicatingClasses() [][]int` | CommunicatingClasses returns the communicating classes (strongly connected components of the transition graph) as a slice of sorted state slices. |
| `func (c *MarkovChain) EntropyRate(pi []float64) float64` | EntropyRate returns the entropy rate (in nats) of the chain under the given stationary distribution pi: H = -Σ_i π_i Σ_j P_ij log P_ij. |
| `func (c *MarkovChain) ExpectedHittingTime(from int, target []int) (float64, error)` | ExpectedHittingTime returns the expected number of steps to first reach the target set, starting from state from (0 if from is already in the target). |
| `func (c *MarkovChain) ExpectedHittingTimes(target []int) ([]float64, error)` | ExpectedHittingTimes returns, for every state i, the expected number of steps to first reach the target set (0 for states already in the target). |
| `func (c *MarkovChain) ExpectedReward(reward []float64) (float64, error)` | ExpectedReward returns the long-run average reward per step under the stationary distribution, given a per-state reward vector: Σ_i π_i r_i. |
| `func (c *MarkovChain) ExpectedStepsToAbsorption() ([]float64, []int, error)` | ExpectedStepsToAbsorption returns, for each transient state, the expected number of steps until the chain is absorbed. |
| `func (c *MarkovChain) ExpectedVisits(from int) ([]float64, []int, error)` | ExpectedVisits returns the expected number of visits to each transient state before absorption, starting from transient state i (given in original… |
| `func (c *MarkovChain) FundamentalMatrix() ([][]float64, []int, error)` | FundamentalMatrix returns the fundamental matrix N = (I - Q)^{-1} of an absorbing chain, indexed by the transient states in ascending order. |
| `func (c *MarkovChain) HittingProbabilities(target []int) ([]float64, error)` | HittingProbabilities returns, for every state i, the probability that the chain starting at i ever reaches the target set. |
| `func (c *MarkovChain) HittingProbability(from int, target []int) (float64, error)` | HittingProbability returns the probability that the chain starting in state from ever reaches the target set. |
| `func (c *MarkovChain) IsAbsorbing() bool` | IsAbsorbing reports whether the chain is an absorbing chain: it has at least one absorbing state and from every state some absorbing state is… |
| `func (c *MarkovChain) IsAbsorbingState(i int) bool` | IsAbsorbingState reports whether state i is absorbing (P`i``i` = 1). |
| `func (c *MarkovChain) IsAperiodic() bool` | IsAperiodic reports whether the chain is aperiodic, meaning every state has period 1. |
| `func (c *MarkovChain) IsClosedClass(states []int) bool` | IsClosedClass reports whether the given set of states is closed: no transition leaves the set. |
| `func (c *MarkovChain) IsErgodic() bool` | IsErgodic reports whether the chain is ergodic (irreducible and aperiodic), which for a finite chain is equivalent to being regular. |
| `func (c *MarkovChain) IsIrreducible() bool` | IsIrreducible reports whether the chain is irreducible: all states form a single communicating class. |
| `func (c *MarkovChain) IsRegular() bool` | IsRegular reports whether the chain is regular: some power P^k has all strictly positive entries. |
| `func (c *MarkovChain) IsReversible(pi []float64, tol float64) bool` | IsReversible reports whether the chain satisfies the detailed-balance equations π_i P_ij = π_j P_ji for all i,j within tol, for the given… |
| `func (c *MarkovChain) KemenyConstant() (float64, error)` | KemenyConstant returns Kemeny's constant of an ergodic chain: the expected number of steps to reach a state drawn from the stationary distribution,… |
| `func (c *MarkovChain) LimitingMatrix(tol float64, iters int) ([][]float64, error)` | LimitingMatrix returns P^k for large k, approximating the limiting matrix whose rows are the stationary distribution (for an ergodic chain). |
| `func (c *MarkovChain) LogProbabilityOfPath(states []int) float64` | LogProbabilityOfPath returns the natural logarithm of ProbabilityOfPath, computed as a sum of logs for numerical stability. |
| `func (c *MarkovChain) Matrix() [][]float64` | Matrix returns a copy of the chain's transition matrix. |
| `func (c *MarkovChain) MeanFirstPassageTime(i, j int) (float64, error)` | MeanFirstPassageTime returns the expected number of steps to first reach state j starting from state i in an ergodic chain. |
| `func (c *MarkovChain) MeanFirstPassageTimes() ([][]float64, error)` | MeanFirstPassageTimes returns the matrix M of mean first-passage times for an ergodic chain: M`i``j` is the expected number of steps to first reach… |
| `func (c *MarkovChain) MeanRecurrenceTimes() ([]float64, error)` | MeanRecurrenceTimes returns the vector of mean recurrence times of an ergodic chain, m_jj = 1/π_j. |
| `func (c *MarkovChain) MixingTime(epsilon float64, maxSteps int) (int, error)` | MixingTime returns the smallest number of steps k such that, from every starting state, the total-variation distance between the k-step distribution… |
| `func (c *MarkovChain) N() int` | N returns the number of states of the chain. |
| `func (c *MarkovChain) NStep(k int) [][]float64` | NStep returns the k-step transition matrix P^k. |
| `func (c *MarkovChain) NStepProb(i, j, k int) float64` | NStepProb returns the probability of being in state j after k steps starting from state i, i.e. |
| `func (c *MarkovChain) NextDistribution(dist []float64) []float64` | NextDistribution returns the distribution after one step given the current distribution dist, i.e. |
| `func (c *MarkovChain) NumClasses() int` | NumClasses returns the number of communicating classes. |
| `func (c *MarkovChain) Period(i int) int` | Period returns the period of state i: the greatest common divisor of the lengths of all closed paths (returns) from i to itself. |
| `func (c *MarkovChain) ProbabilityOfPath(states []int) float64` | ProbabilityOfPath returns the probability of the state sequence conditional on starting in states[0], i.e. |
| `func (c *MarkovChain) ReachabilityMatrix() [][]bool` | ReachabilityMatrix returns the boolean matrix whose `i``j` entry is true iff state j is accessible from state i (transitive closure of the one-step… |
| `func (c *MarkovChain) ReachableFrom(i int) []int` | ReachableFrom returns the sorted list of states reachable from state i (including i itself). |
| `func (c *MarkovChain) RecurrentClasses() [][]int` | RecurrentClasses returns the communicating classes that are recurrent (closed under the dynamics), each as a sorted slice. |
| `func (c *MarkovChain) RecurrentStates() []int` | RecurrentStates returns the sorted list of recurrent states: those belonging to a closed communicating class. |
| `func (c *MarkovChain) Reverse(pi []float64) *MarkovChain` | Reverse returns the time-reversed chain with respect to the distribution pi, whose transition matrix is P*_ij = π_j P_ji / π_i. |
| `func (c *MarkovChain) Reverse01() *MarkovChain` | Reverse01 returns a chain whose transition graph is the edge-reversal of this chain's positive-probability graph, with rows normalized to be… |
| `func (c *MarkovChain) Row(i int) []float64` | Row returns a copy of row i of the transition matrix (the conditional distribution of the next state given current state i). |
| `func (c *MarkovChain) Simulate(start, steps int, rng *rand.Rand) []int` | Simulate returns a trajectory of the chain of the given number of steps (inclusive of the starting state, so the slice has length steps+1), starting… |
| `func (c *MarkovChain) SimulateFrom(init []float64, steps int, rng *rand.Rand) []int` | SimulateFrom returns a trajectory whose starting state is sampled from the initial distribution init. |
| `func (c *MarkovChain) StationaryDistribution() ([]float64, error)` | StationaryDistribution returns the stationary distribution of the chain. |
| `func (c *MarkovChain) StationaryEntropy() (float64, error)` | StationaryEntropy returns the Shannon entropy (in nats) of the chain's stationary distribution. |
| `func (c *MarkovChain) StationaryPower(tol float64, maxIter int) ([]float64, error)` | StationaryPower approximates a stationary distribution by power iteration starting from the uniform distribution, stopping when successive iterates… |
| `func (c *MarkovChain) StationarySolve() ([]float64, error)` | StationarySolve computes a stationary distribution π satisfying πP = π and Σπ = 1 by solving the linear system exactly. |
| `func (c *MarkovChain) StepDistribution(init []float64, k int) []float64` | StepDistribution returns the distribution after k steps starting from the initial distribution init, i.e. |
| `func (c *MarkovChain) SubmatrixQ() ([][]float64, []int)` | SubmatrixQ returns the transient-to-transient block Q of the canonical form (indexed by the transient states in ascending order). |
| `func (c *MarkovChain) SubmatrixR() (r [][]float64, transient, absorbing []int)` | SubmatrixR returns the transient-to-absorbing block R of the canonical form, together with the transient and absorbing index lists. |
| `func (c *MarkovChain) TVDistanceToStationary(init []float64, k int) (float64, error)` | TVDistanceToStationary returns the total-variation distance between the distribution after k steps (starting from init) and the stationary… |
| `func (c *MarkovChain) TransientClasses() [][]int` | TransientClasses returns the communicating classes that are transient (not closed), each as a sorted slice. |
| `func (c *MarkovChain) TransientStates() []int` | TransientStates returns the sorted list of transient states: those not in any closed communicating class. |
| `func (c *MarkovChain) TransitionProb(i, j int) float64` | TransitionProb returns P`i``j`, the one-step probability of moving from state i to state j. |
| `func (c *MarkovChain) VarianceStepsToAbsorption() ([]float64, []int, error)` | VarianceStepsToAbsorption returns the variance of the number of steps to absorption for each transient state, using the standard formula (2N - I)t -… |

</details>

### Constants

`DefaultTol`

### Variables

`ErrNotSquare`, `ErrDimMismatch`, `ErrSingular`, `ErrEmpty`, `ErrNotStochastic`, `ErrNoConvergence`, `ErrNotAbsorbing`, `ErrNotErgodic`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/markov).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
