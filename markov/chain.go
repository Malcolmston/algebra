package markov

import (
	"math"
	"math/rand"
)

// MarkovChain is a finite-state discrete-time Markov chain represented by a
// row-stochastic transition matrix P, where P[i][j] is the probability of
// moving from state i to state j in one step.
type MarkovChain struct {
	p [][]float64
	n int
}

// NewMarkovChain builds a MarkovChain from the row-stochastic matrix p. The
// matrix is copied. It returns ErrNotStochastic if p is not square,
// non-negative, and row-summing to 1 within DefaultTol.
func NewMarkovChain(p [][]float64) (*MarkovChain, error) {
	if !IsSquare(p) {
		return nil, ErrNotSquare
	}
	if !IsStochastic(p, 1e-9) {
		return nil, ErrNotStochastic
	}
	return &MarkovChain{p: CopyMatrix(p), n: len(p)}, nil
}

// NewMarkovChainUnchecked builds a MarkovChain from p without validating that
// it is stochastic. The matrix is copied. It is the caller's responsibility to
// ensure p is a valid transition matrix; use this only when p is known good
// (for example produced by NormalizeRows). It still requires p to be square.
func NewMarkovChainUnchecked(p [][]float64) (*MarkovChain, error) {
	if !IsSquare(p) {
		return nil, ErrNotSquare
	}
	return &MarkovChain{p: CopyMatrix(p), n: len(p)}, nil
}

// NewMarkovChainFromCounts builds a MarkovChain by row-normalizing a matrix of
// transition counts. A row that is all zeros is turned into a self-loop (the
// state becomes absorbing) so the result is always stochastic.
func NewMarkovChainFromCounts(counts [][]float64) (*MarkovChain, error) {
	if !IsSquare(counts) {
		return nil, ErrNotSquare
	}
	p := CopyMatrix(counts)
	for i := range p {
		var s float64
		for _, x := range p[i] {
			if x < 0 {
				return nil, ErrDimMismatch
			}
			s += x
		}
		if s == 0 {
			p[i][i] = 1
		} else {
			for j := range p[i] {
				p[i][j] /= s
			}
		}
	}
	return &MarkovChain{p: p, n: len(p)}, nil
}

// N returns the number of states of the chain.
func (c *MarkovChain) N() int { return c.n }

// Matrix returns a copy of the chain's transition matrix.
func (c *MarkovChain) Matrix() [][]float64 { return CopyMatrix(c.p) }

// TransitionProb returns P[i][j], the one-step probability of moving from state
// i to state j. It returns NaN if an index is out of range.
func (c *MarkovChain) TransitionProb(i, j int) float64 {
	if i < 0 || i >= c.n || j < 0 || j >= c.n {
		return math.NaN()
	}
	return c.p[i][j]
}

// Row returns a copy of row i of the transition matrix (the conditional
// distribution of the next state given current state i).
func (c *MarkovChain) Row(i int) []float64 {
	if i < 0 || i >= c.n {
		return nil
	}
	return CopyVector(c.p[i])
}

// NStep returns the k-step transition matrix P^k. NStep(0) is the identity.
func (c *MarkovChain) NStep(k int) [][]float64 {
	return MatPow(c.p, k)
}

// NStepProb returns the probability of being in state j after k steps starting
// from state i, i.e. (P^k)[i][j].
func (c *MarkovChain) NStepProb(i, j, k int) float64 {
	pk := MatPow(c.p, k)
	if pk == nil || i < 0 || i >= c.n || j < 0 || j >= c.n {
		return math.NaN()
	}
	return pk[i][j]
}

// NextDistribution returns the distribution after one step given the current
// distribution dist, i.e. dist·P. It returns nil if len(dist) != N.
func (c *MarkovChain) NextDistribution(dist []float64) []float64 {
	if len(dist) != c.n {
		return nil
	}
	return VecMat(dist, c.p)
}

// StepDistribution returns the distribution after k steps starting from the
// initial distribution init, i.e. init·P^k. It returns nil if len(init) != N or
// k is negative.
func (c *MarkovChain) StepDistribution(init []float64, k int) []float64 {
	if len(init) != c.n || k < 0 {
		return nil
	}
	d := CopyVector(init)
	for step := 0; step < k; step++ {
		d = VecMat(d, c.p)
	}
	return d
}

// StationaryPower approximates a stationary distribution by power iteration
// starting from the uniform distribution, stopping when successive iterates
// differ by at most tol in L1 norm or after maxIter iterations. It returns
// ErrNoConvergence if the tolerance is not reached. For chains with more than
// one stationary distribution (reducible chains) the result depends on the
// start vector.
func (c *MarkovChain) StationaryPower(tol float64, maxIter int) ([]float64, error) {
	if c.n == 0 {
		return nil, ErrEmpty
	}
	if maxIter <= 0 {
		maxIter = 10000
	}
	v := UniformDistribution(c.n)
	for iter := 0; iter < maxIter; iter++ {
		next := VecMat(v, c.p)
		// Renormalize to guard against tiny drift.
		next = Normalize(next)
		if L1Distance(v, next) <= tol {
			return next, nil
		}
		v = next
	}
	return v, ErrNoConvergence
}

// StationarySolve computes a stationary distribution π satisfying πP = π and
// Σπ = 1 by solving the linear system exactly. For an irreducible chain the
// stationary distribution is unique; for reducible chains this returns one
// particular solution (which may contain negative entries if none is strictly
// positive, in which case interpret with care).
func (c *MarkovChain) StationarySolve() ([]float64, error) {
	n := c.n
	if n == 0 {
		return nil, ErrEmpty
	}
	// Solve (P^T - I) π = 0 with the last equation replaced by Σπ = 1.
	a := make([][]float64, n)
	b := make([]float64, n)
	for i := 0; i < n; i++ {
		a[i] = make([]float64, n)
		for j := 0; j < n; j++ {
			a[i][j] = c.p[j][i] // (P^T)[i][j]
			if i == j {
				a[i][j] -= 1
			}
		}
	}
	// Replace last row with the normalization constraint.
	for j := 0; j < n; j++ {
		a[n-1][j] = 1
	}
	b[n-1] = 1
	pi, err := SolveLinear(a, b)
	if err != nil {
		return nil, err
	}
	return pi, nil
}

// StationaryDistribution returns the stationary distribution of the chain. It
// solves the linear system exactly and, on success, normalizes the result. It
// is the recommended entry point for well-behaved (irreducible) chains.
func (c *MarkovChain) StationaryDistribution() ([]float64, error) {
	pi, err := c.StationarySolve()
	if err != nil {
		return nil, err
	}
	return Normalize(pi), nil
}

// LimitingMatrix returns P^k for large k, approximating the limiting matrix
// whose rows are the stationary distribution (for an ergodic chain). It squares
// the matrix repeatedly until successive powers agree within tol (max-norm) or
// iters doublings are reached.
func (c *MarkovChain) LimitingMatrix(tol float64, iters int) ([][]float64, error) {
	if c.n == 0 {
		return nil, ErrEmpty
	}
	if iters <= 0 {
		iters = 60
	}
	m := CopyMatrix(c.p)
	for i := 0; i < iters; i++ {
		next := MatMul(m, m)
		if MatEqual(m, next, tol) {
			return next, nil
		}
		m = next
	}
	return m, ErrNoConvergence
}

// IsReversible reports whether the chain satisfies the detailed-balance
// equations π_i P_ij = π_j P_ji for all i,j within tol, for the given
// distribution pi.
func (c *MarkovChain) IsReversible(pi []float64, tol float64) bool {
	if len(pi) != c.n {
		return false
	}
	for i := 0; i < c.n; i++ {
		for j := 0; j < c.n; j++ {
			if math.Abs(pi[i]*c.p[i][j]-pi[j]*c.p[j][i]) > tol {
				return false
			}
		}
	}
	return true
}

// Reverse returns the time-reversed chain with respect to the distribution pi,
// whose transition matrix is P*_ij = π_j P_ji / π_i. States with π_i = 0 are
// given a self-loop. It returns nil if len(pi) != N.
func (c *MarkovChain) Reverse(pi []float64) *MarkovChain {
	if len(pi) != c.n {
		return nil
	}
	q := make([][]float64, c.n)
	for i := 0; i < c.n; i++ {
		q[i] = make([]float64, c.n)
		if pi[i] == 0 {
			q[i][i] = 1
			continue
		}
		for j := 0; j < c.n; j++ {
			q[i][j] = pi[j] * c.p[j][i] / pi[i]
		}
	}
	return &MarkovChain{p: NormalizeRows(q), n: c.n}
}

// EntropyRate returns the entropy rate (in nats) of the chain under the given
// stationary distribution pi: H = -Σ_i π_i Σ_j P_ij log P_ij. It returns NaN if
// len(pi) != N.
func (c *MarkovChain) EntropyRate(pi []float64) float64 {
	if len(pi) != c.n {
		return math.NaN()
	}
	var h float64
	for i := 0; i < c.n; i++ {
		var row float64
		for j := 0; j < c.n; j++ {
			if c.p[i][j] > 0 {
				row -= c.p[i][j] * math.Log(c.p[i][j])
			}
		}
		h += pi[i] * row
	}
	return h
}

// KemenyConstant returns Kemeny's constant of an ergodic chain: the expected
// number of steps to reach a state drawn from the stationary distribution,
// starting from any fixed state (the value is independent of the start state).
// It equals Σ_j π_j m_ij for any i, where m are the mean first-passage times.
func (c *MarkovChain) KemenyConstant() (float64, error) {
	pi, err := c.StationaryDistribution()
	if err != nil {
		return 0, err
	}
	m, err := c.MeanFirstPassageTimes()
	if err != nil {
		return 0, err
	}
	// Use row 0; sum over j != 0 of pi_j * m_0j (diagonal m_00 excluded because
	// the standard Kemeny constant sums first-passage times to distinct states).
	var k float64
	for j := 0; j < c.n; j++ {
		if j == 0 {
			continue
		}
		k += pi[j] * m[0][j]
	}
	return k, nil
}

// Simulate returns a trajectory of the chain of the given number of steps
// (inclusive of the starting state, so the slice has length steps+1), starting
// from state start and using rng for the random transitions. It returns nil if
// start is out of range or steps is negative.
func (c *MarkovChain) Simulate(start, steps int, rng *rand.Rand) []int {
	if start < 0 || start >= c.n || steps < 0 || rng == nil {
		return nil
	}
	path := make([]int, steps+1)
	path[0] = start
	cur := start
	for t := 1; t <= steps; t++ {
		cur = sampleRow(c.p[cur], rng)
		path[t] = cur
	}
	return path
}

// SimulateFrom returns a trajectory whose starting state is sampled from the
// initial distribution init. The returned slice has length steps+1. It returns
// nil on dimension mismatch.
func (c *MarkovChain) SimulateFrom(init []float64, steps int, rng *rand.Rand) []int {
	if len(init) != c.n || rng == nil {
		return nil
	}
	start := SampleCategorical(init, rng)
	return c.Simulate(start, steps, rng)
}

// sampleRow draws a state index from the categorical distribution row.
func sampleRow(row []float64, rng *rand.Rand) int {
	u := rng.Float64()
	var cum float64
	for j, p := range row {
		cum += p
		if u < cum {
			return j
		}
	}
	return len(row) - 1
}
