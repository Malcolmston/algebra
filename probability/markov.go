package probability

import "math"

// MarkovChain is a finite-state discrete-time Markov chain represented by a
// row-stochastic transition matrix P, where P[i][j] is the probability of moving
// from state i to state j in one step. Each row is non-negative and sums to one.
type MarkovChain struct {
	// P is the square row-stochastic transition matrix.
	P [][]float64
}

// NewMarkovChain builds a MarkovChain from a transition matrix, validating that
// it is square, non-negative, and row-stochastic (each row sums to one within
// [probabilityTol]). The matrix is copied. It returns an error otherwise.
func NewMarkovChain(p [][]float64) (MarkovChain, error) {
	n := len(p)
	if n == 0 {
		return MarkovChain{}, probabilityErrorf("NewMarkovChain: empty matrix")
	}
	cp := make([][]float64, n)
	for i := range p {
		if len(p[i]) != n {
			return MarkovChain{}, probabilityErrorf("NewMarkovChain: row %d has length %d, want %d", i, len(p[i]), n)
		}
		cp[i] = make([]float64, n)
		rowSum := 0.0
		for j, v := range p[i] {
			if v < 0 || math.IsNaN(v) || math.IsInf(v, 0) {
				return MarkovChain{}, probabilityErrorf("NewMarkovChain: invalid entry %g at (%d,%d)", v, i, j)
			}
			cp[i][j] = v
			rowSum += v
		}
		if probabilityAbs(rowSum-1) > probabilityTol {
			return MarkovChain{}, probabilityErrorf("NewMarkovChain: row %d sums to %g, not 1", i, rowSum)
		}
	}
	return MarkovChain{P: cp}, nil
}

// Size returns the number of states in the chain.
func (m MarkovChain) Size() int { return len(m.P) }

// Step advances a distribution over states by one transition, returning the row
// vector dist·P. dist must have length equal to the number of states. It returns
// an error on a length mismatch.
func (m MarkovChain) Step(dist []float64) ([]float64, error) {
	if len(dist) != len(m.P) {
		return nil, probabilityErrorf("Step: distribution length %d != state count %d", len(dist), len(m.P))
	}
	return probabilityVecMat(dist, m.P), nil
}

// NStep returns the n-step transition matrix P^n, whose (i, j) entry is the
// probability of moving from state i to state j in exactly n steps.
// NStep(0) is the identity. It returns an error for negative n.
func (m MarkovChain) NStep(n int) ([][]float64, error) {
	if n < 0 {
		return nil, probabilityErrorf("NStep: negative n=%d", n)
	}
	return probabilityMatPow(m.P, n), nil
}

// DistributionAfter returns the state distribution after n steps starting from
// the initial distribution, i.e. initial·P^n. initial must have length equal to
// the number of states. It returns an error on a length mismatch or negative n.
func (m MarkovChain) DistributionAfter(initial []float64, n int) ([]float64, error) {
	if len(initial) != len(m.P) {
		return nil, probabilityErrorf("DistributionAfter: distribution length %d != state count %d", len(initial), len(m.P))
	}
	if n < 0 {
		return nil, probabilityErrorf("DistributionAfter: negative n=%d", n)
	}
	pn := probabilityMatPow(m.P, n)
	return probabilityVecMat(initial, pn), nil
}

// StationaryDistribution returns a stationary distribution π satisfying πP = π
// and Σ π = 1. For an irreducible chain the stationary distribution is unique.
// It solves the linear system (P^T - I)π = 0 with a normalization constraint via
// Gaussian elimination. It returns an error if the system is singular (e.g. a
// reducible chain without a unique stationary distribution).
func (m MarkovChain) StationaryDistribution() ([]float64, error) {
	n := len(m.P)
	// Build A where A[j][i] = P[i][j] - (i==j), i.e. (P^T - I).
	a := make([][]float64, n)
	for j := 0; j < n; j++ {
		a[j] = make([]float64, n)
		for i := 0; i < n; i++ {
			a[j][i] = m.P[i][j]
			if i == j {
				a[j][i] -= 1
			}
		}
	}
	// Replace the last equation with the normalization Σ π = 1.
	for i := 0; i < n; i++ {
		a[n-1][i] = 1
	}
	b := make([]float64, n)
	b[n-1] = 1
	x, err := probabilitySolve(a, b)
	if err != nil {
		return nil, probabilityErrorf("StationaryDistribution: %v", err)
	}
	// Guard against tiny negative round-off.
	for i := range x {
		if x[i] < 0 && x[i] > -probabilityTol {
			x[i] = 0
		}
	}
	return x, nil
}

// MeanRecurrenceTimes returns the vector of mean recurrence times of an
// irreducible chain, where entry i is 1/π_i and π is the stationary
// distribution. It returns an error if the stationary distribution cannot be
// computed or has a zero component.
func (m MarkovChain) MeanRecurrenceTimes() ([]float64, error) {
	pi, err := m.StationaryDistribution()
	if err != nil {
		return nil, probabilityErrorf("MeanRecurrenceTimes: %v", err)
	}
	out := make([]float64, len(pi))
	for i, p := range pi {
		if p <= 0 {
			return nil, probabilityErrorf("MeanRecurrenceTimes: zero stationary probability at state %d", i)
		}
		out[i] = 1 / p
	}
	return out, nil
}

// Reachable reports whether state j is reachable from state i in zero or more
// steps (a state is always reachable from itself).
func (m MarkovChain) Reachable(i, j int) bool {
	n := len(m.P)
	if i < 0 || i >= n || j < 0 || j >= n {
		return false
	}
	seen := make([]bool, n)
	stack := []int{i}
	seen[i] = true
	for len(stack) > 0 {
		s := stack[len(stack)-1]
		stack = stack[:len(stack)-1]
		if s == j {
			return true
		}
		for t := 0; t < n; t++ {
			if m.P[s][t] > 0 && !seen[t] {
				seen[t] = true
				stack = append(stack, t)
			}
		}
	}
	return false
}

// IsIrreducible reports whether the chain is irreducible, i.e. every state is
// reachable from every other state.
func (m MarkovChain) IsIrreducible() bool {
	n := len(m.P)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if !m.Reachable(i, j) {
				return false
			}
		}
	}
	return true
}

// IsRegular reports whether the chain is regular (primitive): some power P^k has
// all strictly positive entries. A regular chain is irreducible and aperiodic
// and has a unique limiting distribution. By Wielandt's bound it suffices to
// check powers up to (n-1)^2 + 1.
func (m MarkovChain) IsRegular() bool {
	n := len(m.P)
	limit := (n-1)*(n-1) + 1
	power := probabilityCopyMatrix(m.P)
	for k := 1; k <= limit; k++ {
		allPos := true
		for i := 0; i < n && allPos; i++ {
			for j := 0; j < n; j++ {
				if power[i][j] <= 0 {
					allPos = false
					break
				}
			}
		}
		if allPos {
			return true
		}
		power = probabilityMatMul(power, m.P)
	}
	return false
}

// AbsorbingStates returns the indices of the absorbing states, those i with
// P[i][i] equal to one (and hence no probability of leaving), in ascending
// order.
func (m MarkovChain) AbsorbingStates() []int {
	var out []int
	for i := range m.P {
		if probabilityAbs(m.P[i][i]-1) <= probabilityTol {
			out = append(out, i)
		}
	}
	return out
}

// TransientStates returns the indices of the non-absorbing states in ascending
// order. For an absorbing chain these are exactly the transient states.
func (m MarkovChain) TransientStates() []int {
	var out []int
	for i := range m.P {
		if probabilityAbs(m.P[i][i]-1) > probabilityTol {
			out = append(out, i)
		}
	}
	return out
}

// IsAbsorbing reports whether the chain is an absorbing Markov chain: it has at
// least one absorbing state and every state can reach an absorbing state.
func (m MarkovChain) IsAbsorbing() bool {
	abs := m.AbsorbingStates()
	if len(abs) == 0 {
		return false
	}
	absSet := make(map[int]bool, len(abs))
	for _, a := range abs {
		absSet[a] = true
	}
	for i := range m.P {
		reaches := false
		for _, a := range abs {
			if m.Reachable(i, a) {
				reaches = true
				break
			}
		}
		if !reaches {
			return false
		}
	}
	return true
}

// FundamentalMatrix returns the fundamental matrix N = (I - Q)^{-1} of an
// absorbing chain, where Q is the transient-to-transient submatrix of P. Rows
// and columns are indexed in the order returned by [MarkovChain.TransientStates].
// Entry N[i][j] is the expected number of visits to transient state j before
// absorption, starting from transient state i. It returns an error if the chain
// is not absorbing.
func (m MarkovChain) FundamentalMatrix() ([][]float64, error) {
	trans := m.TransientStates()
	if !m.IsAbsorbing() {
		return nil, probabilityErrorf("FundamentalMatrix: chain is not absorbing")
	}
	t := len(trans)
	// Build I - Q.
	imq := make([][]float64, t)
	for a := 0; a < t; a++ {
		imq[a] = make([]float64, t)
		for b := 0; b < t; b++ {
			imq[a][b] = -m.P[trans[a]][trans[b]]
			if a == b {
				imq[a][b] += 1
			}
		}
	}
	n, err := probabilityInverse(imq)
	if err != nil {
		return nil, probabilityErrorf("FundamentalMatrix: %v", err)
	}
	return n, nil
}

// ExpectedStepsToAbsorption returns, for each state, the expected number of
// steps until the chain is absorbed. The result is a full-length vector indexed
// by state: absorbing states have value zero and transient states hold the
// expected number of steps to absorption (the corresponding row sum of the
// fundamental matrix). It returns an error if the chain is not absorbing.
func (m MarkovChain) ExpectedStepsToAbsorption() ([]float64, error) {
	n, err := m.FundamentalMatrix()
	if err != nil {
		return nil, probabilityErrorf("ExpectedStepsToAbsorption: %v", err)
	}
	trans := m.TransientStates()
	out := make([]float64, len(m.P))
	for a := range trans {
		s := 0.0
		for b := range n[a] {
			s += n[a][b]
		}
		out[trans[a]] = s
	}
	return out, nil
}

// AbsorptionProbabilities returns the matrix B = N·R of absorption
// probabilities, where N is the fundamental matrix and R is the
// transient-to-absorbing submatrix of P. Rows are indexed in the order of
// [MarkovChain.TransientStates] and columns in the order of
// [MarkovChain.AbsorbingStates]; B[i][k] is the probability that a chain started
// in transient state i is eventually absorbed in absorbing state k. It returns
// an error if the chain is not absorbing.
func (m MarkovChain) AbsorptionProbabilities() ([][]float64, error) {
	n, err := m.FundamentalMatrix()
	if err != nil {
		return nil, probabilityErrorf("AbsorptionProbabilities: %v", err)
	}
	trans := m.TransientStates()
	abs := m.AbsorbingStates()
	t := len(trans)
	// Build R (transient-by-absorbing).
	r := make([][]float64, t)
	for a := 0; a < t; a++ {
		r[a] = make([]float64, len(abs))
		for k := range abs {
			r[a][k] = m.P[trans[a]][abs[k]]
		}
	}
	return probabilityMatMul(n, r), nil
}
