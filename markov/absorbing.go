package markov

import (
	"math"
	"sort"
)

// CanonicalForm reorders the states of an absorbing chain into transient states
// followed by absorbing states and returns the reordered transition matrix
// along with the transient and absorbing index lists (in original numbering).
// The reordered matrix has the block structure [[Q, R], [0, I]].
func (c *MarkovChain) CanonicalForm() (reordered [][]float64, transient, absorbing []int) {
	absorbing = c.AbsorbingStates()
	absSet := make(map[int]bool)
	for _, a := range absorbing {
		absSet[a] = true
	}
	for i := 0; i < c.n; i++ {
		if !absSet[i] {
			transient = append(transient, i)
		}
	}
	order := append(append([]int{}, transient...), absorbing...)
	reordered = make([][]float64, c.n)
	for ni, oi := range order {
		reordered[ni] = make([]float64, c.n)
		for nj, oj := range order {
			reordered[ni][nj] = c.p[oi][oj]
		}
	}
	return reordered, transient, absorbing
}

// SubmatrixQ returns the transient-to-transient block Q of the canonical form
// (indexed by the transient states in ascending order).
func (c *MarkovChain) SubmatrixQ() ([][]float64, []int) {
	transient := c.TransientStates()
	t := len(transient)
	q := make([][]float64, t)
	for a, i := range transient {
		q[a] = make([]float64, t)
		for b, j := range transient {
			q[a][b] = c.p[i][j]
		}
	}
	return q, transient
}

// SubmatrixR returns the transient-to-absorbing block R of the canonical form,
// together with the transient and absorbing index lists.
func (c *MarkovChain) SubmatrixR() (r [][]float64, transient, absorbing []int) {
	transient = c.TransientStates()
	absorbing = c.AbsorbingStates()
	r = make([][]float64, len(transient))
	for a, i := range transient {
		r[a] = make([]float64, len(absorbing))
		for b, j := range absorbing {
			r[a][b] = c.p[i][j]
		}
	}
	return r, transient, absorbing
}

// FundamentalMatrix returns the fundamental matrix N = (I - Q)^{-1} of an
// absorbing chain, indexed by the transient states in ascending order. Entry
// N[i][j] is the expected number of visits to transient state j starting from
// transient state i before absorption. It returns ErrNotAbsorbing if the chain
// has no absorbing states, and the transient index list.
func (c *MarkovChain) FundamentalMatrix() ([][]float64, []int, error) {
	if len(c.AbsorbingStates()) == 0 {
		return nil, nil, ErrNotAbsorbing
	}
	q, transient := c.SubmatrixQ()
	t := len(transient)
	if t == 0 {
		return [][]float64{}, transient, nil
	}
	imq := MatSub(Identity(t), q)
	n, err := MatInverse(imq)
	if err != nil {
		return nil, nil, err
	}
	return n, transient, nil
}

// ExpectedStepsToAbsorption returns, for each transient state, the expected
// number of steps until the chain is absorbed. The result is indexed by the
// transient states in ascending order (also returned).
func (c *MarkovChain) ExpectedStepsToAbsorption() ([]float64, []int, error) {
	n, transient, err := c.FundamentalMatrix()
	if err != nil {
		return nil, nil, err
	}
	t := len(transient)
	steps := make([]float64, t)
	for i := 0; i < t; i++ {
		for j := 0; j < t; j++ {
			steps[i] += n[i][j]
		}
	}
	return steps, transient, nil
}

// VarianceStepsToAbsorption returns the variance of the number of steps to
// absorption for each transient state, using the standard formula
// (2N - I)t - t∘t where t is the expected-steps vector and ∘ is elementwise
// product. The result is indexed by the transient states in ascending order.
func (c *MarkovChain) VarianceStepsToAbsorption() ([]float64, []int, error) {
	n, transient, err := c.FundamentalMatrix()
	if err != nil {
		return nil, nil, err
	}
	t := len(transient)
	tv := make([]float64, t) // expected steps
	for i := 0; i < t; i++ {
		for j := 0; j < t; j++ {
			tv[i] += n[i][j]
		}
	}
	// (2N - I) t
	varr := make([]float64, t)
	for i := 0; i < t; i++ {
		var s float64
		for j := 0; j < t; j++ {
			coef := 2 * n[i][j]
			if i == j {
				coef -= 1
			}
			s += coef * tv[j]
		}
		varr[i] = s - tv[i]*tv[i]
	}
	return varr, transient, nil
}

// AbsorptionProbabilities returns the matrix B = N·R of absorption
// probabilities: B[i][k] is the probability that a chain started in transient
// state i is eventually absorbed in absorbing state k. Rows are indexed by the
// transient states and columns by the absorbing states, both in ascending
// order (also returned).
func (c *MarkovChain) AbsorptionProbabilities() (b [][]float64, transient, absorbing []int, err error) {
	n, transient, err := c.FundamentalMatrix()
	if err != nil {
		return nil, nil, nil, err
	}
	r, _, absorbing := c.SubmatrixR()
	b = MatMul(n, r)
	if b == nil && len(transient) == 0 {
		b = [][]float64{}
	}
	return b, transient, absorbing, nil
}

// AbsorptionProbability returns the probability that the chain started in state
// from is eventually absorbed in the absorbing state target. It returns
// (1, nil) if from == target, (0, nil) if from is a different absorbing state,
// and otherwise looks up the value in the absorption-probability matrix.
func (c *MarkovChain) AbsorptionProbability(from, target int) (float64, error) {
	if from < 0 || from >= c.n || target < 0 || target >= c.n {
		return 0, ErrDimMismatch
	}
	if !c.IsAbsorbingState(target) {
		return 0, ErrNotAbsorbing
	}
	if c.IsAbsorbingState(from) {
		if from == target {
			return 1, nil
		}
		return 0, nil
	}
	b, transient, absorbing, err := c.AbsorptionProbabilities()
	if err != nil {
		return 0, err
	}
	ti := indexOf(transient, from)
	ai := indexOf(absorbing, target)
	if ti < 0 || ai < 0 {
		return 0, ErrDimMismatch
	}
	return b[ti][ai], nil
}

// MeanFirstPassageTimes returns the matrix M of mean first-passage times for an
// ergodic chain: M[i][j] is the expected number of steps to first reach state j
// starting from state i (i≠j). The diagonal M[j][j] holds the mean recurrence
// time 1/π_j. It returns ErrNotErgodic if the chain is not ergodic.
func (c *MarkovChain) MeanFirstPassageTimes() ([][]float64, error) {
	if !c.IsErgodic() {
		return nil, ErrNotErgodic
	}
	n := c.n
	m := make([][]float64, n)
	for i := range m {
		m[i] = make([]float64, n)
	}
	// For each target j, solve (I - Q_j) x = 1 where Q_j is P with row/column j
	// deleted (states other than j), giving expected steps to first hit j.
	for j := 0; j < n; j++ {
		// Build reduced system over states != j.
		idx := make([]int, 0, n-1)
		for s := 0; s < n; s++ {
			if s != j {
				idx = append(idx, s)
			}
		}
		k := len(idx)
		a := make([][]float64, k)
		b := make([]float64, k)
		for r := 0; r < k; r++ {
			a[r] = make([]float64, k)
			for ccol := 0; ccol < k; ccol++ {
				val := -c.p[idx[r]][idx[ccol]]
				if r == ccol {
					val += 1
				}
				a[r][ccol] = val
			}
			b[r] = 1
		}
		x, err := SolveLinear(a, b)
		if err != nil {
			return nil, err
		}
		for r, s := range idx {
			m[s][j] = x[r]
		}
	}
	// Diagonal: mean recurrence time = 1 + Σ_k P_jk m_kj = 1/π_j.
	pi, err := c.StationaryDistribution()
	if err != nil {
		return nil, err
	}
	for j := 0; j < n; j++ {
		if pi[j] > 0 {
			m[j][j] = 1 / pi[j]
		}
	}
	return m, nil
}

// MeanRecurrenceTimes returns the vector of mean recurrence times of an ergodic
// chain, m_jj = 1/π_j.
func (c *MarkovChain) MeanRecurrenceTimes() ([]float64, error) {
	pi, err := c.StationaryDistribution()
	if err != nil {
		return nil, err
	}
	r := make([]float64, c.n)
	for j := 0; j < c.n; j++ {
		if pi[j] > 0 {
			r[j] = 1 / pi[j]
		}
	}
	return r, nil
}

// HittingProbabilities returns, for every state i, the probability that the
// chain starting at i ever reaches the target set. States in the target set
// have probability 1, and states from which the target is unreachable have
// probability 0. The remaining probabilities are obtained by solving the
// linear system h_i = Σ_j P_ij h_j exactly. It returns an error if target is
// empty or contains an out-of-range index.
func (c *MarkovChain) HittingProbabilities(target []int) ([]float64, error) {
	inTarget := make([]bool, c.n)
	for _, t := range target {
		if t < 0 || t >= c.n {
			return nil, ErrDimMismatch
		}
		inTarget[t] = true
	}
	if len(target) == 0 {
		return nil, ErrEmpty
	}
	canReach := c.canReachSet(inTarget)
	// Unknown states are those that can reach the target but are not in it.
	idx := make([]int, 0, c.n)
	for s := 0; s < c.n; s++ {
		if canReach[s] && !inTarget[s] {
			idx = append(idx, s)
		}
	}
	k := len(idx)
	h := make([]float64, c.n)
	for _, t := range target {
		h[t] = 1
	}
	if k == 0 {
		return h, nil
	}
	// (I - Q) x = R·1_target, where Q is transitions among non-target states and
	// R·1 is the one-step probability of hitting the target directly.
	a := make([][]float64, k)
	b := make([]float64, k)
	for r, s := range idx {
		a[r] = make([]float64, k)
		var rhs float64
		for j := 0; j < c.n; j++ {
			if inTarget[j] {
				rhs += c.p[s][j]
			}
		}
		for ccol, s2 := range idx {
			val := -c.p[s][s2]
			if r == ccol {
				val += 1
			}
			a[r][ccol] = val
		}
		b[r] = rhs
	}
	x, err := SolveLinear(a, b)
	if err != nil {
		return nil, err
	}
	for r, s := range idx {
		h[s] = x[r]
	}
	return h, nil
}

// ExpectedHittingTimes returns, for every state i, the expected number of steps
// to first reach the target set (0 for states already in the target). A state
// from which the target is not reached with probability 1 has an infinite
// expected hitting time (+Inf). The finite values are obtained by solving
// (I - Q) x = 1 over the states that reach the target almost surely.
func (c *MarkovChain) ExpectedHittingTimes(target []int) ([]float64, error) {
	inTarget := make([]bool, c.n)
	for _, t := range target {
		if t < 0 || t >= c.n {
			return nil, ErrDimMismatch
		}
		inTarget[t] = true
	}
	if len(target) == 0 {
		return nil, ErrEmpty
	}
	// Finite expected hitting time requires the target to be reached almost
	// surely, i.e. hitting probability 1.
	hp, err := c.HittingProbabilities(target)
	if err != nil {
		return nil, err
	}
	idx := make([]int, 0, c.n)
	out := make([]float64, c.n)
	for s := 0; s < c.n; s++ {
		if inTarget[s] {
			continue
		}
		if hp[s] >= 1-1e-12 {
			idx = append(idx, s)
		} else {
			out[s] = math.Inf(1)
		}
	}
	k := len(idx)
	if k == 0 {
		return out, nil
	}
	a := make([][]float64, k)
	b := make([]float64, k)
	for r, s := range idx {
		a[r] = make([]float64, k)
		for ccol, s2 := range idx {
			val := -c.p[s][s2]
			if r == ccol {
				val += 1
			}
			a[r][ccol] = val
		}
		b[r] = 1
	}
	x, err := SolveLinear(a, b)
	if err != nil {
		return nil, err
	}
	for r, s := range idx {
		out[s] = x[r]
	}
	return out, nil
}

// ExpectedVisits returns the expected number of visits to each transient state
// before absorption, starting from transient state i (given in original state
// numbering). It returns ErrNotAbsorbing if there are no absorbing states, and
// the transient index list identifying the columns.
func (c *MarkovChain) ExpectedVisits(from int) ([]float64, []int, error) {
	n, transient, err := c.FundamentalMatrix()
	if err != nil {
		return nil, nil, err
	}
	ti := indexOf(transient, from)
	if ti < 0 {
		return nil, transient, ErrDimMismatch
	}
	return CopyVector(n[ti]), transient, nil
}

// canReachSet returns the set of states that can reach the target set (states
// marked true in inTarget) via a positive-probability path, including the
// target states themselves. It is computed by a reverse breadth-first search.
func (c *MarkovChain) canReachSet(inTarget []bool) []bool {
	canReach := make([]bool, c.n)
	var queue []int
	for s := 0; s < c.n; s++ {
		if inTarget[s] {
			canReach[s] = true
			queue = append(queue, s)
		}
	}
	for len(queue) > 0 {
		v := queue[0]
		queue = queue[1:]
		for u := 0; u < c.n; u++ {
			if !canReach[u] && c.p[u][v] > 0 {
				canReach[u] = true
				queue = append(queue, u)
			}
		}
	}
	return canReach
}

func indexOf(s []int, v int) int {
	i := sort.SearchInts(s, v)
	if i < len(s) && s[i] == v {
		return i
	}
	return -1
}
