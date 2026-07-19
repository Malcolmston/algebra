package markov

import "math"

// ProbabilityOfPath returns the probability of the state sequence conditional on
// starting in states[0], i.e. the product of the one-step transition
// probabilities P[states[t]][states[t+1]]. It returns NaN if any index is out of
// range and 1 for a path of length 0 or 1.
func (c *MarkovChain) ProbabilityOfPath(states []int) float64 {
	prob := 1.0
	for t := 0; t < len(states); t++ {
		if states[t] < 0 || states[t] >= c.n {
			return math.NaN()
		}
		if t+1 < len(states) {
			prob *= c.p[states[t]][states[t+1]]
		}
	}
	return prob
}

// LogProbabilityOfPath returns the natural logarithm of ProbabilityOfPath,
// computed as a sum of logs for numerical stability. A zero-probability
// transition yields -Inf.
func (c *MarkovChain) LogProbabilityOfPath(states []int) float64 {
	var lp float64
	for t := 0; t+1 < len(states); t++ {
		if states[t] < 0 || states[t] >= c.n || states[t+1] < 0 || states[t+1] >= c.n {
			return math.NaN()
		}
		p := c.p[states[t]][states[t+1]]
		if p <= 0 {
			return math.Inf(-1)
		}
		lp += math.Log(p)
	}
	return lp
}

// MeanFirstPassageTime returns the expected number of steps to first reach state
// j starting from state i in an ergodic chain. For i==j it returns the mean
// recurrence time 1/π_j.
func (c *MarkovChain) MeanFirstPassageTime(i, j int) (float64, error) {
	m, err := c.MeanFirstPassageTimes()
	if err != nil {
		return 0, err
	}
	if i < 0 || i >= c.n || j < 0 || j >= c.n {
		return 0, ErrDimMismatch
	}
	return m[i][j], nil
}

// HittingProbability returns the probability that the chain starting in state
// from ever reaches the target set.
func (c *MarkovChain) HittingProbability(from int, target []int) (float64, error) {
	h, err := c.HittingProbabilities(target)
	if err != nil {
		return 0, err
	}
	if from < 0 || from >= c.n {
		return 0, ErrDimMismatch
	}
	return h[from], nil
}

// ExpectedHittingTime returns the expected number of steps to first reach the
// target set, starting from state from (0 if from is already in the target).
func (c *MarkovChain) ExpectedHittingTime(from int, target []int) (float64, error) {
	h, err := c.ExpectedHittingTimes(target)
	if err != nil {
		return 0, err
	}
	if from < 0 || from >= c.n {
		return 0, ErrDimMismatch
	}
	return h[from], nil
}

// StationaryEntropy returns the Shannon entropy (in nats) of the chain's
// stationary distribution.
func (c *MarkovChain) StationaryEntropy() (float64, error) {
	pi, err := c.StationaryDistribution()
	if err != nil {
		return 0, err
	}
	return ShannonEntropy(pi), nil
}

// TVDistanceToStationary returns the total-variation distance between the
// distribution after k steps (starting from init) and the stationary
// distribution. It requires the chain to have a well-defined stationary
// distribution.
func (c *MarkovChain) TVDistanceToStationary(init []float64, k int) (float64, error) {
	if len(init) != c.n {
		return 0, ErrDimMismatch
	}
	pi, err := c.StationaryDistribution()
	if err != nil {
		return 0, err
	}
	dist := c.StepDistribution(init, k)
	return TotalVariationDistance(dist, pi), nil
}

// MixingTime returns the smallest number of steps k such that, from every
// starting state, the total-variation distance between the k-step distribution
// and the stationary distribution is at most epsilon. It searches up to
// maxSteps and returns -1 if the bound is not met within that horizon. The
// chain must be ergodic.
func (c *MarkovChain) MixingTime(epsilon float64, maxSteps int) (int, error) {
	if !c.IsErgodic() {
		return 0, ErrNotErgodic
	}
	pi, err := c.StationaryDistribution()
	if err != nil {
		return 0, err
	}
	if maxSteps <= 0 {
		maxSteps = 100000
	}
	// Iterate the full power, row by row.
	pk := Identity(c.n)
	for k := 0; k <= maxSteps; k++ {
		worst := 0.0
		for i := 0; i < c.n; i++ {
			if d := TotalVariationDistance(pk[i], pi); d > worst {
				worst = d
			}
		}
		if worst <= epsilon {
			return k, nil
		}
		pk = MatMul(pk, c.p)
	}
	return -1, nil
}

// RecurrentClasses returns the communicating classes that are recurrent (closed
// under the dynamics), each as a sorted slice.
func (c *MarkovChain) RecurrentClasses() [][]int {
	var out [][]int
	for _, comp := range c.CommunicatingClasses() {
		if c.IsClosedClass(comp) {
			out = append(out, comp)
		}
	}
	return out
}

// TransientClasses returns the communicating classes that are transient (not
// closed), each as a sorted slice.
func (c *MarkovChain) TransientClasses() [][]int {
	var out [][]int
	for _, comp := range c.CommunicatingClasses() {
		if !c.IsClosedClass(comp) {
			out = append(out, comp)
		}
	}
	return out
}

// NumClasses returns the number of communicating classes.
func (c *MarkovChain) NumClasses() int {
	return len(c.CommunicatingClasses())
}

// ClassifyStates returns a map from every state to its classification label,
// one of "absorbing", "recurrent", or "transient".
func (c *MarkovChain) ClassifyStates() map[int]string {
	out := make(map[int]string, c.n)
	rec := make(map[int]bool)
	for _, s := range c.RecurrentStates() {
		rec[s] = true
	}
	for s := 0; s < c.n; s++ {
		switch {
		case c.IsAbsorbingState(s):
			out[s] = "absorbing"
		case rec[s]:
			out[s] = "recurrent"
		default:
			out[s] = "transient"
		}
	}
	return out
}

// ExpectedReward returns the long-run average reward per step under the
// stationary distribution, given a per-state reward vector: Σ_i π_i r_i.
func (c *MarkovChain) ExpectedReward(reward []float64) (float64, error) {
	if len(reward) != c.n {
		return 0, ErrDimMismatch
	}
	pi, err := c.StationaryDistribution()
	if err != nil {
		return 0, err
	}
	return Dot(pi, reward), nil
}
