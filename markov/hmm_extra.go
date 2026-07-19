package markov

// TransitionProb returns the hidden-state transition probability A[i][j].
func (h *HMM) TransitionProb(i, j int) float64 {
	if i < 0 || i >= h.n || j < 0 || j >= h.n {
		return -1
	}
	return h.A[i][j]
}

// EmissionProb returns the probability B[state][symbol] of emitting symbol in
// the given hidden state.
func (h *HMM) EmissionProb(state, symbol int) float64 {
	if state < 0 || state >= h.n || symbol < 0 || symbol >= h.m {
		return -1
	}
	return h.B[state][symbol]
}

// TransitionMatrix returns a copy of the hidden-state transition matrix A.
func (h *HMM) TransitionMatrix() [][]float64 { return CopyMatrix(h.A) }

// EmissionMatrix returns a copy of the emission matrix B.
func (h *HMM) EmissionMatrix() [][]float64 { return CopyMatrix(h.B) }

// InitialDistribution returns a copy of the initial state distribution Pi.
func (h *HMM) InitialDistribution() []float64 { return CopyVector(h.Pi) }

// Decode is an alias for Viterbi: it returns the most likely hidden-state path
// and its log probability.
func (h *HMM) Decode(obs []int) ([]int, float64, error) {
	return h.Viterbi(obs)
}

// Filter returns the filtering distribution P(x_{T-1} | o_0..o_{T-1}), i.e. the
// posterior over the final hidden state given the whole observed prefix. It is
// the normalized final row of the scaled forward pass.
func (h *HMM) Filter(obs []int) ([]float64, error) {
	alpha, _, _, err := h.ForwardScaled(obs)
	if err != nil {
		return nil, err
	}
	return CopyVector(alpha[len(alpha)-1]), nil
}

// PredictNextStateDistribution returns the distribution over the hidden state at
// time T (one step past the last observation) given obs.
func (h *HMM) PredictNextStateDistribution(obs []int) ([]float64, error) {
	filt, err := h.Filter(obs)
	if err != nil {
		return nil, err
	}
	return VecMat(filt, h.A), nil
}

// PredictNextObservationDistribution returns the predictive distribution over
// the next observed symbol given obs.
func (h *HMM) PredictNextObservationDistribution(obs []int) ([]float64, error) {
	stateDist, err := h.PredictNextStateDistribution(obs)
	if err != nil {
		return nil, err
	}
	out := make([]float64, h.m)
	for k := 0; k < h.m; k++ {
		var s float64
		for i := 0; i < h.n; i++ {
			s += stateDist[i] * h.B[i][k]
		}
		out[k] = s
	}
	return out, nil
}

// ExpectedTransitionCounts returns the matrix of expected numbers of i→j
// transitions over the sequence obs, Σ_t P(x_t=i, x_{t+1}=j | obs), from the
// posterior xi statistics. It is the E-step transition accumulator of
// Baum-Welch for a single sequence.
func (h *HMM) ExpectedTransitionCounts(obs []int) ([][]float64, error) {
	alpha, scale, ll, err := h.ForwardScaled(obs)
	if err != nil {
		return nil, err
	}
	out := make([][]float64, h.n)
	for i := range out {
		out[i] = make([]float64, h.n)
	}
	if len(obs) < 2 || ll == 0 {
		return out, nil
	}
	beta, err := h.BackwardScaled(obs, scale)
	if err != nil {
		return nil, err
	}
	T := len(obs)
	for t := 0; t < T-1; t++ {
		var denom float64
		tmp := make([][]float64, h.n)
		for i := 0; i < h.n; i++ {
			tmp[i] = make([]float64, h.n)
			for j := 0; j < h.n; j++ {
				tmp[i][j] = alpha[t][i] * h.A[i][j] * h.B[j][obs[t+1]] * beta[t+1][j]
				denom += tmp[i][j]
			}
		}
		if denom == 0 {
			continue
		}
		for i := 0; i < h.n; i++ {
			for j := 0; j < h.n; j++ {
				out[i][j] += tmp[i][j] / denom
			}
		}
	}
	return out, nil
}
