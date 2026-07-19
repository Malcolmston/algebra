package markov

import (
	"math"
	"math/rand"
)

// HMM is a discrete hidden Markov model. It has N hidden states and M
// observable symbols, a state-transition matrix A (N×N, row-stochastic), an
// emission matrix B (N×M, row-stochastic) where B[i][k] is the probability of
// emitting symbol k in state i, and an initial state distribution Pi (length N).
type HMM struct {
	A  [][]float64
	B  [][]float64
	Pi []float64
	n  int
	m  int
}

// NewHMM constructs a hidden Markov model from the transition matrix a, the
// emission matrix b, and the initial distribution pi. All inputs are copied and
// validated to be (row-)stochastic within DefaultTol. It returns an error if
// the shapes are inconsistent or the rows are not stochastic.
func NewHMM(a, b [][]float64, pi []float64) (*HMM, error) {
	if !IsSquare(a) {
		return nil, ErrNotSquare
	}
	n := len(a)
	if len(b) != n || len(pi) != n {
		return nil, ErrDimMismatch
	}
	if len(b) == 0 || len(b[0]) == 0 {
		return nil, ErrEmpty
	}
	m := len(b[0])
	for i := range b {
		if len(b[i]) != m {
			return nil, ErrDimMismatch
		}
	}
	if !IsStochastic(a, 1e-9) {
		return nil, ErrNotStochastic
	}
	if !IsProbabilityVector(pi, 1e-9) {
		return nil, ErrNotStochastic
	}
	for i := range b {
		var s float64
		for _, x := range b[i] {
			if x < -1e-9 {
				return nil, ErrNotStochastic
			}
			s += x
		}
		if math.Abs(s-1) > 1e-9 {
			return nil, ErrNotStochastic
		}
	}
	return &HMM{
		A:  CopyMatrix(a),
		B:  CopyMatrix(b),
		Pi: CopyVector(pi),
		n:  n,
		m:  m,
	}, nil
}

// NumStates returns the number of hidden states N.
func (h *HMM) NumStates() int { return h.n }

// NumSymbols returns the number of observable symbols M.
func (h *HMM) NumSymbols() int { return h.m }

// validObs reports whether obs contains only symbols in [0, M).
func (h *HMM) validObs(obs []int) bool {
	for _, o := range obs {
		if o < 0 || o >= h.m {
			return false
		}
	}
	return true
}

// Forward runs the unscaled forward algorithm and returns the alpha trellis
// (T×N), where alpha[t][i] = P(o_0..o_t, x_t=i). For long sequences the values
// underflow; prefer ForwardScaled. It also returns the total likelihood
// P(obs) = Σ_i alpha[T-1][i].
func (h *HMM) Forward(obs []int) (alpha [][]float64, likelihood float64, err error) {
	if len(obs) == 0 {
		return nil, 0, ErrEmpty
	}
	if !h.validObs(obs) {
		return nil, 0, ErrDimMismatch
	}
	T := len(obs)
	alpha = make([][]float64, T)
	alpha[0] = make([]float64, h.n)
	for i := 0; i < h.n; i++ {
		alpha[0][i] = h.Pi[i] * h.B[i][obs[0]]
	}
	for t := 1; t < T; t++ {
		alpha[t] = make([]float64, h.n)
		for j := 0; j < h.n; j++ {
			var s float64
			for i := 0; i < h.n; i++ {
				s += alpha[t-1][i] * h.A[i][j]
			}
			alpha[t][j] = s * h.B[j][obs[t]]
		}
	}
	for i := 0; i < h.n; i++ {
		likelihood += alpha[T-1][i]
	}
	return alpha, likelihood, nil
}

// ForwardScaled runs the numerically stable scaled forward algorithm. It
// returns the scaled alpha trellis (each row sums to 1), the per-step scaling
// factors c (length T), and the log-likelihood log P(obs) = -Σ_t log c[t].
func (h *HMM) ForwardScaled(obs []int) (alpha [][]float64, scale []float64, logLik float64, err error) {
	if len(obs) == 0 {
		return nil, nil, 0, ErrEmpty
	}
	if !h.validObs(obs) {
		return nil, nil, 0, ErrDimMismatch
	}
	T := len(obs)
	alpha = make([][]float64, T)
	scale = make([]float64, T)
	alpha[0] = make([]float64, h.n)
	var s0 float64
	for i := 0; i < h.n; i++ {
		alpha[0][i] = h.Pi[i] * h.B[i][obs[0]]
		s0 += alpha[0][i]
	}
	if s0 == 0 {
		return nil, nil, math.Inf(-1), nil
	}
	scale[0] = s0
	for i := 0; i < h.n; i++ {
		alpha[0][i] /= s0
	}
	logLik = math.Log(s0)
	for t := 1; t < T; t++ {
		alpha[t] = make([]float64, h.n)
		var st float64
		for j := 0; j < h.n; j++ {
			var s float64
			for i := 0; i < h.n; i++ {
				s += alpha[t-1][i] * h.A[i][j]
			}
			alpha[t][j] = s * h.B[j][obs[t]]
			st += alpha[t][j]
		}
		if st == 0 {
			return nil, nil, math.Inf(-1), nil
		}
		scale[t] = st
		for j := 0; j < h.n; j++ {
			alpha[t][j] /= st
		}
		logLik += math.Log(st)
	}
	return alpha, scale, logLik, nil
}

// BackwardScaled runs the scaled backward algorithm using the scaling factors
// produced by ForwardScaled. It returns the scaled beta trellis (T×N).
func (h *HMM) BackwardScaled(obs []int, scale []float64) (beta [][]float64, err error) {
	if len(obs) == 0 {
		return nil, ErrEmpty
	}
	if len(scale) != len(obs) {
		return nil, ErrDimMismatch
	}
	if !h.validObs(obs) {
		return nil, ErrDimMismatch
	}
	T := len(obs)
	beta = make([][]float64, T)
	beta[T-1] = make([]float64, h.n)
	for i := 0; i < h.n; i++ {
		beta[T-1][i] = 1 / scale[T-1]
	}
	for t := T - 2; t >= 0; t-- {
		beta[t] = make([]float64, h.n)
		for i := 0; i < h.n; i++ {
			var s float64
			for j := 0; j < h.n; j++ {
				s += h.A[i][j] * h.B[j][obs[t+1]] * beta[t+1][j]
			}
			beta[t][i] = s / scale[t]
		}
	}
	return beta, nil
}

// Backward runs the unscaled backward algorithm and returns the beta trellis
// (T×N), where beta[t][i] = P(o_{t+1}..o_{T-1} | x_t=i). Prefer BackwardScaled
// for long sequences.
func (h *HMM) Backward(obs []int) (beta [][]float64, err error) {
	if len(obs) == 0 {
		return nil, ErrEmpty
	}
	if !h.validObs(obs) {
		return nil, ErrDimMismatch
	}
	T := len(obs)
	beta = make([][]float64, T)
	beta[T-1] = make([]float64, h.n)
	for i := 0; i < h.n; i++ {
		beta[T-1][i] = 1
	}
	for t := T - 2; t >= 0; t-- {
		beta[t] = make([]float64, h.n)
		for i := 0; i < h.n; i++ {
			var s float64
			for j := 0; j < h.n; j++ {
				s += h.A[i][j] * h.B[j][obs[t+1]] * beta[t+1][j]
			}
			beta[t][i] = s
		}
	}
	return beta, nil
}

// Likelihood returns P(obs), the total probability of the observation sequence
// under the model.
func (h *HMM) Likelihood(obs []int) (float64, error) {
	_, _, logLik, err := h.ForwardScaled(obs)
	if err != nil {
		return 0, err
	}
	return math.Exp(logLik), nil
}

// LogLikelihood returns log P(obs), computed stably via the scaled forward
// algorithm.
func (h *HMM) LogLikelihood(obs []int) (float64, error) {
	_, _, logLik, err := h.ForwardScaled(obs)
	if err != nil {
		return 0, err
	}
	return logLik, nil
}

// Viterbi returns the most likely hidden-state sequence for obs and the log
// probability of that path, using the log-domain Viterbi algorithm.
func (h *HMM) Viterbi(obs []int) (path []int, logProb float64, err error) {
	if len(obs) == 0 {
		return nil, 0, ErrEmpty
	}
	if !h.validObs(obs) {
		return nil, 0, ErrDimMismatch
	}
	T := len(obs)
	logInf := math.Inf(-1)
	delta := make([][]float64, T)
	psi := make([][]int, T)
	delta[0] = make([]float64, h.n)
	psi[0] = make([]int, h.n)
	for i := 0; i < h.n; i++ {
		delta[0][i] = logSafe(h.Pi[i]) + logSafe(h.B[i][obs[0]])
	}
	for t := 1; t < T; t++ {
		delta[t] = make([]float64, h.n)
		psi[t] = make([]int, h.n)
		for j := 0; j < h.n; j++ {
			best := logInf
			arg := 0
			for i := 0; i < h.n; i++ {
				v := delta[t-1][i] + logSafe(h.A[i][j])
				if v > best {
					best = v
					arg = i
				}
			}
			delta[t][j] = best + logSafe(h.B[j][obs[t]])
			psi[t][j] = arg
		}
	}
	// Termination.
	logProb = logInf
	last := 0
	for i := 0; i < h.n; i++ {
		if delta[T-1][i] > logProb {
			logProb = delta[T-1][i]
			last = i
		}
	}
	path = make([]int, T)
	path[T-1] = last
	for t := T - 2; t >= 0; t-- {
		path[t] = psi[t+1][path[t+1]]
	}
	return path, logProb, nil
}

// PosteriorMarginals returns the gamma matrix (T×N) of smoothed state
// posteriors, gamma[t][i] = P(x_t=i | obs), using the scaled forward-backward
// algorithm.
func (h *HMM) PosteriorMarginals(obs []int) ([][]float64, error) {
	alpha, scale, _, err := h.ForwardScaled(obs)
	if err != nil {
		return nil, err
	}
	beta, err := h.BackwardScaled(obs, scale)
	if err != nil {
		return nil, err
	}
	T := len(obs)
	gamma := make([][]float64, T)
	for t := 0; t < T; t++ {
		gamma[t] = make([]float64, h.n)
		var s float64
		for i := 0; i < h.n; i++ {
			gamma[t][i] = alpha[t][i] * beta[t][i]
			s += gamma[t][i]
		}
		if s > 0 {
			for i := 0; i < h.n; i++ {
				gamma[t][i] /= s
			}
		}
	}
	return gamma, nil
}

// MostLikelyStates returns the sequence of individually most probable states,
// argmax_i P(x_t=i | obs) for each t. Unlike Viterbi this maximizes each
// marginal independently and may not be a valid path.
func (h *HMM) MostLikelyStates(obs []int) ([]int, error) {
	gamma, err := h.PosteriorMarginals(obs)
	if err != nil {
		return nil, err
	}
	out := make([]int, len(obs))
	for t := range gamma {
		best := -1.0
		arg := 0
		for i := 0; i < h.n; i++ {
			if gamma[t][i] > best {
				best = gamma[t][i]
				arg = i
			}
		}
		out[t] = arg
	}
	return out, nil
}

// StationaryDistribution returns the stationary distribution of the HMM's
// hidden-state transition matrix A.
func (h *HMM) StationaryDistribution() ([]float64, error) {
	chain, err := NewMarkovChainUnchecked(h.A)
	if err != nil {
		return nil, err
	}
	return chain.StationaryDistribution()
}

// Generate samples a state path and observation sequence of the given length
// from the model using rng. It returns the hidden states and the emitted
// symbols. It returns nil slices if length <= 0 or rng is nil.
func (h *HMM) Generate(length int, rng *rand.Rand) (states, obs []int) {
	if length <= 0 || rng == nil {
		return nil, nil
	}
	states = make([]int, length)
	obs = make([]int, length)
	cur := SampleCategorical(h.Pi, rng)
	for t := 0; t < length; t++ {
		states[t] = cur
		obs[t] = SampleCategorical(h.B[cur], rng)
		cur = SampleCategorical(h.A[cur], rng)
	}
	return states, obs
}

// BaumWelch trains the model on a single observation sequence using the
// scaled Baum-Welch (EM) algorithm. It returns a new trained HMM, the
// per-iteration log-likelihood history, and any error. Iteration stops when the
// log-likelihood improves by less than tol or after maxIter iterations. The
// receiver is not modified.
func (h *HMM) BaumWelch(obs []int, maxIter int, tol float64) (*HMM, []float64, error) {
	return h.BaumWelchMultiple([][]int{obs}, maxIter, tol)
}

// BaumWelchMultiple trains the model on several observation sequences at once,
// accumulating the expected sufficient statistics across all sequences before
// each re-estimation. It returns the trained model and the log-likelihood
// (summed over sequences) history. The receiver is not modified.
func (h *HMM) BaumWelchMultiple(seqs [][]int, maxIter int, tol float64) (*HMM, []float64, error) {
	if len(seqs) == 0 {
		return nil, nil, ErrEmpty
	}
	for _, o := range seqs {
		if len(o) == 0 {
			return nil, nil, ErrEmpty
		}
		if !h.validObs(o) {
			return nil, nil, ErrDimMismatch
		}
	}
	if maxIter <= 0 {
		maxIter = 100
	}
	N, M := h.n, h.m
	// Work on copies.
	A := CopyMatrix(h.A)
	B := CopyMatrix(h.B)
	Pi := CopyVector(h.Pi)
	var history []float64
	prevLL := math.Inf(-1)

	for iter := 0; iter < maxIter; iter++ {
		cur := &HMM{A: A, B: B, Pi: Pi, n: N, m: M}
		// Accumulators.
		piNum := make([]float64, N)
		aNum := make([][]float64, N)
		aDen := make([]float64, N)
		bNum := make([][]float64, N)
		bDen := make([]float64, N)
		for i := 0; i < N; i++ {
			aNum[i] = make([]float64, N)
			bNum[i] = make([]float64, M)
		}
		var totalLL float64

		for _, obs := range seqs {
			alpha, scale, ll, err := cur.ForwardScaled(obs)
			if err != nil {
				return nil, nil, err
			}
			if math.IsInf(ll, -1) {
				continue
			}
			beta, err := cur.BackwardScaled(obs, scale)
			if err != nil {
				return nil, nil, err
			}
			totalLL += ll
			T := len(obs)
			// gamma[t][i]
			gamma := make([][]float64, T)
			for t := 0; t < T; t++ {
				gamma[t] = make([]float64, N)
				var s float64
				for i := 0; i < N; i++ {
					gamma[t][i] = alpha[t][i] * beta[t][i]
					s += gamma[t][i]
				}
				if s > 0 {
					for i := 0; i < N; i++ {
						gamma[t][i] /= s
					}
				}
			}
			// Initial distribution.
			for i := 0; i < N; i++ {
				piNum[i] += gamma[0][i]
			}
			// xi transitions.
			for t := 0; t < T-1; t++ {
				var denom float64
				tmp := make([][]float64, N)
				for i := 0; i < N; i++ {
					tmp[i] = make([]float64, N)
					for j := 0; j < N; j++ {
						tmp[i][j] = alpha[t][i] * A[i][j] * B[j][obs[t+1]] * beta[t+1][j]
						denom += tmp[i][j]
					}
				}
				if denom == 0 {
					continue
				}
				for i := 0; i < N; i++ {
					for j := 0; j < N; j++ {
						xi := tmp[i][j] / denom
						aNum[i][j] += xi
					}
					aDen[i] += gamma[t][i]
				}
			}
			// Emission accumulation over all t.
			for t := 0; t < T; t++ {
				for i := 0; i < N; i++ {
					bNum[i][obs[t]] += gamma[t][i]
					bDen[i] += gamma[t][i]
				}
			}
		}

		history = append(history, totalLL)
		// Re-estimate.
		nSeq := float64(len(seqs))
		newPi := make([]float64, N)
		for i := 0; i < N; i++ {
			newPi[i] = piNum[i] / nSeq
		}
		newA := make([][]float64, N)
		newB := make([][]float64, N)
		for i := 0; i < N; i++ {
			newA[i] = make([]float64, N)
			if aDen[i] > 0 {
				for j := 0; j < N; j++ {
					newA[i][j] = aNum[i][j] / aDen[i]
				}
			} else {
				copy(newA[i], A[i])
			}
			newB[i] = make([]float64, M)
			if bDen[i] > 0 {
				for k := 0; k < M; k++ {
					newB[i][k] = bNum[i][k] / bDen[i]
				}
			} else {
				copy(newB[i], B[i])
			}
		}
		A, B, Pi = NormalizeRows(newA), NormalizeRows(newB), Normalize(newPi)

		if iter > 0 && totalLL-prevLL < tol {
			prevLL = totalLL
			break
		}
		prevLL = totalLL
	}

	trained := &HMM{A: A, B: B, Pi: Pi, n: N, m: M}
	return trained, history, nil
}

// logSafe returns log(x) for x>0 and -Inf for x<=0, avoiding NaN.
func logSafe(x float64) float64 {
	if x <= 0 {
		return math.Inf(-1)
	}
	return math.Log(x)
}
