package markov

import (
	"math"
	"math/rand"
)

// MHResult holds the output of a Metropolis-Hastings run: the chain of samples
// and the fraction of proposed moves that were accepted.
type MHResult struct {
	// Samples is the recorded chain (one entry per retained iteration).
	Samples [][]float64
	// AcceptanceRate is the fraction of proposals accepted over the whole run.
	AcceptanceRate float64
}

// RandomWalkMetropolis runs a symmetric random-walk Metropolis sampler in one
// dimension targeting the distribution whose log-density is logTarget (specified
// up to an additive constant). At each step a candidate x' = x + stepSize·Z is
// proposed with Z standard normal, accepted with probability min(1,
// exp(logTarget(x') - logTarget(x))). It returns n samples starting from x0
// (x0 itself is not included; the first returned sample is after one step).
func RandomWalkMetropolis(logTarget func(float64) float64, x0, stepSize float64, n int, rng *rand.Rand) ([]float64, float64) {
	if n <= 0 || rng == nil || logTarget == nil {
		return nil, 0
	}
	samples := make([]float64, n)
	x := x0
	lp := logTarget(x)
	accepted := 0
	for i := 0; i < n; i++ {
		cand := x + stepSize*rng.NormFloat64()
		lpc := logTarget(cand)
		if math.Log(rng.Float64()) < lpc-lp {
			x = cand
			lp = lpc
			accepted++
		}
		samples[i] = x
	}
	return samples, float64(accepted) / float64(n)
}

// MetropolisHastings1D runs a general (possibly asymmetric) Metropolis-Hastings
// sampler in one dimension. logTarget is the target log-density (up to a
// constant). propose draws a candidate given the current state and rng.
// logProposal returns the proposal log-density log q(to | from). The
// Metropolis-Hastings acceptance ratio includes the proposal correction. It
// returns the chain and the acceptance rate.
func MetropolisHastings1D(
	logTarget func(float64) float64,
	propose func(x float64, rng *rand.Rand) float64,
	logProposal func(to, from float64) float64,
	x0 float64, n int, rng *rand.Rand,
) ([]float64, float64) {
	if n <= 0 || rng == nil || logTarget == nil || propose == nil || logProposal == nil {
		return nil, 0
	}
	samples := make([]float64, n)
	x := x0
	lp := logTarget(x)
	accepted := 0
	for i := 0; i < n; i++ {
		cand := propose(x, rng)
		lpc := logTarget(cand)
		logAlpha := (lpc + logProposal(x, cand)) - (lp + logProposal(cand, x))
		if math.Log(rng.Float64()) < logAlpha {
			x = cand
			lp = lpc
			accepted++
		}
		samples[i] = x
	}
	return samples, float64(accepted) / float64(n)
}

// RandomWalkMetropolisND runs a random-walk Metropolis sampler in d dimensions
// with an isotropic Gaussian proposal of the given step size. logTarget maps a
// point to its target log-density up to a constant. It returns an MHResult
// whose Samples has n rows of length len(x0).
func RandomWalkMetropolisND(logTarget func([]float64) float64, x0 []float64, stepSize float64, n int, rng *rand.Rand) *MHResult {
	if n <= 0 || rng == nil || logTarget == nil || len(x0) == 0 {
		return nil
	}
	d := len(x0)
	samples := make([][]float64, n)
	x := CopyVector(x0)
	lp := logTarget(x)
	accepted := 0
	for i := 0; i < n; i++ {
		cand := make([]float64, d)
		for k := 0; k < d; k++ {
			cand[k] = x[k] + stepSize*rng.NormFloat64()
		}
		lpc := logTarget(cand)
		if math.Log(rng.Float64()) < lpc-lp {
			x = cand
			lp = lpc
			accepted++
		}
		samples[i] = CopyVector(x)
	}
	return &MHResult{Samples: samples, AcceptanceRate: float64(accepted) / float64(n)}
}

// GibbsSampler runs a systematic-scan Gibbs sampler. conditionals[k] draws a new
// value for coordinate k given the current full state x (which the sampler
// updates in place between coordinate draws) and rng. It returns n recorded
// states, each of length len(x0). One recorded state corresponds to a full
// sweep over all coordinates.
func GibbsSampler(conditionals []func(x []float64, rng *rand.Rand) float64, x0 []float64, n int, rng *rand.Rand) [][]float64 {
	d := len(x0)
	if n <= 0 || rng == nil || len(conditionals) != d || d == 0 {
		return nil
	}
	x := CopyVector(x0)
	out := make([][]float64, n)
	for i := 0; i < n; i++ {
		for k := 0; k < d; k++ {
			x[k] = conditionals[k](x, rng)
		}
		out[i] = CopyVector(x)
	}
	return out
}

// Burnin returns the chain with the first burn samples discarded.
func Burnin(chain [][]float64, burn int) [][]float64 {
	if burn <= 0 || burn >= len(chain) {
		if burn >= len(chain) {
			return nil
		}
		return chain
	}
	return chain[burn:]
}

// Thin returns every step-th sample of the chain (step>=1).
func Thin(chain [][]float64, step int) [][]float64 {
	if step <= 1 {
		return chain
	}
	var out [][]float64
	for i := 0; i < len(chain); i += step {
		out = append(out, chain[i])
	}
	return out
}

// SampleMean returns the mean of a scalar chain.
func SampleMean(chain []float64) float64 {
	if len(chain) == 0 {
		return math.NaN()
	}
	var s float64
	for _, x := range chain {
		s += x
	}
	return s / float64(len(chain))
}

// SampleVariance returns the unbiased (n-1) sample variance of a scalar chain.
func SampleVariance(chain []float64) float64 {
	n := len(chain)
	if n < 2 {
		return math.NaN()
	}
	m := SampleMean(chain)
	var s float64
	for _, x := range chain {
		d := x - m
		s += d * d
	}
	return s / float64(n-1)
}

// SampleStdDev returns the square root of the unbiased sample variance.
func SampleStdDev(chain []float64) float64 {
	return math.Sqrt(SampleVariance(chain))
}

// Autocorrelation returns the sample autocorrelation of a scalar chain at the
// given lag (using the biased 1/N normalization). It returns NaN for an empty
// chain and 1 for lag 0.
func Autocorrelation(chain []float64, lag int) float64 {
	n := len(chain)
	if n == 0 || lag < 0 || lag >= n {
		return math.NaN()
	}
	m := SampleMean(chain)
	var denom float64
	for _, x := range chain {
		d := x - m
		denom += d * d
	}
	if denom == 0 {
		return math.NaN()
	}
	var num float64
	for i := 0; i < n-lag; i++ {
		num += (chain[i] - m) * (chain[i+lag] - m)
	}
	return num / denom
}

// AutocorrelationFunction returns the autocorrelation of the chain for lags 0
// through maxLag inclusive.
func AutocorrelationFunction(chain []float64, maxLag int) []float64 {
	if maxLag < 0 {
		return nil
	}
	out := make([]float64, maxLag+1)
	for k := 0; k <= maxLag; k++ {
		out[k] = Autocorrelation(chain, k)
	}
	return out
}

// IntegratedAutocorrelationTime estimates the integrated autocorrelation time
// tau = 1 + 2 Σ_{k>=1} ρ(k), truncating the sum at the first lag where the
// autocorrelation becomes non-positive (the initial-positive-sequence rule) or
// at maxLag. It returns at least 1.
func IntegratedAutocorrelationTime(chain []float64, maxLag int) float64 {
	if len(chain) < 2 {
		return math.NaN()
	}
	if maxLag <= 0 || maxLag >= len(chain) {
		maxLag = len(chain) - 1
	}
	tau := 1.0
	for k := 1; k <= maxLag; k++ {
		r := Autocorrelation(chain, k)
		if r <= 0 || math.IsNaN(r) {
			break
		}
		tau += 2 * r
	}
	return tau
}

// EffectiveSampleSize estimates the effective sample size of a scalar chain as
// N / tau, where tau is the integrated autocorrelation time.
func EffectiveSampleSize(chain []float64) float64 {
	n := len(chain)
	if n < 2 {
		return math.NaN()
	}
	tau := IntegratedAutocorrelationTime(chain, 0)
	if tau <= 0 {
		return float64(n)
	}
	return float64(n) / tau
}

// BatchMeansVariance estimates the variance of the sample mean of a chain by the
// non-overlapping batch-means method with the given number of batches. It
// returns the estimated variance of the mean (not of the chain).
func BatchMeansVariance(chain []float64, numBatches int) float64 {
	n := len(chain)
	if numBatches < 2 || n < numBatches {
		return math.NaN()
	}
	batchSize := n / numBatches
	means := make([]float64, numBatches)
	for b := 0; b < numBatches; b++ {
		var s float64
		for i := 0; i < batchSize; i++ {
			s += chain[b*batchSize+i]
		}
		means[b] = s / float64(batchSize)
	}
	grand := SampleMean(means)
	var s float64
	for _, m := range means {
		d := m - grand
		s += d * d
	}
	// Variance of the overall mean estimate.
	return s / float64(numBatches-1) / float64(numBatches)
}

// GelmanRubin returns the Gelman-Rubin potential scale reduction factor R-hat
// for several scalar chains (all of the same length m). Values close to 1
// indicate convergence. It returns NaN if fewer than two chains are supplied or
// the lengths differ.
func GelmanRubin(chains [][]float64) float64 {
	M := len(chains)
	if M < 2 {
		return math.NaN()
	}
	m := len(chains[0])
	if m < 2 {
		return math.NaN()
	}
	for _, ch := range chains {
		if len(ch) != m {
			return math.NaN()
		}
	}
	chainMeans := make([]float64, M)
	chainVars := make([]float64, M)
	for j, ch := range chains {
		chainMeans[j] = SampleMean(ch)
		chainVars[j] = SampleVariance(ch)
	}
	grand := SampleMean(chainMeans)
	// Between-chain variance B.
	var b float64
	for j := 0; j < M; j++ {
		d := chainMeans[j] - grand
		b += d * d
	}
	b *= float64(m) / float64(M-1)
	// Within-chain variance W.
	var w float64
	for j := 0; j < M; j++ {
		w += chainVars[j]
	}
	w /= float64(M)
	if w == 0 {
		return math.NaN()
	}
	varHat := (float64(m-1)/float64(m))*w + b/float64(m)
	return math.Sqrt(varHat / w)
}

// Column extracts column k from a multidimensional chain (a slice of state
// vectors) as a scalar series, for use with the scalar diagnostics.
func Column(chain [][]float64, k int) []float64 {
	out := make([]float64, 0, len(chain))
	for _, row := range chain {
		if k < len(row) {
			out = append(out, row[k])
		}
	}
	return out
}

// AcceptanceRate returns the fraction of consecutive samples in a scalar chain
// that changed value, an empirical proxy for the Metropolis acceptance rate
// when the true rate was not recorded.
func AcceptanceRate(chain []float64) float64 {
	if len(chain) < 2 {
		return math.NaN()
	}
	changes := 0
	for i := 1; i < len(chain); i++ {
		if chain[i] != chain[i-1] {
			changes++
		}
	}
	return float64(changes) / float64(len(chain)-1)
}
