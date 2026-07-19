package markov

import (
	"fmt"
	"math"
	"math/rand"
	"testing"
)

const tol = 1e-9

func approx(a, b, eps float64) bool { return math.Abs(a-b) <= eps }

// ---------------------------------------------------------------------------
// Linear algebra
// ---------------------------------------------------------------------------

func TestMatMulAndPow(t *testing.T) {
	p := [][]float64{{0.9, 0.1}, {0.5, 0.5}}
	p2 := MatMul(p, p)
	want := [][]float64{{0.86, 0.14}, {0.70, 0.30}}
	if !MatEqual(p2, want, 1e-12) {
		t.Errorf("P^2 = %v, want %v", p2, want)
	}
	if !MatEqual(MatPow(p, 2), want, 1e-12) {
		t.Errorf("MatPow mismatch")
	}
	if !MatEqual(MatPow(p, 0), Identity(2), 1e-12) {
		t.Errorf("MatPow(0) should be identity")
	}
}

func TestSolveInverseDeterminant(t *testing.T) {
	a := [][]float64{{2, 1, 1}, {1, 3, 2}, {1, 0, 0}}
	inv, err := MatInverse(a)
	if err != nil {
		t.Fatal(err)
	}
	prod := MatMul(a, inv)
	if !MatEqual(prod, Identity(3), 1e-9) {
		t.Errorf("A·A^-1 != I: %v", prod)
	}
	det, err := Determinant(a)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(det, -1, 1e-9) {
		t.Errorf("det = %v, want -1", det)
	}
	// Solve A x = b.
	b := []float64{4, 5, 6}
	x, err := SolveLinear(a, b)
	if err != nil {
		t.Fatal(err)
	}
	check := MatVec(a, x)
	if !VecEqual(check, b, 1e-9) {
		t.Errorf("A·x = %v, want %v", check, b)
	}
}

func TestSingular(t *testing.T) {
	a := [][]float64{{1, 2}, {2, 4}}
	if _, err := MatInverse(a); err == nil {
		t.Errorf("expected singular error")
	}
	det, _ := Determinant(a)
	if det != 0 {
		t.Errorf("singular det = %v, want 0", det)
	}
}

func TestNorms(t *testing.T) {
	v := []float64{3, -4}
	if !approx(VecNorm2(v), 5, tol) {
		t.Errorf("norm2 = %v", VecNorm2(v))
	}
	if !approx(VecNorm1(v), 7, tol) {
		t.Errorf("norm1 = %v", VecNorm1(v))
	}
	if !approx(VecNormInf(v), 4, tol) {
		t.Errorf("normInf = %v", VecNormInf(v))
	}
}

// ---------------------------------------------------------------------------
// Stochastic predicates & distances
// ---------------------------------------------------------------------------

func TestStochasticPredicates(t *testing.T) {
	p := [][]float64{{0.9, 0.1}, {0.5, 0.5}}
	if !IsStochastic(p, tol) {
		t.Errorf("should be stochastic")
	}
	ds := [][]float64{{0.5, 0.5}, {0.5, 0.5}}
	if !IsDoublyStochastic(ds, tol) {
		t.Errorf("should be doubly stochastic")
	}
	if IsStochastic([][]float64{{0.9, 0.2}, {0.5, 0.5}}, tol) {
		t.Errorf("row sum 1.1 should fail")
	}
}

func TestDistances(t *testing.T) {
	p := []float64{0.5, 0.5}
	q := []float64{0.25, 0.75}
	if !approx(TotalVariationDistance(p, q), 0.25, tol) {
		t.Errorf("TV = %v", TotalVariationDistance(p, q))
	}
	if !approx(L1Distance(p, q), 0.5, tol) {
		t.Errorf("L1 = %v", L1Distance(p, q))
	}
	wantKL := 0.5*math.Log(2) + 0.5*math.Log(2.0/3.0)
	if !approx(KLDivergence(p, q), wantKL, 1e-12) {
		t.Errorf("KL = %v, want %v", KLDivergence(p, q), wantKL)
	}
	// Hellinger is symmetric and in [0,1].
	h := HellingerDistance(p, q)
	if h < 0 || h > 1 {
		t.Errorf("Hellinger out of range: %v", h)
	}
	if !approx(HellingerDistance(p, p), 0, tol) {
		t.Errorf("Hellinger of identical should be 0")
	}
	if !approx(BhattacharyyaCoefficient(p, p), 1, tol) {
		t.Errorf("BC of identical should be 1")
	}
	if !approx(ShannonEntropy([]float64{0.5, 0.5}), math.Log(2), tol) {
		t.Errorf("entropy of fair coin should be ln2")
	}
}

func TestSoftmaxLogSumExp(t *testing.T) {
	if !approx(LogSumExp([]float64{0, 0}), math.Log(2), 1e-12) {
		t.Errorf("LogSumExp([0,0]) = %v", LogSumExp([]float64{0, 0}))
	}
	sm := Softmax([]float64{0, 0})
	if !VecEqual(sm, []float64{0.5, 0.5}, 1e-12) {
		t.Errorf("softmax = %v", sm)
	}
	if !approx(VecSum(Softmax([]float64{1, 2, 3})), 1, 1e-12) {
		t.Errorf("softmax should sum to 1")
	}
}

// ---------------------------------------------------------------------------
// MarkovChain core
// ---------------------------------------------------------------------------

func twoState(t *testing.T) *MarkovChain {
	t.Helper()
	c, err := NewMarkovChain([][]float64{{0.9, 0.1}, {0.5, 0.5}})
	if err != nil {
		t.Fatal(err)
	}
	return c
}

func TestNStepDistribution(t *testing.T) {
	c := twoState(t)
	d := c.StepDistribution([]float64{1, 0}, 1)
	if !VecEqual(d, []float64{0.9, 0.1}, 1e-12) {
		t.Errorf("1-step dist = %v", d)
	}
	d2 := c.StepDistribution([]float64{1, 0}, 2)
	if !VecEqual(d2, []float64{0.86, 0.14}, 1e-12) {
		t.Errorf("2-step dist = %v", d2)
	}
	if !approx(c.NStepProb(0, 0, 2), 0.86, 1e-12) {
		t.Errorf("NStepProb = %v", c.NStepProb(0, 0, 2))
	}
}

func TestStationaryTwoState(t *testing.T) {
	c := twoState(t)
	pi, err := c.StationaryDistribution()
	if err != nil {
		t.Fatal(err)
	}
	want := []float64{5.0 / 6.0, 1.0 / 6.0}
	if !VecEqual(pi, want, 1e-9) {
		t.Errorf("stationary = %v, want %v", pi, want)
	}
	// π must be a fixed point of P.
	if !VecEqual(c.NextDistribution(pi), pi, 1e-9) {
		t.Errorf("π is not stationary")
	}
	// Power iteration should agree.
	pip, err := c.StationaryPower(1e-14, 100000)
	if err != nil {
		t.Fatal(err)
	}
	if !VecEqual(pip, want, 1e-6) {
		t.Errorf("power stationary = %v, want %v", pip, want)
	}
}

func TestReversibleAndEntropyRate(t *testing.T) {
	c := twoState(t)
	pi, _ := c.StationaryDistribution()
	// Every 2-state ergodic chain is reversible w.r.t. its stationary dist.
	if !c.IsReversible(pi, 1e-9) {
		t.Errorf("2-state chain should be reversible")
	}
	er := c.EntropyRate(pi)
	if !approx(er, 0.3864279, 1e-5) {
		t.Errorf("entropy rate = %v, want ~0.3864", er)
	}
}

func TestMeanFirstPassage(t *testing.T) {
	c := twoState(t)
	m, err := c.MeanFirstPassageTimes()
	if err != nil {
		t.Fatal(err)
	}
	want := [][]float64{{1.2, 10}, {2, 6}}
	if !MatEqual(m, want, 1e-7) {
		t.Errorf("MFPT = %v, want %v", m, want)
	}
	mr, _ := c.MeanRecurrenceTimes()
	if !VecEqual(mr, []float64{1.2, 6}, 1e-7) {
		t.Errorf("mean recurrence = %v", mr)
	}
	k, err := c.KemenyConstant()
	if err != nil {
		t.Fatal(err)
	}
	if !approx(k, 5.0/3.0, 1e-6) {
		t.Errorf("Kemeny = %v, want 1.6667", k)
	}
}

// ---------------------------------------------------------------------------
// Absorbing chains — gambler's ruin on {0,1,2,3,4}
// ---------------------------------------------------------------------------

func gamblersRuin(t *testing.T) *MarkovChain {
	t.Helper()
	p := [][]float64{
		{1, 0, 0, 0, 0},
		{0.5, 0, 0.5, 0, 0},
		{0, 0.5, 0, 0.5, 0},
		{0, 0, 0.5, 0, 0.5},
		{0, 0, 0, 0, 1},
	}
	c, err := NewMarkovChain(p)
	if err != nil {
		t.Fatal(err)
	}
	return c
}

func TestAbsorbing(t *testing.T) {
	c := gamblersRuin(t)
	if !c.IsAbsorbing() {
		t.Errorf("should be absorbing")
	}
	abs := c.AbsorbingStates()
	if len(abs) != 2 || abs[0] != 0 || abs[1] != 4 {
		t.Errorf("absorbing states = %v, want [0 4]", abs)
	}
	trans := c.TransientStates()
	if len(trans) != 3 || trans[0] != 1 || trans[2] != 3 {
		t.Errorf("transient = %v, want [1 2 3]", trans)
	}

	steps, transient, err := c.ExpectedStepsToAbsorption()
	if err != nil {
		t.Fatal(err)
	}
	if !VecEqual(transient2float(transient), []float64{1, 2, 3}, tol) {
		t.Errorf("transient order = %v", transient)
	}
	if !VecEqual(steps, []float64{3, 4, 3}, 1e-9) {
		t.Errorf("expected steps = %v, want [3 4 3]", steps)
	}

	b, _, absorbing, err := c.AbsorptionProbabilities()
	if err != nil {
		t.Fatal(err)
	}
	if absorbing[0] != 0 || absorbing[1] != 4 {
		t.Errorf("absorbing order = %v", absorbing)
	}
	wantB := [][]float64{{0.75, 0.25}, {0.5, 0.5}, {0.25, 0.75}}
	if !MatEqual(b, wantB, 1e-9) {
		t.Errorf("absorption probs = %v, want %v", b, wantB)
	}

	// Convenience accessor.
	p14, err := c.AbsorptionProbability(1, 4)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(p14, 0.25, 1e-9) {
		t.Errorf("P(absorb at 4 | start 1) = %v, want 0.25", p14)
	}

	// Fundamental matrix rows should sum to expected steps.
	n, _, err := c.FundamentalMatrix()
	if err != nil {
		t.Fatal(err)
	}
	for i := range n {
		if !approx(VecSum(n[i]), steps[i], 1e-9) {
			t.Errorf("row %d of N sums to %v, want %v", i, VecSum(n[i]), steps[i])
		}
	}
}

func transient2float(s []int) []float64 {
	out := make([]float64, len(s))
	for i, v := range s {
		out[i] = float64(v)
	}
	return out
}

func TestHittingProbabilities(t *testing.T) {
	c := gamblersRuin(t)
	h, err := c.HittingProbabilities([]int{4})
	if err != nil {
		t.Fatal(err)
	}
	// Probability of ever hitting state 4 = i/4.
	want := []float64{0, 0.25, 0.5, 0.75, 1}
	if !VecEqual(h, want, 1e-9) {
		t.Errorf("hitting probs = %v, want %v", h, want)
	}
	et, err := c.ExpectedHittingTimes([]int{0, 4})
	if err != nil {
		t.Fatal(err)
	}
	if !VecEqual(et, []float64{0, 3, 4, 3, 0}, 1e-9) {
		t.Errorf("expected hitting times = %v", et)
	}
}

// ---------------------------------------------------------------------------
// Classification
// ---------------------------------------------------------------------------

func TestClassification(t *testing.T) {
	p := [][]float64{
		{0, 1, 0, 0},
		{1, 0, 0, 0},
		{0, 0, 0.5, 0.5},
		{0, 0, 0, 1},
	}
	c, err := NewMarkovChain(p)
	if err != nil {
		t.Fatal(err)
	}
	if c.IsIrreducible() {
		t.Errorf("chain should be reducible")
	}
	classes := c.CommunicatingClasses()
	if len(classes) != 3 {
		t.Fatalf("expected 3 classes, got %v", classes)
	}
	// {0,1} closed recurrent, period 2.
	if c.Period(0) != 2 {
		t.Errorf("period of state 0 = %d, want 2", c.Period(0))
	}
	if c.IsAperiodic() {
		t.Errorf("chain should be periodic")
	}
	rec := c.RecurrentStates()
	if !VecEqual(transient2float(rec), []float64{0, 1, 3}, tol) {
		t.Errorf("recurrent = %v, want [0 1 3]", rec)
	}
	trans := c.TransientStates()
	if len(trans) != 1 || trans[0] != 2 {
		t.Errorf("transient = %v, want [2]", trans)
	}
	labels := c.ClassifyStates()
	if labels[3] != "absorbing" || labels[2] != "transient" || labels[0] != "recurrent" {
		t.Errorf("labels = %v", labels)
	}
	if !c.Communicates(0, 1) || c.Communicates(0, 2) {
		t.Errorf("communication relation wrong")
	}
}

func TestRegularErgodic(t *testing.T) {
	c := twoState(t)
	if !c.IsIrreducible() || !c.IsAperiodic() || !c.IsErgodic() || !c.IsRegular() {
		t.Errorf("two-state chain should be ergodic and regular")
	}
	// A pure 2-cycle is irreducible but periodic, hence not regular.
	cyc, _ := NewMarkovChain([][]float64{{0, 1}, {1, 0}})
	if cyc.IsAperiodic() || cyc.IsRegular() {
		t.Errorf("2-cycle should be periodic and not regular")
	}
	if !cyc.IsIrreducible() {
		t.Errorf("2-cycle should be irreducible")
	}
}

// ---------------------------------------------------------------------------
// HMM
// ---------------------------------------------------------------------------

func testHMM(t *testing.T) *HMM {
	t.Helper()
	a := [][]float64{{0.7, 0.3}, {0.4, 0.6}}
	b := [][]float64{{0.9, 0.1}, {0.2, 0.8}}
	pi := []float64{0.6, 0.4}
	h, err := NewHMM(a, b, pi)
	if err != nil {
		t.Fatal(err)
	}
	return h
}

// bruteForceLikelihood enumerates all state paths.
func bruteForceLikelihood(h *HMM, obs []int) float64 {
	N := h.n
	T := len(obs)
	var total float64
	path := make([]int, T)
	var rec func(t int)
	rec = func(t int) {
		if t == T {
			p := h.Pi[path[0]] * h.B[path[0]][obs[0]]
			for k := 1; k < T; k++ {
				p *= h.A[path[k-1]][path[k]] * h.B[path[k]][obs[k]]
			}
			total += p
			return
		}
		for s := 0; s < N; s++ {
			path[t] = s
			rec(t + 1)
		}
	}
	rec(0)
	return total
}

func bruteForceViterbi(h *HMM, obs []int) ([]int, float64) {
	N := h.n
	T := len(obs)
	best := math.Inf(-1)
	var bestPath []int
	path := make([]int, T)
	var rec func(t int)
	rec = func(t int) {
		if t == T {
			lp := logSafe(h.Pi[path[0]]) + logSafe(h.B[path[0]][obs[0]])
			for k := 1; k < T; k++ {
				lp += logSafe(h.A[path[k-1]][path[k]]) + logSafe(h.B[path[k]][obs[k]])
			}
			if lp > best {
				best = lp
				bestPath = append([]int{}, path...)
			}
			return
		}
		for s := 0; s < N; s++ {
			path[t] = s
			rec(t + 1)
		}
	}
	rec(0)
	return bestPath, best
}

func TestHMMForwardBackward(t *testing.T) {
	h := testHMM(t)
	obs := []int{0, 1, 0, 0, 1}

	bf := bruteForceLikelihood(h, obs)
	lik, err := h.Likelihood(obs)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(lik, bf, 1e-10) {
		t.Errorf("Likelihood = %v, brute force = %v", lik, bf)
	}
	ll, _ := h.LogLikelihood(obs)
	if !approx(ll, math.Log(bf), 1e-9) {
		t.Errorf("LogLikelihood = %v, want %v", ll, math.Log(bf))
	}
	// Unscaled forward total must also match.
	_, likU, _ := h.Forward(obs)
	if !approx(likU, bf, 1e-10) {
		t.Errorf("unscaled forward = %v, want %v", likU, bf)
	}
	// Posterior marginals sum to 1 at each time.
	gamma, _ := h.PosteriorMarginals(obs)
	for t0 := range gamma {
		if !approx(VecSum(gamma[t0]), 1, 1e-9) {
			t.Errorf("gamma row %d sums to %v", t0, VecSum(gamma[t0]))
		}
	}
}

func TestHMMViterbi(t *testing.T) {
	h := testHMM(t)
	obs := []int{0, 1, 0, 0, 1, 1, 0}
	path, lp, err := h.Viterbi(obs)
	if err != nil {
		t.Fatal(err)
	}
	bfPath, bfLP := bruteForceViterbi(h, obs)
	if !approx(lp, bfLP, 1e-9) {
		t.Errorf("Viterbi logprob = %v, brute = %v", lp, bfLP)
	}
	for i := range path {
		if path[i] != bfPath[i] {
			t.Errorf("Viterbi path %v != brute %v", path, bfPath)
			break
		}
	}
}

func TestHMMPredict(t *testing.T) {
	h := testHMM(t)
	obs := []int{0, 1, 0}
	pn, err := h.PredictNextObservationDistribution(obs)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(VecSum(pn), 1, 1e-9) {
		t.Errorf("predictive dist sums to %v", VecSum(pn))
	}
	filt, _ := h.Filter(obs)
	if !approx(VecSum(filt), 1, 1e-9) {
		t.Errorf("filter dist sums to %v", VecSum(filt))
	}
}

func TestBaumWelchIncreasesLikelihood(t *testing.T) {
	// Start from a deliberately mediocre model and train on data generated by a
	// crisp model; the log-likelihood must not decrease across iterations.
	rng := rand.New(rand.NewSource(7))
	trueA := [][]float64{{0.8, 0.2}, {0.3, 0.7}}
	trueB := [][]float64{{0.9, 0.1}, {0.15, 0.85}}
	truePi := []float64{0.5, 0.5}
	gen, _ := NewHMM(trueA, trueB, truePi)
	var seqs [][]int
	for i := 0; i < 20; i++ {
		_, o := gen.Generate(30, rng)
		seqs = append(seqs, o)
	}

	initA := [][]float64{{0.6, 0.4}, {0.4, 0.6}}
	initB := [][]float64{{0.7, 0.3}, {0.3, 0.7}}
	initPi := []float64{0.5, 0.5}
	h, _ := NewHMM(initA, initB, initPi)

	trained, history, err := h.BaumWelchMultiple(seqs, 60, 1e-9)
	if err != nil {
		t.Fatal(err)
	}
	for i := 1; i < len(history); i++ {
		if history[i] < history[i-1]-1e-6 {
			t.Errorf("log-likelihood decreased at iter %d: %v -> %v", i, history[i-1], history[i])
		}
	}
	// Trained model should assign higher likelihood than the initial one.
	var initLL, trainLL float64
	for _, o := range seqs {
		l0, _ := h.LogLikelihood(o)
		l1, _ := trained.LogLikelihood(o)
		initLL += l0
		trainLL += l1
	}
	if trainLL < initLL {
		t.Errorf("training did not improve likelihood: %v -> %v", initLL, trainLL)
	}
}

// ---------------------------------------------------------------------------
// MCMC
// ---------------------------------------------------------------------------

func TestRandomWalkMetropolisNormal(t *testing.T) {
	rng := rand.New(rand.NewSource(42))
	logTarget := func(x float64) float64 { return -0.5 * x * x } // standard normal
	samples, accRate := RandomWalkMetropolis(logTarget, 0, 2.4, 200000, rng)
	if accRate <= 0 || accRate >= 1 {
		t.Errorf("acceptance rate = %v out of range", accRate)
	}
	burned := samples[20000:]
	mean := SampleMean(burned)
	varr := SampleVariance(burned)
	if math.Abs(mean) > 0.05 {
		t.Errorf("sample mean = %v, want ~0", mean)
	}
	if math.Abs(varr-1) > 0.1 {
		t.Errorf("sample variance = %v, want ~1", varr)
	}
}

func TestGibbsSampler(t *testing.T) {
	// Sample from a bivariate normal with correlation rho via Gibbs.
	rng := rand.New(rand.NewSource(1))
	rho := 0.6
	conds := []func([]float64, *rand.Rand) float64{
		func(x []float64, r *rand.Rand) float64 { return rho*x[1] + math.Sqrt(1-rho*rho)*r.NormFloat64() },
		func(x []float64, r *rand.Rand) float64 { return rho*x[0] + math.Sqrt(1-rho*rho)*r.NormFloat64() },
	}
	chain := GibbsSampler(conds, []float64{0, 0}, 100000, rng)
	chain = chain[5000:]
	cov := CovarianceMatrix(chain)
	if math.Abs(cov[0][0]-1) > 0.1 || math.Abs(cov[1][1]-1) > 0.1 {
		t.Errorf("marginal variances off: %v", cov)
	}
	if math.Abs(cov[0][1]-rho) > 0.05 {
		t.Errorf("correlation = %v, want %v", cov[0][1], rho)
	}
}

func TestDiagnostics(t *testing.T) {
	if !approx(Autocorrelation([]float64{1, 2, 3, 4}, 0), 1, tol) {
		t.Errorf("lag-0 autocorrelation must be 1")
	}
	// Two well-mixed chains from the same distribution -> R-hat close to 1.
	rng := rand.New(rand.NewSource(3))
	var chains [][]float64
	for c := 0; c < 4; c++ {
		ch := make([]float64, 5000)
		for i := range ch {
			ch[i] = rng.NormFloat64()
		}
		chains = append(chains, ch)
	}
	rhat := GelmanRubin(chains)
	if math.Abs(rhat-1) > 0.05 {
		t.Errorf("R-hat = %v, want ~1", rhat)
	}
	ess := EffectiveSampleSize(chains[0])
	if ess <= 0 || ess > float64(len(chains[0])+1) {
		t.Errorf("ESS out of range: %v", ess)
	}
}

// ---------------------------------------------------------------------------
// Sampling
// ---------------------------------------------------------------------------

func TestSampleCategoricalFrequencies(t *testing.T) {
	rng := rand.New(rand.NewSource(9))
	p := []float64{0.1, 0.3, 0.6}
	counts := make([]int, 3)
	n := 200000
	for i := 0; i < n; i++ {
		counts[SampleCategorical(p, rng)]++
	}
	for k := range p {
		f := float64(counts[k]) / float64(n)
		if math.Abs(f-p[k]) > 0.01 {
			t.Errorf("category %d freq = %v, want %v", k, f, p[k])
		}
	}
}

func TestAliasTable(t *testing.T) {
	rng := rand.New(rand.NewSource(11))
	w := []float64{1, 2, 3, 4}
	at := NewAliasTable(w)
	if at.Len() != 4 {
		t.Fatalf("alias len = %d", at.Len())
	}
	counts := make([]int, 4)
	n := 200000
	for i := 0; i < n; i++ {
		counts[at.Sample(rng)]++
	}
	total := 10.0
	for k := range w {
		f := float64(counts[k]) / float64(n)
		want := w[k] / total
		if math.Abs(f-want) > 0.01 {
			t.Errorf("alias category %d freq = %v, want %v", k, f, want)
		}
	}
}

func TestSimulate(t *testing.T) {
	c := twoState(t)
	rng := rand.New(rand.NewSource(5))
	path := c.Simulate(0, 100000, rng)
	if len(path) != 100001 {
		t.Fatalf("path length = %d", len(path))
	}
	// Empirical occupation should approximate stationary [5/6, 1/6].
	count0 := 0
	for _, s := range path {
		if s == 0 {
			count0++
		}
	}
	frac := float64(count0) / float64(len(path))
	if math.Abs(frac-5.0/6.0) > 0.01 {
		t.Errorf("occupation of state 0 = %v, want ~0.8333", frac)
	}
}

// ---------------------------------------------------------------------------
// Kronecker / extra
// ---------------------------------------------------------------------------

func TestKronecker(t *testing.T) {
	a := [][]float64{{1, 0}, {0, 1}}
	b := [][]float64{{0, 1}, {1, 0}}
	k := KroneckerProduct(a, b)
	if len(k) != 4 || len(k[0]) != 4 {
		t.Fatalf("kron shape wrong: %dx%d", len(k), len(k[0]))
	}
	// The Kronecker product of two stochastic matrices is stochastic.
	if !IsStochastic(k, 1e-12) {
		t.Errorf("kron of stochastic matrices should be stochastic")
	}
}

func TestNewMarkovChainFromCounts(t *testing.T) {
	counts := [][]float64{{2, 2, 0}, {0, 0, 0}, {1, 1, 2}}
	c, err := NewMarkovChainFromCounts(counts)
	if err != nil {
		t.Fatal(err)
	}
	m := c.Matrix()
	if !VecEqual(m[0], []float64{0.5, 0.5, 0}, tol) {
		t.Errorf("row 0 = %v", m[0])
	}
	// Zero row becomes a self-loop (absorbing).
	if !VecEqual(m[1], []float64{0, 1, 0}, tol) {
		t.Errorf("row 1 = %v, want self-loop", m[1])
	}
	if !VecEqual(m[2], []float64{0.25, 0.25, 0.5}, tol) {
		t.Errorf("row 2 = %v", m[2])
	}
}

// ---------------------------------------------------------------------------
// Runnable example
// ---------------------------------------------------------------------------

func ExampleMarkovChain_StationaryDistribution() {
	// A simple weather chain: state 0 = sunny, state 1 = rainy.
	c, _ := NewMarkovChain([][]float64{
		{0.9, 0.1},
		{0.5, 0.5},
	})
	pi, _ := c.StationaryDistribution()
	fmt.Printf("sunny=%.4f rainy=%.4f\n", pi[0], pi[1])
	// Output: sunny=0.8333 rainy=0.1667
}
