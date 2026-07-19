package stochastic

import (
	"fmt"
	"math"
	"math/rand"
	"testing"
)

const tol = 1e-9

func approx(a, b, eps float64) bool {
	return math.Abs(a-b) <= eps
}

// meanOf returns the mean of a slice of float64.
func meanOf(x []float64) float64 {
	s := 0.0
	for _, v := range x {
		s += v
	}
	return s / float64(len(x))
}

func TestUniformSampleRange(t *testing.T) {
	rng := NewRNG(1)
	for i := 0; i < 10000; i++ {
		v := UniformSample(rng, 2, 5)
		if v < 2 || v >= 5 {
			t.Fatalf("uniform out of range: %v", v)
		}
	}
}

func TestDeterministicSeed(t *testing.T) {
	a := BrownianMotion(NewRNG(42), 1, 100)
	b := BrownianMotion(NewRNG(42), 1, 100)
	for i := range a.Values {
		if a.Values[i] != b.Values[i] {
			t.Fatalf("non-deterministic at %d: %v vs %v", i, a.Values[i], b.Values[i])
		}
	}
}

func TestSamplerMeans(t *testing.T) {
	rng := NewRNG(7)
	const n = 200000
	tests := []struct {
		name string
		draw func() float64
		mean float64
		eps  float64
	}{
		{"exponential", func() float64 { return ExponentialSample(rng, 2) }, 0.5, 0.02},
		{"normal", func() float64 { return NormalSample(rng, 3, 2) }, 3, 0.03},
		{"gamma", func() float64 { return GammaSample(rng, 2, 3) }, 6, 0.1},
		{"beta", func() float64 { return BetaSample(rng, 2, 3) }, 0.4, 0.01},
		{"chi2", func() float64 { return ChiSquaredSample(rng, 4) }, 4, 0.05},
		{"weibull", func() float64 { return WeibullSample(rng, 1, 2) }, 2, 0.05},
		{"rayleigh", func() float64 { return RayleighSample(rng, 1) }, math.Sqrt(math.Pi / 2), 0.02},
		{"lognormal", func() float64 { return LogNormalSample(rng, 0, 1) }, math.Exp(0.5), 0.05},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			s := 0.0
			for i := 0; i < n; i++ {
				s += tc.draw()
			}
			m := s / n
			if !approx(m, tc.mean, tc.eps) {
				t.Errorf("%s mean = %v, want ~%v", tc.name, m, tc.mean)
			}
		})
	}
}

func TestPoissonSample(t *testing.T) {
	rng := NewRNG(11)
	for _, lambda := range []float64{0.5, 3, 12, 50, 200} {
		const n = 200000
		s := 0.0
		ss := 0.0
		for i := 0; i < n; i++ {
			k := float64(PoissonSample(rng, lambda))
			s += k
			ss += k * k
		}
		m := s / n
		v := ss/n - m*m
		if !approx(m, lambda, 0.05*lambda+0.05) {
			t.Errorf("lambda=%v mean=%v", lambda, m)
		}
		if !approx(v, lambda, 0.1*lambda+0.2) {
			t.Errorf("lambda=%v var=%v", lambda, v)
		}
	}
}

func TestPoissonPMFCDF(t *testing.T) {
	// Poisson(2): P(0)=e^-2, P(1)=2e^-2, P(2)=2e^-2
	e2 := math.Exp(-2)
	tests := []struct {
		k    int
		want float64
	}{
		{0, e2},
		{1, 2 * e2},
		{2, 2 * e2},
	}
	for _, tc := range tests {
		if got := PoissonPMF(tc.k, 2); !approx(got, tc.want, 1e-12) {
			t.Errorf("PMF(%d)=%v want %v", tc.k, got, tc.want)
		}
	}
	// CDF must be monotone and approach 1.
	if PoissonCDF(0, 2) >= PoissonCDF(3, 2) {
		t.Error("CDF not increasing")
	}
	if got := PoissonCDF(100, 5); !approx(got, 1, 1e-9) {
		t.Errorf("CDF tail = %v", got)
	}
}

func TestPoissonProcessCount(t *testing.T) {
	rng := NewRNG(3)
	const n = 50000
	s := 0
	for i := 0; i < n; i++ {
		s += PoissonProcessCount(rng, 4, 2.5) // mean 10
	}
	m := float64(s) / n
	if !approx(m, 10, 0.1) {
		t.Errorf("process count mean = %v want 10", m)
	}
}

func TestPoissonArrivalsSorted(t *testing.T) {
	rng := NewRNG(5)
	arr := PoissonArrivals(rng, 5, 10)
	for i := 1; i < len(arr); i++ {
		if arr[i] < arr[i-1] {
			t.Fatal("arrivals not sorted")
		}
	}
	if len(arr) == 0 {
		t.Fatal("expected some arrivals")
	}
}

func TestSuperposedPoisson(t *testing.T) {
	a := []float64{0.1, 0.5, 0.9}
	b := []float64{0.2, 0.3, 1.0}
	m := SuperposedPoisson(a, b)
	want := []float64{0.1, 0.2, 0.3, 0.5, 0.9, 1.0}
	if len(m) != len(want) {
		t.Fatalf("len = %d", len(m))
	}
	for i := range want {
		if m[i] != want[i] {
			t.Errorf("merge[%d]=%v want %v", i, m[i], want[i])
		}
	}
}

func TestCompoundPoissonMoments(t *testing.T) {
	rng := NewRNG(9)
	rate, T := 3.0, 4.0
	jump := func(r *rand.Rand) float64 { return NormalSample(r, 2, 1) }
	const n = 40000
	vals := make([]float64, n)
	for i := 0; i < n; i++ {
		vals[i] = CompoundPoisson(rng, rate, T, jump)
	}
	m := meanOf(vals)
	wantMean := CompoundPoissonMean(rate, T, 2) // 24
	if !approx(m, wantMean, 0.3) {
		t.Errorf("compound mean = %v want ~%v", m, wantMean)
	}
	v := SampleVariance(vals)
	wantVar := CompoundPoissonVariance(rate, T, 2, 1) // 3*4*(1+4)=60
	if !approx(v, wantVar, 3) {
		t.Errorf("compound var = %v want ~%v", v, wantVar)
	}
}

func TestGamblersRuin(t *testing.T) {
	tests := []struct {
		i, N     int
		p        float64
		wantRuin float64
	}{
		{5, 10, 0.5, 0.5},
		{1, 4, 0.5, 0.75},
		{3, 4, 0.5, 0.25},
	}
	for _, tc := range tests {
		got := GamblersRuinProbability(tc.i, tc.N, tc.p)
		if !approx(got, tc.wantRuin, 1e-12) {
			t.Errorf("ruin(%d,%d,%v)=%v want %v", tc.i, tc.N, tc.p, got, tc.wantRuin)
		}
	}
	// biased: p=0.6, i=1, N=2. r=q/p=2/3. win=(1-r)/(1-r^2)=(1/3)/(5/9)=3/5. ruin=2/5.
	if got := GamblersRuinProbability(1, 2, 0.6); !approx(got, 0.4, 1e-12) {
		t.Errorf("biased ruin = %v want 0.4", got)
	}
	// duration fair: i(N-i) = 5*5 = 25
	if got := GamblersRuinExpectedDuration(5, 10, 0.5); !approx(got, 25, 1e-12) {
		t.Errorf("duration = %v want 25", got)
	}
}

func TestGamblersRuinSimulation(t *testing.T) {
	rng := NewRNG(13)
	const trials = 20000
	ruin := 0
	for i := 0; i < trials; i++ {
		pos, _ := AbsorbingRandomWalk(rng, 3, 0, 10, 0.5)
		if pos == 0 {
			ruin++
		}
	}
	p := float64(ruin) / trials
	if !approx(p, 0.7, 0.02) {
		t.Errorf("simulated ruin = %v want ~0.7", p)
	}
}

func TestRandomWalkMoments(t *testing.T) {
	if got := RandomWalkMean(100, 0.5); got != 0 {
		t.Errorf("mean = %v", got)
	}
	if got := RandomWalkVariance(100, 0.5); got != 100 {
		t.Errorf("var = %v", got)
	}
	rng := NewRNG(17)
	const trials = 20000
	ends := make([]float64, trials)
	for i := 0; i < trials; i++ {
		w := SimpleRandomWalk(rng, 100, 0.5)
		ends[i] = float64(w[len(w)-1])
	}
	if v := SampleVariance(ends); !approx(v, 100, 5) {
		t.Errorf("endpoint var = %v want ~100", v)
	}
}

func TestBrownianMotionStats(t *testing.T) {
	rng := NewRNG(21)
	const trials = 20000
	finals := make([]float64, trials)
	for i := 0; i < trials; i++ {
		bm := BrownianMotion(rng, 1, 50)
		finals[i] = bm.Final()
		if bm.Times[0] != 0 || !approx(bm.EndTime(), 1, tol) {
			t.Fatal("bad grid")
		}
	}
	if m := meanOf(finals); !approx(m, 0, 0.03) {
		t.Errorf("W(1) mean = %v", m)
	}
	if v := SampleVariance(finals); !approx(v, 1, 0.05) {
		t.Errorf("W(1) var = %v want 1", v)
	}
}

func TestBrownianBridgePinned(t *testing.T) {
	rng := NewRNG(23)
	for i := 0; i < 50; i++ {
		b := BrownianBridge(rng, 1, 4, 2, 40)
		if !approx(b.Initial(), 1, 1e-9) {
			t.Fatalf("bridge start = %v", b.Initial())
		}
		if !approx(b.Final(), 4, 1e-9) {
			t.Fatalf("bridge end = %v", b.Final())
		}
	}
}

func TestGBMMoments(t *testing.T) {
	rng := NewRNG(29)
	s0, mu, sigma, T := 100.0, 0.05, 0.2, 1.0
	const trials = 40000
	finals := make([]float64, trials)
	for i := 0; i < trials; i++ {
		finals[i] = GeometricBrownianMotion(rng, s0, mu, sigma, T, 20).Final()
	}
	m := meanOf(finals)
	wantMean := GBMMean(s0, mu, T)
	if !approx(m, wantMean, 1.5) {
		t.Errorf("GBM mean = %v want %v", m, wantMean)
	}
	v := SampleVariance(finals)
	wantVar := GBMVariance(s0, mu, sigma, T)
	if math.Abs(v-wantVar)/wantVar > 0.1 {
		t.Errorf("GBM var = %v want %v", v, wantVar)
	}
}

func TestOUStationary(t *testing.T) {
	rng := NewRNG(31)
	theta, mu, sigma := 1.0, 5.0, 2.0
	// run long enough to reach stationarity, sample endpoint
	const trials = 30000
	finals := make([]float64, trials)
	for i := 0; i < trials; i++ {
		finals[i] = OrnsteinUhlenbeck(rng, 0, theta, mu, sigma, 10, 200).Final()
	}
	if m := meanOf(finals); !approx(m, mu, 0.05) {
		t.Errorf("OU mean = %v want %v", m, mu)
	}
	wantVar := OUStationaryVariance(theta, sigma) // 4/2 = 2
	if v := SampleVariance(finals); !approx(v, wantVar, 0.1) {
		t.Errorf("OU var = %v want %v", v, wantVar)
	}
}

func TestOUAnalyticFormulas(t *testing.T) {
	// OUMean(x0=0,theta=1,mu=5,t=0)=0; t->inf ~5
	if got := OUMean(0, 1, 5, 0); !approx(got, 0, tol) {
		t.Errorf("OUMean t=0 = %v", got)
	}
	if got := OUMean(0, 1, 5, 100); !approx(got, 5, 1e-9) {
		t.Errorf("OUMean t=inf = %v", got)
	}
	if got := OUStationaryVariance(2, 2); !approx(got, 1, tol) {
		t.Errorf("stationary var = %v want 1", got)
	}
}

func TestEulerMaruyamaGBM(t *testing.T) {
	// dX = mu X dt + sigma X dW should track exact GBM mean.
	rng := NewRNG(37)
	mu, sigma, s0, T := 0.1, 0.3, 50.0, 1.0
	sde := SDE{
		Drift:     func(t, x float64) float64 { return mu * x },
		Diffusion: func(t, x float64) float64 { return sigma * x },
	}
	const trials = 40000
	finals := make([]float64, trials)
	for i := 0; i < trials; i++ {
		finals[i] = EulerMaruyama(rng, sde, s0, T, 100).Final()
	}
	m := meanOf(finals)
	want := GBMMean(s0, mu, T)
	if math.Abs(m-want)/want > 0.02 {
		t.Errorf("EM GBM mean = %v want %v", m, want)
	}
}

func TestMilsteinMatchesExact(t *testing.T) {
	// For GBM the Milstein scheme is strongly convergent; check its mean.
	rng := NewRNG(41)
	mu, sigma, s0, T := 0.05, 0.4, 20.0, 1.0
	sde := SDE{
		Drift:     func(t, x float64) float64 { return mu * x },
		Diffusion: func(t, x float64) float64 { return sigma * x },
	}
	dDiff := func(t, x float64) float64 { return sigma }
	const trials = 30000
	finals := make([]float64, trials)
	for i := 0; i < trials; i++ {
		finals[i] = Milstein(rng, sde, dDiff, s0, T, 50).Final()
	}
	m := meanOf(finals)
	want := GBMMean(s0, mu, T)
	if math.Abs(m-want)/want > 0.03 {
		t.Errorf("Milstein mean = %v want %v", m, want)
	}
}

func TestMonteCarloOption(t *testing.T) {
	// Black-Scholes call price via risk-neutral GBM SDE.
	r, sigma, s0, K, T := 0.05, 0.2, 100.0, 100.0, 1.0
	sde := SDE{
		Drift:     func(t, x float64) float64 { return r * x },
		Diffusion: func(t, x float64) float64 { return sigma * x },
	}
	payoff := func(p Path) float64 {
		v := p.Final() - K
		if v < 0 {
			v = 0
		}
		return math.Exp(-r*T) * v
	}
	mean, se := MonteCarloExpectation(101, sde, s0, T, 50, 40000, payoff)
	// Analytic BS call price ~ 10.4506
	bs := blackScholesCall(s0, K, r, sigma, T)
	if math.Abs(mean-bs) > 4*se+0.3 {
		t.Errorf("MC price = %v (se %v) want ~%v", mean, se, bs)
	}
	if se <= 0 {
		t.Error("standard error should be positive")
	}
}

func blackScholesCall(s, k, r, sigma, T float64) float64 {
	d1 := (math.Log(s/k) + (r+0.5*sigma*sigma)*T) / (sigma * math.Sqrt(T))
	d2 := d1 - sigma*math.Sqrt(T)
	return s*normalCDF(d1) - k*math.Exp(-r*T)*normalCDF(d2)
}

func TestReflectionPrinciple(t *testing.T) {
	// P(max >= a) should match simulated probability.
	rng := NewRNG(43)
	a, T, sigma := 1.0, 1.0, 1.0
	analytic := ReflectionPrincipleMaxProb(a, T, sigma)
	const trials = 30000
	hit := 0
	for i := 0; i < trials; i++ {
		bm := BrownianMotion(rng, T, 400)
		if bm.Max() >= a {
			hit++
		}
	}
	p := float64(hit) / trials
	if math.Abs(p-analytic) > 0.03 {
		t.Errorf("reflection: sim %v analytic %v", p, analytic)
	}
}

func TestInverseGaussian(t *testing.T) {
	// mean and variance identities
	mu, lambda := 2.0, 3.0
	if got := InverseGaussianMean(mu, lambda); got != mu {
		t.Errorf("IG mean = %v", got)
	}
	if got := InverseGaussianVariance(mu, lambda); !approx(got, mu*mu*mu/lambda, tol) {
		t.Errorf("IG var = %v", got)
	}
	// CDF is monotone increasing in x and in [0,1]
	prev := 0.0
	for _, x := range []float64{0.1, 0.5, 1, 2, 5, 20} {
		c := InverseGaussianCDF(x, mu, lambda)
		if c < prev-1e-12 || c < 0 || c > 1 {
			t.Errorf("IG CDF(%v) = %v not monotone/in-range", x, c)
		}
		prev = c
	}
	// numeric integral of PDF over a fine grid should approach the CDF
	sum := 0.0
	dx := 0.001
	for x := dx / 2; x < 5; x += dx {
		sum += InverseGaussianPDF(x, mu, lambda) * dx
	}
	if !approx(sum, InverseGaussianCDF(5, mu, lambda), 0.01) {
		t.Errorf("IG pdf integral = %v vs cdf %v", sum, InverseGaussianCDF(5, mu, lambda))
	}
}

func TestFirstPassageDrift(t *testing.T) {
	// Simulate BM with drift; mean hitting time of level a should be a/drift.
	rng := NewRNG(47)
	a, drift, sigma := 5.0, 1.0, 1.0
	sde := SDE{
		Drift:     func(t, x float64) float64 { return drift },
		Diffusion: func(t, x float64) float64 { return sigma },
	}
	const trials = 8000
	sum := 0.0
	cnt := 0
	for i := 0; i < trials; i++ {
		p := EulerMaruyama(rng, sde, 0, 30, 3000)
		if th, ok := p.TimeToHit(a); ok {
			sum += th
			cnt++
		}
	}
	if cnt == 0 {
		t.Fatal("never hit level")
	}
	m := sum / float64(cnt)
	want := FirstPassageDriftMean(a, drift) // 5
	if math.Abs(m-want) > 0.5 {
		t.Errorf("mean first passage = %v want %v", m, want)
	}
}

func TestPathMethods(t *testing.T) {
	p := NewPath(
		[]float64{0, 1, 2, 3, 4},
		[]float64{0, 2, 1, 5, 3},
	)
	if p.Len() != 5 {
		t.Fatalf("len = %d", p.Len())
	}
	if p.Max() != 5 || p.Min() != 0 {
		t.Errorf("max/min = %v/%v", p.Max(), p.Min())
	}
	if p.ArgMax() != 3 || p.ArgMin() != 0 {
		t.Errorf("argmax/argmin = %d/%d", p.ArgMax(), p.ArgMin())
	}
	if !approx(p.Range(), 5, tol) {
		t.Errorf("range = %v", p.Range())
	}
	// increments: 2,-1,4,-2 ; QV = 4+1+16+4 = 25 ; TV = 2+1+4+2 = 9
	if !approx(p.QuadraticVariation(), 25, tol) {
		t.Errorf("QV = %v", p.QuadraticVariation())
	}
	if !approx(p.TotalVariation(), 9, tol) {
		t.Errorf("TV = %v", p.TotalVariation())
	}
	// max drawdown: peak 2 then down to 1 (dd=1); peak 5 then to 3 (dd=2) -> 2
	if !approx(p.MaxDrawdown(), 2, tol) {
		t.Errorf("drawdown = %v", p.MaxDrawdown())
	}
	// interpolation at t=1.5 between (1,2) and (2,1) -> 1.5
	if !approx(p.At(1.5), 1.5, tol) {
		t.Errorf("At(1.5) = %v", p.At(1.5))
	}
	if !approx(p.At(-5), 0, tol) || !approx(p.At(100), 3, tol) {
		t.Errorf("At clamp failed")
	}
	// first passage to level 5 (up from 0) at index 3, time 3
	if th, ok := p.TimeToHit(5); !ok || !approx(th, 3, tol) {
		t.Errorf("hit = %v ok=%v", th, ok)
	}
	// occupation fraction in [0,2]: values 0,2,1,3? -> 0,2,1 are in [0,2] => 3/5
	if !approx(p.OccupationFraction(0, 2), 0.6, tol) {
		t.Errorf("occ frac = %v", p.OccupationFraction(0, 2))
	}
	// antithetic reflects about start (0): 0,-2,-1,-5,-3
	anti := p.Antithetic()
	if anti.Values[1] != -2 || anti.Values[3] != -5 {
		t.Errorf("antithetic = %v", anti.Values)
	}
}

func TestPathRunning(t *testing.T) {
	p := NewPath([]float64{0, 1, 2, 3}, []float64{1, 3, 2, 5})
	rmax := p.RunningMax()
	wantMax := []float64{1, 3, 3, 5}
	rmin := p.RunningMin()
	wantMin := []float64{1, 1, 1, 1}
	for i := range wantMax {
		if rmax[i] != wantMax[i] || rmin[i] != wantMin[i] {
			t.Errorf("running at %d: max %v min %v", i, rmax[i], rmin[i])
		}
	}
	rmean := p.RunningMean()
	if !approx(rmean[3], (1+3+2+5)/4.0, tol) {
		t.Errorf("running mean = %v", rmean[3])
	}
}

func TestGillespieBirthDeath(t *testing.T) {
	// Immigration-death: dX = k - d*X. Stationary mean = k/d = Poisson(k/d).
	rng := NewRNG(53)
	net := ImmigrationDeathNetwork(10, 1) // stationary mean 10
	const trials = 3000
	sum := 0.0
	for i := 0; i < trials; i++ {
		res := Gillespie(rng, net, []int{0}, 20)
		sum += float64(res.FinalState()[0])
	}
	m := sum / trials
	if !approx(m, 10, 0.5) {
		t.Errorf("immigration-death mean = %v want ~10", m)
	}
}

func TestGillespieConservation(t *testing.T) {
	// SIR conserves total population S+I+R.
	rng := NewRNG(59)
	net := SIRNetwork(0.005, 1.0)
	init := []int{99, 1, 0}
	total := init[0] + init[1] + init[2]
	res := Gillespie(rng, net, init, 20)
	for _, st := range res.States {
		if st[0]+st[1]+st[2] != total {
			t.Fatalf("population not conserved: %v", st)
		}
		if st[0] < 0 || st[1] < 0 || st[2] < 0 {
			t.Fatalf("negative counts: %v", st)
		}
	}
}

func TestMassActionPropensity(t *testing.T) {
	// dimerization A+A -> ... has propensity k*C(x,2) = k*x*(x-1)/2
	prop := MassAction(2.0, []int{2})
	if got := prop([]int{5}); !approx(got, 2*10, tol) { // C(5,2)=10
		t.Errorf("mass action = %v want 20", got)
	}
	// bimolecular A+B has propensity k*x_A*x_B
	prop2 := MassAction(3.0, []int{1, 1})
	if got := prop2([]int{4, 5}); !approx(got, 3*20, tol) {
		t.Errorf("bimolecular = %v want 60", got)
	}
}

func TestTauLeapingMean(t *testing.T) {
	rng := NewRNG(61)
	net := ImmigrationDeathNetwork(20, 1) // stationary mean 20
	const trials = 2000
	sum := 0.0
	for i := 0; i < trials; i++ {
		res := TauLeaping(rng, net, []int{20}, 15, 0.01)
		sum += float64(res.FinalState()[0])
	}
	m := sum / trials
	if !approx(m, 20, 1.0) {
		t.Errorf("tau-leaping mean = %v want ~20", m)
	}
}

func TestEstimateGBMParams(t *testing.T) {
	rng := NewRNG(67)
	mu, sigma := 0.08, 0.25
	p := GeometricBrownianMotion(rng, 100, mu, sigma, 50, 200000)
	emu, esig := EstimateGBMParams(p)
	if math.Abs(esig-sigma) > 0.01 {
		t.Errorf("est sigma = %v want %v", esig, sigma)
	}
	if math.Abs(emu-mu) > 0.03 {
		t.Errorf("est mu = %v want %v", emu, mu)
	}
}

func TestEstimateOUParams(t *testing.T) {
	rng := NewRNG(71)
	theta, mu, sigma := 1.5, 3.0, 0.8
	p := OrnsteinUhlenbeck(rng, 3, theta, mu, sigma, 500, 100000)
	et, em, es := EstimateOUParams(p)
	if math.Abs(et-theta) > 0.15 {
		t.Errorf("est theta = %v want %v", et, theta)
	}
	if math.Abs(em-mu) > 0.1 {
		t.Errorf("est mu = %v want %v", em, mu)
	}
	if math.Abs(es-sigma) > 0.05 {
		t.Errorf("est sigma = %v want %v", es, sigma)
	}
}

func TestCorrelatedBrownianMotion(t *testing.T) {
	rng := NewRNG(73)
	corr := [][]float64{{1, 0.5}, {0.5, 1}}
	const trials = 20000
	var num, d1, d2 float64
	for i := 0; i < trials; i++ {
		paths, err := CorrelatedBrownianMotion(rng, corr, 1, 1)
		if err != nil {
			t.Fatal(err)
		}
		a := paths[0].Final()
		b := paths[1].Final()
		num += a * b
		d1 += a * a
		d2 += b * b
	}
	rho := num / math.Sqrt(d1*d2)
	if math.Abs(rho-0.5) > 0.03 {
		t.Errorf("empirical corr = %v want 0.5", rho)
	}
	// non-PD matrix should error
	if _, err := CorrelatedBrownianMotion(rng, [][]float64{{1, 2}, {2, 1}}, 1, 1); err == nil {
		t.Error("expected error for non-PD matrix")
	}
}

func TestFractionalBrownianMotion(t *testing.T) {
	rng := NewRNG(79)
	// H=0.5 should reduce to standard BM with variance t at endpoint.
	const trials = 4000
	finals := make([]float64, trials)
	for i := 0; i < trials; i++ {
		p, err := FractionalBrownianMotion(rng, 0.5, 1, 8)
		if err != nil {
			t.Fatal(err)
		}
		finals[i] = p.Final()
	}
	if v := SampleVariance(finals); !approx(v, 1, 0.1) {
		t.Errorf("fBm H=0.5 endpoint var = %v want 1", v)
	}
	if _, err := FractionalBrownianMotion(rng, 1.5, 1, 4); err == nil {
		t.Error("expected error for H out of range")
	}
}

func TestSampleStatistics(t *testing.T) {
	x := []float64{2, 4, 4, 4, 5, 5, 7, 9}
	if !approx(SampleMean(x), 5, tol) {
		t.Errorf("mean = %v", SampleMean(x))
	}
	// population var = 4, sample var = 32/7
	if !approx(SampleVariance(x), 32.0/7.0, tol) {
		t.Errorf("var = %v", SampleVariance(x))
	}
	if !approx(Autocorrelation(x, 0), 1, tol) {
		t.Errorf("autocorr lag0 = %v", Autocorrelation(x, 0))
	}
}

func TestDiscreteSample(t *testing.T) {
	rng := NewRNG(83)
	weights := []float64{1, 0, 3} // index 1 impossible
	counts := make([]int, 3)
	const n = 60000
	for i := 0; i < n; i++ {
		counts[DiscreteSample(rng, weights)]++
	}
	if counts[1] != 0 {
		t.Error("zero-weight category was sampled")
	}
	frac0 := float64(counts[0]) / n
	if !approx(frac0, 0.25, 0.02) {
		t.Errorf("category 0 frac = %v want 0.25", frac0)
	}
	if DiscreteSample(rng, []float64{0, 0}) != -1 {
		t.Error("all-zero weights should return -1")
	}
}

func ExampleGeometricBrownianMotion() {
	rng := NewRNG(1)
	path := GeometricBrownianMotion(rng, 100, 0.05, 0.2, 1, 4)
	fmt.Printf("start=%.2f end=%.4f\n", path.Initial(), path.Final())
	// Output: start=100.00 end=107.3004
}

func ExampleGillespie() {
	rng := NewRNG(7)
	net := SIRNetwork(0.01, 1.0)
	res := Gillespie(rng, net, []int{50, 1, 0}, 30)
	final := res.FinalState()
	fmt.Printf("S=%d I=%d R=%d total=%d\n", final[0], final[1], final[2], final[0]+final[1]+final[2])
	// Output: S=48 I=0 R=3 total=51
}
