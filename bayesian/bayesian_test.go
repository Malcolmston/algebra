package bayesian

import (
	"fmt"
	"math"
	"testing"
)

const tol = 1e-9

func approx(a, b, eps float64) bool {
	if math.IsInf(a, 0) || math.IsInf(b, 0) {
		return a == b
	}
	return math.Abs(a-b) <= eps*(1+math.Abs(b))
}

// ------------------------------------------------------------------
// Special functions
// ------------------------------------------------------------------

func TestSpecialFunctions(t *testing.T) {
	tests := []struct {
		name string
		got  float64
		want float64
		eps  float64
	}{
		{"LogGamma(5)=log24", LogGamma(5), math.Log(24), 1e-12},
		{"LogBeta(2,3)", LogBeta(2, 3), math.Log(1.0 / 12.0), 1e-12},
		{"BetaFunc(2,3)", BetaFunc(2, 3), 1.0 / 12.0, 1e-12},
		{"Digamma(1)=-gamma", Digamma(1), -0.5772156649015329, 1e-9},
		{"Digamma(2)=1-gamma", Digamma(2), 1 - 0.5772156649015329, 1e-9},
		{"Trigamma(1)=pi^2/6", Trigamma(1), math.Pi * math.Pi / 6, 1e-8},
		{"RegGammaP(1,x)=1-e^-x", RegularizedGammaP(1, 2), 1 - math.Exp(-2), 1e-10},
		{"RegIncBeta(0.5,1,1)", RegularizedIncompleteBeta(0.5, 1, 1), 0.5, 1e-12},
		{"RegIncBeta(x,2,2)", RegularizedIncompleteBeta(0.5, 2, 2), 0.5, 1e-10},
		{"ErfInv(erf(0.7))", ErfInv(math.Erf(0.7)), 0.7, 1e-10},
		{"StdNormalCDF(0)", StdNormalCDF(0), 0.5, 1e-12},
		{"StdNormalQuantile(0.975)", StdNormalQuantile(0.975), 1.959963984540054, 1e-8},
		{"LogChoose(5,2)=log10", LogChoose(5, 2), math.Log(10), 1e-12},
		{"LogFactorial(4)=log24", LogFactorial(4), math.Log(24), 1e-12},
	}
	for _, tt := range tests {
		if !approx(tt.got, tt.want, tt.eps) {
			t.Errorf("%s: got %.12g want %.12g", tt.name, tt.got, tt.want)
		}
	}
}

func TestInverseSpecials(t *testing.T) {
	// Round trips.
	for _, p := range []float64{0.05, 0.25, 0.5, 0.75, 0.95} {
		x := InverseRegularizedIncompleteBeta(p, 2.5, 3.5)
		back := RegularizedIncompleteBeta(x, 2.5, 3.5)
		if !approx(back, p, 1e-8) {
			t.Errorf("inv beta round trip p=%v got %v", p, back)
		}
		xg := InverseRegularizedGammaP(p, 3.2)
		backg := RegularizedGammaP(3.2, xg)
		if !approx(backg, p, 1e-8) {
			t.Errorf("inv gamma round trip p=%v got %v", p, backg)
		}
	}
}

func TestLogSumExp(t *testing.T) {
	got := LogSumExp([]float64{0, 0, 0})
	if !approx(got, math.Log(3), 1e-12) {
		t.Errorf("LogSumExp got %v", got)
	}
	if !math.IsInf(LogSumExp(nil), -1) {
		t.Errorf("LogSumExp(nil) should be -Inf")
	}
}

// ------------------------------------------------------------------
// Distributions
// ------------------------------------------------------------------

func TestBeta(t *testing.T) {
	d := Beta{Alpha: 2, Beta: 3}
	if !approx(d.Mean(), 0.4, tol) {
		t.Errorf("mean %v", d.Mean())
	}
	if !approx(d.Variance(), 0.04, tol) {
		t.Errorf("var %v", d.Variance())
	}
	if !approx(d.Mode(), 1.0/3.0, tol) {
		t.Errorf("mode %v", d.Mode())
	}
	// PDF at 0.5 for Beta(2,3): x*(1-x)^2 / B(2,3) = 0.5*0.25/(1/12)=1.5
	if !approx(d.PDF(0.5), 1.5, 1e-9) {
		t.Errorf("pdf %v", d.PDF(0.5))
	}
	// CDF and quantile round trip
	q := d.Quantile(0.3)
	if !approx(d.CDF(q), 0.3, 1e-8) {
		t.Errorf("beta cdf/quantile %v", d.CDF(q))
	}
	if !approx(d.MeanLog(), Digamma(2)-Digamma(5), 1e-12) {
		t.Errorf("meanlog")
	}
}

func TestGamma(t *testing.T) {
	d := Gamma{Shape: 3, Rate: 2}
	if !approx(d.Mean(), 1.5, tol) {
		t.Errorf("mean %v", d.Mean())
	}
	if !approx(d.Variance(), 0.75, tol) {
		t.Errorf("var %v", d.Variance())
	}
	if !approx(d.Mode(), 1.0, tol) {
		t.Errorf("mode %v", d.Mode())
	}
	q := d.Quantile(0.4)
	if !approx(d.CDF(q), 0.4, 1e-8) {
		t.Errorf("gamma cdf/quantile %v", d.CDF(q))
	}
	// Exp(rate=2) is Gamma(1,2): CDF at x = 1-e^{-2x}
	e := Gamma{Shape: 1, Rate: 2}
	if !approx(e.CDF(1), 1-math.Exp(-2), 1e-10) {
		t.Errorf("exp cdf %v", e.CDF(1))
	}
}

func TestInverseGamma(t *testing.T) {
	d := InverseGamma{Shape: 3, Scale: 2}
	if !approx(d.Mean(), 1.0, tol) {
		t.Errorf("mean %v", d.Mean())
	}
	if !approx(d.Mode(), 0.5, tol) {
		t.Errorf("mode %v", d.Mode())
	}
	q := d.Quantile(0.6)
	if !approx(d.CDF(q), 0.6, 1e-7) {
		t.Errorf("ig cdf/quantile %v", d.CDF(q))
	}
}

func TestNormal(t *testing.T) {
	d := Normal{Mu: 1, Sigma: 2}
	if !approx(d.CDF(1), 0.5, tol) {
		t.Errorf("cdf %v", d.CDF(1))
	}
	if !approx(d.Quantile(0.975), 1+2*1.959963984540054, 1e-7) {
		t.Errorf("quantile %v", d.Quantile(0.975))
	}
	if !approx(d.PDF(1), 1/(2*math.Sqrt(2*math.Pi)), 1e-12) {
		t.Errorf("pdf %v", d.PDF(1))
	}
}

func TestStudentT(t *testing.T) {
	d := StudentT{Nu: 10, Loc: 0, Scale: 1}
	// Known t_10 0.975 quantile = 2.228138852
	if !approx(d.Quantile(0.975), 2.228138852, 1e-6) {
		t.Errorf("t quantile %v", d.Quantile(0.975))
	}
	if !approx(d.CDF(0), 0.5, tol) {
		t.Errorf("t cdf %v", d.CDF(0))
	}
	// symmetric
	if !approx(d.CDF(1)+d.CDF(-1), 1, 1e-10) {
		t.Errorf("t symmetry")
	}
	if !approx(d.Variance(), 10.0/8.0, 1e-12) {
		t.Errorf("t var %v", d.Variance())
	}
}

func TestDiscrete(t *testing.T) {
	p := Poisson{Lambda: 3}
	if !approx(p.PMF(2), math.Exp(-3)*9/2, 1e-12) {
		t.Errorf("poisson pmf %v", p.PMF(2))
	}
	// Poisson CDF check against direct sum
	var s float64
	for k := 0; k <= 4; k++ {
		s += p.PMF(k)
	}
	if !approx(p.CDF(4), s, 1e-10) {
		t.Errorf("poisson cdf %v vs %v", p.CDF(4), s)
	}
	b := Binomial{N: 10, P: 0.3}
	if !approx(b.PMF(3), 120*math.Pow(0.3, 3)*math.Pow(0.7, 7), 1e-12) {
		t.Errorf("binom pmf %v", b.PMF(3))
	}
	if !approx(b.Mean(), 3, tol) || !approx(b.Variance(), 2.1, tol) {
		t.Errorf("binom moments")
	}
	// Binomial CDF vs direct
	var bs float64
	for k := 0; k <= 3; k++ {
		bs += b.PMF(k)
	}
	if !approx(b.CDF(3), bs, 1e-10) {
		t.Errorf("binom cdf %v vs %v", b.CDF(3), bs)
	}
}

func TestBetaBinomial(t *testing.T) {
	bb := BetaBinomial{N: 5, Alpha: 2, Beta: 2}
	// PMF should sum to 1
	var s float64
	for k := 0; k <= 5; k++ {
		s += bb.PMF(k)
	}
	if !approx(s, 1, 1e-10) {
		t.Errorf("betabinom sum %v", s)
	}
	if !approx(bb.Mean(), 2.5, 1e-10) {
		t.Errorf("betabinom mean %v", bb.Mean())
	}
	// Uniform prior Beta(1,1): PMF is uniform 1/(N+1)
	u := BetaBinomial{N: 4, Alpha: 1, Beta: 1}
	if !approx(u.PMF(2), 1.0/5.0, 1e-12) {
		t.Errorf("uniform betabinom %v", u.PMF(2))
	}
}

func TestNegativeBinomial(t *testing.T) {
	nb := NegativeBinomial{R: 4, P: 0.6}
	var s float64
	for k := 0; k <= 400; k++ {
		s += nb.PMF(k)
	}
	if !approx(s, 1, 1e-8) {
		t.Errorf("negbinom sum %v", s)
	}
	if !approx(nb.Mean(), 4*0.4/0.6, 1e-9) {
		t.Errorf("negbinom mean %v", nb.Mean())
	}
	if !approx(nb.Variance(), 4*0.4/(0.6*0.6), 1e-9) {
		t.Errorf("negbinom var %v", nb.Variance())
	}
}

func TestDirichlet(t *testing.T) {
	d := Dirichlet{Alpha: []float64{1, 2, 3}}
	m := d.Mean()
	want := []float64{1.0 / 6, 2.0 / 6, 3.0 / 6}
	for i := range m {
		if !approx(m[i], want[i], 1e-12) {
			t.Errorf("dir mean[%d] %v", i, m[i])
		}
	}
	// Marginal of component 2 is Beta(3,3)
	mg := d.Marginal(2)
	if !approx(mg.Alpha, 3, tol) || !approx(mg.Beta, 3, tol) {
		t.Errorf("dir marginal %v", mg)
	}
	// covariance symmetry & diagonal matches Variance()
	cov := d.Covariance()
	v := d.Variance()
	for i := range v {
		if !approx(cov[i][i], v[i], 1e-12) {
			t.Errorf("dir cov diag %d", i)
		}
	}
}

// ------------------------------------------------------------------
// Conjugate updates
// ------------------------------------------------------------------

func TestBetaBernoulliUpdate(t *testing.T) {
	post := BetaBernoulliPosterior(UniformBetaPrior(), 7, 3)
	if !approx(post.Alpha, 8, tol) || !approx(post.Beta, 4, tol) {
		t.Errorf("posterior %v", post)
	}
	if !approx(BetaBernoulliPredictiveProb(post), 8.0/12.0, tol) {
		t.Errorf("pred prob")
	}
	// sequential equals batch
	seq := BetaBernoulliSequentialPosterior(UniformBetaPrior(),
		[]bool{true, true, false, true})
	if !approx(seq.Alpha, 4, tol) || !approx(seq.Beta, 2, tol) {
		t.Errorf("seq %v", seq)
	}
}

func TestGammaPoissonUpdate(t *testing.T) {
	prior := Gamma{Shape: 2, Rate: 1}
	post := GammaPoissonPosteriorData(prior, []int{3, 4, 2, 5})
	if !approx(post.Shape, 2+14, tol) || !approx(post.Rate, 1+4, tol) {
		t.Errorf("gp posterior %v", post)
	}
	pred := GammaPoissonPredictive(post)
	if !approx(pred.Mean(), post.Mean(), 1e-9) {
		t.Errorf("gp predictive mean %v vs %v", pred.Mean(), post.Mean())
	}
}

func TestNormalKnownVarUpdate(t *testing.T) {
	prior := Normal{Mu: 0, Sigma: 1}
	data := []float64{2, 3, 4}
	post := NormalKnownVariancePosteriorData(prior, 1, data)
	// precision: prior 1, data 3 -> postvar=1/4, mean = (0 + 9)/4 = 2.25
	if !approx(post.Variance(), 0.25, 1e-12) {
		t.Errorf("nk var %v", post.Variance())
	}
	if !approx(post.Mu, 2.25, 1e-12) {
		t.Errorf("nk mean %v", post.Mu)
	}
	pred := NormalNormalPredictive(post, 1)
	if !approx(pred.Variance(), 1.25, 1e-12) {
		t.Errorf("nk pred var %v", pred.Variance())
	}
}

func TestNormalInverseGamma(t *testing.T) {
	prior := NormalInverseGamma{Mu: 0, Kappa: 1, Alpha: 2, Beta: 2}
	data := []float64{1, 2, 3, 4, 5}
	post := NormalInverseGammaPosterior(prior, data)
	// kappa_n = 6, mu_n = (0 + 15)/6 = 2.5
	if !approx(post.Kappa, 6, tol) {
		t.Errorf("nig kappa %v", post.Kappa)
	}
	if !approx(post.Mu, 2.5, tol) {
		t.Errorf("nig mu %v", post.Mu)
	}
	// alpha_n = 2 + 2.5 = 4.5
	if !approx(post.Alpha, 4.5, tol) {
		t.Errorf("nig alpha %v", post.Alpha)
	}
	// SS about mean = 10 -> 0.5*SS = 5
	// cross term = 0.5 * kappa0*n/kappaN * (xbar-mu0)^2 = 0.5*5/6*9 = 3.75
	// beta_n = 2 + 5 + 3.75 = 10.75
	if !approx(post.Beta, 10.75, tol) {
		t.Errorf("nig beta %v", post.Beta)
	}
	mm := post.MarginalMean()
	if !approx(mm.Nu, 9, tol) || !approx(mm.Loc, 2.5, tol) {
		t.Errorf("nig marginal mean %v", mm)
	}
	pv := post.Predictive()
	if pv.Scale <= mm.Scale {
		t.Errorf("predictive scale should exceed marginal-mean scale")
	}
}

func TestDirichletMultinomialUpdate(t *testing.T) {
	prior := Dirichlet{Alpha: []float64{1, 1, 1}}
	post := DirichletMultinomialPosterior(prior, []int{3, 5, 2})
	want := []float64{4, 6, 3}
	for i := range want {
		if !approx(post.Alpha[i], want[i], tol) {
			t.Errorf("dm posterior[%d] %v", i, post.Alpha[i])
		}
	}
	pred := DirichletMultinomialPredictive(post)
	if !approx(pred[1], 6.0/13.0, 1e-12) {
		t.Errorf("dm predictive %v", pred[1])
	}
}

func TestInverseGammaVarianceUpdate(t *testing.T) {
	prior := InverseGamma{Shape: 1, Scale: 1}
	post := InverseGammaVariancePosterior(prior, 0, []float64{1, -1, 2, -2})
	// n=4 -> shape 1+2=3 ; ss = 1+1+4+4=10 -> scale 1+5=6
	if !approx(post.Shape, 3, tol) || !approx(post.Scale, 6, tol) {
		t.Errorf("ig var posterior %v", post)
	}
}

// ------------------------------------------------------------------
// Credible intervals
// ------------------------------------------------------------------

func TestCredibleIntervals(t *testing.T) {
	d := Normal{Mu: 0, Sigma: 1}
	iv := NormalCredibleInterval(d, 0.95)
	if !approx(iv.Lower, -1.959963984540054, 1e-6) || !approx(iv.Upper, 1.959963984540054, 1e-6) {
		t.Errorf("normal ci %v", iv)
	}
	// HDI of symmetric normal ~ equal tailed
	hdi := HighestDensityInterval(d, 0.95)
	if math.Abs(hdi.Width()-iv.Width()) > 1e-2 {
		t.Errorf("normal hdi width %v vs %v", hdi.Width(), iv.Width())
	}
	// coverage of a beta interval
	b := Beta{Alpha: 5, Beta: 3}
	bi := BetaCredibleInterval(b, 0.9)
	cov := b.CDF(bi.Upper) - b.CDF(bi.Lower)
	if !approx(cov, 0.9, 1e-6) {
		t.Errorf("beta ci coverage %v", cov)
	}
}

func TestProbabilityBetaExceeds(t *testing.T) {
	// Identical betas => P(X>Y)=0.5
	a := Beta{Alpha: 3, Beta: 4}
	p := ProbabilityBetaExceedsBeta(a, a)
	if !approx(p, 0.5, 1e-4) {
		t.Errorf("P(X>Y) identical %v", p)
	}
	// Strongly separated
	hi := Beta{Alpha: 30, Beta: 3}
	lo := Beta{Alpha: 3, Beta: 30}
	p2 := ProbabilityBetaExceedsBeta(hi, lo)
	if p2 < 0.999 {
		t.Errorf("P(X>Y) separated %v", p2)
	}
}

// ------------------------------------------------------------------
// Evidence and Bayes factors
// ------------------------------------------------------------------

func TestMarginalLikelihoodBetaBinomial(t *testing.T) {
	// Uniform prior: marginal likelihood of k in n is 1/(n+1) for all k.
	lm := LogMarginalLikelihoodBetaBinomial(UniformBetaPrior(), 3, 10)
	if !approx(math.Exp(lm), 1.0/11.0, 1e-10) {
		t.Errorf("ml betabinom %v", math.Exp(lm))
	}
}

func TestMarginalLikelihoodDirichlet(t *testing.T) {
	prior := Dirichlet{Alpha: []float64{1, 1}}
	// Two categories, counts (a,b): equals 1/(N+1) * ... check normalization by
	// summing over all count vectors with N=3.
	var s float64
	N := 3
	for a := 0; a <= N; a++ {
		b := N - a
		lm := LogMarginalLikelihoodDirichletMultinomial(prior, []int{a, b})
		s += math.Exp(lm)
	}
	if !approx(s, 1, 1e-9) {
		t.Errorf("dirichlet ml normalization %v", s)
	}
}

func TestMarginalLikelihoodNormalKnownVar(t *testing.T) {
	prior := Normal{Mu: 0, Sigma: 2}
	data := []float64{1.0, -0.5, 0.3}
	lm := LogMarginalLikelihoodNormalKnownVariance(prior, 1, data)
	// Compare against direct multivariate-normal evidence:
	// D ~ N(0, sigma2 I + sigma0^2 * 11^T)
	sigma2 := 1.0
	s0 := 4.0
	n := len(data)
	// covariance matrix
	cov := make([][]float64, n)
	for i := range cov {
		cov[i] = make([]float64, n)
		for j := range cov[i] {
			cov[i][j] = s0
			if i == j {
				cov[i][j] += sigma2
			}
		}
	}
	logdet, quad := gaussianLogDetQuad(cov, data)
	want := -0.5*float64(n)*math.Log(2*math.Pi) - 0.5*logdet - 0.5*quad
	if !approx(lm, want, 1e-8) {
		t.Errorf("normal ml %v want %v", lm, want)
	}
}

// gaussianLogDetQuad computes log|C| and x^T C^{-1} x via Cholesky for testing.
func gaussianLogDetQuad(c [][]float64, x []float64) (float64, float64) {
	n := len(c)
	l := make([][]float64, n)
	for i := range l {
		l[i] = make([]float64, n)
	}
	for i := 0; i < n; i++ {
		for j := 0; j <= i; j++ {
			sum := c[i][j]
			for k := 0; k < j; k++ {
				sum -= l[i][k] * l[j][k]
			}
			if i == j {
				l[i][j] = math.Sqrt(sum)
			} else {
				l[i][j] = sum / l[j][j]
			}
		}
	}
	var logdet float64
	for i := 0; i < n; i++ {
		logdet += 2 * math.Log(l[i][i])
	}
	// solve L y = x
	y := make([]float64, n)
	for i := 0; i < n; i++ {
		sum := x[i]
		for k := 0; k < i; k++ {
			sum -= l[i][k] * y[k]
		}
		y[i] = sum / l[i][i]
	}
	var quad float64
	for i := 0; i < n; i++ {
		quad += y[i] * y[i]
	}
	return logdet, quad
}

func TestBayesFactor(t *testing.T) {
	bf := BayesFactor(math.Log(4), math.Log(2))
	if !approx(bf, 2, 1e-12) {
		t.Errorf("bf %v", bf)
	}
	if BayesFactorInterpretation(2) != "negligible" {
		t.Errorf("interp 2")
	}
	if BayesFactorInterpretation(50) != "strong" {
		t.Errorf("interp 50")
	}
	if BayesFactorInterpretation(1.0/200.0) != "decisive" {
		t.Errorf("interp small")
	}
	probs := PosteriorModelProbabilities([]float64{math.Log(2), math.Log(1)}, []float64{1, 1})
	if !approx(probs[0], 2.0/3.0, 1e-12) || !approx(probs[1], 1.0/3.0, 1e-12) {
		t.Errorf("model probs %v", probs)
	}
}

// ------------------------------------------------------------------
// Naive Bayes
// ------------------------------------------------------------------

func TestGaussianNB(t *testing.T) {
	X := [][]float64{
		{1, 1}, {1.2, 0.9}, {0.9, 1.1},
		{5, 5}, {5.2, 4.9}, {4.8, 5.1},
	}
	y := []int{0, 0, 0, 1, 1, 1}
	clf := NewGaussianNB()
	if err := clf.Fit(X, y); err != nil {
		t.Fatal(err)
	}
	pred, err := clf.Predict([]float64{1, 1})
	if err != nil || pred != 0 {
		t.Errorf("gnb predict got %v err %v", pred, err)
	}
	pred2, _ := clf.Predict([]float64{5, 5})
	if pred2 != 1 {
		t.Errorf("gnb predict2 %v", pred2)
	}
	probs, _ := clf.PredictProba([]float64{1, 1})
	if probs[0] <= probs[1] {
		t.Errorf("gnb proba %v", probs)
	}
	var s float64
	for _, p := range probs {
		s += p
	}
	if !approx(s, 1, 1e-10) {
		t.Errorf("gnb proba sum %v", s)
	}
	acc, _ := Accuracy(clf, X, y)
	if !approx(acc, 1, tol) {
		t.Errorf("gnb accuracy %v", acc)
	}
}

func TestMultinomialNB(t *testing.T) {
	// two "topics" over 3 words
	X := [][]float64{
		{3, 0, 0}, {2, 1, 0}, {0, 0, 3}, {0, 1, 2},
	}
	y := []int{0, 0, 1, 1}
	clf := NewMultinomialNB(1)
	if err := clf.Fit(X, y); err != nil {
		t.Fatal(err)
	}
	p, _ := clf.Predict([]float64{4, 0, 0})
	if p != 0 {
		t.Errorf("mnb predict %v", p)
	}
	p2, _ := clf.Predict([]float64{0, 0, 5})
	if p2 != 1 {
		t.Errorf("mnb predict2 %v", p2)
	}
	probs, _ := clf.PredictProba([]float64{4, 0, 0})
	var s float64
	for _, v := range probs {
		s += v
	}
	if !approx(s, 1, 1e-10) {
		t.Errorf("mnb proba sum %v", s)
	}
}

func TestBernoulliNB(t *testing.T) {
	X := [][]float64{
		{1, 0, 0}, {1, 1, 0}, {0, 0, 1}, {0, 1, 1},
	}
	y := []int{0, 0, 1, 1}
	clf := NewBernoulliNB(1, 0.5)
	if err := clf.Fit(X, y); err != nil {
		t.Fatal(err)
	}
	p, _ := clf.Predict([]float64{1, 0, 0})
	if p != 0 {
		t.Errorf("bnb predict %v", p)
	}
	p2, _ := clf.Predict([]float64{0, 0, 1})
	if p2 != 1 {
		t.Errorf("bnb predict2 %v", p2)
	}
}

func TestNotFitted(t *testing.T) {
	clf := NewGaussianNB()
	if _, err := clf.Predict([]float64{1}); err != ErrNotFitted {
		t.Errorf("expected ErrNotFitted got %v", err)
	}
}

// ------------------------------------------------------------------
// Factors and Bayesian network
// ------------------------------------------------------------------

func TestFactorOps(t *testing.T) {
	// f(A,B) with A,B binary
	f, err := NewFactor([]string{"A", "B"}, []int{2, 2}, []float64{0.1, 0.2, 0.3, 0.4})
	if err != nil {
		t.Fatal(err)
	}
	if !approx(f.Get(map[string]int{"A": 1, "B": 0}), 0.3, tol) {
		t.Errorf("factor get %v", f.Get(map[string]int{"A": 1, "B": 0}))
	}
	m := f.Marginalize("B")
	if !approx(m.Get(map[string]int{"A": 0}), 0.3, tol) || !approx(m.Get(map[string]int{"A": 1}), 0.7, tol) {
		t.Errorf("marginalize %v", m.Table)
	}
	r := f.Reduce(map[string]int{"A": 1})
	if !approx(r.Get(map[string]int{"B": 0}), 0.3, tol) || !approx(r.Get(map[string]int{"B": 1}), 0.4, tol) {
		t.Errorf("reduce %v", r.Table)
	}
	n := f.Normalize()
	if !approx(n.Sum(), 1, tol) {
		t.Errorf("normalize sum %v", n.Sum())
	}
}

func TestFactorMultiply(t *testing.T) {
	fa, _ := NewFactor([]string{"A"}, []int{2}, []float64{0.6, 0.4})
	fb, _ := NewFactor([]string{"A", "B"}, []int{2, 2}, []float64{0.5, 0.5, 0.1, 0.9})
	prod := fa.Multiply(fb)
	// P(A=0,B=1) = 0.6*0.5=0.3
	if !approx(prod.Get(map[string]int{"A": 0, "B": 1}), 0.3, tol) {
		t.Errorf("multiply %v", prod.Get(map[string]int{"A": 0, "B": 1}))
	}
	// P(A=1,B=1) = 0.4*0.9=0.36
	if !approx(prod.Get(map[string]int{"A": 1, "B": 1}), 0.36, tol) {
		t.Errorf("multiply2 %v", prod.Get(map[string]int{"A": 1, "B": 1}))
	}
}

// buildSprinkler constructs the classic Rain/Sprinkler/WetGrass network.
func buildSprinkler(t *testing.T) *BayesianNetwork {
	bn := NewBayesianNetwork()
	must := func(err error) {
		if err != nil {
			t.Fatal(err)
		}
	}
	must(bn.AddVariable("Rain", 2))
	must(bn.AddVariable("Sprinkler", 2))
	must(bn.AddVariable("Wet", 2))
	// P(Rain): 0->false, 1->true
	must(bn.SetCPT("Rain", nil, []float64{0.8, 0.2}))
	// P(Sprinkler | Rain): rows Rain=0,1 ; cols Sprinkler=0,1
	must(bn.SetCPT("Sprinkler", []string{"Rain"}, []float64{
		0.6, 0.4, // Rain=0
		0.99, 0.01, // Rain=1
	}))
	// P(Wet | Sprinkler, Rain): order parents [Sprinkler, Rain]
	must(bn.SetCPT("Wet", []string{"Sprinkler", "Rain"}, []float64{
		1.0, 0.0, // S=0,R=0
		0.2, 0.8, // S=0,R=1
		0.1, 0.9, // S=1,R=0
		0.01, 0.99, // S=1,R=1
	}))
	return bn
}

func TestBayesianNetwork(t *testing.T) {
	bn := buildSprinkler(t)
	// Joint should sum to 1 over all 8 assignments.
	var total float64
	for r := 0; r < 2; r++ {
		for s := 0; s < 2; s++ {
			for w := 0; w < 2; w++ {
				total += bn.JointProbability(map[string]int{"Rain": r, "Sprinkler": s, "Wet": w})
			}
		}
	}
	if !approx(total, 1, 1e-12) {
		t.Errorf("joint sum %v", total)
	}
	// Marginal P(Rain=1)=0.2
	pr, _ := bn.MarginalProbability("Rain", 1)
	if !approx(pr, 0.2, 1e-12) {
		t.Errorf("marginal rain %v", pr)
	}
	// Marginal of Wet: compute by hand
	// P(Wet=1) = sum over r,s P(r)P(s|r)P(Wet=1|s,r)
	pw, _ := bn.MarginalProbability("Wet", 1)
	want := 0.8*0.6*0.0 + 0.8*0.4*0.9 + 0.2*0.99*0.8 + 0.2*0.01*0.99
	if !approx(pw, want, 1e-12) {
		t.Errorf("marginal wet %v want %v", pw, want)
	}
	// Conditional P(Rain=1 | Wet=1) via Bayes
	prw, _ := bn.ConditionalProbability("Rain", 1, map[string]int{"Wet": 1})
	// numerator P(Rain=1,Wet=1) = 0.2*(0.99*0.8+0.01*0.99)
	num := 0.2 * (0.99*0.8 + 0.01*0.99)
	wantC := num / want
	if !approx(prw, wantC, 1e-10) {
		t.Errorf("cond rain|wet %v want %v", prw, wantC)
	}
	// Evidence probability P(Wet=1) matches marginal
	ev, _ := bn.EvidenceProbability(map[string]int{"Wet": 1})
	if !approx(ev, want, 1e-12) {
		t.Errorf("evidence %v want %v", ev, want)
	}
}

// ------------------------------------------------------------------
// Runnable examples
// ------------------------------------------------------------------

func ExampleBetaBernoulliPosterior() {
	prior := UniformBetaPrior()
	post := BetaBernoulliPosterior(prior, 8, 2)
	fmt.Printf("posterior mean = %.3f\n", post.Mean())
	ci := BetaCredibleInterval(post, 0.95)
	fmt.Printf("95%% CI = [%.3f, %.3f]\n", ci.Lower, ci.Upper)
	// Output:
	// posterior mean = 0.750
	// 95% CI = [0.482, 0.940]
}

func ExampleBayesianNetwork() {
	bn := NewBayesianNetwork()
	bn.AddVariable("Rain", 2)
	bn.AddVariable("Wet", 2)
	bn.SetCPT("Rain", nil, []float64{0.8, 0.2})
	bn.SetCPT("Wet", []string{"Rain"}, []float64{
		0.9, 0.1, // Rain=0 -> mostly dry
		0.2, 0.8, // Rain=1 -> mostly wet
	})
	p, _ := bn.ConditionalProbability("Rain", 1, map[string]int{"Wet": 1})
	fmt.Printf("P(Rain=1 | Wet=1) = %.3f\n", p)
	// Output:
	// P(Rain=1 | Wet=1) = 0.667
}

func ExampleGammaPoissonPosterior() {
	prior := Gamma{Shape: 2, Rate: 1}
	post := GammaPoissonPosteriorData(prior, []int{3, 4, 2, 5})
	fmt.Printf("rate posterior mean = %.4f\n", post.Mean())
	// Output:
	// rate posterior mean = 3.2000
}
