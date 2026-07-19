package infogeom

import (
	"errors"
	"fmt"
	"math"
	"testing"
)

// approx reports whether a and b agree to within the absolute tolerance tol.
func approx(a, b, tol float64) bool {
	if math.IsInf(a, 0) || math.IsInf(b, 0) {
		return a == b
	}
	return math.Abs(a-b) <= tol
}

func approxVec(a, b []float64, tol float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !approx(a[i], b[i], tol) {
			return false
		}
	}
	return true
}

func approxMat(a, b [][]float64, tol float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !approxVec(a[i], b[i], tol) {
			return false
		}
	}
	return true
}

func TestKLDivergence(t *testing.T) {
	tests := []struct {
		name string
		p, q []float64
		want float64
	}{
		{"identical", []float64{0.5, 0.5}, []float64{0.5, 0.5}, 0},
		{"half-quarter", []float64{0.5, 0.5}, []float64{0.25, 0.75}, 0.14384103622589042},
		{"skewed", []float64{0.7, 0.3}, []float64{0.5, 0.5}, 0.08228287850505178},
	}
	for _, tt := range tests {
		got, err := KLDivergence(tt.p, tt.q)
		if err != nil {
			t.Fatalf("%s: unexpected error %v", tt.name, err)
		}
		if !approx(got, tt.want, 1e-12) {
			t.Errorf("%s: KLDivergence = %v, want %v", tt.name, got, tt.want)
		}
	}
	// Absolute continuity failure yields +Inf.
	got, err := KLDivergence([]float64{0.5, 0.5}, []float64{1, 0})
	if err != nil {
		t.Fatalf("unexpected error %v", err)
	}
	if !math.IsInf(got, 1) {
		t.Errorf("expected +Inf, got %v", got)
	}
	// Malformed input.
	if _, err := KLDivergence([]float64{0.5, 0.6}, []float64{0.5, 0.5}); !errors.Is(err, ErrNotProb) {
		t.Errorf("expected ErrNotProb, got %v", err)
	}
	if _, err := KLDivergence([]float64{1}, []float64{0.5, 0.5}); !errors.Is(err, ErrDim) {
		t.Errorf("expected ErrDim, got %v", err)
	}
}

func TestEntropyAndCrossEntropy(t *testing.T) {
	p := []float64{0.5, 0.5}
	h, err := Entropy(p)
	if err != nil || !approx(h, ln2, 1e-12) {
		t.Errorf("Entropy = %v (err %v), want ln2", h, err)
	}
	if hb, _ := EntropyBits(p); !approx(hb, 1, 1e-12) {
		t.Errorf("EntropyBits = %v, want 1", hb)
	}
	// Cross entropy H(p,p) equals entropy.
	ce, _ := CrossEntropy(p, p)
	if !approx(ce, h, 1e-12) {
		t.Errorf("CrossEntropy(p,p) = %v, want %v", ce, h)
	}
	// KL = cross entropy - entropy.
	q := []float64{0.25, 0.75}
	ce2, _ := CrossEntropy(p, q)
	kl, _ := KLDivergence(p, q)
	if !approx(ce2-h, kl, 1e-12) {
		t.Errorf("CrossEntropy-Entropy = %v, want KL %v", ce2-h, kl)
	}
}

func TestSymmetricDivergences(t *testing.T) {
	p := []float64{0.5, 0.5}
	q := []float64{0.25, 0.75}
	if jsd, _ := JensenShannonDivergence(p, q); jsd < 0 || jsd > ln2 {
		t.Errorf("JSD out of range: %v", jsd)
	}
	// JSD of disjoint distributions is ln2.
	jsd, _ := JensenShannonDivergence([]float64{1, 0}, []float64{0, 1})
	if !approx(jsd, ln2, 1e-12) {
		t.Errorf("JSD disjoint = %v, want ln2", jsd)
	}
	if jb, _ := JensenShannonDivergenceBits([]float64{1, 0}, []float64{0, 1}); !approx(jb, 1, 1e-12) {
		t.Errorf("JSD bits disjoint = %v, want 1", jb)
	}
	// Jeffreys is symmetric.
	j1, _ := JeffreysDivergence(p, q)
	j2, _ := JeffreysDivergence(q, p)
	if !approx(j1, j2, 1e-12) {
		t.Errorf("Jeffreys not symmetric: %v vs %v", j1, j2)
	}
	// JS distance is the sqrt of the divergence.
	d, _ := JensenShannonDistance(p, q)
	if !approx(d*d, jsdOf(t, p, q), 1e-12) {
		t.Errorf("JS distance squared mismatch")
	}
}

func jsdOf(t *testing.T, p, q []float64) float64 {
	t.Helper()
	v, err := JensenShannonDivergence(p, q)
	if err != nil {
		t.Fatal(err)
	}
	return v
}

func TestOtherDivergences(t *testing.T) {
	p := []float64{0.5, 0.5}
	q := []float64{0.25, 0.75}
	if v, _ := RenyiDivergence(p, q, 2); !approx(v, 0.28768207245178085, 1e-12) {
		t.Errorf("RenyiDivergence(2) = %v", v)
	}
	// Renyi order 1 equals KL.
	kl, _ := KLDivergence(p, q)
	if v, _ := RenyiDivergence(p, q, 1); !approx(v, kl, 1e-12) {
		t.Errorf("RenyiDivergence(1) = %v, want KL %v", v, kl)
	}
	if v, _ := ChiSquaredDivergence(p, q); !approx(v, 1.0/3.0, 1e-12) {
		t.Errorf("ChiSquared = %v", v)
	}
	if v, _ := TotalVariationDistance(p, q); !approx(v, 0.25, 1e-12) {
		t.Errorf("TV = %v", v)
	}
	if v, _ := BhattacharyyaCoefficient(p, q); !approx(v, 0.9659258262890682, 1e-12) {
		t.Errorf("BC = %v", v)
	}
	if v, _ := HellingerDistance(p, q); !approx(v, 0.18459191128251476, 1e-12) {
		t.Errorf("Hellinger = %v", v)
	}
	// Alpha-divergence at alpha=0 equals 4*(1-BC).
	if v, _ := AlphaDivergence(p, q, 0); !approx(v, 0.1362966948437272, 1e-12) {
		t.Errorf("AlphaDivergence(0) = %v", v)
	}
	// Alpha-divergence limit alpha->1 equals KL.
	if v, _ := AlphaDivergence(p, q, 1); !approx(v, kl, 1e-12) {
		t.Errorf("AlphaDivergence(1) = %v, want KL %v", v, kl)
	}
	// Tsallis limit equals KL.
	if v, _ := TsallisDivergence(p, q, 1); !approx(v, kl, 1e-12) {
		t.Errorf("TsallisDivergence(1) = %v", v)
	}
}

func TestFDivergenceMatchesKL(t *testing.T) {
	p := []float64{0.5, 0.5}
	q := []float64{0.25, 0.75}
	// generator f(t) = t ln t reproduces the KL divergence.
	f := func(x float64) float64 {
		if x <= 0 {
			return 0
		}
		return x * math.Log(x)
	}
	got, _ := FDivergence(p, q, f)
	kl, _ := KLDivergence(p, q)
	if !approx(got, kl, 1e-12) {
		t.Errorf("FDivergence = %v, want KL %v", got, kl)
	}
}

func TestBregman(t *testing.T) {
	if v, _ := SquaredEuclideanDivergence([]float64{1, 2}, []float64{0, 0}); !approx(v, 2.5, 1e-12) {
		t.Errorf("SquaredEuclidean = %v", v)
	}
	// negative-entropy Bregman reproduces generalised KL, which on the simplex
	// equals ordinary KL.
	p := []float64{0.5, 0.5}
	q := []float64{0.25, 0.75}
	gen := NegativeEntropyGenerator()
	bd, _ := BregmanDivergence(gen, p, q)
	kl, _ := KLDivergence(p, q)
	if !approx(bd, kl, 1e-10) {
		t.Errorf("Bregman(negentropy) = %v, want KL %v", bd, kl)
	}
	gkl, _ := GeneralizedKLDivergence(p, q)
	if !approx(gkl, kl, 1e-12) {
		t.Errorf("GeneralizedKL = %v, want KL %v", gkl, kl)
	}
	if v, _ := ItakuraSaitoDivergence([]float64{2}, []float64{1}); !approx(v, 0.3068528194400546, 1e-12) {
		t.Errorf("ItakuraSaito = %v", v)
	}
	// centroid is the arithmetic mean.
	c, _ := BregmanCentroid([][]float64{{0, 0}, {2, 4}})
	if !approxVec(c, []float64{1, 2}, 1e-12) {
		t.Errorf("BregmanCentroid = %v", c)
	}
}

func TestFisherInformation(t *testing.T) {
	if v, _ := FisherInformationBernoulli(0.5); !approx(v, 4, 1e-12) {
		t.Errorf("FisherInformationBernoulli(0.5) = %v", v)
	}
	if v, _ := FisherInformationPoisson(2); !approx(v, 0.5, 1e-12) {
		t.Errorf("FisherInformationPoisson(2) = %v", v)
	}
	if v, _ := FisherInformationExponential(2); !approx(v, 0.25, 1e-12) {
		t.Errorf("FisherInformationExponential(2) = %v", v)
	}
	// Gaussian Fisher matrix in (mu, sigma).
	g := Gaussian{Mu: 0, Sigma: 2}
	fm, _ := g.FisherInformationMuSigma()
	want := [][]float64{{0.25, 0}, {0, 0.5}}
	if !approxMat(fm, want, 1e-12) {
		t.Errorf("FisherInformationMuSigma = %v", fm)
	}
	// Categorical reduced Fisher matrix.
	fc, _ := FisherInformationCategorical([]float64{0.2, 0.3, 0.5})
	wantc := [][]float64{{7, 2}, {2, 1.0/0.3 + 2}}
	if !approxMat(fc, wantc, 1e-12) {
		t.Errorf("FisherInformationCategorical = %v", fc)
	}
	if _, err := FisherInformationBernoulli(1); !errors.Is(err, ErrDomain) {
		t.Errorf("expected ErrDomain for p=1, got %v", err)
	}
}

func TestFisherRao(t *testing.T) {
	if d, _ := FisherRaoCategorical([]float64{1, 0}, []float64{0, 1}); !approx(d, math.Pi, 1e-12) {
		t.Errorf("FisherRaoCategorical disjoint = %v, want pi", d)
	}
	if d, _ := FisherRaoCategorical([]float64{0.5, 0.5}, []float64{0.25, 0.75}); !approx(d, 0.5235987755982995, 1e-12) {
		t.Errorf("FisherRaoCategorical = %v", d)
	}
	if d, _ := FisherRaoBernoulli(0, 1); !approx(d, math.Pi, 1e-12) {
		t.Errorf("FisherRaoBernoulli(0,1) = %v, want pi", d)
	}
	if d, _ := FisherRaoBernoulli(0.5, 0.5); !approx(d, 0, 1e-12) {
		t.Errorf("FisherRaoBernoulli identical = %v", d)
	}
	// Gaussian geodesic.
	d, _ := FisherRaoGaussian(Gaussian{0, 1}, Gaussian{1, 1})
	if !approx(d, 0.9802581434685472, 1e-12) {
		t.Errorf("FisherRaoGaussian = %v", d)
	}
	if d, _ := FisherRaoGaussianFixedVariance(0, 1, 1); !approx(d, 1, 1e-12) {
		t.Errorf("FisherRaoGaussianFixedVariance = %v", d)
	}
	if d, _ := FisherRaoGaussianFixedMean(1, math.E); !approx(d, math.Sqrt2, 1e-12) {
		t.Errorf("FisherRaoGaussianFixedMean = %v", d)
	}
	if d, _ := FisherRaoPoisson(1, 4); !approx(d, 2, 1e-12) {
		t.Errorf("FisherRaoPoisson = %v", d)
	}
	if d, _ := FisherRaoExponential(1, math.E); !approx(d, 1, 1e-12) {
		t.Errorf("FisherRaoExponential = %v", d)
	}
	// Scalar geodesic distance for the Poisson metric 1/lambda between 1 and 4
	// should reproduce the closed form 2.
	sd, err := FisherRaoScalar(func(x float64) float64 { return 1 / x }, 1, 4, 200)
	if err != nil || !approx(sd, 2, 1e-6) {
		t.Errorf("FisherRaoScalar = %v (err %v), want 2", sd, err)
	}
}

func TestGaussianKL(t *testing.T) {
	v, err := KLDivergenceGaussian(Gaussian{0, 1}, Gaussian{1, 2})
	if err != nil || !approx(v, 0.4431471805599453, 1e-12) {
		t.Errorf("KLDivergenceGaussian = %v (err %v)", v, err)
	}
	if v, _ := KLDivergenceGaussian(Gaussian{3, 2}, Gaussian{3, 2}); !approx(v, 0, 1e-12) {
		t.Errorf("KLDivergenceGaussian identical = %v", v)
	}
}

func TestExponentialFamily(t *testing.T) {
	f := BernoulliFamily()
	// mean parameter at theta=0 is 0.5.
	if m := f.MeanParameters([]float64{0}); !approxVec(m, []float64{0.5}, 1e-12) {
		t.Errorf("MeanParameters = %v", m)
	}
	// Fisher information at theta=0 is 0.25.
	if fi := f.FisherInformationNatural([]float64{0}); !approxMat(fi, [][]float64{{0.25}}, 1e-12) {
		t.Errorf("FisherInformationNatural = %v", fi)
	}
	// KL between two Bernoulli members via the Bregman form.
	theta2, _ := Bernoulli{P: Sigmoid(1)}.NaturalParameter()
	kl, _ := KLDivergenceExpFamily(f, []float64{0}, []float64{theta2})
	if !approx(kl, 0.12011450695827752, 1e-9) {
		t.Errorf("KLDivergenceExpFamily = %v", kl)
	}
	// Numerical mean parameter matches analytic for the Poisson family.
	pf := PoissonFamily()
	num := ExponentialFamily{LogPartition: pf.LogPartition}
	if !approxVec(num.MeanParameters([]float64{0.5}), pf.MeanParameters([]float64{0.5}), 1e-5) {
		t.Errorf("numerical mean mismatch")
	}
}

func TestNaturalFromExpectation(t *testing.T) {
	f := BernoulliFamily()
	theta, err := NaturalFromExpectation(f, []float64{0.5}, []float64{0.4}, 1e-12, 50)
	if err != nil || !approxVec(theta, []float64{0}, 1e-9) {
		t.Errorf("NaturalFromExpectation = %v (err %v)", theta, err)
	}
	// Recover a non-trivial mean.
	theta, err = NaturalFromExpectation(f, []float64{Sigmoid(1)}, []float64{0}, 1e-12, 50)
	if err != nil || !approxVec(theta, []float64{1}, 1e-8) {
		t.Errorf("NaturalFromExpectation = %v (err %v), want [1]", theta, err)
	}
}

func TestNaturalGradient(t *testing.T) {
	fisher := [][]float64{{2, 0}, {0, 4}}
	grad := []float64{2, 4}
	ng, err := NaturalGradient(fisher, grad)
	if err != nil || !approxVec(ng, []float64{1, 1}, 1e-12) {
		t.Errorf("NaturalGradient = %v (err %v)", ng, err)
	}
	step, _ := NaturalGradientStep([]float64{0, 0}, grad, fisher, 0.5)
	if !approxVec(step, []float64{-0.5, -0.5}, 1e-12) {
		t.Errorf("NaturalGradientStep = %v", step)
	}
	// dual norm sqrt(g^T F^-1 g) = sqrt(2*1 + 4*1) = sqrt 6.
	nn, _ := NaturalGradientNorm(fisher, grad)
	if !approx(nn, math.Sqrt(6), 1e-12) {
		t.Errorf("NaturalGradientNorm = %v", nn)
	}
	// damped gradient with large lambda shrinks toward zero.
	d, _ := DampedNaturalGradient(fisher, grad, 1e6)
	if Norm2(d) > 1e-3 {
		t.Errorf("DampedNaturalGradient too large: %v", d)
	}
}

func TestMirrorDescent(t *testing.T) {
	// zero gradient leaves the distribution unchanged.
	p, err := MirrorDescentStep([]float64{0.5, 0.5}, []float64{0, 0}, 1)
	if err != nil || !approxVec(p, []float64{0.5, 0.5}, 1e-12) {
		t.Errorf("MirrorDescentStep = %v (err %v)", p, err)
	}
	// output is a valid distribution.
	p2, _ := MirrorDescentStep([]float64{0.3, 0.7}, []float64{1, -1}, 0.5)
	if !IsProbabilityVector(p2, 1e-12) {
		t.Errorf("MirrorDescentStep not a distribution: %v", p2)
	}
}

func TestGeodesics(t *testing.T) {
	p := []float64{0.2, 0.8}
	q := []float64{0.6, 0.4}
	// endpoints.
	if m, _ := MixtureGeodesic(p, q, 0); !approxVec(m, p, 1e-12) {
		t.Errorf("MixtureGeodesic t=0 = %v", m)
	}
	if m, _ := MixtureGeodesic(p, q, 1); !approxVec(m, q, 1e-12) {
		t.Errorf("MixtureGeodesic t=1 = %v", m)
	}
	// midpoint of the mixture geodesic is the average.
	if m, _ := MixtureGeodesic(p, q, 0.5); !approxVec(m, []float64{0.4, 0.6}, 1e-12) {
		t.Errorf("MixtureGeodesic midpoint = %v", m)
	}
	// exponential geodesic endpoints.
	if e, _ := ExponentialGeodesic(p, q, 0); !approxVec(e, p, 1e-12) {
		t.Errorf("ExponentialGeodesic t=0 = %v", e)
	}
	if e, _ := ExponentialGeodesic(p, q, 1); !approxVec(e, q, 1e-12) {
		t.Errorf("ExponentialGeodesic t=1 = %v", e)
	}
	// alpha-mixture with alpha=1 equals exponential geodesic.
	am, _ := AlphaMixture(p, q, 0.5, 1)
	eg, _ := ExponentialGeodesic(p, q, 0.5)
	if !approxVec(am, eg, 1e-12) {
		t.Errorf("AlphaMixture(alpha=1) = %v, want %v", am, eg)
	}
	// alpha-mixture with alpha=-1 equals mixture geodesic.
	am2, _ := AlphaMixture(p, q, 0.5, -1)
	mg, _ := MixtureGeodesic(p, q, 0.5)
	if !approxVec(am2, mg, 1e-12) {
		t.Errorf("AlphaMixture(alpha=-1) = %v, want %v", am2, mg)
	}
}

func TestConnectionTensors(t *testing.T) {
	pf := PoissonFamily()
	// third cumulant of Poisson at theta=0 is e^0 = 1.
	tt := CumulantTensor3(pf.HessLogPartition, []float64{0}, 1e-4)
	if !approx(tt[0][0][0], 1, 1e-4) {
		t.Errorf("CumulantTensor3 = %v, want 1", tt[0][0][0])
	}
	// e-connection (alpha=1) coefficients vanish.
	e := AlphaConnectionCoefficients(pf.HessLogPartition, []float64{0}, 1, 1e-4)
	if !approx(e[0][0][0], 0, 1e-6) {
		t.Errorf("e-connection coefficient = %v, want 0", e[0][0][0])
	}
	// m-connection (alpha=-1) coefficients equal the full skewness tensor.
	m := AlphaConnectionCoefficients(pf.HessLogPartition, []float64{0}, -1, 1e-4)
	if !approx(m[0][0][0], 1, 1e-4) {
		t.Errorf("m-connection coefficient = %v, want 1", m[0][0][0])
	}
	if DualAlpha(0.3) != -0.3 {
		t.Errorf("DualAlpha wrong")
	}
}

func TestEscortAndEmbedding(t *testing.T) {
	// escort at alpha=1 is the identity.
	e, _ := EscortDistribution([]float64{0.2, 0.3, 0.5}, 1)
	if !approxVec(e, []float64{0.2, 0.3, 0.5}, 1e-12) {
		t.Errorf("EscortDistribution(1) = %v", e)
	}
	// large alpha concentrates on the mode.
	e2, _ := EscortDistribution([]float64{0.2, 0.3, 0.5}, 20)
	if e2[2] < 0.99 {
		t.Errorf("EscortDistribution(20) not concentrated: %v", e2)
	}
	// alpha-embedding at alpha=1 is the logarithm.
	emb, _ := AlphaEmbedding([]float64{0.5, 0.5}, 1)
	if !approxVec(emb, []float64{math.Log(0.5), math.Log(0.5)}, 1e-12) {
		t.Errorf("AlphaEmbedding(1) = %v", emb)
	}
}

func TestLinearAlgebra(t *testing.T) {
	a := [][]float64{{4, 0}, {0, 2}}
	inv, err := Inverse(a)
	if err != nil || !approxMat(inv, [][]float64{{0.25, 0}, {0, 0.5}}, 1e-12) {
		t.Errorf("Inverse = %v (err %v)", inv, err)
	}
	if d, _ := Determinant([][]float64{{1, 2}, {3, 4}}); !approx(d, -2, 1e-12) {
		t.Errorf("Determinant = %v", d)
	}
	x, err := Solve([][]float64{{2, 1}, {1, 3}}, []float64{3, 4})
	if err != nil || !approxVec(x, []float64{1, 1}, 1e-12) {
		t.Errorf("Solve = %v (err %v)", x, err)
	}
	// Cholesky reconstruction.
	spd := [][]float64{{4, 2}, {2, 3}}
	l, err := Cholesky(spd)
	if err != nil {
		t.Fatal(err)
	}
	lt, _ := Transpose(l)
	recon, _ := MatMul(l, lt)
	if !approxMat(recon, spd, 1e-12) {
		t.Errorf("Cholesky reconstruction = %v", recon)
	}
	if !IsPositiveDefinite(spd) {
		t.Errorf("expected positive definite")
	}
	ld, _ := LogDet(spd)
	det, _ := Determinant(spd)
	if !approx(ld, math.Log(det), 1e-10) {
		t.Errorf("LogDet = %v, want %v", ld, math.Log(det))
	}
	if _, err := Inverse([][]float64{{1, 1}, {1, 1}}); !errors.Is(err, ErrSingular) {
		t.Errorf("expected ErrSingular, got %v", err)
	}
}

func TestSoftmaxAndSigmoid(t *testing.T) {
	if !approx(Sigmoid(0), 0.5, 1e-12) {
		t.Errorf("Sigmoid(0) = %v", Sigmoid(0))
	}
	if !approx(Logit(Sigmoid(1.3)), 1.3, 1e-12) {
		t.Errorf("Logit/Sigmoid roundtrip failed")
	}
	s := Softmax([]float64{0, 0, 0})
	if !approxVec(s, []float64{1.0 / 3, 1.0 / 3, 1.0 / 3}, 1e-12) {
		t.Errorf("Softmax uniform = %v", s)
	}
	if !approx(LogSumExp([]float64{0, 0}), ln2, 1e-12) {
		t.Errorf("LogSumExp = %v", LogSumExp([]float64{0, 0}))
	}
	// numerical stability with large inputs.
	big := Softmax([]float64{1000, 1000})
	if !approxVec(big, []float64{0.5, 0.5}, 1e-12) {
		t.Errorf("Softmax overflow: %v", big)
	}
}

func TestNumericalDerivatives(t *testing.T) {
	f := func(x []float64) float64 { return x[0]*x[0] + 2*x[1]*x[1] }
	g := NumericalGradient(f, []float64{1, 1}, 1e-6)
	if !approxVec(g, []float64{2, 4}, 1e-4) {
		t.Errorf("NumericalGradient = %v", g)
	}
	h := NumericalHessian(f, []float64{1, 1}, 1e-3)
	if !approxMat(h, [][]float64{{2, 0}, {0, 4}}, 1e-3) {
		t.Errorf("NumericalHessian = %v", h)
	}
}

func TestMetricHelpers(t *testing.T) {
	g := [][]float64{{2, 0}, {0, 8}}
	n, _ := MetricNorm(g, []float64{1, 0})
	if !approx(n, math.Sqrt(2), 1e-12) {
		t.Errorf("MetricNorm = %v", n)
	}
	// orthogonal under the metric.
	ang, _ := MetricAngle(g, []float64{1, 0}, []float64{0, 1})
	if !approx(ang, math.Pi/2, 1e-12) {
		t.Errorf("MetricAngle = %v", ang)
	}
}

// ExampleKLDivergence demonstrates computing the Kullback-Leibler divergence
// between two discrete distributions and its symmetric Jensen-Shannon
// counterpart.
func ExampleKLDivergence() {
	p := []float64{0.5, 0.5}
	q := []float64{0.25, 0.75}
	kl, _ := KLDivergence(p, q)
	js, _ := JensenShannonDivergence(p, q)
	fmt.Printf("KL = %.4f nats\nJS = %.4f nats\n", kl, js)
	// Output:
	// KL = 0.1438 nats
	// JS = 0.0338 nats
}

// ExampleFisherRaoCategorical demonstrates the Fisher-Rao geodesic distance
// between two categorical distributions, which for antipodal (disjoint)
// distributions equals pi.
func ExampleFisherRaoCategorical() {
	d, _ := FisherRaoCategorical([]float64{1, 0}, []float64{0, 1})
	fmt.Printf("%.4f\n", d)
	// Output:
	// 3.1416
}
