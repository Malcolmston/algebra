package stats

import (
	"math"
	"testing"
)

const tol = 1e-9

func approx(t *testing.T, got, want, eps float64, name string) {
	t.Helper()
	if math.IsNaN(want) {
		if !math.IsNaN(got) {
			t.Errorf("%s = %v, want NaN", name, got)
		}
		return
	}
	if math.Abs(got-want) > eps {
		t.Errorf("%s = %v, want %v (diff %g)", name, got, want, math.Abs(got-want))
	}
}

func TestDescriptiveBasics(t *testing.T) {
	xs := []float64{2, 4, 4, 4, 5, 5, 7, 9}
	approx(t, Mean(xs), 5, tol, "Mean")
	approx(t, PopVariance(xs), 4, tol, "PopVariance")
	approx(t, PopStdDev(xs), 2, tol, "PopStdDev")
	approx(t, Variance(xs), 32.0/7, tol, "Variance")
	approx(t, StdDev(xs), math.Sqrt(32.0/7), tol, "StdDev")
	approx(t, Sum(xs), 40, tol, "Sum")
	approx(t, Min(xs), 2, tol, "Min")
	approx(t, Max(xs), 9, tol, "Max")
	approx(t, Range(xs), 7, tol, "Range")
	approx(t, Product([]float64{2, 3, 4}), 24, tol, "Product")
}

func TestMedian(t *testing.T) {
	approx(t, Median([]float64{3, 1, 2}), 2, tol, "Median odd")
	approx(t, Median([]float64{4, 2, 1, 3}), 2.5, tol, "Median even")
	approx(t, Median([]float64{5}), 5, tol, "Median single")
	// Input must not be mutated.
	in := []float64{3, 1, 2}
	_ = Median(in)
	if in[0] != 3 || in[1] != 1 || in[2] != 2 {
		t.Errorf("Median mutated its input: %v", in)
	}
}

func TestMode(t *testing.T) {
	got := Mode([]float64{1, 2, 2, 3, 3, 4})
	if len(got) != 2 || got[0] != 2 || got[1] != 3 {
		t.Errorf("Mode multimodal = %v, want [2 3]", got)
	}
	if Mode([]float64{1, 2, 3}) != nil {
		t.Errorf("Mode of all-unique should be nil")
	}
	single := Mode([]float64{4, 4, 1, 2})
	if len(single) != 1 || single[0] != 4 {
		t.Errorf("Mode = %v, want [4]", single)
	}
}

func TestQuantileAndIQR(t *testing.T) {
	xs := []float64{2, 4, 4, 4, 5, 5, 7, 9}
	approx(t, Quantile(xs, 0), 2, tol, "Quantile 0")
	approx(t, Quantile(xs, 1), 9, tol, "Quantile 1")
	approx(t, Quantile(xs, 0.5), Median(xs), tol, "Quantile 0.5 == Median")
	approx(t, Percentile(xs, 25), 4, tol, "Percentile 25")
	approx(t, IQR(xs), 1.5, tol, "IQR")
	if !math.IsNaN(Quantile(xs, 1.5)) {
		t.Errorf("Quantile out of range should be NaN")
	}
}

func TestMeansAndMoments(t *testing.T) {
	approx(t, GeometricMean([]float64{1, 2, 4}), 2, tol, "GeometricMean")
	approx(t, HarmonicMean([]float64{1, 2, 4}), 12.0/7, tol, "HarmonicMean")
	approx(t, WeightedMean([]float64{1, 2, 3}, []float64{1, 1, 2}), 2.25, tol, "WeightedMean")
	approx(t, ZScore(1.96, 0, 1), 1.96, tol, "ZScore")
	xs := []float64{2, 4, 4, 4, 5, 5, 7, 9}
	approx(t, Skewness(xs), 0.65625, tol, "Skewness")
	approx(t, Kurtosis(xs), -0.21875, tol, "Kurtosis")
	// Symmetric data has zero skew.
	approx(t, Skewness([]float64{1, 2, 3, 4, 5}), 0, tol, "Skewness symmetric")
}

func TestCovarianceCorrelation(t *testing.T) {
	a := []float64{1, 2, 3, 4, 5}
	b := []float64{2, 4, 6, 8, 10}
	approx(t, Correlation(a, b), 1, tol, "Correlation perfect positive")
	approx(t, Covariance(a, b), 5, tol, "Covariance")
	c := []float64{10, 8, 6, 4, 2}
	approx(t, Correlation(a, c), -1, tol, "Correlation perfect negative")
	if !math.IsNaN(Correlation(a, []float64{7, 7, 7, 7, 7})) {
		t.Errorf("Correlation with zero-variance series should be NaN")
	}
}

func TestNaNGuards(t *testing.T) {
	if !math.IsNaN(Mean(nil)) {
		t.Errorf("Mean(nil) should be NaN")
	}
	if !math.IsNaN(Variance([]float64{1})) {
		t.Errorf("Variance of single element should be NaN")
	}
	if !math.IsNaN(GeometricMean([]float64{1, -2})) {
		t.Errorf("GeometricMean with non-positive should be NaN")
	}
	if !math.IsNaN(Covariance([]float64{1, 2}, []float64{1})) {
		t.Errorf("Covariance mismatched lengths should be NaN")
	}
}

func TestCombinatorics(t *testing.T) {
	approx(t, Factorial(0), 1, tol, "Factorial 0")
	approx(t, Factorial(5), 120, tol, "Factorial 5")
	approx(t, Factorial(10), 3628800, tol, "Factorial 10")
	if !math.IsNaN(Factorial(-1)) {
		t.Errorf("Factorial(-1) should be NaN")
	}
	approx(t, Choose(5, 2), 10, tol, "Choose 5,2")
	approx(t, Choose(52, 5), 2598960, 0.5, "Choose 52,5")
	approx(t, Choose(5, 0), 1, tol, "Choose n,0")
	approx(t, Choose(5, 6), 0, tol, "Choose k>n")
	approx(t, Perm(5, 2), 20, tol, "Perm 5,2")
	approx(t, Perm(10, 3), 720, tol, "Perm 10,3")
	// Symmetry C(n,k) == C(n,n-k).
	approx(t, Choose(30, 12), Choose(30, 18), 1, "Choose symmetry")
}

func TestNormal(t *testing.T) {
	n := Normal{Mu: 0, Sigma: 1}
	approx(t, n.CDF(0), 0.5, tol, "NormalCDF(0)")
	approx(t, n.CDF(1.96), 0.9750021048, 1e-9, "NormalCDF(1.96)")
	approx(t, n.CDF(-1.96), 0.0249978952, 1e-9, "NormalCDF(-1.96)")
	approx(t, n.CDF(1), 0.8413447461, 1e-9, "NormalCDF(1)")
	approx(t, n.PDF(0), 1/math.Sqrt(2*math.Pi), tol, "NormalPDF(0)")
	approx(t, NormalCDF(0, 0, 1), 0.5, tol, "NormalCDF wrapper")
	approx(t, NormalPDF(0, 0, 1), 1/math.Sqrt(2*math.Pi), tol, "NormalPDF wrapper")
	approx(t, n.Mean(), 0, tol, "Normal Mean")
	approx(t, n.Variance(), 1, tol, "Normal Variance")
	// Quantile is the inverse of CDF.
	approx(t, n.Quantile(0.975), 1.959963985, 1e-7, "Normal Quantile 0.975")
	for _, p := range []float64{0.01, 0.1, 0.25, 0.5, 0.83, 0.99} {
		approx(t, n.CDF(n.Quantile(p)), p, 1e-9, "CDF(Quantile) roundtrip")
	}
	// Shifted/scaled.
	m := Normal{Mu: 10, Sigma: 2}
	approx(t, m.CDF(10), 0.5, tol, "shifted CDF at mean")
	approx(t, m.Quantile(0.5), 10, 1e-9, "shifted Quantile 0.5")
}

func TestBinomial(t *testing.T) {
	b := Binomial{N: 5, P: 0.5}
	// C(5,2)*0.5^5 = 10/32 = 0.3125
	approx(t, b.PMF(2), 0.3125, tol, "Binomial PMF(2;5,0.5)")
	approx(t, b.PMF(0), 1.0/32, tol, "Binomial PMF(0)")
	approx(t, b.PMF(5), 1.0/32, tol, "Binomial PMF(5)")
	approx(t, b.PMF(6), 0, tol, "Binomial PMF out of range")
	approx(t, b.CDF(5), 1, tol, "Binomial CDF(N)")
	approx(t, b.Mean(), 2.5, tol, "Binomial Mean")
	approx(t, b.Variance(), 1.25, tol, "Binomial Variance")
	// PMF sums to 1.
	sum := 0.0
	for k := 0; k <= 5; k++ {
		sum += b.PMF(k)
	}
	approx(t, sum, 1, tol, "Binomial PMF sum")
	// Asymmetric hand calc: C(4,1)*0.2*0.8^3 = 4*0.2*0.512 = 0.4096
	approx(t, Binomial{N: 4, P: 0.2}.PMF(1), 0.4096, tol, "Binomial PMF(1;4,0.2)")
}

func TestPoisson(t *testing.T) {
	p := Poisson{Lambda: 3}
	// e^-3 * 3^2 / 2 = 0.224042
	approx(t, p.PMF(2), math.Exp(-3)*9/2, tol, "Poisson PMF(2;3)")
	approx(t, p.PMF(0), math.Exp(-3), tol, "Poisson PMF(0)")
	approx(t, p.Mean(), 3, tol, "Poisson Mean")
	approx(t, p.Variance(), 3, tol, "Poisson Variance")
	sum := 0.0
	for k := 0; k <= 60; k++ {
		sum += p.PMF(k)
	}
	approx(t, sum, 1, 1e-9, "Poisson PMF sum")
}

func TestUniform(t *testing.T) {
	u := Uniform{A: 0, B: 10}
	approx(t, u.PDF(5), 0.1, tol, "Uniform PDF inside")
	approx(t, u.PDF(-1), 0, tol, "Uniform PDF outside")
	approx(t, u.CDF(2.5), 0.25, tol, "Uniform CDF")
	approx(t, u.Quantile(0.25), 2.5, tol, "Uniform Quantile")
	approx(t, u.Mean(), 5, tol, "Uniform Mean")
	approx(t, u.Variance(), 100.0/12, tol, "Uniform Variance")
}

func TestExponential(t *testing.T) {
	e := Exponential{Lambda: 2}
	approx(t, e.PDF(0), 2, tol, "Exponential PDF(0)")
	approx(t, e.CDF(1), 1-math.Exp(-2), tol, "Exponential CDF(1)")
	approx(t, e.Mean(), 0.5, tol, "Exponential Mean")
	approx(t, e.Variance(), 0.25, tol, "Exponential Variance")
	// Median is ln(2)/lambda.
	approx(t, e.Quantile(0.5), math.Ln2/2, tol, "Exponential Quantile 0.5")
	approx(t, e.CDF(e.Quantile(0.3)), 0.3, 1e-9, "Exponential roundtrip")
}

func TestStudentT(t *testing.T) {
	tt := StudentT{Nu: 10}
	approx(t, tt.CDF(0), 0.5, tol, "StudentT CDF(0)")
	// t_{0.975, 10} = 2.228
	approx(t, tt.CDF(2.228138852), 0.975, 1e-6, "StudentT CDF at 0.975 quantile")
	approx(t, tt.CDF(-2.228138852), 0.025, 1e-6, "StudentT CDF negative")
	approx(t, tt.Mean(), 0, tol, "StudentT Mean")
	approx(t, tt.Variance(), 10.0/8, tol, "StudentT Variance")
	if !math.IsNaN(StudentT{Nu: 1}.Mean()) {
		t.Errorf("StudentT Mean for nu=1 should be NaN")
	}
	// Large nu approaches the normal distribution.
	approx(t, StudentT{Nu: 100000}.CDF(1.96), NormalCDF(1.96, 0, 1), 1e-3, "StudentT -> Normal")
}

func TestChiSquared(t *testing.T) {
	c := ChiSquared{K: 2}
	// For k=2, CDF(x) = 1 - e^{-x/2}.
	approx(t, c.CDF(2), 1-math.Exp(-1), 1e-9, "ChiSquared CDF(2;2)")
	approx(t, c.PDF(0), 0.5, tol, "ChiSquared PDF(0;2)")
	approx(t, c.Mean(), 2, tol, "ChiSquared Mean")
	approx(t, c.Variance(), 4, tol, "ChiSquared Variance")
	// k=1, CDF(3.841) ~ 0.95 (classic critical value).
	approx(t, ChiSquared{K: 1}.CDF(3.841458821), 0.95, 1e-6, "ChiSquared 95% critical")
}

func TestGamma(t *testing.T) {
	g := Gamma{Shape: 2, Scale: 2}
	approx(t, g.Mean(), 4, tol, "Gamma Mean")
	approx(t, g.Variance(), 8, tol, "Gamma Variance")
	// Gamma(shape=1, scale=1/lambda) is Exponential(lambda).
	approx(t, Gamma{Shape: 1, Scale: 0.5}.CDF(1), Exponential{Lambda: 2}.CDF(1), 1e-9, "Gamma == Exponential")
	// Gamma(shape=k/2, scale=2) is ChiSquared(k).
	approx(t, Gamma{Shape: 1, Scale: 2}.CDF(2), ChiSquared{K: 2}.CDF(2), 1e-9, "Gamma == ChiSquared")
}

func TestLinearRegression(t *testing.T) {
	// Perfectly collinear: y = 2x + 1.
	xs := []float64{1, 2, 3, 4, 5}
	ys := []float64{3, 5, 7, 9, 11}
	slope, intercept, r2 := LinearRegression(xs, ys)
	approx(t, slope, 2, tol, "slope")
	approx(t, intercept, 1, tol, "intercept")
	approx(t, r2, 1, tol, "r2")
	// Noisy but strongly linear data still has a sane fit.
	xs2 := []float64{0, 1, 2, 3, 4}
	ys2 := []float64{1, 3, 2, 5, 4}
	s2, i2, r22 := LinearRegression(xs2, ys2)
	approx(t, s2, 0.8, tol, "noisy slope")
	approx(t, i2, 1.4, tol, "noisy intercept")
	if r22 <= 0 || r22 >= 1 {
		t.Errorf("noisy r2 = %v, want in (0,1)", r22)
	}
	// Degenerate: vertical line (zero x-variance).
	_, _, r2v := LinearRegression([]float64{2, 2, 2}, []float64{1, 2, 3})
	if !math.IsNaN(r2v) {
		t.Errorf("vertical regression r2 should be NaN, got %v", r2v)
	}
}

func TestSpecialFunctions(t *testing.T) {
	// gammaLn(n+1) == ln(n!)
	approx(t, gammaLn(6), math.Log(120), tol, "gammaLn(6)")
	// Regularized gamma P endpoints and a known value.
	approx(t, regularizedGammaP(1, 1), 1-math.Exp(-1), 1e-12, "P(1,1)")
	approx(t, regularizedGammaP(2, 0), 0, tol, "P(a,0)")
	approx(t, regularizedGammaP(3, 100), 1, 1e-12, "P(a, large)")
	approx(t, regularizedGammaP(2, 5)+regularizedGammaQ(2, 5), 1, tol, "P+Q=1")
	// Regularized incomplete beta symmetry: I_x(a,b) = 1 - I_{1-x}(b,a).
	approx(t, regularizedIncompleteBeta(2, 3, 0.4),
		1-regularizedIncompleteBeta(3, 2, 0.6), 1e-12, "beta symmetry")
	approx(t, regularizedIncompleteBeta(2, 2, 0.5), 0.5, 1e-12, "beta symmetric point")
}

func TestDistributionPDFsAndEdges(t *testing.T) {
	// StudentT PDF: symmetric, peak at 0; integrates via a coarse check.
	tt := StudentT{Nu: 5}
	approx(t, tt.PDF(1), tt.PDF(-1), tol, "StudentT PDF symmetric")
	if tt.PDF(0) <= tt.PDF(2) {
		t.Errorf("StudentT PDF should peak at 0")
	}
	// Gamma/ChiSquared PDF known point: ChiSquared(k=2) PDF(x)=0.5*e^{-x/2}.
	approx(t, ChiSquared{K: 2}.PDF(2), 0.5*math.Exp(-1), 1e-12, "ChiSquared PDF(2)")
	// Gamma(1, theta) is Exponential(1/theta): PDF(x) = (1/theta) e^{-x/theta}.
	approx(t, Gamma{Shape: 1, Scale: 2}.PDF(1), 0.5*math.Exp(-0.5), 1e-12, "Gamma PDF exp")
	// PDF/CDF below support are zero.
	approx(t, ChiSquared{K: 3}.PDF(-1), 0, tol, "ChiSquared PDF negative")
	approx(t, ChiSquared{K: 3}.CDF(-1), 0, tol, "ChiSquared CDF negative")
	approx(t, Gamma{Shape: 2, Scale: 1}.PDF(-1), 0, tol, "Gamma PDF negative")
	approx(t, Gamma{Shape: 2, Scale: 1}.CDF(-1), 0, tol, "Gamma CDF negative")
	approx(t, Exponential{Lambda: 1}.PDF(-1), 0, tol, "Exponential PDF negative")
	approx(t, Exponential{Lambda: 1}.CDF(-1), 0, tol, "Exponential CDF negative")
	// Origin behaviour of shape-dependent densities.
	approx(t, ChiSquared{K: 2}.PDF(0), 0.5, tol, "ChiSquared PDF(0) k=2")
	if !math.IsInf(ChiSquared{K: 1}.PDF(0), 1) {
		t.Errorf("ChiSquared PDF(0) for k<2 should be +Inf")
	}
	approx(t, ChiSquared{K: 3}.PDF(0), 0, tol, "ChiSquared PDF(0) k>2")
	approx(t, Gamma{Shape: 1, Scale: 2}.PDF(0), 0.5, tol, "Gamma PDF(0) shape=1")
	if !math.IsInf(Gamma{Shape: 0.5, Scale: 1}.PDF(0), 1) {
		t.Errorf("Gamma PDF(0) for shape<1 should be +Inf")
	}
	approx(t, Gamma{Shape: 2, Scale: 1}.PDF(0), 0, tol, "Gamma PDF(0) shape>1")
}

func TestDistributionDegenerateAndQuantiles(t *testing.T) {
	// Binomial with degenerate p.
	approx(t, Binomial{N: 3, P: 0}.PMF(0), 1, tol, "Binomial p=0 PMF(0)")
	approx(t, Binomial{N: 3, P: 0}.PMF(1), 0, tol, "Binomial p=0 PMF(1)")
	approx(t, Binomial{N: 3, P: 1}.PMF(3), 1, tol, "Binomial p=1 PMF(N)")
	approx(t, Binomial{N: 3, P: 1}.PMF(0), 0, tol, "Binomial p=1 PMF(0)")
	approx(t, Binomial{N: 5, P: 0.5}.CDF(-1), 0, tol, "Binomial CDF(<0)")
	approx(t, Poisson{Lambda: 2}.PMF(-1), 0, tol, "Poisson PMF(<0)")
	approx(t, Poisson{Lambda: 2}.CDF(-1), 0, tol, "Poisson CDF(<0)")
	// Uniform/Exponential quantile out of range -> NaN.
	if !math.IsNaN(Uniform{A: 0, B: 1}.Quantile(2)) {
		t.Errorf("Uniform Quantile out of range should be NaN")
	}
	if !math.IsNaN(Exponential{Lambda: 1}.Quantile(-0.1)) {
		t.Errorf("Exponential Quantile out of range should be NaN")
	}
	// StudentT variance branches.
	if !math.IsInf(StudentT{Nu: 1.5}.Variance(), 1) {
		t.Errorf("StudentT Variance for 1<nu<=2 should be +Inf")
	}
	if !math.IsNaN(StudentT{Nu: 0.5}.Variance()) {
		t.Errorf("StudentT Variance for nu<=1 should be NaN")
	}
	// normQuantile boundary values.
	if !math.IsInf(Normal{Mu: 0, Sigma: 1}.Quantile(0), -1) {
		t.Errorf("Normal Quantile(0) should be -Inf")
	}
	if !math.IsInf(Normal{Mu: 0, Sigma: 1}.Quantile(1), 1) {
		t.Errorf("Normal Quantile(1) should be +Inf")
	}
	if !math.IsNaN(Normal{Mu: 0, Sigma: 1}.Quantile(1.5)) {
		t.Errorf("Normal Quantile out of range should be NaN")
	}
	// Extreme-tail quantiles exercise the outer Acklam branches.
	approx(t, Normal{Mu: 0, Sigma: 1}.CDF(Normal{Mu: 0, Sigma: 1}.Quantile(0.001)), 0.001, 1e-9, "far low tail")
	approx(t, Normal{Mu: 0, Sigma: 1}.CDF(Normal{Mu: 0, Sigma: 1}.Quantile(0.999)), 0.999, 1e-9, "far high tail")
	// Combinatorics edge branches.
	if !math.IsInf(Factorial(200), 1) {
		t.Errorf("Factorial(200) should overflow to +Inf")
	}
	approx(t, Perm(5, 0), 1, tol, "Perm k=0")
	if !math.IsNaN(Choose(-1, 0)) {
		t.Errorf("Choose(-1,0) should be NaN")
	}
	if !math.IsNaN(Perm(-1, 0)) {
		t.Errorf("Perm(-1,0) should be NaN")
	}
}
