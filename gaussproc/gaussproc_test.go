package gaussproc

import (
	"fmt"
	"math"
	"math/rand"
	"testing"
)

const tol = 1e-9

func approx(a, b, eps float64) bool { return math.Abs(a-b) <= eps }

func TestKernelValues(t *testing.T) {
	x := []float64{0}
	y := []float64{1}
	tests := []struct {
		name string
		k    Kernel
		x, y []float64
		want float64
		eps  float64
	}{
		{"RBF", NewRBF(1, 1), x, y, math.Exp(-0.5), tol},
		{"RBFsameVar", NewRBF(2.5, 1), x, x, 2.5, tol},
		{"Matern12", NewMatern12(1, 1), x, y, math.Exp(-1), tol},
		{"Matern32", NewMatern32(1, 1), x, y, (1 + math.Sqrt(3)) * math.Exp(-math.Sqrt(3)), tol},
		{"Matern52", NewMatern52(1, 1), x, y, (1 + math.Sqrt(5) + 5.0/3.0) * math.Exp(-math.Sqrt(5)), tol},
		{"RationalQuadratic", NewRationalQuadratic(1, 1, 1), x, y, 1.0 / 1.5, tol},
		{"Periodic", NewPeriodic(1, 1, 1), []float64{0}, []float64{0.5}, math.Exp(-2), tol},
		{"Constant", NewConstant(3.3), x, y, 3.3, tol},
		{"WhiteNoiseEq", NewWhiteNoise(0.7), x, x, 0.7, tol},
		{"WhiteNoiseNeq", NewWhiteNoise(0.7), x, y, 0, tol},
		{"Linear", NewLinear(0, 1, 0), []float64{2}, []float64{3}, 6, tol},
		{"LinearBias", NewLinear(1, 2, 0), []float64{2}, []float64{3}, 1 + 2*6, tol},
		{"Polynomial", NewPolynomial(1, 1, 2), []float64{1, 2}, []float64{1, 1}, 16, tol},
		{"Cosine", NewCosine(1, 1), []float64{0}, []float64{0.5}, math.Cos(math.Pi), tol},
		{"GammaExp2", NewGammaExponential(1, 1, 2), x, y, math.Exp(-1), tol},
		{"GammaExp1", NewGammaExponential(1, 1, 1), x, y, math.Exp(-1), tol},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := tt.k.Eval(tt.x, tt.y)
			if !approx(got, tt.want, tt.eps) {
				t.Fatalf("%s: got %v want %v", tt.name, got, tt.want)
			}
			// symmetry
			if sym := tt.k.Eval(tt.y, tt.x); !approx(sym, got, tt.eps) {
				t.Fatalf("%s: not symmetric: %v vs %v", tt.name, got, sym)
			}
		})
	}
}

func TestGammaExponentialMatchesRBF(t *testing.T) {
	rbf := NewRBF(1, 1.3)
	// gamma-exponential with gamma=2 uses exp(-(r/l)^2); RBF uses
	// exp(-r^2/(2 l^2)). They match if length scales differ by sqrt(2).
	ge := NewGammaExponential(1, 1.3*math.Sqrt(2), 2)
	x, y := []float64{0.2, -1}, []float64{1.1, 0.5}
	if !approx(rbf.Eval(x, y), ge.Eval(x, y), 1e-12) {
		t.Fatalf("mismatch: %v vs %v", rbf.Eval(x, y), ge.Eval(x, y))
	}
}

func TestKernelAlgebra(t *testing.T) {
	a := NewRBF(1, 1)
	b := NewConstant(2)
	x, y := []float64{0}, []float64{1}
	if got, want := Sum(a, b).Eval(x, y), a.Eval(x, y)+2; !approx(got, want, tol) {
		t.Fatalf("sum: %v want %v", got, want)
	}
	if got, want := Product(a, b).Eval(x, y), a.Eval(x, y)*2; !approx(got, want, tol) {
		t.Fatalf("product: %v want %v", got, want)
	}
	if got, want := Scale(3, a).Eval(x, y), 3*a.Eval(x, y); !approx(got, want, tol) {
		t.Fatalf("scale: %v want %v", got, want)
	}
	if got, want := Exponentiate(a, 3).Eval(x, y), math.Pow(a.Eval(x, y), 3); !approx(got, want, tol) {
		t.Fatalf("exponentiate: %v want %v", got, want)
	}
}

func TestGramMatrix(t *testing.T) {
	x := [][]float64{{0}, {1}, {2}}
	k := NewRBF(1, 1)
	g := GramMatrix(k, x)
	if !g.IsSymmetric(1e-12) {
		t.Fatal("Gram matrix not symmetric")
	}
	for i := range x {
		if !approx(g[i][i], 1, tol) {
			t.Fatalf("diagonal %d = %v want 1", i, g[i][i])
		}
	}
	if !approx(g[0][1], math.Exp(-0.5), tol) {
		t.Fatalf("g[0][1] = %v", g[0][1])
	}
	if !IsKernelPSD(k, x, 1e-12) {
		t.Fatal("expected PSD Gram matrix")
	}
}

func TestCholeskyRoundTrip(t *testing.T) {
	a := Matrix{
		{4, 2, 2},
		{2, 5, 1},
		{2, 1, 3},
	}
	l, err := Cholesky(a)
	if err != nil {
		t.Fatal(err)
	}
	recon := MatMul(l, l.Transpose())
	if d := MaxAbsMatrixDiff(recon, a); d > 1e-12 {
		t.Fatalf("reconstruction error %v", d)
	}
	// solve A x = b
	b := []float64{1, 2, 3}
	x, err := CholeskySolve(l, b)
	if err != nil {
		t.Fatal(err)
	}
	if d := MaxAbsDiff(MatVec(a, x), b); d > 1e-10 {
		t.Fatalf("solve residual %v", d)
	}
	// log determinant
	logdet := LogDetFromCholesky(l)
	// det of a computed directly: 4*(5*3-1)-2*(2*3-2)+2*(2-10) = 4*14 -2*4 +2*(-8)=56-8-16=32
	if !approx(math.Exp(logdet), 32, 1e-8) {
		t.Fatalf("det = %v want 32", math.Exp(logdet))
	}
}

func TestCholeskyNotPD(t *testing.T) {
	a := Matrix{{1, 2}, {2, 1}}
	if _, err := Cholesky(a); err == nil {
		t.Fatal("expected non-positive-definite error")
	}
}

func TestInvertSPD(t *testing.T) {
	a := Matrix{{2, 0.5}, {0.5, 1}}
	inv, err := InvertSPD(a)
	if err != nil {
		t.Fatal(err)
	}
	prod := MatMul(a, inv)
	if d := MaxAbsMatrixDiff(prod, Identity(2)); d > 1e-12 {
		t.Fatalf("A A^-1 != I, err %v", d)
	}
}

func TestGPInterpolation(t *testing.T) {
	x := [][]float64{{-2}, {-1}, {0}, {1}, {2}}
	y := []float64{4, 1, 0, 1, 4} // y = x^2
	g := NewGP(NewRBF(1, 1), 1e-8)
	if err := g.Fit(x, y); err != nil {
		t.Fatal(err)
	}
	mean, variance, err := g.Predict(x)
	if err != nil {
		t.Fatal(err)
	}
	for i := range x {
		if !approx(mean[i], y[i], 1e-3) {
			t.Fatalf("mean[%d]=%v want %v", i, mean[i], y[i])
		}
		if variance[i] < -1e-9 || variance[i] > 1e-3 {
			t.Fatalf("variance[%d]=%v expected near 0", i, variance[i])
		}
	}
}

func TestGPPredictBetweenPoints(t *testing.T) {
	x := [][]float64{{0}, {1}, {2}, {3}}
	y := []float64{0, 1, 2, 3} // linear
	g := NewGP(NewRBF(1, 1.5), 1e-6)
	if err := g.Fit(x, y); err != nil {
		t.Fatal(err)
	}
	mean, variance, err := g.Predict([][]float64{{1.5}})
	if err != nil {
		t.Fatal(err)
	}
	if mean[0] < 0.5 || mean[0] > 2.5 {
		t.Fatalf("interpolated mean out of range: %v", mean[0])
	}
	if variance[0] <= 0 {
		t.Fatalf("expected positive predictive variance, got %v", variance[0])
	}
}

func TestLogMarginalLikelihoodOnePoint(t *testing.T) {
	x := [][]float64{{0.5}}
	y := []float64{2}
	k := NewRBF(1, 1)
	noise := 0.1
	got, err := LogMarginalLikelihood(k, x, y, noise)
	if err != nil {
		t.Fatal(err)
	}
	// K00 = 1, variance = 1 + noise (jitter negligible relative to eps)
	want := GaussianLogPDF(2, 0, 1+noise)
	if !approx(got, want, 1e-6) {
		t.Fatalf("LML=%v want %v", got, want)
	}
}

func TestPredictCovarianceDiagonalMatchesVariance(t *testing.T) {
	x := [][]float64{{0}, {1}, {2}}
	y := []float64{1, -1, 0.5}
	g := NewGP(NewMatern32(1, 1), 0.05)
	if err := g.Fit(x, y); err != nil {
		t.Fatal(err)
	}
	xs := [][]float64{{0.5}, {1.5}}
	variance, err := g.PredictVariance(xs)
	if err != nil {
		t.Fatal(err)
	}
	cov, err := g.PredictCovariance(xs)
	if err != nil {
		t.Fatal(err)
	}
	for i := range xs {
		if !approx(cov[i][i], variance[i], 1e-9) {
			t.Fatalf("cov diag %d = %v want %v", i, cov[i][i], variance[i])
		}
	}
	if !cov.IsSymmetric(1e-9) {
		t.Fatal("predictive covariance not symmetric")
	}
}

func TestKernelRidgeRegression(t *testing.T) {
	x := [][]float64{{0}, {1}, {2}, {3}, {4}}
	y := []float64{0, 1, 4, 9, 16} // x^2
	f, err := KernelRidgeRegression(NewRBF(4, 1.5), x, y, 1e-8)
	if err != nil {
		t.Fatal(err)
	}
	for i := range x {
		if got := f.Eval(x[i]); !approx(got, y[i], 1e-2) {
			t.Fatalf("KRR fit at x[%d]: got %v want %v", i, got, y[i])
		}
	}
}

func TestRKHSReproducing(t *testing.T) {
	k := NewRBF(1, 1)
	x0 := []float64{0.3}
	// f = k(x0, .) has <f,f>_H = k(x0,x0)
	f := NewRKHSFunction(k, [][]float64{x0}, []float64{1})
	if !approx(f.SquaredNorm(), k.Eval(x0, x0), tol) {
		t.Fatalf("||f||^2 = %v want %v", f.SquaredNorm(), k.Eval(x0, x0))
	}
	// reproducing property: <f, k(z,.)> = f(z)
	z := []float64{1.2}
	g := NewRKHSFunction(k, [][]float64{z}, []float64{1})
	if !approx(f.InnerProduct(g), f.Eval(z), tol) {
		t.Fatalf("reproducing failed: %v vs %v", f.InnerProduct(g), f.Eval(z))
	}
}

func TestRKHSDistanceAndAlgebra(t *testing.T) {
	k := NewRBF(1, 1)
	pts := [][]float64{{0}, {1}}
	f := NewRKHSFunction(k, pts, []float64{1, 2})
	g := NewRKHSFunction(k, pts, []float64{1, 2})
	if d := RKHSDistance(f, g); d > 1e-12 {
		t.Fatalf("distance to identical function %v", d)
	}
	// (f + (-1)*g) has zero norm
	diff := f.Add(g.Scale(-1))
	if diff.Norm() > 1e-9 {
		t.Fatalf("norm of f-g = %v", diff.Norm())
	}
	// scaling
	h := f.Scale(2)
	if !approx(h.Eval(pts[0]), 2*f.Eval(pts[0]), tol) {
		t.Fatal("scale eval mismatch")
	}
}

func TestSampleReproducible(t *testing.T) {
	x := [][]float64{{0}, {0.5}, {1}}
	k := NewRBF(1, 1)
	s1, err := SamplePrior(k, x, 1e-9, 2, rand.New(rand.NewSource(42)))
	if err != nil {
		t.Fatal(err)
	}
	s2, err := SamplePrior(k, x, 1e-9, 2, rand.New(rand.NewSource(42)))
	if err != nil {
		t.Fatal(err)
	}
	for i := range s1 {
		if MaxAbsDiff(s1[i], s2[i]) > tol {
			t.Fatal("prior samples not reproducible for same seed")
		}
	}
}

func TestPosteriorSampleMeanConverges(t *testing.T) {
	x := [][]float64{{0}, {1}, {2}}
	y := []float64{0.5, -0.5, 1}
	g := NewGP(NewRBF(1, 1), 0.01)
	if err := g.Fit(x, y); err != nil {
		t.Fatal(err)
	}
	xs := [][]float64{{0.5}, {1.5}}
	analMean, err := g.PredictMean(xs)
	if err != nil {
		t.Fatal(err)
	}
	rng := rand.New(rand.NewSource(7))
	const N = 4000
	sum := make([]float64, len(xs))
	samples, err := g.Sample(xs, N, rng)
	if err != nil {
		t.Fatal(err)
	}
	for _, s := range samples {
		for i := range s {
			sum[i] += s[i]
		}
	}
	for i := range xs {
		empMean := sum[i] / N
		if !approx(empMean, analMean[i], 0.05) {
			t.Fatalf("empirical mean %v vs analytic %v", empMean, analMean[i])
		}
	}
}

func TestMaternKernelDispatch(t *testing.T) {
	x, y := []float64{0}, []float64{1}
	if got := MaternKernel(1.5, 1, 1).Eval(x, y); !approx(got, NewMatern32(1, 1).Eval(x, y), tol) {
		t.Fatal("MaternKernel dispatch mismatch")
	}
}

func TestErrorPaths(t *testing.T) {
	g := NewGP(NewRBF(1, 1), 0.1)
	if _, err := g.PredictMean([][]float64{{0}}); err == nil {
		t.Fatal("expected error predicting before fit")
	}
	if err := g.Fit([][]float64{{0}, {1}}, []float64{1}); err == nil {
		t.Fatal("expected dimension mismatch error")
	}
	if err := g.Fit(nil, nil); err == nil {
		t.Fatal("expected empty error")
	}
}

func ExampleGP() {
	// Fit a Gaussian process to noisy observations of a simple function and
	// predict at a new location.
	x := [][]float64{{-1}, {0}, {1}}
	y := []float64{1, 0, 1}

	gp := NewGP(NewRBF(1, 1), 1e-6)
	if err := gp.Fit(x, y); err != nil {
		panic(err)
	}
	mean, variance, err := gp.Predict([][]float64{{0.5}})
	if err != nil {
		panic(err)
	}
	fmt.Printf("mean=%.3f var=%.3f\n", mean[0], variance[0])
	// Output: mean=0.342 var=0.018
}

func ExampleRBF() {
	k := NewRBF(1, 1)
	fmt.Printf("%.5f\n", k.Eval([]float64{0}, []float64{1}))
	// Output: 0.60653
}
