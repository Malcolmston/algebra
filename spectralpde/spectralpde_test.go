package spectralpde

import (
	"fmt"
	"math"
	"testing"
)

func approx(a, b, tol float64) bool { return math.Abs(a-b) <= tol }

func vecApprox(a, b []float64, tol float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if math.Abs(a[i]-b[i]) > tol {
			return false
		}
	}
	return true
}

func TestChebyshevGaussLobattoNodes(t *testing.T) {
	got := ChebyshevGaussLobattoNodes(2)
	want := []float64{1, 0, -1}
	if !vecApprox(got, want, 1e-12) {
		t.Fatalf("CGL(2) = %v, want %v", got, want)
	}
	got4 := ChebyshevGaussLobattoNodes(4)
	r := math.Sqrt2 / 2
	want4 := []float64{1, r, 0, -r, -1}
	if !vecApprox(got4, want4, 1e-12) {
		t.Fatalf("CGL(4) = %v, want %v", got4, want4)
	}
}

func TestChebyshevT(t *testing.T) {
	tests := []struct {
		n    int
		x    float64
		want float64
	}{
		{0, 0.3, 1},
		{1, 0.3, 0.3},
		{2, 0.5, -0.5},
		{3, 0.5, -1},
		{4, 1, 1},
		{5, -1, -1},
	}
	for _, tc := range tests {
		if got := ChebyshevT(tc.n, tc.x); !approx(got, tc.want, 1e-12) {
			t.Errorf("T_%d(%v) = %v, want %v", tc.n, tc.x, got, tc.want)
		}
	}
}

func TestChebyshevDerivativeMatrixPolynomial(t *testing.T) {
	N := 8
	x := ChebyshevGaussLobattoNodes(N)
	f := make([]float64, N+1)
	want := make([]float64, N+1)
	for i, xi := range x {
		f[i] = xi * xi * xi // x^3
		want[i] = 3 * xi * xi
	}
	D := ChebyshevDiffMatrix(N)
	got := MatVec(D, f)
	if !vecApprox(got, want, 1e-10) {
		t.Fatalf("d/dx x^3 mismatch: got %v want %v", got, want)
	}
}

func TestChebyshevSecondDerivativeMatrix(t *testing.T) {
	N := 10
	x := ChebyshevGaussLobattoNodes(N)
	f := make([]float64, N+1)
	want := make([]float64, N+1)
	for i, xi := range x {
		f[i] = math.Cos(xi)
		want[i] = -math.Cos(xi)
	}
	D2 := ChebyshevDiffMatrix2(N)
	got := MatVec(D2, f)
	if !vecApprox(got, want, 1e-6) {
		t.Fatalf("second derivative of cos mismatch")
	}
}

func TestFourierDiffMatrix(t *testing.T) {
	N := 16
	x := FourierNodes(N)
	f := make([]float64, N)
	want := make([]float64, N)
	for i, xi := range x {
		f[i] = math.Sin(xi)
		want[i] = math.Cos(xi)
	}
	D := FourierDiffMatrix(N)
	got := MatVec(D, f)
	if !vecApprox(got, want, 1e-10) {
		t.Fatalf("fourier d/dx sin mismatch: got %v", got)
	}
}

func TestFourierDiffMatrix2(t *testing.T) {
	N := 16
	x := FourierNodes(N)
	f := make([]float64, N)
	want := make([]float64, N)
	for i, xi := range x {
		f[i] = math.Sin(2 * xi)
		want[i] = -4 * math.Sin(2*xi)
	}
	D2 := FourierDiffMatrix2(N)
	got := MatVec(D2, f)
	if !vecApprox(got, want, 1e-9) {
		t.Fatalf("fourier d2/dx2 sin(2x) mismatch: got %v", got)
	}
}

func TestFFTvsDFT(t *testing.T) {
	x := RealToComplex([]float64{1, 2, 3, 4, 5, 6, 7, 8})
	a := FFT(x)
	b := DFT(x)
	for i := range a {
		if math.Abs(real(a[i])-real(b[i])) > 1e-9 || math.Abs(imag(a[i])-imag(b[i])) > 1e-9 {
			t.Fatalf("FFT vs DFT mismatch at %d: %v vs %v", i, a[i], b[i])
		}
	}
}

func TestFFTRoundTrip(t *testing.T) {
	orig := []float64{3, 1, 4, 1, 5, 9, 2, 6}
	back := ComplexReal(IFFT(FFTReal(orig)))
	if !vecApprox(orig, back, 1e-10) {
		t.Fatalf("FFT round trip failed: %v", back)
	}
	// Non-power-of-two path.
	orig3 := []float64{1, 2, 3, 4, 5}
	back3 := ComplexReal(IDFT(DFT(RealToComplex(orig3))))
	if !vecApprox(orig3, back3, 1e-10) {
		t.Fatalf("DFT round trip failed: %v", back3)
	}
}

func TestChebyshevCoefficients(t *testing.T) {
	N := 5
	x := ChebyshevGaussLobattoNodes(N)
	// f = T_3.
	f := make([]float64, N+1)
	for i := range x {
		f[i] = ChebyshevT(3, x[i])
	}
	c := ChebyshevCoefficients(f)
	want := []float64{0, 0, 0, 1, 0, 0}
	if !vecApprox(c, want, 1e-10) {
		t.Fatalf("Chebyshev coeffs of T_3 = %v", c)
	}
	// Round trip.
	back := ChebyshevValuesFromCoeffs(c)
	if !vecApprox(back, f, 1e-10) {
		t.Fatalf("coeff round trip failed")
	}
}

func TestChebyshevInterpolate(t *testing.T) {
	f := func(x float64) float64 { return math.Exp(x) }
	coeffs := ChebyshevFit(f, 20)
	for _, x := range []float64{-0.9, -0.3, 0.0, 0.37, 0.88} {
		got := ClenshawEval(coeffs, x)
		if !approx(got, f(x), 1e-11) {
			t.Errorf("exp interp at %v: got %v want %v", x, got, f(x))
		}
	}
}

func TestChebyshevDifferentiateCoeffs(t *testing.T) {
	// f = T_2, f' = 4 T_1.
	d := ChebyshevDifferentiateCoeffs([]float64{0, 0, 1})
	if !vecApprox(d, []float64{0, 4}, 1e-12) {
		t.Fatalf("d T_2 = %v, want [0 4]", d)
	}
}

func TestChebyshevIntegral(t *testing.T) {
	// x^2 = (T_0 + T_2)/2, integral over [-1,1] = 2/3.
	coeffs := []float64{0.5, 0, 0.5}
	if got := ChebyshevIntegral(coeffs); !approx(got, 2.0/3.0, 1e-12) {
		t.Fatalf("integral x^2 = %v want 0.6667", got)
	}
}

func TestClenshawCurtis(t *testing.T) {
	w := ClenshawCurtisWeights(8)
	if s := Sum(w); !approx(s, 2, 1e-12) {
		t.Fatalf("CC weights sum = %v want 2", s)
	}
	// Integrate x^4 on [-1,1] = 2/5.
	got := ClenshawCurtisIntegrate(func(x float64) float64 { return x * x * x * x }, 8, -1, 1)
	if !approx(got, 0.4, 1e-12) {
		t.Fatalf("CC integral x^4 = %v want 0.4", got)
	}
	// Integrate exp on [0,2].
	e := ClenshawCurtisIntegrate(math.Exp, 24, 0, 2)
	if !approx(e, math.Exp(2)-1, 1e-10) {
		t.Fatalf("CC integral exp = %v", e)
	}
}

func TestFejer(t *testing.T) {
	w := FejerWeights(10)
	if s := Sum(w); !approx(s, 2, 1e-12) {
		t.Fatalf("Fejer weights sum = %v want 2", s)
	}
	got := FejerIntegrate(func(x float64) float64 { return x * x }, 10, -1, 1)
	if !approx(got, 2.0/3.0, 1e-12) {
		t.Fatalf("Fejer integral x^2 = %v", got)
	}
}

func TestGaussLegendre(t *testing.T) {
	nodes, weights := LegendreGaussNodesWeights(2)
	inv := 1 / math.Sqrt(3)
	if !vecApprox(nodes, []float64{-inv, inv}, 1e-12) {
		t.Fatalf("GL(2) nodes = %v", nodes)
	}
	if !vecApprox(weights, []float64{1, 1}, 1e-12) {
		t.Fatalf("GL(2) weights = %v", weights)
	}
	// Integrate x^4 exactly with n=3 points.
	got := GaussLegendreIntegrate(func(x float64) float64 { return x * x * x * x }, 3, -1, 1)
	if !approx(got, 0.4, 1e-12) {
		t.Fatalf("GL integral x^4 = %v", got)
	}
}

func TestLegendreP(t *testing.T) {
	tests := []struct {
		n    int
		x    float64
		want float64
	}{
		{0, 0.4, 1},
		{1, 0.4, 0.4},
		{2, 0.5, -0.125},
		{3, 1, 1},
		{4, 1, 1},
	}
	for _, tc := range tests {
		if got := LegendreP(tc.n, tc.x); !approx(got, tc.want, 1e-12) {
			t.Errorf("P_%d(%v) = %v want %v", tc.n, tc.x, got, tc.want)
		}
	}
}

func TestLegendreProjection(t *testing.T) {
	// Project P_2 -> coefficient vector e_2.
	c := LegendreProjection(func(x float64) float64 { return 1.5*x*x - 0.5 }, 4)
	want := []float64{0, 0, 1, 0, 0}
	if !vecApprox(c, want, 1e-10) {
		t.Fatalf("Legendre projection of P2 = %v", c)
	}
}

func TestLegendreGaussLobatto(t *testing.T) {
	nodes, weights := LegendreGaussLobattoNodesWeights(4)
	// Known LGL(4): -1, -sqrt(3/7), 0, sqrt(3/7), 1.
	s := math.Sqrt(3.0 / 7.0)
	want := []float64{-1, -s, 0, s, 1}
	if !vecApprox(nodes, want, 1e-10) {
		t.Fatalf("LGL nodes = %v want %v", nodes, want)
	}
	if sum := Sum(weights); !approx(sum, 2, 1e-12) {
		t.Fatalf("LGL weights sum = %v want 2", sum)
	}
}

func TestBarycentric(t *testing.T) {
	nodes := []float64{-1, 0, 1}
	values := []float64{1, 0, 1} // x^2
	w := BarycentricWeights(nodes)
	if got := BarycentricInterpolate(nodes, values, w, 0.5); !approx(got, 0.25, 1e-12) {
		t.Fatalf("barycentric x^2 at 0.5 = %v", got)
	}
	// Chebyshev closed-form weights reproduce values at nodes.
	N := 6
	cx := ChebyshevGaussLobattoNodes(N)
	fv := ApplyFunc(func(x float64) float64 { return math.Sin(3 * x) }, cx)
	cw := BarycentricWeightsChebyshev(N)
	for i := range cx {
		got := BarycentricInterpolate(cx, fv, cw, cx[i])
		if !approx(got, fv[i], 1e-12) {
			t.Fatalf("barycentric node reproduction failed")
		}
	}
}

func TestLagrangeAndNewton(t *testing.T) {
	nodes := []float64{0, 1, 2, 3}
	f := func(x float64) float64 { return x*x*x - 2*x + 1 }
	values := ApplyFunc(f, nodes)
	for _, x := range []float64{0.5, 1.5, 2.7} {
		if got := LagrangeInterpolate(nodes, values, x); !approx(got, f(x), 1e-10) {
			t.Errorf("Lagrange at %v = %v want %v", x, got, f(x))
		}
		coef := NewtonDividedDifferences(nodes, values)
		if got := NewtonEval(coef, nodes, x); !approx(got, f(x), 1e-10) {
			t.Errorf("Newton at %v = %v want %v", x, got, f(x))
		}
	}
}

func TestLinearAlgebra(t *testing.T) {
	A := [][]float64{{2, 1}, {1, 3}}
	b := []float64{3, 5}
	x, err := SolveLinearSystem(A, b)
	if err != nil {
		t.Fatal(err)
	}
	// Solution: x=0.8, y=1.4.
	if !vecApprox(x, []float64{0.8, 1.4}, 1e-12) {
		t.Fatalf("solve = %v", x)
	}
	det, _ := Determinant(A)
	if !approx(det, 5, 1e-12) {
		t.Fatalf("det = %v want 5", det)
	}
	inv, _ := Inverse(A)
	prod := MatMul(A, inv)
	if !vecApprox(prod[0], []float64{1, 0}, 1e-12) || !vecApprox(prod[1], []float64{0, 1}, 1e-12) {
		t.Fatalf("inverse wrong: %v", prod)
	}
}

func TestKron(t *testing.T) {
	a := [][]float64{{1, 2}, {3, 4}}
	b := [][]float64{{0, 1}, {1, 0}}
	k := Kron(a, b)
	want := [][]float64{
		{0, 1, 0, 2},
		{1, 0, 2, 0},
		{0, 3, 0, 4},
		{3, 0, 4, 0},
	}
	for i := range want {
		if !vecApprox(k[i], want[i], 1e-12) {
			t.Fatalf("kron row %d = %v want %v", i, k[i], want[i])
		}
	}
}

func TestJacobiEigen(t *testing.T) {
	A := [][]float64{{2, 1}, {1, 2}}
	vals, vecs, err := JacobiEigenSymmetric(A)
	if err != nil {
		t.Fatal(err)
	}
	if !vecApprox(vals, []float64{1, 3}, 1e-10) {
		t.Fatalf("eigenvalues = %v want [1 3]", vals)
	}
	// Verify A v = lambda v for the first eigenpair.
	v0 := []float64{vecs[0][0], vecs[1][0]}
	Av := MatVec(A, v0)
	lv := VectorScale(v0, vals[0])
	if !vecApprox(Av, lv, 1e-9) {
		t.Fatalf("eigenvector check failed: %v vs %v", Av, lv)
	}
}

func TestPoissonSolve1D(t *testing.T) {
	// u'' = f, u(x) = sin(pi x) on [-1,1], u(-1)=u(1)=0.
	exact := func(x float64) float64 { return math.Sin(math.Pi * x) }
	f := func(x float64) float64 { return -math.Pi * math.Pi * math.Sin(math.Pi*x) }
	nodes, u, err := PoissonSolve1D(f, 24, -1, 1, exact(-1), exact(1))
	if err != nil {
		t.Fatal(err)
	}
	want := ApplyFunc(exact, nodes)
	if LinfError(u, want) > 1e-8 {
		t.Fatalf("Poisson1D error = %v", LinfError(u, want))
	}
}

func TestHelmholtzSolve1D(t *testing.T) {
	// u'' + u = f, exact u = x^3 on [0,1], f = 6x + x^3.
	exact := func(x float64) float64 { return x * x * x }
	f := func(x float64) float64 { return 6*x + x*x*x }
	nodes, u, err := HelmholtzSolve1D(f, 1, 16, 0, 1, exact(0), exact(1))
	if err != nil {
		t.Fatal(err)
	}
	want := ApplyFunc(exact, nodes)
	if LinfError(u, want) > 1e-9 {
		t.Fatalf("Helmholtz1D error = %v", LinfError(u, want))
	}
}

func TestPoissonSolve1DFourier(t *testing.T) {
	// u'' = sin(x) on [0,2pi) => u = -sin(x).
	f := func(x float64) float64 { return math.Sin(x) }
	nodes, u := PoissonSolve1DFourier(f, 32, 2*math.Pi)
	want := make([]float64, len(nodes))
	for i, x := range nodes {
		want[i] = -math.Sin(x)
	}
	if LinfError(u, want) > 1e-10 {
		t.Fatalf("Fourier Poisson error = %v", LinfError(u, want))
	}
}

func TestHeatSolveFourier(t *testing.T) {
	N := 32
	x := FourierNodes(N)
	u0 := ApplyFunc(math.Sin, x)
	u := HeatSolveFourier(u0, 1.0, 0.5, 2*math.Pi)
	want := make([]float64, N)
	for i := range x {
		want[i] = math.Exp(-0.5) * math.Sin(x[i])
	}
	if LinfError(u, want) > 1e-10 {
		t.Fatalf("heat fourier error = %v", LinfError(u, want))
	}
}

func TestAdvectionSolveFourier(t *testing.T) {
	N := 32
	x := FourierNodes(N)
	u0 := ApplyFunc(func(v float64) float64 { return math.Sin(v) }, x)
	c := 1.3
	tt := 0.7
	u := AdvectionSolveFourier(u0, c, tt, 2*math.Pi)
	want := make([]float64, N)
	for i := range x {
		want[i] = math.Sin(x[i] - c*tt)
	}
	if LinfError(u, want) > 1e-9 {
		t.Fatalf("advection error = %v", LinfError(u, want))
	}
}

func TestHeatSolveChebyshev(t *testing.T) {
	// u_t = u_xx on [-1,1], u=exp(-(pi/2)^2 t) sin(pi(x+1)/2).
	mode := func(x float64) float64 { return math.Sin(math.Pi * (x + 1) / 2) }
	nu := 1.0
	dt := 0.0005
	steps := 200
	tt := dt * float64(steps)
	nodes, u, err := HeatSolveChebyshev(mode, nu, dt, steps, 24, -1, 1)
	if err != nil {
		t.Fatal(err)
	}
	decay := math.Exp(-nu * (math.Pi / 2) * (math.Pi / 2) * tt)
	want := make([]float64, len(nodes))
	for i, x := range nodes {
		want[i] = decay * mode(x)
	}
	if LinfError(u, want) > 1e-4 {
		t.Fatalf("heat chebyshev error = %v", LinfError(u, want))
	}
}

func TestPoisson2D(t *testing.T) {
	pi := math.Pi
	exact := func(x, y float64) float64 { return math.Sin(pi*x) * math.Sin(pi*y) }
	f := func(x, y float64) float64 { return -2 * pi * pi * math.Sin(pi*x) * math.Sin(pi*y) }
	g := func(x, y float64) float64 { return 0 }
	xn, yn, U, err := Poisson2D(f, g, 16, 16, -1, 1, -1, 1)
	if err != nil {
		t.Fatal(err)
	}
	var maxErr float64
	for i := range xn {
		for j := range yn {
			e := math.Abs(U[i][j] - exact(xn[i], yn[j]))
			if e > maxErr {
				maxErr = e
			}
		}
	}
	if maxErr > 1e-5 {
		t.Fatalf("Poisson2D error = %v", maxErr)
	}
}

func TestPoisson2DFourier(t *testing.T) {
	// f = -2 sin(x) sin(y) => u = sin(x) sin(y).
	f := func(x, y float64) float64 { return -2 * math.Sin(x) * math.Sin(y) }
	xn, yn, U := Poisson2DFourier(f, 16, 16, 2*math.Pi, 2*math.Pi)
	var maxErr float64
	for i := range xn {
		for j := range yn {
			e := math.Abs(U[i][j] - math.Sin(xn[i])*math.Sin(yn[j]))
			if e > maxErr {
				maxErr = e
			}
		}
	}
	if maxErr > 1e-10 {
		t.Fatalf("Poisson2DFourier error = %v", maxErr)
	}
}

func TestDCTRoundTrip(t *testing.T) {
	x := []float64{1, 2, 0.5, -1, 3, 2, 1}
	back := IDCT1(DCT1(x))
	if !vecApprox(back, x, 1e-10) {
		t.Fatalf("DCT1 round trip failed: %v", back)
	}
}

func TestDSTRoundTrip(t *testing.T) {
	x := []float64{1, -2, 3, 0.5, 2}
	back := IDST1(DST1(x))
	if !vecApprox(back, x, 1e-10) {
		t.Fatalf("DST1 round trip failed: %v", back)
	}
}

func TestSpectralConvergenceRate(t *testing.T) {
	// Errors decaying like exp(-2 n).
	ns := []int{4, 8, 12, 16}
	errs := make([]float64, len(ns))
	for i, n := range ns {
		errs[i] = math.Exp(-2 * float64(n))
	}
	r := SpectralConvergenceRate(ns, errs)
	if !approx(r, 2, 1e-6) {
		t.Fatalf("convergence rate = %v want ~2", r)
	}
	if !IsSpectrallyConverging(ns, errs, 1.5) {
		t.Fatal("expected spectral convergence")
	}
}

func TestFourierInterpolate(t *testing.T) {
	N := 16
	x := FourierNodes(N)
	f := func(v float64) float64 { return math.Sin(v) + math.Cos(3*v) }
	vals := ApplyFunc(f, x)
	for _, q := range []float64{0.3, 1.1, 2.7, 5.0} {
		if got := FourierInterpolate(vals, q); !approx(got, f(q), 1e-9) {
			t.Errorf("Fourier interp at %v = %v want %v", q, got, f(q))
		}
	}
}

func TestVectorOps(t *testing.T) {
	a := []float64{3, 4}
	if !approx(Norm2(a), 5, 1e-12) {
		t.Fatalf("norm2 = %v", Norm2(a))
	}
	if !approx(NormInf(a), 4, 1e-12) {
		t.Fatalf("normInf = %v", NormInf(a))
	}
	if !approx(DotProduct(a, []float64{1, 2}), 11, 1e-12) {
		t.Fatalf("dot wrong")
	}
}

func ExampleChebyshevT() {
	fmt.Printf("%.1f\n", ChebyshevT(3, 0.5))
	// Output: -1.0
}

func ExamplePoissonSolve1D() {
	// Solve u'' = -pi^2 sin(pi x) on [-1,1] with u(-1)=u(1)=0; exact sin(pi x).
	f := func(x float64) float64 { return -math.Pi * math.Pi * math.Sin(math.Pi*x) }
	_, u, _ := PoissonSolve1D(f, 20, -1, 1, 0, 0)
	// The value at the centre node (x=0) should be ~sin(0)=0.
	fmt.Printf("%.4f\n", u[10])
	// Output: 0.0000
}
