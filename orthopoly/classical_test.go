package orthopoly

import (
	"math"
	"testing"
)

// orthopolyClassicalClose reports whether a and b agree within an absolute
// tolerance.
func orthopolyClassicalClose(a, b, tol float64) bool {
	if math.IsNaN(a) || math.IsNaN(b) {
		return false
	}
	return math.Abs(a-b) <= tol
}

// orthopolyClassicalSlicesClose reports whether two float slices match
// elementwise within tol.
func orthopolyClassicalSlicesClose(a, b []float64, tol float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !orthopolyClassicalClose(a[i], b[i], tol) {
			return false
		}
	}
	return true
}

func TestLegendreClassical(t *testing.T) {
	cases := []struct {
		name string
		got  float64
		want float64
	}{
		{"P0", LegendreP(0, 0.5), 1},
		{"P1", LegendreP(1, 0.5), 0.5},
		{"P2", LegendreP(2, 0.5), -0.125},
		{"P3", LegendreP(3, 0.5), -0.4375},
		{"P4(1)", LegendreP(4, 1), 1},
		{"P5(-1)", LegendreP(5, -1), -1},
		{"dP2", LegendrePDerivative(2, 0.5), 1.5},
		{"dP3(1)", LegendrePDerivative(3, 1), 6},
		{"dP3(-1)", LegendrePDerivative(3, -1), 6},
		{"dP2(-1)", LegendrePDerivative(2, -1), -3},
		{"shifted P2(0.75)", ShiftedLegendreP(2, 0.75), LegendreP(2, 0.5)},
		{"norm2", LegendrePNorm(2), 0.4},
	}
	for _, c := range cases {
		if !orthopolyClassicalClose(c.got, c.want, 1e-12) {
			t.Errorf("%s = %.15g, want %.15g", c.name, c.got, c.want)
		}
	}
	if got := LegendrePValues(3, 0.5); !orthopolyClassicalSlicesClose(got, []float64{1, 0.5, -0.125, -0.4375}, 1e-12) {
		t.Errorf("LegendrePValues = %v", got)
	}
	// Orthonormality: integral of NormalizedLegendreP(2)^2 over [-1,1] == 1.
	rule := NewGaussLegendreRule(8)
	if got := rule.Integrate(func(x float64) float64 {
		v := NormalizedLegendreP(2, x)
		return v * v
	}); !orthopolyClassicalClose(got, 1, 1e-12) {
		t.Errorf("normalized Legendre norm = %.15g, want 1", got)
	}
}

func TestLegendreZerosClassical(t *testing.T) {
	r := LegendreZeros(2)
	want := []float64{-1 / math.Sqrt(3), 1 / math.Sqrt(3)}
	if !orthopolyClassicalSlicesClose(r, want, 1e-12) {
		t.Errorf("LegendreZeros(2) = %v, want %v", r, want)
	}
	for n := 1; n <= 12; n++ {
		for _, x := range LegendreZeros(n) {
			if !orthopolyClassicalClose(LegendreP(n, x), 0, 1e-11) {
				t.Errorf("P_%d(%.15g) = %.3g, want 0", n, x, LegendreP(n, x))
			}
		}
	}
	var sum float64
	for _, w := range LegendreWeights(10) {
		sum += w
	}
	if !orthopolyClassicalClose(sum, 2, 1e-12) {
		t.Errorf("sum of Gauss-Legendre weights = %.15g, want 2", sum)
	}
}

func TestChebyshevFirstClassical(t *testing.T) {
	cases := []struct {
		name string
		got  float64
		want float64
	}{
		{"T0", ChebyshevFirst(0, 0.3), 1},
		{"T1", ChebyshevFirst(1, 0.3), 0.3},
		{"T2", ChebyshevFirst(2, 0.5), -0.5},
		{"T3", ChebyshevFirst(3, 0.5), -1},
		{"T4(cos)", ChebyshevFirst(4, math.Cos(0.7)), math.Cos(4 * 0.7)},
		{"dT3(1)", ChebyshevFirstDerivative(3, 1), 9},
		{"dT2", ChebyshevFirstDerivative(2, 0.5), 4 * 0.5},
		{"normT0", ChebyshevFirstNorm(0), math.Pi},
		{"normT3", ChebyshevFirstNorm(3), math.Pi / 2},
		{"weight", ChebyshevFirstWeight(0), 1},
		{"shiftedT2", ShiftedChebyshevFirst(2, 0.75), ChebyshevFirst(2, 0.5)},
	}
	for _, c := range cases {
		if !orthopolyClassicalClose(c.got, c.want, 1e-12) {
			t.Errorf("%s = %.15g, want %.15g", c.name, c.got, c.want)
		}
	}
	roots := ChebyshevFirstRoots(2)
	if !orthopolyClassicalSlicesClose(roots, []float64{-math.Sqrt2 / 2, math.Sqrt2 / 2}, 1e-12) {
		t.Errorf("ChebyshevFirstRoots(2) = %v", roots)
	}
	ex := ChebyshevFirstExtrema(2)
	if !orthopolyClassicalSlicesClose(ex, []float64{-1, 0, 1}, 1e-12) {
		t.Errorf("ChebyshevFirstExtrema(2) = %v", ex)
	}
	if got := ClenshawChebyshevFirst([]float64{0, 0, 0, 1}, 0.4); !orthopolyClassicalClose(got, ChebyshevFirst(3, 0.4), 1e-12) {
		t.Errorf("ClenshawChebyshevFirst = %.15g, want %.15g", got, ChebyshevFirst(3, 0.4))
	}
	if got := ChebyshevFirstValues(3, 0.5); !orthopolyClassicalSlicesClose(got, []float64{1, 0.5, -0.5, -1}, 1e-12) {
		t.Errorf("ChebyshevFirstValues = %v", got)
	}
}

func TestChebyshevSecondClassical(t *testing.T) {
	cases := []struct {
		name string
		got  float64
		want float64
	}{
		{"U0", ChebyshevSecond(0, 0.3), 1},
		{"U1", ChebyshevSecond(1, 0.3), 0.6},
		{"U2", ChebyshevSecond(2, 0.5), 0},
		{"U3", ChebyshevSecond(3, 0.5), -1},
		{"dU2(1)", ChebyshevSecondDerivative(2, 1), 8},
		{"dU2(-1)", ChebyshevSecondDerivative(2, -1), -8},
		{"normU", ChebyshevSecondNorm(4), math.Pi / 2},
		{"weightU", ChebyshevSecondWeight(0), 1},
	}
	for _, c := range cases {
		if !orthopolyClassicalClose(c.got, c.want, 1e-12) {
			t.Errorf("%s = %.15g, want %.15g", c.name, c.got, c.want)
		}
	}
	// U_n(cosθ) sinθ == sin((n+1)θ).
	th := 0.9
	if got := ChebyshevSecond(4, math.Cos(th)) * math.Sin(th); !orthopolyClassicalClose(got, math.Sin(5*th), 1e-12) {
		t.Errorf("U4 identity = %.15g, want %.15g", got, math.Sin(5*th))
	}
	roots := ChebyshevSecondRoots(3)
	want := []float64{-math.Sqrt2 / 2, 0, math.Sqrt2 / 2}
	if !orthopolyClassicalSlicesClose(roots, want, 1e-12) {
		t.Errorf("ChebyshevSecondRoots(3) = %v, want %v", roots, want)
	}
}

func TestChebyshevThirdFourthClassical(t *testing.T) {
	// V_n(cosθ) = cos((n+1/2)θ)/cos(θ/2); W_n(cosθ) = sin((n+1/2)θ)/sin(θ/2).
	th := 0.8
	x := math.Cos(th)
	for n := 0; n <= 5; n++ {
		gotV := ChebyshevThird(n, x)
		wantV := math.Cos((float64(n)+0.5)*th) / math.Cos(th/2)
		if !orthopolyClassicalClose(gotV, wantV, 1e-11) {
			t.Errorf("V_%d = %.15g, want %.15g", n, gotV, wantV)
		}
		gotW := ChebyshevFourth(n, x)
		wantW := math.Sin((float64(n)+0.5)*th) / math.Sin(th/2)
		if !orthopolyClassicalClose(gotW, wantW, 1e-11) {
			t.Errorf("W_%d = %.15g, want %.15g", n, gotW, wantW)
		}
	}
}

func TestHermitePhysicistsClassical(t *testing.T) {
	cases := []struct {
		name string
		got  float64
		want float64
	}{
		{"H0", HermiteH(0, 0.5), 1},
		{"H1", HermiteH(1, 0.5), 1},
		{"H2", HermiteH(2, 0.5), -1}, // 4x^2-2
		{"H3", HermiteH(3, 0.5), -5}, // 8x^3-12x
		{"H4", HermiteH(4, 1), -20},  // 16-48+12
		{"dH3", HermiteHDerivative(3, 0.5), 6 * HermiteH(2, 0.5)},
		{"norm2", HermiteHNorm(2), 8 * math.Sqrt(math.Pi)},
	}
	for _, c := range cases {
		if !orthopolyClassicalClose(c.got, c.want, 1e-11) {
			t.Errorf("%s = %.15g, want %.15g", c.name, c.got, c.want)
		}
	}
	r := HermiteHRoots(2)
	if !orthopolyClassicalSlicesClose(r, []float64{-math.Sqrt2 / 2, math.Sqrt2 / 2}, 1e-12) {
		t.Errorf("HermiteHRoots(2) = %v", r)
	}
	// Gauss-Hermite integrates x^2 e^{-x^2} = sqrt(pi)/2.
	rule := NewGaussHermiteRule(5)
	if got := rule.Integrate(func(x float64) float64 { return x * x }); !orthopolyClassicalClose(got, math.Sqrt(math.Pi)/2, 1e-10) {
		t.Errorf("Gauss-Hermite x^2 = %.15g, want %.15g", got, math.Sqrt(math.Pi)/2)
	}
	// HermiteFunction orthonormality: integral psi_2^2 == 1.
	nodes, weights := orthopolyGaussHermiteTestNodes(12)
	var s float64
	for i, x := range nodes {
		p := HermiteFunction(2, x)
		s += weights[i] * math.Exp(x*x) * p * p
	}
	if !orthopolyClassicalClose(s, 1, 1e-9) {
		t.Errorf("HermiteFunction norm = %.15g, want 1", s)
	}
	if got := HermiteHWeight(0); !orthopolyClassicalClose(got, 1, 1e-15) {
		t.Errorf("HermiteHWeight(0) = %.15g, want 1", got)
	}
}

// orthopolyGaussHermiteTestNodes exposes the Hermite nodes and weights for use
// in the orthonormality check above.
func orthopolyGaussHermiteTestNodes(n int) ([]float64, []float64) {
	return HermiteHRoots(n), HermiteHWeights(n)
}

func TestHermiteProbabilistsClassical(t *testing.T) {
	cases := []struct {
		name string
		got  float64
		want float64
	}{
		{"He0", HermiteHe(0, 0.5), 1},
		{"He1", HermiteHe(1, 0.5), 0.5},
		{"He2", HermiteHe(2, 0.5), -0.75},  // x^2-1
		{"He3", HermiteHe(3, 0.5), -1.375}, // x^3-3x
		{"He4", HermiteHe(4, 1), -2},       // x^4-6x^2+3
		{"dHe3", HermiteHeDerivative(3, 0.5), 3 * HermiteHe(2, 0.5)},
		{"norm3", HermiteHeNorm(3), math.Sqrt(2*math.Pi) * 6},
		{"weightHe", HermiteHeWeight(0), 1},
	}
	for _, c := range cases {
		if !orthopolyClassicalClose(c.got, c.want, 1e-11) {
			t.Errorf("%s = %.15g, want %.15g", c.name, c.got, c.want)
		}
	}
	// He_n(x) = 2^{-n/2} H_n(x/sqrt2).
	x := 0.7
	for n := 0; n <= 6; n++ {
		want := math.Pow(2, -float64(n)/2) * HermiteH(n, x/math.Sqrt2)
		if !orthopolyClassicalClose(HermiteHe(n, x), want, 1e-11) {
			t.Errorf("He_%d(%.2f) relation mismatch", n, x)
		}
	}
	r := HermiteHeRoots(2)
	if !orthopolyClassicalSlicesClose(r, []float64{-1, 1}, 1e-12) {
		t.Errorf("HermiteHeRoots(2) = %v", r)
	}
	if got := NormalizedHermiteHe(0, 0); !orthopolyClassicalClose(got, 1/math.Sqrt(math.Sqrt(2*math.Pi)), 1e-12) {
		t.Errorf("NormalizedHermiteHe(0,0) = %.15g", got)
	}
}

func TestLaguerreClassical(t *testing.T) {
	cases := []struct {
		name string
		got  float64
		want float64
	}{
		{"L0", LaguerreL(0, 2), 1},
		{"L1", LaguerreL(1, 2), -1}, // 1-x
		{"L2(1)", LaguerreL(2, 1), -0.5},
		{"L3(0)", LaguerreL(3, 0), 1},
		{"dL2", LaguerreLDerivative(2, 1), -1}, // -L_1^(1)(1) = -(2-1)
		{"norm", LaguerreNorm(4), 1},
		{"weight", LaguerreWeightFunction(0), 1},
	}
	for _, c := range cases {
		if !orthopolyClassicalClose(c.got, c.want, 1e-11) {
			t.Errorf("%s = %.15g, want %.15g", c.name, c.got, c.want)
		}
	}
	r := LaguerreRoots(2)
	want := []float64{2 - math.Sqrt2, 2 + math.Sqrt2}
	if !orthopolyClassicalSlicesClose(r, want, 1e-11) {
		t.Errorf("LaguerreRoots(2) = %v, want %v", r, want)
	}
	// Gauss-Laguerre integrates x^3 e^{-x} = 3! = 6.
	rule := NewGaussLaguerreRule(6)
	if got := rule.Integrate(func(x float64) float64 { return x * x * x }); !orthopolyClassicalClose(got, 6, 1e-9) {
		t.Errorf("Gauss-Laguerre x^3 = %.15g, want 6", got)
	}
	var sum float64
	for _, w := range LaguerreWeights(8) {
		sum += w
	}
	if !orthopolyClassicalClose(sum, 1, 1e-12) {
		t.Errorf("sum of Gauss-Laguerre weights = %.15g, want 1", sum)
	}
	if got := LaguerreLValues(2, 1); !orthopolyClassicalSlicesClose(got, []float64{1, 0, -0.5}, 1e-12) {
		t.Errorf("LaguerreLValues(2,1) = %v", got)
	}
}

func TestGeneralizedLaguerreClassical(t *testing.T) {
	// L_2^(1)(x) = (x^2 - 6x + 6)/2.
	f := func(x float64) float64 { return (x*x - 6*x + 6) / 2 }
	for _, x := range []float64{0, 1, 2.5, 4} {
		if !orthopolyClassicalClose(GeneralizedLaguerreL(2, 1, x), f(x), 1e-11) {
			t.Errorf("L_2^(1)(%.2f) = %.15g, want %.15g", x, GeneralizedLaguerreL(2, 1, x), f(x))
		}
	}
	for n := 0; n <= 6; n++ {
		if !orthopolyClassicalClose(GeneralizedLaguerreL(n, 0, 1.3), LaguerreL(n, 1.3), 1e-12) {
			t.Errorf("generalized/ordinary mismatch at n=%d", n)
		}
	}
	// d/dx L_2^(1) = x-3, at x=2 gives -1.
	if got := GeneralizedLaguerreLDerivative(2, 1, 2); !orthopolyClassicalClose(got, -1, 1e-11) {
		t.Errorf("dL_2^(1)(2) = %.15g, want -1", got)
	}
	// Norm L_2^(1) = Gamma(4)/2! = 3.
	if got := GeneralizedLaguerreNorm(2, 1); !orthopolyClassicalClose(got, 3, 1e-11) {
		t.Errorf("GeneralizedLaguerreNorm(2,1) = %.15g, want 3", got)
	}
	// Weight function.
	if got := GeneralizedLaguerreWeight(2, 1); !orthopolyClassicalClose(got, math.Exp(-1), 1e-12) {
		t.Errorf("GeneralizedLaguerreWeight(2,1) = %.15g", got)
	}
	// Generalized Gauss-Laguerre weights sum to Gamma(alpha+1).
	alpha := 1.5
	var sum float64
	for _, w := range GeneralizedLaguerreWeights(8, alpha) {
		sum += w
	}
	if !orthopolyClassicalClose(sum, math.Gamma(alpha+1), 1e-9) {
		t.Errorf("generalized Laguerre weight sum = %.15g, want %.15g", sum, math.Gamma(alpha+1))
	}
	for _, x := range GeneralizedLaguerreRoots(7, alpha) {
		if !orthopolyClassicalClose(GeneralizedLaguerreL(7, alpha, x), 0, 1e-8) {
			t.Errorf("L_7^(1.5)(%.15g) not a root", x)
		}
	}
	// Rule-based integration: integral x^alpha e^{-x} over [0,inf) = Gamma(alpha+1).
	rule := NewGeneralizedLaguerreRule(6, alpha)
	if got := rule.Integrate(func(x float64) float64 { return 1 }); !orthopolyClassicalClose(got, math.Gamma(alpha+1), 1e-9) {
		t.Errorf("generalized Laguerre rule mass = %.15g, want %.15g", got, math.Gamma(alpha+1))
	}
}

func TestGegenbauerJacobiClassical(t *testing.T) {
	x := 0.4
	for n := 0; n <= 6; n++ {
		if !orthopolyClassicalClose(GegenbauerC(n, 0.5, x), LegendreP(n, x), 1e-11) {
			t.Errorf("C_%d^(1/2) != P_%d", n, n)
		}
		if !orthopolyClassicalClose(GegenbauerC(n, 1, x), ChebyshevSecond(n, x), 1e-11) {
			t.Errorf("C_%d^(1) != U_%d", n, n)
		}
		if !orthopolyClassicalClose(JacobiP(n, 0, 0, x), LegendreP(n, x), 1e-11) {
			t.Errorf("P_%d^(0,0) != Legendre", n)
		}
	}
	// C_2^(alpha) = 2a(a+1)x^2 - a; derivative 4a(a+1)x.
	a := 1.5
	if got := GegenbauerCDerivative(2, a, x); !orthopolyClassicalClose(got, 4*a*(a+1)*x, 1e-11) {
		t.Errorf("GegenbauerCDerivative(2,%.1f,%.1f) = %.15g", a, x, got)
	}
	// Jacobi P_1^(a,b) = (a-b)/2 + (a+b+2)/2 x.
	if got := JacobiP(1, 2, 1, x); !orthopolyClassicalClose(got, 0.5+2.5*x, 1e-12) {
		t.Errorf("JacobiP(1,2,1,x) = %.15g", got)
	}
	// Jacobi derivative: d/dx P_2^(1,1) = 5/2 P_1^(2,2) = 5/2 * 3x = 7.5x.
	if got := JacobiPDerivative(2, 1, 1, x); !orthopolyClassicalClose(got, 7.5*x, 1e-11) {
		t.Errorf("JacobiPDerivative(2,1,1,x) = %.15g, want %.15g", got, 7.5*x)
	}
	gv := GegenbauerCValues(5, 1.2, x)
	for n := 0; n <= 5; n++ {
		if !orthopolyClassicalClose(gv[n], GegenbauerC(n, 1.2, x), 1e-12) {
			t.Errorf("GegenbauerCValues mismatch at %d", n)
		}
	}
	jv := JacobiPValues(5, 1, 2, x)
	for n := 0; n <= 5; n++ {
		if !orthopolyClassicalClose(jv[n], JacobiP(n, 1, 2, x), 1e-12) {
			t.Errorf("JacobiPValues mismatch at %d", n)
		}
	}
}

func TestQuadratureRuleClassical(t *testing.T) {
	rule := NewGaussLegendreRule(4)
	if rule.Order() != 4 {
		t.Errorf("Order = %d, want 4", rule.Order())
	}
	// Exact for degree <= 2n-1 = 7: integral of x^6 over [-1,1] = 2/7.
	if got := rule.Integrate(func(x float64) float64 { return math.Pow(x, 6) }); !orthopolyClassicalClose(got, 2.0/7.0, 1e-12) {
		t.Errorf("Gauss-Legendre x^6 = %.15g, want %.15g", got, 2.0/7.0)
	}
}

func BenchmarkGeneralizedLaguerreRule64(b *testing.B) {
	for i := 0; i < b.N; i++ {
		_ = NewGeneralizedLaguerreRule(64, 0.5)
	}
}
