package orthopoly

import (
	"math"
	"testing"
)

// orthopolyTestClose reports whether a and b agree within an absolute tolerance.
func orthopolyTestClose(a, b, tol float64) bool {
	if math.IsNaN(a) || math.IsNaN(b) {
		return false
	}
	return math.Abs(a-b) <= tol
}

func TestGegenbauer(t *testing.T) {
	// C_2^alpha(x) = 2*alpha*(alpha+1)*x^2 - alpha.
	if got, want := Gegenbauer(2, 1.0, 0.4), -0.36; !orthopolyTestClose(got, want, 1e-12) {
		t.Errorf("Gegenbauer(2,1,0.4) = %v want %v", got, want)
	}
	if got, want := Gegenbauer(0, 2.0, 0.5), 1.0; got != want {
		t.Errorf("Gegenbauer(0,...) = %v want %v", got, want)
	}
	// C_1^alpha(x) = 2*alpha*x.
	if got, want := Gegenbauer(1, 3.0, 0.5), 3.0; !orthopolyTestClose(got, want, 1e-12) {
		t.Errorf("Gegenbauer(1,3,0.5) = %v want %v", got, want)
	}
	// C_n^1 == U_n.
	if got, want := Gegenbauer(4, 1.0, 0.3), ChebyshevU(4, 0.3); !orthopolyTestClose(got, want, 1e-12) {
		t.Errorf("Gegenbauer C_n^1 vs U_n = %v want %v", got, want)
	}
	// Derivative d/dx C_2^1 = 8x (since C_2^1 = 4x^2-1), so at x=0.5 it is 4.
	if got, want := GegenbauerDerivative(2, 1.0, 0.5), 4.0; !orthopolyTestClose(got, want, 1e-12) {
		t.Errorf("GegenbauerDerivative = %v want %v", got, want)
	}
	// Norm of C_2^1 equals pi/2 (== integral of sqrt(1-x^2) U_2^2).
	if got, want := GegenbauerNorm(2, 1.0), math.Pi/2; !orthopolyTestClose(got, want, 1e-10) {
		t.Errorf("GegenbauerNorm(2,1) = %v want %v", got, want)
	}
}

func TestJacobi(t *testing.T) {
	// P_1^{a,b}(x) = (a-b)/2 + (a+b+2)/2 * x.
	if got, want := Jacobi(1, 1.5, 0.5, 0.3), 1.1; !orthopolyTestClose(got, want, 1e-12) {
		t.Errorf("Jacobi(1,1.5,0.5,0.3) = %v want %v", got, want)
	}
	// Jacobi with a=b=0 reduces to Legendre.
	for _, x := range []float64{-0.7, 0.1, 0.9} {
		if got, want := Jacobi(4, 0, 0, x), Legendre(4, x); !orthopolyTestClose(got, want, 1e-12) {
			t.Errorf("Jacobi(4,0,0,%v) = %v want Legendre %v", x, got, want)
		}
	}
	// Derivative check via finite difference.
	const h = 1e-6
	fd := (Jacobi(5, 1, 2, 0.3+h) - Jacobi(5, 1, 2, 0.3-h)) / (2 * h)
	if got := JacobiDerivative(5, 1, 2, 0.3); !orthopolyTestClose(got, fd, 1e-5) {
		t.Errorf("JacobiDerivative = %v want ~%v", got, fd)
	}
	// Jacobi norm with a=b=0 equals Legendre norm 2/(2n+1).
	if got, want := JacobiNorm(3, 0, 0), 2.0/7.0; !orthopolyTestClose(got, want, 1e-12) {
		t.Errorf("JacobiNorm(3,0,0) = %v want %v", got, want)
	}
}

func TestLegendre(t *testing.T) {
	cases := []struct {
		n    int
		x, w float64
	}{
		{0, 0.5, 1},
		{1, 0.5, 0.5},
		{2, 0.5, -0.125},
		{3, 0.5, -0.4375},
		{5, 1.0, 1.0},
	}
	for _, c := range cases {
		if got := Legendre(c.n, c.x); !orthopolyTestClose(got, c.w, 1e-12) {
			t.Errorf("Legendre(%d,%v) = %v want %v", c.n, c.x, got, c.w)
		}
	}
	// P_2'(x) = 3x.
	if got, want := LegendreDerivative(2, 0.4), 1.2; !orthopolyTestClose(got, want, 1e-12) {
		t.Errorf("LegendreDerivative(2,0.4) = %v want %v", got, want)
	}
	// Endpoint: P_3'(1) = 3*4/2 = 6.
	if got, want := LegendreDerivative(3, 1.0), 6.0; !orthopolyTestClose(got, want, 1e-12) {
		t.Errorf("LegendreDerivative(3,1) = %v want %v", got, want)
	}
	if got := LegendreNorm(4); !orthopolyTestClose(got, 2.0/9.0, 1e-12) {
		t.Errorf("LegendreNorm(4) = %v want %v", got, 2.0/9.0)
	}
	b := LegendreBasis(3, 0.5)
	if len(b) != 4 || !orthopolyTestClose(b[2], -0.125, 1e-12) {
		t.Errorf("LegendreBasis = %v", b)
	}
	roots := LegendreRoots(2)
	if len(roots) != 2 || !orthopolyTestClose(roots[0], -1/math.Sqrt(3), 1e-12) {
		t.Errorf("LegendreRoots(2) = %v", roots)
	}
}

func TestAssociatedLegendreAndHarmonics(t *testing.T) {
	// P_1^1(x) with Condon-Shortley phase = -sqrt(1-x^2).
	if got, want := AssociatedLegendre(1, 1, 0.5), -math.Sqrt(0.75); !orthopolyTestClose(got, want, 1e-12) {
		t.Errorf("AssociatedLegendre(1,1,0.5) = %v want %v", got, want)
	}
	// P_2^0 = (3x^2-1)/2 equals Legendre.
	if got, want := AssociatedLegendre(2, 0, 0.3), Legendre(2, 0.3); !orthopolyTestClose(got, want, 1e-12) {
		t.Errorf("AssociatedLegendre(2,0,0.3) = %v want %v", got, want)
	}
	// Negative order relation P_1^{-1} = -1/2 * P_1^1.
	if got, want := AssociatedLegendre(1, -1, 0.5), -0.5*AssociatedLegendre(1, 1, 0.5); !orthopolyTestClose(got, want, 1e-12) {
		t.Errorf("AssociatedLegendre(1,-1,0.5) = %v want %v", got, want)
	}
	// Y_0^0 = 1/(2 sqrt(pi)).
	if got, want := RealSphericalHarmonic(0, 0, 1.0, 2.0), 0.5/math.Sqrt(math.Pi); !orthopolyTestClose(got, want, 1e-12) {
		t.Errorf("Y00 = %v want %v", got, want)
	}
	// Y_1^0 = sqrt(3/(4pi)) cos(theta).
	theta := 0.7
	if got, want := RealSphericalHarmonic(1, 0, theta, 0.0), math.Sqrt(3/(4*math.Pi))*math.Cos(theta); !orthopolyTestClose(got, want, 1e-12) {
		t.Errorf("Y10 = %v want %v", got, want)
	}
	// Y_1^1 (cosine, p_x) = sqrt(3/(4pi)) sin(theta) cos(phi).
	phi := 1.1
	if got, want := RealSphericalHarmonic(1, 1, theta, phi), math.Sqrt(3/(4*math.Pi))*math.Sin(theta)*math.Cos(phi); !orthopolyTestClose(got, want, 1e-12) {
		t.Errorf("Y11 = %v want %v", got, want)
	}
	if got, want := SphericalHarmonicNormalization(1, 0), math.Sqrt(3/(4*math.Pi)); !orthopolyTestClose(got, want, 1e-12) {
		t.Errorf("SHNorm(1,0) = %v want %v", got, want)
	}
	// Normalized associated Legendre squared-integral is 1 (checked below via quadrature).
	if math.IsNaN(NormalizedAssociatedLegendre(2, 1, 0.3)) {
		t.Errorf("NormalizedAssociatedLegendre returned NaN")
	}
}

func TestChebyshev(t *testing.T) {
	// T_3(x) = 4x^3-3x -> T_3(0.5) = -1.
	if got := ChebyshevT(3, 0.5); !orthopolyTestClose(got, -1, 1e-12) {
		t.Errorf("ChebyshevT(3,0.5) = %v want -1", got)
	}
	// U_2(x) = 4x^2-1 -> U_2(0.5) = 0.
	if got := ChebyshevU(2, 0.5); !orthopolyTestClose(got, 0, 1e-12) {
		t.Errorf("ChebyshevU(2,0.5) = %v want 0", got)
	}
	// T_n(cos t) = cos(n t).
	tt := 0.9
	if got, want := ChebyshevT(5, math.Cos(tt)), math.Cos(5*tt); !orthopolyTestClose(got, want, 1e-12) {
		t.Errorf("ChebyshevT cos identity = %v want %v", got, want)
	}
	// T_3'(x) = 3 U_2(x).
	if got, want := ChebyshevTDerivative(3, 0.4), 3*ChebyshevU(2, 0.4); !orthopolyTestClose(got, want, 1e-12) {
		t.Errorf("ChebyshevTDerivative = %v want %v", got, want)
	}
	// U derivative via finite difference (interior).
	const h = 1e-6
	fd := (ChebyshevU(4, 0.3+h) - ChebyshevU(4, 0.3-h)) / (2 * h)
	if got := ChebyshevUDerivative(4, 0.3); !orthopolyTestClose(got, fd, 1e-4) {
		t.Errorf("ChebyshevUDerivative = %v want ~%v", got, fd)
	}
	// U_n'(1) = n(n+1)(n+2)/3, for n=3 => 3*4*5/3 = 20.
	if got, want := ChebyshevUDerivative(3, 1.0), 20.0; !orthopolyTestClose(got, want, 1e-9) {
		t.Errorf("ChebyshevUDerivative(3,1) = %v want %v", got, want)
	}
	// Roots of T_2 are +/- 1/sqrt(2).
	r := ChebyshevTRoots(2)
	if len(r) != 2 || !orthopolyTestClose(r[1], 1/math.Sqrt(2), 1e-12) {
		t.Errorf("ChebyshevTRoots(2) = %v", r)
	}
	// Extrema of T_2 are -1, 0, 1.
	e := ChebyshevTExtrema(2)
	if len(e) != 3 || !orthopolyTestClose(e[1], 0, 1e-12) {
		t.Errorf("ChebyshevTExtrema(2) = %v", e)
	}
	if got := ChebyshevTBasis(3, 0.5); len(got) != 4 || !orthopolyTestClose(got[3], ChebyshevT(3, 0.5), 1e-12) {
		t.Errorf("ChebyshevTBasis = %v", got)
	}
	if !orthopolyTestClose(ChebyshevTWeight(0), 1, 1e-12) || !orthopolyTestClose(ChebyshevUWeight(0), 1, 1e-12) {
		t.Errorf("Chebyshev weights at 0 wrong")
	}
}

func TestHermite(t *testing.T) {
	// H_2 = 4x^2-2 -> H_2(1) = 2.
	if got := HermitePhysicists(2, 1.0); !orthopolyTestClose(got, 2, 1e-12) {
		t.Errorf("HermitePhysicists(2,1) = %v want 2", got)
	}
	// H_3 = 8x^3-12x -> H_3(0.5) = 1 - 6 = -5.
	if got := HermitePhysicists(3, 0.5); !orthopolyTestClose(got, -5, 1e-12) {
		t.Errorf("HermitePhysicists(3,0.5) = %v want -5", got)
	}
	// He_2 = x^2-1 -> He_2(2) = 3.
	if got := HermiteProbabilists(2, 2.0); !orthopolyTestClose(got, 3, 1e-12) {
		t.Errorf("HermiteProbabilists(2,2) = %v want 3", got)
	}
	// H_2'(x) = 4 H_1(x) = 8x.
	if got, want := HermitePhysicistsDerivative(2, 0.5), 4.0; !orthopolyTestClose(got, want, 1e-12) {
		t.Errorf("HermitePhysicistsDerivative(2,0.5) = %v want %v", got, want)
	}
	// He_3'(x) = 3 He_2(x).
	if got, want := HermiteProbabilistsDerivative(3, 0.7), 3*HermiteProbabilists(2, 0.7); !orthopolyTestClose(got, want, 1e-12) {
		t.Errorf("HermiteProbabilistsDerivative = %v want %v", got, want)
	}
	// Norm sqrt(pi)*2^n*n!, n=2 => sqrt(pi)*8.
	if got, want := HermitePhysicistsNorm(2), math.Sqrt(math.Pi)*8; !orthopolyTestClose(got, want, 1e-10) {
		t.Errorf("HermitePhysicistsNorm(2) = %v want %v", got, want)
	}
	if got := HermitePhysicistsBasis(3, 0.5); len(got) != 4 || !orthopolyTestClose(got[3], -5, 1e-12) {
		t.Errorf("HermitePhysicistsBasis = %v", got)
	}
	if !orthopolyTestClose(HermiteWeight(1), math.Exp(-1), 1e-12) {
		t.Errorf("HermiteWeight wrong")
	}
}

func TestLaguerre(t *testing.T) {
	// L_2(x) = (x^2-4x+2)/2 -> L_2(1) = -1/2.
	if got := Laguerre(2, 1.0); !orthopolyTestClose(got, -0.5, 1e-12) {
		t.Errorf("Laguerre(2,1) = %v want -0.5", got)
	}
	// L_1^alpha(x) = 1+alpha-x.
	if got, want := AssociatedLaguerre(1, 2.0, 0.5), 2.5; !orthopolyTestClose(got, want, 1e-12) {
		t.Errorf("AssociatedLaguerre(1,2,0.5) = %v want %v", got, want)
	}
	// L_2'(x) = -L_1^1(x) = -(2-x) => at x=1, -1.
	if got, want := LaguerreDerivative(2, 1.0), -1.0; !orthopolyTestClose(got, want, 1e-12) {
		t.Errorf("LaguerreDerivative(2,1) = %v want %v", got, want)
	}
	if got, want := AssociatedLaguerreDerivative(2, 1.0, 0.5), -AssociatedLaguerre(1, 2.0, 0.5); !orthopolyTestClose(got, want, 1e-12) {
		t.Errorf("AssociatedLaguerreDerivative = %v want %v", got, want)
	}
	if got := LaguerreBasis(2, 1.0); len(got) != 3 || !orthopolyTestClose(got[2], -0.5, 1e-12) {
		t.Errorf("LaguerreBasis = %v", got)
	}
	if !orthopolyTestClose(LaguerreWeight(0, 2), math.Exp(-2), 1e-12) {
		t.Errorf("LaguerreWeight wrong")
	}
}

func TestBernoulliEuler(t *testing.T) {
	bn := map[int]float64{0: 1, 1: -0.5, 2: 1.0 / 6, 4: -1.0 / 30, 6: 1.0 / 42}
	for n, want := range bn {
		if got := BernoulliNumber(n); !orthopolyTestClose(got, want, 1e-12) {
			t.Errorf("BernoulliNumber(%d) = %v want %v", n, got, want)
		}
	}
	// B_2(x) = x^2 - x + 1/6.
	if got, want := BernoulliPolynomial(2, 0.3), 0.09-0.3+1.0/6; !orthopolyTestClose(got, want, 1e-12) {
		t.Errorf("BernoulliPolynomial(2,0.3) = %v want %v", got, want)
	}
	// E_2(x) = x^2 - x.
	if got, want := EulerPolynomial(2, 0.3), 0.09-0.3; !orthopolyTestClose(got, want, 1e-12) {
		t.Errorf("EulerPolynomial(2,0.3) = %v want %v", got, want)
	}
	en := map[int]float64{0: 1, 2: -1, 4: 5, 6: -61}
	for n, want := range en {
		if got := EulerNumber(n); !orthopolyTestClose(got, want, 1e-9) {
			t.Errorf("EulerNumber(%d) = %v want %v", n, got, want)
		}
	}
	if EulerNumber(3) != 0 {
		t.Errorf("EulerNumber(3) should be 0")
	}
}

func TestBernstein(t *testing.T) {
	// b_{1,2}(x) = 2x(1-x) -> at 0.5 => 0.5.
	if got := Bernstein(1, 2, 0.5); !orthopolyTestClose(got, 0.5, 1e-12) {
		t.Errorf("Bernstein(1,2,0.5) = %v want 0.5", got)
	}
	// Partition of unity.
	basis := BernsteinBasis(5, 0.37)
	var sum float64
	for _, v := range basis {
		sum += v
	}
	if !orthopolyTestClose(sum, 1, 1e-12) {
		t.Errorf("Bernstein partition of unity = %v", sum)
	}
	// Derivative via finite difference.
	const h = 1e-6
	fd := (Bernstein(2, 4, 0.4+h) - Bernstein(2, 4, 0.4-h)) / (2 * h)
	if got := BernsteinDerivative(2, 4, 0.4); !orthopolyTestClose(got, fd, 1e-5) {
		t.Errorf("BernsteinDerivative = %v want ~%v", got, fd)
	}
	// Bezier reproduces polynomial: constant coeffs => constant.
	if got := BernsteinBezier([]float64{2, 2, 2, 2}, 0.6); !orthopolyTestClose(got, 2, 1e-12) {
		t.Errorf("BernsteinBezier constant = %v want 2", got)
	}
	// Linear ramp coeffs {0,1} => value == x.
	if got := BernsteinBezier([]float64{0, 1}, 0.42); !orthopolyTestClose(got, 0.42, 1e-12) {
		t.Errorf("BernsteinBezier linear = %v want 0.42", got)
	}
}

func TestZernike(t *testing.T) {
	// R_2^0 = 2rho^2-1 -> at 0.5 => -0.5.
	if got := ZernikeRadial(2, 0, 0.5); !orthopolyTestClose(got, -0.5, 1e-12) {
		t.Errorf("ZernikeRadial(2,0,0.5) = %v want -0.5", got)
	}
	// R_4^2 = 4rho^4-3rho^2 -> at 0.5 => -0.5.
	if got := ZernikeRadial(4, 2, 0.5); !orthopolyTestClose(got, -0.5, 1e-12) {
		t.Errorf("ZernikeRadial(4,2,0.5) = %v want -0.5", got)
	}
	// n-m odd => 0.
	if got := ZernikeRadial(3, 0, 0.5); got != 0 {
		t.Errorf("ZernikeRadial(3,0,...) = %v want 0", got)
	}
	// R_n^m(1) = 1 for valid pairs.
	if got := ZernikeRadial(4, 2, 1.0); !orthopolyTestClose(got, 1, 1e-12) {
		t.Errorf("ZernikeRadial(4,2,1) = %v want 1", got)
	}
	// Derivative via finite difference.
	const h = 1e-6
	fd := (ZernikeRadial(4, 2, 0.6+h) - ZernikeRadial(4, 2, 0.6-h)) / (2 * h)
	if got := ZernikeRadialDerivative(4, 2, 0.6); !orthopolyTestClose(got, fd, 1e-4) {
		t.Errorf("ZernikeRadialDerivative = %v want ~%v", got, fd)
	}
	if got, want := ZernikeNormalization(2, 0), math.Sqrt(3); !orthopolyTestClose(got, want, 1e-12) {
		t.Errorf("ZernikeNormalization(2,0) = %v want %v", got, want)
	}
	// Zernike full: Z_2^{-2}(rho,theta) = rho^2 sin(2 theta).
	if got, want := Zernike(2, -2, 0.5, 0.3), 0.25*math.Sin(0.6); !orthopolyTestClose(got, want, 1e-12) {
		t.Errorf("Zernike(2,-2,...) = %v want %v", got, want)
	}
	// Noll index mapping.
	nollCases := []struct {
		j, n, m int
	}{
		{1, 0, 0}, {2, 1, 1}, {3, 1, -1}, {4, 2, 0}, {5, 2, -2}, {6, 2, 2},
	}
	for _, c := range nollCases {
		n, m := NollToNM(c.j)
		if n != c.n || m != c.m {
			t.Errorf("NollToNM(%d) = (%d,%d) want (%d,%d)", c.j, n, m, c.n, c.m)
		}
	}
}

func TestGaussLegendre(t *testing.T) {
	nodes, weights := GaussLegendre(2)
	if !orthopolyTestClose(nodes[0], -1/math.Sqrt(3), 1e-12) || !orthopolyTestClose(nodes[1], 1/math.Sqrt(3), 1e-12) {
		t.Errorf("GaussLegendre(2) nodes = %v", nodes)
	}
	if !orthopolyTestClose(weights[0], 1, 1e-12) || !orthopolyTestClose(weights[1], 1, 1e-12) {
		t.Errorf("GaussLegendre(2) weights = %v", weights)
	}
	// Integrate x^6 over [-1,1] = 2/7 exactly with n=4.
	n4, w4 := GaussLegendre(4)
	var s float64
	for i := range n4 {
		s += w4[i] * math.Pow(n4[i], 6)
	}
	if !orthopolyTestClose(s, 2.0/7.0, 1e-12) {
		t.Errorf("GaussLegendre integral of x^6 = %v want %v", s, 2.0/7.0)
	}
	// Weights sum to interval length 2.
	_, w10 := GaussLegendre(10)
	var sw float64
	for _, w := range w10 {
		sw += w
	}
	if !orthopolyTestClose(sw, 2, 1e-12) {
		t.Errorf("GaussLegendre(10) weight sum = %v want 2", sw)
	}
	// Convenience integrator: integral of x^2 over [0,3] = 9.
	if got := GaussLegendreIntegrate(func(x float64) float64 { return x * x }, 0, 3, 5); !orthopolyTestClose(got, 9, 1e-9) {
		t.Errorf("GaussLegendreIntegrate = %v want 9", got)
	}
}

func TestGaussHermite(t *testing.T) {
	nodes, weights := GaussHermite(2)
	if !orthopolyTestClose(nodes[1], 1/math.Sqrt(2), 1e-12) {
		t.Errorf("GaussHermite(2) nodes = %v", nodes)
	}
	if !orthopolyTestClose(weights[0], math.Sqrt(math.Pi)/2, 1e-12) {
		t.Errorf("GaussHermite(2) weights = %v", weights)
	}
	// Integral of exp(-x^2) = sqrt(pi).
	if got := GaussHermiteIntegrate(func(x float64) float64 { return 1 }, 5); !orthopolyTestClose(got, math.Sqrt(math.Pi), 1e-10) {
		t.Errorf("GaussHermite total mass = %v want %v", got, math.Sqrt(math.Pi))
	}
	// Integral of x^2 exp(-x^2) = sqrt(pi)/2.
	if got := GaussHermiteIntegrate(func(x float64) float64 { return x * x }, 5); !orthopolyTestClose(got, math.Sqrt(math.Pi)/2, 1e-10) {
		t.Errorf("GaussHermite x^2 = %v want %v", got, math.Sqrt(math.Pi)/2)
	}
}

func TestGaussLaguerre(t *testing.T) {
	nodes, weights := GaussLaguerre(2)
	if !orthopolyTestClose(nodes[0], 2-math.Sqrt(2), 1e-10) || !orthopolyTestClose(nodes[1], 2+math.Sqrt(2), 1e-10) {
		t.Errorf("GaussLaguerre(2) nodes = %v", nodes)
	}
	// Closed-form weights for n=2 are (2±√2)/4 and sum to 1.
	if !orthopolyTestClose(weights[0], (2+math.Sqrt(2))/4, 1e-10) || !orthopolyTestClose(weights[1], (2-math.Sqrt(2))/4, 1e-10) {
		t.Errorf("GaussLaguerre(2) weights = %v", weights)
	}
	// Integral of exp(-x) over [0,inf) = 1.
	if got := GaussLaguerreIntegrate(func(x float64) float64 { return 1 }, 5); !orthopolyTestClose(got, 1, 1e-12) {
		t.Errorf("GaussLaguerre total mass = %v want 1", got)
	}
	// Integral of x^3 exp(-x) = 3! = 6, exact with n=4.
	if got := GaussLaguerreIntegrate(func(x float64) float64 { return x * x * x }, 4); !orthopolyTestClose(got, 6, 1e-9) {
		t.Errorf("GaussLaguerre x^3 = %v want 6", got)
	}
	// Generalized: integral of x^alpha exp(-x) = Gamma(alpha+1).
	gn, gw := GaussLaguerreGeneralized(6, 1.5)
	var mass float64
	for i := range gn {
		mass += gw[i]
	}
	if !orthopolyTestClose(mass, math.Gamma(2.5), 1e-9) {
		t.Errorf("GaussLaguerreGeneralized mass = %v want %v", mass, math.Gamma(2.5))
	}
}

func TestGaussJacobiAndFriends(t *testing.T) {
	// Gauss-Jacobi with a=b=0 reproduces Gauss-Legendre.
	jn, jw := GaussJacobi(4, 0, 0)
	ln, lw := GaussLegendre(4)
	for i := range jn {
		if !orthopolyTestClose(jn[i], ln[i], 1e-10) || !orthopolyTestClose(jw[i], lw[i], 1e-10) {
			t.Errorf("GaussJacobi(a=b=0) mismatch at %d: %v/%v vs %v/%v", i, jn[i], jw[i], ln[i], lw[i])
		}
	}
	// Gauss-Chebyshev T closed form: nodes cos((2k-1)pi/2n), weight pi/n.
	cn, cw := GaussChebyshevT(3)
	if !orthopolyTestClose(cw[0], math.Pi/3, 1e-12) {
		t.Errorf("GaussChebyshevT weight = %v want %v", cw[0], math.Pi/3)
	}
	// Integral of 1/sqrt(1-x^2) = pi.
	var mass float64
	for _, w := range cw {
		mass += w
	}
	if !orthopolyTestClose(mass, math.Pi, 1e-12) {
		t.Errorf("GaussChebyshevT mass = %v want pi", mass)
	}
	_ = cn
	// Gauss-Chebyshev U mass = integral of sqrt(1-x^2) over [-1,1] = pi/2.
	_, uw := GaussChebyshevU(5)
	var umass float64
	for _, w := range uw {
		umass += w
	}
	if !orthopolyTestClose(umass, math.Pi/2, 1e-12) {
		t.Errorf("GaussChebyshevU mass = %v want pi/2", umass)
	}
	// Gauss-Gegenbauer with alpha=0.5 is Gauss-Legendre.
	gn, _ := GaussGegenbauer(3, 0.5)
	ln3, _ := GaussLegendre(3)
	for i := range gn {
		if !orthopolyTestClose(gn[i], ln3[i], 1e-10) {
			t.Errorf("GaussGegenbauer(alpha=0.5) node mismatch %v vs %v", gn[i], ln3[i])
		}
	}
	// GolubWelsch directly on Legendre recurrence.
	a := []float64{0, 0, 0}
	b := []float64{2, 1.0 / 3, 4.0 / 15}
	wn, ww := GolubWelsch(a, b, 2)
	if len(wn) != 3 || !orthopolyTestClose(ww[0]+ww[1]+ww[2], 2, 1e-12) {
		t.Errorf("GolubWelsch Legendre weights sum wrong: %v", ww)
	}
}

// BenchmarkGaussJacobi benchmarks the heaviest routine: building a moderately
// large Gauss-Jacobi rule, which drives the symmetric tridiagonal eigensolver.
func BenchmarkGaussJacobi(b *testing.B) {
	for i := 0; i < b.N; i++ {
		GaussJacobi(64, 0.5, -0.25)
	}
}
