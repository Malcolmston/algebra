package special

import (
	"math"
	"testing"
)

func TestBesselFirstKind(t *testing.T) {
	cases := []struct {
		name string
		got  float64
		want float64
	}{
		{"J0(1)", BesselJ0(1), 0.7651976865579666},
		{"J1(1)", BesselJ1(1), 0.4400505857449335},
		{"Jn(2,1)", BesselJn(2, 1), 0.1149034849319005},
		{"Y0(1)", BesselY0(1), 0.08825696421567697},
		{"Y1(1)", BesselY1(1), -0.7812128213002887},
		{"Yn(2,1)", BesselYn(2, 1), -1.6506826068162548},
	}
	for _, c := range cases {
		if !specialClose(c.got, c.want, 1e-10) {
			t.Errorf("%s = %.15g, want %.15g", c.name, c.got, c.want)
		}
	}
}

func TestModifiedBessel(t *testing.T) {
	cases := []struct {
		name string
		got  float64
		want float64
		tol  float64
	}{
		{"I0(1)", BesselI0(1), 1.2660658777520082, 1e-6},
		{"I1(1)", BesselI1(1), 0.565159103992485, 1e-6},
		{"In(2,1)", BesselIn(2, 1), 0.1357476697670383, 1e-6},
		{"I0(3)", BesselI0(3), 4.880792585865024, 1e-5},
		{"In(3,2)", BesselIn(3, 2), 0.2127399592398527, 1e-6},
		{"K0(1)", BesselK0(1), 0.42102443824070834, 1e-6},
		{"K1(1)", BesselK1(1), 0.6019072301972346, 1e-6},
		{"Kn(2,1)", BesselKn(2, 1), 1.6248388986351774, 1e-5},
		{"K0(3)", BesselK0(3), 0.03473950438627925, 1e-7},
	}
	for _, c := range cases {
		if !specialClose(c.got, c.want, c.tol) {
			t.Errorf("%s = %.15g, want %.15g", c.name, c.got, c.want)
		}
	}
}

func TestModifiedBesselWronskianConsistency(t *testing.T) {
	// I_n(x) K_n'(x) - I_n'(x) K_n(x) = -1/x, checked numerically per order.
	for _, x := range []float64{0.5, 1.0, 2.5, 4.0} {
		for n := 0; n <= 3; n++ {
			got := BesselIn(n, x)*BesselKnPrime(n, x) - BesselInPrime(n, x)*BesselKn(n, x)
			if !specialClose(got, BesselWronskianIK(x), 1e-6) {
				t.Errorf("Wronskian IK n=%d x=%g: got %.12g want %.12g", n, x, got, -1/x)
			}
		}
	}
}

func TestExponentiallyScaled(t *testing.T) {
	for _, x := range []float64{0.3, 1.0, 3.0, 10.0, 40.0} {
		if got, want := BesselI0e(x), BesselI0(x)*math.Exp(-x); math.Abs(x) < 20 && !specialClose(got, want, 1e-6) {
			t.Errorf("I0e(%g) = %.12g, want %.12g", x, got, want)
		}
		if got, want := BesselK0e(x), BesselK0(x)*math.Exp(x); math.Abs(x) < 20 && !specialClose(got, want, 1e-6) {
			t.Errorf("K0e(%g) = %.12g, want %.12g", x, got, want)
		}
	}
	// Scaled forms must stay finite where the unscaled ones overflow.
	if v := BesselI0e(700); math.IsInf(v, 0) || v <= 0 {
		t.Errorf("I0e(700) = %v, want small finite positive", v)
	}
	if v := BesselK0e(700); math.IsInf(v, 0) || v <= 0 {
		t.Errorf("K0e(700) = %v, want small finite positive", v)
	}
}

func TestSphericalBessel(t *testing.T) {
	x := 1.0
	// Closed-form reference values.
	j0 := math.Sin(x) / x
	j1 := math.Sin(x)/(x*x) - math.Cos(x)/x
	j2 := (3/(x*x*x)-1/x)*math.Sin(x) - 3/(x*x)*math.Cos(x)
	y0 := -math.Cos(x) / x
	y1 := -math.Cos(x)/(x*x) - math.Sin(x)/x
	cases := []struct {
		name string
		got  float64
		want float64
	}{
		{"j0(1)", SphericalJ0(1), j0},
		{"j1(1)", SphericalJ1(1), j1},
		{"jn(2,1)", SphericalJn(2, 1), j2},
		{"y0(1)", SphericalY0(1), y0},
		{"y1(1)", SphericalY1(1), y1},
	}
	for _, c := range cases {
		if !specialClose(c.got, c.want, 1e-10) {
			t.Errorf("%s = %.15g, want %.15g", c.name, c.got, c.want)
		}
	}
	// Downward recurrence branch (x < n): compare jn against the recurrence
	// relation j_{n-1}+j_{n+1} = (2n+1)/x j_n.
	for n := 2; n <= 8; n++ {
		lhs := SphericalJn(n-1, 0.5) + SphericalJn(n+1, 0.5)
		rhs := float64(2*n+1) / 0.5 * SphericalJn(n, 0.5)
		if !specialClose(lhs, rhs, 1e-9) {
			t.Errorf("spherical recurrence n=%d: %.12g vs %.12g", n, lhs, rhs)
		}
	}
}

func TestModifiedSpherical(t *testing.T) {
	x := 1.5
	if !specialClose(ModifiedSphericalI0(x), math.Sinh(x)/x, 1e-12) {
		t.Errorf("i0 mismatch")
	}
	if !specialClose(ModifiedSphericalI1(x), (x*math.Cosh(x)-math.Sinh(x))/(x*x), 1e-12) {
		t.Errorf("i1 mismatch")
	}
	if !specialClose(ModifiedSphericalK0(x), math.Pi/2*math.Exp(-x)/x, 1e-12) {
		t.Errorf("k0 mismatch")
	}
	if !specialClose(ModifiedSphericalK1(x), math.Pi/2*math.Exp(-x)*(1+1/x)/x, 1e-12) {
		t.Errorf("k1 mismatch")
	}
}

func TestAiry(t *testing.T) {
	cases := []struct {
		name string
		got  float64
		want float64
		tol  float64
	}{
		{"Ai(0)", AiryAi(0), 0.3550280538878172, 1e-9},
		{"Ai(1)", AiryAi(1), 0.13529241631288147, 1e-8},
		{"Ai(2)", AiryAi(2), 0.03492413042327235, 1e-8},
		{"Ai(-1)", AiryAi(-1), 0.5355608832923521, 1e-8},
		{"Bi(0)", AiryBi(0), 0.6149266274460007, 1e-9},
		{"Bi(1)", AiryBi(1), 1.2074235949528712, 1e-8},
		{"Bi(-1)", AiryBi(-1), 0.10399738949694459, 1e-8},
		{"Ai'(0)", AiryAiPrime(0), -0.2588194037928068, 1e-9},
		{"Ai'(1)", AiryAiPrime(1), -0.1591474412967932, 1e-8},
		{"Bi'(0)", AiryBiPrime(0), 0.4482883573538264, 1e-9},
		{"Bi'(1)", AiryBiPrime(1), 0.9324359333927756, 1e-8},
	}
	for _, c := range cases {
		if !specialClose(c.got, c.want, c.tol) {
			t.Errorf("%s = %.15g, want %.15g", c.name, c.got, c.want)
		}
	}
	// Wronskian identity Ai(x)Bi'(x) - Ai'(x)Bi(x) = 1/pi.
	for _, x := range []float64{-2, -0.5, 0, 0.5, 2} {
		w := AiryAi(x)*AiryBiPrime(x) - AiryAiPrime(x)*AiryBi(x)
		if !specialClose(w, 1/math.Pi, 1e-6) {
			t.Errorf("Airy Wronskian at %g = %.10g, want %.10g", x, w, 1/math.Pi)
		}
	}
}

func TestStruve(t *testing.T) {
	if !specialClose(StruveH0(1), 0.5686566270825707, 1e-5) {
		t.Errorf("H0(1) = %.12g", StruveH0(1))
	}
	if !specialClose(StruveH1(1), 0.19845643737084569, 1e-5) {
		t.Errorf("H1(1) = %.12g", StruveH1(1))
	}
	// Odd/even symmetry.
	if !specialClose(StruveH0(-1), -StruveH0(1), 1e-12) {
		t.Errorf("H0 should be odd")
	}
	if !specialClose(StruveH1(-1), StruveH1(1), 1e-12) {
		t.Errorf("H1 should be even")
	}
}

func TestModifiedStruve(t *testing.T) {
	if !specialClose(StruveL0(1), 0.7102431745671832, 1e-5) {
		t.Errorf("L0(1) = %.12g", StruveL0(1))
	}
	if !specialClose(StruveL1(1), 0.2267643858228386, 1e-5) {
		t.Errorf("L1(1) = %.12g", StruveL1(1))
	}
}

func TestKelvin(t *testing.T) {
	cases := []struct {
		name string
		got  float64
		want float64
	}{
		{"ber(1)", KelvinBer(1), 0.9843817812130869},
		{"bei(1)", KelvinBei(1), 0.2495660400366597},
		{"ker(1)", KelvinKer(1), 0.2867062087283160},
		{"kei(1)", KelvinKei(1), -0.4949946365187199},
	}
	for _, c := range cases {
		if !specialClose(c.got, c.want, 1e-6) {
			t.Errorf("%s = %.15g, want %.15g", c.name, c.got, c.want)
		}
	}
}

func TestBesselRealOrderHalfInteger(t *testing.T) {
	// J_{1/2}(x) = sqrt(2/(pi x)) sin(x); I_{1/2}(x) = sqrt(2/(pi x)) sinh(x).
	for _, x := range []float64{0.25, 1.0, 3.0, 6.0} {
		wantJ := math.Sqrt(2/(math.Pi*x)) * math.Sin(x)
		if !specialClose(BesselJnu(0.5, x), wantJ, 1e-9) {
			t.Errorf("Jnu(0.5,%g) = %.12g, want %.12g", x, BesselJnu(0.5, x), wantJ)
		}
		wantI := math.Sqrt(2/(math.Pi*x)) * math.Sinh(x)
		if !specialClose(BesselInu(0.5, x), wantI, 1e-8) {
			t.Errorf("Inu(0.5,%g) = %.12g, want %.12g", x, BesselInu(0.5, x), wantI)
		}
	}
	// Integer order must agree with the stdlib-backed routines.
	if !specialClose(BesselJnu(0, 2.5), math.J0(2.5), 1e-10) {
		t.Errorf("Jnu(0,2.5) mismatch")
	}
	if !specialClose(BesselJnu(1, 2.5), math.J1(2.5), 1e-10) {
		t.Errorf("Jnu(1,2.5) mismatch")
	}
}

func TestRiccatiBessel(t *testing.T) {
	x := 2.0
	if !specialClose(RiccatiBesselPsi(1, x), x*SphericalJn(1, x), 1e-12) {
		t.Errorf("psi mismatch")
	}
	if !specialClose(RiccatiBesselChi(1, x), -x*SphericalYn(1, x), 1e-12) {
		t.Errorf("chi mismatch")
	}
}

func TestDerivatives(t *testing.T) {
	// Compare analytic derivatives against central finite differences.
	h := 1e-6
	fd := func(f func(float64) float64, x float64) float64 {
		return (f(x+h) - f(x-h)) / (2 * h)
	}
	x := 1.7
	checks := []struct {
		name string
		d    float64
		fd   float64
	}{
		{"J0'", BesselJ0Prime(x), fd(BesselJ0, x)},
		{"J1'", BesselJ1Prime(x), fd(BesselJ1, x)},
		{"Y0'", BesselY0Prime(x), fd(BesselY0, x)},
		{"Y1'", BesselY1Prime(x), fd(BesselY1, x)},
		{"I0'", BesselI0Prime(x), fd(BesselI0, x)},
		{"I1'", BesselI1Prime(x), fd(BesselI1, x)},
		{"K0'", BesselK0Prime(x), fd(BesselK0, x)},
		{"K1'", BesselK1Prime(x), fd(BesselK1, x)},
		{"Jn'", BesselJnPrime(3, x), fd(func(t float64) float64 { return BesselJn(3, t) }, x)},
		{"In'", BesselInPrime(3, x), fd(func(t float64) float64 { return BesselIn(3, t) }, x)},
		{"jn'", SphericalJnPrime(2, x), fd(func(t float64) float64 { return SphericalJn(2, t) }, x)},
		{"yn'", SphericalYnPrime(2, x), fd(func(t float64) float64 { return SphericalYn(2, t) }, x)},
	}
	for _, c := range checks {
		if !specialClose(c.d, c.fd, 1e-5) {
			t.Errorf("%s = %.10g, finite-diff %.10g", c.name, c.d, c.fd)
		}
	}
}

func TestZeros(t *testing.T) {
	// McMahon asymptotic zeros of J0.
	if !specialClose(BesselJZero(0, 5), 14.930917708487786, 1e-4) {
		t.Errorf("J0 5th zero = %.12g", BesselJZero(0, 5))
	}
	if !specialClose(BesselJZero(0, 10), 30.634606468431976, 1e-4) {
		t.Errorf("J0 10th zero = %.12g", BesselJZero(0, 10))
	}
	// Airy zeros (asymptotic).
	if !specialClose(AiryZeroAi(1), -2.338107410459767, 5e-3) {
		t.Errorf("Ai 1st zero = %.12g", AiryZeroAi(1))
	}
	if !specialClose(AiryZeroAi(4), -6.786708090071759, 3e-3) {
		t.Errorf("Ai 4th zero = %.12g", AiryZeroAi(4))
	}
}

func TestWronskians(t *testing.T) {
	x := 3.3
	// J_n Y_n' - J_n' Y_n = 2/(pi x).
	got := BesselJ0(x)*BesselY0Prime(x) - BesselJ0Prime(x)*BesselY0(x)
	if !specialClose(got, BesselWronskianJY(x), 1e-9) {
		t.Errorf("JY Wronskian = %.12g, want %.12g", got, BesselWronskianJY(x))
	}
}

func BenchmarkAiryAi(b *testing.B) {
	var s float64
	for i := 0; i < b.N; i++ {
		s += AiryAi(float64(i%13)*0.5 - 3.0)
	}
	_ = s
}
