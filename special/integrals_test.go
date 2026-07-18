package special

import (
	"math"
	"testing"
)

// close reports whether got is within tol (absolute) of want, or within a
// relative tolerance for large magnitudes.
func specialClose(got, want, tol float64) bool {
	if math.IsNaN(got) || math.IsNaN(want) {
		return math.IsNaN(got) == math.IsNaN(want)
	}
	if math.IsInf(want, 0) {
		return math.IsInf(got, int(math.Copysign(1, want)))
	}
	d := math.Abs(got - want)
	if d <= tol {
		return true
	}
	return d <= tol*math.Max(1, math.Abs(want))
}

func TestErrorFunctions(t *testing.T) {
	cases := []struct {
		name string
		got  float64
		want float64
		tol  float64
	}{
		{"Erf(1)", Erf(1), 0.8427007929497149, 1e-12},
		{"Erfc(1)", Erfc(1), 0.15729920705028513, 1e-12},
		{"Erfcx(0)", Erfcx(0), 1, 1e-12},
		{"Erfcx(1)", Erfcx(1), 0.4275835761558070, 1e-9},
		{"Erfcx(5)", Erfcx(5), 0.1107046377166892, 1e-9},
		{"Erfi(1)", Erfi(1), 1.6504257587975428, 1e-9},
		{"Dawson(1)", Dawson(1), 0.5380795069127684, 1e-9},
		{"Dawson(2)", Dawson(2), 0.3013403889237919, 1e-9},
		{"Dawson(5)", Dawson(5), 0.1021340744242768, 1e-8},
		{"InverseErf(0.8427)", InverseErf(0.8427007929497149), 1, 1e-9},
		{"InverseErfc(0.1573)", InverseErfc(0.15729920705028513), 1, 1e-9},
	}
	for _, c := range cases {
		if !specialClose(c.got, c.want, c.tol) {
			t.Errorf("%s = %.15g, want %.15g", c.name, c.got, c.want)
		}
	}
}

func TestFresnel(t *testing.T) {
	cases := []struct {
		x     float64
		wantS float64
		wantC float64
		tol   float64
	}{
		{0, 0, 0, 1e-12},
		{0.5, 0.06473243285999929, 0.4923442258714464, 1e-9},
		{1, 0.4382591473903548, 0.7798934003768228, 1e-9},
		{2, 0.3434156783636982, 0.4882534060753408, 1e-8},
		{3, 0.4963129989673750, 0.6057207892976856, 1e-7},
	}
	for _, c := range cases {
		s, cc := Fresnel(c.x)
		if !specialClose(s, c.wantS, c.tol) {
			t.Errorf("FresnelS(%g) = %.12g, want %.12g", c.x, s, c.wantS)
		}
		if !specialClose(cc, c.wantC, c.tol) {
			t.Errorf("FresnelC(%g) = %.12g, want %.12g", c.x, cc, c.wantC)
		}
	}
}

func TestExponentialIntegrals(t *testing.T) {
	cases := []struct {
		name string
		got  float64
		want float64
		tol  float64
	}{
		{"Ei(1)", Ei(1), 1.8951178163559368, 1e-10},
		{"Ei(2)", Ei(2), 4.954234356001890, 1e-9},
		{"Ei(10)", Ei(10), 2492.2289762418777, 1e-6},
		{"E1(1)", E1(1), 0.21938393439552029, 1e-11},
		{"E1(2)", E1(2), 0.04890051070806112, 1e-11},
		{"E1(0.5)", E1(0.5), 0.5597735947761608, 1e-11},
		{"En(2,1)", En(2, 1), 0.14849550677592205, 1e-10},
		{"En(3,1)", En(3, 1), 0.10969197582165664, 1e-7},
		{"Li(2)", Li(2), 1.0451637801174927, 1e-10},
		{"Li(10)", Li(10), 6.165599504787280, 1e-8},
	}
	for _, c := range cases {
		if !specialClose(c.got, c.want, c.tol) {
			t.Errorf("%s = %.15g, want %.15g", c.name, c.got, c.want)
		}
	}
}

func TestSineCosineIntegrals(t *testing.T) {
	cases := []struct {
		name string
		got  float64
		want float64
		tol  float64
	}{
		{"Si(1)", Si(1), 0.9460830703671830, 1e-10},
		{"Si(5)", Si(5), 1.5499312449446741, 1e-8},
		{"Si(10)", Si(10), 1.6583475942188740, 1e-7},
		{"Ci(1)", Ci(1), 0.3374039229009681, 1e-10},
		{"Ci(5)", Ci(5), -0.19002974965664387, 1e-8},
		{"Ci(10)", Ci(10), -0.045456433004455, 1e-7},
		{"Shi(1)", Shi(1), 1.0572508753757285, 1e-10},
		{"Chi(1)", Chi(1), 0.8378669409802082, 1e-10},
		{"Shi(2)", Shi(2), 2.5015674333549756, 1e-9},
	}
	for _, c := range cases {
		if !specialClose(c.got, c.want, c.tol) {
			t.Errorf("%s = %.15g, want %.15g", c.name, c.got, c.want)
		}
	}
}

func TestElliptic(t *testing.T) {
	cases := []struct {
		name string
		got  float64
		want float64
		tol  float64
	}{
		{"EllipticK(0)", EllipticK(0), math.Pi / 2, 1e-12},
		{"EllipticK(0.5)", EllipticK(0.5), 1.8540746773013719, 1e-11},
		{"EllipticE(0)", EllipticE(0), math.Pi / 2, 1e-12},
		{"EllipticE(0.5)", EllipticE(0.5), 1.3506438810427978, 1e-11},
		{"EllipticE(1)", EllipticE(1), 1, 1e-12},
		{"EllipticF(pi/2,0.5)", EllipticF(math.Pi/2, 0.5), 1.8540746773013719, 1e-11},
		{"EllipticEInc(pi/2,0.5)", EllipticEInc(math.Pi/2, 0.5), 1.3506438810427978, 1e-11},
		{"EllipticF(pi/4,0.5)", EllipticF(math.Pi/4, 0.5), 0.826017876249246, 1e-10},
		{"EllipticPi(0.5,0)", EllipticPi(0.5, 0), 2.2214414690791831, 1e-10},
	}
	for _, c := range cases {
		if !specialClose(c.got, c.want, c.tol) {
			t.Errorf("%s = %.15g, want %.15g", c.name, c.got, c.want)
		}
	}
}

func TestZetaFamily(t *testing.T) {
	cases := []struct {
		name string
		got  float64
		want float64
		tol  float64
	}{
		{"Zeta(2)", Zeta(2), math.Pi * math.Pi / 6, 1e-10},
		{"Zeta(4)", Zeta(4), math.Pow(math.Pi, 4) / 90, 1e-10},
		{"Zeta(3)", Zeta(3), 1.2020569031595942, 1e-9},
		{"Zeta(0)", Zeta(0), -0.5, 1e-9},
		{"Zeta(-1)", Zeta(-1), -1.0 / 12, 1e-9},
		{"Zeta(6)", Zeta(6), math.Pow(math.Pi, 6) / 945, 1e-9},
		{"DirichletEta(1)", DirichletEta(1), math.Ln2, 1e-12},
		{"DirichletEta(2)", DirichletEta(2), math.Pi * math.Pi / 12, 1e-10},
		{"DirichletBeta(1)", DirichletBeta(1), math.Pi / 4, 1e-9},
		{"DirichletBeta(2)", DirichletBeta(2), 0.9159655941772190, 1e-9},
		{"HurwitzZeta(2,1)", HurwitzZeta(2, 1), math.Pi * math.Pi / 6, 1e-10},
	}
	for _, c := range cases {
		if !specialClose(c.got, c.want, c.tol) {
			t.Errorf("%s = %.15g, want %.15g", c.name, c.got, c.want)
		}
	}
}

func TestPolylog(t *testing.T) {
	cases := []struct {
		name string
		got  float64
		want float64
		tol  float64
	}{
		{"Li2(1)", Li2(1), math.Pi * math.Pi / 6, 1e-12},
		{"Li2(-1)", Li2(-1), -math.Pi * math.Pi / 12, 1e-10},
		{"Li2(0.5)", Li2(0.5), 0.5822405264650125, 1e-10},
		{"Li2(-2)", Li2(-2), -1.4367463668836082, 1e-9},
		{"Li3(1)", Li3(1), 1.2020569031595942, 1e-9},
		{"Li3(0.5)", Li3(0.5), 0.5372131936080402, 1e-9},
		{"Polylog(2,0.5)", Polylog(2, 0.5), 0.5822405264650125, 1e-10},
	}
	for _, c := range cases {
		if !specialClose(c.got, c.want, c.tol) {
			t.Errorf("%s = %.15g, want %.15g", c.name, c.got, c.want)
		}
	}
}

func TestLambertW(t *testing.T) {
	cases := []struct {
		name string
		got  float64
		want float64
		tol  float64
	}{
		{"LambertW(0)", LambertW(0), 0, 1e-12},
		{"LambertW(1)", LambertW(1), 0.5671432904097838, 1e-12},
		{"LambertW(e)", LambertW(math.E), 1, 1e-12},
		{"LambertW(10)", LambertW(10), 1.7455280027406994, 1e-11},
		{"LambertW(-0.2)", LambertW(-0.2), -0.2591711018190738, 1e-10},
		{"LambertWm1(-0.2)", LambertWm1(-0.2), -2.542641357773526, 1e-9},
		{"LambertWm1(-0.1)", LambertWm1(-0.1), -3.577152063957297, 1e-9},
	}
	for _, c := range cases {
		if !specialClose(c.got, c.want, c.tol) {
			t.Errorf("%s = %.15g, want %.15g", c.name, c.got, c.want)
		}
	}
	// Round-trip check: W(x)·e^{W(x)} = x.
	for _, x := range []float64{0.5, 1, 3, 20, 100} {
		w := LambertW(x)
		if !specialClose(w*math.Exp(w), x, 1e-9) {
			t.Errorf("LambertW round-trip failed at x=%g", x)
		}
	}
}

func TestGammaFamily(t *testing.T) {
	cases := []struct {
		name string
		got  float64
		want float64
		tol  float64
	}{
		{"Digamma(1)", Digamma(1), -specialEulerGamma, 1e-11},
		{"Digamma(2)", Digamma(2), 1 - specialEulerGamma, 1e-11},
		{"Digamma(0.5)", Digamma(0.5), -specialEulerGamma - 2*math.Ln2, 1e-11},
		{"Trigamma(1)", Trigamma(1), math.Pi * math.Pi / 6, 1e-10},
		{"Trigamma(2)", Trigamma(2), math.Pi*math.Pi/6 - 1, 1e-10},
		{"Polygamma(1,1)", Polygamma(1, 1), math.Pi * math.Pi / 6, 1e-9},
		{"Polygamma(2,1)", Polygamma(2, 1), -2 * 1.2020569031595942, 1e-8},
		{"Beta(2,3)", Beta(2, 3), 1.0 / 12, 1e-12},
		{"GammaP(2,1)", GammaP(2, 1), 0.2642411176571153, 1e-11},
		{"GammaQ(2,1)", GammaQ(2, 1), 0.7357588823428847, 1e-11},
		{"GammaP(3,5)", GammaP(3, 5), 0.8753479805169189, 1e-10},
		{"BetaInc(2,3,0.5)", BetaInc(2, 3, 0.5), 0.6875, 1e-11},
		{"BetaInc(0.5,0.5,0.5)", BetaInc(0.5, 0.5, 0.5), 0.5, 1e-11},
	}
	for _, c := range cases {
		if !specialClose(c.got, c.want, c.tol) {
			t.Errorf("%s = %.15g, want %.15g", c.name, c.got, c.want)
		}
	}
	// P + Q = 1.
	if !specialClose(GammaP(2.5, 3)+GammaQ(2.5, 3), 1, 1e-12) {
		t.Errorf("GammaP+GammaQ != 1")
	}
}

func TestBessel(t *testing.T) {
	cases := []struct {
		name string
		got  float64
		want float64
		tol  float64
	}{
		{"BesselI0(1)", BesselI0(1), 1.2660658777520084, 1e-7},
		{"BesselI1(1)", BesselI1(1), 0.5651591039924850, 1e-7},
		{"BesselI(2,1)", BesselI(2, 1), 0.1357476697670383, 1e-7},
		{"BesselK0(1)", BesselK0(1), 0.4210244382407083, 1e-7},
		{"BesselK1(1)", BesselK1(1), 0.6019072301972346, 1e-7},
		{"BesselK(2,1)", BesselK(2, 1), 1.6248388986351774, 1e-6},
		{"BesselJ0(1)", BesselJ0(1), 0.7651976865579666, 1e-9},
		{"BesselY0(1)", BesselY0(1), 0.08825696421567697, 1e-9},
	}
	for _, c := range cases {
		if !specialClose(c.got, c.want, c.tol) {
			t.Errorf("%s = %.15g, want %.15g", c.name, c.got, c.want)
		}
	}
}

func TestMisc(t *testing.T) {
	if !specialClose(Sinc(0), 1, 1e-12) {
		t.Errorf("Sinc(0) wrong")
	}
	if !specialClose(SincNorm(1), 0, 1e-12) {
		t.Errorf("SincNorm(1) wrong")
	}
	if !specialClose(Sinc(math.Pi), 0, 1e-12) {
		t.Errorf("Sinc(pi) wrong")
	}
	if !specialClose(ReciprocalGamma(5), 1.0/24, 1e-12) {
		t.Errorf("ReciprocalGamma(5) wrong")
	}
	if !specialClose(Gamma(5), 24, 1e-10) {
		t.Errorf("Gamma(5) wrong")
	}
	if !specialClose(LogGamma(5), math.Log(24), 1e-12) {
		t.Errorf("LogGamma(5) wrong")
	}
	// Struve H0 and H1 reference values (from the defining series).
	if !specialClose(Struve(0, 1), 0.5686566270482879, 1e-9) {
		t.Errorf("Struve(0,1) wrong: %v", Struve(0, 1))
	}
	if !specialClose(Struve(1, 2), 0.6467637282835621, 1e-9) {
		t.Errorf("Struve(1,2) wrong: %v", Struve(1, 2))
	}
}

// BenchmarkEllipticPiInc exercises the heaviest routine (Carlson R_J via R_C).
func BenchmarkEllipticPiInc(b *testing.B) {
	var s float64
	for i := 0; i < b.N; i++ {
		s += EllipticPiInc(0.3, 1.2, 0.5)
	}
	_ = s
}
