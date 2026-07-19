package analyticnt_test

import (
	"fmt"
	"math"
	"math/cmplx"
	"testing"

	ant "github.com/malcolmston/algebra/analyticnt"
)

// approx reports whether got and want agree within an absolute tolerance tol,
// or within a relative tolerance tol when want is large.
func approx(got, want, tol float64) bool {
	d := math.Abs(got - want)
	if d <= tol {
		return true
	}
	return d <= tol*math.Abs(want)
}

func TestPrimeCountingExact(t *testing.T) {
	tests := []struct {
		x    int64
		want int64
	}{
		{0, 0}, {1, 0}, {2, 1}, {10, 4}, {100, 25}, {1000, 168},
		{10000, 1229}, {100000, 9592},
	}
	for _, tc := range tests {
		if got := ant.PrimePi(tc.x); got != tc.want {
			t.Errorf("PrimePi(%d)=%d want %d", tc.x, got, tc.want)
		}
	}
	// Cross-check the three exact algorithms.
	for _, x := range []int64{50, 500, 5000, 50000} {
		v, ok := ant.PrimePiCheck(x)
		if !ok {
			t.Errorf("PrimePiCheck(%d): algorithms disagree", x)
		}
		if v != ant.PrimePi(x) {
			t.Errorf("PrimePiCheck(%d)=%d mismatch", x, v)
		}
		if ant.PrimePiLegendre(x) != ant.PrimePi(x) {
			t.Errorf("PrimePiLegendre(%d) mismatch", x)
		}
		if ant.PrimePiMeissel(x) != ant.PrimePi(x) {
			t.Errorf("PrimePiMeissel(%d) mismatch", x)
		}
	}
}

func TestPrimesAndSieve(t *testing.T) {
	got := ant.PrimesUpTo(20)
	want := []int64{2, 3, 5, 7, 11, 13, 17, 19}
	if len(got) != len(want) {
		t.Fatalf("PrimesUpTo(20) len=%d want %d", len(got), len(want))
	}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("PrimesUpTo(20)[%d]=%d want %d", i, got[i], want[i])
		}
	}
	if !ant.IsPrime(97) || ant.IsPrime(91) || ant.IsPrime(1) {
		t.Errorf("IsPrime failures")
	}
	if ant.NextPrime(13) != 17 || ant.PrevPrime(13) != 11 {
		t.Errorf("NextPrime/PrevPrime failure")
	}
	if ant.NthPrime(1) != 2 || ant.NthPrime(10) != 29 || ant.NthPrime(100) != 541 {
		t.Errorf("NthPrime failure")
	}
	ok, p, k := ant.IsPrimePower(27)
	if !ok || p != 3 || k != 3 {
		t.Errorf("IsPrimePower(27)=(%v,%d,%d)", ok, p, k)
	}
	if ok, _, _ := ant.IsPrimePower(12); ok {
		t.Errorf("IsPrimePower(12) should be false")
	}
}

func TestLiEi(t *testing.T) {
	tests := []struct {
		name string
		got  float64
		want float64
		tol  float64
	}{
		{"Ei(1)", ant.Ei(1), 1.8951178163559368, 1e-9},
		{"Ei(2)", ant.Ei(2), 4.954234356001890, 1e-8},
		{"E1(1)", ant.E1(1), 0.21938393439552029, 1e-9},
		{"E1(2)", ant.E1(2), 0.04890051070806112, 1e-9},
		{"li(2)", ant.LiValueAt2(), 1.0451637801174927, 1e-9},
		{"li(1000)", ant.Li(1000), 177.6096580155163, 1e-6},
		{"Li(1000)", ant.LiOffset(1000), 176.56449421, 1e-5},
	}
	for _, tc := range tests {
		if !approx(tc.got, tc.want, tc.tol) {
			t.Errorf("%s=%.12g want %.12g", tc.name, tc.got, tc.want)
		}
	}
	// li vanishes at the Ramanujan–Soldner constant.
	if math.Abs(ant.Li(ant.SoldnerConstant)) > 1e-9 {
		t.Errorf("Li(SoldnerConstant)=%g want 0", ant.Li(ant.SoldnerConstant))
	}
}

func TestRiemannR(t *testing.T) {
	// Both formulations must agree, and R must beat li at estimating pi.
	for _, x := range []float64{100, 1000, 10000} {
		r1 := ant.RiemannR(x)
		r2 := ant.RiemannRScaled(x)
		if !approx(r1, r2, 3e-3) {
			t.Errorf("RiemannR(%g)=%g vs Scaled %g", x, r1, r2)
		}
	}
	if !approx(ant.RiemannR(1000), 168.359446, 1e-4) {
		t.Errorf("RiemannR(1000)=%g", ant.RiemannR(1000))
	}
	// R(1e4) closer to pi(1e4)=1229 than Li.
	pi := 1229.0
	if math.Abs(ant.RiemannR(1e4)-pi) > math.Abs(ant.LiOffset(1e4)-pi) {
		t.Errorf("RiemannR should beat Li at x=1e4")
	}
}

func TestZetaReal(t *testing.T) {
	tests := []struct {
		s    float64
		want float64
		tol  float64
	}{
		{2, math.Pi * math.Pi / 6, 1e-10},
		{4, math.Pow(math.Pi, 4) / 90, 1e-10},
		{3, 1.2020569031595943, 1e-9},
		{6, 1.0173430619844492, 1e-9},
		{-1, -1.0 / 12, 1e-9},
		{0, -0.5, 1e-9},
		{-3, 1.0 / 120, 1e-8},
		{0.5, -1.4603545088095868, 1e-6},
		{10, 1.0009945751278182, 1e-10},
	}
	for _, tc := range tests {
		if got := ant.Zeta(tc.s); !approx(got, tc.want, tc.tol) {
			t.Errorf("Zeta(%g)=%.12g want %.12g", tc.s, got, tc.want)
		}
	}
	// Trivial zeros.
	if ant.Zeta(-2) != 0 || ant.Zeta(-4) != 0 {
		t.Errorf("even negative integers must give exact zero")
	}
}

func TestZetaComplexAndEta(t *testing.T) {
	// zeta(2+0i) should equal real zeta(2).
	z := ant.ZetaComplex(complex(2, 0))
	if !approx(real(z), math.Pi*math.Pi/6, 1e-9) || math.Abs(imag(z)) > 1e-9 {
		t.Errorf("ZetaComplex(2)=%v", z)
	}
	if got := ant.DirichletEta(1); !approx(got, math.Ln2, 1e-12) {
		t.Errorf("DirichletEta(1)=%g want ln2", got)
	}
	if got := ant.DirichletEta(2); !approx(got, math.Pi*math.Pi/12, 1e-9) {
		t.Errorf("DirichletEta(2)=%g want pi^2/12", got)
	}
	// eta on the critical line is finite.
	e := ant.DirichletEtaComplex(complex(0.5, 10))
	if cmplx.IsNaN(e) {
		t.Errorf("DirichletEtaComplex NaN")
	}
}

func TestZetaZeros(t *testing.T) {
	want := []float64{14.134725141734693, 21.02203963877155, 25.010857580145688,
		30.424876125859513, 32.935061587739189}
	got := ant.ZetaZeros(5)
	if len(got) != 5 {
		t.Fatalf("ZetaZeros(5) len=%d", len(got))
	}
	for i := range want {
		if !approx(got[i], want[i], 1e-4) {
			t.Errorf("zero %d = %.9f want %.9f", i+1, got[i], want[i])
		}
	}
	if !approx(ant.ZetaZero(1), 14.134725, 1e-4) {
		t.Errorf("ZetaZero(1)=%g", ant.ZetaZero(1))
	}
	// HardyZ changes sign across the first zero.
	if ant.HardyZ(14.0) > 0 == (ant.HardyZ(14.2) > 0) {
		t.Errorf("HardyZ should change sign across first zero")
	}
	// Riemann–von Mangoldt smooth count just past the 5th zero (~32.94) should
	// be close to 5.
	if n := ant.RiemannVonMangoldtN(35); n < 4.3 || n > 5.6 {
		t.Errorf("RiemannVonMangoldtN(35)=%g", n)
	}
}

func TestChebyshev(t *testing.T) {
	if !approx(ant.ChebyshevTheta(10), math.Log(210), 1e-9) {
		t.Errorf("theta(10)=%g want ln210", ant.ChebyshevTheta(10))
	}
	psiWant := 3*math.Ln2 + 2*math.Log(3) + math.Log(5) + math.Log(7)
	if !approx(ant.ChebyshevPsi(10), psiWant, 1e-9) {
		t.Errorf("psi(10)=%g want %g", ant.ChebyshevPsi(10), psiWant)
	}
	// psi(x) ~ x.
	if e := ant.PsiError(100000); math.Abs(e) > 5 {
		t.Errorf("PsiError(1e5)=%g unexpectedly large", e)
	}
	if !approx(ant.PrimeHarmonic(10), 1.0/2+1.0/3+1.0/5+1.0/7, 1e-12) {
		t.Errorf("PrimeHarmonic(10) wrong")
	}
}

func TestVonMangoldt(t *testing.T) {
	tests := []struct {
		n    int64
		want float64
	}{
		{1, 0}, {2, math.Ln2}, {4, math.Ln2}, {8, math.Ln2},
		{9, math.Log(3)}, {7, math.Log(7)}, {6, 0}, {12, 0},
	}
	for _, tc := range tests {
		if got := ant.VonMangoldt(tc.n); !approx(got, tc.want, 1e-12) {
			t.Errorf("VonMangoldt(%d)=%g want %g", tc.n, got, tc.want)
		}
	}
}

func TestMobiusMertens(t *testing.T) {
	muTests := []struct {
		n    int64
		want int
	}{
		{1, 1}, {2, -1}, {3, -1}, {4, 0}, {6, 1}, {30, -1}, {12, 0}, {210, 1},
	}
	for _, tc := range muTests {
		if got := ant.MobiusMu(tc.n); got != tc.want {
			t.Errorf("MobiusMu(%d)=%d want %d", tc.n, got, tc.want)
		}
	}
	// Sieve must match direct.
	mu := ant.MobiusSieve(100)
	for n := int64(1); n <= 100; n++ {
		if mu[n] != ant.MobiusMu(n) {
			t.Errorf("MobiusSieve[%d]=%d disagrees", n, mu[n])
		}
	}
	mTests := []struct {
		n    int
		want int64
	}{
		{1, 1}, {10, -1}, {100, 1}, {1000, 2},
	}
	for _, tc := range mTests {
		if got := ant.MertensFunction(tc.n); got != tc.want {
			t.Errorf("MertensFunction(%d)=%d want %d", tc.n, got, tc.want)
		}
	}
}

func TestMultiplicative(t *testing.T) {
	if ant.EulerPhi(36) != 12 || ant.EulerPhi(1) != 1 || ant.EulerPhi(10) != 4 {
		t.Errorf("EulerPhi failure")
	}
	if ant.DivisorSigma(1, 6) != 12 || ant.DivisorSigma(0, 12) != 6 || ant.DivisorSigma(2, 3) != 10 {
		t.Errorf("DivisorSigma failure")
	}
	if ant.DivisorCount(12) != 6 || ant.Omega(12) != 2 || ant.BigOmega(12) != 3 {
		t.Errorf("divisor/omega failure")
	}
	if ant.Radical(360) != 30 || ant.Liouville(12) != -1 || ant.Liouville(4) != 1 {
		t.Errorf("radical/liouville failure")
	}
	if ant.JordanTotient(2, 6) != 24 || ant.TotientSummatory(10) != 32 {
		t.Errorf("jordan/totient-summatory failure")
	}
	if !ant.IsSquareFree(30) || ant.IsSquareFree(12) {
		t.Errorf("IsSquareFree failure")
	}
}

func TestConstantsAndZetaValues(t *testing.T) {
	if !approx(ant.EulerGammaExp(), math.Exp(-ant.EulerGamma), 1e-15) {
		t.Errorf("EulerGammaExp failure")
	}
	if !approx(ant.ZetaAtEvenInteger(1), math.Pi*math.Pi/6, 1e-12) {
		t.Errorf("ZetaAtEvenInteger(1) failure")
	}
	if !approx(ant.PrimeZeta(2), 0.4522474200410655, 1e-9) {
		t.Errorf("PrimeZeta(2)=%g", ant.PrimeZeta(2))
	}
	if !approx(ant.AperyConstant, 1.2020569031595942, 1e-12) {
		t.Errorf("AperyConstant wrong")
	}
	// Mertens' third theorem: product ~ e^{-gamma}/ln x.
	x := 100000.0
	got := ant.MertensProduct(x)
	want := math.Exp(-ant.EulerGamma) / math.Log(x)
	if !approx(got, want, 0.02) {
		t.Errorf("MertensProduct(%g)=%g want ~%g", x, got, want)
	}
}

func TestDirichletCharacters(t *testing.T) {
	chi := ant.LegendreCharacter(7)
	if !chi.IsReal() || chi.IsPrincipal() || chi.Order() != 2 {
		t.Errorf("Legendre character mod 7 properties wrong")
	}
	// (a|7): QRs mod 7 are {1,2,4}.
	wants := map[int64]int{1: 1, 2: 1, 3: -1, 4: 1, 5: -1, 6: -1}
	for a, w := range wants {
		if real(chi.Eval(a)) != float64(w) {
			t.Errorf("chi7(%d)=%v want %d", a, chi.Eval(a), w)
		}
	}
	if ant.LegendreSymbol(2, 7) != 1 || ant.LegendreSymbol(3, 7) != -1 {
		t.Errorf("LegendreSymbol failure")
	}
	if ant.JacobiSymbol(2, 15) != 1 || ant.JacobiSymbol(7, 15) != -1 {
		t.Errorf("JacobiSymbol failure")
	}
	// Gauss sum magnitude for a primitive character mod q is sqrt(q).
	if g := ant.GaussSum(chi); !approx(cmplx.Abs(g), math.Sqrt(7), 1e-9) {
		t.Errorf("GaussSum magnitude=%g want sqrt7", cmplx.Abs(g))
	}
	// Class numbers of imaginary quadratic fields.
	if ant.ClassNumberDirichlet(23) != 3 || ant.ClassNumberDirichlet(163) != 1 {
		t.Errorf("ClassNumberDirichlet failure")
	}
	// Principal character mod 1 gives zeta.
	pc := ant.PrincipalCharacter(1)
	l := ant.DirichletL(complex(2, 0), pc)
	if !approx(real(l), math.Pi*math.Pi/6, 1e-8) {
		t.Errorf("L(2, triv)=%v want zeta(2)", l)
	}
}

func TestLFunctions(t *testing.T) {
	// Non-principal character mod 4: L(2,chi)=Catalan, L(1,chi)=pi/4.
	c4 := ant.CharacterFromPrimitiveRoot(4, 1)
	l2 := ant.DirichletL(complex(2, 0), c4)
	if !approx(real(l2), ant.CatalanConstant, 1e-8) {
		t.Errorf("L(2, chi4)=%v want Catalan", l2)
	}
	if !approx(ant.DirichletBeta(2), ant.CatalanConstant, 1e-9) {
		t.Errorf("DirichletBeta(2)=%g want Catalan", ant.DirichletBeta(2))
	}
	if !approx(ant.DirichletBeta(1), math.Pi/4, 1e-12) {
		t.Errorf("DirichletBeta(1)=%g want pi/4", ant.DirichletBeta(1))
	}
	// Hurwitz zeta reduces to Riemann zeta at a=1.
	if !approx(ant.HurwitzZeta(2, 1), math.Pi*math.Pi/6, 1e-9) {
		t.Errorf("HurwitzZeta(2,1) failure")
	}
	if !approx(ant.DirichletLambda(2), math.Pi*math.Pi/8, 1e-9) {
		t.Errorf("DirichletLambda(2)=%g want pi^2/8", ant.DirichletLambda(2))
	}
}

func TestGapsAndTwins(t *testing.T) {
	if ant.TwinPrimeCount(100) != 8 {
		t.Errorf("TwinPrimeCount(100)=%d want 8", ant.TwinPrimeCount(100))
	}
	twins := ant.TwinPrimesUpTo(20)
	want := [][2]int64{{3, 5}, {5, 7}, {11, 13}, {17, 19}}
	if len(twins) != len(want) {
		t.Fatalf("TwinPrimesUpTo(20) len=%d", len(twins))
	}
	for i := range want {
		if twins[i] != want[i] {
			t.Errorf("twin %d=%v want %v", i, twins[i], want[i])
		}
	}
	if g, _ := ant.MaxPrimeGap(1000); g != 20 {
		t.Errorf("MaxPrimeGap(1000)=%d want 20", g)
	}
	if ant.GoldbachPartitions(100) != 6 {
		t.Errorf("GoldbachPartitions(100)=%d want 6", ant.GoldbachPartitions(100))
	}
	p, q, ok := ant.GoldbachPair(100)
	if !ok || p+q != 100 || !ant.IsPrime(p) || !ant.IsPrime(q) {
		t.Errorf("GoldbachPair(100)=(%d,%d,%v)", p, q, ok)
	}
	// Deterministic random prime with a fixed seed is reproducible.
	if ant.RandomPrimeBelow(100, 42) != ant.RandomPrimeBelow(100, 42) {
		t.Errorf("RandomPrimeBelow not deterministic")
	}
	// Average gap grows like ln x.
	if a := ant.AverageGap(100000); a < 8 || a > 14 {
		t.Errorf("AverageGap(1e5)=%g out of expected range", a)
	}
}

func TestExtras(t *testing.T) {
	if !approx(ant.Bernoulli(2), 1.0/6, 1e-12) || !approx(ant.Bernoulli(4), -1.0/30, 1e-12) {
		t.Errorf("Bernoulli failure")
	}
	if !approx(ant.Digamma(1), -ant.EulerGamma, 1e-8) {
		t.Errorf("Digamma(1)=%g want -gamma", ant.Digamma(1))
	}
	if !approx(ant.Trigamma(1), math.Pi*math.Pi/6, 1e-6) {
		t.Errorf("Trigamma(1)=%g want pi^2/6", ant.Trigamma(1))
	}
	if ant.SquareFreeCount(100) != 61 || ant.SquareFreeCount(1000) != 608 {
		t.Errorf("SquareFreeCount failure")
	}
	if ant.SemiprimeCount(100) != 34 {
		t.Errorf("SemiprimeCount(100)=%d want 34", ant.SemiprimeCount(100))
	}
	if !ant.IsWilsonPrime(5) || !ant.IsWilsonPrime(13) || ant.IsWilsonPrime(7) {
		t.Errorf("IsWilsonPrime failure")
	}
	if ant.GCD(48, 36) != 12 || ant.LCM(4, 6) != 12 || !ant.IsCoprime(9, 10) {
		t.Errorf("gcd/lcm/coprime failure")
	}
	if !ant.IsSophieGermain(23) || !ant.IsSafePrime(47) {
		t.Errorf("SophieGermain/SafePrime failure")
	}
	fac := ant.FactorInt(360)
	if fac[2] != 3 || fac[3] != 2 || fac[5] != 1 {
		t.Errorf("FactorInt(360)=%v", fac)
	}
	// Riemann xi at 1/2 is a known real value.
	if xi := ant.RiemannXiReal(0.5); !approx(real(xi), 0.4971207782, 1e-6) {
		t.Errorf("RiemannXiReal(0.5)=%v", xi)
	}
	// Chebyshev bias: N(3 mod 4) leads N(1 mod 4) up to 1000.
	if ant.ChebyshevBiasCount(1000) <= 0 {
		t.Errorf("ChebyshevBiasCount(1000) should be positive")
	}
}

func TestNthPrimeApproximations(t *testing.T) {
	// Inverting R should land within a few of the true nth prime.
	for _, n := range []int{100, 1000} {
		est := ant.NthPrimeInverseR(n)
		true := float64(ant.NthPrime(n))
		if math.Abs(est-true) > 0.02*true {
			t.Errorf("NthPrimeInverseR(%d)=%g true %g", n, est, true)
		}
	}
}

func TestPanics(t *testing.T) {
	mustPanic := func(name string, f func()) {
		defer func() {
			if recover() == nil {
				t.Errorf("%s did not panic", name)
			}
		}()
		f()
	}
	mustPanic("NthPrime(0)", func() { ant.NthPrime(0) })
	mustPanic("MobiusMu(0)", func() { ant.MobiusMu(0) })
	mustPanic("Li(-1)", func() { ant.Li(-1) })
	mustPanic("EulerPhi(0)", func() { ant.EulerPhi(0) })
}

// Examples -----------------------------------------------------------------

func ExamplePrimePi() {
	fmt.Println(ant.PrimePi(100))
	// Output: 25
}

func ExampleZeta() {
	fmt.Printf("%.6f\n", ant.Zeta(2))
	// Output: 1.644934
}

func ExampleZetaZero() {
	fmt.Printf("%.4f\n", ant.ZetaZero(1))
	// Output: 14.1347
}

func ExampleMobiusMu() {
	fmt.Println(ant.MobiusMu(30), ant.MobiusMu(12))
	// Output: -1 0
}

func ExampleLiOffset() {
	fmt.Printf("%.2f\n", ant.LiOffset(1000))
	// Output: 176.56
}

func ExampleTwinPrimesUpTo() {
	fmt.Println(ant.TwinPrimesUpTo(20))
	// Output: [[3 5] [5 7] [11 13] [17 19]]
}
