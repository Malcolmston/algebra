package bigfloat_test

import (
	"errors"
	"fmt"
	"math"
	"math/big"
	"testing"

	bf "github.com/malcolmston/algebra/bigfloat"
)

const testPrec = 160

// f64 converts a *big.Float to float64.
func f64(x *big.Float) float64 { v, _ := x.Float64(); return v }

// bfv builds a *big.Float from a float64 at test precision.
func bfv(x float64) *big.Float { return bf.FromFloat64(x, testPrec) }

// close reports whether got is within a relative/absolute tolerance of want.
func closeTo(got, want, tol float64) bool {
	return math.Abs(got-want) <= tol*(1+math.Abs(want))
}

func check(t *testing.T, name string, got *big.Float, want, tol float64) {
	t.Helper()
	g := f64(got)
	if !closeTo(g, want, tol) {
		t.Errorf("%s = %.17g, want %.17g (diff %.3g)", name, g, want, g-want)
	}
}

func TestConstants(t *testing.T) {
	const tol = 1e-14
	tests := []struct {
		name string
		got  *big.Float
		want float64
	}{
		{"Pi", bf.Pi(testPrec), math.Pi},
		{"TwoPi", bf.TwoPi(testPrec), 2 * math.Pi},
		{"HalfPi", bf.HalfPi(testPrec), math.Pi / 2},
		{"QuarterPi", bf.QuarterPi(testPrec), math.Pi / 4},
		{"InvPi", bf.InvPi(testPrec), 1 / math.Pi},
		{"TwoInvPi", bf.TwoInvPi(testPrec), 2 / math.Pi},
		{"E", bf.E(testPrec), math.E},
		{"InvE", bf.InvE(testPrec), 1 / math.E},
		{"Ln2", bf.Ln2(testPrec), math.Ln2},
		{"Ln10", bf.Ln10(testPrec), math.Log(10)},
		{"LnPi", bf.LnPi(testPrec), math.Log(math.Pi)},
		{"Ln2Pi", bf.Ln2Pi(testPrec), math.Log(2 * math.Pi)},
		{"Log2E", bf.Log2E(testPrec), math.Log2E},
		{"Log10E", bf.Log10E(testPrec), math.Log10E},
		{"Sqrt2", bf.Sqrt2(testPrec), math.Sqrt2},
		{"Sqrt3", bf.Sqrt3(testPrec), math.Sqrt(3)},
		{"Sqrt5", bf.Sqrt5(testPrec), math.Sqrt(5)},
		{"InvSqrt2", bf.InvSqrt2(testPrec), 1 / math.Sqrt2},
		{"SqrtPi", bf.SqrtPi(testPrec), math.Sqrt(math.Pi)},
		{"SqrtTwoPi", bf.SqrtTwoPi(testPrec), math.Sqrt(2 * math.Pi)},
		{"InvSqrtPi", bf.InvSqrtPi(testPrec), 1 / math.Sqrt(math.Pi)},
		{"EulerGamma", bf.EulerGamma(testPrec), 0.5772156649015329},
		{"Catalan", bf.Catalan(testPrec), 0.915965594177219},
		{"Apery", bf.Apery(testPrec), 1.2020569031595943},
		{"Zeta3", bf.Zeta3(testPrec), 1.2020569031595943},
		{"GoldenRatio", bf.GoldenRatio(testPrec), (1 + math.Sqrt(5)) / 2},
		{"InvGoldenRatio", bf.InvGoldenRatio(testPrec), (math.Sqrt(5) - 1) / 2},
		{"SilverRatio", bf.SilverRatio(testPrec), 1 + math.Sqrt2},
		{"PlasticNumber", bf.PlasticNumber(testPrec), 1.3247179572447458},
		{"Degree", bf.Degree(testPrec), math.Pi / 180},
	}
	for _, tc := range tests {
		check(t, tc.name, tc.got, tc.want, tol)
	}
}

func TestExpLog(t *testing.T) {
	const tol = 1e-14
	xs := []float64{0, 1, -1, 0.5, -0.5, 2, -3.5, 10, -10, 50, -50, 0.001, 123.456}
	for _, x := range xs {
		check(t, fmt.Sprintf("Exp(%g)", x), bf.Exp(bfv(x), testPrec), math.Exp(x), tol)
		check(t, fmt.Sprintf("Expm1(%g)", x), bf.Expm1(bfv(x), testPrec), math.Expm1(x), tol)
		check(t, fmt.Sprintf("Exp2(%g)", x), bf.Exp2(bfv(x), testPrec), math.Exp2(x), tol)
		check(t, fmt.Sprintf("Sinh(%g)", x), bf.Sinh(bfv(x), testPrec), math.Sinh(x), tol)
		check(t, fmt.Sprintf("Cosh(%g)", x), bf.Cosh(bfv(x), testPrec), math.Cosh(x), tol)
		check(t, fmt.Sprintf("Tanh(%g)", x), bf.Tanh(bfv(x), testPrec), math.Tanh(x), tol)
	}
	pos := []float64{1, 2, 0.5, 7, 1e-3, 1e6, 123.456, math.E}
	for _, x := range pos {
		l, err := bf.Log(bfv(x), testPrec)
		if err != nil {
			t.Fatalf("Log(%g): %v", x, err)
		}
		check(t, fmt.Sprintf("Log(%g)", x), l, math.Log(x), tol)
		l2, _ := bf.Log2(bfv(x), testPrec)
		check(t, fmt.Sprintf("Log2(%g)", x), l2, math.Log2(x), tol)
		l10, _ := bf.Log10(bfv(x), testPrec)
		check(t, fmt.Sprintf("Log10(%g)", x), l10, math.Log10(x), tol)
	}
	// Log1p for small arguments.
	for _, x := range []float64{-0.5, 0, 1e-6, 5, 1e3} {
		l, err := bf.Log1p(bfv(x), testPrec)
		if err != nil {
			t.Fatalf("Log1p(%g): %v", x, err)
		}
		check(t, fmt.Sprintf("Log1p(%g)", x), l, math.Log1p(x), tol)
	}
	// Exp10.
	check(t, "Exp10(3)", bf.Exp10(bfv(3), testPrec), 1000, tol)
	check(t, "Exp10(-2)", bf.Exp10(bfv(-2), testPrec), 0.01, tol)
}

func TestPowRoots(t *testing.T) {
	const tol = 1e-14
	cases := []struct{ x, y, want float64 }{
		{2, 10, 1024}, {2, 0.5, math.Sqrt2}, {9, 0.5, 3},
		{2, -3, 0.125}, {5, 0, 1}, {1.5, 3.7, math.Pow(1.5, 3.7)},
		{10, 2.5, math.Pow(10, 2.5)},
	}
	for _, c := range cases {
		got, err := bf.Pow(bfv(c.x), bfv(c.y), testPrec)
		if err != nil {
			t.Fatalf("Pow(%g,%g): %v", c.x, c.y, err)
		}
		check(t, fmt.Sprintf("Pow(%g,%g)", c.x, c.y), got, c.want, tol)
	}
	// Negative base, integer exponent.
	g, err := bf.Pow(bfv(-2), bfv(3), testPrec)
	if err != nil {
		t.Fatal(err)
	}
	check(t, "Pow(-2,3)", g, -8, tol)

	check(t, "PowInt(1.5,20)", bf.PowInt(bfv(1.5), 20, testPrec), math.Pow(1.5, 20), tol)
	check(t, "PowInt(2,-5)", bf.PowInt(bfv(2), -5, testPrec), math.Pow(2, -5), tol)

	sq, err := bf.Sqrt(bfv(2), testPrec)
	if err != nil {
		t.Fatal(err)
	}
	check(t, "Sqrt(2)", sq, math.Sqrt2, tol)
	check(t, "Cbrt(27)", bf.Cbrt(bfv(27), testPrec), 3, tol)
	check(t, "Cbrt(-8)", bf.Cbrt(bfv(-8), testPrec), -2, tol)
	r5, err := bf.Root(bfv(32), 5, testPrec)
	if err != nil {
		t.Fatal(err)
	}
	check(t, "Root(32,5)", r5, 2, tol)
	check(t, "Hypot(3,4)", bf.Hypot(bfv(3), bfv(4), testPrec), 5, tol)
	rs, _ := bf.Rsqrt(bfv(4), testPrec)
	check(t, "Rsqrt(4)", rs, 0.5, tol)
}

func TestTrig(t *testing.T) {
	const tol = 1e-13
	xs := []float64{0, 0.3, 1, 1.2, -0.7, math.Pi / 6, math.Pi / 3, 2, 3, 100, -100, 1000}
	for _, x := range xs {
		check(t, fmt.Sprintf("Sin(%g)", x), bf.Sin(bfv(x), testPrec), math.Sin(x), tol)
		check(t, fmt.Sprintf("Cos(%g)", x), bf.Cos(bfv(x), testPrec), math.Cos(x), tol)
	}
	check(t, "Tan(1.2)", bf.Tan(bfv(1.2), testPrec), math.Tan(1.2), tol)
	check(t, "Cot(1.2)", bf.Cot(bfv(1.2), testPrec), 1/math.Tan(1.2), tol)
	check(t, "Sec(1)", bf.Sec(bfv(1), testPrec), 1/math.Cos(1), tol)
	check(t, "Csc(1)", bf.Csc(bfv(1), testPrec), 1/math.Sin(1), tol)
	s, c := bf.SinCos(bfv(0.7), testPrec)
	check(t, "SinCos.sin", s, math.Sin(0.7), tol)
	check(t, "SinCos.cos", c, math.Cos(0.7), tol)
	check(t, "Sinpi(0.25)", bf.Sinpi(bfv(0.25), testPrec), math.Sin(math.Pi*0.25), tol)
	check(t, "Cospi(0.25)", bf.Cospi(bfv(0.25), testPrec), math.Cos(math.Pi*0.25), tol)
	check(t, "Sinpi(3.5)", bf.Sinpi(bfv(3.5), testPrec), -1, tol)
	check(t, "Sinc(0)", bf.Sinc(bfv(0), testPrec), 1, tol)
	check(t, "Sinc(2)", bf.Sinc(bfv(2), testPrec), math.Sin(2)/2, tol)
	check(t, "Versin(1)", bf.Versin(bfv(1), testPrec), 1-math.Cos(1), tol)
	check(t, "Haversin(1)", bf.Haversin(bfv(1), testPrec), (1-math.Cos(1))/2, tol)
	check(t, "Deg2Rad(180)", bf.Deg2Rad(bfv(180), testPrec), math.Pi, tol)
	check(t, "Rad2Deg(pi)", bf.Rad2Deg(bf.Pi(testPrec), testPrec), 180, tol)
}

func TestInverseTrig(t *testing.T) {
	const tol = 1e-13
	for _, x := range []float64{-3, -1, -0.3, 0, 0.3, 1, 3, 100} {
		check(t, fmt.Sprintf("Atan(%g)", x), bf.Atan(bfv(x), testPrec), math.Atan(x), tol)
	}
	for _, x := range []float64{-0.9, -0.5, 0, 0.25, 0.6, 0.99} {
		as, err := bf.Asin(bfv(x), testPrec)
		if err != nil {
			t.Fatal(err)
		}
		check(t, fmt.Sprintf("Asin(%g)", x), as, math.Asin(x), tol)
		ac, _ := bf.Acos(bfv(x), testPrec)
		check(t, fmt.Sprintf("Acos(%g)", x), ac, math.Acos(x), tol)
	}
	// Boundary.
	as1, _ := bf.Asin(bfv(1), testPrec)
	check(t, "Asin(1)", as1, math.Pi/2, tol)
	cases := []struct {
		name       string
		y, x, want float64
	}{
		{"Atan2(1,1)", 1, 1, math.Atan2(1, 1)},
		{"Atan2(1,-1)", 1, -1, math.Atan2(1, -1)},
		{"Atan2(-1,-1)", -1, -1, math.Atan2(-1, -1)},
		{"Atan2(-1,1)", -1, 1, math.Atan2(-1, 1)},
		{"Atan2(0,-1)", 0, -1, math.Atan2(0, -1)},
		{"Atan2(3,4)", 3, 4, math.Atan2(3, 4)},
	}
	for _, c := range cases {
		check(t, c.name, bf.Atan2(bfv(c.y), bfv(c.x), testPrec), c.want, tol)
	}
	check(t, "Acot(2)", bf.Acot(bfv(2), testPrec), math.Atan(0.5), tol)
}

func TestInverseHyperbolic(t *testing.T) {
	const tol = 1e-13
	for _, x := range []float64{-5, -1, 0, 0.5, 1, 3, 20} {
		check(t, fmt.Sprintf("Asinh(%g)", x), bf.Asinh(bfv(x), testPrec), math.Asinh(x), tol)
	}
	for _, x := range []float64{1, 1.5, 2, 10, 100} {
		ac, err := bf.Acosh(bfv(x), testPrec)
		if err != nil {
			t.Fatal(err)
		}
		check(t, fmt.Sprintf("Acosh(%g)", x), ac, math.Acosh(x), tol)
	}
	for _, x := range []float64{-0.9, -0.5, 0, 0.3, 0.5, 0.99} {
		at, err := bf.Atanh(bfv(x), testPrec)
		if err != nil {
			t.Fatal(err)
		}
		check(t, fmt.Sprintf("Atanh(%g)", x), at, math.Atanh(x), tol)
	}
	acoth, _ := bf.Acoth(bfv(2), testPrec)
	check(t, "Acoth(2)", acoth, math.Atanh(0.5), tol)
	check(t, "Coth(1)", must(bf.Coth(bfv(1), testPrec)), 1/math.Tanh(1), tol)
	check(t, "Sech(1)", bf.Sech(bfv(1), testPrec), 1/math.Cosh(1), tol)
	check(t, "Csch(1)", must(bf.Csch(bfv(1), testPrec)), 1/math.Sinh(1), tol)
}

func must(x *big.Float, err error) *big.Float {
	if err != nil {
		panic(err)
	}
	return x
}

func TestGammaFamily(t *testing.T) {
	const tol = 1e-13
	gammaCases := []struct{ x, want float64 }{
		{1, 1}, {2, 1}, {3, 2}, {5, 24}, {6, 120},
		{0.5, math.Sqrt(math.Pi)}, {1.5, 0.5 * math.Sqrt(math.Pi)},
		{2.5, 1.3293403881791370}, {4.5, 11.631728396567448},
		{-0.5, -2 * math.Sqrt(math.Pi)}, {-1.5, 4.0 / 3.0 * math.Sqrt(math.Pi)},
		{0.1, 9.513507698668732}, {10, 362880},
	}
	for _, c := range gammaCases {
		g, err := bf.Gamma(bfv(c.x), 120)
		if err != nil {
			t.Fatalf("Gamma(%g): %v", c.x, err)
		}
		check(t, fmt.Sprintf("Gamma(%g)", c.x), g, c.want, tol)
	}
	// Lgamma value and sign.
	lg, sign, err := bf.Lgamma(bfv(-0.5), 120)
	if err != nil {
		t.Fatal(err)
	}
	check(t, "Lgamma(-0.5).val", lg, math.Log(2*math.Sqrt(math.Pi)), tol)
	if sign != -1 {
		t.Errorf("Lgamma(-0.5) sign = %d, want -1", sign)
	}
	// Digamma.
	digCases := []struct{ x, want float64 }{
		{1, -0.5772156649015329}, {2, 1 - 0.5772156649015329},
		{0.5, -1.9635100260214235}, {10, 2.2517525890667211},
		{-0.5, 0.036489973978576520},
	}
	for _, c := range digCases {
		d, err := bf.Digamma(bfv(c.x), 120)
		if err != nil {
			t.Fatalf("Digamma(%g): %v", c.x, err)
		}
		check(t, fmt.Sprintf("Digamma(%g)", c.x), d, c.want, tol)
	}
	// Beta and LogBeta.
	b, _ := bf.Beta(bfv(2), bfv(3), 120)
	check(t, "Beta(2,3)", b, 1.0/12.0, tol)
	b2, _ := bf.Beta(bfv(0.5), bfv(0.5), 120)
	check(t, "Beta(0.5,0.5)", b2, math.Pi, tol)
	lb, _ := bf.LogBeta(bfv(3), bfv(4), 120)
	check(t, "LogBeta(3,4)", lb, math.Log(1.0/60.0), tol)
	// Factorials and binomials.
	check(t, "Factorial(10)", bf.Factorial(10, testPrec), 3628800, tol)
	if got := bf.FactorialBig(20).String(); got != "2432902008176640000" {
		t.Errorf("FactorialBig(20) = %s", got)
	}
	check(t, "DoubleFactorial(9)", bf.DoubleFactorial(9, testPrec), 945, tol)
	check(t, "Binomial(10,3)", bf.Binomial(10, 3, testPrec), 120, tol)
	if got := bf.BinomialBig(52, 5).String(); got != "2598960" {
		t.Errorf("BinomialBig(52,5) = %s", got)
	}
	br, _ := bf.BinomialReal(bfv(10), bfv(3), 120)
	check(t, "BinomialReal(10,3)", br, 120, tol)
	check(t, "RisingFactorial(2,4)", bf.RisingFactorial(bfv(2), 4, testPrec), 120, tol)
	check(t, "FallingFactorial(6,3)", bf.FallingFactorial(bfv(6), 3, testPrec), 120, tol)
	// Bernoulli numbers (exact rationals).
	bernCases := []struct {
		n    uint
		want string
	}{
		{0, "1"}, {1, "-1/2"}, {2, "1/6"}, {4, "-1/30"},
		{6, "1/42"}, {8, "-1/30"}, {10, "5/66"}, {12, "-691/2730"},
	}
	for _, c := range bernCases {
		if got := bf.Bernoulli(c.n).RatString(); got != c.want {
			t.Errorf("Bernoulli(%d) = %s, want %s", c.n, got, c.want)
		}
	}
}

func TestSpecial(t *testing.T) {
	const tol = 1e-13
	agm, err := bf.AGM(bfv(1), bfv(2), testPrec)
	if err != nil {
		t.Fatal(err)
	}
	check(t, "AGM(1,2)", agm, 1.4567910310469068, tol)
	agm2, _ := bf.AGM(bfv(24), bfv(6), testPrec)
	check(t, "AGM(24,6)", agm2, 13.458171481725616, tol)

	ekCases := []struct{ k, want float64 }{
		{0, math.Pi / 2}, {0.5, 1.6857503548125961}, {0.9, 2.2805491384227703},
	}
	for _, c := range ekCases {
		ek, err := bf.EllipticK(bfv(c.k), testPrec)
		if err != nil {
			t.Fatal(err)
		}
		check(t, fmt.Sprintf("EllipticK(%g)", c.k), ek, c.want, tol)
	}
	eeCases := []struct{ k, want float64 }{
		{0, math.Pi / 2}, {0.5, 1.4674622093394272}, {0.9, 1.1716970527816142},
	}
	for _, c := range eeCases {
		ee, err := bf.EllipticE(bfv(c.k), testPrec)
		if err != nil {
			t.Fatal(err)
		}
		check(t, fmt.Sprintf("EllipticE(%g)", c.k), ee, c.want, tol)
	}
	lwCases := []struct{ x, want float64 }{
		{0, 0}, {1, 0.5671432904097838}, {math.E, 1}, {10, 1.7455280027406994},
		{-0.2, -0.2591711018190738}, {-0.3, -0.4894022271802149}, {100, 3.3856301402900502},
	}
	for _, c := range lwCases {
		w, err := bf.LambertW(bfv(c.x), testPrec)
		if err != nil {
			t.Fatalf("LambertW(%g): %v", c.x, err)
		}
		check(t, fmt.Sprintf("LambertW(%g)", c.x), w, c.want, tol)
	}
	check(t, "Sigmoid(0)", bf.Sigmoid(bfv(0), testPrec), 0.5, tol)
	check(t, "Sigmoid(2)", bf.Sigmoid(bfv(2), testPrec), 1/(1+math.Exp(-2)), tol)
	lg, _ := bf.Logit(bfv(0.75), testPrec)
	check(t, "Logit(0.75)", lg, math.Log(3), tol)
	check(t, "Softplus(0)", bf.Softplus(bfv(0), testPrec), math.Ln2, tol)
	check(t, "Softplus(10)", bf.Softplus(bfv(10), testPrec), math.Log(1+math.Exp(10)), tol)
	check(t, "LogAddExp", bf.LogAddExp(bfv(1000), bfv(1001), testPrec), 1001+math.Log1p(math.Exp(-1)), tol)
}

func TestUtil(t *testing.T) {
	const tol = 1e-14
	check(t, "Add", bf.Add(bfv(1.5), bfv(2.25), testPrec), 3.75, tol)
	check(t, "Sub", bf.Sub(bfv(1.5), bfv(2.25), testPrec), -0.75, tol)
	check(t, "Mul", bf.Mul(bfv(1.5), bfv(2), testPrec), 3, tol)
	check(t, "Div", bf.Div(bfv(3), bfv(4), testPrec), 0.75, tol)
	check(t, "Recip", bf.Recip(bfv(8), testPrec), 0.125, tol)
	check(t, "FMA", bf.FMA(bfv(2), bfv(3), bfv(4), testPrec), 10, tol)
	check(t, "Square", bf.Square(bfv(7), testPrec), 49, tol)
	check(t, "Cube", bf.Cube(bfv(3), testPrec), 27, tol)
	check(t, "Abs", bf.Abs(bfv(-4), testPrec), 4, tol)
	check(t, "Neg", bf.Neg(bfv(4), testPrec), -4, tol)
	check(t, "Dim", bf.Dim(bfv(5), bfv(2), testPrec), 3, tol)
	check(t, "Dim<0", bf.Dim(bfv(2), bfv(5), testPrec), 0, tol)
	check(t, "Min", bf.Min(bfv(2), bfv(5), testPrec), 2, tol)
	check(t, "Max", bf.Max(bfv(2), bfv(5), testPrec), 5, tol)
	check(t, "Clamp", bf.Clamp(bfv(9), bfv(0), bfv(4), testPrec), 4, tol)
	check(t, "Copysign", bf.Copysign(bfv(3), bfv(-1), testPrec), -3, tol)

	// Rounding.
	round := []struct {
		name string
		got  *big.Float
		want float64
	}{
		{"Floor(2.7)", bf.Floor(bfv(2.7), testPrec), 2},
		{"Floor(-2.3)", bf.Floor(bfv(-2.3), testPrec), -3},
		{"Ceil(2.3)", bf.Ceil(bfv(2.3), testPrec), 3},
		{"Ceil(-2.7)", bf.Ceil(bfv(-2.7), testPrec), -2},
		{"Trunc(2.9)", bf.Trunc(bfv(2.9), testPrec), 2},
		{"Trunc(-2.9)", bf.Trunc(bfv(-2.9), testPrec), -2},
		{"Round(2.5)", bf.Round(bfv(2.5), testPrec), 3},
		{"Round(-2.5)", bf.Round(bfv(-2.5), testPrec), -3},
		{"Round(2.4)", bf.Round(bfv(2.4), testPrec), 2},
		{"Frac(2.75)", bf.Frac(bfv(2.75), testPrec), 0.75},
		{"Mod(7,3)", bf.Mod(bfv(7), bfv(3), testPrec), 1},
		{"Remainder(7,3)", bf.Remainder(bfv(7), bfv(3), testPrec), 1},
		{"Ldexp(3,4)", bf.Ldexp(bfv(3), 4, testPrec), 48},
	}
	for _, c := range round {
		check(t, c.name, c.got, c.want, tol)
	}

	// Predicates and comparisons.
	if !bf.IsZero(bf.Zero(testPrec)) {
		t.Error("IsZero(0) false")
	}
	if !bf.IsInteger(bfv(4)) || bf.IsInteger(bfv(4.5)) {
		t.Error("IsInteger wrong")
	}
	if bf.Sign(bfv(-2)) != -1 || bf.Sign(bfv(3)) != 1 {
		t.Error("Sign wrong")
	}
	if bf.Cmp(bfv(1), bfv(2)) != -1 {
		t.Error("Cmp wrong")
	}
	if !bf.Equal(bfv(2), bfv(2)) {
		t.Error("Equal wrong")
	}
	if !bf.AlmostEqual(bfv(1), bfv(1.0000001), bfv(1e-3)) {
		t.Error("AlmostEqual wrong")
	}
	if !bf.Signbit(bfv(-1)) {
		t.Error("Signbit wrong")
	}
	frac, exp := bf.Frexp(bfv(12), testPrec)
	if exp != 4 || f64(frac) != 0.75 {
		t.Errorf("Frexp(12) = %g, %d", f64(frac), exp)
	}

	// Construction and conversion.
	p, err := bf.Parse("3.14159", testPrec)
	if err != nil {
		t.Fatal(err)
	}
	check(t, "Parse", p, 3.14159, 1e-6)
	if _, err := bf.Parse("not-a-number", testPrec); !errors.Is(err, bf.ErrParse) {
		t.Errorf("Parse bad input err = %v", err)
	}
	check(t, "FromInt", bf.FromInt(42, testPrec), 42, tol)
	check(t, "FromRat", bf.FromRat(big.NewRat(1, 4), testPrec), 0.25, tol)
	check(t, "FromBig", bf.FromBig(big.NewInt(99), testPrec), 99, tol)
	check(t, "MustParse", bf.MustParse("2.5", testPrec), 2.5, tol)
	if v, _ := bf.Float64(bfv(1.25)); v != 1.25 {
		t.Errorf("Float64 = %g", v)
	}
	if s := bf.String(bfv(1.5)); s == "" {
		t.Error("String empty")
	}
	if s := bf.Text(bfv(0.5), 'f', 2); s != "0.50" {
		t.Errorf("Text = %q", s)
	}
}

func TestDomainErrors(t *testing.T) {
	cases := []struct {
		name string
		err  error
		want error
	}{
		{"Log(0)", errOf(func() error { _, e := bf.Log(bfv(0), testPrec); return e }), bf.ErrDomain},
		{"Log(-1)", errOf(func() error { _, e := bf.Log(bfv(-1), testPrec); return e }), bf.ErrDomain},
		{"Sqrt(-1)", errOf(func() error { _, e := bf.Sqrt(bfv(-1), testPrec); return e }), bf.ErrNegative},
		{"Asin(2)", errOf(func() error { _, e := bf.Asin(bfv(2), testPrec); return e }), bf.ErrDomain},
		{"Acos(-2)", errOf(func() error { _, e := bf.Acos(bfv(-2), testPrec); return e }), bf.ErrDomain},
		{"Acosh(0.5)", errOf(func() error { _, e := bf.Acosh(bfv(0.5), testPrec); return e }), bf.ErrDomain},
		{"Atanh(1)", errOf(func() error { _, e := bf.Atanh(bfv(1), testPrec); return e }), bf.ErrDomain},
		{"Gamma(0)", errOf(func() error { _, e := bf.Gamma(bfv(0), testPrec); return e }), bf.ErrPole},
		{"Gamma(-3)", errOf(func() error { _, e := bf.Gamma(bfv(-3), testPrec); return e }), bf.ErrPole},
		{"Digamma(-2)", errOf(func() error { _, e := bf.Digamma(bfv(-2), testPrec); return e }), bf.ErrPole},
		{"Beta(-1,2)", errOf(func() error { _, e := bf.Beta(bfv(-1), bfv(2), testPrec); return e }), bf.ErrDomain},
		{"EllipticK(1)", errOf(func() error { _, e := bf.EllipticK(bfv(1), testPrec); return e }), bf.ErrDomain},
		{"LambertW(-1)", errOf(func() error { _, e := bf.LambertW(bfv(-1), testPrec); return e }), bf.ErrDomain},
		{"Pow(0,-1)", errOf(func() error { _, e := bf.Pow(bfv(0), bfv(-1), testPrec); return e }), bf.ErrDomain},
		{"Root(-4,2)", errOf(func() error { _, e := bf.Root(bfv(-4), 2, testPrec); return e }), bf.ErrDomain},
		{"Logit(1.5)", errOf(func() error { _, e := bf.Logit(bfv(1.5), testPrec); return e }), bf.ErrDomain},
	}
	for _, c := range cases {
		if !errors.Is(c.err, c.want) {
			t.Errorf("%s: err = %v, want %v", c.name, c.err, c.want)
		}
	}
}

func errOf(f func() error) error { return f() }

// TestHighPrecision confirms the constants agree with reference digit strings
// far beyond float64 precision.
func TestHighPrecision(t *testing.T) {
	const p = 400
	refs := []struct{ name, ref string }{
		{"Pi", "3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844609550582231725359408128481117450284102701938521105559644622948954930381964"},
		{"E", "2.718281828459045235360287471352662497757247093699959574966967627724076630353547594571382178525166427427466391932003059921817413596629043572900334295260595630738132328627943490763233829880753195251019011"},
		{"Ln2", "0.693147180559945309417232121458176568075500134360255254120680009493393621969694715605863326996418687542001481020570685733685520235758130557032670751635075961930727570828371435190307038623891673471123350"},
	}
	for _, r := range refs {
		want, _, err := big.ParseFloat(r.ref, 10, p+64, big.ToNearestEven)
		if err != nil {
			t.Fatal(err)
		}
		var got *big.Float
		switch r.name {
		case "Pi":
			got = bf.Pi(p)
		case "E":
			got = bf.E(p)
		case "Ln2":
			got = bf.Ln2(p)
		}
		d := new(big.Float).SetPrec(p+64).Sub(got, want)
		d.Abs(d)
		// Require agreement to at least p-8 bits.
		tol := new(big.Float).SetMantExp(big.NewFloat(1), -(p - 8))
		if d.Cmp(tol) > 0 {
			t.Errorf("%s: only agrees to exponent %d bits", r.name, d.MantExp(nil))
		}
	}
}

func TestPrecisionCarried(t *testing.T) {
	// Results should carry exactly the requested precision.
	for _, p := range []uint{53, 100, 256, 500} {
		if got := bf.Pi(p).Prec(); got != p {
			t.Errorf("Pi(%d).Prec() = %d", p, got)
		}
		if got := bf.Exp(bfv(1), p).Prec(); got != p {
			t.Errorf("Exp precision = %d, want %d", got, p)
		}
	}
}

// ExampleExp computes e to 50 decimal digits by evaluating exp(1) at 200 bits
// of precision.
func ExampleExp() {
	one := bf.FromInt(1, 200)
	e := bf.Exp(one, 200)
	fmt.Println(e.Text('f', 50))
	// Output: 2.71828182845904523536028747135266249775724709369996
}

// ExamplePi prints pi to 30 decimal places.
func ExamplePi() {
	fmt.Println(bf.Pi(160).Text('f', 30))
	// Output: 3.141592653589793238462643383280
}

// ExampleLambertW solves w*exp(w) = 1 for w, the omega constant.
func ExampleLambertW() {
	w, _ := bf.LambertW(bf.FromInt(1, 120), 120)
	fmt.Println(w.Text('f', 20))
	// Output: 0.56714329040978387300
}
