package padic

import (
	"math/big"
	"testing"
)

func p(n int64) *big.Int { return big.NewInt(n) }

func must(x *Padic, err error) *Padic {
	if err != nil {
		panic(err)
	}
	return x
}

func TestValuationInt(t *testing.T) {
	tests := []struct {
		p, n int64
		want int
	}{
		{2, 8, 3},
		{2, 12, 2},
		{3, 54, 3},
		{5, 7, 0},
		{5, 0, -1},
		{7, -49, 2},
	}
	for _, tc := range tests {
		if got := ValuationInt(p(tc.p), p(tc.n)); got != tc.want {
			t.Errorf("ValuationInt(%d,%d)=%d want %d", tc.p, tc.n, got, tc.want)
		}
	}
}

func TestUnitPartInt(t *testing.T) {
	tests := []struct {
		p, n, want int64
	}{
		{2, 8, 1},
		{2, 12, 3},
		{3, 54, 2},
		{5, 7, 7},
		{7, -49, -1},
	}
	for _, tc := range tests {
		if got := UnitPartInt(p(tc.p), p(tc.n)); got.Int64() != tc.want {
			t.Errorf("UnitPartInt(%d,%d)=%v want %d", tc.p, tc.n, got, tc.want)
		}
	}
}

func TestArithmetic(t *testing.T) {
	pr := p(5)
	prec := 6
	// 1/2 + 1/3 = 5/6
	a := must(FromRational(pr, p(1), p(2), prec))
	b := must(FromRational(pr, p(1), p(3), prec))
	sum := must(mustErr(a.Add(b)))
	expect := must(FromRational(pr, p(5), p(6), prec))
	if !sum.Equal(expect) {
		t.Errorf("1/2+1/3: got %v want %v", sum, expect)
	}
	// (1/2)*(2) = 1
	two := FromInt(pr, 2, prec)
	prod := must(mustErr(a.Mul(two)))
	if !prod.IsOne() {
		t.Errorf("(1/2)*2 = %v, want 1", prod)
	}
	// inverse of 3 times 3 == 1
	three := FromInt(pr, 3, prec)
	inv := must(three.Inv())
	back := must(mustErr(inv.Mul(three)))
	if !back.IsOne() {
		t.Errorf("3^-1 * 3 = %v want 1", back)
	}
	// subtraction cancellation: 7 - 7 == 0
	seven := FromInt(pr, 7, prec)
	zero := must(mustErr(seven.Sub(seven)))
	if !zero.IsZero() {
		t.Errorf("7-7 not zero: %v", zero)
	}
	// pow: 2^10 = 1024
	pw := must(FromInt(pr, 2, prec).Pow(10))
	e := FromInt(pr, 1024, prec)
	if !pw.Equal(e) {
		t.Errorf("2^10 = %v want %v", pw, e)
	}
	// negative power: 2^-1 * 2 = 1
	np := must(FromInt(pr, 2, prec).Pow(-1))
	if !must(mustErr(np.Mul(two))).IsOne() {
		t.Errorf("2^-1 * 2 != 1")
	}
}

func mustErr(x *Padic, err error) (*Padic, error) { return x, err }

func TestValuationOfProducts(t *testing.T) {
	pr := p(3)
	prec := 5
	// val(9) = 2, val(1/3) = -1
	nine := FromInt(pr, 9, prec)
	if nine.Valuation() != 2 {
		t.Errorf("val(9 base 3)=%d want 2", nine.Valuation())
	}
	oneThird := must(FromRational(pr, p(1), p(3), prec))
	if oneThird.Valuation() != -1 {
		t.Errorf("val(1/3)=%d want -1", oneThird.Valuation())
	}
	// product val adds
	prod := must(mustErr(nine.Mul(oneThird)))
	if prod.Valuation() != 1 {
		t.Errorf("val(9 * 1/3)=%d want 1", prod.Valuation())
	}
	// abs value |9|_3 = 1/9
	if prod2 := nine.AbsValue(); prod2.Cmp(new(big.Rat).SetFrac(p(1), p(9))) != 0 {
		t.Errorf("|9|_3 = %v want 1/9", prod2)
	}
}

func TestTeichmuller(t *testing.T) {
	tests := []struct {
		p, a int64
		prec int
		want int64
	}{
		{5, 2, 3, 57},
		{5, 1, 3, 1},
		{7, 1, 2, 1},
	}
	for _, tc := range tests {
		got, err := TeichmullerRep(p(tc.p), p(tc.a), tc.prec)
		if err != nil {
			t.Fatalf("TeichmullerRep err: %v", err)
		}
		if got.Int64() != tc.want {
			t.Errorf("TeichmullerRep(%d,%d,%d)=%v want %d", tc.p, tc.a, tc.prec, got, tc.want)
		}
		// verify omega^(p-1) == 1 and omega == a mod p
		if !IsTeichmuller(p(tc.p), got, tc.prec) {
			t.Errorf("result not Teichmuller: %v", got)
		}
	}
	// Teichmuller of a unit powered to (p-1) is 1.
	om := must(Teichmuller(p(5), p(2), 5))
	p1 := must(om.Pow(4)) // p-1 = 4
	if !p1.IsOne() {
		t.Errorf("omega^(p-1) = %v want 1", p1)
	}
}

func TestSqrt(t *testing.T) {
	// sqrt(2) in Q_7
	s := must(FromInt(p(7), 2, 6).Sqrt())
	if !s.Square().Equal(FromInt(p(7), 2, 6)) {
		t.Errorf("sqrt(2)^2 != 2 in Q7")
	}
	// 2 is a square in Q7
	if !FromInt(p(7), 2, 6).IsSquare() {
		t.Errorf("2 should be a square in Q7")
	}
	// 3 is not a square in Q7 (3 is a non-residue mod 7)
	if FromInt(p(7), 3, 6).IsSquare() {
		t.Errorf("3 should not be a square in Q7")
	}
	if _, err := FromInt(p(7), 3, 6).Sqrt(); err == nil {
		t.Errorf("expected error for sqrt(3) in Q7")
	}
	// p=2: 17 is a square (17 == 1 mod 8)
	s2 := must(FromInt(p(2), 17, 8).Sqrt())
	if !s2.Square().Equal(FromInt(p(2), 17, 8)) {
		t.Errorf("sqrt(17)^2 != 17 in Q2")
	}
	// p=2: 3 is not a square (3 mod 8 != 1)
	if FromInt(p(2), 3, 8).IsSquare() {
		t.Errorf("3 should not be square in Q2")
	}
	// odd valuation is not a square
	if FromInt(p(3), 3, 6).IsSquare() {
		t.Errorf("3 (val 1) should not be a square in Q3")
	}
	// both roots
	roots, err := FromInt(p(7), 2, 6).SqrtBothRoots()
	if err != nil || len(roots) != 2 {
		t.Fatalf("SqrtBothRoots failed: %v", err)
	}
	if roots[0].Equal(roots[1]) {
		t.Errorf("both roots identical")
	}
}

func TestSqrtInt(t *testing.T) {
	r, err := SqrtInt(p(7), p(2), 5)
	if err != nil {
		t.Fatal(err)
	}
	sq := new(big.Int).Mod(new(big.Int).Mul(r, r), PPow(p(7), 5))
	if sq.Cmp(p(2)) != 0 {
		t.Errorf("SqrtInt(2)^2 mod 7^5 = %v want 2", sq)
	}
}

func TestHensel(t *testing.T) {
	// x^2 - 2 over 7, root near 3
	f := []*big.Int{p(-2), p(0), p(1)}
	root, err := HenselLift(f, p(3), p(7), 5)
	if err != nil {
		t.Fatal(err)
	}
	val := PolyEval(f, root)
	val.Mod(val, PPow(p(7), 5))
	if val.Sign() != 0 {
		t.Errorf("Hensel root does not satisfy f: f(%v)=%v", root, val)
	}
	// simple roots mod p
	simple := SimpleRootsModP(f, p(7))
	if len(simple) != 2 {
		t.Errorf("expected 2 simple roots of x^2-2 mod 7, got %d", len(simple))
	}
	// PadicRoots returns lifted roots
	roots, err := PadicRoots(f, p(7), 5)
	if err != nil || len(roots) != 2 {
		t.Fatalf("PadicRoots failed: %v (%d)", err, len(roots))
	}
	for _, r := range roots {
		fr := must(EvalPolyPadic(f, r))
		if !fr.IsZero() {
			t.Errorf("PadicRoots root not a root: %v", fr)
		}
	}
}

func TestExpLog(t *testing.T) {
	pr := p(5)
	prec := 6
	// log(1+5) then exp gives back 1+5
	u := FromInt(pr, 6, prec)
	lg, err := u.Log()
	if err != nil {
		t.Fatal(err)
	}
	ex, err := lg.Exp()
	if err != nil {
		t.Fatal(err)
	}
	if !ex.Equal(u) {
		t.Errorf("exp(log(6)) = %v want 6", ex)
	}
	// homomorphism log(ab) = log a + log b
	a := FromInt(pr, 11, prec)
	b := FromInt(pr, 16, prec)
	la := must(a.Log())
	lb := must(b.Log())
	ab := must(mustErr(a.Mul(b)))
	lab := must(ab.Log())
	sumlog := must(mustErr(la.Add(lb)))
	if !lab.IsCloseTo(sumlog, prec-1) {
		t.Errorf("log(ab) != log a + log b: %v vs %v", lab, sumlog)
	}
	// domain errors
	if _, err := FromInt(pr, 2, prec).Log(); err == nil {
		t.Errorf("log of non-1-unit should error")
	}
	if _, err := FromInt(p(2), 2, 6).Exp(); err == nil {
		t.Errorf("exp with val 1 in Q2 should error")
	}
	// exp/log convergence flags
	if !must(FromInt(pr, 5, prec).Add(FromInt(pr, 1, prec))).LogConverges() {
		t.Errorf("1+5 should have converging log")
	}
	if !FromInt(pr, 5, prec).ExpConverges() {
		t.Errorf("5 should have converging exp in Q5")
	}
}

func TestNewtonPolygon(t *testing.T) {
	// x^2 - 6 over p=3: roots have valuation 1/2 each
	np := NewtonPolygonFromInts(p(3), []*big.Int{p(-6), p(0), p(1)})
	if len(np.Vertices) != 2 {
		t.Fatalf("expected 2 vertices, got %v", np.Vertices)
	}
	slopes := np.Slopes()
	if len(slopes) != 1 || slopes[0].Cmp(new(big.Rat).SetFrac(p(-1), p(2))) != 0 {
		t.Errorf("slope = %v want -1/2", slopes)
	}
	rv := np.RootValuations()
	if len(rv) != 2 {
		t.Fatalf("expected 2 root valuations, got %v", rv)
	}
	half := new(big.Rat).SetFrac(p(1), p(2))
	for _, v := range rv {
		if v.Cmp(half) != 0 {
			t.Errorf("root valuation = %v want 1/2", v)
		}
	}
	// Eisenstein-like: x^2 + 3x + 3 over p=3 -> pure, slope -1/2? valuations
	// a0=1(v?) actually val(3)=1,val(3)=1,val(1)=0: points (0,1),(1,1),(2,0)
	np2 := NewtonPolygonFromInts(p(3), []*big.Int{p(3), p(3), p(1)})
	if !np2.IsPure() {
		t.Errorf("x^2+3x+3 should be pure Newton polygon, vertices=%v", np2.Vertices)
	}
	// two distinct slopes: x^2 + 9x + 3? points (0,1),(1,2),(2,0): hull (0,1)-(2,0)
	np3 := NewtonPolygonFromInts(p(3), []*big.Int{p(9), p(3), p(1)})
	// points (0,2),(1,1),(2,0): collinear slope -1 single segment
	if len(np3.Slopes()) != 1 {
		t.Errorf("expected 1 slope, got %v", np3.Slopes())
	}
}

func TestStrassmann(t *testing.T) {
	tests := []struct {
		coeffs []int64
		want   int
	}{
		{[]int64{5, 1, 25, 5}, 1},
		{[]int64{1, 5, 25}, 0},
		{[]int64{5, 5, 1}, 2},
	}
	for _, tc := range tests {
		cs := make([]*big.Int, len(tc.coeffs))
		for i, c := range tc.coeffs {
			cs[i] = p(c)
		}
		got, err := StrassmannBoundInts(p(5), cs)
		if err != nil {
			t.Fatal(err)
		}
		if got != tc.want {
			t.Errorf("StrassmannBound(%v)=%d want %d", tc.coeffs, got, tc.want)
		}
	}
}

func TestExpansion(t *testing.T) {
	// 1/2 in Q_3
	v, digits, err := ExpandRational(p(3), p(1), p(2), 6)
	if err != nil {
		t.Fatal(err)
	}
	if v != 0 {
		t.Errorf("val = %d want 0", v)
	}
	want := []int64{2, 1, 1, 1, 1, 1}
	for i, d := range digits {
		if d.Int64() != want[i] {
			t.Errorf("digit %d = %v want %d", i, d, want[i])
		}
	}
	// reconstruction is congruent to 1/2 mod 3^6
	r := DigitsToRat(p(3), v, digits)
	// r should equal 1/2 mod 3^6, i.e. 2*r - 1 divisible by 3^6
	check := new(big.Int).Sub(new(big.Int).Mul(r.Num(), p(2)), r.Denom())
	check.Mul(check, p(1))
	num := new(big.Int).Mul(r.Num(), p(2))
	num.Sub(num, r.Denom())
	if new(big.Int).Mod(num, PPow(p(3), 6)).Sign() != 0 {
		t.Errorf("reconstruction not congruent to 1/2 mod 3^6")
	}
	// DigitsToPadic round-trips digits
	x := DigitsToPadic(p(3), v, digits)
	orig := must(FromRational(p(3), p(1), p(2), 6))
	if !x.Equal(orig) {
		t.Errorf("DigitsToPadic = %v want %v", x, orig)
	}
}

func TestRationalReconstruction(t *testing.T) {
	pr := p(7)
	orig := must(FromRational(pr, p(3), p(5), 8))
	num, den, err := orig.RationalReconstruction()
	if err != nil {
		t.Fatal(err)
	}
	if num.Int64() != 3 || den.Int64() != 5 {
		t.Errorf("reconstruction = %v/%v want 3/5", num, den)
	}
}

func TestNumberTheory(t *testing.T) {
	if !IsPrime(p(97)) {
		t.Errorf("97 is prime")
	}
	if IsPrime(p(91)) {
		t.Errorf("91 = 7*13 not prime")
	}
	if NextPrime(p(90)).Int64() != 97 {
		t.Errorf("NextPrime(90) = %v want 97", NextPrime(p(90)))
	}
	if LegendreSymbol(p(2), p(7)) != 1 {
		t.Errorf("(2/7) should be 1")
	}
	if LegendreSymbol(p(3), p(7)) != -1 {
		t.Errorf("(3/7) should be -1")
	}
	if JacobiSymbol(p(2), p(15)) != 1 {
		t.Errorf("(2/15) should be 1")
	}
	// SqrtModP
	r, err := SqrtModP(p(2), p(7), nil)
	if err != nil {
		t.Fatal(err)
	}
	if new(big.Int).Mod(new(big.Int).Mul(r, r), p(7)).Cmp(p(2)) != 0 {
		t.Errorf("SqrtModP(2,7)^2 != 2")
	}
	// p = 1 mod 4 case (Tonelli-Shanks main branch): p=13, a=4->2
	r2, err := SqrtModP(p(10), p(13), nil)
	if err != nil {
		t.Fatal(err)
	}
	if new(big.Int).Mod(new(big.Int).Mul(r2, r2), p(13)).Cmp(p(10)) != 0 {
		t.Errorf("SqrtModP(10,13) wrong")
	}
	// InvMod / PowMod
	inv, _ := InvMod(p(3), p(7))
	if inv.Int64() != 5 {
		t.Errorf("3^-1 mod 7 = %v want 5", inv)
	}
	if PowMod(p(2), p(10), p(1000)).Int64() != 24 {
		t.Errorf("2^10 mod 1000 = %v want 24", PowMod(p(2), p(10), p(1000)))
	}
}

func TestErrors(t *testing.T) {
	if _, err := FromRational(p(5), p(1), p(0), 5); err == nil {
		t.Errorf("expected zero-division error")
	}
	if _, err := Zero(p(5), 5).Inv(); err == nil {
		t.Errorf("expected error inverting zero")
	}
	if _, err := New(p(4), 0, p(1), 5); err == nil {
		t.Errorf("expected not-prime error for New with p=4")
	}
	a := FromInt(p(3), 1, 5)
	b := FromInt(p(5), 1, 5)
	if _, err := a.Add(b); err == nil {
		t.Errorf("expected prime-mismatch error")
	}
}

func TestPrecisionTracking(t *testing.T) {
	pr := p(5)
	a := FromInt(pr, 1, 10)
	b := FromInt(pr, 1, 4)
	// adding lower precision reduces result precision
	sum := must(mustErr(a.Add(b)))
	if sum.AbsolutePrecision() != 4 {
		t.Errorf("abs precision = %d want 4", sum.AbsolutePrecision())
	}
	// ReduceTo cannot increase precision
	red := b.ReduceTo(20)
	if red.AbsolutePrecision() != 4 {
		t.Errorf("ReduceTo should not increase precision, got %d", red.AbsolutePrecision())
	}
	red2 := a.ReduceTo(3)
	if red2.AbsolutePrecision() != 3 {
		t.Errorf("ReduceTo(3) = %d", red2.AbsolutePrecision())
	}
}

func TestPolyHelpers(t *testing.T) {
	a := []*big.Int{p(1), p(2), p(3)} // 1 + 2x + 3x^2
	b := []*big.Int{p(0), p(1)}       // x
	sum := PolyAdd(a, b)              // 1 + 3x + 3x^2
	if PolyDegree(sum) != 2 || sum[1].Int64() != 3 {
		t.Errorf("PolyAdd wrong: %v", sum)
	}
	prod := PolyMul(b, b) // x^2
	if PolyDegree(prod) != 2 || prod[2].Int64() != 1 {
		t.Errorf("PolyMul wrong: %v", prod)
	}
	if got := PolyEval(a, p(2)); got.Int64() != 1+4+12 {
		t.Errorf("PolyEval = %v want 17", got)
	}
	// Eisenstein: x^2 + 3x + 3 over p=3
	if !IsEisenstein(p(3), []*big.Int{p(3), p(3), p(1)}) {
		t.Errorf("x^2+3x+3 should be Eisenstein at 3")
	}
	// not Eisenstein: constant term divisible by 9
	if IsEisenstein(p(3), []*big.Int{p(9), p(3), p(1)}) {
		t.Errorf("constant div by 9 is not Eisenstein")
	}
	v, ok := PolyValuation(p(3), []*big.Int{p(9), p(6), p(3)})
	if !ok || v != 1 {
		t.Errorf("PolyValuation = %d,%v want 1,true", v, ok)
	}
}

func TestIntUtils(t *testing.T) {
	if Factorial(5).Int64() != 120 {
		t.Errorf("5! = %v", Factorial(5))
	}
	if ValuationFactorial(p(3), 9) != 4 {
		// v_3(9!) = 3 + 1 = 4
		t.Errorf("v_3(9!) = %d want 4", ValuationFactorial(p(3), 9))
	}
	if Binomial(6, 2).Int64() != 15 {
		t.Errorf("C(6,2) = %v", Binomial(6, 2))
	}
	if MultiplicativeOrderModP(p(2), p(7)) != 3 {
		// 2^3 = 8 = 1 mod 7
		t.Errorf("ord_7(2) = %d want 3", MultiplicativeOrderModP(p(2), p(7)))
	}
	x, m, err := CRTPair(p(2), p(3), p(3), p(5))
	if err != nil || m.Int64() != 15 || x.Int64() != 8 {
		t.Errorf("CRT = %v mod %v want 8 mod 15", x, m)
	}
	pp, k, ok := IsPrimePower(p(27))
	if !ok || pp.Int64() != 3 || k != 3 {
		t.Errorf("IsPrimePower(27) = %v^%d,%v", pp, k, ok)
	}
	if _, _, ok := IsPrimePower(p(12)); ok {
		t.Errorf("12 is not a prime power")
	}
	if ValuationInt64(2, 24) != 3 {
		t.Errorf("v_2(24) = %d want 3", ValuationInt64(2, 24))
	}
}

func TestPadicMethods(t *testing.T) {
	pr := p(5)
	x := FromInt(pr, 7, 6) // 7, val 0
	// MulPow shifts valuation
	y := x.MulPow(2)
	if y.Valuation() != 2 {
		t.Errorf("MulPow val = %d want 2", y.Valuation())
	}
	// AddInt/MulInt
	if !x.AddInt(3).Equal(FromInt(pr, 10, 6)) {
		t.Errorf("7+3 != 10")
	}
	if !x.MulInt(2).Equal(FromInt(pr, 14, 6)) {
		t.Errorf("7*2 != 14")
	}
	// Distance
	d := FromInt(pr, 1, 6).Distance(FromInt(pr, 26, 6)) // 1-26 = -25, val 2 -> |.|=1/25
	if d.Cmp(new(big.Rat).SetFrac(p(1), p(25))) != 0 {
		t.Errorf("distance = %v want 1/25", d)
	}
	// ConstantResidue
	cr, _ := x.ConstantResidue()
	if cr.Int64() != 2 { // 7 mod 5 = 2
		t.Errorf("ConstantResidue = %v want 2", cr)
	}
	// Digit: 7 = 2 + 1*5, so digit0=2, digit1=1
	d0, _ := x.Digit(0)
	d1, _ := x.Digit(1)
	if d0.Int64() != 2 || d1.Int64() != 1 {
		t.Errorf("digits = %v,%v want 2,1", d0, d1)
	}
	// IsRootOfUnity: Teichmuller rep of 2 in Q5
	om := must(Teichmuller(pr, p(2), 6))
	if !om.IsRootOfUnity() {
		t.Errorf("Teichmuller rep should be a root of unity")
	}
	if x.IsRootOfUnity() {
		t.Errorf("7 is not a root of unity in Q5")
	}
	// FromString
	fs := must(FromString(pr, "3/4", 6))
	if !fs.Equal(must(FromRational(pr, p(3), p(4), 6))) {
		t.Errorf("FromString 3/4 mismatch")
	}
}

func TestNewtonPolygonExtras(t *testing.T) {
	// (x - p)(x - p^2) over p=3 => x^2 - (p+p^2) x + p^3: coeffs [27,-12,1]
	np := NewtonPolygonFromInts(p(3), []*big.Int{p(27), p(-12), p(1)})
	// valuations: (0,3),(1,1),(2,0) -> hull vertices all three, slopes -2,-1
	rv := np.RootValuations()
	if len(rv) != 2 {
		t.Fatalf("expected 2 roots, got %v", rv)
	}
	// root valuations should be 1 and 2
	got := map[string]bool{rv[0].RatString(): true, rv[1].RatString(): true}
	if !got["1"] || !got["2"] {
		t.Errorf("root valuations = %v want {1,2}", rv)
	}
	if np.NumRoots() != 2 {
		t.Errorf("NumRoots = %d want 2", np.NumRoots())
	}
	if np.Width() != 2 {
		t.Errorf("Width = %d", np.Width())
	}
}

func TestRootOfUnity(t *testing.T) {
	// omega^4 = 1 in Q_5 for a primitive root
	om, err := RootOfUnity(p(5), 6)
	if err != nil {
		t.Fatal(err)
	}
	if !must(om.Pow(4)).IsOne() {
		t.Errorf("root of unity^4 != 1")
	}
	// order of the underlying primitive root mod 5 is 4
	g := PrimitiveRootModP(p(5))
	if MultiplicativeOrderModP(g, p(5)) != 4 {
		t.Errorf("primitive root order wrong")
	}
}
