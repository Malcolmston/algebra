package ellipticcurves

import (
	"fmt"
	"math"
	"math/big"
	"math/rand"
	"testing"
)

func bi(n int64) *big.Int { return big.NewInt(n) }
func rat(n, d int64) *big.Rat {
	return new(big.Rat).SetFrac(big.NewInt(n), big.NewInt(d))
}

// ---------------------------------------------------------------------------
// Field arithmetic
// ---------------------------------------------------------------------------

func TestModularArithmetic(t *testing.T) {
	p := bi(97)
	tests := []struct {
		name string
		got  *big.Int
		want int64
	}{
		{"add", ModAdd(bi(90), bi(20), p), 13},
		{"sub", ModSub(bi(5), bi(20), p), 82},
		{"mul", ModMul(bi(50), bi(50), p), 2500 % 97},
		{"neg", ModNeg(bi(10), p), 87},
		{"double", ModDouble(bi(60), p), 120 % 97},
		{"square", ModSquare(bi(13), p), 169 % 97},
		{"mod-negative", Mod(bi(-3), p), 94},
	}
	for _, tc := range tests {
		if tc.got.Cmp(bi(tc.want)) != 0 {
			t.Errorf("%s: got %s want %d", tc.name, tc.got, tc.want)
		}
	}
	inv, err := ModInverse(bi(13), p)
	if err != nil {
		t.Fatalf("inverse: %v", err)
	}
	if ModMul(bi(13), inv, p).Cmp(bigOne) != 0 {
		t.Errorf("13 * 13^-1 != 1 mod 97")
	}
	if _, err := ModInverse(bi(0), p); err == nil {
		t.Errorf("expected error inverting 0")
	}
	e, err := ModExp(bi(5), bi(-1), p)
	if err != nil {
		t.Fatalf("modexp neg: %v", err)
	}
	if ModMul(bi(5), e, p).Cmp(bigOne) != 0 {
		t.Errorf("5 * 5^-1 != 1")
	}
}

func TestSquareRootsAndSymbols(t *testing.T) {
	primes := []int64{7, 11, 13, 101, 1009}
	for _, pv := range primes {
		p := bi(pv)
		for a := int64(1); a < pv; a++ {
			leg := LegendreSymbol(bi(a), p)
			jac := JacobiSymbol(bi(a), p)
			kro := KroneckerSymbol(bi(a), p)
			if leg != jac || leg != kro {
				t.Fatalf("symbol mismatch a=%d p=%d: leg=%d jac=%d kro=%d", a, pv, leg, jac, kro)
			}
			if leg == 1 {
				r, err := ModSqrt(bi(a), p)
				if err != nil {
					t.Fatalf("sqrt %d mod %d: %v", a, pv, err)
				}
				if ModSquare(r, p).Cmp(bi(a)) != 0 {
					t.Errorf("bad sqrt %d mod %d", a, pv)
				}
				rt, err := ModSqrtTonelli(bi(a), p)
				if err != nil || ModSquare(rt, p).Cmp(bi(a)) != 0 {
					t.Errorf("tonelli bad sqrt %d mod %d", a, pv)
				}
			}
		}
	}
}

func TestJacobiKroneckerKnown(t *testing.T) {
	tests := []struct {
		a, n int64
		jac  int
		kro  int
	}{
		{2, 15, 1, 1},
		{7, 15, -1, -1},
		{1001, 9907, -1, -1},
		{2, 7, 1, 1},
	}
	for _, tc := range tests {
		if g := JacobiSymbol(bi(tc.a), bi(tc.n)); g != tc.jac {
			t.Errorf("Jacobi(%d,%d)=%d want %d", tc.a, tc.n, g, tc.jac)
		}
		if g := KroneckerSymbol(bi(tc.a), bi(tc.n)); g != tc.kro {
			t.Errorf("Kronecker(%d,%d)=%d want %d", tc.a, tc.n, g, tc.kro)
		}
	}
	if KroneckerSymbol(bi(-1), bi(5)) != 1 {
		t.Errorf("Kronecker(-1,5) want 1")
	}
	if KroneckerSymbol(bi(-1), bi(7)) != -1 {
		t.Errorf("Kronecker(-1,7) want -1")
	}
}

func TestCRT(t *testing.T) {
	x, m, err := CRT([]*big.Int{bi(2), bi(3), bi(2)}, []*big.Int{bi(3), bi(5), bi(7)})
	if err != nil {
		t.Fatalf("CRT: %v", err)
	}
	if x.Cmp(bi(23)) != 0 || m.Cmp(bi(105)) != 0 {
		t.Errorf("CRT got (%s,%s) want (23,105)", x, m)
	}
	for _, mod := range []int64{3, 5, 7} {
		if new(big.Int).Mod(x, bi(mod)).Cmp(bi(map[int64]int64{3: 2, 5: 3, 7: 2}[mod])) != 0 {
			t.Errorf("CRT residue wrong mod %d", mod)
		}
	}
}

func TestFactorization(t *testing.T) {
	tests := []struct {
		n        int64
		divisors int
		totient  int64
	}{
		{12, 6, 4},
		{97, 2, 96},
		{360, 24, 96},
		{1, 1, 1},
	}
	for _, tc := range tests {
		if d := NumDivisors(bi(tc.n)); d != tc.divisors {
			t.Errorf("NumDivisors(%d)=%d want %d", tc.n, d, tc.divisors)
		}
		if len(Divisors(bi(tc.n))) != tc.divisors {
			t.Errorf("len Divisors(%d) mismatch", tc.n)
		}
		if e := EulerTotient(bi(tc.n)); e.Cmp(bi(tc.totient)) != 0 {
			t.Errorf("EulerTotient(%d)=%s want %d", tc.n, e, tc.totient)
		}
	}
	// Factorization reconstructs n.
	for _, n := range []int64{360, 1001, 999983 * 2, 65537} {
		prod := big.NewInt(1)
		for prime, exp := range Factorize(bi(n)) {
			if !prime.ProbablyPrime(20) {
				t.Errorf("non-prime factor %s of %d", prime, n)
			}
			prod.Mul(prod, new(big.Int).Exp(prime, bi(int64(exp)), nil))
		}
		if prod.Cmp(bi(n)) != 0 {
			t.Errorf("Factorize(%d) product %s", n, prod)
		}
	}
}

// ---------------------------------------------------------------------------
// Group law over F_p
// ---------------------------------------------------------------------------

func TestGroupLawFp(t *testing.T) {
	c, err := NewCurveFp(bi(2), bi(3), bi(97))
	if err != nil {
		t.Fatal(err)
	}
	pts := c.EnumeratePoints()
	inf := PointAtInfinityFp()

	// Identity, inverse, and on-curve for all points.
	for _, p := range pts {
		if !c.IsOnCurve(p) {
			t.Fatalf("point %s not on curve", p)
		}
		if !c.PointEqual(c.Add(p, inf), p) || !c.PointEqual(c.Add(inf, p), p) {
			t.Errorf("identity law failed for %s", p)
		}
		if !c.Add(p, c.Neg(p)).Infinity {
			t.Errorf("inverse law failed for %s", p)
		}
	}

	// Associativity and commutativity on a sample.
	rng := rand.New(rand.NewSource(11))
	for i := 0; i < 200; i++ {
		a := pts[rng.Intn(len(pts))]
		b := pts[rng.Intn(len(pts))]
		d := pts[rng.Intn(len(pts))]
		if !c.PointEqual(c.Add(a, b), c.Add(b, a)) {
			t.Fatalf("commutativity failed")
		}
		left := c.Add(c.Add(a, b), d)
		right := c.Add(a, c.Add(b, d))
		if !c.PointEqual(left, right) {
			t.Fatalf("associativity failed: %s %s %s", a, b, d)
		}
	}

	// Doubling matches Add(p,p), and ScalarMul matches repeated addition.
	for _, p := range pts[:iMin(10, len(pts))] {
		if !c.PointEqual(c.Double(p), c.Add(p, p)) {
			t.Errorf("Double != Add(p,p) for %s", p)
		}
		acc := inf
		for k := int64(0); k <= 12; k++ {
			if !c.PointEqual(c.ScalarMul(bi(k), p), acc) {
				t.Errorf("ScalarMul mismatch k=%d p=%s", k, p)
			}
			acc = c.Add(acc, p)
		}
	}
}

func TestScalarMulNegativeAndOrder(t *testing.T) {
	c, _ := NewCurveFp(bi(2), bi(3), bi(97))
	order := c.CurveOrderNaive()
	rng := rand.New(rand.NewSource(5))
	for i := 0; i < 20; i++ {
		p := c.RandomPointFp(rng)
		// order*p = O
		if !c.ScalarMul(order, p).Infinity {
			t.Errorf("order*p != O")
		}
		// (-k)p = -(kp)
		k := bi(int64(rng.Intn(50) + 1))
		if !c.PointEqual(c.ScalarMul(new(big.Int).Neg(k), p), c.Neg(c.ScalarMul(k, p))) {
			t.Errorf("negative scalar mismatch")
		}
	}
}

// ---------------------------------------------------------------------------
// Invariants
// ---------------------------------------------------------------------------

func TestCurveInvariantsFp(t *testing.T) {
	tests := []struct {
		a, b, p int64
		disc    int64
		j       int64
		order   int64
		trace   int64
	}{
		{1, 1, 5, 4, 2, 9, -3},
		{0, 2, 7, 1, 0, 9, -1},
	}
	for _, tc := range tests {
		c, err := NewCurveFp(bi(tc.a), bi(tc.b), bi(tc.p))
		if err != nil {
			t.Fatal(err)
		}
		if c.Discriminant().Cmp(bi(tc.disc)) != 0 {
			t.Errorf("disc(%d,%d,%d)=%s want %d", tc.a, tc.b, tc.p, c.Discriminant(), tc.disc)
		}
		j, err := c.JInvariant()
		if err != nil {
			t.Fatal(err)
		}
		if j.Cmp(bi(tc.j)) != 0 {
			t.Errorf("j=%s want %d", j, tc.j)
		}
		if c.CurveOrderNaive().Cmp(bi(tc.order)) != 0 {
			t.Errorf("order=%s want %d", c.CurveOrderNaive(), tc.order)
		}
		if c.TraceOfFrobenius().Cmp(bi(tc.trace)) != 0 {
			t.Errorf("trace=%s want %d", c.TraceOfFrobenius(), tc.trace)
		}
	}
}

func TestSingularCurveRejected(t *testing.T) {
	// y^2 = x^3 (A=B=0) is singular; also 4A^3+27B^2=0.
	if _, err := NewCurveFp(bi(0), bi(0), bi(101)); err == nil {
		t.Errorf("expected singular rejection")
	}
	// 4*(-3)^3 + 27*2^2 = -108+108 = 0 => singular.
	if _, err := NewCurveFp(bi(-3), bi(2), bi(101)); err == nil {
		t.Errorf("expected singular rejection for cusp")
	}
}

func TestSupersingularAndTwist(t *testing.T) {
	c, _ := NewCurveFp(bi(0), bi(1), bi(5)) // y^2=x^3+1 over F5
	if !c.IsSupersingular() {
		t.Errorf("expected supersingular")
	}
	if c.CurveOrderNaive().Cmp(bi(6)) != 0 {
		t.Errorf("order want 6")
	}
	// Twist order: N + N_twist = 2(p+1).
	c2, _ := NewCurveFp(bi(2), bi(3), bi(101))
	sum := new(big.Int).Add(c2.CurveOrderNaive(), c2.TwistOrder())
	if sum.Cmp(bi(2*(101+1))) != 0 {
		t.Errorf("twist order relation failed: %s", sum)
	}
}

// ---------------------------------------------------------------------------
// Point counting
// ---------------------------------------------------------------------------

func TestPointCountingHasseAndBSGS(t *testing.T) {
	curves := []struct{ a, b, p int64 }{
		{1, 1, 5},
		{2, 3, 97},
		{37, 11, 1009},
		{0, 7, 1013},
	}
	rng := rand.New(rand.NewSource(99))
	for _, cc := range curves {
		c, err := NewCurveFp(bi(cc.a), bi(cc.b), bi(cc.p))
		if err != nil {
			t.Fatal(err)
		}
		naive := c.CurveOrderNaive()
		// Hasse bound.
		lo, hi := c.HasseInterval()
		if naive.Cmp(lo) < 0 || naive.Cmp(hi) > 0 {
			t.Errorf("order %s outside Hasse [%s,%s]", naive, lo, hi)
		}
		// BSGS agrees.
		bsgs, err := c.CurveOrderBSGS(rng)
		if err != nil {
			t.Fatalf("BSGS: %v", err)
		}
		if bsgs.Cmp(naive) != 0 {
			t.Errorf("BSGS %s != naive %s for p=%d", bsgs, naive, cc.p)
		}
		// Trace bound.
		ap := c.TraceOfFrobenius()
		twoRootP := new(big.Int).Mul(bigTwo, IntSqrt(bi(cc.p)))
		twoRootP.Add(twoRootP, bigOne)
		if new(big.Int).Abs(ap).Cmp(twoRootP) > 0 {
			t.Errorf("trace %s exceeds Hasse bound", ap)
		}
	}
}

func TestPointOrderAndDiscreteLog(t *testing.T) {
	c, _ := NewCurveFp(bi(37), bi(11), bi(1009))
	order := c.CurveOrderNaive()
	rng := rand.New(rand.NewSource(3))
	for i := 0; i < 15; i++ {
		p := c.RandomPointFp(rng)
		ob, err := c.PointOrderBSGS(p)
		if err != nil {
			t.Fatalf("PointOrderBSGS: %v", err)
		}
		on := c.PointOrder(p, order)
		if ob.Cmp(on) != 0 {
			t.Errorf("point order BSGS %s != naive %s", ob, on)
		}
		if !c.ScalarMul(ob, p).Infinity {
			t.Errorf("order*p != O")
		}
		if new(big.Int).Mod(order, ob).Sign() != 0 {
			t.Errorf("point order %s does not divide group order %s", ob, order)
		}
		// Discrete log: q = k*p, recover k mod ord.
		k := bi(int64(rng.Intn(1000)))
		q := c.ScalarMul(k, p)
		x, err := c.DiscreteLogBSGS(p, q, ob)
		if err != nil {
			t.Fatalf("DiscreteLogBSGS: %v", err)
		}
		if !c.PointEqual(c.ScalarMul(x, p), q) {
			t.Errorf("discrete log wrong")
		}
	}
}

func TestFindPointOfOrder(t *testing.T) {
	c, _ := NewCurveFp(bi(2), bi(3), bi(97))
	order := c.CurveOrderNaive()
	rng := rand.New(rand.NewSource(7))
	for _, d := range Divisors(order) {
		if d.Cmp(bigOne) == 0 {
			continue
		}
		pt, err := c.FindPointOfOrder(d, rng)
		if err != nil {
			continue // exponent may not admit this order
		}
		if c.PointOrder(pt, order).Cmp(d) != 0 {
			t.Errorf("FindPointOfOrder returned wrong order for %s", d)
		}
	}
}

// ---------------------------------------------------------------------------
// Division polynomials
// ---------------------------------------------------------------------------

func TestDivisionPolynomials(t *testing.T) {
	c, _ := NewCurveFp(bi(2), bi(3), bi(101))
	// psi_3 explicit value at x=1: 3 + 6A + 12B - A^2 = 3+12+36-4 = 47.
	if c.DivisionPolynomial3(bi(1), bi(0)).Cmp(bi(47)) != 0 {
		t.Errorf("psi_3(1) wrong: %s", c.DivisionPolynomial3(bi(1), bi(0)))
	}
	// psi_2 = 2y.
	if c.DivisionPolynomial2(bi(0), bi(9)).Cmp(bi(18)) != 0 {
		t.Errorf("psi_2 wrong")
	}
	// Torsion detection: psi_n(P)=0 iff n*P = O.
	rng := rand.New(rand.NewSource(4))
	order := c.CurveOrderNaive()
	for i := 0; i < 10; i++ {
		p := c.RandomPointFp(rng)
		if p.Y.Sign() == 0 {
			continue
		}
		ord := c.PointOrder(p, order)
		if ord.Cmp(bi(64)) > 0 {
			continue
		}
		n := int(ord.Int64())
		got, err := c.IsNTorsionByDivisionPoly(n, p)
		if err != nil {
			t.Fatalf("division poly: %v", err)
		}
		if !got {
			t.Errorf("psi_%d(P) should vanish for point of order %d", n, n)
		}
		if c.IsNTorsionPoint(bi(int64(n)), p) != got {
			t.Errorf("division-poly torsion disagrees with direct test")
		}
		// A non-multiple must not vanish.
		if n > 2 {
			ng, _ := c.IsNTorsionByDivisionPoly(n-1, p)
			if ng {
				t.Errorf("psi_%d(P) unexpectedly vanished", n-1)
			}
		}
	}
}

// ---------------------------------------------------------------------------
// Weil pairing
// ---------------------------------------------------------------------------

func TestWeilPairing(t *testing.T) {
	c, _ := NewCurveFp(bi(0), bi(2), bi(7)) // full 3-torsion, order 9
	n := bi(3)
	p := PointFp{X: bi(0), Y: bi(4)}
	q := PointFp{X: bi(3), Y: bi(1)}
	if !c.IsOnCurve(p) || !c.IsOnCurve(q) {
		t.Fatal("test points not on curve")
	}
	rng := rand.New(rand.NewSource(1))
	e, err := c.WeilPairing(n, p, q, rng)
	if err != nil {
		t.Fatalf("WeilPairing: %v", err)
	}
	// Known primitive cube root of unity value.
	if e.Cmp(bi(2)) != 0 {
		t.Errorf("e_3(P,Q)=%s want 2", e)
	}
	// n-th root of unity.
	if !c.WeilPairingRootOfUnity(e, n) {
		t.Errorf("pairing value not an n-th root of unity")
	}
	// Nondegeneracy: primitive root.
	if MultiplicativeOrder(e, c.P).Cmp(n) != 0 {
		t.Errorf("pairing not nondegenerate")
	}
	// Bilinearity: e(2P,Q) = e(P,Q)^2.
	e2, _ := c.WeilPairing(n, c.ScalarMul(bi(2), p), q, rand.New(rand.NewSource(2)))
	if e2.Cmp(ModExp2(e, 2, c.P)) != 0 {
		t.Errorf("bilinearity failed: e(2P,Q)=%s e^2=%s", e2, ModExp2(e, 2, c.P))
	}
	// Antisymmetry: e(P,Q)*e(Q,P) = 1.
	eqp, _ := c.WeilPairing(n, q, p, rand.New(rand.NewSource(3)))
	if ModMul(e, eqp, c.P).Cmp(bigOne) != 0 {
		t.Errorf("antisymmetry failed")
	}
	// e(P,P) = 1.
	epp, _ := c.WeilPairing(n, p, p, rand.New(rand.NewSource(4)))
	if epp.Cmp(bigOne) != 0 {
		t.Errorf("e(P,P) != 1")
	}
}

func ModExp2(a *big.Int, e int64, p *big.Int) *big.Int {
	return new(big.Int).Exp(a, bi(e), p)
}

func TestTorsionBasisFp(t *testing.T) {
	c, _ := NewCurveFp(bi(0), bi(2), bi(7))
	rng := rand.New(rand.NewSource(8))
	p, q, ok := c.TorsionBasisFp(bi(3), rng)
	if !ok {
		t.Fatal("expected full 3-torsion basis")
	}
	// Both order 3.
	if !c.ScalarMul(bi(3), p).Infinity || !c.ScalarMul(bi(3), q).Infinity {
		t.Errorf("basis points not 3-torsion")
	}
	v, err := c.WeilPairing(bi(3), p, q, rng)
	if err != nil {
		t.Fatalf("pairing: %v", err)
	}
	if MultiplicativeOrder(v, c.P).Cmp(bi(3)) != 0 {
		t.Errorf("basis not independent")
	}
}

// ---------------------------------------------------------------------------
// Group law and torsion over Q
// ---------------------------------------------------------------------------

func TestGroupLawQ(t *testing.T) {
	c, _ := NewCurveQInt(bi(0), bi(1)) // y^2=x^3+1
	pts, err := c.TorsionSubgroup()
	if err != nil {
		t.Fatal(err)
	}
	inf := PointAtInfinityQ()
	for _, p := range pts {
		if !c.IsOnCurve(p) {
			t.Fatalf("torsion point %s not on curve", p)
		}
		if !c.PointEqual(c.Add(p, inf), p) {
			t.Errorf("identity failed over Q")
		}
		if !c.Add(p, c.Neg(p)).Infinity {
			t.Errorf("inverse failed over Q")
		}
	}
	// Associativity within the finite torsion subgroup.
	for _, a := range pts {
		for _, b := range pts {
			for _, d := range pts {
				l := c.Add(c.Add(a, b), d)
				r := c.Add(a, c.Add(b, d))
				if !c.PointEqual(l, r) {
					t.Fatalf("associativity failed over Q")
				}
			}
		}
	}
}

func TestTorsionSubgroupQ(t *testing.T) {
	tests := []struct {
		a, b  int64
		order int
		expo  int
	}{
		{-1, 0, 4, 2},    // y^2=x^3-x : Z/2 x Z/2
		{0, 1, 6, 6},     // y^2=x^3+1 : Z/6
		{-43, 166, 7, 7}, // 7a curve: Z/7
		{0, -2, 1, 1},    // y^2=x^3-2 : trivial torsion, rank 1
	}
	for _, tc := range tests {
		c, err := NewCurveQInt(bi(tc.a), bi(tc.b))
		if err != nil {
			t.Fatalf("curve (%d,%d): %v", tc.a, tc.b, err)
		}
		ord, err := c.TorsionOrder()
		if err != nil {
			t.Fatal(err)
		}
		if ord != tc.order {
			pts, _ := c.TorsionSubgroup()
			t.Errorf("torsion order (%d,%d)=%d want %d: %v", tc.a, tc.b, ord, tc.order, pts)
		}
		expo, _ := c.TorsionExponent()
		if expo != tc.expo {
			t.Errorf("torsion exponent (%d,%d)=%d want %d", tc.a, tc.b, expo, tc.expo)
		}
		// Every torsion point genuinely has finite order.
		pts, _ := c.TorsionSubgroup()
		for _, p := range pts {
			if c.PointOrderQ(p) == 0 {
				t.Errorf("point %s reported infinite order in torsion set", p)
			}
		}
	}
}

func TestNagellLutzIntegrality(t *testing.T) {
	c, _ := NewCurveQInt(bi(0), bi(1))
	pts, err := c.TorsionSubgroup()
	if err != nil {
		t.Fatal(err)
	}
	for _, p := range pts {
		if !c.IsIntegralPoint(p) {
			t.Errorf("torsion point %s is not integral (Nagell-Lutz)", p)
		}
	}
}

// ---------------------------------------------------------------------------
// Heights and rank heuristics
// ---------------------------------------------------------------------------

func TestHeightsAndRank(t *testing.T) {
	// y^2 = x^3 - 2 has rank 1 generated by (3,5).
	c, _ := NewCurveQInt(bi(0), bi(-2))
	p, err := c.NewPointQ(rat(3, 1), rat(5, 1))
	if err != nil {
		t.Fatalf("point: %v", err)
	}
	h := c.CanonicalHeightApprox(p, 6)
	if h <= 0 {
		t.Errorf("canonical height should be positive, got %g", h)
	}
	// Convergence: successive estimates stabilize.
	h1 := c.CanonicalHeightApprox(p, 5)
	h2 := c.CanonicalHeightApprox(p, 7)
	if math.Abs(h1-h2) > 0.05 {
		t.Errorf("canonical height not converging: %g vs %g", h1, h2)
	}
	// A single non-torsion point gives rank lower bound 1.
	if rb := c.RankLowerBound([]PointQ{p}, 6, 1e-6); rb != 1 {
		t.Errorf("rank lower bound = %d want 1", rb)
	}
	// Point and its double are dependent: rank still 1.
	p2 := c.Double(p)
	if rb := c.RankLowerBound([]PointQ{p, p2}, 6, 1e-2); rb != 1 {
		t.Errorf("dependent points rank = %d want 1", rb)
	}
	if !c.AreIndependent([]PointQ{p}, 6, 1e-6) {
		t.Errorf("single point should be independent")
	}
}

func TestNaiveHeight(t *testing.T) {
	c, _ := NewCurveQInt(bi(0), bi(1))
	p, _ := c.NewPointQ(rat(2, 1), rat(3, 1))
	// x=2 -> height log 2.
	if math.Abs(c.NaiveHeightQ(p)-math.Log(2)) > 1e-12 {
		t.Errorf("naive height wrong: %g", c.NaiveHeightQ(p))
	}
	if c.NaiveHeightQ(PointAtInfinityQ()) != 0 {
		t.Errorf("height of infinity should be 0")
	}
}

// ---------------------------------------------------------------------------
// Isomorphisms and twists
// ---------------------------------------------------------------------------

func TestIsomorphismFp(t *testing.T) {
	c, _ := NewCurveFp(bi(2), bi(3), bi(101))
	u := bi(5)
	scaled, err := c.ScaleCurveFp(u)
	if err != nil {
		t.Fatal(err)
	}
	if !c.IsIsomorphicFp(scaled) {
		t.Errorf("scaled curve should be isomorphic")
	}
	if !c.SameJInvariantFp(scaled) {
		t.Errorf("isomorphic curves share j-invariant")
	}
	// The isomorphism maps points on-curve to points on-curve.
	uu, ok := c.IsomorphismScaleFp(scaled)
	if !ok {
		t.Fatal("expected isomorphism scale")
	}
	rng := rand.New(rand.NewSource(6))
	for i := 0; i < 20; i++ {
		p := c.RandomPointFp(rng)
		img := c.TransformPointFp(uu, p)
		if !scaled.IsOnCurve(img) {
			t.Errorf("transformed point not on scaled curve")
		}
	}
	// Quadratic twist by a non-residue is not isomorphic but shares j.
	if IsQuadraticResidue(bi(2), bi(101)) {
		t.Fatal("expected 2 to be a non-residue mod 101")
	}
	tw, _ := c.QuadraticTwistFp(bi(2))
	if c.IsIsomorphicFp(tw) {
		t.Errorf("quadratic twist should not be isomorphic")
	}
	if !c.IsQuadraticTwistFp(tw) {
		t.Errorf("expected quadratic twist relation")
	}
}

func TestIsomorphismQ(t *testing.T) {
	c, _ := NewCurveQInt(bi(1), bi(1))
	// Scale by u=2: A'=16, B'=64.
	c2, _ := NewCurveQ(rat(16, 1), rat(64, 1))
	if !c.IsIsomorphicQ(c2) {
		t.Errorf("curves should be Q-isomorphic")
	}
	u, ok := c.IsomorphismScaleQ(c2)
	if !ok || u.Cmp(rat(2, 1)) != 0 {
		t.Errorf("isomorphism scale wrong: %v %v", u, ok)
	}
	// Quadratic twist by 3 (non-square): shares j but not isomorphic.
	tw, _ := c.QuadraticTwistQ(rat(3, 1))
	if !c.SameJInvariantQ(tw) {
		t.Errorf("twist should share j-invariant")
	}
	if c.IsIsomorphicQ(tw) {
		t.Errorf("quadratic twist should not be Q-isomorphic")
	}
	if !c.IsQuadraticTwistQ(tw) {
		t.Errorf("expected quadratic twist over Q")
	}
}

// ---------------------------------------------------------------------------
// Reduction modulo p
// ---------------------------------------------------------------------------

func TestReductionModP(t *testing.T) {
	c, _ := NewCurveQInt(bi(0), bi(-2)) // disc = -1728 = -2^6 * 27
	bad, err := c.BadPrimes()
	if err != nil {
		t.Fatal(err)
	}
	badSet := map[string]bool{}
	for _, p := range bad {
		badSet[p.String()] = true
	}
	if !badSet["2"] || !badSet["3"] {
		t.Errorf("bad primes should include 2 and 3, got %v", bad)
	}
	if c.HasGoodReduction(bi(2)) || c.HasGoodReduction(bi(3)) {
		t.Errorf("2 and 3 are bad primes")
	}
	if !c.HasGoodReduction(bi(7)) {
		t.Errorf("7 should be good reduction")
	}
	// a_7 = 7 + 1 - #E~(F_7); the curve is anomalous at 7 with #E=7.
	red, err := c.ReduceModP(bi(7))
	if err != nil {
		t.Fatal(err)
	}
	if red.CurveOrderNaive().Cmp(bi(7)) != 0 {
		t.Errorf("reduced order want 7, got %s", red.CurveOrderNaive())
	}
	ap, _ := c.APAtPrime(bi(7))
	if ap.Cmp(bi(1)) != 0 {
		t.Errorf("a_7 want 1, got %s", ap)
	}
	// A rational point reduces onto the reduced curve.
	p, _ := c.NewPointQ(rat(3, 1), rat(5, 1))
	rp, err := c.ReducePointModP(bi(7), p)
	if err != nil {
		t.Fatal(err)
	}
	if !red.IsOnCurve(rp) {
		t.Errorf("reduced point not on reduced curve")
	}
	// Bad reduction is reported.
	if _, err := c.ReduceModP(bi(2)); err == nil {
		t.Errorf("expected bad reduction at 2")
	}
}

// ---------------------------------------------------------------------------
// Multiplicative order and primitive roots
// ---------------------------------------------------------------------------

func TestMultiplicativeOrder(t *testing.T) {
	p := bi(101)
	g := PrimitiveRoot(p)
	if g == nil {
		t.Fatal("no primitive root found")
	}
	if !IsPrimitiveRoot(g, p) {
		t.Errorf("PrimitiveRoot did not return a primitive root")
	}
	if MultiplicativeOrder(g, p).Cmp(bi(100)) != 0 {
		t.Errorf("primitive root order should be 100")
	}
	// Order of g^2 is 50.
	g2 := ModSquare(g, p)
	if MultiplicativeOrder(g2, p).Cmp(bi(50)) != 0 {
		t.Errorf("order of g^2 should be 50")
	}
}

func iMin(a, b int) int {
	if a < b {
		return a
	}
	return b
}

// ---------------------------------------------------------------------------
// Runnable examples
// ---------------------------------------------------------------------------

func ExampleCurveFp_Order() {
	// The curve y^2 = x^3 + x + 1 over F_5 has nine points.
	c, _ := NewCurveFp(big.NewInt(1), big.NewInt(1), big.NewInt(5))
	fmt.Println(c.Order())
	// Output: 9
}

func ExampleCurveFp_WeilPairing() {
	// Full 3-torsion curve over F_7; the Weil pairing of a basis is a
	// primitive cube root of unity.
	c, _ := NewCurveFp(big.NewInt(0), big.NewInt(2), big.NewInt(7))
	p := PointFp{X: big.NewInt(0), Y: big.NewInt(4)}
	q := PointFp{X: big.NewInt(3), Y: big.NewInt(1)}
	e, _ := c.WeilPairing(big.NewInt(3), p, q, rand.New(rand.NewSource(1)))
	fmt.Println(e)
	// Output: 2
}

func ExampleCurveQ_TorsionSubgroup() {
	// y^2 = x^3 + 1 has rational torsion isomorphic to Z/6Z.
	c, _ := NewCurveQInt(big.NewInt(0), big.NewInt(1))
	order, _ := c.TorsionOrder()
	fmt.Println(order)
	// Output: 6
}
