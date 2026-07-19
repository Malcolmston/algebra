package contfrac

import (
	"math"
	"math/big"
	"reflect"
	"testing"
)

const tol = 1e-9

func eqI64(a, b []int64) bool { return reflect.DeepEqual(a, b) }

func TestFromRational(t *testing.T) {
	tests := []struct {
		p, q int64
		want []int64
	}{
		{415, 93, []int64{4, 2, 6, 7}},
		{3, 2, []int64{1, 2}},
		{1, 1, []int64{1}},
		{0, 1, []int64{0}},
		{7, 1, []int64{7}},
		{1, 3, []int64{0, 3}},
		{-415, 93, []int64{-5, 1, 1, 6, 7}}, // floor convention: a0 = -5
		{10, 4, []int64{2, 2}},              // 5/2
	}
	for _, tc := range tests {
		got := FromRational(tc.p, tc.q)
		if !eqI64(got, tc.want) {
			t.Errorf("FromRational(%d,%d) = %v, want %v", tc.p, tc.q, got, tc.want)
		}
		// Round trip: evaluating the CF must recover the reduced fraction.
		f := got.Frac()
		wp, wq := ReduceFraction(tc.p, tc.q)
		if f.Num != wp || f.Den != wq {
			t.Errorf("round trip %d/%d -> %v", tc.p, tc.q, f)
		}
	}
}

func TestFromRationalBig(t *testing.T) {
	r := big.NewRat(415, 93)
	got := FromRationalBig(r)
	if !eqI64(got, []int64{4, 2, 6, 7}) {
		t.Fatalf("FromRationalBig = %v", got)
	}
	if got.Rat().Cmp(r) != 0 {
		t.Errorf("big round trip failed: %v", got.Rat())
	}
}

func TestConvergents(t *testing.T) {
	// Convergents of [3;7,15,1,292] (pi) are the classical approximations.
	cf := CF{3, 7, 15, 1, 292}
	conv := cf.Convergents()
	want := []Frac{{3, 1}, {22, 7}, {333, 106}, {355, 113}, {103993, 33102}}
	if len(conv) != len(want) {
		t.Fatalf("len convergents = %d", len(conv))
	}
	for i := range want {
		if conv[i] != want[i] {
			t.Errorf("convergent %d = %v, want %v", i, conv[i], want[i])
		}
	}
	// Convergent(k) agrees.
	if cf.Convergent(3) != (Frac{355, 113}) {
		t.Errorf("Convergent(3) = %v", cf.Convergent(3))
	}
	// Numerators and denominators.
	if !eqI64(cf.Numerators(), []int64{3, 22, 333, 355, 103993}) {
		t.Errorf("numerators = %v", cf.Numerators())
	}
	if !eqI64(cf.Denominators(), []int64{1, 7, 106, 113, 33102}) {
		t.Errorf("denominators = %v", cf.Denominators())
	}
}

func TestSemiconvergents(t *testing.T) {
	// [1;2] traces the Stern-Brocot path 1/1, 2/1, 3/2.
	got := CF{1, 2}.Semiconvergents()
	want := []Frac{{1, 1}, {2, 1}, {3, 2}}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("Semiconvergents = %v, want %v", got, want)
	}
}

func TestValueAndFrac(t *testing.T) {
	cf := FromRational(355, 113)
	if math.Abs(cf.Value()-355.0/113.0) > tol {
		t.Errorf("Value = %v", cf.Value())
	}
	if cf.Frac() != (Frac{355, 113}) {
		t.Errorf("Frac = %v", cf.Frac())
	}
	if cf.Rat().Cmp(big.NewRat(355, 113)) != 0 {
		t.Errorf("Rat = %v", cf.Rat())
	}
}

func TestFromFloat(t *testing.T) {
	cf := FromFloat(math.Sqrt2, 12)
	// sqrt(2) = [1;2,2,2,...]
	if cf[0] != 1 {
		t.Errorf("a0 = %d", cf[0])
	}
	for i := 1; i < 8; i++ {
		if cf[i] != 2 {
			t.Errorf("term %d = %d, want 2", i, cf[i])
		}
	}
	// A 12-term truncation is the convergent 19601/13860, good to ~2e-9.
	if math.Abs(cf.Value()-math.Sqrt2) > 1e-7 {
		t.Errorf("value = %v", cf.Value())
	}
}

func TestCanonicalAndEqual(t *testing.T) {
	a := CF{1, 2, 3, 1} // equals [1;2,4]
	b := CF{1, 2, 4}
	if !a.EqualValue(b) {
		t.Errorf("EqualValue failed")
	}
	if !a.Canonical().Equal(b) {
		t.Errorf("Canonical = %v, want %v", a.Canonical(), b)
	}
	if a.Equal(b) {
		t.Errorf("raw Equal should be false")
	}
}

func TestReverseContinuant(t *testing.T) {
	cf := CF{2, 3, 4}
	// Continuant is invariant under reversal.
	if Continuant(cf...) != Continuant(cf.Reverse()...) {
		t.Errorf("continuant not reversal-invariant")
	}
	if Continuant(2, 3, 4) != 30 {
		t.Errorf("Continuant(2,3,4) = %d, want 30", Continuant(2, 3, 4))
	}
	if ContinuantBig([]int64{2, 3, 4}).Int64() != 30 {
		t.Errorf("ContinuantBig mismatch")
	}
}

func TestParseCF(t *testing.T) {
	cf, err := ParseCF("[3; 7, 15, 1]")
	if err != nil {
		t.Fatal(err)
	}
	if !cf.Equal(CF{3, 7, 15, 1}) {
		t.Errorf("ParseCF = %v", cf)
	}
	if (CF{3, 7, 15, 1}).String() != "[3; 7, 15, 1]" {
		t.Errorf("String = %q", CF{3, 7, 15, 1}.String())
	}
}

func TestBestApproximation(t *testing.T) {
	tests := []struct {
		x            float64
		maxDen       int64
		wantP, wantQ int64
	}{
		{math.Pi, 7, 22, 7},
		{math.Pi, 100, 311, 99},
		{math.Pi, 113, 355, 113},
		{math.Pi, 106, 333, 106},
		{math.E, 100, 193, 71},
		{0.5, 10, 1, 2},
		{-math.Pi, 7, -22, 7},
	}
	for _, tc := range tests {
		p, q := BestApproximation(tc.x, tc.maxDen)
		if p != tc.wantP || q != tc.wantQ {
			t.Errorf("BestApproximation(%v,%d) = %d/%d, want %d/%d", tc.x, tc.maxDen, p, q, tc.wantP, tc.wantQ)
		}
	}
}

func TestBestApproximationRat(t *testing.T) {
	got := BestApproximationRat(big.NewRat(103993, 33102), 200)
	if got.Cmp(big.NewRat(355, 113)) != 0 {
		t.Errorf("BestApproximationRat = %v, want 355/113", got)
	}
}

func TestRationalizeAndSmallestDen(t *testing.T) {
	if got := Rationalize(0.333333333, 10); got != (Frac{1, 3}) {
		t.Errorf("Rationalize = %v", got)
	}
	// 22/7 is off by ~1.26e-3, so it is the smallest-denominator fraction
	// within 2e-3 of pi (denominators 1..6 do no better than ~2.5e-2).
	if got := SmallestDenominatorWithin(math.Pi, 2e-3); got != (Frac{22, 7}) {
		t.Errorf("SmallestDenominatorWithin = %v, want 22/7", got)
	}
}

func TestSqrtCF(t *testing.T) {
	tests := []struct {
		n      int64
		a0     int64
		period []int64
	}{
		{2, 1, []int64{2}},
		{3, 1, []int64{1, 2}},
		{7, 2, []int64{1, 1, 1, 4}},
		{13, 3, []int64{1, 1, 1, 1, 6}},
		{23, 4, []int64{1, 3, 1, 8}},
		{4, 2, nil}, // perfect square
	}
	for _, tc := range tests {
		sc := SqrtCF(tc.n)
		if sc.Head[0] != tc.a0 {
			t.Errorf("SqrtCF(%d) a0 = %d, want %d", tc.n, sc.Head[0], tc.a0)
		}
		if !eqI64(sc.Period, tc.period) {
			t.Errorf("SqrtCF(%d) period = %v, want %v", tc.n, sc.Period, tc.period)
		}
		if tc.period != nil {
			if !IsPeriodPalindrome(tc.period) {
				t.Errorf("period of sqrt(%d) not palindromic", tc.n)
			}
			if math.Abs(sc.Value()-math.Sqrt(float64(tc.n))) > 1e-9 {
				t.Errorf("SqrtCF(%d).Value = %v", tc.n, sc.Value())
			}
		}
	}
}

func TestSqrtConvergents(t *testing.T) {
	// Convergents of sqrt(2): 1, 3/2, 7/5, 17/12, 41/29
	conv := SqrtConvergents(2, 5)
	want := []Frac{{1, 1}, {3, 2}, {7, 5}, {17, 12}, {41, 29}}
	if !reflect.DeepEqual(conv, want) {
		t.Errorf("SqrtConvergents(2,5) = %v, want %v", conv, want)
	}
	if SqrtConvergent(2, 3) != (Frac{17, 12}) {
		t.Errorf("SqrtConvergent(2,3) = %v", SqrtConvergent(2, 3))
	}
	if SqrtCFPeriodLength(7) != 4 {
		t.Errorf("period length sqrt(7) = %d", SqrtCFPeriodLength(7))
	}
}

func TestQuadraticSurdCF(t *testing.T) {
	// Golden ratio = [1;(1)]
	g := NewSurd(1, 2, 5).CF()
	if !g.IsPurelyPeriodic() || !eqI64(g.Period, []int64{1}) {
		t.Errorf("golden ratio CF = %v", g)
	}
	// (1+sqrt(3)) = [2;1,2,...]? check round trip via Surd value.
	s := NewSurd(1, 1, 3)
	pc := s.CF()
	if math.Abs(pc.Value()-s.Value()) > 1e-9 {
		t.Errorf("surd CF value mismatch: %v vs %v", pc.Value(), s.Value())
	}
	// Round-trip: PeriodicCF.Surd of golden ratio recovers (1+sqrt5)/2.
	back := g.Surd()
	if math.Abs(back.Value()-GoldenRatio()) > 1e-12 {
		t.Errorf("Surd round trip = %v", back)
	}
	// Surd of sqrt(7) periodic CF recovers sqrt(7).
	if math.Abs(SqrtCF(7).Surd().Value()-math.Sqrt(7)) > 1e-9 {
		t.Errorf("sqrt7 surd = %v", SqrtCF(7).Surd().Value())
	}
}

func TestSurdFloor(t *testing.T) {
	// floor((1+sqrt(2))/1) = 2
	if got := NewSurd(1, 1, 2).Floor(); got != 2 {
		t.Errorf("floor = %d", got)
	}
	// floor(sqrt(2)) = 1
	if got := NewSurd(0, 1, 2).Floor(); got != 1 {
		t.Errorf("floor sqrt2 = %d", got)
	}
	// Conjugate of (1+sqrt(2)) is 1-sqrt(2) ~ -0.414, with floor -1.
	conj := NewSurd(1, 1, 2).Conjugate()
	if math.Abs(conj.Value()-(1-math.Sqrt2)) > 1e-9 {
		t.Errorf("conjugate value = %v, want %v", conj.Value(), 1-math.Sqrt2)
	}
	if conj.Floor() != -1 {
		t.Errorf("conjugate floor = %d, want -1", conj.Floor())
	}
}

func TestPell(t *testing.T) {
	tests := []struct {
		D           int64
		wantX       int64
		wantY       int64
		negSolvable bool
	}{
		{2, 3, 2, true},
		{3, 2, 1, false},
		{5, 9, 4, true},
		{7, 8, 3, false},
		{13, 649, 180, true},
		{61, 1766319049, 226153980, true},
	}
	for _, tc := range tests {
		x, y, ok := PellFundamentalInt(tc.D)
		if !ok {
			t.Errorf("Pell(%d) not ok", tc.D)
			continue
		}
		if x != tc.wantX || y != tc.wantY {
			t.Errorf("Pell(%d) = (%d,%d), want (%d,%d)", tc.D, x, y, tc.wantX, tc.wantY)
		}
		if !IsPellSolution(tc.D, x, y, 1) {
			t.Errorf("Pell(%d) solution invalid", tc.D)
		}
		if PellNegativeSolvable(tc.D) != tc.negSolvable {
			t.Errorf("PellNegativeSolvable(%d) = %v, want %v", tc.D, PellNegativeSolvable(tc.D), tc.negSolvable)
		}
	}
	// No solution for perfect square.
	if _, _, ok := PellFundamental(9); ok {
		t.Errorf("Pell(9) should have no solution")
	}
}

func TestPellSolutions(t *testing.T) {
	sols := PellSolutions(2, 3)
	// (3,2),(17,12),(99,70)
	want := [][2]int64{{3, 2}, {17, 12}, {99, 70}}
	for i, s := range sols {
		if s[0].Int64() != want[i][0] || s[1].Int64() != want[i][1] {
			t.Errorf("solution %d = (%v,%v), want %v", i, s[0], s[1], want[i])
		}
		if s[0].Int64()*s[0].Int64()-2*s[1].Int64()*s[1].Int64() != 1 {
			t.Errorf("solution %d not valid", i)
		}
	}
	nx, ny, ok := PellNegative(2)
	if !ok || nx.Int64() != 1 || ny.Int64() != 1 {
		t.Errorf("PellNegative(2) = (%v,%v,%v)", nx, ny, ok)
	}
	if nx.Int64()*nx.Int64()-2*ny.Int64()*ny.Int64() != -1 {
		t.Errorf("negative Pell invalid")
	}
}

func TestMediantAndSternBrocot(t *testing.T) {
	mp, mq := Mediant(1, 2, 1, 3)
	if mp != 2 || mq != 5 {
		t.Errorf("Mediant = %d/%d, want 2/5", mp, mq)
	}
	tests := []struct {
		p, q int64
		path string
	}{
		{1, 1, ""},
		{2, 1, "R"},
		{1, 2, "L"},
		{3, 2, "RL"},
		{5, 3, "RLR"},
		{2, 5, "LLR"},
	}
	for _, tc := range tests {
		if got := SternBrocotPath(tc.p, tc.q); got != tc.path {
			t.Errorf("SternBrocotPath(%d,%d) = %q, want %q", tc.p, tc.q, got, tc.path)
		}
		if got := SternBrocotFromPath(tc.path); got != (Frac{tc.p, tc.q}) {
			t.Errorf("SternBrocotFromPath(%q) = %v, want %d/%d", tc.path, got, tc.p, tc.q)
		}
	}
	// Parent and children.
	if SternBrocotParent(3, 2) != (Frac{2, 1}) {
		t.Errorf("parent of 3/2 = %v", SternBrocotParent(3, 2))
	}
	l, r := SternBrocotChildren(1, 1)
	if l != (Frac{1, 2}) || r != (Frac{2, 1}) {
		t.Errorf("children of 1/1 = %v, %v", l, r)
	}
	if SternBrocotDepth(5, 3) != 3 {
		t.Errorf("depth 5/3 = %d", SternBrocotDepth(5, 3))
	}
}

func TestPathCFRoundTrip(t *testing.T) {
	tests := []struct {
		path string
		cf   []int64
	}{
		{"", []int64{1}},
		{"R", []int64{2}},
		{"L", []int64{0, 2}},
		{"RL", []int64{1, 2}},
		{"RLR", []int64{1, 1, 2}},
	}
	for _, tc := range tests {
		if got := PathToCF(tc.path); !eqI64(got, tc.cf) {
			t.Errorf("PathToCF(%q) = %v, want %v", tc.path, got, tc.cf)
		}
		if got := CFToPath(CF(tc.cf)); got != tc.path {
			t.Errorf("CFToPath(%v) = %q, want %q", tc.cf, got, tc.path)
		}
	}
	// CFToPath is insensitive to the trailing-1 ambiguity.
	if CFToPath(CF{1, 1, 1}) != CFToPath(CF{1, 2}) {
		t.Errorf("CFToPath ambiguity not handled")
	}
}

func TestFarey(t *testing.T) {
	f5 := FareySequence(5)
	want := []Frac{{0, 1}, {1, 5}, {1, 4}, {1, 3}, {2, 5}, {1, 2}, {3, 5}, {2, 3}, {3, 4}, {4, 5}, {1, 1}}
	if !reflect.DeepEqual(f5, want) {
		t.Errorf("FareySequence(5) = %v", f5)
	}
	if FareyLength(5) != 11 {
		t.Errorf("FareyLength(5) = %d", FareyLength(5))
	}
	if FareyLength(7) != int64(len(FareySequence(7))) {
		t.Errorf("FareyLength inconsistent with FareySequence")
	}
	// Successor / predecessor.
	if FareySuccessor(1, 3, 7) != (Frac{2, 5}) {
		t.Errorf("FareySuccessor(1/3,7) = %v", FareySuccessor(1, 3, 7))
	}
	if FareyPredecessor(1, 2, 5) != (Frac{2, 5}) {
		t.Errorf("FareyPredecessor(1/2,5) = %v", FareyPredecessor(1, 2, 5))
	}
	if FareySuccessor(0, 1, 5) != (Frac{1, 5}) {
		t.Errorf("FareySuccessor(0/1,5) = %v", FareySuccessor(0, 1, 5))
	}
	// Neighbours consistency against the full sequence.
	for i := 1; i < len(f5)-1; i++ {
		prev, next := FareyNeighbors(f5[i].Num, f5[i].Den, 5)
		if prev != f5[i-1] || next != f5[i+1] {
			t.Errorf("neighbours of %v = (%v,%v), want (%v,%v)", f5[i], prev, next, f5[i-1], f5[i+1])
		}
	}
	if !AreFareyNeighbors(1, 3, 2, 5) {
		t.Errorf("1/3 and 2/5 should be Farey neighbours")
	}
	if AreFareyNeighbors(1, 3, 2, 3) {
		t.Errorf("1/3 and 2/3 are not adjacent (bc-ad = 3, not 1)")
	}
}

func TestEgyptian(t *testing.T) {
	tests := []struct {
		p, q int64
		want []int64
	}{
		{5, 6, []int64{2, 3}},
		{4, 13, []int64{4, 18, 468}},
		{3, 7, []int64{3, 11, 231}},
		{1, 2, []int64{2}},
	}
	for _, tc := range tests {
		got, err := EgyptianFraction(tc.p, tc.q)
		if err != nil {
			t.Fatalf("EgyptianFraction(%d,%d): %v", tc.p, tc.q, err)
		}
		if !eqI64(got, tc.want) {
			t.Errorf("EgyptianFraction(%d,%d) = %v, want %v", tc.p, tc.q, got, tc.want)
		}
		// Verify the sum is exact.
		if SumUnitFractions(got).Cmp(big.NewRat(tc.p, tc.q)) != 0 {
			t.Errorf("sum of %v != %d/%d", got, tc.p, tc.q)
		}
	}
	// Improper fraction is an error.
	if _, err := EgyptianFraction(3, 2); err == nil {
		t.Errorf("expected error for improper fraction")
	}
	// Big version matches.
	dens, err := EgyptianFractionBig(big.NewRat(5, 121))
	if err != nil {
		t.Fatal(err)
	}
	sum := new(big.Rat)
	for _, d := range dens {
		sum.Add(sum, new(big.Rat).SetFrac(big.NewInt(1), d))
	}
	if sum.Cmp(big.NewRat(5, 121)) != 0 {
		t.Errorf("EgyptianFractionBig sum = %v", sum)
	}
	if !IsUnitFraction(3, 3) || IsUnitFraction(2, 3) {
		t.Errorf("IsUnitFraction wrong")
	}
}

func TestConstants(t *testing.T) {
	// e continued fraction.
	if !eqI64([]int64(ECF(11)), []int64{2, 1, 2, 1, 1, 4, 1, 1, 6, 1, 1}) {
		t.Errorf("ECF(11) = %v", ECF(11))
	}
	if math.Abs(ECF(20).Value()-math.E) > 1e-9 {
		t.Errorf("ECF value = %v", ECF(20).Value())
	}
	// pi continued fraction.
	if !eqI64([]int64(PiCF(13)), []int64{3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14}) {
		t.Errorf("PiCF(13) = %v", PiCF(13))
	}
	if math.Abs(PiCF(15).Value()-math.Pi) > 1e-12 {
		t.Errorf("PiCF value = %v", PiCF(15).Value())
	}
	// pi convergents include 22/7 and 355/113.
	conv := PiConvergents(4)
	if conv[1] != (Frac{22, 7}) || conv[3] != (Frac{355, 113}) {
		t.Errorf("PiConvergents = %v", conv)
	}
	// Golden ratio.
	if !eqI64([]int64(GoldenRatioCF(5)), []int64{1, 1, 1, 1, 1}) {
		t.Errorf("GoldenRatioCF = %v", GoldenRatioCF(5))
	}
	if math.Abs(GoldenRatio()-1.6180339887) > 1e-9 {
		t.Errorf("GoldenRatio = %v", GoldenRatio())
	}
	if math.Abs(SqrtTwoCF(15).Value()-math.Sqrt2) > 1e-9 {
		t.Errorf("SqrtTwoCF value = %v", SqrtTwoCF(15).Value())
	}
	// PiFloat close to math.Pi.
	pf, _ := PiFloat(20).Float64()
	if math.Abs(pf-math.Pi) > 1e-12 {
		t.Errorf("PiFloat = %v", pf)
	}
}

func TestFracArithmetic(t *testing.T) {
	a := NewFrac(1, 2)
	b := NewFrac(1, 3)
	if a.Add(b) != (Frac{5, 6}) {
		t.Errorf("Add = %v", a.Add(b))
	}
	if a.Sub(b) != (Frac{1, 6}) {
		t.Errorf("Sub = %v", a.Sub(b))
	}
	if a.Mul(b) != (Frac{1, 6}) {
		t.Errorf("Mul = %v", a.Mul(b))
	}
	if a.Div(b) != (Frac{3, 2}) {
		t.Errorf("Div = %v", a.Div(b))
	}
	if a.Inv() != (Frac{2, 1}) {
		t.Errorf("Inv = %v", a.Inv())
	}
	if a.Neg() != (Frac{-1, 2}) {
		t.Errorf("Neg = %v", a.Neg())
	}
	if NewFrac(7, 2).Floor() != 3 || NewFrac(7, 2).Ceil() != 4 || NewFrac(7, 2).Round() != 4 {
		t.Errorf("Floor/Ceil/Round wrong")
	}
	if NewFrac(-7, 2).Floor() != -4 || NewFrac(-7, 2).Round() != -4 {
		t.Errorf("negative Floor/Round wrong")
	}
	if NewFrac(2, 3).Pow(3) != (Frac{8, 27}) {
		t.Errorf("Pow = %v", NewFrac(2, 3).Pow(3))
	}
	if NewFrac(2, 3).Pow(-1) != (Frac{3, 2}) {
		t.Errorf("Pow(-1) = %v", NewFrac(2, 3).Pow(-1))
	}
	if a.Cmp(b) <= 0 || b.Cmp(a) >= 0 || a.Cmp(a) != 0 {
		t.Errorf("Cmp wrong")
	}
	if a.Mediant(b) != (Frac{2, 5}) {
		t.Errorf("Mediant = %v", a.Mediant(b))
	}
}

func TestFracFromFloat(t *testing.T) {
	f, ok := FracFromFloat(0.375)
	if !ok || f != (Frac{3, 8}) {
		t.Errorf("FracFromFloat(0.375) = %v,%v", f, ok)
	}
	// Every finite float64 is a dyadic rational: math.Pi round-trips exactly.
	if pf, ok := FracFromFloat(math.Pi); !ok || pf.Float() != math.Pi {
		t.Errorf("FracFromFloat(pi) = %v,%v", pf, ok)
	}
	// A value too large to fit an int64 numerator is not representable.
	if _, ok := FracFromFloat(1e300); ok {
		t.Errorf("1e300 should not be representable in int64")
	}
}

func TestHelpers(t *testing.T) {
	if GCD(48, 36) != 12 || LCM(4, 6) != 12 {
		t.Errorf("GCD/LCM wrong")
	}
	if Isqrt(1000000) != 1000 || Isqrt(0) != 0 || Isqrt(15) != 3 {
		t.Errorf("Isqrt wrong")
	}
	if !IsPerfectSquare(49) || IsPerfectSquare(50) {
		t.Errorf("IsPerfectSquare wrong")
	}
	if EulerPhi(36) != 12 || EulerPhi(1) != 1 || EulerPhi(13) != 12 {
		t.Errorf("EulerPhi wrong: %d", EulerPhi(36))
	}
	np, nq := ReduceFraction(-10, -4)
	if np != 5 || nq != 2 {
		t.Errorf("ReduceFraction = %d/%d", np, nq)
	}
}
