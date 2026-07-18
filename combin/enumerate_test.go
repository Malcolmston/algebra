package combin

import (
	"math"
	"math/big"
	"testing"
)

// ---- small helpers -----------------------------------------------------------

func bi(n int64) *big.Int { return big.NewInt(n) }

func eqBig(t *testing.T, got *big.Int, want int64, ctx string) {
	t.Helper()
	if got.Cmp(bi(want)) != 0 {
		t.Errorf("%s = %s, want %d", ctx, got.String(), want)
	}
}

func onesFromMask(g uint64) int {
	c := 0
	for g != 0 {
		c += int(g & 1)
		g >>= 1
	}
	return c
}

// ---- partitions --------------------------------------------------------------

func TestPartitionNumber(t *testing.T) {
	want := []int{1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42}
	for n, w := range want {
		if got := PartitionNumber(n); got != w {
			t.Errorf("PartitionNumber(%d) = %d, want %d", n, got, w)
		}
	}
	tbl := PartitionNumberTable(10)
	for n, w := range want {
		if tbl[n] != w {
			t.Errorf("PartitionNumberTable[%d] = %d, want %d", n, tbl[n], w)
		}
	}
}

func TestPartitionEnumerationConsistency(t *testing.T) {
	for n := 0; n <= 12; n++ {
		if got := len(Partitions(n)); got != PartitionNumber(n) {
			t.Errorf("len(Partitions(%d)) = %d, want %d", n, got, PartitionNumber(n))
		}
		// every partition must sum to n and be non-increasing
		for _, p := range Partitions(n) {
			s := 0
			for i, x := range p {
				s += x
				if i > 0 && p[i-1] < x {
					t.Errorf("Partitions(%d): %v not non-increasing", n, p)
				}
			}
			if s != n {
				t.Errorf("Partitions(%d): %v sums to %d", n, p, s)
			}
		}
	}
}

func TestNextPartitionCoversAll(t *testing.T) {
	for n := 1; n <= 12; n++ {
		p := []int{n}
		count := 1
		for {
			var ok bool
			p, ok = NextPartition(p)
			if !ok {
				break
			}
			count++
		}
		if count != PartitionNumber(n) {
			t.Errorf("NextPartition walk for n=%d visited %d, want %d", n, count, PartitionNumber(n))
		}
	}
}

func TestPartitionNumberInto(t *testing.T) {
	for n := 0; n <= 12; n++ {
		for k := 0; k <= n; k++ {
			if got := len(PartitionsInto(n, k)); got != PartitionNumberInto(n, k) {
				t.Errorf("len(PartitionsInto(%d,%d)) = %d, want %d", n, k, got, PartitionNumberInto(n, k))
			}
		}
	}
	if PartitionNumberInto(5, 2) != 2 {
		t.Errorf("PartitionNumberInto(5,2) = %d, want 2", PartitionNumberInto(5, 2))
	}
	if PartitionNumberInto(7, 3) != 4 {
		t.Errorf("PartitionNumberInto(7,3) = %d, want 4", PartitionNumberInto(7, 3))
	}
}

func TestPartitionNumberAtMost(t *testing.T) {
	for n := 0; n <= 15; n++ {
		if PartitionNumberAtMost(n, n) != PartitionNumber(n) {
			t.Errorf("PartitionNumberAtMost(%d,%d) != PartitionNumber(%d)", n, n, n)
		}
	}
}

func TestDistinctEqualsOdd(t *testing.T) {
	for n := 0; n <= 40; n++ {
		if PartitionNumberDistinct(n) != PartitionNumberOdd(n) {
			t.Errorf("Euler theorem fails at n=%d: distinct=%d odd=%d",
				n, PartitionNumberDistinct(n), PartitionNumberOdd(n))
		}
	}
	if PartitionNumberDistinct(6) != 4 {
		t.Errorf("PartitionNumberDistinct(6) = %d, want 4", PartitionNumberDistinct(6))
	}
	if got := len(PartitionsDistinct(6)); got != 4 {
		t.Errorf("len(PartitionsDistinct(6)) = %d, want 4", got)
	}
}

func TestPartitionShape(t *testing.T) {
	p := []int{4, 2, 1}
	c := PartitionConjugate(p)
	wantConj := []int{3, 2, 1, 1}
	if len(c) != len(wantConj) {
		t.Fatalf("conjugate = %v, want %v", c, wantConj)
	}
	for i := range c {
		if c[i] != wantConj[i] {
			t.Fatalf("conjugate = %v, want %v", c, wantConj)
		}
	}
	// conjugation is an involution
	cc := PartitionConjugate(c)
	if len(cc) != len(p) {
		t.Fatalf("double conjugate = %v, want %v", cc, p)
	}
	for i := range p {
		if cc[i] != p[i] {
			t.Fatalf("double conjugate = %v, want %v", cc, p)
		}
	}
	if PartitionRank(p) != 4-3 {
		t.Errorf("PartitionRank(%v) = %d, want 1", p, PartitionRank(p))
	}
	if PartitionDurfeeSquare(p) != 2 {
		t.Errorf("PartitionDurfeeSquare(%v) = %d, want 2", p, PartitionDurfeeSquare(p))
	}
}

func TestPentagonal(t *testing.T) {
	genWant := []int{0, 1, 2, 5, 7, 12, 15, 22, 26, 35}
	for i, w := range genWant {
		if got := GeneralizedPentagonal(i); got != w {
			t.Errorf("GeneralizedPentagonal(%d) = %d, want %d", i, got, w)
		}
	}
	pentWant := []int{0, 1, 5, 12, 22, 35}
	for k, w := range pentWant {
		if got := PentagonalNumber(k); got != w {
			t.Errorf("PentagonalNumber(%d) = %d, want %d", k, got, w)
		}
	}
}

// ---- compositions ------------------------------------------------------------

func TestCompositions(t *testing.T) {
	for n := 0; n <= 10; n++ {
		if got := len(CompositionList(n)); got != CompositionNumber(n) {
			t.Errorf("len(CompositionList(%d)) = %d, want %d", n, got, CompositionNumber(n))
		}
		for k := 0; k <= n; k++ {
			if got := len(CompositionListInto(n, k)); got != CompositionNumberInto(n, k) {
				t.Errorf("len(CompositionListInto(%d,%d)) = %d, want %d",
					n, k, got, CompositionNumberInto(n, k))
			}
		}
	}
	if CompositionNumber(5) != 16 {
		t.Errorf("CompositionNumber(5) = %d, want 16", CompositionNumber(5))
	}
}

func TestWeakCompositions(t *testing.T) {
	for n := 0; n <= 8; n++ {
		for k := 1; k <= 5; k++ {
			if got := len(WeakCompositionList(n, k)); got != WeakCompositionNumber(n, k) {
				t.Errorf("len(WeakCompositionList(%d,%d)) = %d, want %d",
					n, k, got, WeakCompositionNumber(n, k))
			}
		}
	}
}

// ---- derangements / involutions ---------------------------------------------

func TestDerangements(t *testing.T) {
	want := []int{1, 0, 1, 2, 9, 44, 265, 1854}
	for n, w := range want {
		if got := DerangementNumber(n); got != w {
			t.Errorf("DerangementNumber(%d) = %d, want %d", n, got, w)
		}
	}
	for n := 0; n <= 7; n++ {
		if got := len(DerangementList(n)); got != DerangementNumber(n) {
			t.Errorf("len(DerangementList(%d)) = %d, want %d", n, got, DerangementNumber(n))
		}
	}
	if math.Abs(DerangementProbability(12)-math.Exp(-1)) > 1e-6 {
		t.Errorf("DerangementProbability(12) = %v, want ~1/e", DerangementProbability(12))
	}
}

func TestInvolutions(t *testing.T) {
	want := []int{1, 1, 2, 4, 10, 26, 76}
	for n, w := range want {
		if got := InvolutionNumber(n); got != w {
			t.Errorf("InvolutionNumber(%d) = %d, want %d", n, got, w)
		}
	}
	for n := 0; n <= 6; n++ {
		if got := len(InvolutionList(n)); got != InvolutionNumber(n) {
			t.Errorf("len(InvolutionList(%d)) = %d, want %d", n, got, InvolutionNumber(n))
		}
	}
}

// ---- permutations ------------------------------------------------------------

func factI(n int) int {
	r := 1
	for i := 2; i <= n; i++ {
		r *= i
	}
	return r
}

func TestNextPermutationCycle(t *testing.T) {
	for n := 0; n <= 6; n++ {
		a := make([]int, n)
		for i := range a {
			a[i] = i
		}
		count := 1
		for NextPermutation(a) {
			count++
		}
		if count != factI(n) {
			t.Errorf("NextPermutation cycle n=%d visited %d, want %d", n, count, factI(n))
		}
		// after full cycle we wrapped back to ascending
		for i := range a {
			if a[i] != i {
				t.Errorf("NextPermutation did not wrap to sorted for n=%d: %v", n, a)
				break
			}
		}
	}
}

func TestPrevPermutation(t *testing.T) {
	a := []int{0, 1, 2}
	if PrevPermutation(a) {
		t.Errorf("PrevPermutation of lowest returned true")
	}
	// lowest wraps to highest
	want := []int{2, 1, 0}
	for i := range a {
		if a[i] != want[i] {
			t.Fatalf("wrap = %v, want %v", a, want)
		}
	}
	// next then prev is identity
	b := []int{1, 0, 2}
	c := []int{1, 0, 2}
	NextPermutation(b)
	PrevPermutation(b)
	for i := range b {
		if b[i] != c[i] {
			t.Fatalf("next/prev not identity: %v vs %v", b, c)
		}
	}
}

func TestPermutationRankUnrank(t *testing.T) {
	n := 5
	total := factI(n)
	for r := 0; r < total; r++ {
		p := PermutationUnrank(n, big.NewInt(int64(r)))
		got := PermutationRank(p)
		if got.Cmp(big.NewInt(int64(r))) != 0 {
			t.Errorf("rank(unrank(%d)) = %s", r, got.String())
		}
	}
	if PermutationParity([]int{0, 1, 2}) != 1 {
		t.Errorf("parity identity should be +1")
	}
	if PermutationParity([]int{1, 0, 2}) != -1 {
		t.Errorf("parity single swap should be -1")
	}
	if PermutationParity([]int{2, 0, 1}) != 1 {
		t.Errorf("parity 3-cycle should be +1")
	}
}

func TestMultisetPermutations(t *testing.T) {
	got := len(MultisetPermutations([]int{1, 1, 2}))
	if got != 3 {
		t.Errorf("distinct perms of {1,1,2} = %d, want 3", got)
	}
	eqBig(t, MultisetPermutationNumber([]int{2, 1}), 3, "MultisetPermutationNumber([2,1])")
	eqBig(t, MultisetPermutationNumber([]int{1, 1, 1, 1}), 24, "MultisetPermutationNumber(4 singletons)")
	eqBig(t, MultisetPermutationNumber([]int{2, 2, 2}), 90, "MultisetPermutationNumber([2,2,2])")
	if got := len(MultisetPermutationList([]int{2, 1})); got != 3 {
		t.Errorf("MultisetPermutationList([2,1]) len = %d, want 3", got)
	}
	if got := len(AllPermutations([]int{5, 6, 7, 8})); got != 24 {
		t.Errorf("AllPermutations of 4 = %d, want 24", got)
	}
}

func TestFactoradic(t *testing.T) {
	// 5! rank 349 -> known factoradic in a 5-digit factorial base
	digits := FactorialNumberSystem(big.NewInt(349), 5)
	// verify reconstruction
	r := 0
	for i := 0; i < 5; i++ {
		r = r*(5-i) + digits[i]
	}
	if r != 349 {
		t.Errorf("factoradic reconstruction = %d, want 349", r)
	}
}

// ---- combinations / subsets --------------------------------------------------

func TestCombinationList(t *testing.T) {
	for n := 0; n <= 8; n++ {
		for k := 0; k <= n; k++ {
			want := int(combinBinom(n, k).Int64())
			if got := len(CombinationList(n, k)); got != want {
				t.Errorf("len(CombinationList(%d,%d)) = %d, want %d", n, k, got, want)
			}
		}
	}
	// combinations of items
	got := CombinationListOf([]int{10, 20, 30}, 2)
	if len(got) != 3 {
		t.Fatalf("CombinationListOf len = %d, want 3", len(got))
	}
}

func TestCombinationRankUnrank(t *testing.T) {
	n, k := 8, 4
	for _, c := range CombinationList(n, k) {
		rank := CombinationRank(c)
		back := CombinationUnrank(k, rank)
		for i := range c {
			if back[i] != c[i] {
				t.Fatalf("unrank(rank(%v)) = %v", c, back)
			}
		}
	}
}

func TestPowerSet(t *testing.T) {
	for n := 0; n <= 10; n++ {
		want := int(PowerSetSize(n).Int64())
		if got := len(SubsetList(n)); got != want {
			t.Errorf("len(SubsetList(%d)) = %d, want %d", n, got, want)
		}
	}
	if got := len(PowerSet([]int{1, 2, 3})); got != 8 {
		t.Errorf("PowerSet of 3 = %d, want 8", got)
	}
	eqBig(t, MultisetCoefficient(3, 2), 6, "MultisetCoefficient(3,2)")
}

// ---- gray codes --------------------------------------------------------------

func TestGrayCode(t *testing.T) {
	for i := uint64(0); i < 256; i++ {
		g := GrayEncode(i)
		if GrayDecode(g) != i {
			t.Errorf("GrayDecode(GrayEncode(%d)) = %d", i, GrayDecode(g))
		}
		if GrayCodeRank(GrayCodeAt(i)) != i {
			t.Errorf("GrayCodeRank/At round trip failed at %d", i)
		}
	}
	seq := GrayCodeSequence(6)
	if len(seq) != 64 {
		t.Fatalf("GrayCodeSequence(6) len = %d, want 64", len(seq))
	}
	for i := 1; i < len(seq); i++ {
		if onesFromMask(seq[i]^seq[i-1]) != 1 {
			t.Errorf("gray codes %d and %d differ in !=1 bit", i-1, i)
		}
	}
	if onesFromMask(seq[0]^seq[len(seq)-1]) != 1 {
		t.Errorf("gray sequence is not cyclic in one bit")
	}
}

// ---- classical numbers -------------------------------------------------------

func TestCatalanAndNarayana(t *testing.T) {
	catWant := []int64{1, 1, 2, 5, 14, 42, 132, 429, 1430}
	cat := CatalanSequence(8)
	for n, w := range catWant {
		eqBig(t, cat[n], w, "Catalan")
	}
	// Narayana row sums to Catalan
	for n := 1; n <= 8; n++ {
		sum := big.NewInt(0)
		for _, v := range NarayanaTriangleRow(n) {
			sum.Add(sum, v)
		}
		if sum.Cmp(cat[n]) != 0 {
			t.Errorf("Narayana row %d sum = %s, want %s", n, sum.String(), cat[n].String())
		}
	}
	tri := NarayanaTriangle(4)
	if len(tri) != 5 {
		t.Errorf("NarayanaTriangle(4) rows = %d, want 5", len(tri))
	}
}

func TestMotzkin(t *testing.T) {
	want := []int64{1, 1, 2, 4, 9, 21, 51, 127}
	seq := MotzkinSequence(7)
	for n, w := range want {
		eqBig(t, seq[n], w, "Motzkin")
	}
}

func TestEulerianTriangle(t *testing.T) {
	row3 := EulerianRow(3)
	want := []int64{1, 4, 1}
	if len(row3) != 3 {
		t.Fatalf("EulerianRow(3) len = %d", len(row3))
	}
	for i, w := range want {
		eqBig(t, row3[i], w, "EulerianRow3")
	}
	// each row sums to n!
	for n := 1; n <= 7; n++ {
		sum := big.NewInt(0)
		for _, v := range EulerianRow(n) {
			sum.Add(sum, v)
		}
		if sum.Cmp(combinFactorial(n)) != 0 {
			t.Errorf("Eulerian row %d sum = %s, want %s", n, sum.String(), combinFactorial(n).String())
		}
	}
	if len(EulerianTriangle(5)) != 6 {
		t.Errorf("EulerianTriangle(5) rows wrong")
	}
}

func TestBellTriangle(t *testing.T) {
	bellWant := []int64{1, 1, 2, 5, 15, 52, 203}
	tri := BellTriangleFull(6)
	for n, w := range bellWant {
		eqBig(t, tri[n][0], w, "Bell")
	}
}

func TestBernoulli(t *testing.T) {
	type kv struct {
		n    int
		p, q int64
	}
	want := []kv{
		{0, 1, 1}, {1, -1, 2}, {2, 1, 6}, {3, 0, 1}, {4, -1, 30},
		{5, 0, 1}, {6, 1, 42}, {7, 0, 1}, {8, -1, 30}, {10, 5, 66},
	}
	for _, w := range want {
		got := BernoulliNumber(w.n)
		exp := big.NewRat(w.p, w.q)
		if got.Cmp(exp) != 0 {
			t.Errorf("BernoulliNumber(%d) = %s, want %s", w.n, got.String(), exp.String())
		}
	}
	if math.Abs(BernoulliNumberFloat(2)-1.0/6.0) > 1e-12 {
		t.Errorf("BernoulliNumberFloat(2) = %v, want 1/6", BernoulliNumberFloat(2))
	}
	if len(BernoulliSequence(8)) != 9 {
		t.Errorf("BernoulliSequence(8) length wrong")
	}
}

func TestEulerAndZigzag(t *testing.T) {
	eWant := map[int]int64{0: 1, 1: 0, 2: -1, 3: 0, 4: 5, 6: -61, 8: 1385}
	for n, w := range eWant {
		eqBig(t, EulerNumber(n), w, "EulerNumber")
	}
	zigWant := []int64{1, 1, 1, 2, 5, 16, 61, 272, 1385, 7936}
	for n, w := range zigWant {
		eqBig(t, ZigzagNumber(n), w, "Zigzag")
	}
	tanWant := []int64{1, 2, 16, 272, 7936}
	for n, w := range tanWant {
		eqBig(t, TangentNumber(n), w, "Tangent")
	}
	secWant := []int64{1, 1, 5, 61, 1385}
	for n, w := range secWant {
		eqBig(t, SecantNumber(n), w, "Secant")
	}
	if len(EulerNumberSequence(8)) != 9 {
		t.Errorf("EulerNumberSequence(8) length wrong")
	}
}

func TestNecklaceAndLyndon(t *testing.T) {
	eqBig(t, NecklaceNumber(2, 2), 3, "NecklaceNumber(2,2)")
	eqBig(t, NecklaceNumber(3, 2), 4, "NecklaceNumber(3,2)")
	eqBig(t, NecklaceNumber(4, 2), 6, "NecklaceNumber(4,2)")
	eqBig(t, LyndonWordNumber(2, 2), 1, "LyndonWordNumber(2,2)")
	eqBig(t, LyndonWordNumber(3, 2), 2, "LyndonWordNumber(3,2)")
	eqBig(t, LyndonWordNumber(6, 2), 9, "LyndonWordNumber(6,2)")
}

// ---- benchmark of the heaviest routine --------------------------------------

func BenchmarkBernoulliSequence(b *testing.B) {
	for i := 0; i < b.N; i++ {
		_ = BernoulliSequence(80)
	}
}
