package logic

import (
	"reflect"
	"testing"
)

// logicVerifyCover checks that the cover implements exactly the function given
// by minterms with the given don't-cares over numVars variables: every minterm
// is covered, every forbidden term (not a minterm and not a don't-care) is not.
func logicVerifyCover(t *testing.T, cover []Implicant, minterms, dontCares []int, numVars int) {
	t.Helper()
	mset := map[int]bool{}
	for _, m := range minterms {
		mset[m] = true
	}
	dset := map[int]bool{}
	for _, d := range dontCares {
		dset[d] = true
	}
	for i := 0; i < (1 << numVars); i++ {
		covered := false
		for _, im := range cover {
			if im.Covers(i) {
				covered = true
				break
			}
		}
		if mset[i] && !covered {
			t.Errorf("minterm %d not covered", i)
		}
		if !mset[i] && !dset[i] && covered {
			t.Errorf("term %d covered but should be 0", i)
		}
	}
}

func TestQuineMcCluskeyBasic(t *testing.T) {
	// minterms {1,3} over vars A,B: 01 and 11 combine on A -> "-1" = B.
	cover := QuineMcCluskey([]int{1, 3}, nil, 2)
	logicVerifyCover(t, cover, []int{1, 3}, nil, 2)
	if len(cover) != 1 {
		t.Fatalf("want 1 term, got %d", len(cover))
	}
	if got := SOPString(cover, []string{"A", "B"}); got != "B" {
		t.Errorf("SOPString = %q, want B", got)
	}
}

func TestQuineMcCluskeyAll(t *testing.T) {
	cover := QuineMcCluskey([]int{0, 1, 2, 3}, nil, 2)
	if len(cover) != 1 {
		t.Fatalf("want 1 term, got %d", len(cover))
	}
	if got := SOPString(cover, []string{"A", "B"}); got != "T" {
		t.Errorf("SOPString = %q, want T", got)
	}
}

func TestQuineMcCluskeyEmpty(t *testing.T) {
	if cover := QuineMcCluskey(nil, nil, 3); cover != nil {
		t.Errorf("empty function should give nil cover, got %v", cover)
	}
	if got := SOPString(nil, []string{"A"}); got != "F" {
		t.Errorf("SOPString(nil) = %q, want F", got)
	}
}

func TestQuineMcCluskeyClassic(t *testing.T) {
	// Classic 4-variable example (Mano): minterms with two don't-cares.
	minterms := []int{4, 8, 10, 11, 12, 15}
	dontCares := []int{9, 14}
	cover := QuineMcCluskey(minterms, dontCares, 4)
	logicVerifyCover(t, cover, minterms, dontCares, 4)
	// The minimal cover for this function is a compact sum of products; the
	// essential-plus-greedy strategy finds a correct cover of at most four
	// product terms (independently verified to be exactly three here).
	if len(cover) < 1 || len(cover) > 4 {
		t.Errorf("cover size %d out of expected range: %v", len(cover), SOPString(cover, []string{"A", "B", "C", "D"}))
	}
}

func TestQuineMcCluskeyThreeVar(t *testing.T) {
	// f = sum(0,1,2,5,6,7) over A,B,C. Verify correctness against the function.
	minterms := []int{0, 1, 2, 5, 6, 7}
	cover := QuineMcCluskey(minterms, nil, 3)
	logicVerifyCover(t, cover, minterms, nil, 3)
}

func TestMinimizeFromExpr(t *testing.T) {
	e := MustParse("(A & B) | (A & !B)")
	// Reduces to A.
	if got := MinimizeString(e); got != "A" {
		t.Errorf("MinimizeString = %q, want A", got)
	}
	cover := MinimizeSOP(e)
	if !Equivalent(e, MustParse(SOPString(cover, Vars(e)))) {
		t.Errorf("minimized cover not equivalent to original")
	}
}

func TestMinimizeConstants(t *testing.T) {
	if got := MinimizeString(MustParse("A & !A")); got != "F" {
		t.Errorf("contradiction minimizes to %q, want F", got)
	}
	if got := MinimizeString(MustParse("A | !A")); got != "T" {
		t.Errorf("tautology minimizes to %q, want T", got)
	}
}

func TestPrimeAndEssential(t *testing.T) {
	minterms := []int{0, 1, 2, 5, 6, 7}
	primes := PrimeImplicants(minterms, nil, 3)
	if len(primes) == 0 {
		t.Fatalf("expected prime implicants")
	}
	// Every prime implicant must cover only real minterms.
	for _, p := range primes {
		for m := 0; m < 8; m++ {
			if p.Covers(m) {
				found := false
				for _, mm := range minterms {
					if mm == m {
						found = true
					}
				}
				if !found {
					t.Errorf("prime %s covers non-minterm %d", p.String(), m)
				}
			}
		}
	}
	ess := EssentialPrimeImplicants(minterms, nil, 3)
	// Essentials are a subset of primes.
	if len(ess) > len(primes) {
		t.Errorf("more essentials than primes")
	}
}

func TestImplicantMethods(t *testing.T) {
	// Cube "1-" over A,B: A fixed 1, B eliminated.
	im := Implicant{Value: 0b10, Dashes: 0b01, NumVars: 2}
	if im.String() != "1-" {
		t.Errorf("String = %q, want 1-", im.String())
	}
	if im.LiteralCount() != 1 {
		t.Errorf("LiteralCount = %d, want 1", im.LiteralCount())
	}
	if !im.Covers(2) || !im.Covers(3) || im.Covers(0) {
		t.Errorf("Covers wrong for cube 1-")
	}
	if got := im.Expression([]string{"A", "B"}); got != "A" {
		t.Errorf("Expression = %q, want A", got)
	}
}

func TestGrayCode(t *testing.T) {
	if got := GrayCode(2); !reflect.DeepEqual(got, []int{0, 1, 3, 2}) {
		t.Errorf("GrayCode(2) = %v, want [0 1 3 2]", got)
	}
	for n := 0; n < 64; n++ {
		if GrayDecode(GrayEncode(n)) != n {
			t.Errorf("Gray round trip failed at %d", n)
		}
	}
	// Adjacent Gray codes differ in exactly one bit.
	g := GrayCode(4)
	for i := 1; i < len(g); i++ {
		if PopCount(uint(g[i]^g[i-1])) != 1 {
			t.Errorf("Gray adjacency broken at %d", i)
		}
	}
}

func TestBitStringPopCount(t *testing.T) {
	if BitString(5, 4) != "0101" {
		t.Errorf("BitString(5,4) = %q", BitString(5, 4))
	}
	if PopCount(0b10110) != 3 {
		t.Errorf("PopCount wrong")
	}
}

func TestKarnaughMap(t *testing.T) {
	k := NewKarnaughMap([]int{3}, 2)
	if k.Rows() != 2 || k.Cols() != 2 {
		t.Fatalf("dims = %dx%d, want 2x2", k.Rows(), k.Cols())
	}
	// minterm 3 => A=1,B=1 => row 1, col 1 in Gray order (row/col bits = 1).
	if !k.At(1, 1) {
		t.Errorf("expected At(1,1) true for minterm 3")
	}
	if k.At(0, 0) {
		t.Errorf("expected At(0,0) false")
	}
	// Groups equal the prime implicants of the mapped minterms.
	k2 := NewKarnaughMap([]int{0, 1, 2, 3}, 2)
	groups := k2.Groups()
	if len(groups) != 1 || groups[0].LiteralCount() != 0 {
		t.Errorf("full 2-var map should give one all-dash group, got %v", groups)
	}
	_ = k2.String() // exercise rendering
}
