package ntheory

import (
	"math/big"
	"testing"
)

func eqInt64Slice(a, b []int64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func TestContinuedFraction(t *testing.T) {
	cases := []struct {
		p, q int64
		want []int64
	}{
		{415, 93, []int64{4, 2, 6, 7}},
		{3, 2, []int64{1, 2}},
		{2, 3, []int64{0, 1, 2}},
		{7, 1, []int64{7}},
		{0, 5, []int64{0}},
		{4, 2, []int64{2}},
		{415, -93, []int64{-5, 1, 1, 6, 7}},
		{-415, 93, []int64{-5, 1, 1, 6, 7}},
	}
	for _, c := range cases {
		got := ContinuedFraction(c.p, c.q)
		if !eqInt64Slice(got, c.want) {
			t.Errorf("ContinuedFraction(%d, %d) = %v, want %v", c.p, c.q, got, c.want)
		}
		// Round-trip: the reconstructed rational must equal p/q.
		want := new(big.Rat).SetFrac64(c.p, c.q)
		if got := RatFromContinuedFraction(c.want); got.Cmp(want) != 0 {
			t.Errorf("RatFromContinuedFraction(%v) = %s, want %s", c.want, got, want)
		}
	}
}

func TestContinuedFractionPanicsOnZeroDenominator(t *testing.T) {
	defer func() {
		if recover() == nil {
			t.Error("ContinuedFraction(1, 0) did not panic")
		}
	}()
	ContinuedFraction(1, 0)
}

func TestContinuedFractionRat(t *testing.T) {
	cases := []struct {
		num, den int64
		want     []int64
	}{
		{415, 93, []int64{4, 2, 6, 7}},
		{3, 2, []int64{1, 2}},
		{-415, 93, []int64{-5, 1, 1, 6, 7}},
	}
	for _, c := range cases {
		r := new(big.Rat).SetFrac64(c.num, c.den)
		got := ContinuedFractionRat(r)
		if !eqInt64Slice(got, c.want) {
			t.Errorf("ContinuedFractionRat(%d/%d) = %v, want %v", c.num, c.den, got, c.want)
		}
	}
}

func TestConvergents(t *testing.T) {
	cf := []int64{4, 2, 6, 7}
	want := []*big.Rat{
		big.NewRat(4, 1),
		big.NewRat(9, 2),
		big.NewRat(58, 13),
		big.NewRat(415, 93),
	}
	got := Convergents(cf)
	if len(got) != len(want) {
		t.Fatalf("Convergents(%v) returned %d entries, want %d", cf, len(got), len(want))
	}
	for i := range want {
		if got[i].Cmp(want[i]) != 0 {
			t.Errorf("Convergents(%v)[%d] = %s, want %s", cf, i, got[i], want[i])
		}
	}
	if len(Convergents(nil)) != 0 {
		t.Error("Convergents(nil) should be empty")
	}
}

func TestRatFromContinuedFraction(t *testing.T) {
	if r := RatFromContinuedFraction(nil); r.Sign() != 0 {
		t.Errorf("RatFromContinuedFraction(nil) = %s, want 0", r)
	}
	// The last convergent equals the reconstructed value.
	cf := []int64{4, 2, 6, 7}
	if r := RatFromContinuedFraction(cf); r.Cmp(big.NewRat(415, 93)) != 0 {
		t.Errorf("RatFromContinuedFraction(%v) = %s, want 415/93", cf, r)
	}
}

func TestSqrtContinuedFraction(t *testing.T) {
	cases := []struct {
		n      uint64
		a0     int64
		period []int64
	}{
		{2, 1, []int64{2}},
		{3, 1, []int64{1, 2}},
		{7, 2, []int64{1, 1, 1, 4}},
		{23, 4, []int64{1, 3, 1, 8}},
		{13, 3, []int64{1, 1, 1, 1, 6}},
		{4, 2, nil},
		{9, 3, nil},
		{1, 1, nil},
	}
	for _, c := range cases {
		a0, period := SqrtContinuedFraction(c.n)
		if a0 != c.a0 || !eqInt64Slice(period, c.period) {
			t.Errorf("SqrtContinuedFraction(%d) = (%d, %v), want (%d, %v)",
				c.n, a0, period, c.a0, c.period)
		}
	}
}

func TestPellFundamental(t *testing.T) {
	cases := []struct {
		n    uint64
		x, y string
	}{
		{2, "3", "2"},
		{3, "2", "1"},
		{5, "9", "4"},
		{7, "8", "3"},
		{13, "649", "180"},
		{61, "1766319049", "226153980"},
	}
	for _, c := range cases {
		x, y := PellFundamental(c.n)
		if x.String() != c.x || y.String() != c.y {
			t.Errorf("PellFundamental(%d) = (%s, %s), want (%s, %s)",
				c.n, x.String(), y.String(), c.x, c.y)
		}
		// Verify x^2 - n*y^2 == 1 exactly.
		lhs := new(big.Int).Mul(x, x)
		ny2 := new(big.Int).Mul(y, y)
		ny2.Mul(ny2, new(big.Int).SetUint64(c.n))
		lhs.Sub(lhs, ny2)
		if lhs.Cmp(big.NewInt(1)) != 0 {
			t.Errorf("PellFundamental(%d): x^2 - n*y^2 = %s, want 1", c.n, lhs)
		}
	}
}

func TestPellFundamentalPanicsOnSquare(t *testing.T) {
	defer func() {
		if recover() == nil {
			t.Error("PellFundamental(9) did not panic on a perfect square")
		}
	}()
	PellFundamental(9)
}

func BenchmarkSqrtContinuedFraction(b *testing.B) {
	for i := 0; i < b.N; i++ {
		SqrtContinuedFraction(9949) // long period, non-square
	}
}

func BenchmarkConvergents(b *testing.B) {
	_, period := SqrtContinuedFraction(9949)
	cf := make([]int64, 0, len(period)+1)
	cf = append(cf, 99)
	cf = append(cf, period...)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Convergents(cf)
	}
}

func BenchmarkPellFundamental(b *testing.B) {
	for i := 0; i < b.N; i++ {
		PellFundamental(61)
	}
}
