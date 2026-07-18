package crypto

import (
	"math/big"
	"testing"
)

func TestPrimeFactors(t *testing.T) {
	cases := []struct {
		n    int64
		want []int64
	}{
		{360, []int64{2, 2, 2, 3, 3, 5}},
		{97, []int64{97}},
		{1, nil},
		{100, []int64{2, 2, 5, 5}},
		{1024, []int64{2, 2, 2, 2, 2, 2, 2, 2, 2, 2}},
		{13195, []int64{5, 7, 13, 29}},
	}
	for _, c := range cases {
		got := PrimeFactors(bi(c.n))
		if len(got) != len(c.want) {
			t.Errorf("PrimeFactors(%d) len=%d want %d (%v)", c.n, len(got), len(c.want), got)
			continue
		}
		for i := range c.want {
			if got[i].Int64() != c.want[i] {
				t.Errorf("PrimeFactors(%d)[%d]=%d want %d", c.n, i, got[i].Int64(), c.want[i])
			}
		}
		// Product must reconstruct n.
		prod := big.NewInt(1)
		for _, f := range got {
			prod.Mul(prod, f)
		}
		if c.n > 1 && prod.Int64() != c.n {
			t.Errorf("PrimeFactors(%d) product=%d", c.n, prod.Int64())
		}
	}
}

func TestFactorizationSemiprime(t *testing.T) {
	// A semiprime that forces Pollard rho past trial division.
	n := bi(1000003 * 1000033)
	fs := Factorization(n)
	if len(fs) != 2 {
		t.Fatalf("Factorization semiprime got %d factors", len(fs))
	}
	if fs[0].Prime.Int64() != 1000003 || fs[1].Prime.Int64() != 1000033 {
		t.Errorf("Factorization semiprime = %v", fs)
	}
}

func TestPollardRho(t *testing.T) {
	composites := []int64{8051, 10403, 1234567, 999999937 * 2}
	for _, c := range composites {
		n := bi(c)
		d := PollardRho(n)
		if d.Int64() <= 1 || d.Int64() >= c || c%d.Int64() != 0 {
			t.Errorf("PollardRho(%d)=%d not a proper factor", c, d.Int64())
		}
		db := PollardRhoBrent(n)
		if db.Int64() <= 1 || c%db.Int64() != 0 {
			t.Errorf("PollardRhoBrent(%d)=%d not a proper factor", c, db.Int64())
		}
	}
}

func TestEulerTotient(t *testing.T) {
	cases := []struct{ n, want int64 }{
		{1, 1}, {2, 1}, {9, 6}, {10, 4}, {36, 12}, {97, 96}, {100, 40},
	}
	for _, c := range cases {
		if got := EulerTotient(bi(c.n)); got.Int64() != c.want {
			t.Errorf("EulerTotient(%d)=%d want %d", c.n, got.Int64(), c.want)
		}
	}
}

func TestCarmichaelLambda(t *testing.T) {
	cases := []struct{ n, want int64 }{
		{1, 1}, {2, 1}, {4, 2}, {8, 2}, {15, 4}, {21, 6}, {36, 6}, {561, 80},
	}
	for _, c := range cases {
		if got := CarmichaelLambda(bi(c.n)); got.Int64() != c.want {
			t.Errorf("CarmichaelLambda(%d)=%d want %d", c.n, got.Int64(), c.want)
		}
	}
}

func TestSmallestPrimeFactor(t *testing.T) {
	cases := []struct{ n, want int64 }{{15, 3}, {49, 7}, {97, 97}, {1000003 * 1000033, 1000003}}
	for _, c := range cases {
		if got := SmallestPrimeFactor(bi(c.n)); got.Int64() != c.want {
			t.Errorf("SmallestPrimeFactor(%d)=%d want %d", c.n, got.Int64(), c.want)
		}
	}
}
