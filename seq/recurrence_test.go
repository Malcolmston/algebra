package seq

import (
	"math/big"
	"reflect"
	"testing"
)

func TestFibonacci(t *testing.T) {
	// Known values F0..F20 plus a larger checkpoint.
	want := []uint64{0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181, 6765}
	for n, w := range want {
		if got := Fibonacci(n); got != w {
			t.Errorf("Fibonacci(%d) = %d, want %d", n, got, w)
		}
	}
	if got := Fibonacci(50); got != 12586269025 {
		t.Errorf("Fibonacci(50) = %d, want 12586269025", got)
	}
	// Largest exact 64-bit Fibonacci number.
	if got := Fibonacci(93); got != 12200160415121876738 {
		t.Errorf("Fibonacci(93) = %d, want 12200160415121876738", got)
	}
}

func TestFibonacciBig(t *testing.T) {
	// F100 is a standard reference value.
	want, _ := new(big.Int).SetString("354224848179261915075", 10)
	if got := FibonacciBig(100); got.Cmp(want) != 0 {
		t.Errorf("FibonacciBig(100) = %s, want %s", got, want)
	}
	// Consistency with the uint64 fast path across the exact range.
	for n := 0; n <= 93; n++ {
		if FibonacciBig(n).Uint64() != Fibonacci(n) {
			t.Errorf("FibonacciBig/Fibonacci disagree at n=%d", n)
		}
	}
}

func TestFibonacciSequence(t *testing.T) {
	got := FibonacciSequence(10)
	want := []uint64{0, 1, 1, 2, 3, 5, 8, 13, 21, 34}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("FibonacciSequence(10) = %v, want %v", got, want)
	}
}

func TestLucas(t *testing.T) {
	want := []uint64{2, 1, 3, 4, 7, 11, 18, 29, 47, 76, 123, 199, 322}
	for n, w := range want {
		if got := Lucas(n); got != w {
			t.Errorf("Lucas(%d) = %d, want %d", n, got, w)
		}
		if got := LucasBig(n); got.Uint64() != w {
			t.Errorf("LucasBig(%d) = %s, want %d", n, got, w)
		}
	}
	// Lucas relation L(n) = F(n-1) + F(n+1).
	for n := 1; n <= 80; n++ {
		if Lucas(n) != Fibonacci(n-1)+Fibonacci(n+1) {
			t.Errorf("Lucas identity failed at n=%d", n)
		}
	}
}

func TestPell(t *testing.T) {
	want := []uint64{0, 1, 2, 5, 12, 29, 70, 169, 408, 985, 2378}
	for n, w := range want {
		if got := Pell(n); got != w {
			t.Errorf("Pell(%d) = %d, want %d", n, got, w)
		}
		if got := PellBig(n); got.Uint64() != w {
			t.Errorf("PellBig(%d) = %s, want %d", n, got, w)
		}
	}
	if !reflect.DeepEqual(PellSequence(6), []uint64{0, 1, 2, 5, 12, 29}) {
		t.Errorf("PellSequence(6) mismatch")
	}
}

func TestPellLucas(t *testing.T) {
	want := []uint64{2, 2, 6, 14, 34, 82, 198, 478}
	for n, w := range want {
		if got := PellLucas(n); got != w {
			t.Errorf("PellLucas(%d) = %d, want %d", n, got, w)
		}
	}
}

func TestJacobsthal(t *testing.T) {
	want := []uint64{0, 1, 1, 3, 5, 11, 21, 43, 85, 171, 341}
	for n, w := range want {
		if got := Jacobsthal(n); got != w {
			t.Errorf("Jacobsthal(%d) = %d, want %d", n, got, w)
		}
		// Closed form (2^n - (-1)^n)/3.
		var pow uint64 = 1
		for i := 0; i < n; i++ {
			pow *= 2
		}
		var cf uint64
		if n%2 == 0 {
			cf = (pow - 1) / 3
		} else {
			cf = (pow + 1) / 3
		}
		if cf != w {
			t.Errorf("Jacobsthal closed form check failed at n=%d", n)
		}
	}
	if !reflect.DeepEqual(JacobsthalSequence(5), []uint64{0, 1, 1, 3, 5}) {
		t.Errorf("JacobsthalSequence(5) mismatch")
	}
}

func TestJacobsthalLucas(t *testing.T) {
	want := []uint64{2, 1, 5, 7, 17, 31, 65, 127}
	for n, w := range want {
		if got := JacobsthalLucas(n); got != w {
			t.Errorf("JacobsthalLucas(%d) = %d, want %d", n, got, w)
		}
	}
}

func TestTribonacci(t *testing.T) {
	want := []uint64{0, 0, 1, 1, 2, 4, 7, 13, 24, 44, 81, 149, 274}
	for n, w := range want {
		if got := Tribonacci(n); got != w {
			t.Errorf("Tribonacci(%d) = %d, want %d", n, got, w)
		}
	}
	if !reflect.DeepEqual(TribonacciSequence(7), []uint64{0, 0, 1, 1, 2, 4, 7}) {
		t.Errorf("TribonacciSequence(7) mismatch")
	}
}

func TestTetranacci(t *testing.T) {
	want := []uint64{0, 0, 0, 1, 1, 2, 4, 8, 15, 29, 56, 108}
	for n, w := range want {
		if got := Tetranacci(n); got != w {
			t.Errorf("Tetranacci(%d) = %d, want %d", n, got, w)
		}
	}
}

func TestPadovan(t *testing.T) {
	want := []uint64{1, 1, 1, 2, 2, 3, 4, 5, 7, 9, 12, 16, 21, 28, 37}
	for n, w := range want {
		if got := Padovan(n); got != w {
			t.Errorf("Padovan(%d) = %d, want %d", n, got, w)
		}
	}
	if !reflect.DeepEqual(PadovanSequence(6), []uint64{1, 1, 1, 2, 2, 3}) {
		t.Errorf("PadovanSequence(6) mismatch")
	}
}

func TestPerrin(t *testing.T) {
	want := []uint64{3, 0, 2, 3, 2, 5, 5, 7, 10, 12, 17, 22, 29, 39}
	for n, w := range want {
		if got := Perrin(n); got != w {
			t.Errorf("Perrin(%d) = %d, want %d", n, got, w)
		}
	}
	if !reflect.DeepEqual(PerrinSequence(6), []uint64{3, 0, 2, 3, 2, 5}) {
		t.Errorf("PerrinSequence(6) mismatch")
	}
}

func TestLinearRecurrence(t *testing.T) {
	// Fibonacci via the general engine.
	fibC := []int64{1, 1}
	fibI := []int64{0, 1}
	for n := 0; n <= 30; n++ {
		if got := LinearRecurrence(fibC, fibI, n); got != int64(Fibonacci(n)) {
			t.Errorf("LinearRecurrence Fibonacci mismatch at n=%d: %d", n, got)
		}
	}
	// Perrin has a negative-free recurrence a(n)=a(n-2)+a(n-3); express it as
	// order 3 with a leading zero coefficient.
	perC := []int64{0, 1, 1}
	perI := []int64{3, 0, 2}
	for n := 0; n <= 20; n++ {
		if got := LinearRecurrence(perC, perI, n); got != int64(Perrin(n)) {
			t.Errorf("LinearRecurrence Perrin mismatch at n=%d: %d", n, got)
		}
	}
	// Sequence and big variants agree.
	gotSeq := LinearRecurrenceSequence(fibC, fibI, 12)
	wantSeq := []int64{0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89}
	if !reflect.DeepEqual(gotSeq, wantSeq) {
		t.Errorf("LinearRecurrenceSequence mismatch: %v", gotSeq)
	}
	bC := []*big.Int{big.NewInt(1), big.NewInt(1)}
	bI := []*big.Int{big.NewInt(0), big.NewInt(1)}
	if got := LinearRecurrenceBig(bC, bI, 100); got.Cmp(FibonacciBig(100)) != 0 {
		t.Errorf("LinearRecurrenceBig(100) = %s, want F100", got)
	}
}

func TestIsFibonacci(t *testing.T) {
	fibSet := map[uint64]bool{}
	for n := 0; n <= 90; n++ {
		fibSet[Fibonacci(n)] = true
	}
	for x := uint64(0); x <= 200; x++ {
		if IsFibonacci(x) != fibSet[x] {
			t.Errorf("IsFibonacci(%d) = %v, want %v", x, IsFibonacci(x), fibSet[x])
		}
	}
}

func TestZeckendorf(t *testing.T) {
	cases := map[uint64][]uint64{
		0:   {},
		1:   {1},
		10:  {8, 2},
		11:  {8, 3},
		100: {89, 8, 3},
		64:  {55, 8, 1},
	}
	for n, want := range cases {
		got := Zeckendorf(n)
		if len(want) == 0 {
			if len(got) != 0 {
				t.Errorf("Zeckendorf(%d) = %v, want empty", n, got)
			}
			continue
		}
		if !reflect.DeepEqual(got, want) {
			t.Errorf("Zeckendorf(%d) = %v, want %v", n, got, want)
		}
	}
	// Property: terms sum to n, are Fibonacci, strictly decreasing and
	// non-consecutive in index.
	for n := uint64(1); n <= 500; n++ {
		terms := Zeckendorf(n)
		var sum uint64
		for i, tm := range terms {
			if !IsFibonacci(tm) {
				t.Errorf("Zeckendorf(%d) term %d not Fibonacci", n, tm)
			}
			if i > 0 && terms[i-1] <= tm {
				t.Errorf("Zeckendorf(%d) not strictly decreasing", n)
			}
			sum += tm
		}
		if sum != n {
			t.Errorf("Zeckendorf(%d) sums to %d", n, sum)
		}
	}
}

func BenchmarkFibonacciBig(b *testing.B) {
	for i := 0; i < b.N; i++ {
		_ = FibonacciBig(100000)
	}
}
