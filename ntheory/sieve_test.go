package ntheory

import (
	"reflect"
	"testing"
)

// ntheorySieveReference returns primes in [lo, hi] via the independent
// deterministic primality test, used to cross-check the segmented sieve.
func ntheorySieveReference(lo, hi uint64) []uint64 {
	var out []uint64
	if hi < 2 {
		return nil
	}
	if lo < 2 {
		lo = 2
	}
	for n := lo; n <= hi; n++ {
		if IsPrimeU64(n) {
			out = append(out, n)
		}
	}
	return out
}

func TestSieveSegmentedKnownAnswers(t *testing.T) {
	tests := []struct {
		lo, hi uint64
		want   []uint64
	}{
		{0, 1, nil},
		{5, 3, nil},
		{0, 2, []uint64{2}},
		{0, 3, []uint64{2, 3}},
		{0, 10, []uint64{2, 3, 5, 7}},
		{0, 30, []uint64{2, 3, 5, 7, 11, 13, 17, 19, 23, 29}},
		{2, 2, []uint64{2}},
		{3, 5, []uint64{3, 5}},
		{4, 6, []uint64{5}},
		{6, 6, nil},
		{7, 7, []uint64{7}},
		{8, 10, nil},
		{10, 30, []uint64{11, 13, 17, 19, 23, 29}},
		{90, 100, []uint64{97}},
		{100, 110, []uint64{101, 103, 107, 109}},
		{113, 127, []uint64{113, 127}},
	}
	for _, tt := range tests {
		got := SegmentedSieve(tt.lo, tt.hi)
		if !reflect.DeepEqual(got, tt.want) {
			t.Errorf("SegmentedSieve(%d, %d) = %v, want %v", tt.lo, tt.hi, got, tt.want)
		}
	}
}

func TestSieveSegmentedMatchesReference(t *testing.T) {
	// Cover many segment boundaries: ntheorySieveSpan is 245760, so ranges
	// spanning several spans exercise the multi-segment path.
	ranges := []struct{ lo, hi uint64 }{
		{0, 1000},
		{0, 100000},
		{1, 300000},
		{200000, 260000},
		{245000, 246000},
		{999900, 1000100},
	}
	for _, r := range ranges {
		got := SegmentedSieve(r.lo, r.hi)
		want := ntheorySieveReference(r.lo, r.hi)
		if !reflect.DeepEqual(got, want) {
			t.Errorf("SegmentedSieve(%d, %d) disagrees with reference (got %d primes, want %d)",
				r.lo, r.hi, len(got), len(want))
		}
	}
}

func TestSieveSegmentedMatchesPrimesUpTo(t *testing.T) {
	want := PrimesUpTo(5000)
	gotU := SegmentedSieve(0, 5000)
	got := make([]int64, len(gotU))
	for i, p := range gotU {
		got[i] = int64(p)
	}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("SegmentedSieve(0, 5000) disagrees with PrimesUpTo(5000)")
	}
}

func TestSieveHighNarrowWindow(t *testing.T) {
	// A narrow window far above the 32-bit range: memory use stays bounded.
	const lo, hi = uint64(1_000_000_000_000), uint64(1_000_000_000_100)
	got := SegmentedSieve(lo, hi)
	want := ntheorySieveReference(lo, hi)
	if !reflect.DeepEqual(got, want) {
		t.Errorf("SegmentedSieve high window = %v, want %v", got, want)
	}
}

func TestSievePrimesInRange(t *testing.T) {
	tests := []struct {
		lo, hi uint64
		want   []uint64
	}{
		{0, 0, nil},
		{0, 1, nil},
		{0, 2, nil},
		{0, 3, []uint64{2}},
		{2, 12, []uint64{2, 3, 5, 7, 11}},
		{10, 20, []uint64{11, 13, 17, 19}},
		{11, 13, []uint64{11}},
		{14, 14, nil},
	}
	for _, tt := range tests {
		got := PrimesInRange(tt.lo, tt.hi)
		if !reflect.DeepEqual(got, tt.want) {
			t.Errorf("PrimesInRange(%d, %d) = %v, want %v", tt.lo, tt.hi, got, tt.want)
		}
	}
}

func TestSievePrimePiRange(t *testing.T) {
	tests := []struct {
		lo, hi uint64
		want   uint64
	}{
		{0, 1, 0},
		{5, 3, 0},
		{0, 10, 4},
		{0, 100, 25},
		{0, 1000, 168},
		{0, 10000, 1229},
		{900, 1000, 14},
		{2, 2, 1},
		{200000, 260000, 4853}, // spans multiple segments
	}
	for _, tt := range tests {
		if got := PrimePiRange(tt.lo, tt.hi); got != tt.want {
			t.Errorf("PrimePiRange(%d, %d) = %d, want %d", tt.lo, tt.hi, got, tt.want)
		}
	}
}

func TestSievePrimePiRangeMatchesCount(t *testing.T) {
	ranges := []struct{ lo, hi uint64 }{
		{0, 50000},
		{100000, 200000},
		{245000, 250000},
	}
	for _, r := range ranges {
		want := uint64(len(SegmentedSieve(r.lo, r.hi)))
		if got := PrimePiRange(r.lo, r.hi); got != want {
			t.Errorf("PrimePiRange(%d, %d) = %d, want %d", r.lo, r.hi, got, want)
		}
	}
}

func TestSieveNthPrime(t *testing.T) {
	tests := []struct {
		n    uint64
		want uint64
	}{
		{1, 2},
		{2, 3},
		{3, 5},
		{4, 7},
		{5, 11},
		{6, 13},
		{10, 29},
		{25, 97},
		{100, 541},
		{168, 997},
		{1000, 7919},
		{10000, 104729},
	}
	for _, tt := range tests {
		if got := NthPrime(tt.n); got != tt.want {
			t.Errorf("NthPrime(%d) = %d, want %d", tt.n, got, tt.want)
		}
	}
}

func TestSieveNthPrimePanics(t *testing.T) {
	defer func() {
		if recover() == nil {
			t.Errorf("NthPrime(0) did not panic")
		}
	}()
	NthPrime(0)
}

func TestSievePrimeSieveStream(t *testing.T) {
	want := []uint64{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47}
	s := NewPrimeSieve()
	for i, w := range want {
		if got := s.Next(); got != w {
			t.Fatalf("PrimeSieve.Next() call %d = %d, want %d", i+1, got, w)
		}
	}
}

func TestSievePrimeSieveMatchesReference(t *testing.T) {
	const upTo = uint64(300000)
	want := ntheorySieveReference(2, upTo)
	s := NewPrimeSieve()
	for i, w := range want {
		got := s.Next()
		if got != w {
			t.Fatalf("PrimeSieve.Next() index %d = %d, want %d", i, got, w)
		}
	}
	// The next prime after upTo must be strictly greater than upTo.
	if next := s.Next(); next <= upTo {
		t.Fatalf("PrimeSieve.Next() after %d = %d, want > %d", upTo, next, upTo)
	}
}

func TestSievePrimeSieveAgreesWithNthPrime(t *testing.T) {
	s := NewPrimeSieve()
	for n := uint64(1); n <= 500; n++ {
		got := s.Next()
		if want := NthPrime(n); got != want {
			t.Fatalf("streaming prime %d = %d, NthPrime(%d) = %d", n, got, n, want)
		}
	}
}

func BenchmarkSieveSegmentedHigh(b *testing.B) {
	for i := 0; i < b.N; i++ {
		SegmentedSieve(1_000_000_000_000, 1_000_000_100_000)
	}
}

func BenchmarkSievePrimePiRange(b *testing.B) {
	for i := 0; i < b.N; i++ {
		PrimePiRange(0, 1_000_000)
	}
}

func BenchmarkSieveNthPrime(b *testing.B) {
	for i := 0; i < b.N; i++ {
		NthPrime(100_000)
	}
}

func BenchmarkSievePrimeSieveStream(b *testing.B) {
	for i := 0; i < b.N; i++ {
		s := NewPrimeSieve()
		for j := 0; j < 100_000; j++ {
			s.Next()
		}
	}
}
