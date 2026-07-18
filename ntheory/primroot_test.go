package ntheory

import (
	"reflect"
	"testing"
)

func TestCarmichael(t *testing.T) {
	cases := []struct {
		n    uint64
		want uint64
	}{
		{0, 0},
		{1, 1},
		{2, 1},
		{3, 2},
		{4, 2},
		{5, 4},
		{7, 6},
		{8, 2},
		{12, 2},
		{15, 4},
		{16, 4},
		{32, 8},
		{35, 12},
		{100, 20},
		{561, 80},
		{1105, 48},
	}
	for _, c := range cases {
		if got := Carmichael(c.n); got != c.want {
			t.Errorf("Carmichael(%d) = %d, want %d", c.n, got, c.want)
		}
	}
}

func TestMultiplicativeOrder(t *testing.T) {
	cases := []struct {
		a, m   uint64
		want   uint64
		wantOK bool
	}{
		{1, 7, 1, true},
		{2, 7, 3, true},
		{3, 7, 6, true},
		{4, 7, 3, true},
		{2, 5, 4, true},
		{2, 9, 6, true},
		{3, 10, 4, true},
		{10, 3, 1, true}, // 10 ≡ 1 (mod 3)
		{1, 1, 1, true},
		{2, 4, 0, false}, // gcd(2,4) != 1
		{6, 9, 0, false}, // gcd(6,9) != 1
	}
	for _, c := range cases {
		got, ok := MultiplicativeOrder(c.a, c.m)
		if got != c.want || ok != c.wantOK {
			t.Errorf("MultiplicativeOrder(%d, %d) = (%d, %v), want (%d, %v)",
				c.a, c.m, got, ok, c.want, c.wantOK)
		}
	}
}

// TestMultiplicativeOrderMatchesDefinition brute-forces the order and compares.
func TestMultiplicativeOrderMatchesDefinition(t *testing.T) {
	for m := uint64(2); m <= 200; m++ {
		for a := uint64(0); a < m; a++ {
			got, ok := MultiplicativeOrder(a, m)
			if ntheoryGCDU64(a, m) != 1 {
				if ok || got != 0 {
					t.Fatalf("MultiplicativeOrder(%d, %d) = (%d, %v), want (0, false)", a, m, got, ok)
				}
				continue
			}
			// brute-force smallest k>0 with a^k ≡ 1 (mod m)
			var want uint64
			cur := uint64(1) % m
			for k := uint64(1); k <= m; k++ {
				cur = (cur * a) % m
				if cur == 1 {
					want = k
					break
				}
			}
			if !ok || got != want {
				t.Fatalf("MultiplicativeOrder(%d, %d) = (%d, %v), want (%d, true)", a, m, got, ok, want)
			}
		}
	}
}

func TestPrimitiveRoot(t *testing.T) {
	cases := []struct {
		m      uint64
		want   uint64
		wantOK bool
	}{
		{1, 0, true},
		{2, 1, true},
		{3, 2, true},
		{4, 3, true},
		{5, 2, true},
		{6, 5, true},
		{7, 3, true},
		{9, 2, true},
		{10, 3, true},
		{11, 2, true},
		{13, 2, true},
		{14, 3, true},
		{8, 0, false},
		{12, 0, false},
		{15, 0, false},
		{16, 0, false},
		{24, 0, false},
	}
	for _, c := range cases {
		got, ok := PrimitiveRoot(c.m)
		if got != c.want || ok != c.wantOK {
			t.Errorf("PrimitiveRoot(%d) = (%d, %v), want (%d, %v)",
				c.m, got, ok, c.want, c.wantOK)
		}
	}
}

func TestIsPrimitiveRoot(t *testing.T) {
	cases := []struct {
		g, m uint64
		want bool
	}{
		{1, 1, true},
		{3, 7, true},
		{5, 7, true},
		{2, 7, false}, // order 3, not 6
		{2, 5, true},
		{3, 5, true},
		{2, 9, true},
		{3, 4, true},
		{2, 4, false}, // gcd(2,4) != 1
		{5, 6, true},
		{1, 6, false}, // order 1
		{2, 8, false}, // no primitive root modulo 8
	}
	for _, c := range cases {
		if got := IsPrimitiveRoot(c.g, c.m); got != c.want {
			t.Errorf("IsPrimitiveRoot(%d, %d) = %v, want %v", c.g, c.m, got, c.want)
		}
	}
}

func TestPrimitiveRoots(t *testing.T) {
	cases := []struct {
		m    uint64
		want []uint64
	}{
		{1, []uint64{0}},
		{2, []uint64{1}},
		{4, []uint64{3}},
		{5, []uint64{2, 3}},
		{6, []uint64{5}},
		{7, []uint64{3, 5}},
		{9, []uint64{2, 5}},
		{11, []uint64{2, 6, 7, 8}},
		{8, nil},
		{15, nil},
	}
	for _, c := range cases {
		got := PrimitiveRoots(c.m)
		if !reflect.DeepEqual(got, c.want) {
			t.Errorf("PrimitiveRoots(%d) = %v, want %v", c.m, got, c.want)
		}
	}
}

// TestPrimitiveRootsConsistency checks every reported root is primitive, the
// count equals φ(φ(m)), and the slice is ascending.
func TestPrimitiveRootsConsistency(t *testing.T) {
	for m := uint64(1); m <= 100; m++ {
		roots := PrimitiveRoots(m)
		if roots == nil {
			if ntheoryHasPrimitiveRootU64(m) {
				t.Fatalf("PrimitiveRoots(%d) = nil, but a primitive root exists", m)
			}
			continue
		}
		if want := ntheoryPhiU64(ntheoryPhiU64(m)); uint64(len(roots)) != want {
			t.Fatalf("PrimitiveRoots(%d) count = %d, want φ(φ(m)) = %d", m, len(roots), want)
		}
		for i, g := range roots {
			if i > 0 && roots[i-1] >= g {
				t.Fatalf("PrimitiveRoots(%d) not ascending: %v", m, roots)
			}
			if !IsPrimitiveRoot(g, m) {
				t.Fatalf("PrimitiveRoots(%d) reported %d which is not a primitive root", m, g)
			}
		}
	}
}

func BenchmarkCarmichael(b *testing.B) {
	for i := 0; i < b.N; i++ {
		_ = Carmichael(982451653 * 2)
	}
}

func BenchmarkMultiplicativeOrder(b *testing.B) {
	const m = uint64(1000000007) // prime
	for i := 0; i < b.N; i++ {
		_, _ = MultiplicativeOrder(3, m)
	}
}

func BenchmarkPrimitiveRoot(b *testing.B) {
	const m = uint64(1000000007) // prime
	for i := 0; i < b.N; i++ {
		_, _ = PrimitiveRoot(m)
	}
}
