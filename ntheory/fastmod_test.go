package ntheory

import (
	"math/big"
	"math/rand"
	"testing"
)

// bigMulMod computes (a*b) mod m with math/big as an independent oracle.
func bigMulMod(a, b, m uint64) uint64 {
	ba := new(big.Int).SetUint64(a)
	bb := new(big.Int).SetUint64(b)
	bm := new(big.Int).SetUint64(m)
	ba.Mul(ba, bb)
	ba.Mod(ba, bm)
	return ba.Uint64()
}

// bigModPow computes base**exp mod m with math/big as an independent oracle.
func bigModPow(base, exp, m uint64) uint64 {
	r := new(big.Int).Exp(
		new(big.Int).SetUint64(base),
		new(big.Int).SetUint64(exp),
		new(big.Int).SetUint64(m),
	)
	return r.Uint64()
}

func TestMulModU64Table(t *testing.T) {
	cases := []struct {
		a, b, m, want uint64
	}{
		{6, 7, 10, 2},
		{0, 123, 7, 0},
		{123, 0, 7, 0},
		{1, 1, 2, 1},
		{999, 999, 1000, 1},
		{1 << 63, 1 << 63, 1000000007, bigMulMod(1<<63, 1<<63, 1000000007)},
		{^uint64(0), ^uint64(0), 1000000007, bigMulMod(^uint64(0), ^uint64(0), 1000000007)},
		{^uint64(0) - 1, ^uint64(0) - 2, ^uint64(0), bigMulMod(^uint64(0)-1, ^uint64(0)-2, ^uint64(0))},
	}
	for _, c := range cases {
		if got := MulModU64(c.a, c.b, c.m); got != c.want {
			t.Errorf("MulModU64(%d,%d,%d) = %d, want %d", c.a, c.b, c.m, got, c.want)
		}
	}
}

func TestMulModU64Random(t *testing.T) {
	rng := rand.New(rand.NewSource(42))
	for i := 0; i < 20000; i++ {
		a := rng.Uint64()
		b := rng.Uint64()
		m := rng.Uint64()
		if m == 0 {
			m = 1
		}
		if got, want := MulModU64(a, b, m), bigMulMod(a, b, m); got != want {
			t.Fatalf("MulModU64(%d,%d,%d) = %d, want %d", a, b, m, got, want)
		}
	}
}

func TestMulModU64PanicsOnZero(t *testing.T) {
	defer func() {
		if recover() == nil {
			t.Fatal("MulModU64 with m==0 did not panic")
		}
	}()
	MulModU64(1, 2, 0)
}

func TestAddModU64(t *testing.T) {
	cases := []struct {
		a, b, m, want uint64
	}{
		{7, 8, 10, 5},
		{3, 4, 10, 7},
		{0, 0, 5, 0},
		{^uint64(0) - 2, ^uint64(0) - 2, ^uint64(0), func() uint64 {
			m := ^uint64(0)
			return (new(big.Int).Mod(
				new(big.Int).Add(new(big.Int).SetUint64(m-2), new(big.Int).SetUint64(m-2)),
				new(big.Int).SetUint64(m))).Uint64()
		}()},
	}
	for _, c := range cases {
		if got := AddModU64(c.a, c.b, c.m); got != c.want {
			t.Errorf("AddModU64(%d,%d,%d) = %d, want %d", c.a, c.b, c.m, got, c.want)
		}
	}
	// Cross-check against big over random inputs (exercises overflow paths).
	rng := rand.New(rand.NewSource(7))
	for i := 0; i < 20000; i++ {
		a, b := rng.Uint64(), rng.Uint64()
		m := rng.Uint64()
		if m == 0 {
			m = 1
		}
		want := new(big.Int).Mod(
			new(big.Int).Add(new(big.Int).SetUint64(a), new(big.Int).SetUint64(b)),
			new(big.Int).SetUint64(m)).Uint64()
		if got := AddModU64(a, b, m); got != want {
			t.Fatalf("AddModU64(%d,%d,%d) = %d, want %d", a, b, m, got, want)
		}
	}
}

func TestSubModU64(t *testing.T) {
	cases := []struct {
		a, b, m, want uint64
	}{
		{5, 3, 10, 2},
		{3, 5, 10, 8},
		{0, 1, 7, 6},
		{10, 10, 7, 0},
	}
	for _, c := range cases {
		if got := SubModU64(c.a, c.b, c.m); got != c.want {
			t.Errorf("SubModU64(%d,%d,%d) = %d, want %d", c.a, c.b, c.m, got, c.want)
		}
	}
	rng := rand.New(rand.NewSource(11))
	for i := 0; i < 20000; i++ {
		a, b := rng.Uint64(), rng.Uint64()
		m := rng.Uint64()
		if m == 0 {
			m = 1
		}
		want := new(big.Int).Mod(
			new(big.Int).Sub(new(big.Int).SetUint64(a), new(big.Int).SetUint64(b)),
			new(big.Int).SetUint64(m)).Uint64()
		if got := SubModU64(a, b, m); got != want {
			t.Fatalf("SubModU64(%d,%d,%d) = %d, want %d", a, b, m, got, want)
		}
	}
}

func TestModPowU64Table(t *testing.T) {
	cases := []struct {
		base, exp, m, want uint64
	}{
		{2, 10, 1000, 24},
		{3, 0, 7, 1},
		{0, 0, 7, 1},
		{5, 3, 13, 8},
		{7, 100, 1, 0},
		{2, 100, 13, bigModPow(2, 100, 13)},     // odd modulus -> Montgomery path
		{2, 100, 1024, bigModPow(2, 100, 1024)}, // even modulus -> MulModU64 path
	}
	for _, c := range cases {
		if got := ModPowU64(c.base, c.exp, c.m); got != c.want {
			t.Errorf("ModPowU64(%d,%d,%d) = %d, want %d", c.base, c.exp, c.m, got, c.want)
		}
	}
}

func TestModPowU64Random(t *testing.T) {
	rng := rand.New(rand.NewSource(99))
	for i := 0; i < 5000; i++ {
		base := rng.Uint64()
		exp := rng.Uint64()
		m := rng.Uint64()
		if m == 0 {
			m = 1
		}
		if got, want := ModPowU64(base, exp, m), bigModPow(base, exp, m); got != want {
			t.Fatalf("ModPowU64(%d,%d,%d) = %d, want %d", base, exp, m, got, want)
		}
	}
}

func TestMontgomeryRoundTrip(t *testing.T) {
	mods := []uint64{3, 5, 97, 1000000007, (1 << 61) - 1, ^uint64(0)}
	rng := rand.New(rand.NewSource(2024))
	for _, m := range mods {
		mo := NewMontgomery(m)
		if mo.Modulus() != m {
			t.Fatalf("Modulus() = %d, want %d", mo.Modulus(), m)
		}
		for i := 0; i < 2000; i++ {
			a := rng.Uint64() % m
			if got := mo.FromMont(mo.ToMont(a)); got != a {
				t.Fatalf("FromMont(ToMont(%d)) = %d for m=%d", a, got, m)
			}
		}
	}
}

func TestMontgomeryMulMont(t *testing.T) {
	mods := []uint64{97, 1000000007, (1 << 61) - 1, ^uint64(0)}
	rng := rand.New(rand.NewSource(31))
	for _, m := range mods {
		mo := NewMontgomery(m)
		for i := 0; i < 3000; i++ {
			a := rng.Uint64() % m
			b := rng.Uint64() % m
			// MulMont operates in the Montgomery domain; convert back to compare.
			got := mo.FromMont(mo.MulMont(mo.ToMont(a), mo.ToMont(b)))
			if want := bigMulMod(a, b, m); got != want {
				t.Fatalf("MulMont domain product for %d*%d mod %d = %d, want %d", a, b, m, got, want)
			}
		}
	}
}

func TestMontgomeryPowMont(t *testing.T) {
	mods := []uint64{3, 97, 1000000007, (1 << 61) - 1, ^uint64(0)}
	rng := rand.New(rand.NewSource(1234))
	for _, m := range mods {
		mo := NewMontgomery(m)
		for i := 0; i < 2000; i++ {
			base := rng.Uint64()
			exp := rng.Uint64()
			if got, want := mo.PowMont(base, exp), bigModPow(base%m, exp, m); got != want {
				t.Fatalf("PowMont(%d,%d) mod %d = %d, want %d", base, exp, m, got, want)
			}
		}
	}
}

func TestNewMontgomeryPanics(t *testing.T) {
	for _, m := range []uint64{0, 1, 2, 4, 100} {
		func(m uint64) {
			defer func() {
				if recover() == nil {
					t.Fatalf("NewMontgomery(%d) did not panic", m)
				}
			}()
			NewMontgomery(m)
		}(m)
	}
}

func TestBarrettReduce(t *testing.T) {
	cases := []struct {
		m, x, want uint64
	}{
		{1000, 123456, 456},
		{1, 999, 0},
		{2, 999, 1},
		{7, 0, 0},
	}
	for _, c := range cases {
		b := NewBarrett(c.m)
		if got := b.Reduce(c.x); got != c.want {
			t.Errorf("Barrett(%d).Reduce(%d) = %d, want %d", c.m, c.x, got, c.want)
		}
	}
	// Cross-check Reduce and MulMod against big, including even moduli and moduli
	// with the top bit set (shift == 0).
	mods := []uint64{2, 6, 1000, 1 << 32, (1 << 62), (1 << 63) + 1, ^uint64(0)}
	rng := rand.New(rand.NewSource(555))
	for _, m := range mods {
		b := NewBarrett(m)
		if b.Modulus() != m {
			t.Fatalf("Modulus() = %d, want %d", b.Modulus(), m)
		}
		for i := 0; i < 5000; i++ {
			x := rng.Uint64()
			if got, want := b.Reduce(x), new(big.Int).Mod(new(big.Int).SetUint64(x), new(big.Int).SetUint64(m)).Uint64(); got != want {
				t.Fatalf("Barrett(%d).Reduce(%d) = %d, want %d", m, x, got, want)
			}
			a, c := rng.Uint64(), rng.Uint64()
			if got, want := b.MulMod(a, c), bigMulMod(a, c, m); got != want {
				t.Fatalf("Barrett(%d).MulMod(%d,%d) = %d, want %d", m, a, c, got, want)
			}
		}
	}
}

func TestNewBarrettPanicsOnZero(t *testing.T) {
	defer func() {
		if recover() == nil {
			t.Fatal("NewBarrett(0) did not panic")
		}
	}()
	NewBarrett(0)
}

var sinkU64 uint64

func BenchmarkMulModU64(b *testing.B) {
	const m = 1000000007
	x := uint64(123456789)
	for i := 0; i < b.N; i++ {
		x = MulModU64(x, 2862933555777941757, m)
	}
	sinkU64 = x
}

func BenchmarkModPowU64Odd(b *testing.B) {
	const m = 1000000007
	var r uint64
	for i := 0; i < b.N; i++ {
		r = ModPowU64(2862933555777941757, uint64(i)|1, m)
	}
	sinkU64 = r
}

func BenchmarkModPowU64Even(b *testing.B) {
	const m = 1 << 40
	var r uint64
	for i := 0; i < b.N; i++ {
		r = ModPowU64(2862933555777941757, uint64(i)|1, m)
	}
	sinkU64 = r
}

func BenchmarkMontgomeryPowMont(b *testing.B) {
	mo := NewMontgomery(1000000007)
	var r uint64
	for i := 0; i < b.N; i++ {
		r = mo.PowMont(2862933555777941757, uint64(i)|1)
	}
	sinkU64 = r
}

func BenchmarkBarrettMulMod(b *testing.B) {
	ba := NewBarrett(1 << 40)
	x := uint64(123456789)
	for i := 0; i < b.N; i++ {
		x = ba.MulMod(x, 2862933555777941757)
	}
	sinkU64 = x
}
