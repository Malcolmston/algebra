package discretemath

import (
	"math/bits"
	"strconv"
	"testing"
)

// discretemathSampleWords is a deterministic spread of test words covering
// edges, single bits, alternating patterns and mixed values.
var discretemathSampleWords = []uint64{
	0, 1, 2, 3, 0xff, 0x100, 0xdeadbeef, 0xcafebabe,
	0x5555555555555555, 0xaaaaaaaaaaaaaaaa, 0x0f0f0f0f0f0f0f0f,
	0x8000000000000000, 0xffffffffffffffff, 0x123456789abcdef0,
	1 << 17, 1 << 40, 1<<63 | 1, 0x7fffffffffffffff,
}

func TestPopCountAgainstStdlib(t *testing.T) {
	for _, x := range discretemathSampleWords {
		if got, want := PopCount(x), bits.OnesCount64(x); got != want {
			t.Errorf("PopCount(%#x) = %d, want %d", x, got, want)
		}
	}
}

func TestParity(t *testing.T) {
	for _, x := range discretemathSampleWords {
		if got, want := Parity(x), bits.OnesCount64(x)&1; got != want {
			t.Errorf("Parity(%#x) = %d, want %d", x, got, want)
		}
	}
}

func TestHammingDistance(t *testing.T) {
	cases := []struct {
		a, b uint64
		want int
	}{
		{0, 0, 0},
		{0, 0xffffffffffffffff, 64},
		{0b1010, 0b0101, 4},
		{0x0f, 0x00, 4},
		{0xdeadbeef, 0xdeadbeef, 0},
	}
	for _, c := range cases {
		if got := HammingDistance(c.a, c.b); got != c.want {
			t.Errorf("HammingDistance(%#x, %#x) = %d, want %d", c.a, c.b, got, c.want)
		}
	}
}

func TestHammingDistanceBytes(t *testing.T) {
	a := []byte{0x00, 0xff, 0x0f}
	b := []byte{0xff, 0xff, 0xf0}
	got, err := HammingDistanceBytes(a, b)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if got != 16 {
		t.Errorf("HammingDistanceBytes = %d, want 16", got)
	}
	if _, err := HammingDistanceBytes([]byte{1}, []byte{1, 2}); err == nil {
		t.Error("expected length-mismatch error")
	}
}

func TestGrayRoundTrip(t *testing.T) {
	// Successive Gray codes must differ by exactly one bit.
	for i := uint64(0); i < 1024; i++ {
		g := GrayEncode(i)
		if got := GrayDecode(g); got != i {
			t.Fatalf("GrayDecode(GrayEncode(%d)) = %d", i, got)
		}
		if i > 0 {
			if d := HammingDistance(GrayEncode(i-1), g); d != 1 {
				t.Fatalf("Gray codes %d and %d differ in %d bits, want 1", i-1, i, d)
			}
		}
	}
	// Known small Gray codes.
	known := []uint64{0, 1, 3, 2, 6, 7, 5, 4}
	for i, want := range known {
		if got := GrayEncode(uint64(i)); got != want {
			t.Errorf("GrayEncode(%d) = %d, want %d", i, got, want)
		}
	}
}

func TestMorton2DRoundTrip(t *testing.T) {
	coords := []uint32{0, 1, 2, 3, 0xff, 0xffff, 0x12345678, 0xffffffff, 0xdeadbeef}
	for _, x := range coords {
		for _, y := range coords {
			m := MortonEncode2D(x, y)
			gx, gy := MortonDecode2D(m)
			if gx != x || gy != y {
				t.Fatalf("Morton2D round trip (%#x,%#x) -> %#x -> (%#x,%#x)", x, y, m, gx, gy)
			}
		}
	}
	// Known interleave: x=0b11 (bits 0,1), y=0 -> 0b0101 = 5.
	if got := MortonEncode2D(0b11, 0); got != 0b0101 {
		t.Errorf("MortonEncode2D(3,0) = %#b, want 0101", got)
	}
	// x=0, y=0b11 -> 0b1010 = 10.
	if got := MortonEncode2D(0, 0b11); got != 0b1010 {
		t.Errorf("MortonEncode2D(0,3) = %#b, want 1010", got)
	}
}

func TestMorton3DRoundTrip(t *testing.T) {
	coords := []uint32{0, 1, 2, 7, 0xff, 0x1fffff, 0x155555, 0xabcde}
	for _, x := range coords {
		for _, y := range coords {
			for _, z := range coords {
				m := MortonEncode3D(x, y, z)
				gx, gy, gz := MortonDecode3D(m)
				if gx != x || gy != y || gz != z {
					t.Fatalf("Morton3D round trip (%#x,%#x,%#x) -> (%#x,%#x,%#x)", x, y, z, gx, gy, gz)
				}
			}
		}
	}
}

func TestReverseBits(t *testing.T) {
	for _, x := range discretemathSampleWords {
		if got, want := ReverseBits(x), bits.Reverse64(x); got != want {
			t.Errorf("ReverseBits(%#x) = %#x, want %#x", x, got, want)
		}
	}
	for _, x := range []uint32{0, 1, 0x80000000, 0xdeadbeef, 0xffffffff} {
		if got, want := ReverseBits32(x), bits.Reverse32(x); got != want {
			t.Errorf("ReverseBits32(%#x) = %#x, want %#x", x, got, want)
		}
	}
}

func TestRotate(t *testing.T) {
	for _, x := range discretemathSampleWords {
		for _, k := range []int{0, 1, 7, 31, 63, 64, 65, -1, -63, -64} {
			if got, want := RotateLeft(x, k), bits.RotateLeft64(x, k); got != want {
				t.Errorf("RotateLeft(%#x, %d) = %#x, want %#x", x, k, got, want)
			}
			if got, want := RotateRight(x, k), bits.RotateLeft64(x, -k); got != want {
				t.Errorf("RotateRight(%#x, %d) = %#x, want %#x", x, k, got, want)
			}
		}
	}
}

func TestZeroCounts(t *testing.T) {
	for _, x := range discretemathSampleWords {
		if got, want := TrailingZeros(x), bits.TrailingZeros64(x); got != want {
			t.Errorf("TrailingZeros(%#x) = %d, want %d", x, got, want)
		}
		if got, want := LeadingZeros(x), bits.LeadingZeros64(x); got != want {
			t.Errorf("LeadingZeros(%#x) = %d, want %d", x, got, want)
		}
	}
}

func TestLog2AndPowers(t *testing.T) {
	for _, c := range []struct {
		x           uint64
		floor, ceil int
		prev, next  uint64
		isPow2      bool
	}{
		{1, 0, 0, 1, 1, true},
		{2, 1, 1, 2, 2, true},
		{3, 1, 2, 2, 4, false},
		{4, 2, 2, 4, 4, true},
		{7, 2, 3, 4, 8, false},
		{8, 3, 3, 8, 8, true},
		{1000, 9, 10, 512, 1024, false},
		{1 << 40, 40, 40, 1 << 40, 1 << 40, true},
	} {
		if got := Log2Floor(c.x); got != c.floor {
			t.Errorf("Log2Floor(%d) = %d, want %d", c.x, got, c.floor)
		}
		if got := Log2Ceil(c.x); got != c.ceil {
			t.Errorf("Log2Ceil(%d) = %d, want %d", c.x, got, c.ceil)
		}
		if got := PrevPowerOfTwo(c.x); got != c.prev {
			t.Errorf("PrevPowerOfTwo(%d) = %d, want %d", c.x, got, c.prev)
		}
		if got := NextPowerOfTwo(c.x); got != c.next {
			t.Errorf("NextPowerOfTwo(%d) = %d, want %d", c.x, got, c.next)
		}
		if got := IsPowerOfTwo(c.x); got != c.isPow2 {
			t.Errorf("IsPowerOfTwo(%d) = %v, want %v", c.x, got, c.isPow2)
		}
	}
	if Log2Floor(0) != -1 || Log2Ceil(0) != -1 {
		t.Error("Log2 of zero should be -1")
	}
	if IsPowerOfTwo(0) {
		t.Error("0 is not a power of two")
	}
}

func TestSingleBitOps(t *testing.T) {
	var x uint64
	x = SetBit(x, 5)
	if x != 32 || !TestBit(x, 5) {
		t.Errorf("SetBit/TestBit failed: %d", x)
	}
	x = ToggleBit(x, 5)
	if x != 0 || TestBit(x, 5) {
		t.Errorf("ToggleBit failed: %d", x)
	}
	x = SetBit(SetBit(0, 0), 3) // 0b1001 = 9
	if x != 9 {
		t.Errorf("expected 9, got %d", x)
	}
	x = ClearBit(x, 0)
	if x != 8 {
		t.Errorf("ClearBit failed: %d", x)
	}
	if LowestSetBit(0b101100) != 0b100 {
		t.Error("LowestSetBit failed")
	}
	if ClearLowestSetBit(0b101100) != 0b101000 {
		t.Error("ClearLowestSetBit failed")
	}
}

func TestBitString(t *testing.T) {
	cases := []struct {
		x     uint64
		width int
		want  string
	}{
		{0, 0, "0"},
		{5, 0, "101"},
		{5, 8, "00000101"},
		{255, 8, "11111111"},
		{1, 4, "0001"},
	}
	for _, c := range cases {
		if got := BitString(c.x, c.width); got != c.want {
			t.Errorf("BitString(%d, %d) = %q, want %q", c.x, c.width, got, c.want)
		}
		// Cross-check against strconv for the natural width.
		if c.width == 0 {
			if got := BitString(c.x, 0); got != strconv.FormatUint(c.x, 2) {
				t.Errorf("BitString(%d,0) = %q, want %q", c.x, got, strconv.FormatUint(c.x, 2))
			}
		}
	}
}

// BenchmarkLevenshteinDistance benchmarks the heaviest routine in the package
// (a full O(n*m) dynamic program) on a pair of moderately long strings.
func BenchmarkLevenshteinDistance(b *testing.B) {
	s1 := "the quick brown fox jumps over the lazy dog near the riverbank"
	s2 := "a quick browne foxx jump over some lazy dogs by the river bank!!"
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = LevenshteinDistance(s1, s2)
	}
}
