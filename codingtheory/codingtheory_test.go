package codingtheory

import (
	"fmt"
	"math"
	"math/rand"
	"reflect"
	"testing"
)

func eqInts(a, b []int) bool { return reflect.DeepEqual(a, b) }

func approx(a, b, tol float64) bool { return math.Abs(a-b) <= tol }

// --- GF(2) polynomial arithmetic ---

func TestGF2Poly(t *testing.T) {
	if got := GF2PolyMul(0b110, 0b11); got != 0b1010 {
		t.Errorf("GF2PolyMul: got %b want 1010", got)
	}
	if got := GF2PolyDegree(0b10011); got != 4 {
		t.Errorf("GF2PolyDegree: got %d want 4", got)
	}
	q, r := GF2PolyDivMod(0b101011, 0b11)
	// x^5+x^3+x+1 divided by x+1
	if GF2PolyMul(q, 0b11)^r != 0b101011 {
		t.Errorf("GF2PolyDivMod reconstruction failed: q=%b r=%b", q, r)
	}
	if got := GF2PolyGCD(0b1010, 0b110); got == 0 {
		t.Errorf("GF2PolyGCD returned 0")
	}
	if got := GF2PolyString(0b10011); got != "x^4 + x + 1" {
		t.Errorf("GF2PolyString: got %q", got)
	}
}

// --- GF(2^m) field ---

func TestFieldGF8(t *testing.T) {
	f := NewGF8()
	if f.Size() != 8 || f.Order() != 7 {
		t.Fatalf("size/order wrong")
	}
	// alpha^3 = alpha + 1 for x^3+x+1
	if got := f.Exp(3); got != 0b011 {
		t.Errorf("alpha^3: got %b want 011", got)
	}
	// multiplication/inverse consistency
	for a := 1; a < 8; a++ {
		if f.Mul(a, f.Inv(a)) != 1 {
			t.Errorf("Inv failed for %d", a)
		}
		for b := 1; b < 8; b++ {
			if f.Div(f.Mul(a, b), b) != a {
				t.Errorf("Div/Mul inconsistent a=%d b=%d", a, b)
			}
		}
	}
	// log/exp inverse
	for a := 1; a < 8; a++ {
		if f.Exp(f.Log(a)) != a {
			t.Errorf("Exp(Log) failed for %d", a)
		}
	}
}

func TestFieldGF256(t *testing.T) {
	f := NewGF256()
	if f.Poly() != 0x11d {
		t.Fatalf("poly wrong")
	}
	// exhaustive mul/div check on a sample
	for _, a := range []int{1, 2, 7, 100, 200, 255} {
		if f.Mul(a, f.Inv(a)) != 1 {
			t.Errorf("Inv failed for %d", a)
		}
		if f.Pow(a, 255) != 1 {
			t.Errorf("a^255 != 1 for %d", a)
		}
	}
	// trace is in GF(2)
	for a := 0; a < 256; a++ {
		if tr := f.Trace(a); tr != 0 && tr != 1 {
			t.Fatalf("Trace not in GF(2) for %d: %d", a, tr)
		}
	}
}

func TestNonPrimitiveRejected(t *testing.T) {
	// x^4+x^3+x^2+x+1 divides x^5-1 and is irreducible but NOT primitive.
	if _, err := NewField(4, 0b11111); err != ErrNotPrimitive {
		t.Errorf("expected ErrNotPrimitive, got %v", err)
	}
}

func TestMinimalPoly(t *testing.T) {
	f := NewGF16()
	// minimal poly of alpha in GF(16) with x^4+x+1 is x^4+x+1 = 0b10011
	if got := f.MinimalPoly(1); got != 0b10011 {
		t.Errorf("MinimalPoly(alpha): got %b want 10011", got)
	}
}

// --- bit helpers ---

func TestBitHelpers(t *testing.T) {
	if HammingWeight([]int{1, 0, 1, 1}) != 3 {
		t.Error("HammingWeight")
	}
	if HammingDistance([]int{1, 0, 1}, []int{0, 0, 1}) != 1 {
		t.Error("HammingDistance")
	}
	if HammingWeightUint(0b1011) != 3 {
		t.Error("HammingWeightUint")
	}
	if BitsToUint([]int{1, 0, 1}) != 5 {
		t.Error("BitsToUint")
	}
	if !eqInts(UintToBits(5, 4), []int{0, 1, 0, 1}) {
		t.Error("UintToBits")
	}
	if RankGF2([][]int{{1, 0, 1}, {0, 1, 1}, {1, 1, 0}}) != 2 {
		t.Error("RankGF2")
	}
}

// --- Hamming(7,4) ---

func TestHamming74(t *testing.T) {
	for m := 0; m < 16; m++ {
		data := UintToBits(uint64(m), 4)
		cw, err := Hamming74Encode(data)
		if err != nil {
			t.Fatal(err)
		}
		// zero syndrome for clean codeword
		s, _ := Hamming74Syndrome(cw)
		if HammingWeight(s) != 0 {
			t.Fatalf("clean codeword has nonzero syndrome for m=%d", m)
		}
		// single error in every position is corrected
		for pos := 0; pos < 7; pos++ {
			recv := append([]int(nil), cw...)
			recv[pos] ^= 1
			dec, corrected, ep, err := Hamming74Decode(recv)
			if err != nil {
				t.Fatal(err)
			}
			if ep != pos {
				t.Errorf("m=%d pos=%d: reported error position %d", m, pos, ep)
			}
			if !eqInts(dec, data) {
				t.Errorf("m=%d pos=%d: data mismatch", m, pos)
			}
			if !eqInts(corrected, cw) {
				t.Errorf("m=%d pos=%d: codeword not corrected", m, pos)
			}
		}
	}
}

func TestExtendedHamming74DoubleError(t *testing.T) {
	data := []int{1, 0, 1, 1}
	cw, _ := ExtendedHamming74Encode(data)
	// double error should be flagged
	cw[0] ^= 1
	cw[1] ^= 1
	_, dbl, err := ExtendedHamming74Decode(cw)
	if err != nil {
		t.Fatal(err)
	}
	if !dbl {
		t.Error("expected double error detection")
	}
}

// --- general Hamming code ---

func TestHammingCodeGeneral(t *testing.T) {
	c, err := NewHammingCode(4) // [15,11]
	if err != nil {
		t.Fatal(err)
	}
	if c.N() != 15 || c.K() != 11 {
		t.Fatalf("dimensions wrong: n=%d k=%d", c.N(), c.K())
	}
	lc, err := c.LinearCode()
	if err != nil {
		t.Fatal(err)
	}
	msg := make([]int, 11)
	for i := range msg {
		msg[i] = (i * 7) % 2
	}
	cw, err := lc.Encode(msg)
	if err != nil {
		t.Fatal(err)
	}
	// inject single error at each position, verify syndrome position
	for pos := 0; pos < 15; pos++ {
		recv := append([]int(nil), cw...)
		recv[pos] ^= 1
		corrected, ep, _ := c.Decode(recv)
		if ep != pos {
			t.Errorf("pos %d: reported %d", pos, ep)
		}
		if !eqInts(corrected, cw) {
			t.Errorf("pos %d not corrected", pos)
		}
	}
}

// --- linear code, syndrome decoding ---

func TestLinearCodeSyndromeDecode(t *testing.T) {
	c, err := NewLinearCode(Hamming74Generator())
	if err != nil {
		t.Fatal(err)
	}
	if c.MinimumDistance() != 3 {
		t.Errorf("Hamming(7,4) min distance: got %d want 3", c.MinimumDistance())
	}
	table := c.SyndromeTable()
	msg := []int{1, 1, 0, 1}
	cw, _ := c.Encode(msg)
	recv := append([]int(nil), cw...)
	recv[2] ^= 1
	dec, _ := c.DecodeSyndrome(recv, table)
	if !eqInts(dec, cw) {
		t.Errorf("syndrome decode failed: got %v want %v", dec, cw)
	}
}

// --- cyclic code ---

func TestCyclicCode(t *testing.T) {
	// Hamming(7,4) as cyclic code with g(x)=x^3+x+1 = 0b1011
	c, err := NewCyclicCode(7, 0b1011)
	if err != nil {
		t.Fatal(err)
	}
	if c.K() != 4 {
		t.Fatalf("k wrong: %d", c.K())
	}
	msg := []int{1, 0, 1, 1}
	cw, _ := c.Encode(msg)
	if !c.IsCodeword(cw) {
		t.Error("encoded word is not a codeword")
	}
	// non-divisor generator rejected
	if _, err := NewCyclicCode(7, 0b101); err == nil {
		t.Error("expected error for non-dividing generator")
	}
}

// --- Golay ---

func TestGolay(t *testing.T) {
	g := NewGolayCode()
	rng := rand.New(rand.NewSource(1))
	for trial := 0; trial < 500; trial++ {
		msg := make([]int, 12)
		for i := range msg {
			msg[i] = rng.Intn(2)
		}
		cw, _ := g.Encode(msg)
		if !g.Cyclic().IsCodeword(cw) {
			t.Fatal("golay encode not codeword")
		}
		ne := rng.Intn(4) // up to 3 errors
		recv := append([]int(nil), cw...)
		perm := rng.Perm(23)
		for i := 0; i < ne; i++ {
			recv[perm[i]] ^= 1
		}
		dec, n, err := g.DecodeMessage(recv)
		if err != nil {
			t.Fatal(err)
		}
		if !eqInts(dec, msg) {
			t.Fatalf("golay decode mismatch with %d errors", ne)
		}
		if n != ne {
			t.Fatalf("golay error count: got %d want %d", n, ne)
		}
	}
}

// --- Reed-Solomon ---

func TestReedSolomon(t *testing.T) {
	f := NewGF256()
	for _, fcr := range []int{0, 1} {
		rs, err := NewReedSolomonFCR(f, 8, fcr) // t=4
		if err != nil {
			t.Fatal(err)
		}
		rng := rand.New(rand.NewSource(int64(fcr + 1)))
		for trial := 0; trial < 400; trial++ {
			msg := make([]int, rs.K())
			for i := range msg {
				msg[i] = rng.Intn(256)
			}
			cw, _ := rs.Encode(msg)
			if !rs.IsCodeword(cw) {
				t.Fatal("rs encode not codeword")
			}
			ne := rng.Intn(rs.Correction() + 1)
			recv := append([]int(nil), cw...)
			perm := rng.Perm(rs.N())
			for i := 0; i < ne; i++ {
				recv[perm[i]] ^= 1 + rng.Intn(255)
			}
			dec, n, err := rs.DecodeMessage(recv)
			if err != nil {
				t.Fatalf("fcr=%d decode error with %d errors: %v", fcr, ne, err)
			}
			if !eqInts(dec, msg) {
				t.Fatalf("fcr=%d rs decode mismatch with %d errors", fcr, ne)
			}
			_ = n
		}
	}
}

// --- BCH ---

func TestBCHGenerators(t *testing.T) {
	f := NewGF16()
	tests := []struct {
		t       int
		wantGen int
		wantK   int
	}{
		{1, 0b10011, 11},
		{2, 0b111010001, 7},
		{3, 0b10100110111, 5},
	}
	for _, tc := range tests {
		bch, err := NewBCHCode(f, tc.t)
		if err != nil {
			t.Fatal(err)
		}
		if bch.Generator() != tc.wantGen {
			t.Errorf("t=%d generator: got %b want %b", tc.t, bch.Generator(), tc.wantGen)
		}
		if bch.K() != tc.wantK {
			t.Errorf("t=%d k: got %d want %d", tc.t, bch.K(), tc.wantK)
		}
	}
}

func TestBCHDecode(t *testing.T) {
	f := NewGF16()
	bch, _ := NewBCHCode(f, 3) // corrects 3 errors
	rng := rand.New(rand.NewSource(7))
	for trial := 0; trial < 500; trial++ {
		msg := make([]int, bch.K())
		for i := range msg {
			msg[i] = rng.Intn(2)
		}
		cw, _ := bch.Encode(msg)
		ne := rng.Intn(4)
		recv := append([]int(nil), cw...)
		perm := rng.Perm(bch.N())
		for i := 0; i < ne; i++ {
			recv[perm[i]] ^= 1
		}
		dec, _, err := bch.DecodeMessage(recv)
		if err != nil {
			t.Fatalf("bch decode error with %d errors: %v", ne, err)
		}
		if !eqInts(dec, msg) {
			t.Fatalf("bch decode mismatch with %d errors", ne)
		}
	}
}

// --- CRC ---

func TestCRCCheckValues(t *testing.T) {
	data := []byte("123456789")
	tests := []struct {
		name string
		crc  CRC
		want uint64
	}{
		{"CRC8", CRC8(), 0xF4},
		{"CRC16CCITTFalse", CRC16CCITTFalse(), 0x29B1},
		{"CRC16XModem", CRC16XModem(), 0x31C3},
		{"CRC16IBM", CRC16IBM(), 0xBB3D},
		{"CRC32", CRC32(), 0xCBF43926},
		{"CRC32C", CRC32C(), 0xE3069283},
	}
	for _, tc := range tests {
		if got := tc.crc.Checksum(data); got != tc.want {
			t.Errorf("%s: got %X want %X", tc.name, got, tc.want)
		}
	}
}

func TestCRCAppendVerify(t *testing.T) {
	c := CRC32()
	data := []byte("hello world")
	framed := c.Append(data)
	// recomputing the CRC over data must match the appended value
	sum := c.Checksum(data)
	if !c.Verify(data, sum) {
		t.Error("verify failed")
	}
	if len(framed) != len(data)+4 {
		t.Error("append length wrong")
	}
}

// --- simple codes ---

func TestRepetition(t *testing.T) {
	msg := []int{1, 0, 1}
	enc, _ := RepetitionEncode(msg, 3)
	if len(enc) != 9 {
		t.Fatal("repetition length")
	}
	enc[0] ^= 1 // one error in first block, majority still recovers
	dec, _ := RepetitionDecode(enc, 3)
	if !eqInts(dec, msg) {
		t.Errorf("repetition decode: got %v want %v", dec, msg)
	}
}

func TestParity(t *testing.T) {
	enc, _ := SingleParityEncode([]int{1, 0, 1})
	if !SingleParityCheck(enc) {
		t.Error("parity check should pass on clean word")
	}
	enc[0] ^= 1
	_, detected, _ := SingleParityDecode(enc)
	if !detected {
		t.Error("single error should be detected")
	}
}

// --- Hadamard / Walsh ---

func TestWalshHadamardTransform(t *testing.T) {
	// forward then inverse is identity
	a := []float64{1, 0, 1, 0, 0, 1, 1, 0}
	back := InverseWalshHadamardTransform(WalshHadamardTransform(a))
	for i := range a {
		if !approx(a[i], back[i], 1e-9) {
			t.Fatalf("WHT inverse mismatch at %d", i)
		}
	}
}

func TestHadamardCode(t *testing.T) {
	// [16,4,8] Hadamard code corrects up to 3 errors
	rng := rand.New(rand.NewSource(3))
	for trial := 0; trial < 200; trial++ {
		msg := make([]int, 4)
		for i := range msg {
			msg[i] = rng.Intn(2)
		}
		cw, _ := HadamardEncode(msg)
		recv := append([]int(nil), cw...)
		ne := rng.Intn(4)
		perm := rng.Perm(len(recv))
		for i := 0; i < ne; i++ {
			recv[perm[i]] ^= 1
		}
		dec, _ := HadamardDecode(recv)
		if !eqInts(dec, msg) {
			t.Fatalf("hadamard decode mismatch with %d errors: got %v want %v", ne, dec, msg)
		}
	}
}

func TestWalshOrthogonality(t *testing.T) {
	// distinct Walsh codes are orthogonal (equal number of agreements/disagreements)
	a := WalshCode(3, 2)
	b := WalshCode(3, 5)
	agree := 0
	for i := range a {
		if a[i] == b[i] {
			agree++
		}
	}
	if agree != len(a)/2 {
		t.Errorf("Walsh codes not orthogonal: agree=%d", agree)
	}
}

// --- convolutional ---

func TestConvolutionalEncode(t *testing.T) {
	c, _ := NewConvCodeOctal(3, []int{7, 5})
	cw, _ := c.Encode([]int{1, 0, 1, 1})
	want := []int{1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1}
	if !eqInts(cw, want) {
		t.Errorf("conv encode: got %v want %v", cw, want)
	}
}

func TestViterbi(t *testing.T) {
	c, _ := NewConvCodeOctal(3, []int{7, 5})
	rng := rand.New(rand.NewSource(11))
	for trial := 0; trial < 500; trial++ {
		L := 6 + rng.Intn(10)
		msg := make([]int, L)
		for i := range msg {
			msg[i] = rng.Intn(2)
		}
		enc, _ := c.Encode(msg)
		recv := append([]int(nil), enc...)
		// flip up to 2 bits (within free distance 5, correctable)
		ne := rng.Intn(3)
		perm := rng.Perm(len(recv))
		for i := 0; i < ne; i++ {
			recv[perm[i]] ^= 1
		}
		dec, err := c.HardViterbiDecode(recv)
		if err != nil {
			t.Fatal(err)
		}
		if !eqInts(dec, msg) {
			t.Fatalf("viterbi mismatch with %d errors", ne)
		}
	}
}

// --- Gray codes ---

func TestGrayCode(t *testing.T) {
	for n := uint64(0); n < 256; n++ {
		g := BinaryToGray(n)
		if GrayToBinary(g) != n {
			t.Fatalf("gray roundtrip failed for %d", n)
		}
		if n > 0 {
			// consecutive gray codes differ in exactly one bit
			if HammingDistanceUint(BinaryToGray(n-1), g) != 1 {
				t.Fatalf("gray adjacency failed at %d", n)
			}
		}
	}
}

func TestNaryGray(t *testing.T) {
	radix := 4
	prev := []int(nil)
	for v := 0; v < 64; v++ {
		digits := []int{(v / 16) % 4, (v / 4) % 4, v % 4}
		g := NaryToGray(digits, radix)
		if !eqInts(GrayToNary(g, radix), digits) {
			t.Fatalf("nary gray roundtrip failed at %d", v)
		}
		if prev != nil {
			diff := 0
			for i := range g {
				if g[i] != prev[i] {
					diff++
				}
			}
			if diff != 1 {
				t.Fatalf("nary gray adjacency failed at %d", v)
			}
		}
		prev = g
	}
}

// --- interleaver ---

func TestInterleaver(t *testing.T) {
	data := make([]int, 12)
	for i := range data {
		data[i] = i
	}
	il, _ := BlockInterleave(data, 3, 4)
	back, _ := BlockDeinterleave(il, 3, 4)
	if !eqInts(back, data) {
		t.Error("block interleaver roundtrip failed")
	}
	perm := InterleavePermutation(3, 4)
	if inv := InvertPermutation(perm); len(inv) != 12 {
		t.Error("permutation inverse wrong length")
	}
}

// --- LDPC ---

func TestLDPCBitFlip(t *testing.T) {
	h, err := GallagerLDPC(12, 3, 4, 42)
	if err != nil {
		t.Fatal(err)
	}
	code, err := NewLDPCCode(h)
	if err != nil {
		t.Fatal(err)
	}
	if !code.IsRegular() {
		t.Error("Gallager code should be regular")
	}
	// all-zero word is a codeword; a single error should be correctable
	zero := make([]int, code.N())
	corrected := 0
	for pos := 0; pos < code.N(); pos++ {
		recv := append([]int(nil), zero...)
		recv[pos] ^= 1
		dec, _, ok, _ := code.BitFlipDecode(recv, 20)
		if ok && eqInts(dec, zero) {
			corrected++
		}
	}
	if corrected == 0 {
		t.Error("bit-flip decoder corrected no single errors")
	}
}

// --- Huffman ---

func TestHuffman(t *testing.T) {
	weights := map[int]int{0: 5, 1: 9, 2: 12, 3: 13, 4: 16, 5: 45}
	root, err := BuildHuffmanTree(weights)
	if err != nil {
		t.Fatal(err)
	}
	codes := HuffmanCodes(root)
	// highest-frequency symbol gets the shortest code
	if len(codes[5]) != 1 {
		t.Errorf("expected 1-bit code for most frequent symbol, got %d bits", len(codes[5]))
	}
	data := []int{5, 5, 2, 3, 0, 1, 4, 5, 2, 1, 0}
	bits, _ := HuffmanEncode(data, codes)
	dec, err := HuffmanDecode(bits, root)
	if err != nil {
		t.Fatal(err)
	}
	if !eqInts(dec, data) {
		t.Error("huffman roundtrip failed")
	}
	if avg := HuffmanAverageLength(codes, weights); !approx(avg, 2.24, 1e-9) {
		t.Errorf("average length: got %v want 2.24", avg)
	}
}

// --- arithmetic coding ---

func TestArithmeticCoding(t *testing.T) {
	freqs := []int{3, 3, 2, 1}
	rng := rand.New(rand.NewSource(99))
	for trial := 0; trial < 1000; trial++ {
		L := 1 + rng.Intn(40)
		msg := make([]int, L)
		for i := range msg {
			msg[i] = rng.Intn(len(freqs))
		}
		enc, err := ArithmeticEncode(msg, freqs)
		if err != nil {
			t.Fatal(err)
		}
		dec, err := ArithmeticDecode(enc, freqs, L)
		if err != nil {
			t.Fatal(err)
		}
		if !eqInts(dec, msg) {
			t.Fatalf("arithmetic roundtrip failed: %v -> %v", msg, dec)
		}
	}
	if !approx(ShannonCodeLength(0.25), 2, 1e-12) {
		t.Error("ShannonCodeLength(0.25) should be 2")
	}
}

// --- runnable example ---

func ExampleReedSolomon() {
	f := NewGF256()
	rs, _ := NewReedSolomon(f, 4) // 4 parity symbols, corrects 2 errors
	msg := make([]int, rs.K())
	for i := range msg {
		msg[i] = i % 7
	}
	code, _ := rs.Encode(msg)
	// corrupt two symbols
	code[3] ^= 0x55
	code[10] ^= 0xAA
	decoded, nerr, _ := rs.DecodeMessage(code)
	fmt.Println(nerr, eqInts(decoded, msg))
	// Output: 2 true
}

func ExampleHamming74Decode() {
	data := []int{1, 0, 1, 1}
	cw, _ := Hamming74Encode(data)
	cw[2] ^= 1 // flip one bit
	recovered, _, pos, _ := Hamming74Decode(cw)
	fmt.Println(pos, recovered)
	// Output: 2 [1 0 1 1]
}
