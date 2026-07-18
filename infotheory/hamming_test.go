package infotheory

import "testing"

func TestHamming74RoundTrip(t *testing.T) {
	for nibble := byte(0); nibble < 16; nibble++ {
		code := Hamming74Encode(nibble)
		data, pos, corrected := Hamming74Decode(code)
		if corrected || pos != 0 {
			t.Errorf("clean codeword for %d reported correction (pos=%d)", nibble, pos)
		}
		if data != nibble {
			t.Errorf("round trip: encoded %d decoded %d", nibble, data)
		}
	}
}

func TestHamming74SingleErrorCorrection(t *testing.T) {
	for nibble := byte(0); nibble < 16; nibble++ {
		code := Hamming74Encode(nibble)
		for bit := 0; bit < 7; bit++ {
			corrupted := code ^ (1 << bit)
			data, pos, corrected := Hamming74Decode(corrupted)
			if !corrected {
				t.Errorf("nibble %d bit %d: error not corrected", nibble, bit)
			}
			if pos != bit+1 {
				t.Errorf("nibble %d bit %d: reported pos %d, want %d", nibble, bit, pos, bit+1)
			}
			if data != nibble {
				t.Errorf("nibble %d bit %d: recovered %d, want %d", nibble, bit, data, nibble)
			}
		}
	}
}

func TestHamming74MinimumDistance(t *testing.T) {
	// Every pair of distinct Hamming(7,4) codewords differs in >= 3 bits.
	var codewords [16]byte
	for n := byte(0); n < 16; n++ {
		codewords[n] = Hamming74Encode(n)
	}
	min := 7
	for i := 0; i < 16; i++ {
		for j := i + 1; j < 16; j++ {
			d := HammingDistanceBits(uint64(codewords[i]), uint64(codewords[j]))
			if d < min {
				min = d
			}
		}
	}
	if min != 3 {
		t.Errorf("minimum distance = %d, want 3", min)
	}
}

func TestHammingWeightAndDistance(t *testing.T) {
	if got := HammingWeight(0b1011); got != 3 {
		t.Errorf("HammingWeight(0b1011) = %d, want 3", got)
	}
	if got := HammingDistanceBits(0b1010, 0b0110); got != 2 {
		t.Errorf("HammingDistanceBits = %d, want 2", got)
	}
	if got := Parity(0b111); got != 1 {
		t.Errorf("Parity(0b111) = %d, want 1", got)
	}
	if got := Parity(0b1001); got != 0 {
		t.Errorf("Parity(0b1001) = %d, want 0", got)
	}
	d, err := HammingDistanceBytes([]byte{0xFF, 0x00}, []byte{0x0F, 0x01})
	if err != nil || d != 5 {
		t.Errorf("HammingDistanceBytes = %d (err %v), want 5", d, err)
	}
	if _, err := HammingDistanceBytes([]byte{1}, []byte{1, 2}); err == nil {
		t.Error("expected length mismatch error")
	}
	ds, err := HammingDistanceStrings("karolin", "kathrin")
	if err != nil || ds != 3 {
		t.Errorf("HammingDistanceStrings = %d (err %v), want 3", ds, err)
	}
}
