package infotheory

import (
	"reflect"
	"testing"
)

func TestHuffmanAverageLength(t *testing.T) {
	// Classic example: optimal average code length is 2.2 bits and the code
	// lengths saturate the Kraft equality.
	probs := map[string]float64{"a": 0.4, "b": 0.2, "c": 0.2, "d": 0.1, "e": 0.1}
	codes, err := HuffmanCodes(probs)
	if err != nil {
		t.Fatal(err)
	}
	if len(codes) != 5 {
		t.Fatalf("expected 5 codewords, got %d", len(codes))
	}
	avg := AverageCodeLength(codes, probs)
	if !approx(avg, 2.2, 1e-9) {
		t.Errorf("average code length = %v, want 2.2", avg)
	}
	// Optimal code length must be within [H, H+1).
	h := Entropy([]float64{0.4, 0.2, 0.2, 0.1, 0.1})
	if avg < h-1e-9 || avg >= h+1 {
		t.Errorf("avg length %v not in [H=%v, H+1)", avg, h)
	}
	// Codewords must be prefix-free and complete.
	cw := make([]string, 0, len(codes))
	lengths := make([]int, 0, len(codes))
	for _, c := range codes {
		cw = append(cw, c)
		lengths = append(lengths, len(c))
	}
	if !IsPrefixFree(cw) {
		t.Error("Huffman codewords are not prefix-free")
	}
	if !KraftEquality(lengths, 2) {
		t.Error("Huffman codeword lengths do not satisfy Kraft equality")
	}
}

func TestHuffmanEncodeDecodeRoundTrip(t *testing.T) {
	probs := map[string]float64{"a": 0.5, "b": 0.25, "c": 0.125, "d": 0.125}
	root, err := BuildHuffmanTree(probs)
	if err != nil {
		t.Fatal(err)
	}
	codes := root.Codewords()
	msg := []string{"a", "b", "a", "c", "d", "a", "a", "b"}
	bits, err := HuffmanEncode(codes, msg)
	if err != nil {
		t.Fatal(err)
	}
	got, err := HuffmanDecode(root, bits)
	if err != nil {
		t.Fatal(err)
	}
	if !reflect.DeepEqual(got, msg) {
		t.Errorf("round trip = %v, want %v", got, msg)
	}
	// For this dyadic distribution Huffman attains the entropy exactly.
	avg := AverageCodeLength(codes, probs)
	if !approx(avg, 1.75, 1e-9) {
		t.Errorf("avg length = %v, want 1.75", avg)
	}
}

func TestHuffmanSingleSymbol(t *testing.T) {
	root, err := BuildHuffmanTree(map[string]float64{"x": 1})
	if err != nil {
		t.Fatal(err)
	}
	codes := root.Codewords()
	if codes["x"] != "0" {
		t.Errorf("single symbol code = %q, want \"0\"", codes["x"])
	}
	bits, _ := HuffmanEncode(codes, []string{"x", "x", "x"})
	got, err := HuffmanDecode(root, bits)
	if err != nil {
		t.Fatal(err)
	}
	if !reflect.DeepEqual(got, []string{"x", "x", "x"}) {
		t.Errorf("single symbol round trip = %v", got)
	}
}

func TestHuffmanErrors(t *testing.T) {
	if _, err := BuildHuffmanTree(map[string]float64{}); err == nil {
		t.Error("expected error for empty weights")
	}
	codes := map[string]string{"a": "0"}
	if _, err := HuffmanEncode(codes, []string{"b"}); err == nil {
		t.Error("expected error encoding unknown symbol")
	}
	root, _ := BuildHuffmanTree(map[string]float64{"a": 0.5, "b": 0.5})
	if _, err := HuffmanDecode(root, "0x1"); err == nil {
		t.Error("expected error decoding invalid bit")
	}
}

func TestShannonFanoCodes(t *testing.T) {
	probs := map[string]float64{"a": 0.4, "b": 0.2, "c": 0.2, "d": 0.1, "e": 0.1}
	codes, err := ShannonFanoCodes(probs)
	if err != nil {
		t.Fatal(err)
	}
	cw := make([]string, 0, len(codes))
	for _, c := range codes {
		cw = append(cw, c)
	}
	if !IsPrefixFree(cw) {
		t.Errorf("Shannon-Fano codes not prefix-free: %v", codes)
	}
	// Shannon-Fano is a valid code so its average length is at least the
	// entropy and (being suboptimal) no shorter than Huffman's.
	avg := AverageCodeLength(codes, probs)
	huff, _ := HuffmanCodes(probs)
	if avg < AverageCodeLength(huff, probs)-1e-9 {
		t.Errorf("Shannon-Fano avg %v shorter than Huffman", avg)
	}
}
