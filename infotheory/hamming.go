package infotheory

import (
	"errors"
	"math/bits"
)

// ErrHamming is returned by the Hamming distance helpers when their operands
// have mismatched lengths.
var ErrHamming = errors.New("infotheory: mismatched lengths")

// Hamming74Encode encodes the four data bits in the low nibble of data (bit 0
// is d1, bit 1 is d2, bit 2 is d3, bit 3 is d4) into a 7-bit Hamming(7,4)
// codeword returned in the low seven bits of a byte. The bit positions 1..7
// are laid out as p1, p2, d1, p4, d2, d3, d4 with parity bits
// p1 = d1^d2^d4, p2 = d1^d3^d4 and p4 = d2^d3^d4. Only the low four bits of
// data are used.
func Hamming74Encode(data byte) byte {
	d1 := (data >> 0) & 1
	d2 := (data >> 1) & 1
	d3 := (data >> 2) & 1
	d4 := (data >> 3) & 1
	p1 := d1 ^ d2 ^ d4
	p2 := d1 ^ d3 ^ d4
	p4 := d2 ^ d3 ^ d4
	var c byte
	c |= p1 << 0 // position 1
	c |= p2 << 1 // position 2
	c |= d1 << 2 // position 3
	c |= p4 << 3 // position 4
	c |= d2 << 4 // position 5
	c |= d3 << 5 // position 6
	c |= d4 << 6 // position 7
	return c
}

// Hamming74Syndrome returns the error-position syndrome of a 7-bit Hamming(7,4)
// codeword held in the low seven bits of code. A result of zero indicates no
// detected error; a non-zero result in 1..7 is the 1-indexed bit position at
// which a single-bit error is located.
func Hamming74Syndrome(code byte) int {
	b := func(pos uint) byte { return (code >> pos) & 1 }
	// Positions are 1-indexed; bit i of code is position i+1.
	s1 := b(0) ^ b(2) ^ b(4) ^ b(6) // positions 1,3,5,7
	s2 := b(1) ^ b(2) ^ b(5) ^ b(6) // positions 2,3,6,7
	s4 := b(3) ^ b(4) ^ b(5) ^ b(6) // positions 4,5,6,7
	return int(s1) + int(s2)<<1 + int(s4)<<2
}

// Hamming74Decode decodes a 7-bit Hamming(7,4) codeword held in the low seven
// bits of code, correcting any single-bit error. It returns the recovered four
// data bits in the low nibble, the 1-indexed position of the corrected error
// (0 if no error was detected), and a flag reporting whether a correction was
// applied. Because Hamming(7,4) has minimum distance three it corrects any
// single-bit error but cannot reliably detect two-bit errors.
func Hamming74Decode(code byte) (data byte, errorPos int, corrected bool) {
	syn := Hamming74Syndrome(code)
	if syn != 0 {
		code ^= 1 << (syn - 1) // flip the erroneous bit
		errorPos = syn
		corrected = true
	}
	d1 := (code >> 2) & 1 // position 3
	d2 := (code >> 4) & 1 // position 5
	d3 := (code >> 5) & 1 // position 6
	d4 := (code >> 6) & 1 // position 7
	data = d1<<0 | d2<<1 | d3<<2 | d4<<3
	return data, errorPos, corrected
}

// HammingWeight returns the Hamming weight of x, the number of 1-bits in its
// binary representation.
func HammingWeight(x uint64) int {
	return bits.OnesCount64(x)
}

// HammingDistanceBits returns the Hamming distance between a and b, the number
// of bit positions in which their binary representations differ.
func HammingDistanceBits(a, b uint64) int {
	return bits.OnesCount64(a ^ b)
}

// Parity returns the parity bit of x: 1 if x has an odd number of 1-bits and 0
// otherwise (even parity).
func Parity(x uint64) int {
	return bits.OnesCount64(x) & 1
}

// HammingDistanceBytes returns the bit-level Hamming distance between two equal
// length byte slices, the total number of differing bits across all bytes. It
// returns ErrHamming if the slices differ in length.
func HammingDistanceBytes(a, b []byte) (int, error) {
	if len(a) != len(b) {
		return 0, ErrHamming
	}
	var d int
	for i := range a {
		d += bits.OnesCount8(a[i] ^ b[i])
	}
	return d, nil
}

// HammingDistanceStrings returns the symbol-level Hamming distance between two
// equal-length strings, the number of positions at which the corresponding
// bytes differ. It returns ErrHamming if the strings differ in length.
func HammingDistanceStrings(a, b string) (int, error) {
	if len(a) != len(b) {
		return 0, ErrHamming
	}
	var d int
	for i := 0; i < len(a); i++ {
		if a[i] != b[i] {
			d++
		}
	}
	return d, nil
}
