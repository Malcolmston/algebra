package codingtheory

import "math"

// This file implements a static-model integer arithmetic coder using 32-bit
// range arithmetic with E1/E2/E3 renormalisation, together with cumulative
// frequency helpers.

const (
	acPrecision    = 32
	acWhole        = uint64(1) << acPrecision
	acHalf         = acWhole >> 1
	acQuarter      = acWhole >> 2
	acThreeQuarter = 3 * acQuarter
	acMask         = acWhole - 1
)

// CumulativeFreqs returns the cumulative-frequency table of a frequency slice:
// a slice of length len(freqs)+1 where entry i is the sum of freqs[:i]. The
// final entry is the total frequency.
func CumulativeFreqs(freqs []int) []int {
	cum := make([]int, len(freqs)+1)
	for i, f := range freqs {
		cum[i+1] = cum[i] + f
	}
	return cum
}

// validFreqs reports whether every frequency is positive and the total fits the
// coder's precision budget.
func validFreqs(freqs []int) bool {
	if len(freqs) == 0 {
		return false
	}
	total := 0
	for _, f := range freqs {
		if f <= 0 {
			return false
		}
		total += f
	}
	return total > 0 && uint64(total) < acQuarter
}

// ArithmeticEncode encodes a sequence of symbols (each an index into freqs)
// under the static frequency model freqs, returning the compressed bit slice.
// Frequencies must be positive and sum to less than 2^30. It returns ErrHuffman
// for an invalid model or out-of-range symbol.
func ArithmeticEncode(symbols []int, freqs []int) ([]int, error) {
	if !validFreqs(freqs) {
		return nil, ErrHuffman
	}
	cum := CumulativeFreqs(freqs)
	total := uint64(cum[len(cum)-1])
	low, high := uint64(0), acMask
	pending := 0
	var out []int
	emit := func(bit int) {
		out = append(out, bit)
		for ; pending > 0; pending-- {
			out = append(out, bit^1)
		}
	}
	for _, s := range symbols {
		if s < 0 || s >= len(freqs) {
			return nil, ErrHuffman
		}
		rng := high - low + 1
		high = low + rng*uint64(cum[s+1])/total - 1
		low = low + rng*uint64(cum[s])/total
		for {
			if high < acHalf {
				emit(0)
			} else if low >= acHalf {
				emit(1)
				low -= acHalf
				high -= acHalf
			} else if low >= acQuarter && high < acThreeQuarter {
				pending++
				low -= acQuarter
				high -= acQuarter
			} else {
				break
			}
			low = (low << 1) & acMask
			high = ((high << 1) | 1) & acMask
		}
	}
	// flush
	pending++
	if low < acQuarter {
		emit(0)
	} else {
		emit(1)
	}
	return out, nil
}

// ArithmeticDecode decodes exactly count symbols from a bit slice produced by
// ArithmeticEncode under the same frequency model.
func ArithmeticDecode(bits []int, freqs []int, count int) ([]int, error) {
	if !validFreqs(freqs) || count < 0 {
		return nil, ErrHuffman
	}
	cum := CumulativeFreqs(freqs)
	total := uint64(cum[len(cum)-1])
	pos := 0
	nextBit := func() uint64 {
		var b uint64
		if pos < len(bits) {
			b = uint64(bits[pos] & 1)
		}
		pos++
		return b
	}
	var value uint64
	for i := 0; i < acPrecision; i++ {
		value = (value << 1) | nextBit()
	}
	low, high := uint64(0), acMask
	out := make([]int, 0, count)
	for k := 0; k < count; k++ {
		rng := high - low + 1
		scaled := ((value-low+1)*total - 1) / rng
		// locate symbol s with cum[s] <= scaled < cum[s+1]
		s := 0
		for s < len(freqs) && uint64(cum[s+1]) <= scaled {
			s++
		}
		if s >= len(freqs) {
			return nil, ErrHuffman
		}
		out = append(out, s)
		high = low + rng*uint64(cum[s+1])/total - 1
		low = low + rng*uint64(cum[s])/total
		for {
			if high < acHalf {
				// nothing
			} else if low >= acHalf {
				value -= acHalf
				low -= acHalf
				high -= acHalf
			} else if low >= acQuarter && high < acThreeQuarter {
				value -= acQuarter
				low -= acQuarter
				high -= acQuarter
			} else {
				break
			}
			low = (low << 1) & acMask
			high = ((high << 1) | 1) & acMask
			value = ((value << 1) | nextBit()) & acMask
		}
	}
	return out, nil
}

// ShannonCodeLength returns the ideal (fractional) number of bits to encode a
// symbol of probability p under an optimal arithmetic coder, -log2(p). It
// returns 0 for p<=0.
func ShannonCodeLength(p float64) float64 {
	if p <= 0 {
		return 0
	}
	return -math.Log2(p)
}
