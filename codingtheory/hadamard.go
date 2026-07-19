package codingtheory

// This file implements Sylvester-Hadamard matrices, the fast Walsh-Hadamard
// transform, Walsh spreading codes, and the binary Hadamard code together with
// its maximum-likelihood decoder.

// IsPowerOfTwo reports whether n is a positive power of two.
func IsPowerOfTwo(n int) bool { return n > 0 && n&(n-1) == 0 }

// log2Exact returns the base-2 logarithm of a power of two, or -1 otherwise.
func log2Exact(n int) int {
	if !IsPowerOfTwo(n) {
		return -1
	}
	k := 0
	for n > 1 {
		n >>= 1
		k++
	}
	return k
}

// HadamardMatrix returns the 2^k-by-2^k Sylvester-Hadamard matrix with entries
// +1 and -1, built recursively from H_1 = [[1]] and the Kronecker doubling
// H_{2n} = [[H,H],[H,-H]].
func HadamardMatrix(k int) [][]int {
	n := 1 << uint(k)
	h := make([][]int, n)
	for i := range h {
		h[i] = make([]int, n)
		for j := range h[i] {
			// entry = (-1)^{popcount(i & j)}
			if HammingWeightUint(uint64(i&j))&1 == 0 {
				h[i][j] = 1
			} else {
				h[i][j] = -1
			}
		}
	}
	return h
}

// WalshHadamardTransform returns the natural-ordered fast Walsh-Hadamard
// transform of a, whose length must be a power of two. The transform is its own
// inverse up to the scale factor len(a) (see InverseWalshHadamardTransform).
func WalshHadamardTransform(a []float64) []float64 {
	n := len(a)
	out := append([]float64(nil), a...)
	for span := 1; span < n; span <<= 1 {
		for i := 0; i < n; i += span << 1 {
			for j := i; j < i+span; j++ {
				x, y := out[j], out[j+span]
				out[j] = x + y
				out[j+span] = x - y
			}
		}
	}
	return out
}

// InverseWalshHadamardTransform returns the inverse Walsh-Hadamard transform,
// i.e. the forward transform divided by len(a).
func InverseWalshHadamardTransform(a []float64) []float64 {
	out := WalshHadamardTransform(a)
	n := float64(len(a))
	for i := range out {
		out[i] /= n
	}
	return out
}

// WalshHadamardTransformInt is the integer-valued forward Walsh-Hadamard
// transform for a slice whose length is a power of two.
func WalshHadamardTransformInt(a []int) []int {
	n := len(a)
	out := append([]int(nil), a...)
	for span := 1; span < n; span <<= 1 {
		for i := 0; i < n; i += span << 1 {
			for j := i; j < i+span; j++ {
				x, y := out[j], out[j+span]
				out[j] = x + y
				out[j+span] = x - y
			}
		}
	}
	return out
}

// WalshCode returns the index-th Walsh code of length 2^k as a bit slice, taken
// from the rows of the Sylvester-Hadamard matrix (a +1 entry maps to bit 0 and
// -1 to bit 1). Distinct Walsh codes of the same length are mutually orthogonal.
func WalshCode(k, index int) []int {
	n := 1 << uint(k)
	out := make([]int, n)
	for j := 0; j < n; j++ {
		if HammingWeightUint(uint64(index&j))&1 != 0 {
			out[j] = 1
		}
	}
	return out
}

// HadamardEncode encodes an m-bit message into a length-2^m binary Hadamard
// code word: the bit at position x equals the GF(2) inner product of the
// message with the binary expansion of x. The code is linear [2^m, m, 2^{m-1}].
func HadamardEncode(msg []int) ([]int, error) {
	if !validBits(msg) {
		return nil, ErrBit
	}
	m := len(msg)
	if m == 0 || m > 20 {
		return nil, ErrFieldParam
	}
	n := 1 << uint(m)
	out := make([]int, n)
	for x := 0; x < n; x++ {
		bit := 0
		for i := 0; i < m; i++ {
			if x&(1<<uint(i)) != 0 {
				bit ^= msg[i] & 1
			}
		}
		out[x] = bit
	}
	return out, nil
}

// HadamardDecode maximum-likelihood decodes a length-2^m Hadamard code word by
// mapping bits to +/-1 and taking the fast Walsh-Hadamard transform: the index
// of the largest-magnitude coefficient is the message. It returns the recovered
// m message bits. The word length must be a power of two.
func HadamardDecode(word []int) ([]int, error) {
	m := log2Exact(len(word))
	if m < 1 {
		return nil, ErrLength
	}
	if !validBits(word) {
		return nil, ErrBit
	}
	signal := make([]float64, len(word))
	for i, b := range word {
		if b&1 == 0 {
			signal[i] = 1
		} else {
			signal[i] = -1
		}
	}
	tr := WalshHadamardTransform(signal)
	best, bestVal := 0, -1.0
	for i, v := range tr {
		av := v
		if av < 0 {
			av = -av
		}
		if av > bestVal {
			bestVal = av
			best = i
		}
	}
	// best is the message value; unpack its low m bits little-endian.
	msg := make([]int, m)
	for i := 0; i < m; i++ {
		msg[i] = (best >> uint(i)) & 1
	}
	return msg, nil
}
