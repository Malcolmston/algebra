package infotheory

import "math"

// KraftSum returns the Kraft sum sum_i base^{-l_i} of a code whose codeword
// lengths are l_i, over an alphabet of the given base (2 for binary). Lengths
// that are non-positive are treated as zero-length and contribute base^0 = 1.
// The base must be at least two.
func KraftSum(lengths []int, base int) float64 {
	b := float64(base)
	var s float64
	for _, l := range lengths {
		s += math.Pow(b, -float64(l))
	}
	return s
}

// KraftInequality reports whether the codeword lengths satisfy the Kraft
// inequality sum_i base^{-l_i} <= 1, the necessary and sufficient condition for
// a prefix code with those lengths to exist over an alphabet of the given base.
// A small tolerance absorbs floating-point round-off.
func KraftInequality(lengths []int, base int) bool {
	return KraftSum(lengths, base) <= 1+1e-9
}

// KraftEquality reports whether the codeword lengths make the Kraft sum equal to
// one within a small tolerance, the condition characterising a complete (no
// wasted codeword) prefix code such as an optimal Huffman code.
func KraftEquality(lengths []int, base int) bool {
	return math.Abs(KraftSum(lengths, base)-1) <= 1e-9
}

// McMillanInequality reports whether the codeword lengths satisfy the
// McMillan inequality sum_i base^{-l_i} <= 1. By McMillan's theorem this is the
// same numerical condition as the Kraft inequality, but it is the necessary and
// sufficient condition for a uniquely decodable (not merely prefix) code with
// those lengths to exist.
func McMillanInequality(lengths []int, base int) bool {
	return KraftInequality(lengths, base)
}

// IsPrefixFree reports whether the set of binary codewords is prefix-free: no
// codeword is a prefix of another. A prefix-free code is uniquely and
// instantaneously decodable. Duplicate codewords make the code not prefix-free.
func IsPrefixFree(codewords []string) bool {
	for i := 0; i < len(codewords); i++ {
		for j := 0; j < len(codewords); j++ {
			if i == j {
				continue
			}
			a, b := codewords[i], codewords[j]
			if len(a) <= len(b) && b[:len(a)] == a {
				return false
			}
		}
	}
	return true
}

// FanoInequalityBound returns the Fano inequality's upper bound on the
// conditional entropy H(X|Y) in bits, given the probability of error errorProb
// of an estimator of X and the size alphabet of the alphabet of X:
// H_b(P_e) + P_e log2(|X|-1). Any decoder achieving error probability errorProb
// must have conditional entropy no greater than this bound. errorProb must lie
// in [0,1] and alphabet must be at least one.
func FanoInequalityBound(errorProb float64, alphabet int) float64 {
	bound := BinaryEntropy(errorProb)
	if alphabet > 1 {
		bound += errorProb * infotheoryLog2(float64(alphabet-1))
	}
	return bound
}

// CodeRate returns the rate k/n of a block code that maps k information symbols
// to n coded symbols. It returns zero for a non-positive n.
func CodeRate(k, n int) float64 {
	if n <= 0 {
		return 0
	}
	return float64(k) / float64(n)
}

// EntropyRate returns the entropy rate H(p)/L in bits per symbol of a source
// with per-block entropy H(p) emitting blocks of length L symbols. It returns
// zero for a non-positive L.
func EntropyRate(p []float64, blockLength int) float64 {
	if blockLength <= 0 {
		return 0
	}
	return Entropy(p) / float64(blockLength)
}

// CodingEfficiency returns the efficiency of a code, the ratio of the source
// entropy in bits to the average codeword length: H / avgLength. A value of one
// indicates an optimal code that attains the entropy bound. It returns zero for
// a non-positive average length.
func CodingEfficiency(entropyBits, avgLength float64) float64 {
	if avgLength <= 0 {
		return 0
	}
	return entropyBits / avgLength
}
