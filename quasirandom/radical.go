package quasirandom

import "math"

// RadicalInverse returns the radical inverse phi_base(n): the value obtained by
// reflecting the base-b digits of n about the radix point, mapping the integer
// with digits d0 d1 d2 ... to the fraction 0.d0 d1 d2 ... in base b. The result
// lies in [0,1). It returns an error when base < 2.
func RadicalInverse(base int, n uint64) (float64, error) {
	if base < 2 {
		return 0, ErrBadBase
	}
	b := float64(base)
	inv := 1 / b
	factor := inv
	var result float64
	ub := uint64(base)
	for n > 0 {
		d := n % ub
		result += float64(d) * factor
		factor *= inv
		n /= ub
	}
	return result, nil
}

// RadicalInverseBase2 returns the radical inverse in base two using pure
// bit reversal, faster than the general routine and exact for indices below
// 2^52.
func RadicalInverseBase2(n uint64) float64 {
	var result float64
	factor := 0.5
	for n > 0 {
		result += float64(n&1) * factor
		factor *= 0.5
		n >>= 1
	}
	return result
}

// ScrambledRadicalInverse returns the radical inverse of n in the given base
// with each digit d replaced by perm[d] before reflection. The permutation must
// be a bijection of {0,...,base-1}; supplying the identity reproduces
// RadicalInverse. It returns an error when base < 2 or perm has the wrong
// length. A permutation fixing zero keeps the result in [0,1); a permutation
// moving zero adds the geometric tail perm[0]/(base-1) so the result still lies
// in [0,1).
func ScrambledRadicalInverse(base int, n uint64, perm []int) (float64, error) {
	if base < 2 {
		return 0, ErrBadBase
	}
	if len(perm) != base {
		return 0, ErrDimension
	}
	b := float64(base)
	inv := 1 / b
	factor := inv
	ub := uint64(base)
	var result float64
	// Reflect the digits actually present in n.
	digitsSeen := 0
	m := n
	for m > 0 {
		d := int(m % ub)
		result += float64(perm[d]) * factor
		factor *= inv
		m /= ub
		digitsSeen++
	}
	// If perm[0] != 0 the (infinitely many) leading zero digits of n each
	// contribute perm[0]*base^-k for k beyond digitsSeen; sum the geometric
	// tail so the mapping remains a bijection onto a subset of [0,1).
	if perm[0] != 0 {
		tail := float64(perm[0]) * factor / (b - 1)
		result += tail
	}
	_ = n
	return result, nil
}

// VanDerCorput returns the n-th term (zero-based) of the van der Corput
// sequence in the given base, which is exactly RadicalInverse(base, n).
func VanDerCorput(base int, n uint64) (float64, error) {
	return RadicalInverse(base, n)
}

// VanDerCorputSequence returns the first count terms of the van der Corput
// sequence in the given base, starting from index zero (whose value is 0).
func VanDerCorputSequence(base, count int) ([]float64, error) {
	if base < 2 {
		return nil, ErrBadBase
	}
	if count < 0 {
		return nil, ErrNonPositive
	}
	out := make([]float64, count)
	for i := 0; i < count; i++ {
		v, err := RadicalInverse(base, uint64(i))
		if err != nil {
			return nil, err
		}
		out[i] = v
	}
	return out, nil
}

// VanDerCorputSequenceOffset returns count terms of the van der Corput sequence
// in the given base beginning at index start rather than zero.
func VanDerCorputSequenceOffset(base, start, count int) ([]float64, error) {
	if base < 2 {
		return nil, ErrBadBase
	}
	if count < 0 || start < 0 {
		return nil, ErrNonPositive
	}
	out := make([]float64, count)
	for i := 0; i < count; i++ {
		v, err := RadicalInverse(base, uint64(start+i))
		if err != nil {
			return nil, err
		}
		out[i] = v
	}
	return out, nil
}

// RadicalInverseDigits returns the reflected base-b digit vector of n, ordered
// from the most significant original digit to the least; interpreting the slice
// as 0.d0 d1 d2 ... reproduces RadicalInverse. It returns an error when base<2.
func RadicalInverseDigits(base int, n uint64) ([]int, error) {
	return Digits(n, base)
}

// FoldedRadicalInverse returns the folded (Hammersley–Kronecker) radical
// inverse, in which digit position k contributes (d_k + k) mod base rather than
// d_k. The folding removes the correlation artefacts of the plain radical
// inverse for large bases. It returns an error when base < 2.
func FoldedRadicalInverse(base int, n uint64) (float64, error) {
	if base < 2 {
		return 0, ErrBadBase
	}
	b := float64(base)
	inv := 1 / b
	factor := inv
	ub := uint64(base)
	var result float64
	offset := 0
	// Fold across enough digit positions to reach the working precision.
	for factor > 1e-17 {
		d := int(n % ub)
		result += float64((d+offset)%base) * factor
		factor *= inv
		n /= ub
		offset++
	}
	return result, nil
}

// radicalInverseInverse recovers the smallest index whose radical inverse in
// the given base rounds to x within width digits. It is used by tests and by
// the leaped Halton helper.
func radicalInverseIndex(base int, x float64, width int) uint64 {
	b := uint64(base)
	var n uint64
	var p uint64 = 1
	for i := 0; i < width; i++ {
		x *= float64(base)
		d := uint64(math.Floor(x + 1e-9))
		if d >= b {
			d = b - 1
		}
		x -= float64(d)
		n += d * p
		p *= b
	}
	return n
}
