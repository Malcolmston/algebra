package quasirandom

// Halton returns the n-th point (zero-based) of the Halton sequence in the
// given dimension, using the first dim prime numbers as the coordinate bases.
// The returned point lies in [0,1)^dim. It returns an error when dim < 1.
func Halton(dim int, n uint64) ([]float64, error) {
	if dim < 1 {
		return nil, ErrDimension
	}
	bases, err := PrimeBases(dim)
	if err != nil {
		return nil, err
	}
	return HaltonWithBases(bases, n)
}

// HaltonWithBases returns the n-th Halton point using the supplied coordinate
// bases, one per dimension. The bases should be pairwise coprime for the point
// set to be well distributed. It returns an error when no bases are given or a
// base is smaller than two.
func HaltonWithBases(bases []int, n uint64) ([]float64, error) {
	if len(bases) == 0 {
		return nil, ErrDimension
	}
	out := make([]float64, len(bases))
	for i, b := range bases {
		v, err := RadicalInverse(b, n)
		if err != nil {
			return nil, err
		}
		out[i] = v
	}
	return out, nil
}

// HaltonCoordinate returns coordinate k (zero-based) of the n-th Halton point in
// the given dimension, the radical inverse of n in the (k+1)-th prime. It
// returns an error when dim < 1 or k is out of range.
func HaltonCoordinate(dim, k int, n uint64) (float64, error) {
	if dim < 1 {
		return 0, ErrDimension
	}
	if k < 0 || k >= dim {
		return 0, ErrDimension
	}
	p, err := Prime(k + 1)
	if err != nil {
		return 0, err
	}
	return RadicalInverse(p, n)
}

// HaltonSequence returns the first count points of the dim-dimensional Halton
// sequence starting from index zero as a slice of points.
func HaltonSequence(dim, count int) ([][]float64, error) {
	if dim < 1 {
		return nil, ErrDimension
	}
	if count < 0 {
		return nil, ErrNonPositive
	}
	bases, err := PrimeBases(dim)
	if err != nil {
		return nil, err
	}
	out := make([][]float64, count)
	for i := 0; i < count; i++ {
		p, err := HaltonWithBases(bases, uint64(i))
		if err != nil {
			return nil, err
		}
		out[i] = p
	}
	return out, nil
}

// HaltonSequenceOffset returns count Halton points beginning at index start.
// Skipping a leading burn-in of a few dozen points is a common way to reduce
// the correlation between the low-index terms of high-dimensional Halton
// sequences.
func HaltonSequenceOffset(dim, start, count int) ([][]float64, error) {
	if dim < 1 {
		return nil, ErrDimension
	}
	if count < 0 || start < 0 {
		return nil, ErrNonPositive
	}
	bases, err := PrimeBases(dim)
	if err != nil {
		return nil, err
	}
	out := make([][]float64, count)
	for i := 0; i < count; i++ {
		p, err := HaltonWithBases(bases, uint64(start+i))
		if err != nil {
			return nil, err
		}
		out[i] = p
	}
	return out, nil
}

// LeapedHalton returns the n-th point of a leaped Halton sequence, which takes
// every leap-th point of the ordinary sequence (index n*leap). Leaping by a
// value coprime to the product of the bases is a simple decorrelation device.
// It returns an error when dim < 1 or leap < 1.
func LeapedHalton(dim int, leap, n uint64) ([]float64, error) {
	if dim < 1 {
		return nil, ErrDimension
	}
	if leap < 1 {
		return nil, ErrNonPositive
	}
	return Halton(dim, n*leap)
}

// ScrambledHalton returns the n-th Halton point in which each coordinate is a
// permutation-scrambled radical inverse. perms must supply one permutation per
// dimension, each a bijection of {0,...,base-1} for the corresponding prime
// base. It returns an error when the shapes disagree.
func ScrambledHalton(perms [][]int, n uint64) ([]float64, error) {
	dim := len(perms)
	if dim < 1 {
		return nil, ErrDimension
	}
	bases, err := PrimeBases(dim)
	if err != nil {
		return nil, err
	}
	out := make([]float64, dim)
	for i := 0; i < dim; i++ {
		v, err := ScrambledRadicalInverse(bases[i], n, perms[i])
		if err != nil {
			return nil, err
		}
		out[i] = v
	}
	return out, nil
}

// ReverseHalton returns the n-th point of the reverse (Faure–Tezuka) Halton
// sequence, using the deterministic per-base scrambling permutations returned
// by FaurePermutation. This variant markedly improves the uniformity of the
// early points in moderate to high dimensions. It returns an error when dim<1.
func ReverseHalton(dim int, n uint64) ([]float64, error) {
	if dim < 1 {
		return nil, ErrDimension
	}
	bases, err := PrimeBases(dim)
	if err != nil {
		return nil, err
	}
	out := make([]float64, dim)
	for i := 0; i < dim; i++ {
		perm := FaurePermutation(bases[i])
		v, err := ScrambledRadicalInverse(bases[i], n, perm)
		if err != nil {
			return nil, err
		}
		out[i] = v
	}
	return out, nil
}
