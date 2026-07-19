package quasirandom

// Hammersley returns the i-th point (zero-based, 0 <= i < total) of the
// dim-dimensional Hammersley point set of total points. The first coordinate is
// i/total and the remaining dim-1 coordinates are radical inverses in the first
// dim-1 primes, giving a finite set whose star discrepancy is asymptotically
// smaller than the corresponding Halton sequence. It returns an error when
// dim < 1, total < 1 or i is out of range.
func Hammersley(dim, i, total int) ([]float64, error) {
	if dim < 1 {
		return nil, ErrDimension
	}
	if total < 1 {
		return nil, ErrNonPositive
	}
	if i < 0 || i >= total {
		return nil, ErrDimension
	}
	out := make([]float64, dim)
	out[0] = float64(i) / float64(total)
	if dim > 1 {
		bases, err := PrimeBases(dim - 1)
		if err != nil {
			return nil, err
		}
		for k := 1; k < dim; k++ {
			v, err := RadicalInverse(bases[k-1], uint64(i))
			if err != nil {
				return nil, err
			}
			out[k] = v
		}
	}
	return out, nil
}

// HammersleySet returns the complete dim-dimensional Hammersley point set of
// total points as a slice of points.
func HammersleySet(dim, total int) ([][]float64, error) {
	if dim < 1 {
		return nil, ErrDimension
	}
	if total < 1 {
		return nil, ErrNonPositive
	}
	out := make([][]float64, total)
	for i := 0; i < total; i++ {
		p, err := Hammersley(dim, i, total)
		if err != nil {
			return nil, err
		}
		out[i] = p
	}
	return out, nil
}

// HammersleyWithBases returns the i-th Hammersley point whose trailing
// coordinates use the supplied bases (one per trailing dimension, so the point
// has len(bases)+1 coordinates). The leading coordinate remains i/total.
func HammersleyWithBases(bases []int, i, total int) ([]float64, error) {
	if total < 1 {
		return nil, ErrNonPositive
	}
	if i < 0 || i >= total {
		return nil, ErrDimension
	}
	out := make([]float64, len(bases)+1)
	out[0] = float64(i) / float64(total)
	for k, b := range bases {
		v, err := RadicalInverse(b, uint64(i))
		if err != nil {
			return nil, err
		}
		out[k+1] = v
	}
	return out, nil
}
