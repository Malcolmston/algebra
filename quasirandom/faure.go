package quasirandom

// FaureBase returns the base used by the Faure sequence in the given dimension:
// the smallest prime greater than or equal to the dimension (and at least two).
// It returns an error when dim < 1.
func FaureBase(dim int) (int, error) {
	if dim < 1 {
		return 0, ErrDimension
	}
	if dim < 2 {
		return 2, nil
	}
	return NextPrimeGE(dim), nil
}

// binomialMod returns Pascal's triangle of binomial coefficients modulo base,
// as a size-by-size lower-triangular table with tab[k][i] = C(k,i) mod base.
func binomialMod(size, base int) [][]int {
	tab := make([][]int, size)
	for k := 0; k < size; k++ {
		tab[k] = make([]int, size)
		tab[k][0] = 1 % base
		for i := 1; i <= k; i++ {
			tab[k][i] = (tab[k-1][i-1] + tab[k-1][i]) % base
		}
	}
	return tab
}

// FaureGeneratorMatrix returns the power-th power of the upper-triangular Pascal
// generator matrix modulo base, truncated to size rows and columns. Entry (i,j)
// equals C(j,i) * power^(j-i) mod base for j >= i and zero below the diagonal,
// which is the generating matrix applied to coordinate index power of a Faure
// sequence in the given base. It returns an error when base < 2 or size < 0.
func FaureGeneratorMatrix(base, power, size int) ([][]int, error) {
	if base < 2 {
		return nil, ErrBadBase
	}
	if size < 0 {
		return nil, ErrNonPositive
	}
	binom := binomialMod(size, base)
	// Precompute powers of `power` modulo base.
	pw := make([]int, size+1)
	pw[0] = 1 % base
	pm := ((power % base) + base) % base
	for e := 1; e <= size; e++ {
		pw[e] = (pw[e-1] * pm) % base
	}
	m := make([][]int, size)
	for i := 0; i < size; i++ {
		m[i] = make([]int, size)
		for j := i; j < size; j++ {
			m[i][j] = (binom[j][i] * pw[j-i]) % base
		}
	}
	return m, nil
}

// faureCoordinate transforms the base-b digit vector a (least significant
// first) into the output digit vector for coordinate index d, computing
// y_i = sum_{k>=i} C(k,i) d^{k-i} a_k (mod base), then reflecting to a fraction.
func faureCoordinate(base, d int, a []int, binom [][]int) float64 {
	m := len(a)
	// Powers of d modulo base up to m.
	pw := make([]int, m+1)
	pw[0] = 1 % base
	dm := d % base
	for e := 1; e <= m; e++ {
		pw[e] = (pw[e-1] * dm) % base
	}
	inv := 1 / float64(base)
	factor := inv
	var result float64
	for i := 0; i < m; i++ {
		y := 0
		for k := i; k < m; k++ {
			y = (y + binom[k][i]*pw[k-i]*a[k]) % base
		}
		result += float64(y) * factor
		factor *= inv
	}
	return result
}

// Faure returns the n-th point (zero-based) of the Faure sequence in the given
// dimension. All coordinates share the base returned by FaureBase and are the
// images of n under successive powers of the Pascal generator matrix. The point
// lies in [0,1)^dim. It returns an error when dim < 1.
func Faure(dim int, n uint64) ([]float64, error) {
	if dim < 1 {
		return nil, ErrDimension
	}
	base, err := FaureBase(dim)
	if err != nil {
		return nil, err
	}
	return FaureWithBase(dim, base, n)
}

// FaureWithBase returns the n-th Faure point using an explicit prime base, which
// must be at least as large as dim for the construction to be a genuine
// (t,s)-sequence. It returns an error when dim < 1 or base < 2.
func FaureWithBase(dim, base int, n uint64) ([]float64, error) {
	if dim < 1 {
		return nil, ErrDimension
	}
	if base < 2 {
		return nil, ErrBadBase
	}
	a, err := Digits(n, base)
	if err != nil {
		return nil, err
	}
	if len(a) == 0 {
		return make([]float64, dim), nil
	}
	binom := binomialMod(len(a), base)
	out := make([]float64, dim)
	for d := 0; d < dim; d++ {
		out[d] = faureCoordinate(base, d, a, binom)
	}
	return out, nil
}

// FaureSequence returns the first count points of the dim-dimensional Faure
// sequence beginning at index zero.
func FaureSequence(dim, count int) ([][]float64, error) {
	if dim < 1 {
		return nil, ErrDimension
	}
	if count < 0 {
		return nil, ErrNonPositive
	}
	base, err := FaureBase(dim)
	if err != nil {
		return nil, err
	}
	out := make([][]float64, count)
	for i := 0; i < count; i++ {
		p, err := FaureWithBase(dim, base, uint64(i))
		if err != nil {
			return nil, err
		}
		out[i] = p
	}
	return out, nil
}

// FaureSequenceOffset returns count Faure points starting at index start.
func FaureSequenceOffset(dim, start, count int) ([][]float64, error) {
	if dim < 1 {
		return nil, ErrDimension
	}
	if count < 0 || start < 0 {
		return nil, ErrNonPositive
	}
	base, err := FaureBase(dim)
	if err != nil {
		return nil, err
	}
	out := make([][]float64, count)
	for i := 0; i < count; i++ {
		p, err := FaureWithBase(dim, base, uint64(start+i))
		if err != nil {
			return nil, err
		}
		out[i] = p
	}
	return out, nil
}
