package wavelet

// DWT2DResult holds the four subbands produced by one level of a separable
// two-dimensional discrete wavelet transform. The transform is applied first
// along rows (horizontally) and then along columns (vertically). Each subband
// is an (R/2)-by-(C/2) matrix for an R-by-C input.
type DWT2DResult struct {
	// LL is the approximation subband (low-pass along both axes).
	LL [][]float64
	// LH is the horizontal-detail subband (low-pass rows, high-pass columns).
	LH [][]float64
	// HL is the vertical-detail subband (high-pass rows, low-pass columns).
	HL [][]float64
	// HH is the diagonal-detail subband (high-pass along both axes).
	HH [][]float64
}

// waveletColumn returns column c of matrix m as a slice.
func waveletColumn(m [][]float64, c int) []float64 {
	out := make([]float64, len(m))
	for r := range m {
		out[r] = m[r][c]
	}
	return out
}

// waveletSetColumn writes v into column c of matrix m.
func waveletSetColumn(m [][]float64, c int, v []float64) {
	for r := range m {
		m[r][c] = v[r]
	}
}

// waveletNewMatrix allocates a rows-by-cols matrix of zeros.
func waveletNewMatrix(rows, cols int) [][]float64 {
	m := make([][]float64, rows)
	for r := range m {
		m[r] = make([]float64, cols)
	}
	return m
}

// waveletCheckRect verifies that image is a non-empty rectangular matrix with
// even dimensions, panicking otherwise.
func waveletCheckRect(image [][]float64) (rows, cols int) {
	rows = len(image)
	if rows == 0 {
		panic("wavelet: 2D transform requires a non-empty image")
	}
	cols = len(image[0])
	for _, row := range image {
		if len(row) != cols {
			panic("wavelet: 2D transform requires a rectangular image")
		}
	}
	if rows%2 != 0 || cols%2 != 0 || rows < 2 || cols < 2 {
		panic("wavelet: 2D transform requires even dimensions of at least 2")
	}
	return rows, cols
}

// DWT2D computes one level of the separable two-dimensional discrete wavelet
// transform of image with wavelet w, using periodic boundary handling on both
// axes. It returns the four subbands. It panics if image is not a rectangular
// matrix with even dimensions of at least two.
func DWT2D(image [][]float64, w Wavelet) *DWT2DResult {
	rows, cols := waveletCheckRect(image)
	halfC := cols / 2
	halfR := rows / 2

	// Rows pass: split each row into low (A) and high (D) halves.
	a := waveletNewMatrix(rows, halfC)
	d := waveletNewMatrix(rows, halfC)
	for r := 0; r < rows; r++ {
		lo, hi := wavelet1DForward(image[r], w.lo, w.hi)
		copy(a[r], lo)
		copy(d[r], hi)
	}

	res := &DWT2DResult{
		LL: waveletNewMatrix(halfR, halfC),
		LH: waveletNewMatrix(halfR, halfC),
		HL: waveletNewMatrix(halfR, halfC),
		HH: waveletNewMatrix(halfR, halfC),
	}
	// Columns pass on A -> LL (approx) and LH (detail).
	for c := 0; c < halfC; c++ {
		lo, hi := wavelet1DForward(waveletColumn(a, c), w.lo, w.hi)
		waveletSetColumn(res.LL, c, lo)
		waveletSetColumn(res.LH, c, hi)
	}
	// Columns pass on D -> HL (approx) and HH (detail).
	for c := 0; c < halfC; c++ {
		lo, hi := wavelet1DForward(waveletColumn(d, c), w.lo, w.hi)
		waveletSetColumn(res.HL, c, lo)
		waveletSetColumn(res.HH, c, hi)
	}
	return res
}

// IDWT2D reconstructs an image from its four subbands using wavelet w,
// inverting a single level of [DWT2D] to within floating-point rounding. It
// panics if the four subbands are not equally sized rectangular matrices.
func IDWT2D(r *DWT2DResult, w Wavelet) [][]float64 {
	halfR := len(r.LL)
	if halfR == 0 {
		panic("wavelet: IDWT2D requires non-empty subbands")
	}
	halfC := len(r.LL[0])
	for _, sb := range [][][]float64{r.LL, r.LH, r.HL, r.HH} {
		if len(sb) != halfR {
			panic("wavelet: IDWT2D subbands differ in height")
		}
		for _, row := range sb {
			if len(row) != halfC {
				panic("wavelet: IDWT2D subbands differ in width")
			}
		}
	}
	rows := halfR * 2
	cols := halfC * 2

	// Invert the column pass: rebuild A from (LL, LH) and D from (HL, HH).
	a := waveletNewMatrix(rows, halfC)
	d := waveletNewMatrix(rows, halfC)
	for c := 0; c < halfC; c++ {
		waveletSetColumn(a, c, wavelet1DInverse(waveletColumn(r.LL, c), waveletColumn(r.LH, c), w.lo, w.hi))
		waveletSetColumn(d, c, wavelet1DInverse(waveletColumn(r.HL, c), waveletColumn(r.HH, c), w.lo, w.hi))
	}

	// Invert the row pass: rebuild each full row from (A, D).
	out := waveletNewMatrix(rows, cols)
	for row := 0; row < rows; row++ {
		out[row] = wavelet1DInverse(a[row], d[row], w.lo, w.hi)
	}
	return out
}

// Dims returns the height and width of each subband in the result.
func (r *DWT2DResult) Dims() (rows, cols int) {
	rows = len(r.LL)
	if rows > 0 {
		cols = len(r.LL[0])
	}
	return rows, cols
}

// Subbands returns the four subbands in the fixed order LL, LH, HL, HH.
func (r *DWT2DResult) Subbands() (ll, lh, hl, hh [][]float64) {
	return r.LL, r.LH, r.HL, r.HH
}

// Energy returns the total energy (sum of squares) across all four subbands.
// For an orthogonal wavelet this equals the energy of the original image.
func (r *DWT2DResult) Energy() float64 {
	var e float64
	for _, sb := range [][][]float64{r.LL, r.LH, r.HL, r.HH} {
		for _, row := range sb {
			e += Energy(row)
		}
	}
	return e
}
