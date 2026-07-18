package transform

import (
	"math"
	"math/cmplx"
)

// IsPow2 reports whether n is a positive power of two.
func IsPow2(n int) bool {
	return n > 0 && n&(n-1) == 0
}

// NextPow2 returns the smallest power of two that is greater than or equal to
// n. For n <= 1 it returns 1.
func NextPow2(n int) int {
	if n <= 1 {
		return 1
	}
	p := 1
	for p < n {
		p <<= 1
	}
	return p
}

// transformFFTInPlace performs an in-place radix-2 Cooley-Tukey FFT on x whose
// length must be a power of two. When inverse is true the forward twiddle
// factors are conjugated; the caller is responsible for the 1/n scaling.
func transformFFTInPlace(x []complex128, inverse bool) {
	n := len(x)
	// Bit-reversal permutation.
	for i, j := 1, 0; i < n; i++ {
		bit := n >> 1
		for ; j&bit != 0; bit >>= 1 {
			j ^= bit
		}
		j ^= bit
		if i < j {
			x[i], x[j] = x[j], x[i]
		}
	}
	for length := 2; length <= n; length <<= 1 {
		ang := 2 * math.Pi / float64(length)
		if !inverse {
			ang = -ang
		}
		wlen := cmplx.Rect(1, ang)
		half := length >> 1
		for i := 0; i < n; i += length {
			w := complex(1, 0)
			for k := 0; k < half; k++ {
				u := x[i+k]
				v := x[i+k+half] * w
				x[i+k] = u + v
				x[i+k+half] = u - v
				w *= wlen
			}
		}
	}
}

// DFT computes the discrete Fourier transform of x directly in O(n^2) time,
// returning X[k] = sum_n x[n] exp(-2*pi*i*k*n/N). It works for any length and
// serves as the reference implementation for the faster routines.
func DFT(x []complex128) []complex128 {
	n := len(x)
	out := make([]complex128, n)
	for k := 0; k < n; k++ {
		var sum complex128
		for j := 0; j < n; j++ {
			ang := -2 * math.Pi * float64(k) * float64(j) / float64(n)
			sum += x[j] * cmplx.Rect(1, ang)
		}
		out[k] = sum
	}
	return out
}

// IDFT computes the inverse discrete Fourier transform of X directly,
// returning x[n] = (1/N) sum_k X[k] exp(+2*pi*i*k*n/N). It works for any
// length.
func IDFT(X []complex128) []complex128 {
	n := len(X)
	out := make([]complex128, n)
	for j := 0; j < n; j++ {
		var sum complex128
		for k := 0; k < n; k++ {
			ang := 2 * math.Pi * float64(k) * float64(j) / float64(n)
			sum += X[k] * cmplx.Rect(1, ang)
		}
		out[j] = sum / complex(float64(n), 0)
	}
	return out
}

// FFT computes the discrete Fourier transform of x using the radix-2
// Cooley-Tukey algorithm. The length of x must be a power of two; FFT panics
// otherwise. Use [FFTAny] to transform sequences of arbitrary length.
func FFT(x []complex128) []complex128 {
	if !IsPow2(len(x)) {
		panic("transform: FFT requires a power-of-two length; use FFTAny")
	}
	out := make([]complex128, len(x))
	copy(out, x)
	transformFFTInPlace(out, false)
	return out
}

// IFFT computes the inverse radix-2 FFT of X, including the 1/N scaling. The
// length of X must be a power of two; IFFT panics otherwise.
func IFFT(X []complex128) []complex128 {
	if !IsPow2(len(X)) {
		panic("transform: IFFT requires a power-of-two length; use IFFTAny")
	}
	out := make([]complex128, len(X))
	copy(out, X)
	transformFFTInPlace(out, true)
	scale := complex(1/float64(len(out)), 0)
	for i := range out {
		out[i] *= scale
	}
	return out
}

// FFTAny computes the discrete Fourier transform of x for any length. When the
// length is a power of two it uses the radix-2 [FFT]; otherwise it uses
// [Bluestein]'s algorithm, which runs in O(n log n) time.
func FFTAny(x []complex128) []complex128 {
	if IsPow2(len(x)) {
		return FFT(x)
	}
	return Bluestein(x)
}

// IFFTAny computes the inverse discrete Fourier transform of X for any length,
// including the 1/N scaling.
func IFFTAny(X []complex128) []complex128 {
	n := len(X)
	if n == 0 {
		return []complex128{}
	}
	if IsPow2(n) {
		return IFFT(X)
	}
	conj := make([]complex128, n)
	for i, v := range X {
		conj[i] = cmplx.Conj(v)
	}
	y := FFTAny(conj)
	scale := complex(1/float64(n), 0)
	out := make([]complex128, n)
	for i, v := range y {
		out[i] = cmplx.Conj(v) * scale
	}
	return out
}

// FFTReal computes the discrete Fourier transform of a real-valued signal x,
// returning the full complex spectrum of the same length. It accepts any
// length.
func FFTReal(x []float64) []complex128 {
	c := make([]complex128, len(x))
	for i, v := range x {
		c[i] = complex(v, 0)
	}
	return FFTAny(c)
}

// RFFT computes the discrete Fourier transform of a real-valued signal x and
// returns only the non-redundant first floor(N/2)+1 frequency bins. The
// remaining bins are the complex conjugates of these by Hermitian symmetry.
// Use [IRFFT] to invert the result.
func RFFT(x []float64) []complex128 {
	full := FFTReal(x)
	m := len(x)/2 + 1
	out := make([]complex128, m)
	copy(out, full[:m])
	return out
}

// IRFFT reconstructs a real signal of length n from the floor(n/2)+1
// half-spectrum X produced by [RFFT]. The original length n must be supplied
// because it cannot be recovered from the half-spectrum alone.
func IRFFT(X []complex128, n int) []float64 {
	full := make([]complex128, n)
	half := n/2 + 1
	if len(X) < half {
		half = len(X)
	}
	for k := 0; k < half; k++ {
		full[k] = X[k]
	}
	for k := 1; k < n-half+1; k++ {
		full[n-k] = cmplx.Conj(X[k])
	}
	inv := IFFTAny(full)
	out := make([]float64, n)
	for i := range inv {
		out[i] = real(inv[i])
	}
	return out
}

// ZeroPad returns a copy of x extended with zeros to length n. If n is smaller
// than len(x) the input is truncated.
func ZeroPad(x []complex128, n int) []complex128 {
	out := make([]complex128, n)
	copy(out, x)
	return out
}

// ZeroPadReal returns a copy of the real slice x extended with zeros to length
// n. If n is smaller than len(x) the input is truncated.
func ZeroPadReal(x []float64, n int) []float64 {
	out := make([]float64, n)
	copy(out, x)
	return out
}

// FFTShift rearranges a spectrum so that the zero-frequency component is moved
// to the center of the slice, matching the conventional plotting order. It is
// inverted by [IFFTShift].
func FFTShift(x []complex128) []complex128 {
	n := len(x)
	out := make([]complex128, n)
	s := n / 2
	for i := 0; i < n; i++ {
		out[(i+s)%n] = x[i]
	}
	return out
}

// IFFTShift is the inverse of [FFTShift]; it moves the zero-frequency
// component from the center back to the start of the slice.
func IFFTShift(x []complex128) []complex128 {
	n := len(x)
	out := make([]complex128, n)
	s := n - n/2
	for i := 0; i < n; i++ {
		out[(i+s)%n] = x[i]
	}
	return out
}

// FFTFreq returns the n sample frequencies corresponding to the bins of an
// n-point DFT taken with sample spacing d (in seconds). The result follows the
// standard convention: bins 0..(n-1)/2 hold the non-negative frequencies and
// the remainder hold the negative frequencies.
func FFTFreq(n int, d float64) []float64 {
	out := make([]float64, n)
	if n == 0 {
		return out
	}
	val := 1.0 / (float64(n) * d)
	half := (n - 1) / 2
	for i := 0; i <= half; i++ {
		out[i] = float64(i) * val
	}
	for i := half + 1; i < n; i++ {
		out[i] = float64(i-n) * val
	}
	return out
}

// RFFTFreq returns the floor(n/2)+1 non-negative sample frequencies
// corresponding to the bins produced by [RFFT] for an n-point signal with
// sample spacing d (in seconds).
func RFFTFreq(n int, d float64) []float64 {
	m := n/2 + 1
	out := make([]float64, m)
	val := 1.0 / (float64(n) * d)
	for i := 0; i < m; i++ {
		out[i] = float64(i) * val
	}
	return out
}

// FFT2D computes the two-dimensional discrete Fourier transform of a matrix by
// transforming every row and then every column. The matrix must be
// rectangular; each dimension may have any length.
func FFT2D(m [][]complex128) [][]complex128 {
	return transformFFT2D(m, false)
}

// IFFT2D computes the inverse two-dimensional discrete Fourier transform,
// including the 1/(rows*cols) scaling. The matrix must be rectangular.
func IFFT2D(m [][]complex128) [][]complex128 {
	return transformFFT2D(m, true)
}

func transformFFT2D(m [][]complex128, inverse bool) [][]complex128 {
	rows := len(m)
	if rows == 0 {
		return [][]complex128{}
	}
	cols := len(m[0])
	fwd := FFTAny
	if inverse {
		fwd = IFFTAny
	}
	out := make([][]complex128, rows)
	for i := range m {
		out[i] = fwd(m[i])
	}
	col := make([]complex128, rows)
	for j := 0; j < cols; j++ {
		for i := 0; i < rows; i++ {
			col[i] = out[i][j]
		}
		tc := fwd(col)
		for i := 0; i < rows; i++ {
			out[i][j] = tc[i]
		}
	}
	return out
}

// DFTMatrix returns the n-by-n discrete Fourier transform matrix W whose entry
// W[j][k] equals exp(-2*pi*i*j*k/n). Multiplying this matrix by a column
// vector reproduces [DFT].
func DFTMatrix(n int) [][]complex128 {
	w := make([][]complex128, n)
	for j := 0; j < n; j++ {
		w[j] = make([]complex128, n)
		for k := 0; k < n; k++ {
			ang := -2 * math.Pi * float64(j) * float64(k) / float64(n)
			w[j][k] = cmplx.Rect(1, ang)
		}
	}
	return w
}

// Magnitude returns the element-wise magnitudes (absolute values) of a complex
// spectrum.
func Magnitude(X []complex128) []float64 {
	out := make([]float64, len(X))
	for i, v := range X {
		out[i] = cmplx.Abs(v)
	}
	return out
}

// Phase returns the element-wise phase angles, in radians in the range
// (-pi, pi], of a complex spectrum.
func Phase(X []complex128) []float64 {
	out := make([]float64, len(X))
	for i, v := range X {
		out[i] = cmplx.Phase(v)
	}
	return out
}

// PowerSpectrum returns the element-wise squared magnitudes |X[k]|^2 of a
// complex spectrum.
func PowerSpectrum(X []complex128) []float64 {
	out := make([]float64, len(X))
	for i, v := range X {
		m := cmplx.Abs(v)
		out[i] = m * m
	}
	return out
}

// Goertzel computes a single DFT bin X[k] of the length-n signal x using the
// Goertzel algorithm, which is more efficient than a full DFT when only a few
// bins are required. The result equals DFT(x)[k].
func Goertzel(x []float64, k, n int) complex128 {
	w := 2 * math.Pi * float64(k) / float64(n)
	coeff := 2 * math.Cos(w)
	var s0, s1, s2 float64
	for i := 0; i < n; i++ {
		s0 = x[i] + coeff*s1 - s2
		s2 = s1
		s1 = s0
	}
	// Combine the final states and apply the group-delay phase correction so
	// the result matches DFT(x)[k] exactly.
	y := complex(s1, 0) - complex(s2, 0)*cmplx.Rect(1, -w)
	return y * cmplx.Rect(1, -w*float64(n-1))
}

// GoertzelPower returns the power (squared magnitude) of the signal x at the
// frequency targetFreq given the sampleRate, using the Goertzel recurrence.
// Unlike [Goertzel] the target frequency need not fall on an exact DFT bin.
func GoertzelPower(x []float64, targetFreq, sampleRate float64) float64 {
	n := len(x)
	w := 2 * math.Pi * targetFreq / sampleRate
	coeff := 2 * math.Cos(w)
	var s1, s2 float64
	for i := 0; i < n; i++ {
		s0 := x[i] + coeff*s1 - s2
		s2 = s1
		s1 = s0
	}
	return s1*s1 + s2*s2 - coeff*s1*s2
}

// FFTPlan holds precomputed twiddle factors for repeated radix-2 transforms of
// a fixed power-of-two length. Reusing a plan avoids recomputing the twiddle
// table on every call and is convenient when transforming many blocks of the
// same size. A plan is safe for concurrent use because Forward and Inverse do
// not mutate it.
type FFTPlan struct {
	n  int
	tw []complex128 // tw[j] = exp(-2*pi*i*j/n), j = 0..n/2-1
}

// NewFFTPlan creates a reusable transform plan for signals of length n, which
// must be a power of two. NewFFTPlan panics for other lengths.
func NewFFTPlan(n int) *FFTPlan {
	if !IsPow2(n) {
		panic("transform: NewFFTPlan requires a power-of-two length")
	}
	tw := make([]complex128, n/2)
	for j := range tw {
		tw[j] = cmplx.Rect(1, -2*math.Pi*float64(j)/float64(n))
	}
	return &FFTPlan{n: n, tw: tw}
}

// Len returns the fixed transform length the plan was created for.
func (p *FFTPlan) Len() int { return p.n }

func (p *FFTPlan) transform(x []complex128, inverse bool) {
	n := p.n
	for i, j := 1, 0; i < n; i++ {
		bit := n >> 1
		for ; j&bit != 0; bit >>= 1 {
			j ^= bit
		}
		j ^= bit
		if i < j {
			x[i], x[j] = x[j], x[i]
		}
	}
	for length := 2; length <= n; length <<= 1 {
		half := length >> 1
		step := n / length
		for i := 0; i < n; i += length {
			for k := 0; k < half; k++ {
				tw := p.tw[k*step]
				if inverse {
					tw = cmplx.Conj(tw)
				}
				u := x[i+k]
				v := x[i+k+half] * tw
				x[i+k] = u + v
				x[i+k+half] = u - v
			}
		}
	}
}

// Forward computes the discrete Fourier transform of x using the plan. The
// length of x must equal the plan length. The input is not modified.
func (p *FFTPlan) Forward(x []complex128) []complex128 {
	if len(x) != p.n {
		panic("transform: FFTPlan.Forward length mismatch")
	}
	out := make([]complex128, p.n)
	copy(out, x)
	p.transform(out, false)
	return out
}

// Inverse computes the inverse discrete Fourier transform of X using the plan,
// including the 1/N scaling. The length of X must equal the plan length. The
// input is not modified.
func (p *FFTPlan) Inverse(X []complex128) []complex128 {
	if len(X) != p.n {
		panic("transform: FFTPlan.Inverse length mismatch")
	}
	out := make([]complex128, p.n)
	copy(out, X)
	p.transform(out, true)
	scale := complex(1/float64(p.n), 0)
	for i := range out {
		out[i] *= scale
	}
	return out
}
