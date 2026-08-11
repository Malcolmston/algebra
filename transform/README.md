# transform — implements integral and discrete signal transforms in pure Go using only the standard

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/transform.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/transform)

Package transform implements integral and discrete signal transforms in pure
Go using only the standard library.

The package collects the transforms most often needed in numerical analysis,
digital-signal processing and applied mathematics and provides correct,
self-contained implementations of each.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/transform"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func ApplyWindow(x, w []float64) []float64` | ApplyWindow returns the element-wise product of the signal x and the window w. |
| `func AutoCorrelate(x []float64) []float64` | AutoCorrelate returns the full autocorrelation of x, i.e. |
| `func Bartlett(n int) []float64` | Bartlett returns the n-point Bartlett (triangular) window, which tapers linearly to zero at both ends. |
| `func Blackman(n int) []float64` | Blackman returns the n-point Blackman window, w`i` = 0.42 - 0.5*cos(2*pi*i/(n-1)) + 0.08*cos(4*pi*i/(n-1)). |
| `func Bluestein(x []complex128) []complex128` | Bluestein computes the discrete Fourier transform of x for an arbitrary length by expressing it as a convolution (the chirp-z transform with a = 1, w… |
| `func ChirpZTransform(x []complex128, m int, w, a complex128) []complex128` | ChirpZTransform computes the chirp-z transform of x, sampling the Z-transform along a spiral contour. |
| `func CircularConvolve(a, b []float64) []float64` | CircularConvolve returns the circular (cyclic) convolution of a and b. |
| `func Convolve(a, b []float64) []float64` | Convolve returns the full linear convolution of the real sequences a and b, computed directly. |
| `func ConvolveComplex(a, b []complex128) []complex128` | ConvolveComplex returns the full linear convolution of the complex sequences a and b, computed directly. |
| `func ConvolveFFT(a, b []float64) []float64` | ConvolveFFT returns the full linear convolution of the real sequences a and b using the FFT (multiplication in the frequency domain). |
| `func Correlate(a, b []float64) []float64` | Correlate returns the full cross-correlation of a and b, defined as the convolution of a with the time-reversed b. |
| `func CrossCorrelateFFT(a, b []float64) []float64` | CrossCorrelateFFT returns the full cross-correlation of a and b computed via the FFT. |
| `func DCT(x []float64) []float64` | DCT computes the unnormalized type-II discrete cosine transform of x: |
| `func DCT1(x []float64) []float64` | DCT1 computes the unnormalized type-I discrete cosine transform of x, which requires at least two samples: |
| `func DCT4(x []float64) []float64` | DCT4 computes the unnormalized type-IV discrete cosine transform of x: |
| `func DFT(x []complex128) []complex128` | DFT computes the discrete Fourier transform of x directly in O(n^2) time, returning X`k` = sum_n x`n` exp(-2*pi*i*k*n/N). |
| `func DFTMatrix(n int) [][]complex128` | DFTMatrix returns the n-by-n discrete Fourier transform matrix W whose entry W`j``k` equals exp(-2*pi*i*j*k/n). |
| `func DST(x []float64) []float64` | DST computes the unnormalized type-I discrete sine transform of x: |
| `func DTFT(x []float64, omega float64) complex128` | DTFT evaluates the discrete-time Fourier transform of the real sequence x at the angular frequency omega (radians/sample), returning sum_{n=0}^{N-1}… |
| `func DTFTComplex(x []complex128, omega float64) complex128` | DTFTComplex evaluates the discrete-time Fourier transform of the complex sequence x at the angular frequency omega (radians/sample). |
| `func DTFTSample(x []float64, omegas []float64) []complex128` | DTFTSample evaluates the discrete-time Fourier transform of the real sequence x at each of the supplied angular frequencies, returning one complex… |
| `func Envelope(x []float64) []float64` | Envelope returns the amplitude envelope of the real sequence x, i.e. |
| `func FFT(x []complex128) []complex128` | FFT computes the discrete Fourier transform of x using the radix-2 Cooley-Tukey algorithm. |
| `func FFT2D(m [][]complex128) [][]complex128` | FFT2D computes the two-dimensional discrete Fourier transform of a matrix by transforming every row and then every column. |
| `func FFTAny(x []complex128) []complex128` | FFTAny computes the discrete Fourier transform of x for any length. |
| `func FFTFreq(n int, d float64) []float64` | FFTFreq returns the n sample frequencies corresponding to the bins of an n-point DFT taken with sample spacing d (in seconds). |
| `func FFTReal(x []float64) []complex128` | FFTReal computes the discrete Fourier transform of a real-valued signal x, returning the full complex spectrum of the same length. |
| `func FFTShift(x []complex128) []complex128` | FFTShift rearranges a spectrum so that the zero-frequency component is moved to the center of the slice, matching the conventional plotting order. |
| `func Goertzel(x []float64, k, n int) complex128` | Goertzel computes a single DFT bin X`k` of the length-n signal x using the Goertzel algorithm, which is more efficient than a full DFT when only a… |
| `func GoertzelPower(x []float64, targetFreq, sampleRate float64) float64` | GoertzelPower returns the power (squared magnitude) of the signal x at the frequency targetFreq given the sampleRate, using the Goertzel recurrence. |
| `func Hamming(n int) []float64` | Hamming returns the n-point Hamming window, w`i` = 0.54 - 0.46*cos(2*pi*i/(n-1)). |
| `func Hann(n int) []float64` | Hann returns the n-point Hann (raised-cosine) window, w`i` = 0.5*(1 - cos(2*pi*i/(n-1))). |
| `func Hilbert(x []float64) []complex128` | Hilbert returns the analytic signal of the real sequence x, whose real part is x and whose imaginary part is the Hilbert transform of x. |
| `func HilbertTransform(x []float64) []float64` | HilbertTransform returns the Hilbert transform of the real sequence x, i.e. |
| `func IDCT(X []float64) []float64` | IDCT computes the inverse of the type-II transform produced by `DCT` (it is the scaled type-III transform): |
| `func IDFT(X []complex128) []complex128` | IDFT computes the inverse discrete Fourier transform of X directly, returning x`n` = (1/N) sum_k X`k` exp(+2*pi*i*k*n/N). |
| `func IDST(X []float64) []float64` | IDST computes the inverse of the type-I sine transform produced by `DST`. |
| `func IFFT(X []complex128) []complex128` | IFFT computes the inverse radix-2 FFT of X, including the 1/N scaling. |
| `func IFFT2D(m [][]complex128) [][]complex128` | IFFT2D computes the inverse two-dimensional discrete Fourier transform, including the 1/(rows*cols) scaling. |
| `func IFFTAny(X []complex128) []complex128` | IFFTAny computes the inverse discrete Fourier transform of X for any length, including the 1/N scaling. |
| `func IFFTShift(x []complex128) []complex128` | IFFTShift is the inverse of `FFTShift`; it moves the zero-frequency component from the center back to the start of the slice. |
| `func IRFFT(X []complex128, n int) []float64` | IRFFT reconstructs a real signal of length n from the floor(n/2)+1 half-spectrum X produced by `RFFT`. |
| `func InstantaneousFrequency(x []float64, sampleRate float64) []float64` | InstantaneousFrequency returns the instantaneous frequency, in the same units as sampleRate (typically Hz), of the real sequence x. |
| `func InstantaneousPhase(x []float64) []float64` | InstantaneousPhase returns the unwrapped instantaneous phase, in radians, of the real sequence x, obtained as the argument of its analytic signal… |
| `func InverseLaplaceEuler(F func(complex128) complex128, t float64, m int) float64` | InverseLaplaceEuler approximates the inverse Laplace transform f(t) using the Fourier-series method with Euler summation (Abate and Whitt). |
| `func InverseLaplaceStehfest(F func(float64) float64, t float64, n int) float64` | InverseLaplaceStehfest approximates the inverse Laplace transform f(t) using the Gaver-Stehfest algorithm, which evaluates the real-argument… |
| `func InverseLaplaceTalbot(F func(complex128) complex128, t float64, m int) float64` | InverseLaplaceTalbot approximates the inverse Laplace transform f(t) of the transform F, supplied as a callable, using the fixed Talbot method of… |
| `func InverseZTransform(X func(complex128) complex128, n int, radius float64, m int) []complex128` | InverseZTransform recovers the first n samples of a causal sequence from its Z-transform X(z), which is supplied as a callable. |
| `func IsPow2(n int) bool` | IsPow2 reports whether n is a positive power of two. |
| `func Laplace(f func(float64) float64, s complex128, upper float64, n int) complex128` | Laplace numerically approximates the Laplace transform |
| `func Magnitude(X []complex128) []float64` | Magnitude returns the element-wise magnitudes (absolute values) of a complex spectrum. |
| `func NextPow2(n int) int` | NextPow2 returns the smallest power of two that is greater than or equal to n. |
| `func Periodogram(x []float64) []float64` | Periodogram returns the periodogram power-spectral-density estimate of the real signal x, defined as (1/N) \|DFT(x)`k`\|^2 for k = 0 .. |
| `func Phase(X []complex128) []float64` | Phase returns the element-wise phase angles, in radians in the range (-pi, pi], of a complex spectrum. |
| `func PhaseUnwrap(phase []float64) []float64` | PhaseUnwrap returns a copy of the phase sequence with discontinuities larger than pi removed by adding integer multiples of 2*pi, producing a… |
| `func PowerSpectrum(X []complex128) []float64` | PowerSpectrum returns the element-wise squared magnitudes \|X`k`\|^2 of a complex spectrum. |
| `func RFFT(x []float64) []complex128` | RFFT computes the discrete Fourier transform of a real-valued signal x and returns only the non-redundant first floor(N/2)+1 frequency bins. |
| `func RFFTFreq(n int, d float64) []float64` | RFFTFreq returns the floor(n/2)+1 non-negative sample frequencies corresponding to the bins produced by `RFFT` for an n-point signal with sample… |
| `func SampleDTFT(x []float64, m int) []complex128` | SampleDTFT samples the discrete-time Fourier transform of the real sequence x at m frequencies equally spaced over [0, 2*pi), returning |
| `func StehfestCoefficients(n int) []float64` | StehfestCoefficients returns the n Gaver-Stehfest weights V_1..V_n (indexed from zero) used by `InverseLaplaceStehfest`. |
| `func Welch(n int) []float64` | Welch returns the n-point Welch (parabolic) window, w`i` = 1 - ((i-(n-1)/2)/((n-1)/2))^2. |
| `func ZTransform(x []float64, z complex128) complex128` | ZTransform evaluates the one-sided Z-transform of the real sequence x at the complex point z, returning sum_{n=0}^{N-1} x`n` z^{-n}. |
| `func ZTransformComplex(x []complex128, z complex128) complex128` | ZTransformComplex evaluates the one-sided Z-transform of the complex sequence x at the point z, returning sum_{n=0}^{N-1} x`n` z^{-n}. |
| `func ZeroPad(x []complex128, n int) []complex128` | ZeroPad returns a copy of x extended with zeros to length n. |
| `func ZeroPadReal(x []float64, n int) []float64` | ZeroPadReal returns a copy of the real slice x extended with zeros to length n. |

### Types

| Type | What it is |
| --- | --- |
| `FFTPlan` | FFTPlan holds precomputed twiddle factors for repeated radix-2 transforms of a fixed power-of-two length. |

<details>
<summary><code>FFTPlan</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewFFTPlan(n int) *FFTPlan` | NewFFTPlan creates a reusable transform plan for signals of length n, which must be a power of two. |
| `func (p *FFTPlan) Forward(x []complex128) []complex128` | Forward computes the discrete Fourier transform of x using the plan. |
| `func (p *FFTPlan) Inverse(X []complex128) []complex128` | Inverse computes the inverse discrete Fourier transform of X using the plan, including the 1/N scaling. |
| `func (p *FFTPlan) Len() int` | Len returns the fixed transform length the plan was created for. |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/transform).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
