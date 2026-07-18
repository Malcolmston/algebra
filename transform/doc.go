// Package transform implements integral and discrete signal transforms in
// pure Go using only the standard library.
//
// The package collects the transforms most often needed in numerical
// analysis, digital-signal processing and applied mathematics and provides
// correct, self-contained implementations of each.
//
// The following families are covered:
//
//   - Discrete Fourier transforms: the direct O(n^2) [DFT] and [IDFT], the
//     radix-2 Cooley-Tukey [FFT] and [IFFT], the mixed-length [FFTAny] and
//     [IFFTAny] (which fall back to Bluestein's algorithm for non-power-of-two
//     lengths), the real-input transforms [RFFT] and [IRFFT], two-dimensional
//     [FFT2D]/[IFFT2D], the reusable [FFTPlan], and the single-bin [Goertzel]
//     algorithm;
//   - spectral helpers: [FFTShift], [FFTFreq], [Magnitude], [Phase],
//     [PowerSpectrum], [Periodogram] and the analysis windows [Hann],
//     [Hamming], [Blackman], [Bartlett] and [Welch];
//   - the discrete cosine and sine transforms [DCT], [IDCT], [DCT1], [DCT4],
//     [DST] and [IDST];
//   - fast linear, circular and cross correlation via [Convolve],
//     [ConvolveFFT], [CircularConvolve], [Correlate] and [AutoCorrelate];
//   - the numerical Laplace transform [Laplace] together with three classical
//     inversion methods: the fixed Talbot method [InverseLaplaceTalbot], the
//     Gaver-Stehfest method [InverseLaplaceStehfest] and the Fourier-series
//     Euler method [InverseLaplaceEuler];
//   - the [ZTransform], its numerical inverse [InverseZTransform], the
//     chirp-z transform [ChirpZTransform] and [Bluestein]'s algorithm;
//   - the [Hilbert] analytic signal together with the derived [Envelope],
//     [InstantaneousPhase] and [InstantaneousFrequency], and the
//     discrete-time Fourier transform samplers [DTFT] and [SampleDTFT].
//
// Every routine is deterministic and depends only on the Go standard library.
package transform
