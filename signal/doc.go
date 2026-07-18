// Package signal is a small, dependency-free digital-signal-processing
// toolkit written entirely with the Go standard library (only the math
// package). It operates on plain []float64 sample vectors and []complex128
// spectra and is completely standalone: it does not import the parent algebra
// package.
//
// The package is organised into the following groups.
//
// # Window functions
//
// Fixed-length analysis windows: [Rectangular], [Hann], [Hamming],
// [Blackman], [BlackmanHarris], [Nuttall], [FlatTop], [Bartlett], [Welch],
// [Cosine], [Gaussian], [Kaiser] and [Tukey], together with the modified
// Bessel function [BesselI0], the [KaiserBeta] design helper and the
// element-wise [ApplyWindow].
//
// # FIR filter design
//
// Linear-phase windowed-sinc finite-impulse-response filters:
// [FIRLowpass], [FIRHighpass], [FIRBandpass] and [FIRBandstop] (each using a
// Hamming window by default) plus the explicit-window variants
// [FIRLowpassWin], [FIRHighpassWin], [FIRBandpassWin] and [FIRBandstopWin].
// Supporting routines are the normalized [Sinc], the reusable [FIRFilter]
// state object, [ApplyFIR], [FIRFrequencyResponse] and [FIRGroupDelay].
//
// # IIR (biquad) filter design
//
// Second-order recursive sections following the Robert Bristow-Johnson audio
// EQ cookbook: [BiquadLowpass], [BiquadHighpass], [BiquadBandpass],
// [BiquadNotch], [BiquadAllpass], [BiquadPeaking], [BiquadLowShelf] and
// [BiquadHighShelf]. Higher-order maximally-flat filters are built as cascades
// of second-order sections by [ButterworthLowpass] and [ButterworthHighpass]
// and applied with [FilterSOS]. The [Biquad] type carries the coefficients and
// per-sample state.
//
// # Convolution and correlation
//
// [Convolve], [ConvolveSame], [ConvolveValid], [CrossCorrelate] and
// [AutoCorrelate].
//
// # Resampling
//
// [Upsample], [Downsample], [Decimate], [Interpolate] and [ResampleLinear].
//
// # Spectral analysis
//
// The direct [DFT] and [IDFT], the spectral helpers [Magnitude], [Phase] and
// [PowerSpectrum], the [Periodogram] and Welch-averaged [WelchPSD] power
// spectral density estimators and the [FrequencyBins] helper.
//
// # Time-domain analysis
//
// [MovingAverage], [MovingAverageCentered], [ExponentialMovingAverage],
// [CumulativeSum], [Diff], [RMS], [Energy] and [ZeroPad].
//
// Frequencies are expressed throughout as fractions of the Nyquist rate:
// a normalized cutoff of 1.0 corresponds to half the sampling rate. Functions
// that take an explicit sampling rate fs interpret it in hertz.
package signal
