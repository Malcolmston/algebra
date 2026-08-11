# signal — small, dependency-free digital-signal-processing toolkit written entirely with the Go

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/signal.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/signal)

Package signal is a small, dependency-free digital-signal-processing toolkit
written entirely with the Go standard library (only the math package). It
operates on plain []float64 sample vectors and []complex128 spectra and is
completely standalone: it does not import the parent algebra package.

The package is organised into the following groups.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/signal"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func ApplyFIR(taps, x []float64) []float64` | ApplyFIR filters the signal x through the FIR filter defined by taps and returns the full linear convolution, of length len(x)+len(taps)-1. |
| `func ApplyWindow(x, w []float64) []float64` | ApplyWindow returns the element-wise product of the signal x and the window w. |
| `func AutoCorrelate(a []float64) []float64` | AutoCorrelate returns the full autocorrelation of a, equal to CrossCorrelate(a, a). |
| `func Bartlett(n int) []float64` | Bartlett returns the length-n Bartlett (triangular) window with zero-valued endpoints, w`k` = 1 - \|(k - (n-1)/2)/((n-1)/2)\|. |
| `func BesselI0(x float64) float64` | BesselI0 returns the modified Bessel function of the first kind of order zero, I0(x), evaluated by its always-convergent power series Σ… |
| `func Blackman(n int) []float64` | Blackman returns the length-n Blackman window, w`k` = 0.42 - 0.5·cos(2π·k/(n-1)) + 0.08·cos(4π·k/(n-1)). |
| `func BlackmanHarris(n int) []float64` | BlackmanHarris returns the length-n four-term Blackman-Harris window using the coefficients 0.35875, 0.48829, 0.14128 and 0.01168. |
| `func Convolve(a, b []float64) []float64` | Convolve returns the full discrete linear convolution of a and b, (a∗b)`k` = Σ a`i`·b`k-i`, as a new slice of length len(a)+len(b)-1. |
| `func ConvolveSame(a, b []float64) []float64` | ConvolveSame returns the central part of the full convolution of a and b that is the same length as a. |
| `func ConvolveValid(a, b []float64) []float64` | ConvolveValid returns only those elements of the convolution of a and b that are computed without any zero-padding of the shorter sequence, i.e. |
| `func Cosine(n int) []float64` | Cosine returns the length-n cosine (sine) window, w`k` = sin(π·k/(n-1)). |
| `func CrossCorrelate(a, b []float64) []float64` | CrossCorrelate returns the full cross-correlation of a and b as a slice of length len(a)+len(b)-1. |
| `func CumulativeSum(x []float64) []float64` | CumulativeSum returns the running (prefix) sum of x, where out`i` is the sum of x[0…i]. |
| `func DFT(x []float64) []complex128` | DFT returns the discrete Fourier transform of the real signal x, computed directly as X`k` = Σ x`n`·exp(-j2πkn/N) for k = 0…N-1. |
| `func Decimate(x []float64, factor int) []float64` | Decimate reduces the sampling rate of x by an integer factor after applying a linear-phase FIR low-pass anti-aliasing filter with cutoff 1/factor of… |
| `func Diff(x []float64) []float64` | Diff returns the first-order forward difference of x, out`i` = x[i+1] - x`i`. |
| `func Downsample(x []float64, factor int) []float64` | Downsample decreases the sampling rate of x by an integer factor by keeping every factor-th sample starting from index 0 (decimation without… |
| `func Energy(x []float64) float64` | Energy returns the total energy of x, Σ x`i`². |
| `func ExponentialMovingAverage(x []float64, alpha float64) []float64` | ExponentialMovingAverage returns the exponentially-weighted moving average of x with smoothing factor alpha in [0, 1]: out[0] = x[0] and out`i` =… |
| `func FIRBandpass(numtaps int, low, high float64) []float64` | FIRBandpass designs a length-numtaps linear-phase band-pass FIR filter passing normalized frequencies between low and high using a Hamming window. |
| `func FIRBandpassWin(numtaps int, low, high float64, w []float64) []float64` | FIRBandpassWin designs a length-numtaps linear-phase band-pass FIR filter that passes normalized frequencies between low and high (both in fractions… |
| `func FIRBandstop(numtaps int, low, high float64) []float64` | FIRBandstop designs a length-numtaps linear-phase band-stop FIR filter rejecting normalized frequencies between low and high using a Hamming window. |
| `func FIRBandstopWin(numtaps int, low, high float64, w []float64) []float64` | FIRBandstopWin designs a length-numtaps linear-phase band-stop (notch) FIR filter that rejects normalized frequencies between low and high using the… |
| `func FIRFrequencyResponse(taps, freqs []float64) []complex128` | FIRFrequencyResponse evaluates the complex frequency response H(ω) of the FIR filter taps at each normalized frequency in freqs, where a value of… |
| `func FIRGroupDelay(taps, freqs []float64) []float64` | FIRGroupDelay returns the group delay, in samples, of the FIR filter taps at each normalized frequency in freqs (1.0 corresponds to Nyquist). |
| `func FIRHighpass(numtaps int, cutoff float64) []float64` | FIRHighpass designs a length-numtaps linear-phase high-pass FIR filter with the given normalized cutoff using a Hamming window. |
| `func FIRHighpassWin(numtaps int, cutoff float64, w []float64) []float64` | FIRHighpassWin designs a length-numtaps linear-phase high-pass FIR filter with the given normalized cutoff using the supplied window w. |
| `func FIRLowpass(numtaps int, cutoff float64) []float64` | FIRLowpass designs a length-numtaps linear-phase low-pass FIR filter with the given normalized cutoff (in fractions of Nyquist) using a Hamming… |
| `func FIRLowpassWin(numtaps int, cutoff float64, w []float64) []float64` | FIRLowpassWin designs a length-numtaps linear-phase low-pass FIR filter with the given normalized cutoff (in fractions of Nyquist, 0 < cutoff < 1)… |
| `func FilterSOS(sos []Biquad, x []float64) []float64` | FilterSOS filters the signal x through a cascade of second-order sections sos and returns a new slice of the same length. |
| `func FlatTop(n int) []float64` | FlatTop returns the length-n flat-top window using the SR785 five-term coefficients. |
| `func FrequencyBins(n int, fs float64) []float64` | FrequencyBins returns the n bin-centre frequencies of an n-point DFT for a signal sampled at fs hertz, i.e. |
| `func Gaussian(n int, sigma float64) []float64` | Gaussian returns the length-n Gaussian window with normalized standard deviation sigma (0 < sigma <= 0.5 is typical), w`k` = exp(-½·((k -… |
| `func Hamming(n int) []float64` | Hamming returns the length-n Hamming window, w`k` = 0.54 - 0.46·cos(2π·k/(n-1)). |
| `func Hann(n int) []float64` | Hann returns the length-n Hann (raised-cosine) window, w`k` = 0.5·(1 - cos(2π·k/(n-1))). |
| `func IDFT(X []complex128) []complex128` | IDFT returns the inverse discrete Fourier transform of the spectrum X, x`n` = (1/N)·Σ X`k`·exp(+j2πkn/N), as a complex slice of the same length. |
| `func Interpolate(x []float64, factor int) []float64` | Interpolate increases the sampling rate of x by an integer factor, zero-stuffing and then applying a linear-phase FIR low-pass anti-imaging filter of… |
| `func Kaiser(n int, beta float64) []float64` | Kaiser returns the length-n Kaiser window with shape parameter beta, w`k` = I0(beta·√(1 - (2k/(n-1) - 1)²)) / I0(beta), where I0 is `BesselI0`. |
| `func KaiserBeta(atten float64) float64` | KaiserBeta returns the Kaiser shape parameter beta that achieves a stop-band attenuation of atten decibels, using Kaiser's empirical formula. |
| `func Magnitude(X []complex128) []float64` | Magnitude returns the element-wise magnitude \|X`k`\| of a complex spectrum X. The result has the same length as X. |
| `func MovingAverage(x []float64, w int) []float64` | MovingAverage returns the causal (trailing) simple moving average of x with window length w: out`i` is the mean of x[max(0,i-w+1) … i]. |
| `func MovingAverageCentered(x []float64, w int) []float64` | MovingAverageCentered returns the centered simple moving average of x with window length w. |
| `func Nuttall(n int) []float64` | Nuttall returns the length-n four-term Nuttall window using the coefficients 0.355768, 0.487396, 0.144232 and 0.012604. |
| `func Periodogram(x []float64, window []float64, fs float64) (freqs, psd []float64)` | Periodogram estimates the one-sided power spectral density of the real signal x sampled at fs hertz. |
| `func Phase(X []complex128) []float64` | Phase returns the element-wise phase angle arg(X`k`), in radians in the range (-π, π], of a complex spectrum X. The result has the same length as X. |
| `func PowerSpectrum(X []complex128) []float64` | PowerSpectrum returns the element-wise squared magnitude \|X`k`\|² of a complex spectrum X. The result has the same length as X. |
| `func RMS(x []float64) float64` | RMS returns the root-mean-square value of x, √(Σ x`i`²/N). |
| `func Rectangular(n int) []float64` | Rectangular returns the length-n rectangular (boxcar) window, whose samples are all 1. |
| `func ResampleLinear(x []float64, numOut int) []float64` | ResampleLinear resamples x to exactly numOut samples using linear interpolation over the original sample positions. |
| `func Sinc(x float64) float64` | Sinc returns the normalized cardinal sine, sinc(x) = sin(πx)/(πx), with the removable singularity sinc(0) = 1 handled exactly. |
| `func Tukey(n int, alpha float64) []float64` | Tukey returns the length-n Tukey (tapered-cosine) window with taper fraction alpha in [0, 1]. |
| `func Upsample(x []float64, factor int) []float64` | Upsample increases the sampling rate of x by an integer factor by inserting factor-1 zeros between each input sample (zero stuffing). |
| `func Welch(n int) []float64` | Welch returns the length-n Welch (parabolic) window, w`k` = 1 - ((k - (n-1)/2)/((n-1)/2))². |
| `func WelchPSD(x []float64, segLen, overlap int, fs float64) (freqs, psd []float64)` | WelchPSD estimates the one-sided power spectral density of x sampled at fs hertz using Welch's method: x is split into overlapping segments of length… |
| `func ZeroPad(x []float64, n int) []float64` | ZeroPad returns a copy of x extended with trailing zeros to total length n. |

### Types

| Type | What it is |
| --- | --- |
| `Biquad` | Biquad is a second-order recursive (IIR) filter section realising the transfer function |
| `FIRFilter` | FIRFilter is a reusable streaming finite-impulse-response filter. |

<details>
<summary><code>Biquad</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BiquadAllpass(f0, fs, q float64) Biquad` | BiquadAllpass returns a second-order all-pass `Biquad` centred on f0 hertz with sampling rate fs hertz and quality factor q, following the RBJ… |
| `func BiquadBandpass(f0, fs, q float64) Biquad` | BiquadBandpass returns a second-order band-pass `Biquad` centred on f0 hertz with sampling rate fs hertz and quality factor q, normalized to unit (0… |
| `func BiquadHighShelf(f0, fs, q, gainDB float64) Biquad` | BiquadHighShelf returns a second-order high-shelf `Biquad` with shelf midpoint f0 hertz, sampling rate fs hertz, quality factor q and shelf gain… |
| `func BiquadHighpass(f0, fs, q float64) Biquad` | BiquadHighpass returns a second-order high-pass `Biquad` with corner frequency f0 hertz, sampling rate fs hertz and quality factor q, following the… |
| `func BiquadLowShelf(f0, fs, q, gainDB float64) Biquad` | BiquadLowShelf returns a second-order low-shelf `Biquad` with shelf midpoint f0 hertz, sampling rate fs hertz, quality factor q and shelf gain gainDB… |
| `func BiquadLowpass(f0, fs, q float64) Biquad` | BiquadLowpass returns a second-order low-pass `Biquad` with corner frequency f0 hertz, sampling rate fs hertz and quality factor q, following the RBJ… |
| `func BiquadNotch(f0, fs, q float64) Biquad` | BiquadNotch returns a second-order band-reject (notch) `Biquad` centred on f0 hertz with sampling rate fs hertz and quality factor q, following the… |
| `func BiquadPeaking(f0, fs, q, gainDB float64) Biquad` | BiquadPeaking returns a second-order peaking-EQ `Biquad` centred on f0 hertz with sampling rate fs hertz, quality factor q and a boost or cut of… |
| `func ButterworthHighpass(order int, f0, fs float64) []Biquad` | ButterworthHighpass designs an order-order Butterworth high-pass filter with corner frequency f0 hertz and sampling rate fs hertz, returned as a… |
| `func ButterworthLowpass(order int, f0, fs float64) []Biquad` | ButterworthLowpass designs an order-order Butterworth low-pass filter with corner frequency f0 hertz and sampling rate fs hertz, returned as a… |
| `func (bq *Biquad) FrequencyResponse(freqs []float64) []complex128` | FrequencyResponse evaluates the section's complex frequency response H(ω) at each normalized frequency in freqs, where 1.0 corresponds to Nyquist… |
| `func (bq *Biquad) Process(x float64) float64` | Process filters one input sample x through the section and returns the output sample, advancing the internal state. |
| `func (bq *Biquad) ProcessBlock(x []float64) []float64` | ProcessBlock filters an entire block x through the section, returning a new slice of the same length, with state carried across calls. |
| `func (bq *Biquad) Reset()` | Reset clears the section's internal state variables so that subsequent filtering starts from rest. |

</details>

<details>
<summary><code>FIRFilter</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewFIRFilter(taps []float64) *FIRFilter` | NewFIRFilter returns a streaming `FIRFilter` initialized with a copy of taps and a zeroed delay line. |
| `func (f *FIRFilter) Process(x []float64) []float64` | Process filters an entire block x through the streaming filter, returning a new slice of the same length. |
| `func (f *FIRFilter) ProcessSample(x float64) float64` | ProcessSample feeds a single input sample through the filter and returns the corresponding output sample, updating the internal delay line. |
| `func (f *FIRFilter) Reset()` | Reset clears the filter's delay line, returning it to the state of a freshly constructed filter without reallocating. |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/signal).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
