# wavelet — provides discrete wavelet transforms and related multiresolution tools implemented

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/wavelet.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/wavelet)

Package wavelet provides discrete wavelet transforms and related
multiresolution tools implemented with the Go standard library only.

- Wavelet families (see wavelet.go): the orthonormal `Haar` wavelet and the
Daubechies family via `Daubechies`, with the common shortcuts `DB2` and `DB4`.
A `Wavelet` exposes its analysis and synthesis filter banks through
`Wavelet.DecLo`, `Wavelet.DecHi`, `Wavelet.RecLo` and `Wavelet.RecHi`. -
Single-level transforms (see dwt.go): the forward `DWT` and its inverse
`IDWT`, together with the `Coefficients` container. - Multiresolution analysis
(see mra.go): multi-level decomposition `WaveDec` and reconstruction
`WaveRec`, the `Decomposition` type, and projection of each scale back to the
original resolution via `Decomposition.MRAComponents`. - Wavelet packets (see
packet.go): the full binary `PacketTree` produced by `WaveletPacket`, with
perfect reconstruction and the Coifman-Wickerhauser cost measures
`ShannonEntropy` and `LogEnergy`. - Thresholding and denoising (see
threshold.go): `SoftThreshold`, `HardThreshold`, the VisuShrink
`UniversalThreshold`, robust noise estimation via `EstimateNoiseSigma`, and
the end-to-end `Denoise`. - Two-dimensional transforms (see dwt2d.go): a
single level `DWT2D` and its inverse `IDWT2D` returning the `DWT2DResult`
subbands. - Signal utilities (see signal.go): boundary extensions, linear
`Convolve`, `Downsample`/`Upsample`, and small numeric helpers.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/wavelet"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func Convolve(x, h []float64) []float64` | Convolve returns the full linear convolution of x and h, a slice of length len(x)+len(h)-1. |
| `func DWT(signal []float64, w Wavelet) (approx, detail []float64)` | DWT computes one level of the forward discrete wavelet transform of signal using the wavelet w with periodic boundary handling. |
| `func Denoise(signal []float64, w Wavelet, levels int, mode ThresholdMode) []float64` | Denoise removes noise from signal by wavelet shrinkage (the VisuShrink procedure). |
| `func Downsample(x []float64, factor int) []float64` | Downsample returns every factor-th sample of x starting at index 0. |
| `func Energy(x []float64) float64` | Energy returns the sum of squares of the elements of x. |
| `func EstimateNoiseSigma(detail []float64) float64` | EstimateNoiseSigma estimates the standard deviation of white Gaussian noise from a set of finest-scale detail coefficients using the robust median… |
| `func HardThreshold(x, lambda float64) float64` | HardThreshold applies the hard-thresholding rule to x with threshold lambda: x is kept when \|x\| > lambda and set to zero otherwise. |
| `func IDWT(approx, detail []float64, w Wavelet) []float64` | IDWT reconstructs a signal from its approximation and detail coefficients using the wavelet w, inverting a single level of `DWT`. |
| `func IDWT2D(r *DWT2DResult, w Wavelet) [][]float64` | IDWT2D reconstructs an image from its four subbands using wavelet w, inverting a single level of `DWT2D` to within floating-point rounding. |
| `func IsPowerOfTwo(n int) bool` | IsPowerOfTwo reports whether n is a positive power of two. |
| `func L2Norm(x []float64) float64` | L2Norm returns the Euclidean (l2) norm of x, sqrt(sum x_i^2). |
| `func LogEnergy(x []float64) float64` | LogEnergy returns the log-energy cost measure sum_i ln(x_i^2) taken over the non-zero coefficients of x, another additive Coifman-Wickerhauser cost… |
| `func MaxAbs(x []float64) float64` | MaxAbs returns the largest absolute value among the elements of x, or 0 for an empty slice. |
| `func MaxLevel(n int) int` | MaxLevel returns the maximum number of dyadic decomposition levels that a signal of length n admits under periodic boundary handling, which is the… |
| `func MeanAbsoluteDeviation(x []float64) float64` | MeanAbsoluteDeviation returns the median absolute deviation of x about its own median, median(\|x_i - median(x)\|). |
| `func Median(x []float64) float64` | Median returns the median of x. |
| `func NextPowerOfTwo(n int) int` | NextPowerOfTwo returns the smallest power of two that is greater than or equal to n. |
| `func PeriodicExtend(x []float64, pad int) []float64` | PeriodicExtend returns a copy of x extended by pad samples on each side using periodic (wrap-around) boundary handling, so the result has length… |
| `func ShannonEntropy(x []float64) float64` | ShannonEntropy returns the normalized Shannon entropy of x, a Coifman-Wickerhauser cost measure of coefficient concentration. |
| `func SoftThreshold(x, lambda float64) float64` | SoftThreshold applies the soft-thresholding (shrinkage) rule to x with threshold lambda: sign(x) * max(\|x\| - lambda, 0). |
| `func SymmetricExtend(x []float64, pad int) []float64` | SymmetricExtend returns a copy of x extended by pad samples on each side using whole-sample symmetric (mirror) boundary handling, so the result has… |
| `func Threshold(x, lambda float64, mode ThresholdMode) float64` | Threshold applies the thresholding rule selected by mode to x with threshold lambda. |
| `func ThresholdSlice(xs []float64, lambda float64, mode ThresholdMode) []float64` | ThresholdSlice returns a new slice with the thresholding rule selected by mode applied elementwise to xs with threshold lambda. |
| `func UniversalThreshold(n int, sigma float64) float64` | UniversalThreshold returns the Donoho-Johnstone universal (VisuShrink) threshold sigma * sqrt(2 * ln(n)) for a signal of length n and noise standard… |
| `func Upsample(x []float64, factor int) []float64` | Upsample returns x with factor-1 zeros inserted after each sample, producing a slice of length len(x)*factor. |
| `func WaveRec(d *Decomposition) []float64` | WaveRec reconstructs the original signal from a `Decomposition`, inverting `WaveDec` to within floating-point rounding. |
| `func ZeroExtend(x []float64, pad int) []float64` | ZeroExtend returns a copy of x extended by pad zero samples on each side, so the result has length len(x)+2*pad. |
| `func ZeroPadToPowerOfTwo(x []float64) []float64` | ZeroPadToPowerOfTwo returns a copy of x zero-padded on the right up to the next power-of-two length. |

### Types

| Type | What it is |
| --- | --- |
| `Coefficients` | Coefficients holds the approximation and detail coefficients produced by a single level of the discrete wavelet transform. |
| `DWT2DResult` | DWT2DResult holds the four subbands produced by one level of a separable two-dimensional discrete wavelet transform. |
| `Decomposition` | Decomposition is the result of a multi-level (multiresolution) discrete wavelet transform. |
| `PacketTree` | PacketTree is a full (balanced) wavelet packet decomposition. |
| `ThresholdMode` | ThresholdMode selects the shrinkage rule applied to wavelet coefficients. |
| `Wavelet` | Wavelet is an orthonormal discrete wavelet, described by its analysis (decomposition) and synthesis (reconstruction) filter banks. |

<details>
<summary><code>Coefficients</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Transform(signal []float64, w Wavelet) Coefficients` | Transform computes one level of the forward transform of signal with wavelet w and returns the result as a `Coefficients` value. |
| `func (c Coefficients) Energy() float64` | Energy returns the total energy of the coefficients, the sum of squares of the approximation and detail coefficients. |
| `func (c Coefficients) Inverse(w Wavelet) []float64` | Inverse reconstructs the signal from the coefficients using wavelet w, inverting `Transform`. |
| `func (c Coefficients) Length() int` | Length returns the reconstructed signal length, 2*len(c.Approx). |

</details>

<details>
<summary><code>DWT2DResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DWT2D(image [][]float64, w Wavelet) *DWT2DResult` | DWT2D computes one level of the separable two-dimensional discrete wavelet transform of image with wavelet w, using periodic boundary handling on… |
| `func (r *DWT2DResult) Dims() (rows, cols int)` | Dims returns the height and width of each subband in the result. |
| `func (r *DWT2DResult) Energy() float64` | Energy returns the total energy (sum of squares) across all four subbands. |
| `func (r *DWT2DResult) Subbands() (ll, lh, hl, hh [][]float64)` | Subbands returns the four subbands in the fixed order LL, LH, HL, HH. |

</details>

<details>
<summary><code>Decomposition</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func WaveDec(signal []float64, w Wavelet, levels int) (*Decomposition, error)` | WaveDec computes a multi-level forward wavelet transform of signal with wavelet w, cascading the low-pass branch for the requested number of levels. |
| `func (d *Decomposition) ApproxComponent() []float64` | ApproxComponent projects the coarsest approximation back to the original signal resolution, yielding the smooth multiresolution component of the… |
| `func (d *Decomposition) DetailAt(level int) []float64` | DetailAt returns the detail coefficients captured at the given level, where level 1 is the finest scale and d.Levels() is the coarsest. |
| `func (d *Decomposition) DetailComponents() [][]float64` | DetailComponents projects the detail coefficients of each level back to the original signal resolution. |
| `func (d *Decomposition) Energy() float64` | Energy returns the total energy (sum of squares) of every stored coefficient. |
| `func (d *Decomposition) Flatten() []float64` | Flatten returns all coefficients concatenated into a single slice in the order [approx, detail_coarsest, ..., detail_finest], matching the common… |
| `func (d *Decomposition) Levels() int` | Levels returns the number of decomposition levels stored in d. |
| `func (d *Decomposition) MRAComponents() (approx []float64, details [][]float64)` | MRAComponents returns the additive multiresolution components of the signal: the full-length smooth approximation and one full-length detail… |
| `func (d *Decomposition) Reconstruct() []float64` | Reconstruct rebuilds the original signal from the stored coefficients by cascading the inverse transform from the coarsest to the finest level. |

</details>

<details>
<summary><code>PacketTree</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func WaveletPacket(signal []float64, w Wavelet, levels int) (*PacketTree, error)` | WaveletPacket computes a full wavelet packet decomposition of signal with wavelet w down to the requested number of levels. |
| `func (t *PacketTree) Leaves() [][]float64` | Leaves returns copies of every node at the deepest level, left to right. |
| `func (t *PacketTree) Levels() int` | Levels returns the depth of the packet tree. |
| `func (t *PacketTree) Node(level, index int) []float64` | Node returns a copy of the coefficients stored at the given level and index. |
| `func (t *PacketTree) NodeCount(level int) int` | NodeCount returns the number of nodes at the given level, which is 2^level. |
| `func (t *PacketTree) Reconstruct() []float64` | Reconstruct rebuilds the original signal from the leaf nodes of the tree by inverting the packet transform level by level, achieving perfect… |

</details>

<details>
<summary><code>Wavelet</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DB2() Wavelet` | DB2 returns the Daubechies wavelet with two vanishing moments (four taps). |
| `func DB4() Wavelet` | DB4 returns the Daubechies wavelet with four vanishing moments (eight taps). |
| `func Daubechies(p int) (Wavelet, error)` | Daubechies returns the Daubechies wavelet with p vanishing moments (a filter of 2p taps). |
| `func Haar() Wavelet` | Haar returns the Haar wavelet, the orthonormal Daubechies wavelet with one vanishing moment (equivalently db1). |
| `func (w Wavelet) DecHi() []float64` | DecHi returns a copy of the analysis (decomposition) high-pass filter. |
| `func (w Wavelet) DecLo() []float64` | DecLo returns a copy of the analysis (decomposition) low-pass filter. |
| `func (w Wavelet) IsOrthogonal(tol float64) bool` | IsOrthogonal reports whether the scaling filter satisfies the orthonormality conditions of an orthogonal wavelet to within tol: unit norm and… |
| `func (w Wavelet) Name() string` | Name returns the short name of the wavelet family, such as "haar" or "db4". |
| `func (w Wavelet) RecHi() []float64` | RecHi returns a copy of the synthesis (reconstruction) high-pass filter. |
| `func (w Wavelet) RecLo() []float64` | RecLo returns a copy of the synthesis (reconstruction) low-pass filter. |
| `func (w Wavelet) Taps() int` | Taps returns the number of filter coefficients (the filter length). |
| `func (w Wavelet) VanishingMoments() int` | VanishingMoments returns the number of vanishing moments of the wavelet. |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/wavelet).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
