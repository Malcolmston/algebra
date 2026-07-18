// Package wavelet provides discrete wavelet transforms and related
// multiresolution tools implemented with the Go standard library only.
//
// The package is organised into a few areas:
//
//   - Wavelet families (see wavelet.go): the orthonormal [Haar] wavelet and
//     the Daubechies family via [Daubechies], with the common shortcuts [DB2]
//     and [DB4]. A [Wavelet] exposes its analysis and synthesis filter banks
//     through [Wavelet.DecLo], [Wavelet.DecHi], [Wavelet.RecLo] and
//     [Wavelet.RecHi].
//   - Single-level transforms (see dwt.go): the forward [DWT] and its inverse
//     [IDWT], together with the [Coefficients] container.
//   - Multiresolution analysis (see mra.go): multi-level decomposition
//     [WaveDec] and reconstruction [WaveRec], the [Decomposition] type, and
//     projection of each scale back to the original resolution via
//     [Decomposition.MRAComponents].
//   - Wavelet packets (see packet.go): the full binary [PacketTree] produced by
//     [WaveletPacket], with perfect reconstruction and the Coifman-Wickerhauser
//     cost measures [ShannonEntropy] and [LogEnergy].
//   - Thresholding and denoising (see threshold.go): [SoftThreshold],
//     [HardThreshold], the VisuShrink [UniversalThreshold], robust noise
//     estimation via [EstimateNoiseSigma], and the end-to-end [Denoise].
//   - Two-dimensional transforms (see dwt2d.go): a single level [DWT2D] and its
//     inverse [IDWT2D] returning the [DWT2DResult] subbands.
//   - Signal utilities (see signal.go): boundary extensions, linear
//     [Convolve], [Downsample]/[Upsample], and small numeric helpers.
//
// All transforms use periodic (circular) boundary handling, so every subband
// has exactly half the length of its parent and the analysis operator is
// orthogonal. Reconstruction is therefore the exact transpose of analysis and
// achieves perfect reconstruction to within floating-point rounding. Because
// the filters are orthonormal the transform preserves energy (Parseval's
// relation): the sum of squares of the approximation and detail coefficients
// equals the sum of squares of the input.
//
// The package is deterministic and depends only on the math, sort and errors
// standard-library packages.
package wavelet
