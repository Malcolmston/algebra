# randommatrix — implements the core objects of random matrix theory (RMT) using only the Go standard

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/randommatrix.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/randommatrix)

Package randommatrix implements the core objects of random matrix theory (RMT)
using only the Go standard library.

Random matrix theory studies the statistical behaviour of the eigenvalues of
matrices whose entries are drawn at random. Despite the randomness of the
individual entries, the collective spectrum obeys sharp, universal laws. This
package provides seeded generators for the classical Gaussian and Wishart
ensembles together with the analytic machinery used to describe their spectra.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/randommatrix"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func BinomialBig(n, k int) *big.Int` | BinomialBig returns the exact binomial coefficient C(n, k) as a big integer. |
| `func BinomialCoefficient(n, k int) float64` | BinomialCoefficient returns C(n, k) as a float64, computed exactly with big integers and then converted. |
| `func CatalanBig(n int) *big.Int` | CatalanBig returns the n-th Catalan number C_n = C(2n, n)/(n+1) as an exact big integer for n >= 0. |
| `func CatalanNumber(n int) float64` | CatalanNumber returns the n-th Catalan number as a float64. |
| `func CentralSpectralMoment(eigs []float64, k int) float64` | CentralSpectralMoment returns the k-th central moment (1/n) sum (lambda_i - mean)^k of the eigenvalue sample. |
| `func ConditionNumberSymmetric(m *Matrix) (float64, error)` | ConditionNumberSymmetric returns the spectral condition number \|lambda_max\| / \|lambda_min\| of a real symmetric matrix. |
| `func ConsecutiveSpacingRatios(eigs []float64) []float64` | ConsecutiveSpacingRatios returns the restricted ratios r_i = min(s_i, s_{i-1}) / max(s_i, s_{i-1}) in [0, 1], which need no spectral unfolding. |
| `func DeterminantSymmetric(m *Matrix) (float64, error)` | DeterminantSymmetric returns the determinant of a real symmetric matrix as the product of its eigenvalues. |
| `func EigenvaluesHermitian(m *CMatrix) ([]float64, error)` | EigenvaluesHermitian returns the eigenvalues of a Hermitian matrix in ascending order. |
| `func EigenvaluesSymmetric(m *Matrix) ([]float64, error)` | EigenvaluesSymmetric returns the eigenvalues of a real symmetric matrix in ascending order. |
| `func EmpiricalCDF(data []float64, x float64) float64` | EmpiricalCDF returns the value of the empirical cumulative distribution function of the sample at x, i.e. |
| `func EmpiricalCauchy(eigs []float64, z complex128) complex128` | EmpiricalCauchy returns the Cauchy transform (1/n) sum 1/(z - lambda_i), the negative of the Stieltjes transform. |
| `func EmpiricalStieltjes(eigs []float64, z complex128) complex128` | EmpiricalStieltjes returns the Stieltjes transform (1/n) sum 1/(lambda_i - z) of the empirical spectral distribution at complex argument z. |
| `func FactorialBig(n int) *big.Int` | FactorialBig returns n! as an exact big integer for n >= 0 and nil for n < 0. |
| `func FreeConvolutionMoments(momentsA, momentsB []float64) []float64` | FreeConvolutionMoments returns the moments of the free additive convolution of two distributions given their moment sequences. |
| `func FreeConvolutionSemicircleRadius(r1, r2 float64) float64` | FreeConvolutionSemicircleRadius returns the radius of the semicircle law obtained as the free additive convolution of two semicircle laws of the… |
| `func FreeConvolutionSemicircleVariance(var1, var2 float64) float64` | FreeConvolutionSemicircleVariance returns the variance of the free additive convolution of two semicircle laws, which is the sum of their variances. |
| `func FreeCumulantsFromMoments(moments []float64) []float64` | FreeCumulantsFromMoments returns the free cumulants kappa_1..kappa_N implied by the moments m_1..m_N by inverting the free moment-cumulant relation. |
| `func GOEMeanRatio() float64` | GOEMeanRatio returns the surmise value 0.5307 of the mean restricted spacing ratio for the GOE. |
| `func GOESpacingCDF(s float64) float64` | GOESpacingCDF returns the GOE spacing distribution 1 - exp(-pi s^2/4). |
| `func GOESpacingDensity(s float64) float64` | GOESpacingDensity returns the GOE (beta = 1) Wigner surmise density (pi/2) s exp(-pi s^2/4). |
| `func GSEMeanRatio() float64` | GSEMeanRatio returns the surmise value 0.6744 of the mean restricted spacing ratio for the GSE. |
| `func GSESpacingCDF(s float64) float64` | GSESpacingCDF returns the GSE spacing distribution, integrated numerically. |
| `func GSESpacingDensity(s float64) float64` | GSESpacingDensity returns the GSE (beta = 4) Wigner surmise density (2^18/(3^6 pi^3)) s^4 exp(-64 s^2/(9 pi)). |
| `func GUEMeanRatio() float64` | GUEMeanRatio returns the surmise value 0.5996 of the mean restricted spacing ratio for the GUE. |
| `func GUESpacingCDF(s float64) float64` | GUESpacingCDF returns the GUE spacing distribution, integrated numerically. |
| `func GUESpacingDensity(s float64) float64` | GUESpacingDensity returns the GUE (beta = 2) Wigner surmise density (32/pi^2) s^2 exp(-4 s^2/pi). |
| `func GaussianVector(n int, mean, std float64, seed int64) []float64` | GaussianVector returns n independent N(mean, std^2) samples. |
| `func KolmogorovSmirnovStatistic(data []float64, cdf func(float64) float64) float64` | KolmogorovSmirnovStatistic returns the maximum absolute difference between the empirical CDF of the sample and the reference CDF, evaluated at the… |
| `func LargestEigenvalue(m *Matrix) (float64, error)` | LargestEigenvalue returns the algebraically largest eigenvalue of a real symmetric matrix. |
| `func MarchenkoPasturAtomMass(c float64) float64` | MarchenkoPasturAtomMass returns the mass max(0, 1 - 1/c) of the atom at zero of the Marchenko-Pastur law with aspect ratio c. |
| `func MarchenkoPasturCDF(x, c, sigma2 float64) float64` | MarchenkoPasturCDF returns the cumulative distribution function of the Marchenko-Pastur law with aspect ratio c and variance sigma2 at x, including… |
| `func MarchenkoPasturDensity(x, c, sigma2 float64) float64` | MarchenkoPasturDensity returns the continuous part of the Marchenko-Pastur density with aspect ratio c and variance sigma2 at x. |
| `func MarchenkoPasturHasAtom(c float64) bool` | MarchenkoPasturHasAtom reports whether the Marchenko-Pastur law with aspect ratio c places an atom at zero (which happens exactly when c > 1). |
| `func MarchenkoPasturMean(c, sigma2 float64) float64` | MarchenkoPasturMean returns the mean sigma2 of the Marchenko-Pastur law. |
| `func MarchenkoPasturMoment(k int, c, sigma2 float64) float64` | MarchenkoPasturMoment returns the k-th moment of the Marchenko-Pastur law with aspect ratio c and variance sigma2. |
| `func MarchenkoPasturStieltjes(z complex128, c float64) complex128` | MarchenkoPasturStieltjes returns the Stieltjes transform of the standard (sigma2 = 1) Marchenko-Pastur law with ratio c at complex argument z, using… |
| `func MarchenkoPasturSupport(c, sigma2 float64) (lo, hi float64)` | MarchenkoPasturSupport returns the edges lambda_minus and lambda_plus of the Marchenko-Pastur law with aspect ratio c and per-entry variance sigma2. |
| `func MarchenkoPasturVariance(c, sigma2 float64) float64` | MarchenkoPasturVariance returns the variance c*sigma2^2 of the Marchenko-Pastur law. |
| `func Max(xs []float64) float64` | Max returns the maximum of a non-empty slice; it returns NaN for the empty slice. |
| `func Mean(xs []float64) float64` | Mean returns the arithmetic mean of xs, or NaN for the empty slice. |
| `func MeanConsecutiveRatio(eigs []float64) float64` | MeanConsecutiveRatio returns the mean of the restricted spacing ratios of the eigenvalue sample. |
| `func MeanSpacing(eigs []float64) float64` | MeanSpacing returns the mean nearest-neighbour spacing of the eigenvalue sample. |
| `func Min(xs []float64) float64` | Min returns the minimum of a non-empty slice; it returns NaN for the empty slice. |
| `func MomentsFromFreeCumulants(kappa []float64) []float64` | MomentsFromFreeCumulants returns the moments m_1..m_N implied by the free cumulants kappa_1..kappa_N via the free moment-cumulant relation m_n = sum… |
| `func NearestNeighbourSpacings(eigs []float64) []float64` | NearestNeighbourSpacings returns the consecutive gaps of the sorted eigenvalue sample. |
| `func NonCrossingPartitionsCount(n int) float64` | NonCrossingPartitionsCount returns the number of non-crossing partitions of an n-element set, which equals the n-th Catalan number. |
| `func NormalizedSpacings(eigs []float64) []float64` | NormalizedSpacings returns the nearest-neighbour spacings rescaled by their mean, a crude global unfolding that produces unit mean spacing. |
| `func PoissonMeanRatio() float64` | PoissonMeanRatio returns the mean restricted spacing ratio 2 ln 2 - 1 for a Poisson spectrum. |
| `func PoissonRatioDensity(r float64) float64` | PoissonRatioDensity returns the surmise density 1/(1+r)^2 of the unrestricted consecutive-spacing ratio for a Poisson spectrum. |
| `func PoissonSpacingCDF(s float64) float64` | PoissonSpacingCDF returns the cumulative distribution 1 - exp(-s) of Poisson level spacings. |
| `func PoissonSpacingDensity(s float64) float64` | PoissonSpacingDensity returns the density exp(-s) of nearest-neighbour level spacings for an uncorrelated (Poisson) spectrum with unit mean spacing. |
| `func RTransformSeries(kappa []float64) []float64` | RTransformSeries returns the coefficients of the R-transform power series R(z) = sum_{n>=1} kappa_n z^(n-1) from the free cumulants kappa_1..kappa_N.… |
| `func RatioSurmiseDensity(r, beta float64) float64` | RatioSurmiseDensity returns the Atas et al. |
| `func RescaleLargestEigenvalue(lambdaMax float64, n int) float64` | RescaleLargestEigenvalue maps the largest eigenvalue lambdaMax of an n-by-n Gaussian ensemble to the Tracy-Widom scale via (lambdaMax - 2 sqrt(n)) *… |
| `func SampleGOE(n int, seed int64, scale float64) ([]float64, error)` | SampleGOE returns the eigenvalues of GOE(n)/scale in ascending order. |
| `func SampleGUE(n int, seed int64, scale float64) ([]float64, error)` | SampleGUE returns the eigenvalues of GUE(n)/scale in ascending order. |
| `func SampleWishartReal(n, p int, seed int64) ([]float64, error)` | SampleWishartReal returns the eigenvalues of the real sample covariance matrix (1/p) G Gᵀ in ascending order. |
| `func SemicircleCDF(x, radius float64) float64` | SemicircleCDF returns the cumulative distribution function of the semicircle law of radius R evaluated at x. |
| `func SemicircleDensity(x, radius float64) float64` | SemicircleDensity returns the Wigner semicircle density (2/(pi R^2)) sqrt(R^2 - x^2) at x for radius R; it is zero outside [-R, R]. |
| `func SemicircleFreeCumulants(variance float64, N int) []float64` | SemicircleFreeCumulants returns the free cumulants of the semicircle law of the given variance up to order N. Only the second cumulant is non-zero… |
| `func SemicircleMoment(k int, radius float64) float64` | SemicircleMoment returns the k-th moment of the semicircle law of radius R. Odd moments vanish; the (2m)-th moment equals the m-th Catalan number… |
| `func SemicirclePeak(radius float64) float64` | SemicirclePeak returns the maximum value 2/(pi R) of the semicircle density, attained at x = 0. |
| `func SemicircleRadius(variance float64) float64` | SemicircleRadius returns the radius R = 2*sqrt(variance) of the Wigner semicircle law whose second moment (variance) equals the given value. |
| `func SemicircleStieltjes(z complex128, radius float64) complex128` | SemicircleStieltjes returns the Stieltjes transform m(z) of the semicircle law of radius R, the value of the integral of 1/(x-z) against the density. |
| `func SemicircleSupport(radius float64) (lo, hi float64)` | SemicircleSupport returns the endpoints (-R, R) of the semicircle law of radius R. |
| `func SemicircleVariance(radius float64) float64` | SemicircleVariance returns the variance R^2/4 of the semicircle law of radius R. |
| `func SmallestEigenvalue(m *Matrix) (float64, error)` | SmallestEigenvalue returns the algebraically smallest eigenvalue of a real symmetric matrix. |
| `func SortSpacings(spacings []float64) []float64` | SortSpacings returns the ascending-sorted spacings, convenient for building empirical spacing distributions. |
| `func SortedCopy(xs []float64) []float64` | SortedCopy returns a sorted (ascending) copy of xs, leaving the input intact. |
| `func SpacingRatios(eigs []float64) []float64` | SpacingRatios returns the ratios s_{i+1}/s_i of consecutive nearest-neighbour spacings of the eigenvalue sample. |
| `func SpectralGaps(eigs []float64) []float64` | SpectralGaps returns the consecutive differences of the sorted eigenvalues. |
| `func SpectralMean(eigs []float64) float64` | SpectralMean returns the mean of the eigenvalue sample. |
| `func SpectralMoment(eigs []float64, k int) float64` | SpectralMoment returns the k-th empirical spectral moment (1/n) sum lambda_i^k of the eigenvalue sample eigs. |
| `func SpectralMoments(eigs []float64, kmax int) []float64` | SpectralMoments returns the empirical spectral moments of orders 1..kmax. |
| `func SpectralQuantile(eigs []float64, q float64) float64` | SpectralQuantile returns the q-quantile (0 <= q <= 1) of the eigenvalue sample using linear interpolation between order statistics. |
| `func SpectralRadius(m *Matrix) (float64, error)` | SpectralRadius returns the largest absolute eigenvalue of a real symmetric matrix. |
| `func SpectralVariance(eigs []float64) float64` | SpectralVariance returns the population variance of the eigenvalue sample. |
| `func TracyWidomApproxParams(beta int) (k, theta, shift float64, ok bool)` | TracyWidomApproxParams returns the shifted-gamma approximation parameters (shape k, scale theta, location shift) for the Tracy-Widom distribution of… |
| `func TracyWidomCDF(s float64, beta int) float64` | TracyWidomCDF returns the Chiani shifted-gamma approximation of the Tracy-Widom cumulative distribution function of Dyson index beta at s. |
| `func TracyWidomDensity(s float64, beta int) float64` | TracyWidomDensity returns the Chiani shifted-gamma approximation of the Tracy-Widom probability density of Dyson index beta at s. |
| `func TracyWidomEdgeCenter(n int) float64` | TracyWidomEdgeCenter returns the soft-edge centring 2*sqrt(n) for the largest eigenvalue of an n-by-n Gaussian ensemble in the normalisation where… |
| `func TracyWidomEdgeScale(n int) float64` | TracyWidomEdgeScale returns the soft-edge scaling n^(-1/6), so that (lambda_max - 2 sqrt(n)) / scale converges to the Tracy-Widom law. |
| `func TracyWidomExcessKurtosis(beta int) float64` | TracyWidomExcessKurtosis returns the excess kurtosis 6/k of the shifted-gamma approximation to the Tracy-Widom distribution of Dyson index beta. |
| `func TracyWidomMean(beta int) float64` | TracyWidomMean returns the tabulated mean of the Tracy-Widom distribution of Dyson index beta. |
| `func TracyWidomSkewness(beta int) float64` | TracyWidomSkewness returns the skewness 2/sqrt(k) of the shifted-gamma approximation to the Tracy-Widom distribution of Dyson index beta. |
| `func TracyWidomStdDev(beta int) float64` | TracyWidomStdDev returns the standard deviation of the Tracy-Widom distribution of Dyson index beta. |
| `func TracyWidomSupportedBeta(beta int) bool` | TracyWidomSupportedBeta reports whether beta is one of the supported Dyson indices 1, 2 or 4. |
| `func TracyWidomVariance(beta int) float64` | TracyWidomVariance returns the variance k*theta^2 of the shifted-gamma approximation, which matches the tabulated Tracy-Widom variance. |
| `func Unfold(eigs []float64, cdf func(float64) float64) []float64` | Unfold maps each eigenvalue lambda to n*cdf(lambda), where n is the number of eigenvalues and cdf is the integrated limiting spectral density. |
| `func UnfoldedSpacings(eigs []float64, cdf func(float64) float64) []float64` | UnfoldedSpacings returns the nearest-neighbour spacings of the eigenvalues after unfolding with the given integrated spectral density. |
| `func WignerSurmise(s, beta float64) float64` | WignerSurmise returns the Wigner surmise density for Dyson index beta at spacing s, normalised to unit total mass and unit mean spacing. |
| `func WignerSurmiseCDF(s, beta float64) float64` | WignerSurmiseCDF returns the cumulative distribution of the Wigner surmise for Dyson index beta at spacing s, integrated numerically. |
| `func WishartTracyWidomCenter(n, p int) float64` | WishartTracyWidomCenter returns the Johnstone centring (sqrt(n-1) + sqrt(p))^2 for the largest eigenvalue of a real Wishart matrix built from an… |
| `func WishartTracyWidomScale(n, p int) float64` | WishartTracyWidomScale returns the Johnstone scaling (sqrt(n-1) + sqrt(p)) (1/sqrt(n-1) + 1/sqrt(p))^(1/3) for the largest eigenvalue of a real… |

### Types

| Type | What it is |
| --- | --- |
| `CMatrix` | CMatrix is a dense, row-major matrix of complex128 values. |
| `EigenSym` | EigenSym holds the eigen-decomposition of a real symmetric matrix. |
| `Histogram` | Histogram is a fixed-width histogram of a real sample. |
| `Matrix` | Matrix is a dense, row-major matrix of float64 values. |

<details>
<summary><code>CMatrix</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func CIdentity(n int) *CMatrix` | CIdentity returns the n-by-n complex identity matrix. |
| `func GSE(n int, seed int64) *CMatrix` | GSE returns a sample from the Gaussian Symplectic Ensemble as a 2n-by-2n Hermitian, self-dual complex matrix built from an n-by-n array of… |
| `func GUE(n int, seed int64) *CMatrix` | GUE returns a sample from the Gaussian Unitary Ensemble: a Hermitian n-by-n matrix whose diagonal entries are real N(0,1) and whose off-diagonal… |
| `func GaussianCMatrix(rows, cols int, std float64, seed int64) *CMatrix` | GaussianCMatrix returns a rows-by-cols complex matrix whose real and imaginary parts are independent N(0, std^2) samples. |
| `func GinibreComplex(n int, seed int64) *CMatrix` | GinibreComplex returns an n-by-n complex matrix whose entries have real and imaginary parts drawn independently from N(0,1/2), so that each entry has… |
| `func LaguerreUnitary(n, p int, seed int64) *CMatrix` | LaguerreUnitary is a synonym for WishartComplex. |
| `func NewCMatrix(rows, cols int) *CMatrix` | NewCMatrix returns a rows-by-cols complex matrix of zeros. |
| `func NewCMatrixFrom(rows, cols int, data []complex128) *CMatrix` | NewCMatrixFrom builds a complex matrix from a flat row-major slice, which is copied. |
| `func WishartComplex(n, p int, seed int64) *CMatrix` | WishartComplex returns the complex Wishart matrix G G†, where G is an n-by-p complex matrix with real and imaginary parts drawn independently from… |
| `func (m *CMatrix) At(i, j int) complex128` | At returns the entry at row i, column j. |
| `func (m *CMatrix) Clone() *CMatrix` | Clone returns an independent deep copy. |
| `func (m *CMatrix) Cols() int` | Cols returns the number of columns. |
| `func (m *CMatrix) ConjugateTranspose() *CMatrix` | ConjugateTranspose returns the conjugate (Hermitian) transpose m†. |
| `func (m *CMatrix) FrobeniusNorm() float64` | FrobeniusNorm returns the Frobenius norm of the complex matrix. |
| `func (m *CMatrix) Hermitianize() *CMatrix` | Hermitianize returns (m + m†)/2, the Hermitian part of a square matrix. |
| `func (m *CMatrix) ImagPart() *Matrix` | ImagPart returns the matrix of imaginary parts. |
| `func (m *CMatrix) IsHermitian(tol float64) bool` | IsHermitian reports whether m is square and Hermitian to within tol. |
| `func (m *CMatrix) IsSquare() bool` | IsSquare reports whether the matrix is square. |
| `func (m *CMatrix) Minus(other *CMatrix) (*CMatrix, error)` | Minus returns the entrywise difference m - other. |
| `func (m *CMatrix) Mul(other *CMatrix) (*CMatrix, error)` | Mul returns the matrix product m * other. |
| `func (m *CMatrix) Plus(other *CMatrix) (*CMatrix, error)` | Plus returns the entrywise sum m + other. |
| `func (m *CMatrix) RealPart() *Matrix` | RealPart returns the matrix of real parts. |
| `func (m *CMatrix) Rows() int` | Rows returns the number of rows. |
| `func (m *CMatrix) Scale(s complex128) *CMatrix` | Scale multiplies every entry by the complex scalar s. |
| `func (m *CMatrix) Set(i, j int, v complex128)` | Set assigns v to the entry at row i, column j. |
| `func (m *CMatrix) Trace() complex128` | Trace returns the sum of the diagonal entries. |

</details>

<details>
<summary><code>EigenSym</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func EigSymmetric(m *Matrix) (*EigenSym, error)` | EigSymmetric computes the full eigen-decomposition of a real symmetric matrix using the cyclic Jacobi algorithm. |

</details>

<details>
<summary><code>Histogram</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func EmpiricalSpectralDensity(eigs []float64, bins int) *Histogram` | EmpiricalSpectralDensity returns a histogram of the eigenvalue sample using the given number of bins over the observed spectral range. |
| `func NewHistogram(data []float64, bins int, lo, hi float64) *Histogram` | NewHistogram bins the data into the given number of equal-width bins spanning [lo, hi]. |
| `func NewHistogramAuto(data []float64, bins int) *Histogram` | NewHistogramAuto bins the data into the given number of bins spanning its observed range [min, max]. |
| `func (h *Histogram) BinCenters() []float64` | BinCenters returns the midpoint of each bin. |
| `func (h *Histogram) BinWidth() float64` | BinWidth returns the common width of the bins. |
| `func (h *Histogram) Density() []float64` | Density returns the empirical probability density: Counts`i` divided by (Total * BinWidth). |
| `func (h *Histogram) ModeBin() int` | ModeBin returns the index of the bin with the most samples. |
| `func (h *Histogram) NumBins() int` | NumBins returns the number of bins. |
| `func (h *Histogram) Probabilities() []float64` | Probabilities returns the fraction of samples in each bin. |

</details>

<details>
<summary><code>Matrix</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Diag(d []float64) *Matrix` | Diag returns a square diagonal matrix whose diagonal is d. |
| `func GOE(n int, seed int64) *Matrix` | GOE returns a sample from the Gaussian Orthogonal Ensemble: a real symmetric n-by-n matrix whose off-diagonal entries are independent N(0,1) and… |
| `func GaussianMatrix(rows, cols int, mean, std float64, seed int64) *Matrix` | GaussianMatrix returns a rows-by-cols matrix of independent N(mean, std^2) entries. |
| `func GinibreReal(n int, seed int64) *Matrix` | GinibreReal returns an n-by-n matrix of independent standard normal real entries (the real Ginibre ensemble). |
| `func Identity(n int) *Matrix` | Identity returns the n-by-n identity matrix. |
| `func LaguerreOrthogonal(n, p int, seed int64) *Matrix` | LaguerreOrthogonal is a synonym for WishartReal, emphasising its role as the beta = 1 Laguerre ensemble. |
| `func NewMatrix(rows, cols int) *Matrix` | NewMatrix returns a rows-by-cols matrix whose entries are all zero. |
| `func NewMatrixFrom(rows, cols int, data []float64) *Matrix` | NewMatrixFrom builds a rows-by-cols matrix from a flat, row-major slice. |
| `func NewMatrixFromRows(rows [][]float64) *Matrix` | NewMatrixFromRows builds a matrix from a slice of equal-length rows. |
| `func SampleCovarianceReal(n, p int, seed int64) *Matrix` | SampleCovarianceReal returns (1/p) G Gᵀ where G is n-by-p standard normal. |
| `func WignerMatrix(n int, offStd, diagStd float64, seed int64) *Matrix` | WignerMatrix returns a real symmetric n-by-n Wigner matrix with independent off-diagonal entries of standard deviation offStd and diagonal entries of… |
| `func WishartReal(n, p int, seed int64) *Matrix` | WishartReal returns the real Wishart matrix G Gᵀ, where G is an n-by-p matrix of independent standard normal entries. |
| `func Zeros(rows, cols int) *Matrix` | Zeros returns a rows-by-cols zero matrix; it is an alias for NewMatrix. |
| `func (m *Matrix) Add(i, j int, v float64)` | Add increments the entry at row i, column j by v. |
| `func (m *Matrix) At(i, j int) float64` | At returns the entry at row i, column j. |
| `func (m *Matrix) Clone() *Matrix` | Clone returns an independent deep copy of the matrix. |
| `func (m *Matrix) Col(j int) []float64` | Col returns a copy of column j. |
| `func (m *Matrix) Cols() int` | Cols returns the number of columns. |
| `func (m *Matrix) Data() []float64` | Data returns a copy of the underlying row-major slice. |
| `func (m *Matrix) Diagonal() []float64` | Diagonal returns a copy of the main diagonal. |
| `func (m *Matrix) Equals(other *Matrix, tol float64) bool` | Equals reports whether every entry of m and other agree to within tol. |
| `func (m *Matrix) FrobeniusNorm() float64` | FrobeniusNorm returns the Frobenius norm, sqrt of the sum of squared entries. |
| `func (m *Matrix) IsSquare() bool` | IsSquare reports whether the matrix is square. |
| `func (m *Matrix) IsSymmetric(tol float64) bool` | IsSymmetric reports whether m is square and symmetric to within tol. |
| `func (m *Matrix) MaxAbs() float64` | MaxAbs returns the largest absolute value of any entry. |
| `func (m *Matrix) Minus(other *Matrix) (*Matrix, error)` | Minus returns the entrywise difference m - other. |
| `func (m *Matrix) Mul(other *Matrix) (*Matrix, error)` | Mul returns the matrix product m * other. |
| `func (m *Matrix) MulVec(v []float64) ([]float64, error)` | MulVec returns the matrix-vector product m * v. |
| `func (m *Matrix) Plus(other *Matrix) (*Matrix, error)` | Plus returns the entrywise sum m + other. |
| `func (m *Matrix) Row(i int) []float64` | Row returns a copy of row i. |
| `func (m *Matrix) Rows() int` | Rows returns the number of rows. |
| `func (m *Matrix) Scale(s float64) *Matrix` | Scale returns a new matrix equal to m with every entry multiplied by s. |
| `func (m *Matrix) Set(i, j int, v float64)` | Set assigns v to the entry at row i, column j. |
| `func (m *Matrix) String() string` | String renders the matrix with entries formatted to three decimals. |
| `func (m *Matrix) SubMatrix(r0, r1, c0, c1 int) *Matrix` | SubMatrix returns the r0:r1, c0:c1 block (half-open ranges). |
| `func (m *Matrix) Symmetrize() *Matrix` | Symmetrize returns (m + mᵀ)/2, the symmetric part of a square matrix. |
| `func (m *Matrix) Trace() float64` | Trace returns the sum of the diagonal entries. |
| `func (m *Matrix) Transpose() *Matrix` | Transpose returns the transpose of the matrix. |

</details>

### Constants

`BetaOrthogonal`, `BetaUnitary`, `BetaSymplectic`

### Variables

`ErrDimensionMismatch`, `ErrNotSquare`, `ErrNotSymmetric`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/randommatrix).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
