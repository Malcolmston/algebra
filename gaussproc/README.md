# gaussproc — provides a self-contained toolkit for Gaussian process (GP) modelling

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/gaussproc.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/gaussproc)

Package gaussproc provides a self-contained toolkit for Gaussian process (GP)
modelling, kernel methods and the small amount of numerical linear algebra
they require, implemented using only the Go standard library.

A `Kernel` maps a pair of real vectors to a covariance value. The package
ships the standard library of stationary and non-stationary kernels: the
squared-exponential / radial-basis-function kernel (`RBF`); the Matérn family
for smoothness parameters 1/2, 3/2 and 5/2 (`Matern12`, `Matern32`,
`Matern52`); the `RationalQuadratic` kernel; the exponentially-periodic
(`Periodic`) kernel; the `GammaExponential` and `Cosine` kernels; the
non-stationary `Linear` and `Polynomial` dot-product kernels; and the
`Constant` and `WhiteNoise` kernels. Kernels compose through an algebra of
sums (`Sum`, `SumKernel`), products (`Product`, `ProductKernel`) and scalar
scaling (`Scale`, `ScaledKernel`); arbitrary functions become kernels via
`KernelFunc`.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/gaussproc"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func BackSubstitution(u Matrix, b []float64) ([]float64, error)` | BackSubstitution solves the upper-triangular system U·x = b and returns x. |
| `func CholeskySolve(l Matrix, b []float64) ([]float64, error)` | CholeskySolve solves the symmetric positive-definite system A·x = b given the lower-triangular Cholesky factor L of A (as returned by `Cholesky`). |
| `func DeterminantFromCholesky(l Matrix) float64` | DeterminantFromCholesky returns det A computed from the Cholesky factor L of A, as the square of the product of the diagonal of L. |
| `func Distance(a, b []float64) float64` | Distance returns the Euclidean distance between vectors a and b. |
| `func Dot(a, b []float64) float64` | Dot returns the Euclidean inner product of vectors a and b. |
| `func ForwardSubstitution(l Matrix, b []float64) ([]float64, error)` | ForwardSubstitution solves the lower-triangular system L·x = b and returns x. |
| `func GaussianLogPDF(x, mean, variance float64) float64` | GaussianLogPDF returns the log density of a univariate normal distribution with the given mean and variance evaluated at x. |
| `func IsKernelPSD(kernel Kernel, x [][]float64, jitter float64) bool` | IsKernelPSD reports whether the Gram matrix of kernel over the data set x is numerically positive definite, tested by attempting a Cholesky… |
| `func KernelDiagonal(kernel Kernel, x [][]float64) []float64` | KernelDiagonal returns the vector of self-covariances kernel.Eval(x`i`, x`i`). |
| `func KernelMatrixRank(kernel Kernel, x [][]float64, tol float64) int` | KernelMatrixRank estimates the numerical rank of the Gram matrix of kernel over x by counting Cholesky pivots that exceed tol, adding progressively… |
| `func L1Distance(a, b []float64) float64` | L1Distance returns the Manhattan (L1) distance between vectors a and b. |
| `func LogDetFromCholesky(l Matrix) float64` | LogDetFromCholesky returns log(det A) computed from the lower-triangular Cholesky factor L of A, as 2·Σ log L`i``i`. |
| `func LogDetSPD(a Matrix) (float64, error)` | LogDetSPD returns log(det A) for a symmetric positive-definite matrix a. |
| `func LogMarginalLikelihood(kernel Kernel, x [][]float64, y []float64, noiseVar float64) (float64, error)` | LogMarginalLikelihood returns the exact log marginal likelihood of training data (x, y) under a zero-mean Gaussian process with the given kernel and… |
| `func MatVec(a Matrix, x []float64) []float64` | MatVec returns the matrix-vector product A·x. |
| `func MaxAbsDiff(a, b []float64) float64` | MaxAbsDiff returns the largest absolute difference between corresponding entries of a and b. |
| `func MaxAbsMatrixDiff(a, b Matrix) float64` | MaxAbsMatrixDiff returns the largest absolute difference between corresponding entries of a and b. |
| `func Norm(a []float64) float64` | Norm returns the Euclidean (L2) norm of vector a. |
| `func PosteriorMean(kernel Kernel, x [][]float64, y []float64, noiseVar float64, xstar [][]float64) ([]float64, error)` | PosteriorMean returns the Gaussian-process posterior predictive mean at the test inputs xstar for training data (x, y), the given kernel and… |
| `func PosteriorVariance(kernel Kernel, x [][]float64, y []float64, noiseVar float64, xstar [][]float64) ([]float64, error)` | PosteriorVariance returns the Gaussian-process posterior predictive variance at the test inputs xstar for training data (x, y). |
| `func Predict(kernel Kernel, x [][]float64, y []float64, noiseVar float64, xstar [][]float64) (mean, variance []float64, err error)` | Predict returns both the posterior predictive mean and variance at xstar for training data (x, y). |
| `func QuadraticForm(a Matrix, x []float64) float64` | QuadraticForm returns the scalar xᵀ·A·x. |
| `func RKHSDistance(f, g RKHSFunction) float64` | RKHSDistance returns the RKHS norm of the difference f-g, computed as √(⟨f,f⟩ - 2⟨f,g⟩ + ⟨g,g⟩). |
| `func RKHSInnerProduct(kernel Kernel, xa [][]float64, a []float64, xb [][]float64, b []float64) float64` | RKHSInnerProduct returns the RKHS inner product of two finite kernel expansions with weights a at centres xa and weights b at centres xb, all under… |
| `func RKHSNorm(kernel Kernel, x [][]float64, a []float64) float64` | RKHSNorm returns the RKHS norm of the finite kernel expansion with weights a at centres x under the given kernel. |
| `func RepresenterCoefficients(kernel Kernel, x [][]float64, y []float64, lambda float64) ([]float64, error)` | RepresenterCoefficients returns the coefficient vector w solving (K + lambda·I)·w = y for the kernel ridge problem, without wrapping it in an… |
| `func SampleMultivariateNormal(mean []float64, cov Matrix, jitter float64, rng *rand.Rand) ([]float64, error)` | SampleMultivariateNormal returns one draw from the multivariate normal distribution with the given mean vector and covariance matrix, using the… |
| `func SamplePrior(kernel Kernel, x [][]float64, jitter float64, sampleCount int, rng *rand.Rand) ([][]float64, error)` | SamplePrior draws sampleCount joint samples from the Gaussian-process prior with the given kernel evaluated at the inputs x, using the supplied… |
| `func SolveSPD(a Matrix, b []float64) ([]float64, error)` | SolveSPD solves the symmetric positive-definite system A·x = b by Cholesky factorisation. |
| `func SquaredDistance(a, b []float64) float64` | SquaredDistance returns the squared Euclidean distance between a and b. |
| `func SquaredNorm(a []float64) float64` | SquaredNorm returns the squared Euclidean norm of vector a. |
| `func StandardNormalPDF(x float64) float64` | StandardNormalPDF returns the density of the standard normal distribution at x. |
| `func VecAdd(a, b []float64) []float64` | VecAdd returns the element-wise sum a+b as a new slice. |
| `func VecCopy(a []float64) []float64` | VecCopy returns an independent copy of vector a. |
| `func VecScale(s float64, a []float64) []float64` | VecScale returns the vector a scaled by the scalar s as a new slice. |
| `func VecShift(c float64, a []float64) []float64` | VecShift returns a new slice with the scalar c added to every entry of a. |
| `func VecSub(a, b []float64) []float64` | VecSub returns the element-wise difference a-b as a new slice. |
| `func VectorsEqual(a, b []float64) bool` | VectorsEqual reports whether a and b have equal length and identical entries. |
| `func ZeroMean(x []float64) float64` | ZeroMean is the `MeanFunc` that is identically zero. |

### Types

| Type | What it is |
| --- | --- |
| `Constant` | Constant is the kernel k(x, y) = Value for all inputs. |
| `Cosine` | Cosine is the kernel k(x, y) = Variance·cos(2π·r/Period) with r = ‖x-y‖. |
| `GP` | GP is an exact Gaussian-process regression model. |
| `GammaExponential` | GammaExponential is the kernel k(x, y) = Variance·exp(-(r/ℓ)^Gamma) with r = ‖x-y‖ and ℓ = LengthScale. |
| `Kernel` | Kernel is a positive-semidefinite covariance function on pairs of real vectors of equal dimension. |
| `KernelFunc` | KernelFunc adapts an ordinary function to the `Kernel` interface. |
| `Linear` | Linear is the non-stationary dot-product kernel k(x, y) = Bias + Variance·(x-c)·(y-c) with c = Offset broadcast to every coordinate. |
| `Matern12` | Matern12 is the Matérn kernel with smoothness ν = 1/2, equal to the exponential (Ornstein-Uhlenbeck) kernel k(x, y) = Variance·exp(-r/LengthScale)… |
| `Matern32` | Matern32 is the Matérn kernel with smoothness ν = 3/2, k(x, y) = Variance·(1+√3·r/ℓ)·exp(-√3·r/ℓ) with r = ‖x-y‖ and ℓ =… |
| `Matern52` | Matern52 is the Matérn kernel with smoothness ν = 5/2, k(x, y) = Variance·(1+√5·r/ℓ+5r²/(3ℓ²))·exp(-√5·r/ℓ) with r = ‖x-y‖ and… |
| `Matrix` | Matrix is a dense, row-major real matrix stored as a slice of rows. |
| `MeanFunc` | MeanFunc is a deterministic prior mean function m(x) for a Gaussian process. |
| `Periodic` | Periodic is the exp-sine-squared kernel k(x, y) = Variance·exp(-2·sin²(π·r/Period)/ℓ²) with r = ‖x-y‖ and ℓ = LengthScale. |
| `Polynomial` | Polynomial is the kernel k(x, y) = (Alpha·x·y + Coef0)^Degree. |
| `ProductKernel` | ProductKernel is the pointwise product of its component kernels: k(x, y) = Π Kernels`i`.Eval(x, y). |
| `RBF` | RBF is the squared-exponential (radial-basis-function) kernel k(x, y) = Variance·exp(-‖x-y‖²/(2·LengthScale²)). |
| `RKHSFunction` | RKHSFunction represents an element of the reproducing-kernel Hilbert space (RKHS) of a positive-definite kernel as a finite expansion f(·) = Σ… |
| `RationalQuadratic` | RationalQuadratic is the kernel k(x, y) = Variance·(1+r²/(2·Alpha·ℓ²))^(-Alpha) with r = ‖x-y‖ and ℓ = LengthScale. |
| `ScaledKernel` | ScaledKernel multiplies a base kernel by a non-negative scalar: k(x, y) = Factor·Kernel.Eval(x, y). |
| `SumKernel` | SumKernel is the pointwise sum of its component kernels: k(x, y) = Σ Kernels`i`.Eval(x, y). |
| `WhiteNoise` | WhiteNoise is the kernel k(x, y) = Variance when x and y are identical vectors and 0 otherwise. |

<details>
<summary><code>Constant</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewConstant(value float64) Constant` | NewConstant returns a `Constant` kernel with the given value. |
| `func (k Constant) Eval(x, y []float64) float64` | Eval returns the constant value irrespective of x and y. |
| `func (k Constant) String() string` | String returns a human-readable description of the kernel. |

</details>

<details>
<summary><code>Cosine</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewCosine(variance, period float64) Cosine` | NewCosine returns a `Cosine` kernel with the given variance and period. |
| `func (k Cosine) Eval(x, y []float64) float64` | Eval returns the cosine covariance of x and y. |
| `func (k Cosine) String() string` | String returns a human-readable description of the kernel. |

</details>

<details>
<summary><code>GP</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewGP(kernel Kernel, noiseVariance float64) *GP` | NewGP returns a Gaussian-process model with the given kernel and observation noise variance, a zero prior mean and a small default jitter. |
| `func (g *GP) Alpha() []float64` | Alpha returns the weight vector (K+σ²I)^{-1}(y-m) computed during `GP.Fit`. |
| `func (g *GP) CholeskyFactor() Matrix` | CholeskyFactor returns the lower-triangular Cholesky factor of the noisy Gram matrix computed during `GP.Fit`. |
| `func (g *GP) Fit(x [][]float64, y []float64) error` | Fit conditions the model on training inputs x and targets y. |
| `func (g *GP) Fitted() bool` | Fitted reports whether the model has been conditioned on data. |
| `func (g *GP) LogMarginalLikelihood() (float64, error)` | LogMarginalLikelihood returns the exact log marginal likelihood log p(y \| X) = -½(y-m)ᵀα - Σ log L`i``i` - (n/2)·log(2π) of the fitted model,… |
| `func (g *GP) NumTraining() int` | NumTraining returns the number of training points the model was fitted on. |
| `func (g *GP) Predict(xstar [][]float64) (mean, variance []float64, err error)` | Predict returns both the posterior predictive mean and variance at the test inputs xstar. |
| `func (g *GP) PredictCovariance(xstar [][]float64) (Matrix, error)` | PredictCovariance returns the full posterior predictive covariance matrix over the test inputs xstar (excluding observation noise). |
| `func (g *GP) PredictMean(xstar [][]float64) ([]float64, error)` | PredictMean returns the posterior predictive mean at each test input in xstar. |
| `func (g *GP) PredictStd(xstar [][]float64) ([]float64, error)` | PredictStd returns the posterior predictive standard deviation (the square root of `GP.PredictVariance`) at each test input. |
| `func (g *GP) PredictVariance(xstar [][]float64) ([]float64, error)` | PredictVariance returns the posterior predictive variance at each test input in xstar (the marginal variance, ignoring observation noise). |
| `func (g *GP) Sample(xstar [][]float64, sampleCount int, rng *rand.Rand) ([][]float64, error)` | Sample draws sampleCount independent joint samples from the posterior over the test inputs xstar, using the supplied random source. |
| `func (g *GP) WithJitter(j float64) *GP` | WithJitter sets the diagonal jitter and returns the receiver for chaining. |
| `func (g *GP) WithMean(m MeanFunc) *GP` | WithMean sets the prior mean function and returns the receiver for chaining. |

</details>

<details>
<summary><code>GammaExponential</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewGammaExponential(variance, lengthScale, gamma float64) GammaExponential` | NewGammaExponential returns a `GammaExponential` kernel with the given variance, length scale and exponent gamma. |
| `func (k GammaExponential) Eval(x, y []float64) float64` | Eval returns the gamma-exponential covariance of x and y. |
| `func (k GammaExponential) String() string` | String returns a human-readable description of the kernel. |

</details>

<details>
<summary><code>Kernel</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AddKernels(kernels ...Kernel) Kernel` | AddKernels returns the pointwise sum of any number of kernels. |
| `func Exponentiate(base Kernel, p int) Kernel` | Exponentiate returns the kernel raised to a positive integer power p, implemented as a `ProductKernel` with p copies of base. |
| `func MaternKernel(nu, variance, lengthScale float64) Kernel` | MaternKernel returns the `Kernel` of the Matérn family for the smoothness parameter nu, which must be one of 0.5, 1.5 or 2.5. |
| `func MultiplyKernels(kernels ...Kernel) Kernel` | MultiplyKernels returns the pointwise product of any number of kernels. |
| `func Product(a, b Kernel) Kernel` | Product returns a `Kernel` equal to the pointwise product a·b. |
| `func Scale(factor float64, base Kernel) Kernel` | Scale returns a `Kernel` equal to the base kernel multiplied by the scalar factor. |
| `func Sum(a, b Kernel) Kernel` | Sum returns a `Kernel` equal to the pointwise sum a+b. |

</details>

<details>
<summary><code>KernelFunc</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (f KernelFunc) Eval(x, y []float64) float64` | Eval reports k(x, y) by calling the underlying function. |

</details>

<details>
<summary><code>Linear</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewLinear(bias, variance, offset float64) Linear` | NewLinear returns a `Linear` kernel with the given bias, slope variance and offset. |
| `func (k Linear) Eval(x, y []float64) float64` | Eval returns the linear covariance of x and y. |
| `func (k Linear) String() string` | String returns a human-readable description of the kernel. |

</details>

<details>
<summary><code>Matern12</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Exponential(variance, lengthScale float64) Matern12` | Exponential returns the exponential (Ornstein-Uhlenbeck) kernel, which is the Matérn kernel with ν = 1/2. |
| `func NewMatern12(variance, lengthScale float64) Matern12` | NewMatern12 returns a `Matern12` kernel with the given variance and length scale. |
| `func (k Matern12) Eval(x, y []float64) float64` | Eval returns the Matérn-1/2 covariance of x and y. |
| `func (k Matern12) String() string` | String returns a human-readable description of the kernel. |

</details>

<details>
<summary><code>Matern32</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewMatern32(variance, lengthScale float64) Matern32` | NewMatern32 returns a `Matern32` kernel with the given variance and length scale. |
| `func (k Matern32) Eval(x, y []float64) float64` | Eval returns the Matérn-3/2 covariance of x and y. |
| `func (k Matern32) String() string` | String returns a human-readable description of the kernel. |

</details>

<details>
<summary><code>Matern52</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewMatern52(variance, lengthScale float64) Matern52` | NewMatern52 returns a `Matern52` kernel with the given variance and length scale. |
| `func (k Matern52) Eval(x, y []float64) float64` | Eval returns the Matérn-5/2 covariance of x and y. |
| `func (k Matern52) String() string` | String returns a human-readable description of the kernel. |

</details>

<details>
<summary><code>Matrix</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AddJitter(k Matrix, jitter float64) Matrix` | AddJitter returns a copy of square matrix k with the small positive value jitter added to each diagonal entry, improving numerical conditioning for… |
| `func AddToDiagonal(m Matrix, c float64) Matrix` | AddToDiagonal returns a copy of square matrix m with the scalar c added to every main-diagonal entry. |
| `func Cholesky(a Matrix) (Matrix, error)` | Cholesky computes the lower-triangular Cholesky factor L of a symmetric positive-definite matrix a, so that a = L·Lᵀ. |
| `func CholeskySolveMatrix(l Matrix, b Matrix) (Matrix, error)` | CholeskySolveMatrix solves A·X = B column by column given the Cholesky factor L of A and a right-hand-side matrix B, returning X. |
| `func CrossGramMatrix(kernel Kernel, x1, x2 [][]float64) Matrix` | CrossGramMatrix returns the len(x1)-by-len(x2) matrix K with K`i``j` = kernel.Eval(x1`i`, x2`j`). |
| `func Diag(d []float64) Matrix` | Diag returns a square matrix with the entries of d on its main diagonal and zeros elsewhere. |
| `func GramMatrix(kernel Kernel, x [][]float64) Matrix` | GramMatrix returns the symmetric n-by-n Gram matrix K with K`i``j` = kernel.Eval(x`i`, x`j`) for the data set x (a slice of vectors). |
| `func Identity(n int) Matrix` | Identity returns the n-by-n identity matrix. |
| `func InvertSPD(a Matrix) (Matrix, error)` | InvertSPD returns the inverse of a symmetric positive-definite matrix a using its Cholesky factorisation. |
| `func MatAdd(a, b Matrix) Matrix` | MatAdd returns the element-wise sum a+b. |
| `func MatMul(a, b Matrix) Matrix` | MatMul returns the matrix product A·B. It panics if the number of columns of A differs from the number of rows of B. |
| `func MatScale(s float64, m Matrix) Matrix` | MatScale returns matrix m with every entry multiplied by the scalar s. |
| `func MatSub(a, b Matrix) Matrix` | MatSub returns the element-wise difference a-b. |
| `func NewMatrix(rows, cols int) Matrix` | NewMatrix returns a rows-by-cols matrix with all entries set to zero. |
| `func NoisyGramMatrix(kernel Kernel, x [][]float64, noiseVar float64) Matrix` | NoisyGramMatrix returns GramMatrix(kernel, x) with the scalar noiseVar added to every diagonal entry, that is K + noiseVar·I. This is the covariance… |
| `func OuterProduct(a, b []float64) Matrix` | OuterProduct returns the rank-one matrix a·bᵀ. |
| `func (m Matrix) At(i, j int) float64` | At returns the entry in row i and column j of m. |
| `func (m Matrix) Clone() Matrix` | Clone returns an independent deep copy of m. |
| `func (m Matrix) Cols() int` | Cols returns the number of columns of m, or zero if m has no rows. |
| `func (m Matrix) Diagonal() []float64` | Diagonal returns the main diagonal of m as a new slice. |
| `func (m Matrix) IsSquare() bool` | IsSquare reports whether m has equal row and column counts. |
| `func (m Matrix) IsSymmetric(tol float64) bool` | IsSymmetric reports whether m is square and symmetric within the absolute tolerance tol. |
| `func (m Matrix) Rows() int` | Rows returns the number of rows of m. |
| `func (m Matrix) Set(i, j int, v float64)` | Set assigns v to the entry in row i and column j of m. |
| `func (m Matrix) Trace() float64` | Trace returns the sum of the main-diagonal entries of square matrix m. |
| `func (m Matrix) Transpose() Matrix` | Transpose returns a new matrix that is the transpose of m. |

</details>

<details>
<summary><code>MeanFunc</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ConstantMean(c float64) MeanFunc` | ConstantMean returns a `MeanFunc` that is identically c. |

</details>

<details>
<summary><code>Periodic</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewPeriodic(variance, lengthScale, period float64) Periodic` | NewPeriodic returns a `Periodic` kernel with the given variance, length scale and period. |
| `func (k Periodic) Eval(x, y []float64) float64` | Eval returns the periodic covariance of x and y. |
| `func (k Periodic) String() string` | String returns a human-readable description of the kernel. |

</details>

<details>
<summary><code>Polynomial</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewPolynomial(alpha, coef0 float64, degree int) Polynomial` | NewPolynomial returns a `Polynomial` kernel with the given scale, constant offset and degree. |
| `func (k Polynomial) Eval(x, y []float64) float64` | Eval returns the polynomial covariance of x and y. |
| `func (k Polynomial) String() string` | String returns a human-readable description of the kernel. |

</details>

<details>
<summary><code>ProductKernel</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewProductKernel(kernels ...Kernel) ProductKernel` | NewProductKernel returns a `ProductKernel` over the given components. |
| `func (k ProductKernel) Eval(x, y []float64) float64` | Eval returns the product of the component covariances of x and y. |
| `func (k ProductKernel) String() string` | String returns a human-readable description of the kernel. |

</details>

<details>
<summary><code>RBF</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewRBF(variance, lengthScale float64) RBF` | NewRBF returns an `RBF` kernel with the given signal variance and length scale. |
| `func SquaredExponential(variance, lengthScale float64) RBF` | SquaredExponential is an alias constructor for the `RBF` kernel. |
| `func (k RBF) Eval(x, y []float64) float64` | Eval returns the squared-exponential covariance of x and y. |
| `func (k RBF) String() string` | String returns a human-readable description of the kernel. |

</details>

<details>
<summary><code>RKHSFunction</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func KernelRidgeRegression(kernel Kernel, x [][]float64, y []float64, lambda float64) (RKHSFunction, error)` | KernelRidgeRegression fits a kernel ridge regression model to training data (x, y) with regularisation strength lambda ≥ 0, returning the solution… |
| `func NewRKHSFunction(kernel Kernel, points [][]float64, weights []float64) RKHSFunction` | NewRKHSFunction returns an `RKHSFunction` with the given kernel, centres and weights. |
| `func (f RKHSFunction) Add(g RKHSFunction) RKHSFunction` | Add returns the RKHS function f+g by concatenating their expansions. |
| `func (f RKHSFunction) Eval(x []float64) float64` | Eval returns the value f(x) = Σ Weights`i`·Kernel(Points`i`, x). |
| `func (f RKHSFunction) EvalBatch(xs [][]float64) []float64` | EvalBatch returns f evaluated at each input in xs. |
| `func (f RKHSFunction) InnerProduct(g RKHSFunction) float64` | InnerProduct returns the RKHS inner product ⟨f, g⟩_H, which for finite expansions equals Σ_i Σ_j… |
| `func (f RKHSFunction) Norm() float64` | Norm returns the RKHS norm ‖f‖_H = √⟨f, f⟩_H. |
| `func (f RKHSFunction) Scale(c float64) RKHSFunction` | Scale returns the RKHS function c·f, sharing the same centres. |
| `func (f RKHSFunction) SquaredNorm() float64` | SquaredNorm returns the squared RKHS norm ⟨f, f⟩_H. |

</details>

<details>
<summary><code>RationalQuadratic</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewRationalQuadratic(variance, lengthScale, alpha float64) RationalQuadratic` | NewRationalQuadratic returns a `RationalQuadratic` kernel with the given variance, length scale and shape parameter alpha. |
| `func (k RationalQuadratic) Eval(x, y []float64) float64` | Eval returns the rational-quadratic covariance of x and y. |
| `func (k RationalQuadratic) String() string` | String returns a human-readable description of the kernel. |

</details>

<details>
<summary><code>ScaledKernel</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewScaledKernel(factor float64, base Kernel) ScaledKernel` | NewScaledKernel returns a `ScaledKernel` scaling base by factor. |
| `func (k ScaledKernel) Eval(x, y []float64) float64` | Eval returns the scaled covariance of x and y. |
| `func (k ScaledKernel) String() string` | String returns a human-readable description of the kernel. |

</details>

<details>
<summary><code>SumKernel</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewSumKernel(kernels ...Kernel) SumKernel` | NewSumKernel returns a `SumKernel` over the given components. |
| `func (k SumKernel) Eval(x, y []float64) float64` | Eval returns the sum of the component covariances of x and y. |
| `func (k SumKernel) String() string` | String returns a human-readable description of the kernel. |

</details>

<details>
<summary><code>WhiteNoise</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewWhiteNoise(variance float64) WhiteNoise` | NewWhiteNoise returns a `WhiteNoise` kernel with the given variance. |
| `func (k WhiteNoise) Eval(x, y []float64) float64` | Eval returns the noise variance when x equals y and 0 otherwise. |
| `func (k WhiteNoise) String() string` | String returns a human-readable description of the kernel. |

</details>

### Variables

`ErrDimensionMismatch`, `ErrEmpty`, `ErrNoKernel`, `ErrNonFinite`, `ErrNotPositiveDefinite`, `ErrSingular`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/gaussproc).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
