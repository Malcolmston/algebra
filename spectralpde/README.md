# spectralpde — implements spectral methods for the numerical solution of differential equations

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/spectralpde.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/spectralpde)

Package spectralpde implements spectral methods for the numerical solution of
differential equations.

The package provides Fourier and Chebyshev collocation machinery: spectral
differentiation matrices, Chebyshev-Gauss-Lobatto and periodic Fourier grids,
discrete cosine/sine and (fast) Fourier transforms, barycentric interpolation,
Clenshaw-Curtis and Gauss quadrature, and end-to-end solvers for the Poisson,
Helmholtz, heat and advection-diffusion equations in one and two space
dimensions.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/spectralpde"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AXPY(a float64, x, y []float64) []float64` | AXPY returns a*x + y. |
| `func AdvectionDiffusionSolveFourier(u0 []float64, c, nu, t, L float64) []float64` | AdvectionDiffusionSolveFourier evolves the periodic equation u_t + c*u_x = nu*u_xx on [0, L) exactly in Fourier space from u0 to time t. |
| `func AdvectionSolveFourier(u0 []float64, c, t, L float64) []float64` | AdvectionSolveFourier evolves the periodic linear advection equation u_t + c*u_x = 0 on [0, L) exactly to time t. |
| `func AffineMap(x, a, b float64) float64` | AffineMap returns the affine transformation of x from the reference interval [-1, 1] to the interval [a, b]. |
| `func AffineMapInverse(y, a, b float64) float64` | AffineMapInverse maps y from [a, b] back to the reference interval [-1, 1]. |
| `func ApplyFunc(f func(float64) float64, x []float64) []float64` | ApplyFunc evaluates f at each entry of x and returns the resulting vector. |
| `func BarycentricInterpolate(nodes, values, w []float64, x float64) float64` | BarycentricInterpolate evaluates the interpolating polynomial through (nodes, values) at x, using the barycentric weights w (see BarycentricWeights). |
| `func BarycentricInterpolateVec(nodes, values, w, xs []float64) []float64` | BarycentricInterpolateVec evaluates the barycentric interpolant at every point of xs. |
| `func BarycentricWeights(nodes []float64) []float64` | BarycentricWeights returns the barycentric interpolation weights w_j = 1 / prod_{k != j} (nodes`j` - nodes`k`) for an arbitrary set of distinct nodes. |
| `func BarycentricWeightsChebyshev(N int) []float64` | BarycentricWeightsChebyshev returns the closed-form barycentric weights for the N+1 Chebyshev-Gauss-Lobatto nodes: w_j = (-1)^j * delta_j with… |
| `func BuildChebyshevLaplacian1D(N int, a, b float64) [][]float64` | BuildChebyshevLaplacian1D returns the second-derivative (Laplacian in 1-D) Chebyshev collocation operator on [a, b], i.e. |
| `func BuildLaplacian2D(Nx, Ny int, ax, bx, ay, by float64) [][]float64` | BuildLaplacian2D returns the 2-D Chebyshev collocation Laplacian on the tensor-product grid of [ax, bx] x [ay, by] with Nx+1 nodes in x and Ny+1… |
| `func ChebyshevBarycentricInterpolate(values []float64, x float64) float64` | ChebyshevBarycentricInterpolate evaluates the interpolant of values sampled at the Chebyshev-Gauss-Lobatto nodes at x, using the closed-form… |
| `func ChebyshevCoefficients(values []float64) []float64` | ChebyshevCoefficients returns the discrete Chebyshev transform of the values f(x_j) sampled at the N+1 Chebyshev-Gauss-Lobatto nodes, where N =… |
| `func ChebyshevDiffMatrix(N int) [][]float64` | ChebyshevDiffMatrix returns the (N+1)-by-(N+1) Chebyshev collocation differentiation matrix D on the Chebyshev-Gauss-Lobatto nodes. |
| `func ChebyshevDiffMatrix2(N int) [][]float64` | ChebyshevDiffMatrix2 returns the second-order Chebyshev differentiation matrix D^2. |
| `func ChebyshevDiffMatrixInterval(N int, a, b float64) [][]float64` | ChebyshevDiffMatrixInterval returns the first-order Chebyshev differentiation matrix scaled for the physical interval [a, b]. |
| `func ChebyshevDiffMatrixOrder(N, m int) [][]float64` | ChebyshevDiffMatrixOrder returns the m-th order Chebyshev differentiation matrix, computed as the m-th matrix power of the first-order matrix. |
| `func ChebyshevDifferentiateCoeffs(coeffs []float64) []float64` | ChebyshevDifferentiateCoeffs returns the Chebyshev coefficients of the derivative of the series represented by coeffs. |
| `func ChebyshevEvalInterval(coeffs []float64, y, a, b float64) float64` | ChebyshevEvalInterval evaluates a reference-variable Chebyshev series at the physical point y in [a, b]. |
| `func ChebyshevFit(f func(float64) float64, N int) []float64` | ChebyshevFit samples f at the N+1 Chebyshev-Gauss-Lobatto nodes and returns the resulting Chebyshev coefficients of the interpolant. |
| `func ChebyshevFitInterval(f func(float64) float64, N int, a, b float64) []float64` | ChebyshevFitInterval samples f at the CGL nodes mapped to [a, b] and returns the Chebyshev coefficients of the interpolant in the reference variable. |
| `func ChebyshevGaussLobattoAngles(N int) []float64` | ChebyshevGaussLobattoAngles returns the angles theta_j = pi*j/N whose cosines are the Chebyshev-Gauss-Lobatto nodes. |
| `func ChebyshevGaussLobattoNodes(N int) []float64` | ChebyshevGaussLobattoNodes returns the N+1 Chebyshev-Gauss-Lobatto nodes x_j = cos(pi*j/N) on [-1, 1], ordered from +1 (j=0) down to -1 (j=N). |
| `func ChebyshevGaussLobattoNodesInterval(N int, a, b float64) []float64` | ChebyshevGaussLobattoNodesInterval returns the CGL nodes mapped to [a, b]. |
| `func ChebyshevGaussLobattoWeights(N int) []float64` | ChebyshevGaussLobattoWeights returns the quadrature weights associated with the Chebyshev-Gauss-Lobatto nodes for the weight function 1/sqrt(1-x^2). |
| `func ChebyshevGaussNodes(N int) []float64` | ChebyshevGaussNodes returns the N Chebyshev-Gauss nodes, i.e. |
| `func ChebyshevGaussWeights(N int) []float64` | ChebyshevGaussWeights returns the Chebyshev-Gauss quadrature weights, all equal to pi/N, for the weight function 1/sqrt(1-x^2). |
| `func ChebyshevIntegral(coeffs []float64) float64` | ChebyshevIntegral returns the definite integral over [-1, 1] of the Chebyshev series represented by coeffs. |
| `func ChebyshevIntegrateCoeffs(coeffs []float64) []float64` | ChebyshevIntegrateCoeffs returns the Chebyshev coefficients of an antiderivative of the series represented by coeffs. |
| `func ChebyshevInterpolantError(f func(float64) float64, N int, eval []float64) float64` | ChebyshevInterpolantError returns the maximum absolute difference between f and its degree-N Chebyshev interpolant, sampled at the given evaluation… |
| `func ChebyshevInterpolate(values []float64, x float64) float64` | ChebyshevInterpolate evaluates, at x, the polynomial interpolant of values sampled at the N+1 Chebyshev-Gauss-Lobatto nodes. |
| `func ChebyshevL2Projection(f func(float64) float64, n int) []float64` | ChebyshevL2Projection returns the first n+1 Chebyshev coefficients of the L2-orthogonal (Galerkin) projection of f with respect to the weight… |
| `func ChebyshevRoots(n int) []float64` | ChebyshevRoots returns the n roots of T_n (the Chebyshev-Gauss nodes). |
| `func ChebyshevT(n int, x float64) float64` | ChebyshevT evaluates the Chebyshev polynomial of the first kind T_n at x using the stable three-term recurrence. |
| `func ChebyshevTDerivative(n int, x float64) float64` | ChebyshevTDerivative evaluates the derivative of T_n at x, using the identity T_n'(x) = n*U_{n-1}(x). |
| `func ChebyshevTValues(n int, xs []float64) []float64` | ChebyshevTValues evaluates T_n at each point of xs. |
| `func ChebyshevTruncationError(coeffs []float64, keep int) float64` | ChebyshevTruncationError estimates the interpolation error of a Chebyshev series by the sum of absolute values of the tail coefficients beyond index… |
| `func ChebyshevU(n int, x float64) float64` | ChebyshevU evaluates the Chebyshev polynomial of the second kind U_n at x. |
| `func ChebyshevValuesFromCoeffs(coeffs []float64) []float64` | ChebyshevValuesFromCoeffs reconstructs the nodal values at the N+1 Chebyshev-Gauss-Lobatto nodes (N = len(coeffs)-1) from Chebyshev coefficients. |
| `func ChebyshevVandermonde(nodes []float64, degree int) [][]float64` | ChebyshevVandermonde returns the (len(nodes)) x (degree+1) matrix whose (i, j) entry is T_j(nodes`i`). |
| `func ClenshavEval(coeffs []float64, x float64) float64` | ClenshavEval evaluates the Chebyshev series sum_k coeffs`k`*T_k(x) using Clenshaw's recurrence. |
| `func ClenshawCurtisIntegrate(f func(float64) float64, N int, a, b float64) float64` | ClenshawCurtisIntegrate approximates the integral of f over [a, b] using an (N+1)-point Clenshaw-Curtis rule. |
| `func ClenshawCurtisNodes(N int) []float64` | ClenshawCurtisNodes returns the N+1 Clenshaw-Curtis nodes, i.e. |
| `func ClenshawCurtisWeights(N int) []float64` | ClenshawCurtisWeights returns the N+1 Clenshaw-Curtis quadrature weights on [-1, 1] for the weight function 1 (Trefethen's clencurt algorithm). |
| `func ClenshawEval(coeffs []float64, x float64) float64` | ClenshawEval evaluates the Chebyshev series sum_k coeffs`k`*T_k(x) using Clenshaw's backward recurrence, which is numerically stable on [-1, 1]. |
| `func ClenshawEvalDerivative(coeffs []float64, x float64) float64` | ClenshawEvalDerivative evaluates the derivative of the Chebyshev series sum_k coeffs`k`*T_k at x. |
| `func CoefficientDecayRate(coeffs []float64) float64` | CoefficientDecayRate estimates the exponential decay rate of the magnitudes of Chebyshev (or other) spectral coefficients, \|c_k\| ~ C*exp(-r*k). |
| `func ComplexAbs(x []complex128) []float64` | ComplexAbs returns the magnitudes of x. |
| `func ComplexImag(x []complex128) []float64` | ComplexImag returns the imaginary parts of x. |
| `func ComplexReal(x []complex128) []float64` | ComplexReal returns the real parts of x. |
| `func ConvergenceOrder(hs, errs []float64) float64` | ConvergenceOrder estimates the algebraic order p in err ~ C*h^p by a least-squares fit of log(err) against log(h). |
| `func CopyMatrix(a [][]float64) [][]float64` | CopyMatrix returns a deep copy of a [][]float64. |
| `func DCT1(x []float64) []float64` | DCT1 computes the type-I discrete cosine transform of the N+1 input samples (N = len(x)-1): X_k = 0.5*(x_0 + (-1)^k x_N) + sum_{n=1}^{N-1} x_n… |
| `func DFT(x []complex128) []complex128` | DFT computes the (unnormalized) discrete Fourier transform X_k = sum_j x_j exp(-2*pi*i*j*k/N) directly in O(N^2). |
| `func DST1(x []float64) []float64` | DST1 computes the type-I discrete sine transform of the N-1 interior samples (indexing n, k = 1..N-1): X_k = sum_{n=1}^{N-1} x_{n-1} sin(pi*n*k/N). |
| `func Determinant(a [][]float64) (float64, error)` | Determinant returns the determinant of a square matrix. |
| `func Diag(v []float64) [][]float64` | Diag returns a square matrix with the given vector on its main diagonal. |
| `func DifferentiateNodal(D [][]float64, u []float64) []float64` | DifferentiateNodal applies a differentiation matrix D to nodal values u, returning D*u. |
| `func DotProduct(x, y []float64) float64` | DotProduct returns the Euclidean inner product of x and y. |
| `func FFT(x []complex128) []complex128` | FFT computes the discrete Fourier transform. |
| `func FFTReal(x []float64) []complex128` | FFTReal computes the DFT of a real-valued signal. |
| `func FejerIntegrate(f func(float64) float64, n int, a, b float64) float64` | FejerIntegrate approximates the integral of f over [a, b] with an n-point Fejér first rule. |
| `func FejerNodes(n int) []float64` | FejerNodes returns the n nodes of Fejér's first quadrature rule, the Chebyshev-Gauss nodes x_k = cos((k+1/2)*pi/n). |
| `func FejerWeights(n int) []float64` | FejerWeights returns the n weights of Fejér's first quadrature rule on [-1, 1] for the weight function 1. |
| `func FourierCoefficients(values []float64) []complex128` | FourierCoefficients returns the complex Fourier coefficients c_k = (1/N) sum_j f_j exp(-i*k*x_j) of the values sampled on the periodic Fourier grid,… |
| `func FourierDiffMatrix(N int) [][]float64` | FourierDiffMatrix returns the N-by-N first-order Fourier differentiation matrix on the periodic grid x_j = 2*pi*j/N. It is implemented for even N,… |
| `func FourierDiffMatrix2(N int) [][]float64` | FourierDiffMatrix2 returns the N-by-N second-order Fourier differentiation matrix on the periodic grid, for even N. |
| `func FourierDiffMatrixOrder(N, m int) [][]float64` | FourierDiffMatrixOrder returns the m-th order Fourier differentiation matrix. |
| `func FourierDifferentiate(values []float64) []float64` | FourierDifferentiate returns the first derivative of the trigonometric interpolant at the grid nodes. |
| `func FourierDifferentiate2(values []float64) []float64` | FourierDifferentiate2 returns the second derivative of the trigonometric interpolant at the grid nodes. |
| `func FourierDifferentiateOrder(values []float64, m int) []float64` | FourierDifferentiateOrder returns the m-th derivative of the trigonometric interpolant, evaluated at the periodic grid nodes. |
| `func FourierInterpolate(values []float64, x float64) float64` | FourierInterpolate evaluates, at x, the band-limited trigonometric interpolant of values sampled at the N periodic nodes x_j = 2*pi*j/N. |
| `func FourierNodes(N int) []float64` | FourierNodes returns the N equispaced periodic nodes x_j = 2*pi*j/N on [0, 2*pi) for j = 0..N-1. |
| `func FourierNodesInterval(N int, a, b float64) []float64` | FourierNodesInterval returns N equispaced periodic nodes on [a, b). |
| `func FourierWavenumbers(N int) []float64` | FourierWavenumbers returns the signed wavenumbers associated with the length-N periodic FFT ordering: 0, 1, ..., N/2-1, -N/2, ..., -1 for even N. |
| `func FourierWeights(N int) []float64` | FourierWeights returns the trapezoidal (spectrally accurate for periodic functions) quadrature weights 2*pi/N on the Fourier grid. |
| `func FrobeniusNorm(a [][]float64) float64` | FrobeniusNorm returns the Frobenius norm of a. |
| `func GalerkinProjectChebyshev(f func(float64) float64, n int) []float64` | GalerkinProjectChebyshev is a convenience wrapper for ChebyshevL2Projection. |
| `func GalerkinProjectLegendre(f func(float64) float64, n int) []float64` | GalerkinProjectLegendre returns the first n+1 Legendre coefficients of the L2-orthogonal projection of f (see LegendreProjection). |
| `func GaussChebyshevIntegrate(f func(float64) float64, n int) float64` | GaussChebyshevIntegrate approximates the integral over [-1, 1] of f(x) with respect to the Chebyshev weight 1/sqrt(1-x^2), using an n-point… |
| `func GaussLegendreIntegrate(f func(float64) float64, n int, a, b float64) float64` | GaussLegendreIntegrate approximates the integral of f over [a, b] using an n-point Gauss-Legendre rule. |
| `func HeatSolveChebyshev(u0 func(float64) float64, nu, dt float64, steps, N int, a, b float64) (nodes, u []float64, err error)` | HeatSolveChebyshev integrates u_t = nu*u_xx on [a, b] with homogeneous Dirichlet boundary conditions from the initial function u0 to time steps*dt,… |
| `func HeatSolveFourier(u0 []float64, nu, t, L float64) []float64` | HeatSolveFourier evolves the periodic heat equation u_t = nu*u_xx on [0, L) exactly in Fourier space from the initial data u0 (sampled on the Fourier… |
| `func HeatStepChebyshevCN(uInterior []float64, nu, dt float64, N int, a, b float64) ([]float64, error)` | HeatStepChebyshevCN performs one Crank-Nicolson time step of size dt for the heat equation u_t = nu*u_xx on [a, b] with homogeneous Dirichlet… |
| `func Helmholtz2D(f, g func(x, y float64) float64, k2 float64, Nx, Ny int, ax, bx, ay, by float64) (xn, yn []float64, U [][]float64, err error)` | Helmholtz2D solves u_xx + u_yy + k2*u = f on [ax, bx] x [ay, by] with Dirichlet data g, by Chebyshev collocation. |
| `func HelmholtzSolve1D(f func(float64) float64, k2 float64, N int, a, b, ua, ub float64) (nodes, u []float64, err error)` | HelmholtzSolve1D solves the Dirichlet problem u” + k2*u = f on [a, b] with u(a) = ua and u(b) = ub, using Chebyshev collocation. |
| `func IDCT1(x []float64) []float64` | IDCT1 computes the inverse of DCT1, satisfying IDCT1(DCT1(x)) = x. |
| `func IDFT(x []complex128) []complex128` | IDFT computes the inverse discrete Fourier transform x_j = (1/N) sum_k X_k exp(2*pi*i*j*k/N) directly in O(N^2). |
| `func IDST1(x []float64) []float64` | IDST1 computes the inverse of DST1, satisfying IDST1(DST1(x)) = x. |
| `func IFFT(x []complex128) []complex128` | IFFT computes the inverse discrete Fourier transform, matching IDFT. |
| `func Identity(n int) [][]float64` | Identity returns the n-by-n identity matrix as [][]float64. |
| `func IntegrateFunction(f func(float64) float64, nodes, weights []float64) float64` | IntegrateFunction evaluates f at the given nodes and contracts with the weights. |
| `func IntegrateWeighted(weights, values []float64) float64` | IntegrateWeighted returns the dot product of quadrature weights and nodal function values. |
| `func InterpolationError(f func(float64) float64, nodes []float64, eval []float64) float64` | InterpolationError returns the maximum absolute difference between f and its interpolant (through (nodes, f(nodes))) sampled at the given evaluation… |
| `func IntervalScale(a, b float64) float64` | IntervalScale returns d/dx of the reference-to-physical map, i.e. |
| `func Inverse(a [][]float64) ([][]float64, error)` | Inverse returns the inverse of a square, non-singular matrix. |
| `func IsSpectrallyConverging(ns []int, errs []float64, rate float64) bool` | IsSpectrallyConverging reports whether the errors decay at least as fast as exp(-rate*n) over the tested resolutions, using SpectralConvergenceRate. |
| `func JacobiEigenSymmetric(a [][]float64) (vals []float64, vecs [][]float64, err error)` | JacobiEigenSymmetric computes all eigenvalues and eigenvectors of a real symmetric matrix using the cyclic Jacobi method. |
| `func Kron(a, b [][]float64) [][]float64` | Kron returns the Kronecker product of a and b. |
| `func L2Error(a, b []float64) float64` | L2Error returns the Euclidean norm of the difference a-b. |
| `func LUDecompose(a [][]float64) (lu [][]float64, piv []int, sign float64, err error)` | LUDecompose computes an in-place-style LU factorization with partial pivoting. |
| `func LUSolve(lu [][]float64, piv []int, b []float64) []float64` | LUSolve solves A*x = b given the LU factorization and pivot from LUDecompose. |
| `func LagrangeBasis(nodes []float64, j int, x float64) float64` | LagrangeBasis evaluates the j-th Lagrange cardinal polynomial for the given nodes at x. |
| `func LagrangeInterpolate(nodes, values []float64, x float64) float64` | LagrangeInterpolate evaluates the Lagrange interpolating polynomial through (nodes, values) at x directly from the cardinal-polynomial definition. |
| `func LegendreGaussLobattoNodes(N int) []float64` | LegendreGaussLobattoNodes returns just the N+1 Legendre-Gauss-Lobatto nodes. |
| `func LegendreGaussLobattoNodesWeights(N int) ([]float64, []float64)` | LegendreGaussLobattoNodesWeights returns the N+1 Legendre-Gauss-Lobatto nodes and weights on [-1, 1] (endpoints included). |
| `func LegendreGaussLobattoWeights(N int) []float64` | LegendreGaussLobattoWeights returns just the N+1 Legendre-Gauss-Lobatto weights. |
| `func LegendreGaussNodes(n int) []float64` | LegendreGaussNodes returns just the n Gauss-Legendre nodes on [-1, 1]. |
| `func LegendreGaussNodesWeights(n int) ([]float64, []float64)` | LegendreGaussNodesWeights returns the n Gauss-Legendre nodes and weights on [-1, 1], computed by Newton iteration on P_n. |
| `func LegendreGaussWeights(n int) []float64` | LegendreGaussWeights returns just the n Gauss-Legendre weights on [-1, 1]. |
| `func LegendreNormSquared(n int) float64` | LegendreNormSquared returns the squared L2([-1,1]) norm of P_n, 2/(2n+1). |
| `func LegendreP(n int, x float64) float64` | LegendreP evaluates the Legendre polynomial P_n at x via the three-term recurrence. |
| `func LegendrePDerivative(n int, x float64) float64` | LegendrePDerivative evaluates the derivative P_n'(x). |
| `func LegendrePSecondDerivative(n int, x float64) float64` | LegendrePSecondDerivative evaluates P_n”(x) using the Legendre differential equation. |
| `func LegendrePValues(n int, xs []float64) []float64` | LegendrePValues evaluates P_n at every point of xs. |
| `func LegendreProjection(f func(float64) float64, n int) []float64` | LegendreProjection returns the first n+1 Legendre coefficients of f, i.e. |
| `func LegendreSeriesEval(coeffs []float64, x float64) float64` | LegendreSeriesEval evaluates the Legendre series sum_k coeffs`k`*P_k at x. |
| `func LegendreVandermonde(nodes []float64, degree int) [][]float64` | LegendreVandermonde returns the len(nodes) x (degree+1) matrix with entries P_j(nodes`i`). |
| `func LinSpace(a, b float64, n int) []float64` | LinSpace returns n equally spaced points from a to b inclusive. |
| `func LinfError(a, b []float64) float64` | LinfError returns the l-infinity norm of the difference a-b. |
| `func MapFromInterval(nodes []float64, a, b float64) []float64` | MapFromInterval maps a slice of nodes on [a, b] back to [-1, 1]. |
| `func MapToInterval(nodes []float64, a, b float64) []float64` | MapToInterval maps a slice of reference nodes on [-1, 1] to [a, b]. |
| `func MatAdd(a, b [][]float64) [][]float64` | MatAdd returns the elementwise sum a+b. |
| `func MatMaxAbs(a [][]float64) float64` | MatMaxAbs returns the largest absolute entry of a. |
| `func MatMul(a, b [][]float64) [][]float64` | MatMul returns the matrix product a*b. |
| `func MatPow(a [][]float64, p int) [][]float64` | MatPow returns a raised to the non-negative integer power p (a must be square). |
| `func MatScale(a [][]float64, s float64) [][]float64` | MatScale returns s*a. |
| `func MatSub(a, b [][]float64) [][]float64` | MatSub returns the elementwise difference a-b. |
| `func MatVec(a [][]float64, x []float64) []float64` | MatVec returns the matrix-vector product a*x. |
| `func MaxAbs(x []float64) float64` | MaxAbs returns the largest absolute value in x. |
| `func Mean(x []float64) float64` | Mean returns the arithmetic mean of x. |
| `func NewtonDividedDifferences(nodes, values []float64) []float64` | NewtonDividedDifferences returns the Newton divided-difference coefficients for the data (nodes, values). |
| `func NewtonEval(coef, nodes []float64, x float64) float64` | NewtonEval evaluates the Newton form of the interpolating polynomial with the given divided-difference coefficients and nodes at x. |
| `func Norm1(x []float64) float64` | Norm1 returns the l1 norm of x. |
| `func Norm2(x []float64) float64` | Norm2 returns the Euclidean (l2) norm of x. |
| `func NormInf(x []float64) float64` | NormInf returns the l-infinity (maximum absolute) norm of x. |
| `func NormP(x []float64, p float64) float64` | NormP returns the discrete l^p norm of x for p >= 1. |
| `func Poisson2D(f, g func(x, y float64) float64, Nx, Ny int, ax, bx, ay, by float64) (xn, yn []float64, U [][]float64, err error)` | Poisson2D solves the Dirichlet problem u_xx + u_yy = f on the rectangle [ax, bx] x [ay, by] with boundary values given by g(x, y). |
| `func Poisson2DFourier(f func(x, y float64) float64, Nx, Ny int, Lx, Ly float64) (xn, yn []float64, U [][]float64)` | Poisson2DFourier solves the doubly periodic problem u_xx + u_yy = f on [0, Lx) x [0, Ly) with a Fourier spectral method. |
| `func PoissonSolve1D(f func(float64) float64, N int, a, b, ua, ub float64) (nodes, u []float64, err error)` | PoissonSolve1D solves the Dirichlet problem u” = f on [a, b] with u(a) = ua and u(b) = ub, using an (N+1)-node Chebyshev collocation method. |
| `func PoissonSolve1DFourier(f func(float64) float64, N int, L float64) (nodes, u []float64)` | PoissonSolve1DFourier solves the periodic problem u” = f on [0, L) using a Fourier spectral method. |
| `func PolynomialInterpolate(nodes, values []float64, x float64) float64` | PolynomialInterpolate builds barycentric weights for the given nodes and evaluates the interpolant at x. |
| `func RMSError(a, b []float64) float64` | RMSError returns the root-mean-square difference of a and b. |
| `func RealToComplex(x []float64) []complex128` | RealToComplex returns a complex slice whose real parts are x and imaginary parts are zero. |
| `func RelativeL2Error(a, b []float64) float64` | RelativeL2Error returns \|\|a-b\|\|_2 / \|\|b\|\|_2, guarding against division by zero. |
| `func Reverse(x []float64) []float64` | Reverse returns x with its entries in reverse order. |
| `func SimpsonIntegrate(f func(float64) float64, n int, a, b float64) float64` | SimpsonIntegrate approximates the integral of f over [a, b] with the composite Simpson rule on n subintervals (n must be even and >= 2). |
| `func SolveLinearSystem(a [][]float64, b []float64) ([]float64, error)` | SolveLinearSystem solves A*x = b for a square, non-singular A. |
| `func SolveMatrix(a, b [][]float64) ([][]float64, error)` | SolveMatrix solves A*X = B for X, where B (and X) have multiple columns. |
| `func SpectralConvergenceRate(ns []int, errs []float64) float64` | SpectralConvergenceRate estimates the exponential decay rate r in err ~ C*exp(-r*n) by a least-squares fit of log(err) against n. |
| `func SpectralDerivativeChebyshev(values []float64, a, b float64) []float64` | SpectralDerivativeChebyshev returns the nodal derivative values of a function sampled at the Chebyshev-Gauss-Lobatto nodes on [a, b]. |
| `func Sum(x []float64) float64` | Sum returns the sum of the entries of x. |
| `func Trace(a [][]float64) float64` | Trace returns the sum of the diagonal entries of a square matrix. |
| `func Transpose(a [][]float64) [][]float64` | Transpose returns the transpose of a. |
| `func TrapezoidalIntegrate(f func(float64) float64, n int, a, b float64) float64` | TrapezoidalIntegrate approximates the integral of f over [a, b] with the composite trapezoidal rule on n subintervals. |
| `func TrapezoidalWeights(n int, a, b float64) []float64` | TrapezoidalWeights returns the composite trapezoidal weights for n+1 equally spaced points on [a, b]. |
| `func UniformNodes(n int, a, b float64) []float64` | UniformNodes returns n+1 equally spaced nodes on [a, b]. |
| `func VectorAdd(x, y []float64) []float64` | VectorAdd returns x+y. |
| `func VectorCopy(x []float64) []float64` | VectorCopy returns a copy of x. |
| `func VectorFill(n int, v float64) []float64` | VectorFill returns a vector of length n with every entry set to v. |
| `func VectorHadamard(x, y []float64) []float64` | VectorHadamard returns the elementwise product of x and y. |
| `func VectorScale(x []float64, s float64) []float64` | VectorScale returns s*x. |
| `func VectorSub(x, y []float64) []float64` | VectorSub returns x-y. |
| `func VectorZeros(n int) []float64` | VectorZeros returns a zero vector of length n. |
| `func Zeros(r, c int) [][]float64` | Zeros returns an r-by-c matrix of zeros as a plain [][]float64. |

### Types

| Type | What it is |
| --- | --- |
| `Matrix` | Matrix is a dense, row-major real matrix used by the linear-algebra helpers. |

<details>
<summary><code>Matrix</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewMatrix(r, c int) *Matrix` | NewMatrix allocates an r-by-c zero matrix. |
| `func NewMatrixFrom(data [][]float64) *Matrix` | NewMatrixFrom wraps an existing row-major slice as a Matrix. |
| `func (m *Matrix) At(i, j int) float64` | At returns the entry at row i, column j. |
| `func (m *Matrix) Clone() *Matrix` | Clone returns a deep copy of the matrix. |
| `func (m *Matrix) Set(i, j int, v float64)` | Set assigns v to the entry at row i, column j. |

</details>

### Variables

`ErrDimensionMismatch`, `ErrInvalidArgument`, `ErrSingularMatrix`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/spectralpde).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
