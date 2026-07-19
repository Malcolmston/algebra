// Package gaussproc provides a self-contained toolkit for Gaussian process
// (GP) modelling, kernel methods and the small amount of numerical linear
// algebra they require, implemented using only the Go standard library.
//
// # Kernels
//
// A [Kernel] maps a pair of real vectors to a covariance value. The package
// ships the standard library of stationary and non-stationary kernels:
// the squared-exponential / radial-basis-function kernel ([RBF]); the Matérn
// family for smoothness parameters 1/2, 3/2 and 5/2 ([Matern12], [Matern32],
// [Matern52]); the [RationalQuadratic] kernel; the exponentially-periodic
// ([Periodic]) kernel; the [GammaExponential] and [Cosine] kernels; the
// non-stationary [Linear] and [Polynomial] dot-product kernels; and the
// [Constant] and [WhiteNoise] kernels. Kernels compose through an algebra of
// sums ([Sum], [SumKernel]), products ([Product], [ProductKernel]) and scalar
// scaling ([Scale], [ScaledKernel]); arbitrary functions become kernels via
// [KernelFunc].
//
// # Gram matrices
//
// [GramMatrix] evaluates a kernel over a data set to form the symmetric
// covariance matrix; [CrossGramMatrix] does the same for two data sets and
// [KernelDiagonal] returns only the self-covariances. [NoisyGramMatrix] adds
// homoscedastic observation noise to the diagonal, and [AddJitter] adds a small
// stabilising ridge.
//
// # Regression
//
// [GP] performs exact Gaussian-process regression. After [GP.Fit] the model can
// return the predictive mean and variance ([GP.Predict]), the full predictive
// covariance ([GP.PredictCovariance]), the exact log marginal likelihood
// ([GP.LogMarginalLikelihood]) and joint samples from the posterior
// ([GP.Sample]). Prior samples are drawn with [SamplePrior]. The stateless
// helpers [PosteriorMean], [PosteriorVariance] and [Predict] expose the same
// computation without constructing a model. A [MeanFunc] supplies a
// non-zero prior mean via [ConstantMean] or [ZeroMean].
//
// # RKHS utilities
//
// Every positive-definite kernel induces a reproducing-kernel Hilbert space. An
// [RKHSFunction] represents a finite kernel expansion; it can be evaluated
// ([RKHSFunction.Eval]) and combined through the Hilbert-space inner product
// ([RKHSFunction.InnerProduct]), norm ([RKHSFunction.Norm]) and distance
// ([RKHSDistance]). [KernelRidgeRegression] fits such a function by penalised
// least squares.
//
// # Linear algebra
//
// The regression and RKHS routines are built on a compact set of dense
// linear-algebra helpers: the symmetric positive-definite Cholesky
// factorisation ([Cholesky]) with triangular solves ([ForwardSubstitution],
// [BackSubstitution], [CholeskySolve]), log-determinants ([LogDetFromCholesky])
// and general vector and [Matrix] operations.
//
// All randomised routines take a *math/rand.Rand supplied by the caller, so
// results are fully reproducible from a seed; every other routine is
// deterministic and depends only on the standard library.
package gaussproc
