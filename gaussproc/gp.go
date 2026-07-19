package gaussproc

import (
	"math"
	"math/rand"
)

// MeanFunc is a deterministic prior mean function m(x) for a Gaussian process.
type MeanFunc func(x []float64) float64

// ZeroMean is the [MeanFunc] that is identically zero.
func ZeroMean(x []float64) float64 { return 0 }

// ConstantMean returns a [MeanFunc] that is identically c.
func ConstantMean(c float64) MeanFunc {
	return func(x []float64) float64 { return c }
}

// GP is an exact Gaussian-process regression model. Construct one with [NewGP],
// then supply training data with [GP.Fit] before calling the prediction and
// likelihood methods. The zero value is not usable; use [NewGP].
type GP struct {
	// Kernel is the covariance function.
	Kernel Kernel
	// NoiseVariance is the variance of the i.i.d. Gaussian observation noise.
	NoiseVariance float64
	// Mean is the prior mean function; if nil, [ZeroMean] is used.
	Mean MeanFunc
	// Jitter is added to the diagonal for numerical stability.
	Jitter float64

	x     [][]float64
	y     []float64
	l     Matrix    // Cholesky factor of (K + noise·I)
	alpha []float64 // (K + noise·I)^{-1} (y - m)
	fit   bool
}

// NewGP returns a Gaussian-process model with the given kernel and observation
// noise variance, a zero prior mean and a small default jitter.
func NewGP(kernel Kernel, noiseVariance float64) *GP {
	return &GP{
		Kernel:        kernel,
		NoiseVariance: noiseVariance,
		Mean:          ZeroMean,
		Jitter:        1e-10,
	}
}

// WithMean sets the prior mean function and returns the receiver for chaining.
func (g *GP) WithMean(m MeanFunc) *GP {
	g.Mean = m
	return g
}

// WithJitter sets the diagonal jitter and returns the receiver for chaining.
func (g *GP) WithJitter(j float64) *GP {
	g.Jitter = j
	return g
}

func (g *GP) meanFn() MeanFunc {
	if g.Mean == nil {
		return ZeroMean
	}
	return g.Mean
}

// Fit conditions the model on training inputs x and targets y. It factorises
// the noisy Gram matrix once so that subsequent predictions are cheap. It
// returns an error if the shapes are inconsistent or the covariance matrix is
// not positive definite.
func (g *GP) Fit(x [][]float64, y []float64) error {
	if len(x) == 0 {
		return ErrEmpty
	}
	if len(x) != len(y) {
		return ErrDimensionMismatch
	}
	k := NoisyGramMatrix(g.Kernel, x, g.NoiseVariance+g.Jitter)
	l, err := Cholesky(k)
	if err != nil {
		return err
	}
	m := g.meanFn()
	resid := make([]float64, len(y))
	for i := range y {
		resid[i] = y[i] - m(x[i])
	}
	alpha, err := CholeskySolve(l, resid)
	if err != nil {
		return err
	}
	g.x = x
	g.y = y
	g.l = l
	g.alpha = alpha
	g.fit = true
	return nil
}

// Fitted reports whether the model has been conditioned on data.
func (g *GP) Fitted() bool { return g.fit }

// NumTraining returns the number of training points the model was fitted on.
func (g *GP) NumTraining() int { return len(g.x) }

// CholeskyFactor returns the lower-triangular Cholesky factor of the noisy Gram
// matrix computed during [GP.Fit]. The returned matrix is a live reference;
// clone it before mutating.
func (g *GP) CholeskyFactor() Matrix { return g.l }

// Alpha returns the weight vector (K+σ²I)^{-1}(y-m) computed during [GP.Fit].
func (g *GP) Alpha() []float64 { return g.alpha }

// PredictMean returns the posterior predictive mean at each test input in
// xstar. It returns an error if the model has not been fitted.
func (g *GP) PredictMean(xstar [][]float64) ([]float64, error) {
	if !g.fit {
		return nil, ErrEmpty
	}
	m := g.meanFn()
	out := make([]float64, len(xstar))
	for i, xs := range xstar {
		var s float64
		for j, xt := range g.x {
			s += g.Kernel.Eval(xs, xt) * g.alpha[j]
		}
		out[i] = m(xs) + s
	}
	return out, nil
}

// PredictVariance returns the posterior predictive variance at each test input
// in xstar (the marginal variance, ignoring observation noise). Small negative
// values arising from round-off are clamped to zero.
func (g *GP) PredictVariance(xstar [][]float64) ([]float64, error) {
	if !g.fit {
		return nil, ErrEmpty
	}
	out := make([]float64, len(xstar))
	for i, xs := range xstar {
		ks := make([]float64, len(g.x))
		for j, xt := range g.x {
			ks[j] = g.Kernel.Eval(xs, xt)
		}
		v, err := ForwardSubstitution(g.l, ks)
		if err != nil {
			return nil, err
		}
		prior := g.Kernel.Eval(xs, xs)
		val := prior - SquaredNorm(v)
		if val < 0 {
			val = 0
		}
		out[i] = val
	}
	return out, nil
}

// Predict returns both the posterior predictive mean and variance at the test
// inputs xstar.
func (g *GP) Predict(xstar [][]float64) (mean, variance []float64, err error) {
	mean, err = g.PredictMean(xstar)
	if err != nil {
		return nil, nil, err
	}
	variance, err = g.PredictVariance(xstar)
	if err != nil {
		return nil, nil, err
	}
	return mean, variance, nil
}

// PredictStd returns the posterior predictive standard deviation (the square
// root of [GP.PredictVariance]) at each test input.
func (g *GP) PredictStd(xstar [][]float64) ([]float64, error) {
	v, err := g.PredictVariance(xstar)
	if err != nil {
		return nil, err
	}
	out := make([]float64, len(v))
	for i := range v {
		out[i] = math.Sqrt(v[i])
	}
	return out, nil
}

// PredictCovariance returns the full posterior predictive covariance matrix
// over the test inputs xstar (excluding observation noise).
func (g *GP) PredictCovariance(xstar [][]float64) (Matrix, error) {
	if !g.fit {
		return nil, ErrEmpty
	}
	kss := GramMatrix(g.Kernel, xstar)
	ksx := CrossGramMatrix(g.Kernel, xstar, g.x) // m x n
	// v = L^{-1} K_xs  (n x m), where K_xs = ksx^T
	kxs := ksx.Transpose() // n x m
	v, err := forwardSubstMatrix(g.l, kxs)
	if err != nil {
		return nil, err
	}
	// cov = Kss - v^T v
	vt := v.Transpose()
	reduction := MatMul(vt, v)
	cov := MatSub(kss, reduction)
	return cov, nil
}

// forwardSubstMatrix solves L·X = B column by column for lower-triangular L.
func forwardSubstMatrix(l Matrix, b Matrix) (Matrix, error) {
	n := l.Rows()
	if b.Rows() != n {
		return nil, ErrDimensionMismatch
	}
	cols := b.Cols()
	x := NewMatrix(n, cols)
	col := make([]float64, n)
	for j := 0; j < cols; j++ {
		for i := 0; i < n; i++ {
			col[i] = b[i][j]
		}
		sol, err := ForwardSubstitution(l, col)
		if err != nil {
			return nil, err
		}
		for i := 0; i < n; i++ {
			x[i][j] = sol[i]
		}
	}
	return x, nil
}

// LogMarginalLikelihood returns the exact log marginal likelihood
// log p(y | X) = -½(y-m)ᵀα - Σ log L[i][i] - (n/2)·log(2π) of the fitted
// model, where L is the Cholesky factor of the noisy Gram matrix. It returns an
// error if the model has not been fitted.
func (g *GP) LogMarginalLikelihood() (float64, error) {
	if !g.fit {
		return 0, ErrEmpty
	}
	m := g.meanFn()
	var data float64
	for i := range g.y {
		data += (g.y[i] - m(g.x[i])) * g.alpha[i]
	}
	var logdet float64
	for i := 0; i < g.l.Rows(); i++ {
		logdet += math.Log(g.l[i][i])
	}
	n := float64(len(g.y))
	return -0.5*data - logdet - 0.5*n*math.Log(2*math.Pi), nil
}

// Sample draws sampleCount independent joint samples from the posterior over
// the test inputs xstar, using the supplied random source. Each returned slice
// holds one function realisation evaluated at xstar. A small jitter is added to
// the covariance before its Cholesky factorisation for numerical stability.
func (g *GP) Sample(xstar [][]float64, sampleCount int, rng *rand.Rand) ([][]float64, error) {
	mean, err := g.PredictMean(xstar)
	if err != nil {
		return nil, err
	}
	cov, err := g.PredictCovariance(xstar)
	if err != nil {
		return nil, err
	}
	out := make([][]float64, sampleCount)
	for s := 0; s < sampleCount; s++ {
		sample, err := SampleMultivariateNormal(mean, cov, g.Jitter+1e-9, rng)
		if err != nil {
			return nil, err
		}
		out[s] = sample
	}
	return out, nil
}

// PosteriorMean returns the Gaussian-process posterior predictive mean at the
// test inputs xstar for training data (x, y), the given kernel and observation
// noise variance, using a zero prior mean. It is a stateless convenience
// wrapper around [GP].
func PosteriorMean(kernel Kernel, x [][]float64, y []float64, noiseVar float64, xstar [][]float64) ([]float64, error) {
	g := NewGP(kernel, noiseVar)
	if err := g.Fit(x, y); err != nil {
		return nil, err
	}
	return g.PredictMean(xstar)
}

// PosteriorVariance returns the Gaussian-process posterior predictive variance
// at the test inputs xstar for training data (x, y). It is a stateless
// convenience wrapper around [GP].
func PosteriorVariance(kernel Kernel, x [][]float64, y []float64, noiseVar float64, xstar [][]float64) ([]float64, error) {
	g := NewGP(kernel, noiseVar)
	if err := g.Fit(x, y); err != nil {
		return nil, err
	}
	return g.PredictVariance(xstar)
}

// Predict returns both the posterior predictive mean and variance at xstar for
// training data (x, y). It is a stateless convenience wrapper around [GP].
func Predict(kernel Kernel, x [][]float64, y []float64, noiseVar float64, xstar [][]float64) (mean, variance []float64, err error) {
	g := NewGP(kernel, noiseVar)
	if err := g.Fit(x, y); err != nil {
		return nil, nil, err
	}
	return g.Predict(xstar)
}

// LogMarginalLikelihood returns the exact log marginal likelihood of training
// data (x, y) under a zero-mean Gaussian process with the given kernel and
// noise variance. It is a stateless convenience wrapper around [GP].
func LogMarginalLikelihood(kernel Kernel, x [][]float64, y []float64, noiseVar float64) (float64, error) {
	g := NewGP(kernel, noiseVar)
	if err := g.Fit(x, y); err != nil {
		return 0, err
	}
	return g.LogMarginalLikelihood()
}

// SamplePrior draws sampleCount joint samples from the Gaussian-process prior
// with the given kernel evaluated at the inputs x, using the supplied random
// source. A stabilising jitter is added to the Gram matrix diagonal before
// factorisation.
func SamplePrior(kernel Kernel, x [][]float64, jitter float64, sampleCount int, rng *rand.Rand) ([][]float64, error) {
	cov := AddJitter(GramMatrix(kernel, x), jitter)
	mean := make([]float64, len(x))
	out := make([][]float64, sampleCount)
	for s := 0; s < sampleCount; s++ {
		sample, err := SampleMultivariateNormal(mean, cov, 0, rng)
		if err != nil {
			return nil, err
		}
		out[s] = sample
	}
	return out, nil
}

// SampleMultivariateNormal returns one draw from the multivariate normal
// distribution with the given mean vector and covariance matrix, using the
// supplied random source. The scalar jitter is added to the covariance diagonal
// before its Cholesky factorisation to guarantee positive definiteness.
func SampleMultivariateNormal(mean []float64, cov Matrix, jitter float64, rng *rand.Rand) ([]float64, error) {
	n := len(mean)
	if cov.Rows() != n || cov.Cols() != n {
		return nil, ErrDimensionMismatch
	}
	c := cov
	if jitter != 0 {
		c = AddToDiagonal(cov, jitter)
	}
	l, err := Cholesky(c)
	if err != nil {
		return nil, err
	}
	z := make([]float64, n)
	for i := range z {
		z[i] = rng.NormFloat64()
	}
	// sample = mean + L z
	out := make([]float64, n)
	for i := 0; i < n; i++ {
		s := mean[i]
		for j := 0; j <= i; j++ {
			s += l[i][j] * z[j]
		}
		out[i] = s
	}
	return out, nil
}

// StandardNormalPDF returns the density of the standard normal distribution at
// x.
func StandardNormalPDF(x float64) float64 {
	return math.Exp(-0.5*x*x) / math.Sqrt(2*math.Pi)
}

// GaussianLogPDF returns the log density of a univariate normal distribution
// with the given mean and variance evaluated at x. Variance must be positive.
func GaussianLogPDF(x, mean, variance float64) float64 {
	d := x - mean
	return -0.5*math.Log(2*math.Pi*variance) - d*d/(2*variance)
}
