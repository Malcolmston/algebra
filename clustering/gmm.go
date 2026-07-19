package clustering

import (
	"math"
	"math/rand"
)

// CovarianceType selects the structure imposed on the component covariance
// matrices of a Gaussian mixture model.
type CovarianceType int

const (
	// FullCovariance allows each component its own full covariance matrix.
	FullCovariance CovarianceType = iota
	// DiagonalCovariance restricts each component to a diagonal covariance.
	DiagonalCovariance
	// SphericalCovariance restricts each component to an isotropic (scalar times
	// identity) covariance.
	SphericalCovariance
)

// GaussianComponent describes a single multivariate normal component of a
// mixture: its mixing weight, mean and covariance matrix.
type GaussianComponent struct {
	Weight     float64
	Mean       []float64
	Covariance [][]float64
}

// GMM is a fitted Gaussian mixture model.
type GMM struct {
	// Components holds the K fitted components.
	Components []GaussianComponent
	// LogLikelihood is the final data log-likelihood.
	LogLikelihood float64
	// Iterations is the number of EM iterations performed.
	Iterations int
	// Converged reports whether EM stopped due to log-likelihood convergence.
	Converged bool
	// CovType is the covariance structure used during fitting.
	CovType CovarianceType

	invCov [][][]float64
	logDet []float64
}

// GMMOptions configures Gaussian-mixture EM fitting.
type GMMOptions struct {
	// MaxIter is the maximum number of EM iterations (default 100).
	MaxIter int
	// Tol is the log-likelihood improvement tolerance for convergence (default
	// 1e-4).
	Tol float64
	// Reg is a small value added to the covariance diagonal for numerical
	// stability (default 1e-6).
	Reg float64
	// CovType selects the covariance structure (default FullCovariance).
	CovType CovarianceType
	// Rand is the random source used for k-means++ initialisation.
	Rand *rand.Rand
}

func (o GMMOptions) withDefaults() GMMOptions {
	if o.MaxIter <= 0 {
		o.MaxIter = 100
	}
	if o.Tol <= 0 {
		o.Tol = 1e-4
	}
	if o.Reg <= 0 {
		o.Reg = 1e-6
	}
	if o.Rand == nil {
		o.Rand = rand.New(rand.NewSource(1))
	}
	return o
}

// FitGMM fits a Gaussian mixture model with k components to data using the EM
// algorithm and default options.
func FitGMM(data [][]float64, k int) (*GMM, error) {
	return FitGMMWithOptions(data, k, GMMOptions{})
}

// FitGMMWithOptions fits a Gaussian mixture model with k components using the
// supplied options. Means are initialised with k-means, and the EM iterations
// alternate expectation (responsibilities) and maximisation (parameter update)
// steps until the log-likelihood converges or the iteration cap is reached.
func FitGMMWithOptions(data [][]float64, k int, opts GMMOptions) (*GMM, error) {
	if len(data) == 0 {
		return nil, ErrEmptyData
	}
	if k <= 0 || k > len(data) {
		return nil, ErrInvalidK
	}
	opts = opts.withDefaults()
	n := len(data)
	dim := len(data[0])

	// Initialise means via k-means for stable starting points.
	km, err := KMeansWithOptions(data, k, KMeansOptions{NInit: 3, Rand: opts.Rand})
	if err != nil {
		return nil, err
	}
	comps := make([]GaussianComponent, k)
	counts := make([]int, k)
	for _, l := range km.Labels {
		counts[l]++
	}
	globalCov := CovarianceMatrix(data)
	if globalCov == nil {
		globalCov = Identity(dim)
	}
	for c := 0; c < k; c++ {
		comps[c].Weight = float64(counts[c]) / float64(n)
		if comps[c].Weight == 0 {
			comps[c].Weight = 1e-3
		}
		comps[c].Mean = CloneVector(km.Centroids[c])
		comps[c].Covariance = CloneMatrix(globalCov)
		regularizeCov(comps[c].Covariance, opts.Reg)
	}

	g := &GMM{Components: comps, CovType: opts.CovType}
	prevLL := math.Inf(-1)
	resp := make([][]float64, n)
	for i := range resp {
		resp[i] = make([]float64, k)
	}
	iter := 0
	for iter = 0; iter < opts.MaxIter; iter++ {
		if err := g.precompute(); err != nil {
			return nil, err
		}
		ll := g.eStep(data, resp)
		g.mStep(data, resp, opts)
		if ll-prevLL <= opts.Tol && iter > 0 {
			prevLL = ll
			iter++
			g.Converged = true
			break
		}
		prevLL = ll
	}
	if err := g.precompute(); err != nil {
		return nil, err
	}
	g.LogLikelihood = g.eStep(data, resp)
	g.Iterations = iter
	return g, nil
}

func regularizeCov(cov [][]float64, reg float64) {
	for i := range cov {
		cov[i][i] += reg
	}
}

// precompute caches the inverse and log-determinant of each component
// covariance.
func (g *GMM) precompute() error {
	k := len(g.Components)
	g.invCov = make([][][]float64, k)
	g.logDet = make([]float64, k)
	for c := 0; c < k; c++ {
		cov := g.Components[c].Covariance
		inv, err := Inverse(cov)
		if err != nil {
			// Add jitter and retry.
			jit := CloneMatrix(cov)
			regularizeCov(jit, 1e-6)
			inv, err = Inverse(jit)
			if err != nil {
				return err
			}
			cov = jit
			g.Components[c].Covariance = jit
		}
		det, derr := Determinant(cov)
		if derr != nil {
			return derr
		}
		if det <= 0 {
			det = 1e-300
		}
		g.invCov[c] = inv
		g.logDet[c] = math.Log(det)
	}
	return nil
}

// logGaussian returns the log-density of the multivariate normal for component c
// evaluated at x.
func (g *GMM) logGaussian(x []float64, c int) float64 {
	dim := len(x)
	md := mahalanobisSq(x, g.Components[c].Mean, g.invCov[c])
	return -0.5 * (float64(dim)*math.Log(2*math.Pi) + g.logDet[c] + md)
}

func mahalanobisSq(x, mean []float64, invCov [][]float64) float64 {
	dim := len(x)
	diff := make([]float64, dim)
	for i := 0; i < dim; i++ {
		diff[i] = x[i] - mean[i]
	}
	var s float64
	for i := 0; i < dim; i++ {
		var row float64
		for j := 0; j < dim; j++ {
			row += invCov[i][j] * diff[j]
		}
		s += diff[i] * row
	}
	if s < 0 {
		s = 0
	}
	return s
}

// eStep computes responsibilities and returns the total data log-likelihood.
func (g *GMM) eStep(data [][]float64, resp [][]float64) float64 {
	k := len(g.Components)
	var totalLL float64
	logw := make([]float64, k)
	for c := 0; c < k; c++ {
		logw[c] = math.Log(math.Max(g.Components[c].Weight, 1e-300))
	}
	for i, x := range data {
		logProbs := make([]float64, k)
		maxLP := math.Inf(-1)
		for c := 0; c < k; c++ {
			logProbs[c] = logw[c] + g.logGaussian(x, c)
			if logProbs[c] > maxLP {
				maxLP = logProbs[c]
			}
		}
		var sum float64
		for c := 0; c < k; c++ {
			sum += math.Exp(logProbs[c] - maxLP)
		}
		logSum := maxLP + math.Log(sum)
		totalLL += logSum
		for c := 0; c < k; c++ {
			resp[i][c] = math.Exp(logProbs[c] - logSum)
		}
	}
	return totalLL
}

// mStep updates weights, means and covariances from the responsibilities.
func (g *GMM) mStep(data [][]float64, resp [][]float64, opts GMMOptions) {
	n := len(data)
	dim := len(data[0])
	k := len(g.Components)
	nk := make([]float64, k)
	for c := 0; c < k; c++ {
		for i := 0; i < n; i++ {
			nk[c] += resp[i][c]
		}
	}
	for c := 0; c < k; c++ {
		if nk[c] < 1e-300 {
			nk[c] = 1e-300
		}
		// Weight.
		g.Components[c].Weight = nk[c] / float64(n)
		// Mean.
		mean := make([]float64, dim)
		for i := 0; i < n; i++ {
			r := resp[i][c]
			for j := 0; j < dim; j++ {
				mean[j] += r * data[i][j]
			}
		}
		for j := 0; j < dim; j++ {
			mean[j] /= nk[c]
		}
		g.Components[c].Mean = mean
		// Covariance.
		cov := Zeros(dim, dim)
		for i := 0; i < n; i++ {
			r := resp[i][c]
			diff := make([]float64, dim)
			for j := 0; j < dim; j++ {
				diff[j] = data[i][j] - mean[j]
			}
			for a := 0; a < dim; a++ {
				for b := a; b < dim; b++ {
					cov[a][b] += r * diff[a] * diff[b]
				}
			}
		}
		for a := 0; a < dim; a++ {
			for b := a; b < dim; b++ {
				cov[a][b] /= nk[c]
				cov[b][a] = cov[a][b]
			}
		}
		applyCovType(cov, opts.CovType)
		regularizeCov(cov, opts.Reg)
		g.Components[c].Covariance = cov
	}
}

func applyCovType(cov [][]float64, t CovarianceType) {
	dim := len(cov)
	switch t {
	case DiagonalCovariance:
		for a := 0; a < dim; a++ {
			for b := 0; b < dim; b++ {
				if a != b {
					cov[a][b] = 0
				}
			}
		}
	case SphericalCovariance:
		var avg float64
		for a := 0; a < dim; a++ {
			avg += cov[a][a]
		}
		avg /= float64(dim)
		for a := 0; a < dim; a++ {
			for b := 0; b < dim; b++ {
				if a == b {
					cov[a][b] = avg
				} else {
					cov[a][b] = 0
				}
			}
		}
	}
}

// PredictProba returns the posterior probability (responsibility) of each
// component for every row of newData.
func (g *GMM) PredictProba(newData [][]float64) [][]float64 {
	if g.invCov == nil {
		_ = g.precompute()
	}
	k := len(g.Components)
	out := make([][]float64, len(newData))
	logw := make([]float64, k)
	for c := 0; c < k; c++ {
		logw[c] = math.Log(math.Max(g.Components[c].Weight, 1e-300))
	}
	for i, x := range newData {
		out[i] = make([]float64, k)
		logProbs := make([]float64, k)
		maxLP := math.Inf(-1)
		for c := 0; c < k; c++ {
			logProbs[c] = logw[c] + g.logGaussian(x, c)
			if logProbs[c] > maxLP {
				maxLP = logProbs[c]
			}
		}
		var sum float64
		for c := 0; c < k; c++ {
			sum += math.Exp(logProbs[c] - maxLP)
		}
		logSum := maxLP + math.Log(sum)
		for c := 0; c < k; c++ {
			out[i][c] = math.Exp(logProbs[c] - logSum)
		}
	}
	return out
}

// Predict returns the most probable component for each row of newData.
func (g *GMM) Predict(newData [][]float64) []int {
	proba := g.PredictProba(newData)
	labels := make([]int, len(newData))
	for i, p := range proba {
		best := 0
		for c := 1; c < len(p); c++ {
			if p[c] > p[best] {
				best = c
			}
		}
		labels[i] = best
	}
	return labels
}

// ScoreSamples returns the per-sample log-likelihood of newData under the fitted
// mixture.
func (g *GMM) ScoreSamples(newData [][]float64) []float64 {
	if g.invCov == nil {
		_ = g.precompute()
	}
	k := len(g.Components)
	logw := make([]float64, k)
	for c := 0; c < k; c++ {
		logw[c] = math.Log(math.Max(g.Components[c].Weight, 1e-300))
	}
	out := make([]float64, len(newData))
	for i, x := range newData {
		maxLP := math.Inf(-1)
		lp := make([]float64, k)
		for c := 0; c < k; c++ {
			lp[c] = logw[c] + g.logGaussian(x, c)
			if lp[c] > maxLP {
				maxLP = lp[c]
			}
		}
		var sum float64
		for c := 0; c < k; c++ {
			sum += math.Exp(lp[c] - maxLP)
		}
		out[i] = maxLP + math.Log(sum)
	}
	return out
}

// Score returns the mean per-sample log-likelihood of newData.
func (g *GMM) Score(newData [][]float64) float64 {
	s := g.ScoreSamples(newData)
	if len(s) == 0 {
		return 0
	}
	var sum float64
	for _, v := range s {
		sum += v
	}
	return sum / float64(len(s))
}

// NumParameters returns the number of free parameters of the fitted mixture,
// used by the AIC and BIC criteria.
func (g *GMM) NumParameters() int {
	k := len(g.Components)
	if k == 0 {
		return 0
	}
	dim := len(g.Components[0].Mean)
	meanParams := k * dim
	weightParams := k - 1
	var covParams int
	switch g.CovType {
	case DiagonalCovariance:
		covParams = k * dim
	case SphericalCovariance:
		covParams = k
	default:
		covParams = k * dim * (dim + 1) / 2
	}
	return meanParams + weightParams + covParams
}

// AIC returns the Akaike information criterion of the fitted mixture on data.
func (g *GMM) AIC(data [][]float64) float64 {
	ll := g.LogLikelihood
	if len(data) > 0 {
		var sum float64
		for _, v := range g.ScoreSamples(data) {
			sum += v
		}
		ll = sum
	}
	return 2*float64(g.NumParameters()) - 2*ll
}

// BIC returns the Bayesian information criterion of the fitted mixture on data.
func (g *GMM) BIC(data [][]float64) float64 {
	n := len(data)
	var ll float64
	for _, v := range g.ScoreSamples(data) {
		ll += v
	}
	return float64(g.NumParameters())*math.Log(float64(n)) - 2*ll
}

// MultivariateNormalPDF returns the density of a multivariate normal
// distribution with the given mean and covariance evaluated at x.
func MultivariateNormalPDF(x, mean []float64, cov [][]float64) (float64, error) {
	inv, err := Inverse(cov)
	if err != nil {
		return 0, err
	}
	det, err := Determinant(cov)
	if err != nil {
		return 0, err
	}
	if det <= 0 {
		return 0, ErrSingularMatrix
	}
	dim := len(x)
	md := mahalanobisSq(x, mean, inv)
	norm := math.Pow(2*math.Pi, -float64(dim)/2) * math.Pow(det, -0.5)
	return norm * math.Exp(-0.5*md), nil
}
