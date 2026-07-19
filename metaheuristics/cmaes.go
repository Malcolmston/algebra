package metaheuristics

import (
	"math"
	"sort"
)

// CMAESConfig configures [RunCMAES], the covariance matrix adaptation evolution
// strategy.
type CMAESConfig struct {
	// Bounds is the search box; the initial mean defaults to its center.
	Bounds Bounds
	// Mean is the initial distribution mean; if nil, the box center is used.
	Mean []float64
	// Sigma is the initial global step size (standard deviation). If zero, a
	// value of 0.3 times the mean box width is used.
	Sigma float64
	// Lambda is the population (offspring) size; if zero, the default
	// 4 + floor(3 ln N) is used.
	Lambda int
	// MaxIterations is the number of generations.
	MaxIterations int
	// Tolerance stops the run early when sigma times the largest axis length
	// falls below it. Zero disables the check.
	Tolerance float64
	// RecordHistory enables per-generation best-value recording.
	RecordHistory bool
}

// DefaultCMAESConfig returns a reasonable configuration for the given box.
func DefaultCMAESConfig(b Bounds) CMAESConfig {
	return CMAESConfig{
		Bounds:        b,
		Sigma:         0,
		Lambda:        0,
		MaxIterations: 300,
		Tolerance:     1e-12,
	}
}

func (c CMAESConfig) validate() error {
	if !c.Bounds.Valid() {
		return ErrEmptyBounds
	}
	if c.MaxIterations <= 0 {
		return ErrInvalidConfig
	}
	return nil
}

// JacobiEigen computes the eigenvalues and orthonormal eigenvectors of the
// symmetric matrix a using the cyclic Jacobi rotation algorithm. It returns the
// eigenvalues and a matrix whose columns are the corresponding eigenvectors.
// The input is not modified.
func JacobiEigen(a [][]float64) (values []float64, vectors [][]float64) {
	n := len(a)
	// Work on a copy.
	m := make([][]float64, n)
	for i := range m {
		m[i] = make([]float64, n)
		copy(m[i], a[i])
	}
	v := identityMatrix(n)
	for sweep := 0; sweep < 100; sweep++ {
		off := 0.0
		for p := 0; p < n; p++ {
			for q := p + 1; q < n; q++ {
				off += m[p][q] * m[p][q]
			}
		}
		if off < 1e-30 {
			break
		}
		for p := 0; p < n; p++ {
			for q := p + 1; q < n; q++ {
				if math.Abs(m[p][q]) < 1e-300 {
					continue
				}
				theta := (m[q][q] - m[p][p]) / (2 * m[p][q])
				t := sign(theta) / (math.Abs(theta) + math.Sqrt(theta*theta+1))
				if theta == 0 {
					t = 1
				}
				cs := 1 / math.Sqrt(t*t+1)
				sn := t * cs
				for k := 0; k < n; k++ {
					mkp := m[k][p]
					mkq := m[k][q]
					m[k][p] = cs*mkp - sn*mkq
					m[k][q] = sn*mkp + cs*mkq
				}
				for k := 0; k < n; k++ {
					mpk := m[p][k]
					mqk := m[q][k]
					m[p][k] = cs*mpk - sn*mqk
					m[q][k] = sn*mpk + cs*mqk
				}
				for k := 0; k < n; k++ {
					vkp := v[k][p]
					vkq := v[k][q]
					v[k][p] = cs*vkp - sn*vkq
					v[k][q] = sn*vkp + cs*vkq
				}
			}
		}
	}
	values = make([]float64, n)
	for i := 0; i < n; i++ {
		values[i] = m[i][i]
	}
	return values, v
}

func identityMatrix(n int) [][]float64 {
	m := make([][]float64, n)
	for i := range m {
		m[i] = make([]float64, n)
		m[i][i] = 1
	}
	return m
}

func sign(x float64) float64 {
	if x < 0 {
		return -1
	}
	return 1
}

// RunCMAES minimizes f over R^n using CMA-ES, clipping candidate samples to the
// box for evaluation. It is deterministic given rng and returns the best point
// found.
func RunCMAES(f ObjectiveFunc, cfg CMAESConfig, rng *RNG) (Result, error) {
	if err := cfg.validate(); err != nil {
		return Result{}, err
	}
	n := cfg.Bounds.Dim()
	xmean := cfg.Mean
	if xmean == nil {
		xmean = cfg.Bounds.Center()
	} else {
		xmean = VecCopy(xmean)
	}
	if len(xmean) != n {
		return Result{}, ErrDimMismatch
	}
	sigma := cfg.Sigma
	if sigma <= 0 {
		w := cfg.Bounds.Width()
		sigma = 0.3 * w[0]
		if sigma <= 0 {
			sigma = 0.3
		}
	}

	lambda := cfg.Lambda
	if lambda <= 0 {
		lambda = 4 + int(math.Floor(3*math.Log(float64(n))))
	}
	if lambda < 4 {
		lambda = 4
	}
	mu := lambda / 2

	weights := make([]float64, mu)
	for i := 0; i < mu; i++ {
		weights[i] = math.Log(float64(mu)+0.5) - math.Log(float64(i+1))
	}
	sumW := 0.0
	for _, w := range weights {
		sumW += w
	}
	for i := range weights {
		weights[i] /= sumW
	}
	sumW2 := 0.0
	for _, w := range weights {
		sumW2 += w * w
	}
	mueff := 1 / sumW2

	nf := float64(n)
	cc := (4 + mueff/nf) / (nf + 4 + 2*mueff/nf)
	cs := (mueff + 2) / (nf + mueff + 5)
	c1 := 2 / ((nf+1.3)*(nf+1.3) + mueff)
	cmu := math.Min(1-c1, 2*(mueff-2+1/mueff)/((nf+2)*(nf+2)+mueff))
	damps := 1 + 2*math.Max(0, math.Sqrt((mueff-1)/(nf+1))-1) + cs
	chiN := math.Sqrt(nf) * (1 - 1/(4*nf) + 1/(21*nf*nf))

	pc := make([]float64, n)
	ps := make([]float64, n)
	C := identityMatrix(n)
	B := identityMatrix(n)
	D := make([]float64, n)
	for i := range D {
		D[i] = 1
	}
	invSqrtC := identityMatrix(n)

	best := VecCopy(xmean)
	bestF := f(cfg.Bounds.Clip(best))
	evals := 1
	res := Result{}

	eigenEvery := int(float64(lambda) / (c1 + cmu) / nf / 10)
	if eigenEvery < 1 {
		eigenEvery = 1
	}
	lastEigen := 0

	type sample struct {
		x []float64
		z []float64
		f float64
	}

	iter := 0
	for ; iter < cfg.MaxIterations; iter++ {
		samples := make([]sample, lambda)
		for k := 0; k < lambda; k++ {
			z := rng.GaussianVec(n, 0, 1)
			// y = B * (D .* z)
			y := make([]float64, n)
			for i := 0; i < n; i++ {
				s := 0.0
				for j := 0; j < n; j++ {
					s += B[i][j] * (D[j] * z[j])
				}
				y[i] = s
			}
			x := make([]float64, n)
			for i := 0; i < n; i++ {
				x[i] = xmean[i] + sigma*y[i]
			}
			fx := f(cfg.Bounds.Clip(x))
			evals++
			samples[k] = sample{x: x, z: y, f: fx}
		}
		sort.Slice(samples, func(a, b int) bool { return samples[a].f < samples[b].f })
		if samples[0].f < bestF {
			bestF = samples[0].f
			best = cfg.Bounds.Clip(samples[0].x)
		}

		xold := VecCopy(xmean)
		for i := 0; i < n; i++ {
			s := 0.0
			for k := 0; k < mu; k++ {
				s += weights[k] * samples[k].x[i]
			}
			xmean[i] = s
		}

		// ps update
		diff := make([]float64, n)
		for i := 0; i < n; i++ {
			diff[i] = (xmean[i] - xold[i]) / sigma
		}
		invCdiff := matVec(invSqrtC, diff)
		csFac := math.Sqrt(cs * (2 - cs) * mueff)
		for i := 0; i < n; i++ {
			ps[i] = (1-cs)*ps[i] + csFac*invCdiff[i]
		}
		psNorm := VecNorm(ps)
		hsig := 0.0
		denom := math.Sqrt(1 - math.Pow(1-cs, 2*float64(iter+1)))
		if psNorm/denom/chiN < 1.4+2/(nf+1) {
			hsig = 1
		}

		ccFac := math.Sqrt(cc * (2 - cc) * mueff)
		for i := 0; i < n; i++ {
			pc[i] = (1-cc)*pc[i] + hsig*ccFac*diff[i]
		}

		// Covariance update.
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				rankOne := pc[i] * pc[j]
				rankMu := 0.0
				for k := 0; k < mu; k++ {
					yi := (samples[k].x[i] - xold[i]) / sigma
					yj := (samples[k].x[j] - xold[j]) / sigma
					rankMu += weights[k] * yi * yj
				}
				C[i][j] = (1-c1-cmu)*C[i][j] +
					c1*(rankOne+(1-hsig)*cc*(2-cc)*C[i][j]) +
					cmu*rankMu
			}
		}

		// Step-size update.
		sigma *= math.Exp((cs / damps) * (psNorm/chiN - 1))

		// Eigen-decomposition (occasionally).
		if iter-lastEigen >= eigenEvery {
			lastEigen = iter
			// enforce symmetry
			for i := 0; i < n; i++ {
				for j := i + 1; j < n; j++ {
					C[j][i] = C[i][j]
				}
			}
			vals, vecs := JacobiEigen(C)
			B = vecs
			for i := 0; i < n; i++ {
				if vals[i] < 1e-30 {
					vals[i] = 1e-30
				}
				D[i] = math.Sqrt(vals[i])
			}
			invSqrtC = computeInvSqrt(B, vals)
		}

		if cfg.RecordHistory {
			res.History = append(res.History, bestF)
		}
		if cfg.Tolerance > 0 {
			maxD := 0.0
			for _, d := range D {
				if d > maxD {
					maxD = d
				}
			}
			if sigma*maxD < cfg.Tolerance {
				iter++
				res.Converged = true
				break
			}
		}
	}
	res.X = best
	res.F = bestF
	res.Iterations = iter
	res.Evaluations = evals
	return res, nil
}

// matVec returns the matrix-vector product m*x.
func matVec(m [][]float64, x []float64) []float64 {
	n := len(m)
	out := make([]float64, n)
	for i := 0; i < n; i++ {
		s := 0.0
		for j := 0; j < len(x); j++ {
			s += m[i][j] * x[j]
		}
		out[i] = s
	}
	return out
}

// computeInvSqrt returns C^{-1/2} = B diag(1/sqrt(vals)) B^T given eigenvectors
// B (as columns) and eigenvalues vals.
func computeInvSqrt(B [][]float64, vals []float64) [][]float64 {
	n := len(B)
	out := make([][]float64, n)
	for i := range out {
		out[i] = make([]float64, n)
	}
	inv := make([]float64, n)
	for k := 0; k < n; k++ {
		v := vals[k]
		if v < 1e-30 {
			v = 1e-30
		}
		inv[k] = 1 / math.Sqrt(v)
	}
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			s := 0.0
			for k := 0; k < n; k++ {
				s += B[i][k] * inv[k] * B[j][k]
			}
			out[i][j] = s
		}
	}
	return out
}
