package infogeom

import "math"

// Sigmoid returns the logistic sigmoid 1/(1+e^{-x}), computed in a way that
// avoids overflow for large-magnitude x.
func Sigmoid(x float64) float64 {
	if x >= 0 {
		return 1 / (1 + math.Exp(-x))
	}
	e := math.Exp(x)
	return e / (1 + e)
}

// LogSigmoid returns ln sigmoid(x) = -ln(1+e^{-x}), computed stably.
func LogSigmoid(x float64) float64 {
	return -logOnePlusExp(-x)
}

// Logit returns the log-odds ln(p/(1-p)), the inverse of Sigmoid. It returns
// -Inf at p=0 and +Inf at p=1.
func Logit(p float64) float64 {
	return math.Log(p / (1 - p))
}

// logOnePlusExp returns ln(1+e^x) computed without overflow.
func logOnePlusExp(x float64) float64 {
	if x > 0 {
		return x + math.Log1p(math.Exp(-x))
	}
	return math.Log1p(math.Exp(x))
}

// LogSumExp returns ln sum_i e^{x_i}, computed by factoring out the maximum to
// avoid overflow. For an empty slice it returns -Inf.
func LogSumExp(x []float64) float64 {
	if len(x) == 0 {
		return math.Inf(-1)
	}
	m := x[0]
	for _, xi := range x[1:] {
		if xi > m {
			m = xi
		}
	}
	if math.IsInf(m, -1) {
		return m
	}
	var s float64
	for _, xi := range x {
		s += math.Exp(xi - m)
	}
	return m + math.Log(s)
}

// Softmax returns the softmax (Gibbs) distribution with logits x, a probability
// vector proportional to e^{x_i}. It is computed stably via LogSumExp.
func Softmax(x []float64) []float64 {
	out := make([]float64, len(x))
	if len(x) == 0 {
		return out
	}
	lse := LogSumExp(x)
	for i, xi := range x {
		out[i] = math.Exp(xi - lse)
	}
	return out
}

// LogSoftmax returns the elementwise logarithm of Softmax(x), x_i - LogSumExp(x).
func LogSoftmax(x []float64) []float64 {
	out := make([]float64, len(x))
	if len(x) == 0 {
		return out
	}
	lse := LogSumExp(x)
	for i, xi := range x {
		out[i] = xi - lse
	}
	return out
}

// logFactorial returns ln(k!) via the log-gamma function for k >= 0.
func logFactorial(k int) float64 {
	if k < 2 {
		return 0
	}
	lg, _ := math.Lgamma(float64(k) + 1)
	return lg
}

// ExponentialFamily describes a (curved or flat) exponential family through its
// log-partition function and, optionally, the gradient and Hessian of that
// function. In natural coordinates theta the density has the form
// h(x) exp( <theta, T(x)> - A(theta) ), the gradient of A is the expectation
// parameter (mean of the sufficient statistic T) and the Hessian of A is the
// Fisher information matrix.
type ExponentialFamily struct {
	// LogPartition evaluates the cumulant function A(theta).
	LogPartition func(theta []float64) float64
	// GradLogPartition, when non-nil, returns the expectation parameters
	// eta = grad A(theta). When nil the package falls back to numerical
	// differentiation of LogPartition.
	GradLogPartition func(theta []float64) []float64
	// HessLogPartition, when non-nil, returns the Fisher information matrix
	// in natural coordinates, the Hessian of A. When nil the package falls
	// back to numerical differentiation.
	HessLogPartition func(theta []float64) [][]float64
	// SufficientStat, when non-nil, maps an observation to its sufficient
	// statistic T(x).
	SufficientStat func(x float64) []float64
	// LogBaseMeasure, when non-nil, evaluates ln h(x).
	LogBaseMeasure func(x float64) float64
}

// MeanParameters returns the expectation parameters eta = grad A(theta) of the
// family at the natural parameter theta, using the analytic gradient when
// available and a central finite difference otherwise.
func (f ExponentialFamily) MeanParameters(theta []float64) []float64 {
	if f.GradLogPartition != nil {
		return f.GradLogPartition(theta)
	}
	return NumericalGradient(f.LogPartition, theta, 1e-6)
}

// FisherInformationNatural returns the Fisher information matrix in natural
// coordinates, the Hessian of A(theta), using the analytic Hessian when
// available and a finite-difference approximation otherwise.
func (f ExponentialFamily) FisherInformationNatural(theta []float64) [][]float64 {
	if f.HessLogPartition != nil {
		return f.HessLogPartition(theta)
	}
	return NumericalHessian(f.LogPartition, theta, 1e-4)
}

// LogDensity returns the log-density ln p(x; theta) = ln h(x) + <theta, T(x)>
// - A(theta). It returns ErrDim when SufficientStat is nil or produces a
// statistic whose length differs from theta.
func (f ExponentialFamily) LogDensity(theta []float64, x float64) (float64, error) {
	if f.SufficientStat == nil {
		return 0, ErrDim
	}
	t := f.SufficientStat(x)
	inner, err := Dot(theta, t)
	if err != nil {
		return 0, err
	}
	logh := 0.0
	if f.LogBaseMeasure != nil {
		logh = f.LogBaseMeasure(x)
	}
	return logh + inner - f.LogPartition(theta), nil
}

// Density returns the density p(x; theta) = exp(LogDensity).
func (f ExponentialFamily) Density(theta []float64, x float64) (float64, error) {
	ld, err := f.LogDensity(theta, x)
	if err != nil {
		return 0, err
	}
	return math.Exp(ld), nil
}

// LegendreDual returns the value of the convex conjugate (Legendre transform)
// A*(eta) = <theta, eta> - A(theta) evaluated at the dual pair (theta, eta)
// related by eta = grad A(theta). For an exponential family A* is the negative
// entropy expressed in expectation coordinates. It returns ErrDim on a length
// mismatch.
func LegendreDual(a func(theta []float64) float64, theta, eta []float64) (float64, error) {
	inner, err := Dot(theta, eta)
	if err != nil {
		return 0, err
	}
	return inner - a(theta), nil
}

// KLDivergenceExpFamily returns the Kullback-Leibler divergence D(theta1 ||
// theta2) between two members of the same exponential family, expressed as the
// Bregman divergence of the log-partition function,
//
//	A(theta2) - A(theta1) - <grad A(theta1), theta2 - theta1>,
//
// which equals D(p_{theta1} || p_{theta2}). It returns ErrDim on a length
// mismatch.
func KLDivergenceExpFamily(f ExponentialFamily, theta1, theta2 []float64) (float64, error) {
	if len(theta1) != len(theta2) || len(theta1) == 0 {
		return 0, ErrDim
	}
	eta1 := f.MeanParameters(theta1)
	diff, err := SubVectors(theta2, theta1)
	if err != nil {
		return 0, err
	}
	inner, err := Dot(eta1, diff)
	if err != nil {
		return 0, err
	}
	return f.LogPartition(theta2) - f.LogPartition(theta1) - inner, nil
}

// BernoulliFamily returns the ExponentialFamily representation of the Bernoulli
// distribution: sufficient statistic x, log-partition ln(1+e^theta), mean
// sigmoid(theta) and Fisher information sigmoid(theta)(1-sigmoid(theta)).
func BernoulliFamily() ExponentialFamily {
	return ExponentialFamily{
		LogPartition: func(theta []float64) float64 { return logOnePlusExp(theta[0]) },
		GradLogPartition: func(theta []float64) []float64 {
			return []float64{Sigmoid(theta[0])}
		},
		HessLogPartition: func(theta []float64) [][]float64 {
			s := Sigmoid(theta[0])
			return [][]float64{{s * (1 - s)}}
		},
		SufficientStat: func(x float64) []float64 { return []float64{x} },
		LogBaseMeasure: func(x float64) float64 { return 0 },
	}
}

// PoissonFamily returns the ExponentialFamily representation of the Poisson
// distribution with sufficient statistic x, log-partition e^theta, mean
// e^theta and base measure 1/x!.
func PoissonFamily() ExponentialFamily {
	return ExponentialFamily{
		LogPartition: func(theta []float64) float64 { return math.Exp(theta[0]) },
		GradLogPartition: func(theta []float64) []float64 {
			return []float64{math.Exp(theta[0])}
		},
		HessLogPartition: func(theta []float64) [][]float64 {
			return [][]float64{{math.Exp(theta[0])}}
		},
		SufficientStat: func(x float64) []float64 { return []float64{x} },
		LogBaseMeasure: func(x float64) float64 { return -logFactorial(int(x)) },
	}
}

// GaussianKnownVarianceFamily returns the ExponentialFamily representation of
// the normal distribution with known variance sigma^2 and unknown mean, whose
// natural parameter is theta = mu/sigma^2, log-partition sigma^2 theta^2 / 2,
// mean sigma^2 theta and Fisher information sigma^2. It returns ErrDomain when
// sigma is not positive.
func GaussianKnownVarianceFamily(sigma float64) (ExponentialFamily, error) {
	if sigma <= 0 {
		return ExponentialFamily{}, ErrDomain
	}
	s2 := sigma * sigma
	return ExponentialFamily{
		LogPartition: func(theta []float64) float64 { return 0.5 * s2 * theta[0] * theta[0] },
		GradLogPartition: func(theta []float64) []float64 {
			return []float64{s2 * theta[0]}
		},
		HessLogPartition: func(theta []float64) [][]float64 {
			return [][]float64{{s2}}
		},
		SufficientStat: func(x float64) []float64 { return []float64{x} },
		LogBaseMeasure: func(x float64) float64 {
			return -0.5*x*x/s2 - 0.5*math.Log(2*math.Pi*s2)
		},
	}, nil
}

// CategoricalFamily returns the ExponentialFamily representation of a
// categorical distribution over k outcomes with the last outcome as reference,
// using the k-1 natural parameters theta_i = ln(p_i/p_k). The log-partition is
// ln(1 + sum e^{theta_i}). It returns ErrDim when k < 2.
func CategoricalFamily(k int) (ExponentialFamily, error) {
	if k < 2 {
		return ExponentialFamily{}, ErrDim
	}
	logPart := func(theta []float64) float64 {
		ext := make([]float64, len(theta)+1)
		copy(ext, theta)
		return LogSumExp(ext)
	}
	grad := func(theta []float64) []float64 {
		ext := make([]float64, len(theta)+1)
		copy(ext, theta)
		p := Softmax(ext)
		return p[:len(theta)]
	}
	hess := func(theta []float64) [][]float64 {
		ext := make([]float64, len(theta)+1)
		copy(ext, theta)
		p := Softmax(ext)
		n := len(theta)
		h := make([][]float64, n)
		for i := 0; i < n; i++ {
			h[i] = make([]float64, n)
			for j := 0; j < n; j++ {
				if i == j {
					h[i][j] = p[i] * (1 - p[i])
				} else {
					h[i][j] = -p[i] * p[j]
				}
			}
		}
		return h
	}
	return ExponentialFamily{
		LogPartition:     logPart,
		GradLogPartition: grad,
		HessLogPartition: hess,
	}, nil
}
