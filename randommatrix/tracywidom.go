package randommatrix

import "math"

// twParams holds the shifted-gamma approximation parameters of a Tracy-Widom
// distribution: the gamma shape k, scale theta, and the literature mean used to
// fix the location shift.
type twParams struct {
	k, theta, mean float64
}

// twTable stores the Chiani (2014) shape and scale together with the literature
// mean for the Dyson indices 1, 2 and 4.
var twTable = map[int]twParams{
	1: {k: 46.44604884387853, theta: 0.18605402503251536, mean: -1.2065335745820},
	2: {k: 79.6594870748196, theta: 0.10103655775856243, mean: -1.7710868074},
	4: {k: 146.02119327131056, theta: 0.059544295001525084, mean: -2.3068848932},
}

// TracyWidomSupportedBeta reports whether beta is one of the supported Dyson
// indices 1, 2 or 4.
func TracyWidomSupportedBeta(beta int) bool {
	_, ok := twTable[beta]
	return ok
}

// TracyWidomApproxParams returns the shifted-gamma approximation parameters
// (shape k, scale theta, location shift) for the Tracy-Widom distribution of
// Dyson index beta. The approximation, due to Chiani (2014), models
// TW_beta ~ theta * Gamma(k) - shift and is accurate to about 0.01 in
// distribution. It returns ok = false for an unsupported beta.
func TracyWidomApproxParams(beta int) (k, theta, shift float64, ok bool) {
	p, ok := twTable[beta]
	if !ok {
		return 0, 0, 0, false
	}
	shift = p.k*p.theta - p.mean
	return p.k, p.theta, shift, true
}

// regGammaP returns the regularized lower incomplete gamma function P(a, x).
func regGammaP(a, x float64) float64 {
	if x < 0 || a <= 0 {
		return math.NaN()
	}
	if x == 0 {
		return 0
	}
	lg, _ := math.Lgamma(a)
	if x < a+1 {
		ap := a
		sum := 1 / a
		del := sum
		for n := 0; n < 500; n++ {
			ap++
			del *= x / ap
			sum += del
			if math.Abs(del) < math.Abs(sum)*1e-16 {
				break
			}
		}
		return sum * math.Exp(-x+a*math.Log(x)-lg)
	}
	// Lentz continued fraction for Q(a, x).
	const tiny = 1e-300
	b := x + 1 - a
	c := 1 / tiny
	d := 1 / b
	h := d
	for i := 1; i < 500; i++ {
		an := -float64(i) * (float64(i) - a)
		b += 2
		d = an*d + b
		if math.Abs(d) < tiny {
			d = tiny
		}
		c = b + an/c
		if math.Abs(c) < tiny {
			c = tiny
		}
		d = 1 / d
		del := d * c
		h *= del
		if math.Abs(del-1) < 1e-16 {
			break
		}
	}
	q := math.Exp(-x+a*math.Log(x)-lg) * h
	return 1 - q
}

// gammaPDF returns the density of the Gamma(shape=k, scale=1) distribution at x.
func gammaPDF(x, k float64) float64 {
	if x <= 0 {
		return 0
	}
	lg, _ := math.Lgamma(k)
	return math.Exp((k-1)*math.Log(x) - x - lg)
}

// TracyWidomCDF returns the Chiani shifted-gamma approximation of the
// Tracy-Widom cumulative distribution function of Dyson index beta at s. It
// returns NaN for an unsupported beta.
func TracyWidomCDF(s float64, beta int) float64 {
	k, theta, shift, ok := TracyWidomApproxParams(beta)
	if !ok {
		return math.NaN()
	}
	x := (s + shift) / theta
	if x <= 0 {
		return 0
	}
	return regGammaP(k, x)
}

// TracyWidomDensity returns the Chiani shifted-gamma approximation of the
// Tracy-Widom probability density of Dyson index beta at s.
func TracyWidomDensity(s float64, beta int) float64 {
	k, theta, shift, ok := TracyWidomApproxParams(beta)
	if !ok {
		return math.NaN()
	}
	x := (s + shift) / theta
	return gammaPDF(x, k) / theta
}

// TracyWidomMean returns the tabulated mean of the Tracy-Widom distribution of
// Dyson index beta.
func TracyWidomMean(beta int) float64 {
	p, ok := twTable[beta]
	if !ok {
		return math.NaN()
	}
	return p.mean
}

// TracyWidomVariance returns the variance k*theta^2 of the shifted-gamma
// approximation, which matches the tabulated Tracy-Widom variance.
func TracyWidomVariance(beta int) float64 {
	p, ok := twTable[beta]
	if !ok {
		return math.NaN()
	}
	return p.k * p.theta * p.theta
}

// TracyWidomStdDev returns the standard deviation of the Tracy-Widom
// distribution of Dyson index beta.
func TracyWidomStdDev(beta int) float64 {
	return math.Sqrt(TracyWidomVariance(beta))
}

// TracyWidomSkewness returns the skewness 2/sqrt(k) of the shifted-gamma
// approximation to the Tracy-Widom distribution of Dyson index beta.
func TracyWidomSkewness(beta int) float64 {
	p, ok := twTable[beta]
	if !ok {
		return math.NaN()
	}
	return 2 / math.Sqrt(p.k)
}

// TracyWidomExcessKurtosis returns the excess kurtosis 6/k of the shifted-gamma
// approximation to the Tracy-Widom distribution of Dyson index beta.
func TracyWidomExcessKurtosis(beta int) float64 {
	p, ok := twTable[beta]
	if !ok {
		return math.NaN()
	}
	return 6 / p.k
}

// TracyWidomEdgeCenter returns the soft-edge centring 2*sqrt(n) for the largest
// eigenvalue of an n-by-n Gaussian ensemble in the normalisation where the bulk
// fills [-2 sqrt(n), 2 sqrt(n)].
func TracyWidomEdgeCenter(n int) float64 {
	return 2 * math.Sqrt(float64(n))
}

// TracyWidomEdgeScale returns the soft-edge scaling n^(-1/6), so that
// (lambda_max - 2 sqrt(n)) / scale converges to the Tracy-Widom law.
func TracyWidomEdgeScale(n int) float64 {
	return math.Pow(float64(n), -1.0/6.0)
}

// RescaleLargestEigenvalue maps the largest eigenvalue lambdaMax of an n-by-n
// Gaussian ensemble to the Tracy-Widom scale via
// (lambdaMax - 2 sqrt(n)) * n^(1/6).
func RescaleLargestEigenvalue(lambdaMax float64, n int) float64 {
	return (lambdaMax - TracyWidomEdgeCenter(n)) / TracyWidomEdgeScale(n)
}

// WishartTracyWidomCenter returns the Johnstone centring
// (sqrt(n-1) + sqrt(p))^2 for the largest eigenvalue of a real Wishart matrix
// built from an n-by-p Gaussian data matrix.
func WishartTracyWidomCenter(n, p int) float64 {
	a := math.Sqrt(float64(n - 1))
	b := math.Sqrt(float64(p))
	return (a + b) * (a + b)
}

// WishartTracyWidomScale returns the Johnstone scaling
// (sqrt(n-1) + sqrt(p)) (1/sqrt(n-1) + 1/sqrt(p))^(1/3) for the largest
// eigenvalue of a real Wishart matrix, so that (lambda_max - center)/scale
// converges to the Tracy-Widom law with beta = 1.
func WishartTracyWidomScale(n, p int) float64 {
	a := math.Sqrt(float64(n - 1))
	b := math.Sqrt(float64(p))
	return (a + b) * math.Cbrt(1/a+1/b)
}
