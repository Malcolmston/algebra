package stochastic

import (
	"errors"
	"math"
	"math/rand"
)

// BrownianIncrements returns n independent Brownian increments over time steps
// of size dt, each distributed as N(0, dt).
func BrownianIncrements(rng *rand.Rand, dt float64, n int) []float64 {
	if n < 0 {
		n = 0
	}
	sd := math.Sqrt(dt)
	out := make([]float64, n)
	for i := range out {
		out[i] = sd * rng.NormFloat64()
	}
	return out
}

// BrownianMotion returns a standard Brownian motion (Wiener process) on [0, T]
// sampled on nSteps+1 equally spaced grid points, starting at 0.
func BrownianMotion(rng *rand.Rand, T float64, nSteps int) Path {
	return BrownianMotionWithDrift(rng, 0, 1, T, nSteps)
}

// BrownianMotionScaled returns a Brownian motion with volatility sigma (that
// is, variance sigma^2*t) on [0, T], starting at 0.
func BrownianMotionScaled(rng *rand.Rand, sigma, T float64, nSteps int) Path {
	return BrownianMotionWithDrift(rng, 0, sigma, T, nSteps)
}

// BrownianMotionWithDrift returns a Brownian motion with drift mu and
// volatility sigma, X(t) = mu*t + sigma*W(t), on [0, T] sampled on nSteps+1
// grid points, starting at 0.
func BrownianMotionWithDrift(rng *rand.Rand, mu, sigma, T float64, nSteps int) Path {
	if nSteps < 1 {
		nSteps = 1
	}
	dt := T / float64(nSteps)
	sd := math.Sqrt(dt)
	times := make([]float64, nSteps+1)
	vals := make([]float64, nSteps+1)
	x := 0.0
	for i := 1; i <= nSteps; i++ {
		x += mu*dt + sigma*sd*rng.NormFloat64()
		times[i] = float64(i) * dt
		vals[i] = x
	}
	return Path{Times: times, Values: vals}
}

// BrownianBridge returns a Brownian bridge on [0, T] pinned to a at time 0 and
// b at time T, sampled on nSteps+1 grid points. It is constructed from a
// standard Brownian motion W by X(t) = a + (b-a)*t/T + (W(t) - (t/T)*W(T)).
func BrownianBridge(rng *rand.Rand, a, b, T float64, nSteps int) Path {
	if nSteps < 1 {
		nSteps = 1
	}
	w := BrownianMotion(rng, T, nSteps)
	wT := w.Final()
	times := make([]float64, nSteps+1)
	vals := make([]float64, nSteps+1)
	for i := 0; i <= nSteps; i++ {
		t := w.Times[i]
		times[i] = t
		vals[i] = a + (b-a)*t/T + (w.Values[i] - (t/T)*wT)
	}
	return Path{Times: times, Values: vals}
}

// StandardBrownianBridge returns a standard Brownian bridge pinned to 0 at both
// endpoints of [0, 1].
func StandardBrownianBridge(rng *rand.Rand, nSteps int) Path {
	return BrownianBridge(rng, 0, 0, 1, nSteps)
}

// GeometricBrownianMotion returns an exact simulation of geometric Brownian
// motion dS = mu*S dt + sigma*S dW on [0, T] starting at S0, sampled on
// nSteps+1 grid points. Each step uses the exact log-normal transition
// S(t+dt) = S(t)*exp((mu - sigma^2/2)dt + sigma*sqrt(dt)*Z).
func GeometricBrownianMotion(rng *rand.Rand, s0, mu, sigma, T float64, nSteps int) Path {
	if nSteps < 1 {
		nSteps = 1
	}
	dt := T / float64(nSteps)
	sd := sigma * math.Sqrt(dt)
	drift := (mu - 0.5*sigma*sigma) * dt
	times := make([]float64, nSteps+1)
	vals := make([]float64, nSteps+1)
	s := s0
	vals[0] = s0
	for i := 1; i <= nSteps; i++ {
		s *= math.Exp(drift + sd*rng.NormFloat64())
		times[i] = float64(i) * dt
		vals[i] = s
	}
	return Path{Times: times, Values: vals}
}

// GBMMean returns the analytic mean S0*exp(mu*t) of geometric Brownian motion
// at time t.
func GBMMean(s0, mu, t float64) float64 {
	return s0 * math.Exp(mu*t)
}

// GBMVariance returns the analytic variance of geometric Brownian motion at
// time t: S0^2*exp(2 mu t)*(exp(sigma^2 t) - 1).
func GBMVariance(s0, mu, sigma, t float64) float64 {
	return s0 * s0 * math.Exp(2*mu*t) * (math.Exp(sigma*sigma*t) - 1)
}

// OrnsteinUhlenbeck returns an exact simulation of the Ornstein-Uhlenbeck
// process dX = theta*(mu - X) dt + sigma dW on [0, T] starting at x0, sampled
// on nSteps+1 grid points. Each step uses the exact Gaussian transition with
// mean mu + (X - mu)*exp(-theta*dt) and variance sigma^2*(1 - exp(-2 theta
// dt))/(2 theta).
func OrnsteinUhlenbeck(rng *rand.Rand, x0, theta, mu, sigma, T float64, nSteps int) Path {
	if nSteps < 1 {
		nSteps = 1
	}
	dt := T / float64(nSteps)
	ed := math.Exp(-theta * dt)
	var sd float64
	if theta > 0 {
		sd = sigma * math.Sqrt((1-ed*ed)/(2*theta))
	} else {
		sd = sigma * math.Sqrt(dt)
	}
	times := make([]float64, nSteps+1)
	vals := make([]float64, nSteps+1)
	x := x0
	vals[0] = x0
	for i := 1; i <= nSteps; i++ {
		x = mu + (x-mu)*ed + sd*rng.NormFloat64()
		times[i] = float64(i) * dt
		vals[i] = x
	}
	return Path{Times: times, Values: vals}
}

// VasicekProcess is an alias for OrnsteinUhlenbeck using interest-rate
// terminology: dr = a*(b - r) dt + sigma dW.
func VasicekProcess(rng *rand.Rand, r0, a, b, sigma, T float64, nSteps int) Path {
	return OrnsteinUhlenbeck(rng, r0, a, b, sigma, T, nSteps)
}

// OUMean returns the analytic mean mu + (x0 - mu)*exp(-theta*t) of the
// Ornstein-Uhlenbeck process at time t.
func OUMean(x0, theta, mu, t float64) float64 {
	return mu + (x0-mu)*math.Exp(-theta*t)
}

// OUVariance returns the analytic variance sigma^2*(1 - exp(-2 theta t))/(2
// theta) of the Ornstein-Uhlenbeck process at time t.
func OUVariance(theta, sigma, t float64) float64 {
	if theta <= 0 {
		return sigma * sigma * t
	}
	return sigma * sigma * (1 - math.Exp(-2*theta*t)) / (2 * theta)
}

// OUStationaryVariance returns the stationary variance sigma^2/(2 theta) of the
// Ornstein-Uhlenbeck process.
func OUStationaryVariance(theta, sigma float64) float64 {
	if theta <= 0 {
		return math.Inf(1)
	}
	return sigma * sigma / (2 * theta)
}

// CIRProcess simulates the Cox-Ingersoll-Ross process dX = a*(b - X) dt +
// sigma*sqrt(X) dW on [0, T] starting at x0 using a full-truncation Euler
// scheme that keeps the state non-negative. The result is sampled on nSteps+1
// grid points.
func CIRProcess(rng *rand.Rand, x0, a, b, sigma, T float64, nSteps int) Path {
	if nSteps < 1 {
		nSteps = 1
	}
	dt := T / float64(nSteps)
	sd := math.Sqrt(dt)
	times := make([]float64, nSteps+1)
	vals := make([]float64, nSteps+1)
	x := x0
	vals[0] = x0
	for i := 1; i <= nSteps; i++ {
		xp := x
		if xp < 0 {
			xp = 0
		}
		x = x + a*(b-xp)*dt + sigma*math.Sqrt(xp)*sd*rng.NormFloat64()
		if x < 0 {
			x = 0
		}
		times[i] = float64(i) * dt
		vals[i] = x
	}
	return Path{Times: times, Values: vals}
}

// BrownianCovariance returns the covariance min(s, t) of standard Brownian
// motion at times s and t.
func BrownianCovariance(s, t float64) float64 {
	return math.Min(s, t)
}

// BrownianBridgeCovariance returns the covariance s*(T-t)/T (for s <= t) of a
// standard Brownian bridge on [0, T].
func BrownianBridgeCovariance(s, t, T float64) float64 {
	if s > t {
		s, t = t, s
	}
	return s * (T - t) / T
}

// MultiBrownianMotion returns dim independent standard Brownian motions on
// [0, T] sampled on the same grid of nSteps+1 points.
func MultiBrownianMotion(rng *rand.Rand, dim int, T float64, nSteps int) []Path {
	out := make([]Path, dim)
	for d := 0; d < dim; d++ {
		out[d] = BrownianMotion(rng, T, nSteps)
	}
	return out
}

// CorrelatedBrownianMotion returns len(corr) Brownian motions on [0, T] with
// the given instantaneous correlation matrix, sampled on nSteps+1 grid points.
// The correlation matrix must be symmetric positive definite; otherwise an
// error is returned.
func CorrelatedBrownianMotion(rng *rand.Rand, corr [][]float64, T float64, nSteps int) ([]Path, error) {
	d := len(corr)
	l, err := cholesky(corr)
	if err != nil {
		return nil, err
	}
	if nSteps < 1 {
		nSteps = 1
	}
	dt := T / float64(nSteps)
	sd := math.Sqrt(dt)
	paths := make([]Path, d)
	for k := 0; k < d; k++ {
		paths[k].Times = make([]float64, nSteps+1)
		paths[k].Values = make([]float64, nSteps+1)
	}
	x := make([]float64, d)
	z := make([]float64, d)
	for i := 1; i <= nSteps; i++ {
		for k := 0; k < d; k++ {
			z[k] = rng.NormFloat64()
		}
		for k := 0; k < d; k++ {
			corrNoise := 0.0
			for j := 0; j <= k; j++ {
				corrNoise += l[k][j] * z[j]
			}
			x[k] += sd * corrNoise
			paths[k].Times[i] = float64(i) * dt
			paths[k].Values[i] = x[k]
		}
	}
	return paths, nil
}

// FractionalBrownianMotion returns a fractional Brownian motion with Hurst
// exponent H (0 < H < 1) on [0, T] sampled on nSteps+1 grid points, using the
// Cholesky factorization of its covariance matrix. H = 1/2 recovers standard
// Brownian motion.
func FractionalBrownianMotion(rng *rand.Rand, hurst, T float64, nSteps int) (Path, error) {
	if hurst <= 0 || hurst >= 1 {
		return Path{}, errors.New("stochastic: Hurst exponent must lie in (0,1)")
	}
	if nSteps < 1 {
		nSteps = 1
	}
	dt := T / float64(nSteps)
	// covariance of fBm at times t_i = i*dt, i = 1..nSteps
	cov := make([][]float64, nSteps)
	for i := 0; i < nSteps; i++ {
		cov[i] = make([]float64, nSteps)
		s := float64(i+1) * dt
		for j := 0; j < nSteps; j++ {
			u := float64(j+1) * dt
			cov[i][j] = 0.5 * (math.Pow(s, 2*hurst) + math.Pow(u, 2*hurst) - math.Pow(math.Abs(s-u), 2*hurst))
		}
	}
	l, err := cholesky(cov)
	if err != nil {
		return Path{}, err
	}
	z := NormalVector(rng, nSteps)
	times := make([]float64, nSteps+1)
	vals := make([]float64, nSteps+1)
	for i := 0; i < nSteps; i++ {
		s := 0.0
		for j := 0; j <= i; j++ {
			s += l[i][j] * z[j]
		}
		times[i+1] = float64(i+1) * dt
		vals[i+1] = s
	}
	return Path{Times: times, Values: vals}, nil
}

// cholesky returns the lower-triangular Cholesky factor L of a symmetric
// positive-definite matrix a, so that L*L^T = a.
func cholesky(a [][]float64) ([][]float64, error) {
	n := len(a)
	l := make([][]float64, n)
	for i := range l {
		if len(a[i]) != n {
			return nil, errors.New("stochastic: matrix must be square")
		}
		l[i] = make([]float64, n)
	}
	for i := 0; i < n; i++ {
		for j := 0; j <= i; j++ {
			sum := a[i][j]
			for k := 0; k < j; k++ {
				sum -= l[i][k] * l[j][k]
			}
			if i == j {
				if sum <= 0 {
					return nil, errors.New("stochastic: matrix is not positive definite")
				}
				l[i][j] = math.Sqrt(sum)
			} else {
				l[i][j] = sum / l[j][j]
			}
		}
	}
	return l, nil
}
