package stochastic

import (
	"math"
	"math/rand"
)

// SimpleRandomWalk returns the integer positions of an n-step simple random
// walk that starts at 0 and takes a +1 step with probability p and a -1 step
// otherwise. The result has length n+1 (including the start).
func SimpleRandomWalk(rng *rand.Rand, n int, p float64) []int {
	if n < 0 {
		n = 0
	}
	pos := make([]int, n+1)
	x := 0
	pos[0] = 0
	for i := 1; i <= n; i++ {
		if rng.Float64() < p {
			x++
		} else {
			x--
		}
		pos[i] = x
	}
	return pos
}

// SimpleRandomWalkPath returns a simple random walk as a Path with unit time
// steps.
func SimpleRandomWalkPath(rng *rand.Rand, n int, p float64) Path {
	pos := SimpleRandomWalk(rng, n, p)
	times := make([]float64, len(pos))
	vals := make([]float64, len(pos))
	for i, v := range pos {
		times[i] = float64(i)
		vals[i] = float64(v)
	}
	return Path{Times: times, Values: vals}
}

// SymmetricRandomWalk returns a simple symmetric random walk (p = 1/2).
func SymmetricRandomWalk(rng *rand.Rand, n int) []int {
	return SimpleRandomWalk(rng, n, 0.5)
}

// LazyRandomWalk returns an n-step lazy random walk that stays put with
// probability stay, steps +1 with probability (1-stay)*p and -1 otherwise.
func LazyRandomWalk(rng *rand.Rand, n int, p, stay float64) []int {
	if n < 0 {
		n = 0
	}
	pos := make([]int, n+1)
	x := 0
	for i := 1; i <= n; i++ {
		u := rng.Float64()
		if u < stay {
			// stay
		} else if u < stay+(1-stay)*p {
			x++
		} else {
			x--
		}
		pos[i] = x
	}
	return pos
}

// GaussianRandomWalk returns an n-step random walk with independent normal
// increments of mean mu and standard deviation sigma, as a Path with unit time
// steps starting at x0.
func GaussianRandomWalk(rng *rand.Rand, n int, x0, mu, sigma float64) Path {
	if n < 0 {
		n = 0
	}
	times := make([]float64, n+1)
	vals := make([]float64, n+1)
	x := x0
	vals[0] = x
	for i := 1; i <= n; i++ {
		x += mu + sigma*rng.NormFloat64()
		times[i] = float64(i)
		vals[i] = x
	}
	return Path{Times: times, Values: vals}
}

// RandomWalk2D returns the x and y coordinates of an n-step random walk on the
// integer plane, where each step moves one unit in one of the four axis
// directions chosen uniformly. Each returned slice has length n+1.
func RandomWalk2D(rng *rand.Rand, n int) (xs, ys []int) {
	if n < 0 {
		n = 0
	}
	xs = make([]int, n+1)
	ys = make([]int, n+1)
	x, y := 0, 0
	for i := 1; i <= n; i++ {
		switch rng.Intn(4) {
		case 0:
			x++
		case 1:
			x--
		case 2:
			y++
		default:
			y--
		}
		xs[i] = x
		ys[i] = y
	}
	return xs, ys
}

// ReflectingRandomWalk returns an n-step simple random walk confined to the
// closed interval [lo, hi] by reflecting boundaries: a step that would leave the
// interval is suppressed and the walk stays in place.
func ReflectingRandomWalk(rng *rand.Rand, n int, p float64, lo, hi int) []int {
	if n < 0 {
		n = 0
	}
	pos := make([]int, n+1)
	x := 0
	if x < lo {
		x = lo
	}
	if x > hi {
		x = hi
	}
	pos[0] = x
	for i := 1; i <= n; i++ {
		step := -1
		if rng.Float64() < p {
			step = 1
		}
		nx := x + step
		if nx >= lo && nx <= hi {
			x = nx
		}
		pos[i] = x
	}
	return pos
}

// AbsorbingRandomWalk simulates a simple random walk that starts at start,
// steps +1 with probability p, and stops when it first reaches either boundary
// lo or hi. It returns the absorbing position reached and the number of steps
// taken.
func AbsorbingRandomWalk(rng *rand.Rand, start, lo, hi int, p float64) (position, steps int) {
	x := start
	for x > lo && x < hi {
		if rng.Float64() < p {
			x++
		} else {
			x--
		}
		steps++
	}
	return x, steps
}

// ContinuousTimeRandomWalk returns the sample path of a continuous-time random
// walk that waits an Exponential(rate) time between jumps and adds a jump drawn
// from jump(rng), running until time T. The path is piecewise constant and its
// grid consists of the jump times plus the endpoints 0 and T.
func ContinuousTimeRandomWalk(rng *rand.Rand, rate, T float64, jump func(*rand.Rand) float64) Path {
	times := []float64{0}
	vals := []float64{0}
	t := 0.0
	x := 0.0
	for {
		t += ExponentialSample(rng, rate)
		if t > T {
			break
		}
		x += jump(rng)
		times = append(times, t)
		vals = append(vals, x)
	}
	times = append(times, T)
	vals = append(vals, x)
	return Path{Times: times, Values: vals}
}

// LevyFlight returns an n-step random walk in one dimension whose step sizes are
// symmetric with heavy Pareto-distributed magnitudes of tail index alpha
// (0 < alpha < 2), as a Path with unit time steps.
func LevyFlight(rng *rand.Rand, n int, alpha float64) Path {
	if n < 0 {
		n = 0
	}
	times := make([]float64, n+1)
	vals := make([]float64, n+1)
	x := 0.0
	for i := 1; i <= n; i++ {
		mag := ParetoSample(rng, 1, alpha)
		x += float64(RademacherSample(rng)) * mag
		times[i] = float64(i)
		vals[i] = x
	}
	return Path{Times: times, Values: vals}
}

// RandomWalkMean returns the analytic mean n*(2p-1) of the endpoint of an
// n-step simple random walk with up-probability p.
func RandomWalkMean(n int, p float64) float64 {
	return float64(n) * (2*p - 1)
}

// RandomWalkVariance returns the analytic variance 4*n*p*(1-p) of the endpoint
// of an n-step simple random walk with up-probability p.
func RandomWalkVariance(n int, p float64) float64 {
	return 4 * float64(n) * p * (1 - p)
}

// GamblersRuinProbability returns the probability that a simple random walk
// starting at i, absorbing at 0 and N, and stepping up with probability p,
// reaches 0 before N (the gambler's ruin probability).
func GamblersRuinProbability(i, N int, p float64) float64 {
	if i <= 0 {
		return 1
	}
	if i >= N {
		return 0
	}
	if math.Abs(p-0.5) < 1e-15 {
		return float64(N-i) / float64(N)
	}
	q := 1 - p
	r := q / p
	// probability of winning (reaching N) is (1 - r^i)/(1 - r^N)
	win := (1 - math.Pow(r, float64(i))) / (1 - math.Pow(r, float64(N)))
	return 1 - win
}

// GamblersRuinWinProbability returns the probability that the walk of
// GamblersRuinProbability reaches N before 0.
func GamblersRuinWinProbability(i, N int, p float64) float64 {
	return 1 - GamblersRuinProbability(i, N, p)
}

// GamblersRuinExpectedDuration returns the expected number of steps until a
// simple random walk starting at i, absorbing at 0 and N and stepping up with
// probability p, is absorbed.
func GamblersRuinExpectedDuration(i, N int, p float64) float64 {
	if i <= 0 || i >= N {
		return 0
	}
	if math.Abs(p-0.5) < 1e-15 {
		return float64(i) * float64(N-i)
	}
	q := 1 - p
	r := q / p
	fi := float64(i)
	fN := float64(N)
	return fi/(q-p) - fN/(q-p)*(1-math.Pow(r, fi))/(1-math.Pow(r, fN))
}
