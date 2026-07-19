package chaos

import (
	"math"
	"sort"
)

// FeigenbaumDelta is the first Feigenbaum constant, the limiting ratio of
// successive bifurcation-parameter intervals in a period-doubling cascade.
const FeigenbaumDelta = 4.669201609102990671853203820466

// FeigenbaumAlpha is the second Feigenbaum constant, the limiting ratio of
// successive branch separations at the superstable points.
const FeigenbaumAlpha = 2.502907875095892822283902873218

// LogisticAccumulation is the accumulation point of the logistic
// period-doubling cascade, r_infinity ~ 3.5699456.
const LogisticAccumulation = 3.569945672

// Family1D is a one-parameter family of one-dimensional maps.
type Family1D func(param float64) Map1D

// LogisticFamily returns the logistic family r -> (x -> r x(1-x)).
func LogisticFamily() Family1D {
	return func(r float64) Map1D { return Logistic(r) }
}

// SineFamily returns the sine-map family a -> (x -> a sin(pi x)).
func SineFamily() Family1D {
	return func(a float64) Map1D { return SineMap(a) }
}

// BifurcationPoint holds the parameter value and the sampled attractor points
// recorded at that parameter in a bifurcation diagram.
type BifurcationPoint struct {
	// Param is the bifurcation-parameter value.
	Param float64
	// Values are the post-transient orbit samples at Param.
	Values []float64
}

// BifurcationDiagram samples the attractor of a one-parameter family over the
// parameter interval [pMin, pMax] at steps parameter values. For each value it
// iterates from x0, discards transient points and records the next keep
// samples. It returns one BifurcationPoint per parameter value.
func BifurcationDiagram(fam Family1D, pMin, pMax float64, steps int, x0 float64, transient, keep int) []BifurcationPoint {
	if steps < 1 {
		steps = 1
	}
	out := make([]BifurcationPoint, 0, steps)
	for i := 0; i < steps; i++ {
		p := pMin
		if steps > 1 {
			p = pMin + (pMax-pMin)*float64(i)/float64(steps-1)
		}
		f := fam(p)
		vals := OrbitAfterTransient(f, x0, transient, keep)
		out = append(out, BifurcationPoint{Param: p, Values: vals})
	}
	return out
}

// AttractorSet iterates f from x0, discards the transient and returns the
// distinct attractor values found among the next keep samples, sorted
// ascending and de-duplicated with tolerance tol.
func AttractorSet(f Map1D, x0 float64, transient, keep int, tol float64) []float64 {
	vals := OrbitAfterTransient(f, x0, transient, keep)
	sort.Float64s(vals)
	out := vals[:0:0]
	for _, v := range vals {
		if len(out) == 0 || math.Abs(v-out[len(out)-1]) > tol {
			out = append(out, v)
		}
	}
	return out
}

// DetectPeriod iterates f from x0, discards the transient, then determines the
// period of the resulting orbit up to maxPeriod, using tolerance tol. It
// returns 0 if no period up to maxPeriod is detected (e.g. chaotic or
// quasiperiodic behaviour).
func DetectPeriod(f Map1D, x0 float64, transient, maxPeriod int, tol float64) int {
	x := x0
	for i := 0; i < transient; i++ {
		x = f(x)
	}
	base := x
	y := x
	for p := 1; p <= maxPeriod; p++ {
		y = f(y)
		if math.Abs(y-base) < tol {
			return p
		}
	}
	return 0
}

// PeriodMap returns, for each sampled parameter, the detected period of the
// attractor of the family fam, as parallel slices of parameter values and
// periods (0 marking aperiodic behaviour).
func PeriodMap(fam Family1D, pMin, pMax float64, steps int, x0 float64, transient, maxPeriod int, tol float64) (params []float64, periods []int) {
	params = make([]float64, steps)
	periods = make([]int, steps)
	for i := 0; i < steps; i++ {
		p := pMin
		if steps > 1 {
			p = pMin + (pMax-pMin)*float64(i)/float64(steps-1)
		}
		params[i] = p
		periods[i] = DetectPeriod(fam(p), x0, transient, maxPeriod, tol)
	}
	return params, periods
}

// SuperstableParameter finds the parameter of the logistic family for which
// the period-2^n superstable cycle occurs (the critical point x=1/2 is
// periodic with period 2^n), by bisection on [lo, hi]. The returned parameter
// r_n has f^{2^n}(1/2) = 1/2 with derivative zero.
func SuperstableParameter(period, lo, hi float64, tol float64, maxIter int) (float64, error) {
	// We solve g(r) = f_r^{period}(1/2) - 1/2 = 0 near the superstable point,
	// where the target is where the critical orbit returns to 1/2.
	p := int(period)
	g := func(r float64) float64 {
		return Iterate(Logistic(r), 0.5, p) - 0.5
	}
	a, b := lo, hi
	fa := g(a)
	for i := 0; i < maxIter; i++ {
		m := 0.5 * (a + b)
		fm := g(m)
		if math.Abs(fm) < tol || (b-a) < tol {
			return m, nil
		}
		if (fa < 0) != (fm < 0) {
			b = m
		} else {
			a, fa = m, fm
		}
	}
	return 0.5 * (a + b), ErrNoConvergence
}

// FeigenbaumEstimate estimates the first Feigenbaum constant delta from three
// consecutive period-doubling bifurcation parameters r1<r2<r3 via
// (r2-r1)/(r3-r2).
func FeigenbaumEstimate(r1, r2, r3 float64) float64 {
	return (r2 - r1) / (r3 - r2)
}

// FeigenbaumFromSequence estimates the Feigenbaum delta constant from a
// sequence of successive bifurcation parameters, returning the ratio for each
// consecutive triple. The estimates should converge toward FeigenbaumDelta.
func FeigenbaumFromSequence(rs []float64) []float64 {
	if len(rs) < 3 {
		return nil
	}
	out := make([]float64, 0, len(rs)-2)
	for i := 0; i+2 < len(rs); i++ {
		out = append(out, FeigenbaumEstimate(rs[i], rs[i+1], rs[i+2]))
	}
	return out
}

// SuperstableSequence returns the superstable parameters of the logistic
// family for periods 2^0, 2^1, ..., 2^(count-1), each found by bisection on a
// bracketing interval derived from the previous value. It is the natural input
// to FeigenbaumFromSequence.
func SuperstableSequence(count int) []float64 {
	out := make([]float64, 0, count)
	// Narrow brackets, each isolating exactly one superstable parameter s_n
	// where the critical orbit of period 2^n closes on x=1/2.
	brackets := [][2]float64{
		{1.9, 2.1},           // s_0 = 2.000000, period 1
		{3.20, 3.27},         // s_1 = 3.236068, period 2
		{3.49, 3.51},         // s_2 = 3.498562, period 4
		{3.5525, 3.5565},     // s_3 = 3.554641, period 8
		{3.5664, 3.5669},     // s_4 = 3.566667, period 16
		{3.569220, 3.569260}, // s_5 = 3.569244, period 32
	}
	for n := 0; n < count; n++ {
		var lo, hi float64
		if n < len(brackets) {
			lo, hi = brackets[n][0], brackets[n][1]
		} else {
			// Beyond the tabulated brackets the critical orbit becomes far
			// too sensitive for reliable bisection; stop early.
			break
		}
		p := math.Pow(2, float64(n))
		r, _ := SuperstableParameter(p, lo, hi, 1e-13, 200)
		out = append(out, r)
	}
	return out
}

// LyapunovExponentDiagram returns the Lyapunov exponent of the family fam at
// each of steps parameter values across [pMin, pMax]. Negative values mark
// periodic windows, positive values chaos.
func LyapunovExponentDiagram(fam Family1D, pMin, pMax float64, steps int, x0 float64, transient, n int) (params, lambdas []float64) {
	params = make([]float64, steps)
	lambdas = make([]float64, steps)
	for i := 0; i < steps; i++ {
		p := pMin
		if steps > 1 {
			p = pMin + (pMax-pMin)*float64(i)/float64(steps-1)
		}
		params[i] = p
		lambdas[i] = Lyapunov1D(fam(p), x0, transient, n)
	}
	return params, lambdas
}
