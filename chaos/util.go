package chaos

import (
	"math"
	"math/rand"
	"sort"
)

// Linspace returns n evenly spaced values from a to b inclusive.
func Linspace(a, b float64, n int) []float64 {
	if n <= 0 {
		return nil
	}
	if n == 1 {
		return []float64{a}
	}
	out := make([]float64, n)
	for i := 0; i < n; i++ {
		out[i] = a + (b-a)*float64(i)/float64(n-1)
	}
	return out
}

// Logspace returns n values geometrically spaced between a and b (both > 0).
func Logspace(a, b float64, n int) []float64 {
	if n <= 0 || a <= 0 || b <= 0 {
		return nil
	}
	if n == 1 {
		return []float64{a}
	}
	out := make([]float64, n)
	la, lb := math.Log(a), math.Log(b)
	for i := 0; i < n; i++ {
		out[i] = math.Exp(la + (lb-la)*float64(i)/float64(n-1))
	}
	return out
}

// Mean returns the arithmetic mean of xs, or NaN when empty.
func Mean(xs []float64) float64 {
	if len(xs) == 0 {
		return math.NaN()
	}
	var s float64
	for _, x := range xs {
		s += x
	}
	return s / float64(len(xs))
}

// Variance returns the population variance of xs.
func Variance(xs []float64) float64 {
	if len(xs) == 0 {
		return math.NaN()
	}
	m := Mean(xs)
	var s float64
	for _, x := range xs {
		d := x - m
		s += d * d
	}
	return s / float64(len(xs))
}

// StdDev returns the population standard deviation of xs.
func StdDev(xs []float64) float64 {
	return math.Sqrt(Variance(xs))
}

// Autocorrelation returns the sample autocorrelation of xs at the given lag.
func Autocorrelation(xs []float64, lag int) float64 {
	n := len(xs)
	if lag < 0 || lag >= n {
		return math.NaN()
	}
	m := Mean(xs)
	var num, den float64
	for i := 0; i < n; i++ {
		d := xs[i] - m
		den += d * d
	}
	for i := 0; i+lag < n; i++ {
		num += (xs[i] - m) * (xs[i+lag] - m)
	}
	if den == 0 {
		return 0
	}
	return num / den
}

// Histogram bins the values xs into bins equal-width buckets over the range
// [lo, hi] and returns the bucket counts and their centres.
func Histogram(xs []float64, lo, hi float64, bins int) (counts []int, centers []float64) {
	if bins <= 0 || hi <= lo {
		return nil, nil
	}
	counts = make([]int, bins)
	centers = make([]float64, bins)
	w := (hi - lo) / float64(bins)
	for i := 0; i < bins; i++ {
		centers[i] = lo + (float64(i)+0.5)*w
	}
	for _, x := range xs {
		if x < lo || x > hi {
			continue
		}
		k := int((x - lo) / w)
		if k == bins {
			k = bins - 1
		}
		if k >= 0 && k < bins {
			counts[k]++
		}
	}
	return counts, centers
}

// InvariantDensity estimates the invariant density of the one-dimensional map
// f by iterating from x0, discarding the transient and forming a normalised
// histogram of the orbit over [lo, hi]. It returns bin centres and probability
// densities.
func InvariantDensity(f Map1D, x0, lo, hi float64, bins, transient, n int) (centers, density []float64) {
	orbit := OrbitAfterTransient(f, x0, transient, n)
	counts, centers := Histogram(orbit, lo, hi, bins)
	density = make([]float64, len(counts))
	w := (hi - lo) / float64(bins)
	total := 0
	for _, c := range counts {
		total += c
	}
	if total == 0 {
		return centers, density
	}
	for i, c := range counts {
		density[i] = float64(c) / (float64(total) * w)
	}
	return centers, density
}

// ShannonEntropy returns the Shannon entropy (in nats) of a discrete
// probability distribution given by counts; the counts are normalised
// internally.
func ShannonEntropy(counts []int) float64 {
	total := 0
	for _, c := range counts {
		total += c
	}
	if total == 0 {
		return 0
	}
	var h float64
	for _, c := range counts {
		if c > 0 {
			p := float64(c) / float64(total)
			h -= p * math.Log(p)
		}
	}
	return h
}

// RandomInitialConditions returns m random points uniformly distributed in the
// box [lo, hi]^dim, drawn from the supplied random source for reproducibility.
func RandomInitialConditions(rng *rand.Rand, m, dim int, lo, hi float64) []Vec {
	out := make([]Vec, m)
	for i := range out {
		v := make(Vec, dim)
		for j := 0; j < dim; j++ {
			v[j] = lo + (hi-lo)*rng.Float64()
		}
		out[i] = v
	}
	return out
}

// Recurrence reports whether two states are within tolerance tol in Euclidean
// distance, the basic predicate used in recurrence analysis.
func Recurrence(a, b Vec, tol float64) bool {
	return a.Distance(b) < tol
}

// RecurrenceRate returns the fraction of point pairs in the trajectory that lie
// within tol of each other (excluding the diagonal), a scalar summary of a
// recurrence plot.
func RecurrenceRate(points []Vec, tol float64) float64 {
	n := len(points)
	if n < 2 {
		return 0
	}
	var count int
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if points[i].Distance(points[j]) < tol {
				count++
			}
		}
	}
	return float64(count) / float64(n*(n-1)/2)
}

// Median returns the median of xs (which is copied and sorted internally).
func Median(xs []float64) float64 {
	if len(xs) == 0 {
		return math.NaN()
	}
	c := append([]float64(nil), xs...)
	sort.Float64s(c)
	n := len(c)
	if n%2 == 1 {
		return c[n/2]
	}
	return 0.5 * (c[n/2-1] + c[n/2])
}

// WrapAngle reduces an angle to the interval [0, 2pi).
func WrapAngle(x float64) float64 {
	y := math.Mod(x, 2*math.Pi)
	if y < 0 {
		y += 2 * math.Pi
	}
	return y
}

// WrapUnit reduces x to the unit interval [0, 1).
func WrapUnit(x float64) float64 {
	return x - math.Floor(x)
}

// LogisticAnalyticLyapunov returns the exact Lyapunov exponent of the fully
// chaotic logistic map at r=4, which equals log 2.
func LogisticAnalyticLyapunov() float64 {
	return math.Ln2
}

// WindingNumber estimates the winding (rotation) number of the circle map with
// parameters omega and k, as the average increment of the lifted angle per
// iteration over n steps after a transient.
func WindingNumber(omega, k, x0 float64, transient, n int) float64 {
	lift := func(x float64) float64 {
		return x + omega - k/(2*math.Pi)*math.Sin(2*math.Pi*x)
	}
	x := x0
	for i := 0; i < transient; i++ {
		x = WrapUnit(lift(x))
	}
	var total float64
	for i := 0; i < n; i++ {
		nx := lift(x)
		total += nx - x
		x = WrapUnit(nx)
	}
	return total / float64(n)
}
