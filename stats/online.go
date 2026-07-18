package stats

import "math"

// Accumulator is a constant-memory, single-pass streaming replacement for the
// multi-pass descriptive functions ([Mean], [Variance], [Skewness],
// [Kurtosis] and friends). It maintains a running count, mean and the central
// moments M2, M3 and M4 using the numerically stable Welford/Terriberry
// recurrences, together with a running minimum and maximum, all updated in the
// same pass.
//
// The zero value is an empty accumulator ready for use. After construction no
// method allocates, none sorts, and each observation is visited exactly once,
// making Accumulator suitable for arbitrarily large or unbounded streams and,
// via [Accumulator.Merge], for map-reduce over shards processed in parallel.
//
// The summary statistics are numerically consistent with the corresponding
// batch functions in this package: [Accumulator.Variance] matches [Variance],
// [Accumulator.Skewness] matches [Skewness], and so on. Undefined results are
// reported as NaN under the same conditions the batch functions use.
type Accumulator struct {
	n          int
	mean       float64
	m2, m3, m4 float64
	min, max   float64
}

// Push incorporates a single observation x into the accumulator, updating the
// count, mean, central moments M2, M3 and M4, and the running minimum and
// maximum in one pass. It uses the numerically stable Welford/Terriberry
// recurrences, which avoid the cancellation error of the naive
// sum-of-squares approach.
func (a *Accumulator) Push(x float64) {
	if a.n == 0 {
		a.min = x
		a.max = x
	} else {
		if x < a.min {
			a.min = x
		}
		if x > a.max {
			a.max = x
		}
	}

	n1 := float64(a.n)
	a.n++
	n := float64(a.n)

	delta := x - a.mean
	deltaN := delta / n
	deltaN2 := deltaN * deltaN
	term1 := delta * deltaN * n1

	a.mean += deltaN
	// M4 and M3 must be updated using the old M2 and M3, so update them
	// before M2 (and M4 before M3).
	a.m4 += term1*deltaN2*(n*n-3*n+3) + 6*deltaN2*a.m2 - 4*deltaN*a.m3
	a.m3 += term1*deltaN*(n-2) - 3*deltaN*a.m2
	a.m2 += term1
}

// PushAll incorporates every observation in xs by calling [Accumulator.Push]
// on each in order. It allocates nothing.
func (a *Accumulator) PushAll(xs []float64) {
	for _, x := range xs {
		a.Push(x)
	}
}

// Count returns the number of observations pushed so far.
func (a *Accumulator) Count() int { return a.n }

// Sum returns the sum of all observations, computed as the mean times the
// count. The sum of an empty accumulator is 0, matching [Sum].
func (a *Accumulator) Sum() float64 { return a.mean * float64(a.n) }

// Mean returns the arithmetic mean of the observations, or NaN if none have
// been pushed. It matches [Mean].
func (a *Accumulator) Mean() float64 {
	if a.n < 1 {
		return math.NaN()
	}
	return a.mean
}

// Variance returns the unbiased sample variance using the n-1 (Bessel)
// denominator, or NaN if fewer than two observations have been pushed. It
// matches [Variance].
func (a *Accumulator) Variance() float64 {
	if a.n < 2 {
		return math.NaN()
	}
	return a.m2 / float64(a.n-1)
}

// PopVariance returns the population variance using the n denominator, or NaN
// if no observations have been pushed. It matches [PopVariance].
func (a *Accumulator) PopVariance() float64 {
	if a.n < 1 {
		return math.NaN()
	}
	return a.m2 / float64(a.n)
}

// StdDev returns the sample standard deviation (the square root of
// [Accumulator.Variance]), or NaN if fewer than two observations have been
// pushed. It matches [StdDev].
func (a *Accumulator) StdDev() float64 { return math.Sqrt(a.Variance()) }

// PopStdDev returns the population standard deviation (the square root of
// [Accumulator.PopVariance]), or NaN if no observations have been pushed. It
// matches [PopStdDev].
func (a *Accumulator) PopStdDev() float64 { return math.Sqrt(a.PopVariance()) }

// Skewness returns the population coefficient of skewness (the third
// standardized moment), computed from the running moments as
// sqrt(n)*M3/M2^1.5. It returns NaN if fewer than two observations have been
// pushed or the variance is zero. It matches [Skewness].
func (a *Accumulator) Skewness() float64 {
	if a.n < 2 || a.m2 == 0 {
		return math.NaN()
	}
	return math.Sqrt(float64(a.n)) * a.m3 / math.Pow(a.m2, 1.5)
}

// Kurtosis returns the population excess kurtosis (the fourth standardized
// moment minus 3, so a normal distribution has kurtosis 0), computed from the
// running moments as n*M4/M2^2-3. It returns NaN if fewer than two
// observations have been pushed or the variance is zero. It matches
// [Kurtosis].
func (a *Accumulator) Kurtosis() float64 {
	if a.n < 2 || a.m2 == 0 {
		return math.NaN()
	}
	return float64(a.n)*a.m4/(a.m2*a.m2) - 3
}

// Min returns the smallest observation pushed so far, or NaN if none have been
// pushed. It matches [Min].
func (a *Accumulator) Min() float64 {
	if a.n < 1 {
		return math.NaN()
	}
	return a.min
}

// Max returns the largest observation pushed so far, or NaN if none have been
// pushed. It matches [Max].
func (a *Accumulator) Max() float64 {
	if a.n < 1 {
		return math.NaN()
	}
	return a.max
}

// Merge folds the observations summarized by o into the receiver as if they
// had been pushed directly, using the pairwise-moment merge formulas of Chan,
// Golub and LeVeque (extended to the third and fourth moments by Terriberry).
// This enables a map-reduce style computation in which disjoint shards of a
// stream are summarized by separate accumulators in parallel and then combined
// exactly, with the same result as a single sequential pass.
func (a *Accumulator) Merge(o Accumulator) {
	if o.n == 0 {
		return
	}
	if a.n == 0 {
		*a = o
		return
	}

	na := float64(a.n)
	nb := float64(o.n)
	n := na + nb

	delta := o.mean - a.mean
	delta2 := delta * delta
	delta3 := delta2 * delta
	delta4 := delta2 * delta2

	m2 := a.m2 + o.m2 + delta2*na*nb/n
	m3 := a.m3 + o.m3 +
		delta3*na*nb*(na-nb)/(n*n) +
		3*delta*(na*o.m2-nb*a.m2)/n
	m4 := a.m4 + o.m4 +
		delta4*na*nb*(na*na-na*nb+nb*nb)/(n*n*n) +
		6*delta2*(na*na*o.m2+nb*nb*a.m2)/(n*n) +
		4*delta*(na*o.m3-nb*a.m3)/n

	a.mean += delta * nb / n
	a.m2 = m2
	a.m3 = m3
	a.m4 = m4
	a.n = int(n)

	if o.min < a.min {
		a.min = o.min
	}
	if o.max > a.max {
		a.max = o.max
	}
}

// Reset returns the accumulator to its empty zero-value state so it can be
// reused without allocating.
func (a *Accumulator) Reset() { *a = Accumulator{} }

// CovAccumulator is a constant-memory, single-pass streaming replacement for
// the multi-pass bivariate functions [Covariance], [Correlation] and
// [LinearRegression]. It maintains the running means of x and y, their central
// second moments, and the co-moment C = sum((x-meanX)(y-meanY)) using the
// numerically stable Welford recurrences.
//
// The zero value is an empty accumulator ready for use. After construction no
// method allocates and each observation is visited exactly once, making
// CovAccumulator suitable for large or unbounded paired streams and, via
// [CovAccumulator.Merge], for map-reduce over shards. Its results are
// numerically consistent with the corresponding batch functions.
type CovAccumulator struct {
	n            int
	meanX, meanY float64
	m2x, m2y     float64
	c            float64
}

// Push incorporates a single paired observation (x, y), updating the running
// means, the second moments of x and y, and the co-moment C in one pass using
// the numerically stable Welford recurrences.
func (a *CovAccumulator) Push(x, y float64) {
	a.n++
	n := float64(a.n)

	dx := x - a.meanX
	dy := y - a.meanY
	a.meanX += dx / n
	a.meanY += dy / n
	// Use the deviations against the pre-update mean (dx, dy) paired with the
	// post-update mean for the second factor, the standard Welford form.
	a.c += dx * (y - a.meanY)
	a.m2x += dx * (x - a.meanX)
	a.m2y += dy * (y - a.meanY)
}

// PushAll incorporates every paired observation in xs and ys by calling
// [CovAccumulator.Push] on each pair in order. It processes min(len(xs),
// len(ys)) pairs and allocates nothing.
func (a *CovAccumulator) PushAll(xs, ys []float64) {
	n := len(xs)
	if len(ys) < n {
		n = len(ys)
	}
	for i := 0; i < n; i++ {
		a.Push(xs[i], ys[i])
	}
}

// Count returns the number of paired observations pushed so far.
func (a *CovAccumulator) Count() int { return a.n }

// Covariance returns the unbiased sample covariance between x and y using the
// n-1 denominator, or NaN if fewer than two pairs have been pushed. It matches
// [Covariance].
func (a *CovAccumulator) Covariance() float64 {
	if a.n < 2 {
		return math.NaN()
	}
	return a.c / float64(a.n-1)
}

// Correlation returns the Pearson product-moment correlation coefficient
// between x and y, a value in [-1, 1]. It returns NaN if fewer than two pairs
// have been pushed or either variable has zero variance. It matches
// [Correlation].
func (a *CovAccumulator) Correlation() float64 {
	if a.n < 2 {
		return math.NaN()
	}
	den := math.Sqrt(a.m2x * a.m2y)
	if den == 0 {
		return math.NaN()
	}
	return a.c / den
}

// Slope returns the slope of the ordinary-least-squares line y = slope*x +
// intercept fitted to the pushed pairs, equal to the slope reported by
// [LinearRegression]. It returns NaN if fewer than two pairs have been pushed
// or x has zero variance.
func (a *CovAccumulator) Slope() float64 {
	if a.n < 2 || a.m2x == 0 {
		return math.NaN()
	}
	return a.c / a.m2x
}

// Intercept returns the intercept of the ordinary-least-squares line y =
// slope*x + intercept fitted to the pushed pairs, equal to the intercept
// reported by [LinearRegression]. It returns NaN if fewer than two pairs have
// been pushed or x has zero variance.
func (a *CovAccumulator) Intercept() float64 {
	if a.n < 2 || a.m2x == 0 {
		return math.NaN()
	}
	return a.meanY - (a.c/a.m2x)*a.meanX
}

// Merge folds the paired observations summarized by o into the receiver as if
// they had been pushed directly, using the pairwise co-moment merge formula.
// This enables map-reduce over disjoint shards summarized in parallel, with
// the same result as a single sequential pass.
func (a *CovAccumulator) Merge(o CovAccumulator) {
	if o.n == 0 {
		return
	}
	if a.n == 0 {
		*a = o
		return
	}

	na := float64(a.n)
	nb := float64(o.n)
	n := na + nb

	dx := o.meanX - a.meanX
	dy := o.meanY - a.meanY

	a.c += o.c + dx*dy*na*nb/n
	a.m2x += o.m2x + dx*dx*na*nb/n
	a.m2y += o.m2y + dy*dy*na*nb/n
	a.meanX += dx * nb / n
	a.meanY += dy * nb / n
	a.n = int(n)
}

// Reset returns the accumulator to its empty zero-value state so it can be
// reused without allocating.
func (a *CovAccumulator) Reset() { *a = CovAccumulator{} }
