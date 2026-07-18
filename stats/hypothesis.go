package stats

import (
	"math"
	"sort"
)

// TestResult is the common outcome of the classical significance tests in this
// file. Statistic is the test statistic (a t, F, chi-squared or z value
// depending on the test), DF is the degrees of freedom of the reference
// distribution (for tests with two degrees-of-freedom parameters it holds the
// numerator/between value, documented per function), and PValue is the p-value
// under the null hypothesis. Degenerate inputs yield a result whose fields are
// all NaN.
type TestResult struct {
	Statistic float64
	DF        float64
	PValue    float64
}

// statsNaNResult returns a TestResult whose fields are all NaN, used to signal
// degenerate input to any of the tests in this file.
func statsNaNResult() TestResult {
	nan := math.NaN()
	return TestResult{Statistic: nan, DF: nan, PValue: nan}
}

// statsTwoSidedTP returns the two-sided p-value for a Student's t statistic t
// on df degrees of freedom, computed as 2*(1-StudentT{df}.CDF(|t|)).
func statsTwoSidedTP(t, df float64) float64 {
	return 2 * (1 - StudentT{Nu: df}.CDF(math.Abs(t)))
}

// OneSampleTTest performs a two-sided one-sample Student's t-test of the null
// hypothesis that the population mean of xs equals mu0. The statistic is
// t=(mean-mu0)/(s/sqrt(n)) with s the sample standard deviation, the degrees of
// freedom are n-1, and the p-value is 2*(1-StudentT{df}.CDF(|t|)). It returns a
// NaN result when xs has fewer than two elements or zero sample variance.
func OneSampleTTest(xs []float64, mu0 float64) TestResult {
	n := len(xs)
	if n < 2 {
		return statsNaNResult()
	}
	sd := StdDev(xs)
	df := float64(n - 1)
	if sd == 0 || math.IsNaN(sd) {
		return statsNaNResult()
	}
	se := sd / math.Sqrt(float64(n))
	t := (Mean(xs) - mu0) / se
	return TestResult{Statistic: t, DF: df, PValue: statsTwoSidedTP(t, df)}
}

// TwoSampleTTest performs a two-sided two-sample t-test of the null hypothesis
// that xs and ys share a common mean. When equalVar is true it uses the
// pooled-variance (Student's) test with n1+n2-2 degrees of freedom; otherwise
// it uses Welch's test with the Welch–Satterthwaite degrees of freedom. The
// p-value is 2*(1-StudentT{df}.CDF(|t|)). It returns a NaN result when either
// sample has fewer than two elements or the standard error is zero.
func TwoSampleTTest(xs, ys []float64, equalVar bool) TestResult {
	n1, n2 := len(xs), len(ys)
	if n1 < 2 || n2 < 2 {
		return statsNaNResult()
	}
	m1, m2 := Mean(xs), Mean(ys)
	v1, v2 := Variance(xs), Variance(ys)
	f1, f2 := float64(n1), float64(n2)
	var t, df float64
	if equalVar {
		sp2 := ((f1-1)*v1 + (f2-1)*v2) / (f1 + f2 - 2)
		se := math.Sqrt(sp2 * (1/f1 + 1/f2))
		if se == 0 {
			return statsNaNResult()
		}
		t = (m1 - m2) / se
		df = f1 + f2 - 2
	} else {
		se := math.Sqrt(v1/f1 + v2/f2)
		if se == 0 {
			return statsNaNResult()
		}
		t = (m1 - m2) / se
		num := (v1/f1 + v2/f2) * (v1/f1 + v2/f2)
		den := (v1/f1)*(v1/f1)/(f1-1) + (v2/f2)*(v2/f2)/(f2-1)
		df = num / den
	}
	return TestResult{Statistic: t, DF: df, PValue: statsTwoSidedTP(t, df)}
}

// PairedTTest performs a two-sided paired-sample t-test of the null hypothesis
// that the mean of the paired differences xs[i]-ys[i] is zero. It is defined as
// OneSampleTTest applied to the differences with mu0=0, so the degrees of
// freedom are n-1. It returns a NaN result when the slices differ in length,
// have fewer than two pairs, or the differences have zero variance.
func PairedTTest(xs, ys []float64) TestResult {
	n := len(xs)
	if n != len(ys) || n < 2 {
		return statsNaNResult()
	}
	diff := make([]float64, n)
	for i := range xs {
		diff[i] = xs[i] - ys[i]
	}
	return OneSampleTTest(diff, 0)
}

// ChiSquareGoodnessOfFit performs Pearson's chi-squared goodness-of-fit test
// comparing observed counts against expected counts. The statistic is
// X2=sum((o-e)^2/e), the degrees of freedom are k-1 for k categories, and the
// p-value is 1-ChiSquared{df}.CDF(X2). It returns a NaN result when the slices
// differ in length, have fewer than two categories, or any expected count is
// non-positive.
func ChiSquareGoodnessOfFit(observed, expected []float64) TestResult {
	k := len(observed)
	if k != len(expected) || k < 2 {
		return statsNaNResult()
	}
	var x2 float64
	for i := 0; i < k; i++ {
		if expected[i] <= 0 {
			return statsNaNResult()
		}
		d := observed[i] - expected[i]
		x2 += d * d / expected[i]
	}
	df := float64(k - 1)
	return TestResult{Statistic: x2, DF: df, PValue: 1 - ChiSquared{K: df}.CDF(x2)}
}

// ChiSquareIndependence performs Pearson's chi-squared test of independence on
// a contingency table given as a slice of equal-length rows of counts. Expected
// cell counts are formed from the row and column totals as rowTotal*colTotal/
// grandTotal, the statistic is X2=sum((o-e)^2/e), and the degrees of freedom
// are (r-1)(c-1). The p-value is 1-ChiSquared{df}.CDF(X2). It returns a NaN
// result when the table is smaller than 2x2, ragged, or has a zero grand total.
func ChiSquareIndependence(table [][]float64) TestResult {
	rows := len(table)
	if rows < 2 {
		return statsNaNResult()
	}
	cols := len(table[0])
	if cols < 2 {
		return statsNaNResult()
	}
	rowSum := make([]float64, rows)
	colSum := make([]float64, cols)
	var total float64
	for i := 0; i < rows; i++ {
		if len(table[i]) != cols {
			return statsNaNResult()
		}
		for j := 0; j < cols; j++ {
			v := table[i][j]
			rowSum[i] += v
			colSum[j] += v
			total += v
		}
	}
	if total == 0 {
		return statsNaNResult()
	}
	var x2 float64
	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			exp := rowSum[i] * colSum[j] / total
			if exp <= 0 {
				continue
			}
			d := table[i][j] - exp
			x2 += d * d / exp
		}
	}
	df := float64((rows - 1) * (cols - 1))
	return TestResult{Statistic: x2, DF: df, PValue: 1 - ChiSquared{K: df}.CDF(x2)}
}

// OneWayANOVA performs a one-way analysis of variance testing the null
// hypothesis that all supplied groups share a common mean. It computes the
// between-groups and within-groups sums of squares, F=(SSB/dfB)/(SSW/dfW) with
// dfB=k-1 and dfW=N-k, and the upper-tail p-value 1-FDist{dfB,dfW}.CDF(F). The
// returned DF field holds the numerator (between-groups) degrees of freedom.
// It returns a NaN result when fewer than two non-empty groups are supplied,
// the total size does not exceed the group count, or the within-groups sum of
// squares is zero.
func OneWayANOVA(groups ...[]float64) TestResult {
	var k, ntot int
	var grand float64
	for _, g := range groups {
		if len(g) == 0 {
			continue
		}
		k++
		ntot += len(g)
		grand += Sum(g)
	}
	if k < 2 || ntot <= k {
		return statsNaNResult()
	}
	grand /= float64(ntot)
	var ssBetween, ssWithin float64
	for _, g := range groups {
		if len(g) == 0 {
			continue
		}
		gm := Mean(g)
		d := gm - grand
		ssBetween += float64(len(g)) * d * d
		for _, x := range g {
			e := x - gm
			ssWithin += e * e
		}
	}
	dfB := float64(k - 1)
	dfW := float64(ntot - k)
	if ssWithin == 0 {
		return statsNaNResult()
	}
	f := (ssBetween / dfB) / (ssWithin / dfW)
	return TestResult{Statistic: f, DF: dfB, PValue: 1 - FDist{D1: dfB, D2: dfW}.CDF(f)}
}

// statsRankPair pairs a value with the group it belongs to (0 for the first
// sample, 1 for the second) so that the two samples can be ranked jointly.
type statsRankPair struct {
	value float64
	group int
}

// MannWhitneyU performs the two-sided Mann–Whitney U test (equivalently the
// Wilcoxon rank-sum test) of the null hypothesis that xs and ys come from the
// same distribution. Both samples are ranked jointly using average ranks for
// ties; U is the smaller of the two group U values. The p-value uses the normal
// approximation with a continuity correction and a tie correction to the
// variance: z=(U-mu_U)/sigma_U and PValue=2*(1-Normal{0,1}.CDF(|z|)). The
// returned DF field is NaN because the normal approximation has no degrees of
// freedom. Ranking uses a stable sort of index pairs, making the result
// deterministic. It returns a NaN result when either sample is empty or the
// approximate variance is non-positive (for example, when every observation is
// tied).
func MannWhitneyU(xs, ys []float64) TestResult {
	n1, n2 := len(xs), len(ys)
	if n1 == 0 || n2 == 0 {
		return statsNaNResult()
	}
	pairs := make([]statsRankPair, 0, n1+n2)
	for _, v := range xs {
		pairs = append(pairs, statsRankPair{value: v, group: 0})
	}
	for _, v := range ys {
		pairs = append(pairs, statsRankPair{value: v, group: 1})
	}

	// Stable sort of index pairs: sort indices into pairs by value, keeping the
	// original order for equal values so ranking is fully deterministic.
	n := n1 + n2
	idx := make([]int, n)
	for i := range idx {
		idx[i] = i
	}
	sort.SliceStable(idx, func(a, b int) bool {
		return pairs[idx[a]].value < pairs[idx[b]].value
	})

	// Assign average ranks over runs of equal values, accumulating the tie
	// correction term sum(t^3 - t) and the rank sum of the first sample.
	var r1, tieTerm float64
	for i := 0; i < n; {
		j := i
		for j < n && pairs[idx[j]].value == pairs[idx[i]].value {
			j++
		}
		avg := float64(i+j+1) / 2 // mean of ranks (i+1)..j, 1-based
		for m := i; m < j; m++ {
			if pairs[idx[m]].group == 0 {
				r1 += avg
			}
		}
		t := float64(j - i)
		tieTerm += t*t*t - t
		i = j
	}

	f1, f2 := float64(n1), float64(n2)
	u1 := r1 - f1*(f1+1)/2
	u2 := f1*f2 - u1
	u := math.Min(u1, u2)

	nn := f1 + f2
	muU := f1 * f2 / 2
	varU := f1 * f2 / 12 * ((nn + 1) - tieTerm/(nn*(nn-1)))
	if varU <= 0 {
		return statsNaNResult()
	}
	// Continuity correction of 0.5 toward the mean, applied to the two-sided
	// deviation |U - mu_U|.
	z := (math.Abs(u-muU) - 0.5) / math.Sqrt(varU)
	if z < 0 {
		z = 0
	}
	p := 2 * (1 - Normal{Mu: 0, Sigma: 1}.CDF(z))
	return TestResult{Statistic: u, DF: math.NaN(), PValue: p}
}
