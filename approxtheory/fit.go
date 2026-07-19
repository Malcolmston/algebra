package approxtheory

import "math"

// PolyFitResult holds the result of a polynomial fit: the monomial
// coefficients (ascending), the root-mean-square residual, the maximum
// absolute residual and the coefficient of determination R^2.
type PolyFitResult struct {
	Coeffs []float64
	RMS    float64
	MaxAbs float64
	R2     float64
}

// Eval evaluates the fitted polynomial at x.
func (r *PolyFitResult) Eval(x float64) float64 { return Polyval(r.Coeffs, x) }

// PolyFit computes the ordinary least-squares polynomial of the given degree
// through the data (xs, ys) by solving the normal equations. It returns the
// coefficients and residual statistics.
func PolyFit(xs, ys []float64, degree int) (*PolyFitResult, error) {
	return WeightedPolyFit(xs, ys, nil, degree)
}

// WeightedPolyFit computes the weighted least-squares polynomial of the given
// degree. When weights is nil every point has unit weight; otherwise weights
// must match the length of xs and ys and are treated as multipliers on the
// squared residuals.
func WeightedPolyFit(xs, ys, weights []float64, degree int) (*PolyFitResult, error) {
	n := len(xs)
	if n == 0 {
		return nil, ErrEmptyInput
	}
	if len(ys) != n || (weights != nil && len(weights) != n) {
		return nil, ErrDimensionMismatch
	}
	if degree < 0 {
		degree = 0
	}
	p := degree + 1
	if p > n {
		return nil, ErrDimensionMismatch
	}
	// Build normal equations A^T W A c = A^T W y using power sums.
	// Precompute powers of x up to 2*degree.
	sums := make([]float64, 2*degree+1)
	rhs := make([]float64, p)
	for i := 0; i < n; i++ {
		w := 1.0
		if weights != nil {
			w = weights[i]
		}
		xp := 1.0
		powers := make([]float64, 2*degree+1)
		for k := 0; k <= 2*degree; k++ {
			powers[k] = xp
			xp *= xs[i]
		}
		for k := 0; k <= 2*degree; k++ {
			sums[k] += w * powers[k]
		}
		for r := 0; r < p; r++ {
			rhs[r] += w * powers[r] * ys[i]
		}
	}
	A := make([][]float64, p)
	for r := 0; r < p; r++ {
		row := make([]float64, p)
		for c := 0; c < p; c++ {
			row[c] = sums[r+c]
		}
		A[r] = row
	}
	coeffs, err := solveLinear(A, rhs)
	if err != nil {
		return nil, err
	}
	res := &PolyFitResult{Coeffs: coeffs}
	fillFitStats(res, xs, ys)
	return res, nil
}

// fillFitStats computes RMS, MaxAbs and R2 for a fitted polynomial.
func fillFitStats(res *PolyFitResult, xs, ys []float64) {
	n := len(xs)
	var mean float64
	for _, y := range ys {
		mean += y
	}
	mean /= float64(n)
	var ssRes, ssTot, maxAbs float64
	for i := 0; i < n; i++ {
		r := ys[i] - Polyval(res.Coeffs, xs[i])
		ssRes += r * r
		ssTot += (ys[i] - mean) * (ys[i] - mean)
		if a := math.Abs(r); a > maxAbs {
			maxAbs = a
		}
	}
	res.RMS = math.Sqrt(ssRes / float64(n))
	res.MaxAbs = maxAbs
	if ssTot > 0 {
		res.R2 = 1 - ssRes/ssTot
	} else {
		res.R2 = 1
	}
}

// RSquared returns the coefficient of determination for the fit of a monomial
// polynomial to the data (xs, ys).
func RSquared(coeffs, xs, ys []float64) float64 {
	n := len(ys)
	if n == 0 {
		return 0
	}
	var mean float64
	for _, y := range ys {
		mean += y
	}
	mean /= float64(n)
	var ssRes, ssTot float64
	for i := 0; i < n; i++ {
		r := ys[i] - Polyval(coeffs, xs[i])
		ssRes += r * r
		ssTot += (ys[i] - mean) * (ys[i] - mean)
	}
	if ssTot == 0 {
		return 1
	}
	return 1 - ssRes/ssTot
}

// ChebLeastSquares fits a Chebyshev series of the given degree to sampled data
// (xs, ys) on [a, b] by least squares in the Chebyshev basis. It is more
// numerically stable than a monomial fit for higher degrees.
func ChebLeastSquares(xs, ys []float64, degree int, a, b float64) (*ChebSeries, error) {
	n := len(xs)
	if n == 0 {
		return nil, ErrEmptyInput
	}
	if len(ys) != n {
		return nil, ErrDimensionMismatch
	}
	if degree < 0 {
		degree = 0
	}
	p := degree + 1
	// Design matrix rows are Chebyshev values at mapped points.
	At := make([][]float64, p) // p x p normal matrix
	for r := range At {
		At[r] = make([]float64, p)
	}
	rhs := make([]float64, p)
	for i := 0; i < n; i++ {
		t := (2*xs[i] - (a + b)) / (b - a)
		tv := ChebTValues(degree, t)
		for r := 0; r < p; r++ {
			rhs[r] += tv[r] * ys[i]
			for c := 0; c < p; c++ {
				At[r][c] += tv[r] * tv[c]
			}
		}
	}
	coeffs, err := solveLinear(At, rhs)
	if err != nil {
		return nil, err
	}
	return &ChebSeries{Coeffs: coeffs, A: a, B: b}, nil
}

// DiscreteMinimaxPoly computes a degree-n polynomial minimising the maximum
// absolute residual over the discrete data set (xs, ys) using a single-point
// exchange (Remez-type) algorithm. It returns the monomial coefficients and
// the achieved maximum error.
func DiscreteMinimaxPoly(xs, ys []float64, n, maxIter int) ([]float64, float64, error) {
	N := len(xs)
	if N == 0 {
		return nil, 0, ErrEmptyInput
	}
	if len(ys) != N {
		return nil, 0, ErrDimensionMismatch
	}
	if n < 0 {
		n = 0
	}
	m := n + 2
	if N < m {
		return nil, 0, ErrDimensionMismatch
	}
	if maxIter <= 0 {
		maxIter = 100
	}
	// Initial reference: m points spread across the sorted index range.
	order := sortedIndex(xs)
	ref := make([]int, m)
	for i := 0; i < m; i++ {
		ref[i] = order[i*(N-1)/(m-1)]
	}
	var coeffs []float64
	var maxErr float64
	for it := 0; it < maxIter; it++ {
		c, _, err := discreteRemezSolve(xs, ys, ref, n)
		if err != nil {
			return nil, 0, err
		}
		coeffs = c
		// Find the data point of maximum residual.
		worst := -1
		maxErr = 0
		var worstResid float64
		for i := 0; i < N; i++ {
			r := ys[i] - Polyval(coeffs, xs[i])
			if math.Abs(r) > maxErr {
				maxErr = math.Abs(r)
				worst = i
				worstResid = r
			}
		}
		// If the worst point is already in the reference (within the leveled
		// error) we are done.
		inRef := false
		for _, idx := range ref {
			if idx == worst {
				inRef = true
				break
			}
		}
		if inRef || worst < 0 {
			break
		}
		// Exchange: replace the reference point nearest to worst that keeps the
		// sign alternation, i.e. the closest reference index in xs.
		pos := 0
		bestDist := math.Inf(1)
		for j, idx := range ref {
			d := math.Abs(xs[idx] - xs[worst])
			if d < bestDist {
				bestDist = d
				pos = j
			}
		}
		// Preserve alternation: keep the sign of residual at that reference
		// slot matching worstResid's sign.
		refResid := ys[ref[pos]] - Polyval(coeffs, xs[ref[pos]])
		if refResid*worstResid < 0 {
			// choose an adjacent slot instead to keep alternation
			if pos+1 < m {
				pos++
			} else if pos-1 >= 0 {
				pos--
			}
		}
		ref[pos] = worst
		sortRefByX(ref, xs)
	}
	return coeffs, maxErr, nil
}

// discreteRemezSolve solves the alternating linear system on m = n+2 discrete
// reference indices.
func discreteRemezSolve(xs, ys []float64, ref []int, n int) ([]float64, float64, error) {
	m := len(ref)
	A := make([][]float64, m)
	rhs := make([]float64, m)
	for i := 0; i < m; i++ {
		row := make([]float64, m)
		xp := 1.0
		for j := 0; j <= n; j++ {
			row[j] = xp
			xp *= xs[ref[i]]
		}
		sign := 1.0
		if i%2 == 1 {
			sign = -1.0
		}
		row[n+1] = sign
		A[i] = row
		rhs[i] = ys[ref[i]]
	}
	sol, err := solveLinear(A, rhs)
	if err != nil {
		return nil, 0, err
	}
	c := make([]float64, n+1)
	copy(c, sol[:n+1])
	return c, sol[n+1], nil
}

// sortedIndex returns indices that sort xs ascending.
func sortedIndex(xs []float64) []int {
	idx := make([]int, len(xs))
	for i := range idx {
		idx[i] = i
	}
	// simple insertion sort keeps dependencies minimal
	for i := 1; i < len(idx); i++ {
		j := i
		for j > 0 && xs[idx[j-1]] > xs[idx[j]] {
			idx[j-1], idx[j] = idx[j], idx[j-1]
			j--
		}
	}
	return idx
}

// sortRefByX sorts a slice of indices so that xs[ref] is ascending.
func sortRefByX(ref []int, xs []float64) {
	for i := 1; i < len(ref); i++ {
		j := i
		for j > 0 && xs[ref[j-1]] > xs[ref[j]] {
			ref[j-1], ref[j] = ref[j], ref[j-1]
			j--
		}
	}
}
