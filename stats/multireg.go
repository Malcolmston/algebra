package stats

import "math"

// solveDense-style dense linear algebra and multivariate modeling for the
// stats package. Everything here is implemented on top of the standard library
// only, using a small internal Gaussian-elimination solver shared by the
// ordinary-least-squares and ridge regressions.

// statsSolveDense solves the dense linear system a·x = b, where a is an n×n
// matrix in row-major order and b has length n, using Gaussian elimination
// with partial pivoting. It reports ok=false when the system is singular (a
// pivot is effectively zero). The elimination order is deterministic: at each
// column the pivot row is the one with the strictly largest absolute value,
// with ties resolved in favor of the earliest row. Both a and b are modified
// in place, so callers must pass scratch copies they no longer need.
func statsSolveDense(a [][]float64, b []float64) ([]float64, bool) {
	n := len(b)
	for col := 0; col < n; col++ {
		// Partial pivot: pick the row with the largest magnitude in this
		// column. Using strict ">" keeps the earliest row on ties, which
		// makes the elimination order deterministic.
		piv := col
		best := math.Abs(a[col][col])
		for r := col + 1; r < n; r++ {
			if v := math.Abs(a[r][col]); v > best {
				best, piv = v, r
			}
		}
		if best < 1e-300 {
			return nil, false
		}
		if piv != col {
			a[col], a[piv] = a[piv], a[col]
			b[col], b[piv] = b[piv], b[col]
		}
		pivVal := a[col][col]
		for r := col + 1; r < n; r++ {
			f := a[r][col] / pivVal
			if f == 0 {
				continue
			}
			for c := col; c < n; c++ {
				a[r][c] -= f * a[col][c]
			}
			b[r] -= f * b[col]
		}
	}
	// Back-substitution.
	x := make([]float64, n)
	for r := n - 1; r >= 0; r-- {
		s := b[r]
		for c := r + 1; c < n; c++ {
			s -= a[r][c] * x[c]
		}
		x[r] = s / a[r][r]
	}
	return x, true
}

// statsMRBuildDesign returns the (optionally intercept-augmented) design matrix
// for the row-major observation matrix X together with its column count. Each
// row of X holds one observation's predictors. When intercept is true a leading
// column of ones is prepended so the intercept is coefficient zero. It reports
// ok=false on a dimension mismatch (len(X) != len(y) or ragged rows), an empty
// or degenerate problem, or too few rows to identify the coefficients
// (fewer rows than columns).
func statsMRBuildDesign(X [][]float64, y []float64, intercept bool) (design [][]float64, cols int, ok bool) {
	n := len(X)
	if n == 0 || n != len(y) {
		return nil, 0, false
	}
	p := len(X[0])
	cols = p
	if intercept {
		cols = p + 1
	}
	if cols == 0 || n < cols {
		return nil, 0, false
	}
	design = make([][]float64, n)
	for i := 0; i < n; i++ {
		if len(X[i]) != p {
			return nil, 0, false
		}
		row := make([]float64, cols)
		off := 0
		if intercept {
			row[0] = 1
			off = 1
		}
		copy(row[off:], X[i])
		design[i] = row
	}
	return design, cols, true
}

// statsMRNormal forms the normal-equation matrices XᵀX and Xᵀy for the given
// design matrix in a single pass over the rows. XᵀX is symmetric, so only its
// upper triangle is accumulated and then mirrored into the lower triangle; this
// roughly halves the inner-product work and is the performance technique this
// file relies on for the regression fits.
func statsMRNormal(design [][]float64, y []float64, cols int) (xtx [][]float64, xty []float64) {
	xtx = make([][]float64, cols)
	for i := range xtx {
		xtx[i] = make([]float64, cols)
	}
	xty = make([]float64, cols)
	for r, row := range design {
		yr := y[r]
		for i := 0; i < cols; i++ {
			ri := row[i]
			xty[i] += ri * yr
			// Accumulate the upper triangle only.
			for j := i; j < cols; j++ {
				xtx[i][j] += ri * row[j]
			}
		}
	}
	// Mirror the upper triangle into the lower triangle.
	for i := 0; i < cols; i++ {
		for j := 0; j < i; j++ {
			xtx[i][j] = xtx[j][i]
		}
	}
	return xtx, xty
}

// statsMRR2 returns the coefficient of determination R² of the fitted
// coefficients on the training data: 1 - SSres/SStot. A constant response
// (SStot == 0) yields 1 for an exact fit and 0 otherwise.
func statsMRR2(design [][]float64, y, b []float64) float64 {
	my := Mean(y)
	var ssRes, ssTot float64
	for r, row := range design {
		pred := 0.0
		for j, v := range row {
			pred += v * b[j]
		}
		e := y[r] - pred
		ssRes += e * e
		d := y[r] - my
		ssTot += d * d
	}
	if ssTot == 0 {
		if ssRes == 0 {
			return 1
		}
		return 0
	}
	return 1 - ssRes/ssTot
}

// MultipleLinearRegression fits an ordinary-least-squares linear model
// predicting y from the predictors in X and returns the fitted coefficients
// together with the coefficient of determination R².
//
// X is row-major: X[i] holds observation i's predictors, and y[i] is its
// response. When intercept is true a leading column of ones is prepended, so
// coeffs[0] is the intercept and coeffs[1:] align with the columns of X;
// otherwise coeffs aligns directly with the columns of X. The fit solves the
// normal equations (XᵀX)·b = Xᵀy by Gaussian elimination with partial
// pivoting, and R² = 1 - SSres/SStot.
//
// It returns (nil, NaN) on a dimension mismatch (len(X) != len(y) or ragged
// rows), too few rows to identify the coefficients, or a singular XᵀX (for
// example perfectly collinear predictors).
func MultipleLinearRegression(X [][]float64, y []float64, intercept bool) (coeffs []float64, r2 float64) {
	design, cols, ok := statsMRBuildDesign(X, y, intercept)
	if !ok {
		return nil, math.NaN()
	}
	xtx, xty := statsMRNormal(design, y, cols)
	b, ok := statsSolveDense(xtx, xty)
	if !ok {
		return nil, math.NaN()
	}
	return b, statsMRR2(design, y, b)
}

// RidgeRegression fits an L2-penalized (ridge) linear model, solving
// (XᵀX + λI)·b = Xᵀy, and returns the fitted coefficients. X follows the same
// row-major convention as [MultipleLinearRegression]. The penalty lambda must
// be non-negative; lambda == 0 reduces exactly to ordinary least squares. When
// intercept is true the intercept column is left unpenalized (its diagonal
// entry of the added λI is zeroed), which is the conventional treatment.
//
// It returns nil on a negative or NaN lambda, a dimension mismatch, too few
// rows, or a singular coefficient matrix.
func RidgeRegression(X [][]float64, y []float64, lambda float64, intercept bool) []float64 {
	if lambda < 0 || math.IsNaN(lambda) {
		return nil
	}
	design, cols, ok := statsMRBuildDesign(X, y, intercept)
	if !ok {
		return nil
	}
	xtx, xty := statsMRNormal(design, y, cols)
	if lambda != 0 {
		for i := 0; i < cols; i++ {
			if intercept && i == 0 {
				// Do not penalize the intercept column.
				continue
			}
			xtx[i][i] += lambda
		}
	}
	b, ok := statsSolveDense(xtx, xty)
	if !ok {
		return nil
	}
	return b
}

// Predict returns the linear prediction for a single predictor row x given a
// coefficient vector, honoring the same intercept convention as the fitting
// functions. When intercept is true the prediction is coeffs[0] + Σ
// coeffs[1+j]·x[j]; otherwise it is Σ coeffs[j]·x[j]. It returns NaN if x is
// not consistent with coeffs (its length must equal len(coeffs), or
// len(coeffs)-1 when intercept is true).
func Predict(coeffs []float64, x []float64, intercept bool) float64 {
	off := 0
	sum := 0.0
	if intercept {
		if len(coeffs) == 0 {
			return math.NaN()
		}
		sum = coeffs[0]
		off = 1
	}
	if len(x) != len(coeffs)-off {
		return math.NaN()
	}
	for j, v := range x {
		sum += coeffs[off+j] * v
	}
	return sum
}

// statsMRColumnsValid reports whether cols describes a well-formed set of
// variables for a covariance or correlation matrix: at least one column, every
// column of the same length n, and n >= 2 so the unbiased (n-1) estimators are
// defined.
func statsMRColumnsValid(cols [][]float64) bool {
	k := len(cols)
	if k == 0 {
		return false
	}
	n := len(cols[0])
	if n < 2 {
		return false
	}
	for _, c := range cols {
		if len(c) != n {
			return false
		}
	}
	return true
}

// statsMRNaNMatrix returns a k×k matrix whose every entry is NaN.
func statsMRNaNMatrix(k int) [][]float64 {
	out := make([][]float64, k)
	for i := range out {
		out[i] = make([]float64, k)
		for j := range out[i] {
			out[i][j] = math.NaN()
		}
	}
	return out
}

// CovarianceMatrix returns the symmetric k×k unbiased (n-1) sample covariance
// matrix for k variables, where cols[j] holds the samples of variable j and all
// columns share the same length n. The (i, j) entry equals the sample
// covariance of cols[i] and cols[j], so the diagonal equals [Variance] of the
// corresponding column.
//
// It returns nil if there are no columns. If the columns differ in length or
// have fewer than two samples, it returns a k×k matrix whose entries are all
// NaN.
func CovarianceMatrix(cols [][]float64) [][]float64 {
	k := len(cols)
	if k == 0 {
		return nil
	}
	if !statsMRColumnsValid(cols) {
		return statsMRNaNMatrix(k)
	}
	out := make([][]float64, k)
	for i := range out {
		out[i] = make([]float64, k)
	}
	for i := 0; i < k; i++ {
		out[i][i] = Variance(cols[i])
		for j := i + 1; j < k; j++ {
			c := Covariance(cols[i], cols[j])
			out[i][j] = c
			out[j][i] = c
		}
	}
	return out
}

// CorrelationMatrix returns the symmetric k×k Pearson correlation matrix for k
// variables, where cols[j] holds the samples of variable j. The diagonal is 1
// and each off-diagonal (i, j) entry equals [Correlation] of cols[i] and
// cols[j] (which is NaN when either variable has zero variance).
//
// It returns nil if there are no columns. If the columns differ in length or
// have fewer than two samples, it returns a k×k matrix whose entries are all
// NaN.
func CorrelationMatrix(cols [][]float64) [][]float64 {
	k := len(cols)
	if k == 0 {
		return nil
	}
	if !statsMRColumnsValid(cols) {
		return statsMRNaNMatrix(k)
	}
	out := make([][]float64, k)
	for i := range out {
		out[i] = make([]float64, k)
	}
	for i := 0; i < k; i++ {
		out[i][i] = 1
		for j := i + 1; j < k; j++ {
			c := Correlation(cols[i], cols[j])
			out[i][j] = c
			out[j][i] = c
		}
	}
	return out
}
