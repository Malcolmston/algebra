package timeseries

import "math"

// mean returns the arithmetic mean of x, or NaN for an empty slice.
func mean(x []float64) float64 {
	if len(x) == 0 {
		return math.NaN()
	}
	var s float64
	for _, v := range x {
		s += v
	}
	return s / float64(len(x))
}

// sumf returns the sum of x.
func sumf(x []float64) float64 {
	var s float64
	for _, v := range x {
		s += v
	}
	return s
}

// copyf returns a fresh copy of x.
func copyf(x []float64) []float64 {
	out := make([]float64, len(x))
	copy(out, x)
	return out
}

// solveLinear solves A·z = b for z using Gaussian elimination with partial
// pivoting. A is modified in place. It returns the solution and true, or nil
// and false if the system is singular. A must be square with len(A)==len(b).
func solveLinear(A [][]float64, b []float64) ([]float64, bool) {
	n := len(A)
	if n == 0 || len(b) != n {
		return nil, false
	}
	// Work on copies so callers keep their data.
	M := make([][]float64, n)
	for i := range A {
		if len(A[i]) != n {
			return nil, false
		}
		M[i] = make([]float64, n+1)
		copy(M[i], A[i])
		M[i][n] = b[i]
	}
	for col := 0; col < n; col++ {
		// Partial pivot.
		piv := col
		best := math.Abs(M[col][col])
		for r := col + 1; r < n; r++ {
			if a := math.Abs(M[r][col]); a > best {
				best = a
				piv = r
			}
		}
		if best == 0 {
			return nil, false
		}
		M[col], M[piv] = M[piv], M[col]
		// Eliminate.
		for r := 0; r < n; r++ {
			if r == col {
				continue
			}
			f := M[r][col] / M[col][col]
			if f == 0 {
				continue
			}
			for c := col; c <= n; c++ {
				M[r][c] -= f * M[col][c]
			}
		}
	}
	z := make([]float64, n)
	for i := 0; i < n; i++ {
		if M[i][i] == 0 {
			return nil, false
		}
		z[i] = M[i][n] / M[i][i]
	}
	return z, true
}

// leastSquares solves the ordinary least-squares problem min ||X·b - y|| via
// the normal equations XᵀX b = Xᵀy. X is n×p (n rows, p columns). It returns
// the coefficient vector of length p and true, or nil and false if the normal
// matrix is singular.
func leastSquares(X [][]float64, y []float64) ([]float64, bool) {
	n := len(X)
	if n == 0 || len(y) != n {
		return nil, false
	}
	p := len(X[0])
	if p == 0 {
		return nil, false
	}
	// Form XᵀX (p×p) and Xᵀy (p).
	ata := make([][]float64, p)
	for i := range ata {
		ata[i] = make([]float64, p)
	}
	aty := make([]float64, p)
	for r := 0; r < n; r++ {
		row := X[r]
		if len(row) != p {
			return nil, false
		}
		yr := y[r]
		for i := 0; i < p; i++ {
			aty[i] += row[i] * yr
			for j := 0; j < p; j++ {
				ata[i][j] += row[i] * row[j]
			}
		}
	}
	return solveLinear(ata, aty)
}

// levinson runs the Levinson–Durbin recursion on the autocovariance sequence
// gamma (gamma[0..p]) and returns the AR(p) coefficients phi[0..p-1] (so that
// x_t = phi_0 x_{t-1} + … + phi_{p-1} x_{t-p} + e_t) together with the white
// noise variance. If gamma[0] is zero it returns zero coefficients.
func levinson(gamma []float64, p int) ([]float64, float64) {
	phi := make([]float64, p)
	if p <= 0 {
		v := 0.0
		if len(gamma) > 0 {
			v = gamma[0]
		}
		return phi, v
	}
	v := gamma[0]
	if v == 0 {
		return phi, 0
	}
	prev := make([]float64, p)
	for k := 1; k <= p; k++ {
		acc := gamma[k]
		for j := 1; j < k; j++ {
			acc -= phi[j-1] * gamma[k-j]
		}
		refl := acc / v
		copy(prev, phi)
		phi[k-1] = refl
		for j := 1; j < k; j++ {
			phi[j-1] = prev[j-1] - refl*prev[k-1-j]
		}
		v *= 1 - refl*refl
		if v <= 0 {
			v = 0
			break
		}
	}
	return phi, v
}

// invertMatrix returns the inverse of the square matrix A using Gauss–Jordan
// elimination with partial pivoting, or nil and false if A is singular. A is
// not modified.
func invertMatrix(A [][]float64) ([][]float64, bool) {
	n := len(A)
	if n == 0 {
		return nil, false
	}
	// Augmented [A | I].
	M := make([][]float64, n)
	for i := range A {
		if len(A[i]) != n {
			return nil, false
		}
		M[i] = make([]float64, 2*n)
		copy(M[i], A[i])
		M[i][n+i] = 1
	}
	for col := 0; col < n; col++ {
		piv := col
		best := math.Abs(M[col][col])
		for r := col + 1; r < n; r++ {
			if a := math.Abs(M[r][col]); a > best {
				best = a
				piv = r
			}
		}
		if best == 0 {
			return nil, false
		}
		M[col], M[piv] = M[piv], M[col]
		pv := M[col][col]
		for c := 0; c < 2*n; c++ {
			M[col][c] /= pv
		}
		for r := 0; r < n; r++ {
			if r == col {
				continue
			}
			f := M[r][col]
			if f == 0 {
				continue
			}
			for c := 0; c < 2*n; c++ {
				M[r][c] -= f * M[col][c]
			}
		}
	}
	inv := make([][]float64, n)
	for i := 0; i < n; i++ {
		inv[i] = make([]float64, n)
		copy(inv[i], M[i][n:])
	}
	return inv, true
}

// olsStats fits y = X·β by ordinary least squares and returns the coefficient
// vector, the standard errors of the coefficients, the residual variance, and
// true on success. X is n×p with n > p.
func olsStats(X [][]float64, y []float64) (beta, se []float64, sigma2 float64, ok bool) {
	n := len(X)
	if n == 0 || len(y) != n {
		return nil, nil, 0, false
	}
	p := len(X[0])
	if p == 0 || n <= p {
		return nil, nil, 0, false
	}
	ata := make([][]float64, p)
	for i := range ata {
		ata[i] = make([]float64, p)
	}
	aty := make([]float64, p)
	for r := 0; r < n; r++ {
		row := X[r]
		if len(row) != p {
			return nil, nil, 0, false
		}
		for i := 0; i < p; i++ {
			aty[i] += row[i] * y[r]
			for j := 0; j < p; j++ {
				ata[i][j] += row[i] * row[j]
			}
		}
	}
	inv, okInv := invertMatrix(ata)
	if !okInv {
		return nil, nil, 0, false
	}
	beta = make([]float64, p)
	for i := 0; i < p; i++ {
		var s float64
		for j := 0; j < p; j++ {
			s += inv[i][j] * aty[j]
		}
		beta[i] = s
	}
	var sse float64
	for r := 0; r < n; r++ {
		pred := 0.0
		for j := 0; j < p; j++ {
			pred += X[r][j] * beta[j]
		}
		e := y[r] - pred
		sse += e * e
	}
	sigma2 = sse / float64(n-p)
	se = make([]float64, p)
	for i := 0; i < p; i++ {
		v := sigma2 * inv[i][i]
		if v < 0 {
			v = 0
		}
		se[i] = math.Sqrt(v)
	}
	return beta, se, sigma2, true
}

// approxEqual reports whether a and b are within tol.
func approxEqual(a, b, tol float64) bool {
	if math.IsNaN(a) && math.IsNaN(b) {
		return true
	}
	return math.Abs(a-b) <= tol
}
