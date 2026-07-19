package splines

import "errors"

// errSingular is returned by the internal solvers on a numerically singular
// system.
var errSingular = errors.New("splines: singular system")

// SolveTridiagonal solves the tridiagonal linear system whose sub-, main- and
// super-diagonals are a, b and c and whose right-hand side is d, using the
// Thomas algorithm. a[0] and c[n-1] are ignored. The four slices must share a
// common length n >= 1. The returned slice is the solution vector; the inputs
// are not modified.
func SolveTridiagonal(a, b, c, d []float64) ([]float64, error) {
	n := len(b)
	if n == 0 {
		return nil, ErrEmpty
	}
	if len(a) != n || len(c) != n || len(d) != n {
		return nil, ErrLenMismatch
	}
	cp := make([]float64, n)
	dp := make([]float64, n)
	if b[0] == 0 {
		return nil, errSingular
	}
	cp[0] = c[0] / b[0]
	dp[0] = d[0] / b[0]
	for i := 1; i < n; i++ {
		m := b[i] - a[i]*cp[i-1]
		if m == 0 {
			return nil, errSingular
		}
		cp[i] = c[i] / m
		dp[i] = (d[i] - a[i]*dp[i-1]) / m
	}
	x := make([]float64, n)
	x[n-1] = dp[n-1]
	for i := n - 2; i >= 0; i-- {
		x[i] = dp[i] - cp[i]*x[i+1]
	}
	return x, nil
}

// SolveCyclicTridiagonal solves a tridiagonal system that additionally has the
// corner entries alpha (top-right, coupling the last unknown into the first
// row) and beta (bottom-left, coupling the first unknown into the last row).
// Such "periodic" systems arise from periodic splines. It uses the
// Sherman-Morrison formula on top of [SolveTridiagonal]. All of a, b, c and d
// must have the same length n >= 2.
func SolveCyclicTridiagonal(a, b, c, d []float64, alpha, beta float64) ([]float64, error) {
	n := len(b)
	if n < 2 {
		return nil, ErrTooFewPoints
	}
	if len(a) != n || len(c) != n || len(d) != n {
		return nil, ErrLenMismatch
	}
	gamma := -b[0]
	if gamma == 0 {
		gamma = 1
	}
	bb := make([]float64, n)
	copy(bb, b)
	bb[0] = b[0] - gamma
	bb[n-1] = b[n-1] - alpha*beta/gamma

	y, err := SolveTridiagonal(a, bb, c, d)
	if err != nil {
		return nil, err
	}
	u := make([]float64, n)
	u[0] = gamma
	u[n-1] = alpha
	z, err := SolveTridiagonal(a, bb, c, u)
	if err != nil {
		return nil, err
	}
	// fact = (v . y) / (1 + v . z), with v = (1,0,...,0,beta/gamma).
	vy := y[0] + beta/gamma*y[n-1]
	vz := z[0] + beta/gamma*z[n-1]
	fact := vy / (1 + vz)
	x := make([]float64, n)
	for i := 0; i < n; i++ {
		x[i] = y[i] - fact*z[i]
	}
	return x, nil
}
