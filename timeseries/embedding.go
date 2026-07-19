package timeseries

import "math"

// LagMatrix builds a design matrix of lagged values for autoregressive modeling.
// Row i (for i from 0) contains [x[t−1], x[t−2], …, x[t−p]] where t = p+i, so the
// matrix has len(x)−p rows and p columns. The aligned targets x[t] are returned
// as the second result. It returns nil, nil if p < 1 or p ≥ len(x).
func LagMatrix(x []float64, p int) ([][]float64, []float64) {
	n := len(x)
	if p < 1 || p >= n {
		return nil, nil
	}
	rows := n - p
	X := make([][]float64, rows)
	y := make([]float64, rows)
	for t := p; t < n; t++ {
		row := make([]float64, p)
		for j := 1; j <= p; j++ {
			row[j-1] = x[t-j]
		}
		X[t-p] = row
		y[t-p] = x[t]
	}
	return X, y
}

// Embed constructs a time-delay (Takens) embedding of the series with embedding
// dimension m and delay tau. Row i is [x[i], x[i+tau], …, x[i+(m−1)tau]]; the
// matrix has len(x)−(m−1)·tau rows and m columns. It returns nil if the
// arguments are invalid or the series is too short.
func Embed(x []float64, m, tau int) [][]float64 {
	if m < 1 || tau < 1 {
		return nil
	}
	span := (m - 1) * tau
	n := len(x)
	if span >= n {
		return nil
	}
	rows := n - span
	out := make([][]float64, rows)
	for i := 0; i < rows; i++ {
		row := make([]float64, m)
		for j := 0; j < m; j++ {
			row[j] = x[i+j*tau]
		}
		out[i] = row
	}
	return out
}

// TimeDelayEmbedding is an alias for [Embed] with the conventional argument
// order (dimension then delay), returning the delay-coordinate matrix.
func TimeDelayEmbedding(x []float64, dimension, delay int) [][]float64 {
	return Embed(x, dimension, delay)
}

// HankelMatrix builds the Hankel matrix of the series with the given number of
// rows: element (i,j) is x[i+j]. The matrix has len(x)−rows+1 columns. It
// returns nil if rows < 1 or rows > len(x).
func HankelMatrix(x []float64, rows int) [][]float64 {
	n := len(x)
	if rows < 1 || rows > n {
		return nil
	}
	cols := n - rows + 1
	out := make([][]float64, rows)
	for i := 0; i < rows; i++ {
		out[i] = make([]float64, cols)
		for j := 0; j < cols; j++ {
			out[i][j] = x[i+j]
		}
	}
	return out
}

// ToeplitzMatrix builds the symmetric Toeplitz matrix whose first row and
// column are c: element (i,j) is c[|i−j|]. It returns nil for an empty c.
func ToeplitzMatrix(c []float64) [][]float64 {
	n := len(c)
	if n == 0 {
		return nil
	}
	out := make([][]float64, n)
	for i := 0; i < n; i++ {
		out[i] = make([]float64, n)
		for j := 0; j < n; j++ {
			d := i - j
			if d < 0 {
				d = -d
			}
			out[i][j] = c[d]
		}
	}
	return out
}

// AutocorrelationMatrix returns the p×p symmetric Toeplitz autocorrelation
// matrix of the series, with element (i,j) equal to the sample autocorrelation
// at lag |i−j|. It is the coefficient matrix of the Yule–Walker equations.
func AutocorrelationMatrix(x []float64, p int) [][]float64 {
	if p < 1 {
		return nil
	}
	rho := AutoCorrelation(x, p-1)
	return ToeplitzMatrix(rho)
}

// SlidingWindows returns the list of consecutive full windows of length w over
// the series, each stepped by step positions. Each window is a fresh copy. It
// returns nil if w < 1, step < 1 or w > len(x).
func SlidingWindows(x []float64, w, step int) [][]float64 {
	n := len(x)
	if w < 1 || step < 1 || w > n {
		return nil
	}
	var out [][]float64
	for start := 0; start+w <= n; start += step {
		win := make([]float64, w)
		copy(win, x[start:start+w])
		out = append(out, win)
	}
	return out
}

// TakensThetaAutoMI returns the lag at which the sample autocorrelation of the
// series first drops to or below 1/e, a common heuristic for choosing the delay
// in a time-delay embedding. It searches lags 1..maxLag and returns 0 if the
// threshold is never crossed.
func TakensThetaAutoMI(x []float64, maxLag int) int {
	if maxLag < 1 || maxLag >= len(x) {
		if len(x) > 1 {
			maxLag = len(x) - 1
		} else {
			return 0
		}
	}
	acf := AutoCorrelation(x, maxLag)
	thr := 1 / math.E
	for k := 1; k <= maxLag; k++ {
		if acf[k] <= thr {
			return k
		}
	}
	return 0
}

// FlattenMatrix returns a row-major flattening of a matrix into a single slice,
// a convenience for feeding embedded matrices to routines that expect flat data.
func FlattenMatrix(m [][]float64) []float64 {
	var out []float64
	for _, row := range m {
		out = append(out, row...)
	}
	return out
}
