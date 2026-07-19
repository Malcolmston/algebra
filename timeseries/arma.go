package timeseries

import (
	"errors"
	"math"
)

// InnovationsAlgorithm runs the innovations algorithm on the autocovariance
// sequence gamma (length ≥ maxlag+1) and returns the innovation coefficients
// theta as a lower-triangular matrix (theta[n][j] for 1 ≤ j ≤ n) together with
// the one-step prediction error variances v[0..maxlag]. These coefficients are
// the basis of the moving-average innovations estimator.
func InnovationsAlgorithm(gamma []float64, maxlag int) ([][]float64, []float64) {
	if maxlag < 0 || len(gamma) < maxlag+1 {
		return nil, nil
	}
	theta := make([][]float64, maxlag+1)
	for i := range theta {
		theta[i] = make([]float64, maxlag+1)
	}
	v := make([]float64, maxlag+1)
	v[0] = gamma[0]
	for n := 1; n <= maxlag; n++ {
		for k := 0; k < n; k++ {
			acc := gamma[n-k]
			for j := 0; j < k; j++ {
				acc -= theta[k][k-j] * theta[n][n-j] * v[j]
			}
			if v[k] != 0 {
				theta[n][n-k] = acc / v[k]
			}
		}
		vn := gamma[0]
		for j := 0; j < n; j++ {
			vn -= theta[n][n-j] * theta[n][n-j] * v[j]
		}
		v[n] = vn
	}
	return theta, v
}

// MAModel is a fitted moving-average model x_t = μ + e_t + Σ Theta[i]·e_{t−i−1},
// with innovation variance Sigma2.
type MAModel struct {
	Order  int       // MA order q
	Theta  []float64 // MA coefficients, length q
	Mean   float64   // series mean μ
	Sigma2 float64   // innovation variance
}

// MAFit estimates an MA(q) model using the innovations algorithm applied to the
// sample autocovariances of the series. It returns an error if q < 1 or the
// series is too short.
func MAFit(x []float64, q int) (*MAModel, error) {
	if q < 1 {
		return nil, errors.New("timeseries: MA order must be >= 1")
	}
	if len(x) < q+1 {
		return nil, errors.New("timeseries: series too short for the requested MA order")
	}
	gamma := AutoCovariance(x, q)
	theta, v := InnovationsAlgorithm(gamma, q)
	coeffs := make([]float64, q)
	for i := 1; i <= q; i++ {
		coeffs[i-1] = theta[q][i]
	}
	return &MAModel{
		Order:  q,
		Theta:  coeffs,
		Mean:   mean(x),
		Sigma2: v[q],
	}, nil
}

// Residuals returns the estimated innovations of the MA model reconstructed
// from the series by inverting the moving-average recursion:
// e_t = (x_t − μ) − Σ Theta[i]·e_{t−i−1}, with pre-sample innovations taken as
// zero.
func (m *MAModel) Residuals(x []float64) []float64 {
	n := len(x)
	e := make([]float64, n)
	for t := 0; t < n; t++ {
		val := x[t] - m.Mean
		for i := 0; i < m.Order; i++ {
			if t-1-i >= 0 {
				val -= m.Theta[i] * e[t-1-i]
			}
		}
		e[t] = val
	}
	return e
}

// Forecast returns the h-step-ahead forecasts of the MA(q) model given the
// observed history x. Because an MA(q) process is uncorrelated beyond lag q,
// forecasts revert to the mean after q steps.
func (m *MAModel) Forecast(x []float64, h int) []float64 {
	if h < 1 {
		return nil
	}
	e := m.Residuals(x)
	n := len(e)
	out := make([]float64, h)
	for k := 1; k <= h; k++ {
		val := m.Mean
		for i := 0; i < m.Order; i++ {
			// e_{n+k-1-i} is known only when k-1-i < 0, i.e. it refers to a
			// realized innovation at or before time n-1.
			lag := k - 1 - i
			if lag < 0 {
				idx := n + lag
				if idx >= 0 {
					val += m.Theta[i] * e[idx]
				}
			}
		}
		out[k-1] = val
	}
	return out
}

// ARMAModel is a fitted autoregressive moving-average model
// x_t − μ = Σ Phi[i]·(x_{t−i−1} − μ) + e_t + Σ Theta[j]·e_{t−j−1}.
type ARMAModel struct {
	P      int       // AR order
	Q      int       // MA order
	Phi    []float64 // AR coefficients, length P
	Theta  []float64 // MA coefficients, length Q
	Mean   float64   // series mean μ
	Sigma2 float64   // innovation variance
}

// ARMAFit fits an ARMA(p,q) model with the two-stage Hannan–Rissanen
// procedure: a long autoregression estimates the innovations, then x_t is
// regressed on its own lags and the estimated innovation lags. It returns an
// error for invalid orders or an insufficiently long series.
func ARMAFit(x []float64, p, q int) (*ARMAModel, error) {
	if p < 0 || q < 0 || p+q < 1 {
		return nil, errors.New("timeseries: ARMA needs p+q >= 1 and non-negative orders")
	}
	n := len(x)
	mu := mean(x)
	w := make([]float64, n)
	for i, v := range x {
		w[i] = v - mu
	}
	// Stage 1: high-order AR to obtain residuals (innovation proxies).
	m := p + q + int(math.Ceil(math.Sqrt(float64(n))))
	if m < p+q+1 {
		m = p + q + 1
	}
	if m > n/2 {
		m = n / 2
	}
	if m < 1 {
		m = 1
	}
	if n < m+1 {
		return nil, errors.New("timeseries: series too short for ARMA fit")
	}
	arGamma := AutoCovariance(w, m)
	arPhi, _ := levinson(arGamma, m)
	resid := make([]float64, n)
	for t := 0; t < n; t++ {
		if t < m {
			resid[t] = 0
			continue
		}
		pred := 0.0
		for j := 0; j < m; j++ {
			pred += arPhi[j] * w[t-1-j]
		}
		resid[t] = w[t] - pred
	}
	// Stage 2: regress w_t on p own lags and q residual lags.
	start := m + max(p, q)
	if start >= n {
		start = m
	}
	rows := n - start
	if rows < p+q+1 {
		return nil, errors.New("timeseries: not enough rows for ARMA regression")
	}
	X := make([][]float64, rows)
	y := make([]float64, rows)
	for t := start; t < n; t++ {
		row := make([]float64, p+q)
		for j := 0; j < p; j++ {
			row[j] = w[t-1-j]
		}
		for j := 0; j < q; j++ {
			row[p+j] = resid[t-1-j]
		}
		X[t-start] = row
		y[t-start] = w[t]
	}
	beta, ok := leastSquares(X, y)
	if !ok {
		return nil, errors.New("timeseries: singular ARMA regression")
	}
	phi := make([]float64, p)
	theta := make([]float64, q)
	copy(phi, beta[:p])
	copy(theta, beta[p:])
	model := &ARMAModel{P: p, Q: q, Phi: phi, Theta: theta, Mean: mu}
	// Innovation variance from model residuals.
	model.Sigma2 = variance0(model.Residuals(x)[start:])
	return model, nil
}

// variance0 returns the mean of squares of the (assumed zero-mean) residuals.
func variance0(e []float64) float64 {
	if len(e) == 0 {
		return math.NaN()
	}
	var s float64
	for _, v := range e {
		s += v * v
	}
	return s / float64(len(e))
}

// Residuals returns the estimated innovations of the ARMA model reconstructed
// from the series by inverting the recursion, using zero pre-sample values.
func (m *ARMAModel) Residuals(x []float64) []float64 {
	n := len(x)
	e := make([]float64, n)
	for t := 0; t < n; t++ {
		wt := x[t] - m.Mean
		pred := 0.0
		for j := 0; j < m.P; j++ {
			if t-1-j >= 0 {
				pred += m.Phi[j] * (x[t-1-j] - m.Mean)
			}
		}
		for j := 0; j < m.Q; j++ {
			if t-1-j >= 0 {
				pred += m.Theta[j] * e[t-1-j]
			}
		}
		e[t] = wt - pred
	}
	return e
}

// Forecast returns the h-step-ahead forecasts of the ARMA model given the
// observed history x, feeding forecasts back for the AR part and dropping
// future (zero-expectation) innovations for the MA part.
func (m *ARMAModel) Forecast(x []float64, h int) []float64 {
	if h < 1 {
		return nil
	}
	e := m.Residuals(x)
	n := len(x)
	// Extend the working (deviation) series with forecasts.
	w := make([]float64, n+h)
	for i := 0; i < n; i++ {
		w[i] = x[i] - m.Mean
	}
	ee := make([]float64, n+h)
	copy(ee, e)
	out := make([]float64, h)
	for k := 0; k < h; k++ {
		t := n + k
		pred := 0.0
		for j := 0; j < m.P; j++ {
			if t-1-j >= 0 {
				pred += m.Phi[j] * w[t-1-j]
			}
		}
		for j := 0; j < m.Q; j++ {
			if t-1-j >= 0 && t-1-j < n {
				pred += m.Theta[j] * ee[t-1-j]
			}
		}
		w[t] = pred
		out[k] = pred + m.Mean
	}
	return out
}

// ARIMAModel is a fitted ARIMA(p,d,q) model: after d-fold differencing the
// series is modeled as ARMA(p,q). The stored tail retains the last d level
// values needed to integrate forecasts back to the original scale.
type ARIMAModel struct {
	P, D, Q int
	ARMA    *ARMAModel
	tail    []float64 // last D+... original values for integration
	lastVal []float64 // the final D anchor values per differencing level
}

// ARIMAFit fits an ARIMA(p,d,q) model by differencing the series d times and
// fitting an ARMA(p,q) model to the result. It returns an error for invalid
// orders or an insufficiently long series.
func ARIMAFit(x []float64, p, d, q int) (*ARIMAModel, error) {
	if d < 0 {
		return nil, errors.New("timeseries: ARIMA d must be >= 0")
	}
	if len(x) <= d {
		return nil, errors.New("timeseries: series too short for the requested differencing")
	}
	// Keep the last value at each differencing level for integration.
	anchors := make([]float64, d)
	cur := copyf(x)
	for i := 0; i < d; i++ {
		anchors[i] = cur[len(cur)-1]
		cur = Diff(cur)
	}
	arma, err := ARMAFit(cur, p, q)
	if err != nil {
		return nil, err
	}
	return &ARIMAModel{P: p, D: d, Q: q, ARMA: arma, lastVal: anchors, tail: copyf(cur)}, nil
}

// Forecast returns the h-step-ahead forecasts of the ARIMA model on the
// original scale, integrating the differenced ARMA forecasts back up through
// the d retained anchor values.
func (m *ARIMAModel) Forecast(h int) []float64 {
	if h < 1 {
		return nil
	}
	fc := m.ARMA.Forecast(m.tail, h)
	// Integrate d times. At each level, the anchor is the last observed value
	// of the (d-1-k)-times differenced series.
	for level := m.D - 1; level >= 0; level-- {
		anchor := m.lastVal[level]
		acc := anchor
		integrated := make([]float64, h)
		for i := 0; i < h; i++ {
			acc += fc[i]
			integrated[i] = acc
		}
		fc = integrated
	}
	return fc
}

// ARMAToMA returns the first n coefficients ψ_1,…,ψ_n of the MA(∞)
// representation x_t = Σ ψ_j e_{t−j} implied by the AR coefficients phi and MA
// coefficients theta (ψ_0 = 1 is omitted; out[j−1] = ψ_j).
func ARMAToMA(phi, theta []float64, n int) []float64 {
	if n < 1 {
		return nil
	}
	psi := make([]float64, n+1)
	psi[0] = 1
	for j := 1; j <= n; j++ {
		var v float64
		if j-1 < len(theta) {
			v = theta[j-1]
		}
		for i := 0; i < len(phi); i++ {
			if j-1-i >= 0 {
				v += phi[i] * psi[j-1-i]
			}
		}
		psi[j] = v
	}
	return psi[1:]
}

// ARMAToAR returns the first n coefficients π_1,…,π_n of the AR(∞)
// representation Σ π_j x_{t−j} = e_t implied by the AR coefficients phi and MA
// coefficients theta (out[j−1] = π_j).
func ARMAToAR(phi, theta []float64, n int) []float64 {
	if n < 1 {
		return nil
	}
	// Recurrence π_k = φ_k + θ_k − Σ_{j=1}^{k−1} θ_j·π_{k−j}, derived from
	// Π(B)·Θ(B) = Φ(B) with the package's sign convention
	// x_t = Σ φ_i x_{t−i} + Σ θ_j e_{t−j} + e_t and
	// x_t = Σ π_k x_{t−k} + e_t.
	pi := make([]float64, n+1) // pi[k] = π_k, pi[0] unused
	for k := 1; k <= n; k++ {
		var v float64
		if k-1 < len(phi) {
			v += phi[k-1]
		}
		if k-1 < len(theta) {
			v += theta[k-1]
		}
		for j := 1; j < k; j++ {
			if j-1 < len(theta) {
				v -= theta[j-1] * pi[k-j]
			}
		}
		pi[k] = v
	}
	return pi[1:]
}
