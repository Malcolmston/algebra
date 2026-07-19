package timeseries

import (
	"errors"
	"math"
)

// ARModel is a fitted autoregressive model of a given order. The mean-adjusted
// form is x_t − μ = Σ_{i=1}^{p} Phi[i−1]·(x_{t−i} − μ) + e_t, with Intercept =
// μ·(1 − ΣPhi) so that x_t = Intercept + Σ Phi[i]·x_{t−i} + e_t.
type ARModel struct {
	Order     int       // AR order p
	Phi       []float64 // AR coefficients, length p
	Mean      float64   // series mean μ
	Intercept float64   // constant term c
	Sigma2    float64   // white-noise (innovation) variance
}

// LevinsonDurbin solves the Yule–Walker equations for an AR(p) model from the
// autocovariance sequence gamma (which must have length at least p+1). It
// returns the p AR coefficients and the innovation variance.
func LevinsonDurbin(gamma []float64, p int) ([]float64, float64) {
	if p < 0 || len(gamma) < p+1 {
		return nil, math.NaN()
	}
	phi, v := levinson(gamma, p)
	return phi, v
}

// ReflectionCoefficients returns the partial autocorrelation (reflection /
// PARCOR) coefficients k_1,…,k_p produced by the Levinson–Durbin recursion on
// the series, one per order from 1 to p.
func ReflectionCoefficients(x []float64, p int) []float64 {
	if p < 1 {
		return nil
	}
	gamma := AutoCovariance(x, p)
	out := make([]float64, p)
	phi := make([]float64, p)
	prev := make([]float64, p)
	v := gamma[0]
	for k := 1; k <= p; k++ {
		acc := gamma[k]
		for j := 1; j < k; j++ {
			acc -= phi[j-1] * gamma[k-j]
		}
		var refl float64
		if v != 0 {
			refl = acc / v
		}
		copy(prev, phi)
		phi[k-1] = refl
		for j := 1; j < k; j++ {
			phi[j-1] = prev[j-1] - refl*prev[k-1-j]
		}
		v *= 1 - refl*refl
		out[k-1] = refl
	}
	return out
}

// YuleWalker fits an AR(p) model to the series by solving the Yule–Walker
// equations via the Levinson–Durbin recursion. It returns an error if p < 1 or
// the series is shorter than p+1.
func YuleWalker(x []float64, p int) (*ARModel, error) {
	if p < 1 {
		return nil, errors.New("timeseries: YuleWalker order must be >= 1")
	}
	if len(x) < p+1 {
		return nil, errors.New("timeseries: series too short for the requested AR order")
	}
	gamma := AutoCovariance(x, p)
	phi, v := levinson(gamma, p)
	mu := mean(x)
	var sumPhi float64
	for _, c := range phi {
		sumPhi += c
	}
	return &ARModel{
		Order:     p,
		Phi:       phi,
		Mean:      mu,
		Intercept: mu * (1 - sumPhi),
		Sigma2:    v,
	}, nil
}

// ARFitLeastSquares fits an AR(p) model by ordinary least squares regression of
// x_t on its p lagged values and an intercept. It returns an error if p < 1 or
// there are too few usable rows.
func ARFitLeastSquares(x []float64, p int) (*ARModel, error) {
	if p < 1 {
		return nil, errors.New("timeseries: AR order must be >= 1")
	}
	n := len(x)
	if n < p+2 {
		return nil, errors.New("timeseries: series too short for least-squares AR")
	}
	rows := n - p
	X := make([][]float64, rows)
	y := make([]float64, rows)
	for t := p; t < n; t++ {
		row := make([]float64, p+1)
		row[0] = 1 // intercept
		for j := 1; j <= p; j++ {
			row[j] = x[t-j]
		}
		X[t-p] = row
		y[t-p] = x[t]
	}
	beta, ok := leastSquares(X, y)
	if !ok {
		return nil, errors.New("timeseries: singular least-squares system")
	}
	phi := make([]float64, p)
	copy(phi, beta[1:])
	// Innovation variance from residuals.
	var sse float64
	for i := 0; i < rows; i++ {
		pred := beta[0]
		for j := 1; j <= p; j++ {
			pred += beta[j] * X[i][j]
		}
		e := y[i] - pred
		sse += e * e
	}
	sigma2 := sse / float64(rows)
	var sumPhi float64
	for _, c := range phi {
		sumPhi += c
	}
	mu := mean(x)
	return &ARModel{
		Order:     p,
		Phi:       phi,
		Mean:      mu,
		Intercept: beta[0],
		Sigma2:    sigma2,
	}, nil
}

// BurgAR fits an AR(p) model using Burg's method, which minimizes the sum of
// forward and backward prediction errors and is well conditioned for short
// series. It returns an error if p < 1 or the series is too short.
func BurgAR(x []float64, p int) (*ARModel, error) {
	if p < 1 {
		return nil, errors.New("timeseries: Burg order must be >= 1")
	}
	n := len(x)
	if n < p+1 {
		return nil, errors.New("timeseries: series too short for Burg AR")
	}
	// Work on the mean-removed series.
	mu := mean(x)
	f := make([]float64, n)
	b := make([]float64, n)
	for i, v := range x {
		f[i] = v - mu
		b[i] = v - mu
	}
	a := make([]float64, p+1)
	a[0] = 1
	var dk float64
	for i := 0; i < n; i++ {
		dk += 2 * f[i] * f[i]
	}
	dk -= f[0]*f[0] + b[n-1]*b[n-1]

	for m := 0; m < p; m++ {
		// Reflection coefficient.
		var num float64
		for i := m + 1; i < n; i++ {
			num += f[i] * b[i-1]
		}
		var k float64
		if dk != 0 {
			k = -2 * num / dk
		}
		// Update AR coefficients.
		prev := make([]float64, m+2)
		copy(prev, a[:m+2])
		for i := 1; i <= m+1; i++ {
			a[i] = prev[i] + k*prev[m+1-i]
		}
		// Update forward/backward errors.
		for i := n - 1; i >= m+1; i-- {
			fi := f[i]
			f[i] = fi + k*b[i-1]
			b[i] = b[i-1] + k*fi
		}
		dk = (1-k*k)*dk - f[m+1]*f[m+1] - b[n-1]*b[n-1]
	}
	phi := make([]float64, p)
	for i := 1; i <= p; i++ {
		phi[i-1] = -a[i]
	}
	// Innovation variance from residuals of the fitted model.
	var sse float64
	cnt := 0
	for t := p; t < n; t++ {
		pred := 0.0
		for j := 0; j < p; j++ {
			pred += phi[j] * (x[t-1-j] - mu)
		}
		e := (x[t] - mu) - pred
		sse += e * e
		cnt++
	}
	sigma2 := math.NaN()
	if cnt > 0 {
		sigma2 = sse / float64(cnt)
	}
	var sumPhi float64
	for _, c := range phi {
		sumPhi += c
	}
	return &ARModel{
		Order:     p,
		Phi:       phi,
		Mean:      mu,
		Intercept: mu * (1 - sumPhi),
		Sigma2:    sigma2,
	}, nil
}

// Predict returns the in-sample one-step-ahead predictions of the AR model for
// the series x. The first Order elements, which lack sufficient history, are
// set to NaN.
func (m *ARModel) Predict(x []float64) []float64 {
	n := len(x)
	out := make([]float64, n)
	for t := 0; t < n; t++ {
		if t < m.Order {
			out[t] = math.NaN()
			continue
		}
		pred := m.Intercept
		for j := 0; j < m.Order; j++ {
			pred += m.Phi[j] * x[t-1-j]
		}
		out[t] = pred
	}
	return out
}

// Residuals returns the in-sample one-step-ahead residuals x_t − x̂_t of the AR
// model, with NaN for the first Order positions.
func (m *ARModel) Residuals(x []float64) []float64 {
	pred := m.Predict(x)
	out := make([]float64, len(x))
	for t := range x {
		if math.IsNaN(pred[t]) {
			out[t] = math.NaN()
		} else {
			out[t] = x[t] - pred[t]
		}
	}
	return out
}

// Forecast returns the h-step-ahead recursive forecasts of the AR model given
// the observed history x, feeding predicted values back in as needed.
func (m *ARModel) Forecast(x []float64, h int) []float64 {
	if h < 1 {
		return nil
	}
	p := m.Order
	// Buffer of the most recent p values, extended with forecasts.
	hist := make([]float64, 0, p+h)
	start := len(x) - p
	if start < 0 {
		start = 0
	}
	hist = append(hist, x[start:]...)
	out := make([]float64, h)
	for k := 0; k < h; k++ {
		pred := m.Intercept
		for j := 0; j < p; j++ {
			idx := len(hist) - 1 - j
			if idx < 0 {
				pred += m.Phi[j] * m.Mean
			} else {
				pred += m.Phi[j] * hist[idx]
			}
		}
		out[k] = pred
		hist = append(hist, pred)
	}
	return out
}
