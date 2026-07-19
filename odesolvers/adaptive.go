package odesolvers

import (
	"math"
)

// AdaptiveOptions configures the embedded-Runge-Kutta step-size controller.
type AdaptiveOptions struct {
	// RelTol and AbsTol are the relative and absolute error tolerances used to
	// form the per-component scale atol + rtol*max(|y|,|ynext|).
	RelTol float64
	AbsTol float64
	// HInit is the initial step-size guess; when zero an automatic guess is
	// made. HMin and HMax bound the step magnitude (HMax == 0 means unbounded).
	HInit float64
	HMin  float64
	HMax  float64
	// MaxSteps limits the number of accepted-or-rejected steps.
	MaxSteps int
	// Safety, MinScale and MaxScale parameterize the step-update factor
	// h_new = h * clamp(Safety * err^(-1/(q+1)), MinScale, MaxScale).
	Safety   float64
	MinScale float64
	MaxScale float64
	// Dense, when true, records endpoint derivatives to enable Hermite
	// dense output in the returned Solution.
	Dense bool
}

// DefaultAdaptiveOptions returns a reasonable default configuration:
// RelTol=1e-6, AbsTol=1e-9, automatic initial step, up to 100000 steps.
func DefaultAdaptiveOptions() AdaptiveOptions {
	return AdaptiveOptions{
		RelTol:   1e-6,
		AbsTol:   1e-9,
		MaxSteps: 100000,
		Safety:   0.9,
		MinScale: 0.2,
		MaxScale: 5.0,
	}
}

// withDefaults fills unset fields of the options with sensible values.
func (o AdaptiveOptions) withDefaults() AdaptiveOptions {
	if o.RelTol <= 0 {
		o.RelTol = 1e-6
	}
	if o.AbsTol <= 0 {
		o.AbsTol = 1e-9
	}
	if o.MaxSteps <= 0 {
		o.MaxSteps = 100000
	}
	if o.Safety <= 0 {
		o.Safety = 0.9
	}
	if o.MinScale <= 0 {
		o.MinScale = 0.2
	}
	if o.MaxScale <= 0 {
		o.MaxScale = 5.0
	}
	return o
}

// initialStep produces an automatic first-step guess following the heuristic of
// Hairer, Norsett and Wanner.
func initialStep(f Field, t0 float64, y0 []float64, dir float64, order int, o AdaptiveOptions) float64 {
	f0 := f(t0, y0)
	n := len(y0)
	scale := make([]float64, n)
	for i := range scale {
		scale[i] = o.AbsTol + o.RelTol*math.Abs(y0[i])
	}
	d0 := WeightedRMSNorm(y0, scale)
	d1 := WeightedRMSNorm(f0, scale)
	var h0 float64
	if d0 < 1e-5 || d1 < 1e-5 {
		h0 = 1e-6
	} else {
		h0 = 0.01 * d0 / d1
	}
	y1 := AXPY(dir*h0, f0, y0)
	f1 := f(t0+dir*h0, y1)
	df := make([]float64, n)
	for i := range df {
		df[i] = (f1[i] - f0[i]) / scale[i]
	}
	d2 := Norm2(df) / math.Sqrt(float64(n)) / h0
	var h1 float64
	m := math.Max(d1, d2)
	if m <= 1e-15 {
		h1 = math.Max(1e-6, h0*1e-3)
	} else {
		h1 = math.Pow(0.01/m, 1.0/float64(order+1))
	}
	h := math.Min(100*h0, h1)
	if o.HMax > 0 && h > o.HMax {
		h = o.HMax
	}
	return dir * h
}

// IntegrateAdaptive integrates y' = f(t, y) from t0 to tEnd with the embedded
// tableau bt and automatic step-size control. It returns the trajectory and,
// on failure to reach tEnd, a non-nil error (one of [ErrMaxSteps] or
// [ErrStepTooSmall]) alongside the partial Solution.
func IntegrateAdaptive(f Field, bt *ButcherTableau, t0 float64, y0 []float64, tEnd float64, opts AdaptiveOptions) (*Solution, error) {
	if !bt.IsEmbedded() {
		return nil, ErrInvalidInput
	}
	o := opts.withDefaults()
	dir := 1.0
	if tEnd < t0 {
		dir = -1.0
	}
	// The estimator order is the lower of the two embedded orders.
	q := bt.Order
	if bt.OrderStar < q {
		q = bt.OrderStar
	}

	sol := &Solution{Method: bt.Name}
	y := Clone(y0)
	t := t0
	if o.Dense {
		sol.pushWithDeriv(t, y, f(t, y))
	} else {
		sol.push(t, y)
	}

	var h float64
	if o.HInit > 0 {
		h = dir * o.HInit
	} else {
		h = initialStep(f, t0, y0, dir, q, o)
	}
	if o.HMax > 0 && math.Abs(h) > o.HMax {
		h = dir * o.HMax
	}

	minStep := o.HMin
	if minStep <= 0 {
		// A floor a few machine epsilons above the resolution of the time
		// variable, so genuine step stalls are caught without rejecting the
		// (perfectly fine) small steps needed near a distant endpoint.
		scale := math.Max(1, math.Max(math.Abs(t0), math.Abs(tEnd)))
		minStep = 16 * 2.220446049250313e-16 * scale
	}

	steps := 0
	for (dir > 0 && t < tEnd) || (dir < 0 && t > tEnd) {
		if steps >= o.MaxSteps {
			return sol, ErrMaxSteps
		}
		steps++
		// Do not overshoot the endpoint.
		if dir > 0 && t+h > tEnd {
			h = tEnd - t
		}
		if dir < 0 && t+h < tEnd {
			h = tEnd - t
		}

		ynext, errv := bt.StepEmbedded(f, t, y, h)
		// Scaled error norm.
		n := len(y)
		scale := make([]float64, n)
		for i := 0; i < n; i++ {
			scale[i] = o.AbsTol + o.RelTol*math.Max(math.Abs(y[i]), math.Abs(ynext[i]))
		}
		errNorm := WeightedRMSNorm(errv, scale)

		if errNorm <= 1 {
			// Accept.
			t += h
			y = ynext
			if o.Dense {
				sol.pushWithDeriv(t, y, f(t, y))
			} else {
				sol.push(t, y)
			}
			sol.Accepted++
			// Snap exactly onto the endpoint to avoid tiny trailing steps.
			if (dir > 0 && t >= tEnd) || (dir < 0 && t <= tEnd) {
				break
			}
		} else {
			sol.Rejected++
		}

		// Step-size update.
		var factor float64
		if errNorm == 0 {
			factor = o.MaxScale
		} else {
			factor = o.Safety * math.Pow(errNorm, -1.0/float64(q+1))
			factor = math.Max(o.MinScale, math.Min(o.MaxScale, factor))
		}
		h *= factor
		if o.HMax > 0 && math.Abs(h) > o.HMax {
			h = dir * o.HMax
		}
		if math.Abs(h) < minStep {
			return sol, ErrStepTooSmall
		}
	}
	return sol, nil
}

// SolveHeunEuler integrates adaptively with the Heun-Euler 2(1) pair.
func SolveHeunEuler(f Field, t0 float64, y0 []float64, tEnd float64, opts AdaptiveOptions) (*Solution, error) {
	return IntegrateAdaptive(f, HeunEulerTableau(), t0, y0, tEnd, opts)
}

// SolveBogackiShampine integrates adaptively with the Bogacki-Shampine 3(2)
// pair (ode23).
func SolveBogackiShampine(f Field, t0 float64, y0 []float64, tEnd float64, opts AdaptiveOptions) (*Solution, error) {
	return IntegrateAdaptive(f, BogackiShampineTableau(), t0, y0, tEnd, opts)
}

// SolveRKF45 integrates adaptively with the Runge-Kutta-Fehlberg 4(5) pair.
func SolveRKF45(f Field, t0 float64, y0 []float64, tEnd float64, opts AdaptiveOptions) (*Solution, error) {
	return IntegrateAdaptive(f, FehlbergTableau(), t0, y0, tEnd, opts)
}

// SolveCashKarp integrates adaptively with the Cash-Karp 5(4) pair.
func SolveCashKarp(f Field, t0 float64, y0 []float64, tEnd float64, opts AdaptiveOptions) (*Solution, error) {
	return IntegrateAdaptive(f, CashKarpTableau(), t0, y0, tEnd, opts)
}

// SolveDOPRI5 integrates adaptively with the Dormand-Prince 5(4) pair (ode45).
func SolveDOPRI5(f Field, t0 float64, y0 []float64, tEnd float64, opts AdaptiveOptions) (*Solution, error) {
	return IntegrateAdaptive(f, DormandPrinceTableau(), t0, y0, tEnd, opts)
}

// SolveRKF78 integrates adaptively with the high-order Fehlberg 7(8) pair.
func SolveRKF78(f Field, t0 float64, y0 []float64, tEnd float64, opts AdaptiveOptions) (*Solution, error) {
	return IntegrateAdaptive(f, Fehlberg78Tableau(), t0, y0, tEnd, opts)
}
