package reliability

import (
	"errors"
	"math"
	"sort"
)

// KMEstimate holds the result of a Kaplan–Meier product-limit survival
// estimate. The slices are parallel and ordered by ascending event time; only
// times at which at least one failure occurred are recorded.
type KMEstimate struct {
	Times    []float64 // distinct failure times
	AtRisk   []int     // number at risk just before each time
	Events   []int     // number of failures at each time
	Censored []int     // number censored in the interval ending at each time
	Survival []float64 // estimated survival S(t) at and after each time
	Variance []float64 // Greenwood variance of S(t)
	StdErr   []float64 // standard error of S(t)
}

// NAEstimate holds the result of a Nelson–Aalen cumulative-hazard estimate.
// The slices are parallel and ordered by ascending event time.
type NAEstimate struct {
	Times     []float64 // distinct failure times
	AtRisk    []int     // number at risk just before each time
	Events    []int     // number of failures at each time
	CumHazard []float64 // estimated cumulative hazard H(t)
	Variance  []float64 // variance of H(t)
}

// eventGrid summarizes right-censored data at each distinct failure time.
type eventGrid struct {
	times    []float64
	deaths   []int
	censored []int
	atRisk   []int
}

// buildEventGrid pairs times with event indicators (true = observed failure,
// false = right-censored), then tabulates the number of deaths, censorings and
// units at risk at each distinct failure time.
func buildEventGrid(times []float64, events []bool) (*eventGrid, error) {
	if len(times) != len(events) {
		return nil, errors.New("reliability: times and events must have equal length")
	}
	if len(times) == 0 {
		return nil, errors.New("reliability: at least one observation is required")
	}
	type obs struct {
		t     float64
		event bool
	}
	data := make([]obs, len(times))
	for i := range times {
		if times[i] < 0 || math.IsNaN(times[i]) {
			return nil, errors.New("reliability: observation times must be non-negative")
		}
		data[i] = obs{times[i], events[i]}
	}
	sort.Slice(data, func(i, j int) bool { return data[i].t < data[j].t })
	n := len(data)

	var g eventGrid
	i := 0
	remaining := n
	for i < n {
		t := data[i].t
		d, c := 0, 0
		j := i
		for j < n && data[j].t == t {
			if data[j].event {
				d++
			} else {
				c++
			}
			j++
		}
		if d > 0 {
			g.times = append(g.times, t)
			g.deaths = append(g.deaths, d)
			g.censored = append(g.censored, c)
			g.atRisk = append(g.atRisk, remaining)
		}
		remaining -= (j - i)
		i = j
	}
	if len(g.times) == 0 {
		return nil, errors.New("reliability: data contains no observed failures")
	}
	return &g, nil
}

// KaplanMeier computes the Kaplan–Meier product-limit estimator of the survival
// function from right-censored data. events[i] is true when times[i] is an
// observed failure and false when it is a right-censoring time. The Greenwood
// formula supplies the variance of each survival estimate.
func KaplanMeier(times []float64, events []bool) (*KMEstimate, error) {
	g, err := buildEventGrid(times, events)
	if err != nil {
		return nil, err
	}
	m := len(g.times)
	est := &KMEstimate{
		Times:    append([]float64(nil), g.times...),
		AtRisk:   append([]int(nil), g.atRisk...),
		Events:   append([]int(nil), g.deaths...),
		Censored: append([]int(nil), g.censored...),
		Survival: make([]float64, m),
		Variance: make([]float64, m),
		StdErr:   make([]float64, m),
	}
	s := 1.0
	greenwood := 0.0
	for i := 0; i < m; i++ {
		n := float64(g.atRisk[i])
		d := float64(g.deaths[i])
		s *= 1 - d/n
		if n-d > 0 {
			greenwood += d / (n * (n - d))
		} else {
			greenwood += 0
		}
		est.Survival[i] = s
		v := s * s * greenwood
		est.Variance[i] = v
		est.StdErr[i] = math.Sqrt(v)
	}
	return est, nil
}

// SurvivalAt returns the Kaplan–Meier survival estimate S(t). It is a
// right-continuous step function equal to 1 before the first failure time and
// holding its most recent value between failures.
func (e *KMEstimate) SurvivalAt(t float64) float64 {
	s := 1.0
	for i, ti := range e.Times {
		if ti <= t {
			s = e.Survival[i]
		} else {
			break
		}
	}
	return s
}

// CumulativeHazardAt returns the Kaplan–Meier-derived cumulative hazard
// -ln S(t).
func (e *KMEstimate) CumulativeHazardAt(t float64) float64 {
	s := e.SurvivalAt(t)
	if s <= 0 {
		return math.Inf(1)
	}
	return -math.Log(s)
}

// MedianSurvival returns the smallest failure time at which the estimated
// survival drops to 0.5 or below. If survival never reaches 0.5 (heavy
// censoring) it returns NaN, signalling that the median is not reached.
func (e *KMEstimate) MedianSurvival() float64 {
	return e.Quantile(0.5)
}

// Quantile returns the smallest failure time at which the estimated survival
// drops to 1-p or below, i.e. the p-quantile of the lifetime distribution. It
// returns NaN when the quantile is not reached.
func (e *KMEstimate) Quantile(p float64) float64 {
	if p <= 0 || p >= 1 {
		return math.NaN()
	}
	threshold := 1 - p
	for i, s := range e.Survival {
		if s <= threshold {
			return e.Times[i]
		}
	}
	return math.NaN()
}

// RestrictedMeanSurvival returns the restricted mean survival time over [0,tau],
// the area under the Kaplan–Meier curve up to the horizon tau. It equals the
// expected lifetime truncated at tau.
func (e *KMEstimate) RestrictedMeanSurvival(tau float64) float64 {
	if tau <= 0 {
		return 0
	}
	area := 0.0
	prev := 0.0
	cur := 1.0
	for i, ti := range e.Times {
		upper := ti
		if upper > tau {
			upper = tau
		}
		if upper > prev {
			area += cur * (upper - prev)
			prev = upper
		}
		if ti >= tau {
			return area
		}
		cur = e.Survival[i]
	}
	if prev < tau {
		area += cur * (tau - prev)
	}
	return area
}

// NelsonAalen computes the Nelson–Aalen estimator of the cumulative hazard
// function from right-censored data, together with its variance. events[i] is
// true for an observed failure and false for a right-censoring time.
func NelsonAalen(times []float64, events []bool) (*NAEstimate, error) {
	g, err := buildEventGrid(times, events)
	if err != nil {
		return nil, err
	}
	m := len(g.times)
	est := &NAEstimate{
		Times:     append([]float64(nil), g.times...),
		AtRisk:    append([]int(nil), g.atRisk...),
		Events:    append([]int(nil), g.deaths...),
		CumHazard: make([]float64, m),
		Variance:  make([]float64, m),
	}
	h := 0.0
	v := 0.0
	for i := 0; i < m; i++ {
		n := float64(g.atRisk[i])
		d := float64(g.deaths[i])
		h += d / n
		v += d / (n * n)
		est.CumHazard[i] = h
		est.Variance[i] = v
	}
	return est, nil
}

// CumulativeHazardAt returns the Nelson–Aalen cumulative hazard estimate H(t),
// a right-continuous step function equal to 0 before the first failure.
func (e *NAEstimate) CumulativeHazardAt(t float64) float64 {
	h := 0.0
	for i, ti := range e.Times {
		if ti <= t {
			h = e.CumHazard[i]
		} else {
			break
		}
	}
	return h
}

// SurvivalAt returns the survival estimate exp(-H(t)) implied by the
// Nelson–Aalen cumulative hazard.
func (e *NAEstimate) SurvivalAt(t float64) float64 {
	return math.Exp(-e.CumulativeHazardAt(t))
}

// LogRankResult holds the outcome of a two-sample log-rank test comparing the
// survival of group 1 against group 2.
type LogRankResult struct {
	Observed1   float64 // observed failures in group 1
	Expected1   float64 // expected failures in group 1 under H0
	Variance    float64 // variance of (Observed1-Expected1)
	ChiSquare   float64 // log-rank chi-square statistic (1 d.f.)
	PValue      float64 // two-sided p-value
	HazardRatio float64 // one-step hazard-ratio estimate, group 1 vs group 2
}

// LogRankTest performs the two-sample log-rank test on right-censored survival
// data. Each group is supplied as parallel time/event slices where an event
// value of true marks an observed failure. The returned chi-square statistic
// has one degree of freedom.
func LogRankTest(times1 []float64, events1 []bool, times2 []float64, events2 []bool) (*LogRankResult, error) {
	if len(times1) != len(events1) || len(times2) != len(events2) {
		return nil, errors.New("reliability: times and events must have equal length in each group")
	}
	if len(times1) == 0 || len(times2) == 0 {
		return nil, errors.New("reliability: both groups must be non-empty")
	}
	type obs struct {
		t     float64
		event bool
		grp   int
	}
	var all []obs
	for i := range times1 {
		if times1[i] < 0 {
			return nil, errors.New("reliability: times must be non-negative")
		}
		all = append(all, obs{times1[i], events1[i], 1})
	}
	for i := range times2 {
		if times2[i] < 0 {
			return nil, errors.New("reliability: times must be non-negative")
		}
		all = append(all, obs{times2[i], events2[i], 2})
	}
	sort.Slice(all, func(i, j int) bool { return all[i].t < all[j].t })

	n := len(all)
	n1 := len(times1)
	nTotal := n
	// at risk in group 1 and overall are decremented as we advance through time
	atRisk1 := n1
	atRiskTot := nTotal

	var o1, e1, vSum float64
	i := 0
	for i < n {
		t := all[i].t
		d1, dTot, leave1, leaveTot := 0, 0, 0, 0
		j := i
		for j < n && all[j].t == t {
			leaveTot++
			if all[j].grp == 1 {
				leave1++
			}
			if all[j].event {
				dTot++
				if all[j].grp == 1 {
					d1++
				}
			}
			j++
		}
		if dTot > 0 {
			nt := float64(atRiskTot)
			n1t := float64(atRisk1)
			dt := float64(dTot)
			o1 += float64(d1)
			e1 += dt * n1t / nt
			if nt > 1 {
				vSum += dt * (nt - dt) * n1t * (nt - n1t) / (nt * nt * (nt - 1))
			}
		}
		atRisk1 -= leave1
		atRiskTot -= leaveTot
		i = j
	}

	res := &LogRankResult{Observed1: o1, Expected1: e1, Variance: vSum}
	if vSum > 0 {
		res.ChiSquare = (o1 - e1) * (o1 - e1) / vSum
		res.PValue = math.Erfc(math.Sqrt(res.ChiSquare / 2))
		res.HazardRatio = math.Exp((o1 - e1) / vSum)
	} else {
		res.ChiSquare = 0
		res.PValue = 1
		res.HazardRatio = math.NaN()
	}
	return res, nil
}

// HazardRatio returns the ratio of two hazard rates h1/h2, the constant
// proportional-hazards multiplier of group 1 relative to group 2.
func HazardRatio(h1, h2 float64) float64 {
	if h2 <= 0 || h1 < 0 {
		return math.NaN()
	}
	return h1 / h2
}

// ProportionalHazardsReliability returns the reliability of an item whose
// hazard is a constant multiple hr of a baseline hazard: R(t)=R0(t)^{hr},
// where r0 is the baseline reliability at t. This is the survival form of the
// Cox proportional-hazards model.
func ProportionalHazardsReliability(r0, hr float64) float64 {
	if r0 < 0 || r0 > 1 || hr < 0 {
		return math.NaN()
	}
	return math.Pow(r0, hr)
}

// ChiSquarePValue1DF returns the upper-tail probability P(X>=x) for a
// chi-square random variable with one degree of freedom, useful for
// interpreting a log-rank statistic.
func ChiSquarePValue1DF(x float64) float64 {
	if x < 0 {
		return math.NaN()
	}
	return math.Erfc(math.Sqrt(x / 2))
}
