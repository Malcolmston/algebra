package stochastic

import (
	"math"
	"sort"
)

// Path represents a sample path of a continuous-time process, stored as equal
// length slices of Times and Values sampled on a (usually regular) grid. Times
// are assumed strictly increasing.
type Path struct {
	Times  []float64
	Values []float64
}

// NewPath constructs a Path from the given times and values. The slices are
// copied so later mutation of the arguments does not affect the path.
func NewPath(times, values []float64) Path {
	t := append([]float64(nil), times...)
	v := append([]float64(nil), values...)
	if len(t) > len(v) {
		t = t[:len(v)]
	} else if len(v) > len(t) {
		v = v[:len(t)]
	}
	return Path{Times: t, Values: v}
}

// Len returns the number of points in the path.
func (p Path) Len() int { return len(p.Values) }

// Empty reports whether the path has no points.
func (p Path) Empty() bool { return len(p.Values) == 0 }

// Initial returns the first value of the path, or 0 for an empty path.
func (p Path) Initial() float64 {
	if len(p.Values) == 0 {
		return 0
	}
	return p.Values[0]
}

// Final returns the last value of the path, or 0 for an empty path.
func (p Path) Final() float64 {
	if len(p.Values) == 0 {
		return 0
	}
	return p.Values[len(p.Values)-1]
}

// StartTime returns the first time of the path, or 0 for an empty path.
func (p Path) StartTime() float64 {
	if len(p.Times) == 0 {
		return 0
	}
	return p.Times[0]
}

// EndTime returns the last time of the path, or 0 for an empty path.
func (p Path) EndTime() float64 {
	if len(p.Times) == 0 {
		return 0
	}
	return p.Times[len(p.Times)-1]
}

// Max returns the maximum value along the path.
func (p Path) Max() float64 {
	if len(p.Values) == 0 {
		return math.NaN()
	}
	m := p.Values[0]
	for _, v := range p.Values[1:] {
		if v > m {
			m = v
		}
	}
	return m
}

// Min returns the minimum value along the path.
func (p Path) Min() float64 {
	if len(p.Values) == 0 {
		return math.NaN()
	}
	m := p.Values[0]
	for _, v := range p.Values[1:] {
		if v < m {
			m = v
		}
	}
	return m
}

// ArgMax returns the index of the first maximal value along the path.
func (p Path) ArgMax() int {
	if len(p.Values) == 0 {
		return -1
	}
	idx := 0
	for i, v := range p.Values {
		if v > p.Values[idx] {
			idx = i
		}
	}
	return idx
}

// ArgMin returns the index of the first minimal value along the path.
func (p Path) ArgMin() int {
	if len(p.Values) == 0 {
		return -1
	}
	idx := 0
	for i, v := range p.Values {
		if v < p.Values[idx] {
			idx = i
		}
	}
	return idx
}

// Range returns the difference between the maximum and minimum values.
func (p Path) Range() float64 { return p.Max() - p.Min() }

// Mean returns the arithmetic mean of the path values.
func (p Path) Mean() float64 {
	if len(p.Values) == 0 {
		return math.NaN()
	}
	s := 0.0
	for _, v := range p.Values {
		s += v
	}
	return s / float64(len(p.Values))
}

// Variance returns the sample variance (with Bessel's correction) of the path
// values.
func (p Path) Variance() float64 {
	n := len(p.Values)
	if n < 2 {
		return 0
	}
	m := p.Mean()
	s := 0.0
	for _, v := range p.Values {
		d := v - m
		s += d * d
	}
	return s / float64(n-1)
}

// StdDev returns the sample standard deviation of the path values.
func (p Path) StdDev() float64 { return math.Sqrt(p.Variance()) }

// Increments returns the successive differences of the path values. The result
// has length Len()-1.
func (p Path) Increments() []float64 {
	if len(p.Values) < 2 {
		return nil
	}
	d := make([]float64, len(p.Values)-1)
	for i := 1; i < len(p.Values); i++ {
		d[i-1] = p.Values[i] - p.Values[i-1]
	}
	return d
}

// QuadraticVariation returns the sum of squared increments, an estimate of the
// accumulated quadratic variation of the path.
func (p Path) QuadraticVariation() float64 {
	s := 0.0
	for i := 1; i < len(p.Values); i++ {
		d := p.Values[i] - p.Values[i-1]
		s += d * d
	}
	return s
}

// TotalVariation returns the sum of the absolute increments of the path.
func (p Path) TotalVariation() float64 {
	s := 0.0
	for i := 1; i < len(p.Values); i++ {
		s += math.Abs(p.Values[i] - p.Values[i-1])
	}
	return s
}

// RunningMax returns the running maximum of the path values.
func (p Path) RunningMax() []float64 {
	if len(p.Values) == 0 {
		return nil
	}
	out := make([]float64, len(p.Values))
	m := p.Values[0]
	for i, v := range p.Values {
		if v > m {
			m = v
		}
		out[i] = m
	}
	return out
}

// RunningMin returns the running minimum of the path values.
func (p Path) RunningMin() []float64 {
	if len(p.Values) == 0 {
		return nil
	}
	out := make([]float64, len(p.Values))
	m := p.Values[0]
	for i, v := range p.Values {
		if v < m {
			m = v
		}
		out[i] = m
	}
	return out
}

// RunningMean returns the running arithmetic mean of the path values.
func (p Path) RunningMean() []float64 {
	if len(p.Values) == 0 {
		return nil
	}
	out := make([]float64, len(p.Values))
	s := 0.0
	for i, v := range p.Values {
		s += v
		out[i] = s / float64(i+1)
	}
	return out
}

// MaxDrawdown returns the largest peak-to-trough drop of the path, defined as
// the maximum over t of running-max minus value. It is non-negative.
func (p Path) MaxDrawdown() float64 {
	if len(p.Values) == 0 {
		return 0
	}
	peak := p.Values[0]
	dd := 0.0
	for _, v := range p.Values {
		if v > peak {
			peak = v
		}
		if peak-v > dd {
			dd = peak - v
		}
	}
	return dd
}

// IndexToHit returns the first index at which the path reaches or crosses the
// given level, measured relative to the starting value. It returns ok=false if
// the level is never reached.
func (p Path) IndexToHit(level float64) (int, bool) {
	if len(p.Values) == 0 {
		return -1, false
	}
	up := level >= p.Values[0]
	for i, v := range p.Values {
		if up && v >= level {
			return i, true
		}
		if !up && v <= level {
			return i, true
		}
	}
	return -1, false
}

// TimeToHit returns the time of the first index at which the path reaches or
// crosses the given level. It returns ok=false if the level is never reached.
func (p Path) TimeToHit(level float64) (float64, bool) {
	i, ok := p.IndexToHit(level)
	if !ok {
		return math.NaN(), false
	}
	return p.Times[i], true
}

// CrossesLevel reports whether the path ever reaches or crosses the given
// level relative to its starting value.
func (p Path) CrossesLevel(level float64) bool {
	_, ok := p.IndexToHit(level)
	return ok
}

// OccupationFraction returns the fraction of grid points whose value lies in
// the closed interval [lo, hi].
func (p Path) OccupationFraction(lo, hi float64) float64 {
	if len(p.Values) == 0 {
		return 0
	}
	c := 0
	for _, v := range p.Values {
		if v >= lo && v <= hi {
			c++
		}
	}
	return float64(c) / float64(len(p.Values))
}

// At returns the value of the path at time t using linear interpolation between
// grid points. Times outside the range are clamped to the endpoints.
func (p Path) At(t float64) float64 {
	n := len(p.Times)
	if n == 0 {
		return math.NaN()
	}
	if t <= p.Times[0] {
		return p.Values[0]
	}
	if t >= p.Times[n-1] {
		return p.Values[n-1]
	}
	i := sort.SearchFloat64s(p.Times, t)
	if i < len(p.Times) && p.Times[i] == t {
		return p.Values[i]
	}
	// p.Times[i-1] < t < p.Times[i]
	t0, t1 := p.Times[i-1], p.Times[i]
	v0, v1 := p.Values[i-1], p.Values[i]
	w := (t - t0) / (t1 - t0)
	return v0 + w*(v1-v0)
}

// Copy returns a deep copy of the path.
func (p Path) Copy() Path {
	return Path{
		Times:  append([]float64(nil), p.Times...),
		Values: append([]float64(nil), p.Values...),
	}
}

// Scale returns a new path whose values are multiplied by c.
func (p Path) Scale(c float64) Path {
	out := p.Copy()
	for i := range out.Values {
		out.Values[i] *= c
	}
	return out
}

// Shift returns a new path whose values are shifted by the additive constant c.
func (p Path) Shift(c float64) Path {
	out := p.Copy()
	for i := range out.Values {
		out.Values[i] += c
	}
	return out
}

// Antithetic returns the path reflected about its starting value, that is the
// path with value v(0) - (v(t) - v(0)) = 2*v(0) - v(t). This is the antithetic
// variate used for variance reduction.
func (p Path) Antithetic() Path {
	out := p.Copy()
	if len(out.Values) == 0 {
		return out
	}
	v0 := out.Values[0]
	for i := range out.Values {
		out.Values[i] = 2*v0 - out.Values[i]
	}
	return out
}
