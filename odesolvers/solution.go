package odesolvers

import (
	"math"
	"sort"
)

// Solution stores the discrete trajectory produced by an integrator.
//
// T holds the (strictly increasing, or strictly decreasing for backward
// integration) sample times. Y[i] is the state vector at time T[i]. When Deriv
// is non-nil, Deriv[i] holds f(T[i], Y[i]); the presence of derivatives enables
// cubic Hermite dense output through [Solution.At].
type Solution struct {
	T     []float64
	Y     [][]float64
	Deriv [][]float64
	// Accepted and Rejected count the accepted and rejected steps taken by an
	// adaptive integrator. They are left at zero by the fixed-step methods.
	Accepted int
	Rejected int
	// Method names the integrator that produced the solution.
	Method string
}

// newSolution returns an empty Solution seeded with the initial condition.
func newSolution(method string, t0 float64, y0 []float64) *Solution {
	return &Solution{
		T:      []float64{t0},
		Y:      [][]float64{Clone(y0)},
		Method: method,
	}
}

// push appends a sample (t, y) to the solution.
func (s *Solution) push(t float64, y []float64) {
	s.T = append(s.T, t)
	s.Y = append(s.Y, Clone(y))
}

// pushWithDeriv appends a sample (t, y) with its derivative dy.
func (s *Solution) pushWithDeriv(t float64, y, dy []float64) {
	s.T = append(s.T, t)
	s.Y = append(s.Y, Clone(y))
	if s.Deriv == nil {
		s.Deriv = [][]float64{}
	}
	s.Deriv = append(s.Deriv, Clone(dy))
}

// Len returns the number of stored samples.
func (s *Solution) Len() int { return len(s.T) }

// Steps returns the number of integration steps, one fewer than the number of
// stored samples (or zero when empty).
func (s *Solution) Steps() int {
	if len(s.T) == 0 {
		return 0
	}
	return len(s.T) - 1
}

// Dim returns the dimension of the state vectors, or 0 when the solution is
// empty.
func (s *Solution) Dim() int {
	if len(s.Y) == 0 {
		return 0
	}
	return len(s.Y[0])
}

// Initial returns a copy of the first stored state.
func (s *Solution) Initial() []float64 { return Clone(s.Y[0]) }

// Final returns a copy of the last stored state.
func (s *Solution) Final() []float64 { return Clone(s.Y[len(s.Y)-1]) }

// FinalTime returns the last stored time.
func (s *Solution) FinalTime() float64 { return s.T[len(s.T)-1] }

// InitialTime returns the first stored time.
func (s *Solution) InitialTime() float64 { return s.T[0] }

// Component returns the time series of state component j across all samples.
func (s *Solution) Component(j int) []float64 {
	out := make([]float64, len(s.Y))
	for i := range s.Y {
		out[i] = s.Y[i][j]
	}
	return out
}

// Times returns a copy of the sample-time slice.
func (s *Solution) Times() []float64 { return Clone(s.T) }

// isForward reports whether the solution advances in increasing time.
func (s *Solution) isForward() bool {
	return len(s.T) < 2 || s.T[len(s.T)-1] >= s.T[0]
}

// locate returns the index i such that t lies in the bracket [T[i], T[i+1]]
// (respecting the integration direction), clamped to the valid range.
func (s *Solution) locate(t float64) int {
	n := len(s.T)
	if n < 2 {
		return 0
	}
	if s.isForward() {
		i := sort.Search(n, func(k int) bool { return s.T[k] > t }) - 1
		if i < 0 {
			i = 0
		}
		if i > n-2 {
			i = n - 2
		}
		return i
	}
	// Decreasing times: T is descending.
	i := sort.Search(n, func(k int) bool { return s.T[k] < t }) - 1
	if i < 0 {
		i = 0
	}
	if i > n-2 {
		i = n - 2
	}
	return i
}

// At returns the interpolated state at time t. When derivative information is
// available cubic Hermite interpolation is used; otherwise linear interpolation
// is used. Times outside the covered range are extrapolated from the nearest
// bracket.
func (s *Solution) At(t float64) []float64 {
	if len(s.T) == 0 {
		return nil
	}
	if len(s.T) == 1 {
		return Clone(s.Y[0])
	}
	i := s.locate(t)
	t0, t1 := s.T[i], s.T[i+1]
	y0, y1 := s.Y[i], s.Y[i+1]
	h := t1 - t0
	if h == 0 {
		return Clone(y0)
	}
	theta := (t - t0) / h
	if s.Deriv != nil && i < len(s.Deriv)-1 {
		return HermiteInterpolate(theta, h, y0, y1, s.Deriv[i], s.Deriv[i+1])
	}
	// Linear interpolation.
	out := make([]float64, len(y0))
	for k := range y0 {
		out[k] = y0[k] + theta*(y1[k]-y0[k])
	}
	return out
}

// Interpolate is a synonym for [Solution.At].
func (s *Solution) Interpolate(t float64) []float64 { return s.At(t) }

// Sample evaluates the solution at n evenly spaced times over its covered
// interval and returns the times together with the interpolated states.
func (s *Solution) Sample(n int) ([]float64, [][]float64) {
	ts := Linspace(s.InitialTime(), s.FinalTime(), n)
	ys := make([][]float64, len(ts))
	for i, t := range ts {
		ys[i] = s.At(t)
	}
	return ts, ys
}

// HermiteInterpolate evaluates the cubic Hermite interpolant on a single step.
//
// Given the endpoint states y0, y1 and the endpoint derivatives f0, f1 of a
// step of length h, and a normalized position theta in [0, 1], it returns the
// interpolated state y(t0 + theta*h). The interpolant matches both the values
// and the derivatives at the endpoints and is therefore third-order accurate.
func HermiteInterpolate(theta, h float64, y0, y1, f0, f1 []float64) []float64 {
	t2 := theta * theta
	t3 := t2 * theta
	h00 := 2*t3 - 3*t2 + 1
	h10 := t3 - 2*t2 + theta
	h01 := -2*t3 + 3*t2
	h11 := t3 - t2
	out := make([]float64, len(y0))
	for k := range y0 {
		out[k] = h00*y0[k] + h10*h*f0[k] + h01*y1[k] + h11*h*f1[k]
	}
	return out
}

// MaxComponentError returns the maximum absolute difference between the
// components of a and b, a convenience for comparing trajectories in tests.
func MaxComponentError(a, b []float64) float64 {
	mustSameLen(a, b)
	var m float64
	for i := range a {
		if d := math.Abs(a[i] - b[i]); d > m {
			m = d
		}
	}
	return m
}
