package odesolvers

import "math"

// EventFunc is a scalar event (switching) function g(t, y). An event occurs
// where g changes sign along the trajectory.
type EventFunc func(t float64, y []float64) float64

// Event records a detected zero crossing of an [EventFunc].
type Event struct {
	// T is the time of the crossing and Y the interpolated state there.
	T float64
	Y []float64
	// Direction is +1 when g crosses from negative to positive and -1 for the
	// opposite crossing.
	Direction int
	// Interval is the index i such that the crossing lies in [T[i], T[i+1]] of
	// the scanned solution.
	Interval int
}

// ScanEvents scans a computed Solution for zero crossings of the event function
// g. The direction argument filters crossings: +1 keeps only upward
// (negative-to-positive) crossings, -1 keeps only downward crossings and 0
// keeps both. Each detected crossing is refined within its bracket by the
// secant/bisection hybrid [RefineEvent] operating on the dense output of sol.
func ScanEvents(sol *Solution, g EventFunc, direction int) []Event {
	var events []Event
	if sol.Len() < 2 {
		return events
	}
	gPrev := g(sol.T[0], sol.Y[0])
	for i := 0; i < sol.Len()-1; i++ {
		t1 := sol.T[i+1]
		gNext := g(t1, sol.Y[i+1])
		if gPrev == 0 {
			// Exact zero at the left endpoint.
			dir := crossingDirection(gPrev, gNext)
			if dir != 0 && (direction == 0 || dir == direction) {
				events = append(events, Event{T: sol.T[i], Y: Clone(sol.Y[i]), Direction: dir, Interval: i})
			}
		} else if gPrev*gNext < 0 {
			dir := crossingDirection(gPrev, gNext)
			if direction == 0 || dir == direction {
				tc := RefineEvent(sol, g, sol.T[i], t1)
				events = append(events, Event{T: tc, Y: sol.At(tc), Direction: dir, Interval: i})
			}
		}
		gPrev = gNext
	}
	return events
}

// crossingDirection returns +1 for a negative-to-positive crossing between the
// values a and b, -1 for the opposite, and 0 when there is no sign change.
func crossingDirection(a, b float64) int {
	switch {
	case a < 0 && b >= 0:
		return 1
	case a > 0 && b <= 0:
		return -1
	case a == 0 && b > 0:
		return 1
	case a == 0 && b < 0:
		return -1
	default:
		return 0
	}
}

// RefineEvent locates the crossing time of g within the bracket [ta, tb], where
// g(ta) and g(tb) are assumed to have opposite signs. It combines bisection
// with the secant method on the solution's dense output and returns the refined
// crossing time.
func RefineEvent(sol *Solution, g EventFunc, ta, tb float64) float64 {
	ga := g(ta, sol.At(ta))
	gb := g(tb, sol.At(tb))
	if ga == 0 {
		return ta
	}
	if gb == 0 {
		return tb
	}
	for iter := 0; iter < 100; iter++ {
		// Secant estimate.
		var tm float64
		if gb != ga {
			tm = tb - gb*(tb-ta)/(gb-ga)
		} else {
			tm = 0.5 * (ta + tb)
		}
		// Fall back to bisection when the secant step leaves the bracket.
		lo, hi := math.Min(ta, tb), math.Max(ta, tb)
		if tm <= lo || tm >= hi {
			tm = 0.5 * (ta + tb)
		}
		gm := g(tm, sol.At(tm))
		if math.Abs(gm) < 1e-14 || math.Abs(tb-ta) < 1e-14*(1+math.Abs(tm)) {
			return tm
		}
		if (ga < 0) == (gm < 0) {
			ta, ga = tm, gm
		} else {
			tb, gb = tm, gm
		}
	}
	return 0.5 * (ta + tb)
}

// FirstEvent returns the earliest event of g along sol in the given direction,
// or ok == false when none occur. "Earliest" respects the integration
// direction of sol.
func FirstEvent(sol *Solution, g EventFunc, direction int) (Event, bool) {
	events := ScanEvents(sol, g, direction)
	if len(events) == 0 {
		return Event{}, false
	}
	return events[0], true
}

// CountEvents returns the number of zero crossings of g along sol in the given
// direction.
func CountEvents(sol *Solution, g EventFunc, direction int) int {
	return len(ScanEvents(sol, g, direction))
}

// IntegrateUntilEvent integrates y' = f(t, y) with the fixed-step RK4 method
// from t0 with step h and stops at the first zero crossing of g in the given
// direction, or at tMax. It returns the Solution truncated at the crossing and,
// when a crossing was found, the corresponding Event with found == true.
func IntegrateUntilEvent(f Field, g EventFunc, direction int, t0 float64, y0 []float64, tMax, h float64) (*Solution, Event, bool) {
	sol := &Solution{Method: "RK4+event"}
	n, step := stepCount(t0, tMax, h)
	y := Clone(y0)
	t := t0
	sol.pushWithDeriv(t, y, f(t, y))
	gPrev := g(t, y)
	for i := 0; i < n; i++ {
		ynext := RK4Step(f, t, y, step)
		tNext := t0 + float64(i+1)*step
		gNext := g(tNext, ynext)
		y = ynext
		t = tNext
		sol.pushWithDeriv(t, y, f(t, y))
		dir := crossingDirection(gPrev, gNext)
		if dir != 0 && (direction == 0 || dir == direction) {
			tc := RefineEvent(sol, g, sol.T[len(sol.T)-2], t)
			ev := Event{T: tc, Y: sol.At(tc), Direction: dir, Interval: len(sol.T) - 2}
			return sol, ev, true
		}
		gPrev = gNext
	}
	return sol, Event{}, false
}
