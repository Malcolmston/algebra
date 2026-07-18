package interval

import "math"

// Bisect splits x at its midpoint and returns the closed left and right halves,
// which share the midpoint as a common endpoint. An empty interval yields two
// empty intervals.
func Bisect(x Interval) (left, right Interval) {
	if x.IsEmpty() {
		return Empty(), Empty()
	}
	m := x.Midpoint()
	if math.IsNaN(m) || math.IsInf(m, 0) {
		return x, Empty()
	}
	return Interval{Lo: x.Lo, Hi: m}, Interval{Lo: m, Hi: x.Hi}
}

// Newton runs the interval Newton iteration to enclose a root of f within the
// starting interval x. The function df must return an enclosure of the
// derivative of f over any input interval; both f and df must be inclusion
// isotonic interval extensions of the real function and its derivative.
//
// At each step the Newton operator N(X) = m - f(m)/f'(X), evaluated at the
// midpoint m of X, is intersected with X to discard root-free regions. When the
// derivative enclosure straddles 0 the step falls back to bisection so progress
// continues. Iteration stops when the enclosure narrows to within tol or after
// maxIter steps.
//
// The second return value is true only when existence and uniqueness of a root
// in the returned interval were rigorously verified, which happens when N(X) is
// contained in X while f'(X) excludes 0. A false result with a non-empty
// interval means the enclosure is a plausible but unverified bracket; an empty
// interval proves f has no root in the original x.
func Newton(f, df func(Interval) Interval, x Interval, tol float64, maxIter int) (Interval, bool) {
	verified := false
	for i := 0; i < maxIter; i++ {
		if x.IsEmpty() {
			return Empty(), false
		}
		m := x.Midpoint()
		if math.IsNaN(m) || math.IsInf(m, 0) {
			return x, verified
		}
		dfx := df(x)
		if dfx.IsEmpty() || dfx.Contains(0) {
			// Derivative may vanish: take a bisection step, keeping the half
			// whose function enclosure still admits a root.
			left, right := Bisect(x)
			lRoot := !f(left).IsEmpty() && f(left).Contains(0)
			rRoot := !f(right).IsEmpty() && f(right).Contains(0)
			switch {
			case lRoot:
				x = left
			case rRoot:
				x = right
			default:
				return Empty(), false
			}
			if x.Width() <= tol {
				return x, false
			}
			continue
		}
		n := Point(m).Sub(f(Point(m)).Div(dfx))
		if x.ContainsInterval(n) {
			verified = true
		}
		next := x.Intersect(n)
		if next.IsEmpty() {
			return Empty(), false
		}
		if next.Width() <= tol || next.Equal(x) {
			return next, verified
		}
		x = next
	}
	return x, verified
}
