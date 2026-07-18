package interval

import (
	"math"
	"strconv"
)

// Interval is a closed real interval [Lo, Hi] interpreted as the set of real
// numbers x with Lo <= x <= Hi. The zero value [Interval]{} is the degenerate
// interval [0, 0] containing only zero. Use [Empty] for the empty set and
// [Entire] for the whole real line.
type Interval struct {
	// Lo is the (inclusive) lower bound of the interval.
	Lo float64
	// Hi is the (inclusive) upper bound of the interval.
	Hi float64
}

// intervalPred returns the largest representable float64 strictly less than x,
// i.e. x rounded one ULP toward negative infinity. It is used to round lower
// bounds outward.
func intervalPred(x float64) float64 { return math.Nextafter(x, math.Inf(-1)) }

// intervalSucc returns the smallest representable float64 strictly greater than
// x, i.e. x rounded one ULP toward positive infinity. It is used to round upper
// bounds outward.
func intervalSucc(x float64) float64 { return math.Nextafter(x, math.Inf(1)) }

// intervalDownN rounds x outward toward -Inf by n ULPs.
func intervalDownN(x float64, n int) float64 {
	for i := 0; i < n; i++ {
		x = intervalPred(x)
	}
	return x
}

// intervalUpN rounds x outward toward +Inf by n ULPs.
func intervalUpN(x float64, n int) float64 {
	for i := 0; i < n; i++ {
		x = intervalSucc(x)
	}
	return x
}

// New returns the interval [lo, hi]. If lo > hi (and neither is NaN) the bounds
// are swapped so the result is always a well-formed interval. If either bound
// is NaN the result is [Empty].
func New(lo, hi float64) Interval {
	if math.IsNaN(lo) || math.IsNaN(hi) {
		return Empty()
	}
	if lo > hi {
		lo, hi = hi, lo
	}
	return Interval{Lo: lo, Hi: hi}
}

// Point returns the degenerate interval [x, x] containing only the single value
// x. If x is NaN the result is [Empty].
func Point(x float64) Interval {
	if math.IsNaN(x) {
		return Empty()
	}
	return Interval{Lo: x, Hi: x}
}

// Empty returns the empty interval, which contains no real numbers. It is
// encoded with Lo = +Inf and Hi = -Inf so that any Lo > Hi test recognises it.
func Empty() Interval { return Interval{Lo: math.Inf(1), Hi: math.Inf(-1)} }

// Entire returns the interval [-Inf, +Inf] representing the whole real line.
func Entire() Interval { return Interval{Lo: math.Inf(-1), Hi: math.Inf(1)} }

// IsEmpty reports whether the interval is empty, i.e. contains no real numbers.
func (a Interval) IsEmpty() bool {
	return math.IsNaN(a.Lo) || math.IsNaN(a.Hi) || a.Lo > a.Hi
}

// IsPoint reports whether the interval is degenerate, containing exactly one
// finite value (Lo == Hi and finite).
func (a Interval) IsPoint() bool {
	return !a.IsEmpty() && a.Lo == a.Hi && !math.IsInf(a.Lo, 0)
}

// IsEntire reports whether the interval is the whole real line [-Inf, +Inf].
func (a Interval) IsEntire() bool {
	return math.IsInf(a.Lo, -1) && math.IsInf(a.Hi, 1)
}

// IsBounded reports whether both bounds of the interval are finite. The empty
// interval is considered bounded.
func (a Interval) IsBounded() bool {
	if a.IsEmpty() {
		return true
	}
	return !math.IsInf(a.Lo, 0) && !math.IsInf(a.Hi, 0)
}

// Width returns Hi - Lo, the length of the interval, rounded up so the result
// never underestimates the true width. An empty interval has width 0; an
// unbounded interval has width +Inf.
func (a Interval) Width() float64 {
	if a.IsEmpty() {
		return 0
	}
	return intervalSucc(a.Hi - a.Lo)
}

// Radius returns half the width of the interval, rounded up. An empty interval
// has radius 0.
func (a Interval) Radius() float64 {
	if a.IsEmpty() {
		return 0
	}
	return intervalSucc((a.Hi - a.Lo) / 2)
}

// Midpoint returns a representable point near the centre (Lo + Hi) / 2. For an
// empty interval it returns NaN; for [Entire] it returns 0; for a half-bounded
// interval it returns the appropriate signed infinity.
func (a Interval) Midpoint() float64 {
	switch {
	case a.IsEmpty():
		return math.NaN()
	case a.IsEntire():
		return 0
	case math.IsInf(a.Lo, -1):
		return math.Inf(-1)
	case math.IsInf(a.Hi, 1):
		return math.Inf(1)
	default:
		return a.Lo/2 + a.Hi/2
	}
}

// Mag returns the magnitude of the interval, the supremum of |x| over all x in
// the interval (also called the norm). An empty interval returns NaN.
func (a Interval) Mag() float64 {
	if a.IsEmpty() {
		return math.NaN()
	}
	return math.Max(math.Abs(a.Lo), math.Abs(a.Hi))
}

// Mig returns the mignitude of the interval, the infimum of |x| over all x in
// the interval. It is 0 when the interval contains 0. An empty interval returns
// NaN.
func (a Interval) Mig() float64 {
	if a.IsEmpty() {
		return math.NaN()
	}
	if a.Lo <= 0 && a.Hi >= 0 {
		return 0
	}
	return math.Min(math.Abs(a.Lo), math.Abs(a.Hi))
}

// Contains reports whether the real number x lies within the interval. NaN is
// never contained.
func (a Interval) Contains(x float64) bool {
	if a.IsEmpty() || math.IsNaN(x) {
		return false
	}
	return a.Lo <= x && x <= a.Hi
}

// ContainsInterval reports whether b is a subset of a (every point of b lies in
// a). The empty interval is a subset of every interval.
func (a Interval) ContainsInterval(b Interval) bool {
	if b.IsEmpty() {
		return true
	}
	if a.IsEmpty() {
		return false
	}
	return a.Lo <= b.Lo && b.Hi <= a.Hi
}

// Overlaps reports whether a and b share at least one common point.
func (a Interval) Overlaps(b Interval) bool {
	if a.IsEmpty() || b.IsEmpty() {
		return false
	}
	return a.Lo <= b.Hi && b.Lo <= a.Hi
}

// Equal reports whether a and b denote exactly the same set. Two empty
// intervals are equal regardless of their encoded bounds.
func (a Interval) Equal(b Interval) bool {
	if a.IsEmpty() || b.IsEmpty() {
		return a.IsEmpty() && b.IsEmpty()
	}
	return a.Lo == b.Lo && a.Hi == b.Hi
}

// Intersect returns the intersection of a and b, the largest interval contained
// in both. If they do not overlap the result is [Empty].
func (a Interval) Intersect(b Interval) Interval {
	if a.IsEmpty() || b.IsEmpty() {
		return Empty()
	}
	lo := math.Max(a.Lo, b.Lo)
	hi := math.Min(a.Hi, b.Hi)
	if lo > hi {
		return Empty()
	}
	return Interval{Lo: lo, Hi: hi}
}

// Hull returns the convex hull (union enclosure) of a and b, the smallest
// interval containing both. The hull with an empty interval is the other
// interval.
func (a Interval) Hull(b Interval) Interval {
	if a.IsEmpty() {
		return b
	}
	if b.IsEmpty() {
		return a
	}
	return Interval{Lo: math.Min(a.Lo, b.Lo), Hi: math.Max(a.Hi, b.Hi)}
}

// String formats the interval as "[Lo, Hi]" using the shortest decimal that
// round-trips each bound. The empty interval formats as "[empty]".
func (a Interval) String() string {
	if a.IsEmpty() {
		return "[empty]"
	}
	lo := strconv.FormatFloat(a.Lo, 'g', -1, 64)
	hi := strconv.FormatFloat(a.Hi, 'g', -1, 64)
	return "[" + lo + ", " + hi + "]"
}
