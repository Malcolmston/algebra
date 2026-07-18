package interval

import "math"

// Add returns an enclosure of a + b, with the lower bound rounded down and the
// upper bound rounded up so the result rigorously contains every sum x + y for
// x in a and y in b.
func (a Interval) Add(b Interval) Interval {
	if a.IsEmpty() || b.IsEmpty() {
		return Empty()
	}
	return Interval{Lo: intervalPred(a.Lo + b.Lo), Hi: intervalSucc(a.Hi + b.Hi)}
}

// Sub returns an enclosure of a - b, outward rounded.
func (a Interval) Sub(b Interval) Interval {
	if a.IsEmpty() || b.IsEmpty() {
		return Empty()
	}
	return Interval{Lo: intervalPred(a.Lo - b.Hi), Hi: intervalSucc(a.Hi - b.Lo)}
}

// Neg returns -a, the reflection of the interval through the origin. It is
// exact and requires no rounding.
func (a Interval) Neg() Interval {
	if a.IsEmpty() {
		return Empty()
	}
	return Interval{Lo: -a.Hi, Hi: -a.Lo}
}

// Mul returns an enclosure of a * b, outward rounded. It evaluates the four
// endpoint products and takes their extremes, then inflates by one ULP to
// cover the rounding of the selected product.
func (a Interval) Mul(b Interval) Interval {
	if a.IsEmpty() || b.IsEmpty() {
		return Empty()
	}
	p1 := a.Lo * b.Lo
	p2 := a.Lo * b.Hi
	p3 := a.Hi * b.Lo
	p4 := a.Hi * b.Hi
	lo := math.Min(math.Min(p1, p2), math.Min(p3, p4))
	hi := math.Max(math.Max(p1, p2), math.Max(p3, p4))
	// Products involving 0 * Inf yield NaN; sanitise by falling back to the
	// widest safe bound.
	if math.IsNaN(lo) {
		lo = math.Inf(-1)
	}
	if math.IsNaN(hi) {
		hi = math.Inf(1)
	}
	return Interval{Lo: intervalPred(lo), Hi: intervalSucc(hi)}
}

// Recip returns an enclosure of 1 / a. If a contains 0 in its interior the
// result is [Entire]; if a is exactly [0, 0] the result is [Empty].
func (a Interval) Recip() Interval {
	if a.IsEmpty() {
		return Empty()
	}
	switch {
	case a.Lo == 0 && a.Hi == 0:
		return Empty()
	case a.Lo > 0 || a.Hi < 0:
		// 0 not in interval: reciprocal is monotone decreasing.
		return Interval{Lo: intervalPred(1 / a.Hi), Hi: intervalSucc(1 / a.Lo)}
	case a.Hi == 0:
		// a = [lo, 0], lo < 0: 1/a = (-Inf, 1/lo].
		return Interval{Lo: math.Inf(-1), Hi: intervalSucc(1 / a.Lo)}
	case a.Lo == 0:
		// a = [0, hi], hi > 0: 1/a = [1/hi, +Inf).
		return Interval{Lo: intervalPred(1 / a.Hi), Hi: math.Inf(1)}
	default:
		// 0 strictly inside: reciprocal is (-Inf, +Inf).
		return Entire()
	}
}

// Div returns an enclosure of a / b, outward rounded. Division by an interval
// containing 0 follows the same conventions as [Interval.Recip].
func (a Interval) Div(b Interval) Interval {
	if a.IsEmpty() || b.IsEmpty() {
		return Empty()
	}
	// When b does not contain 0 the direct endpoint formula is tighter than
	// multiplying by the reciprocal.
	if b.Lo > 0 || b.Hi < 0 {
		q1 := a.Lo / b.Lo
		q2 := a.Lo / b.Hi
		q3 := a.Hi / b.Lo
		q4 := a.Hi / b.Hi
		lo := math.Min(math.Min(q1, q2), math.Min(q3, q4))
		hi := math.Max(math.Max(q1, q2), math.Max(q3, q4))
		return Interval{Lo: intervalPred(lo), Hi: intervalSucc(hi)}
	}
	return a.Mul(b.Recip())
}

// AddFloat returns an enclosure of a + f, outward rounded.
func (a Interval) AddFloat(f float64) Interval { return a.Add(Point(f)) }

// SubFloat returns an enclosure of a - f, outward rounded.
func (a Interval) SubFloat(f float64) Interval { return a.Sub(Point(f)) }

// MulFloat returns an enclosure of a * f, outward rounded.
func (a Interval) MulFloat(f float64) Interval { return a.Mul(Point(f)) }

// DivFloat returns an enclosure of a / f, outward rounded. Division by 0
// follows [Interval.Div].
func (a Interval) DivFloat(f float64) Interval { return a.Div(Point(f)) }

// Abs returns the interval of absolute values {|x| : x in a}. If a straddles 0
// the lower bound is 0.
func (a Interval) Abs() Interval {
	if a.IsEmpty() {
		return Empty()
	}
	return Interval{Lo: a.Mig(), Hi: a.Mag()}
}

// Min returns the elementwise minimum enclosure {min(x, y) : x in a, y in b}.
func (a Interval) Min(b Interval) Interval {
	if a.IsEmpty() || b.IsEmpty() {
		return Empty()
	}
	return Interval{Lo: math.Min(a.Lo, b.Lo), Hi: math.Min(a.Hi, b.Hi)}
}

// Max returns the elementwise maximum enclosure {max(x, y) : x in a, y in b}.
func (a Interval) Max(b Interval) Interval {
	if a.IsEmpty() || b.IsEmpty() {
		return Empty()
	}
	return Interval{Lo: math.Max(a.Lo, b.Lo), Hi: math.Max(a.Hi, b.Hi)}
}

// Square returns an enclosure of a*a. It is tighter than [Interval.Mul] with a
// itself because it exploits x^2 >= 0: the lower bound is the squared
// mignitude and the upper bound the squared magnitude.
func (a Interval) Square() Interval {
	if a.IsEmpty() {
		return Empty()
	}
	mig := a.Mig()
	mag := a.Mag()
	lo := intervalPred(mig * mig)
	if lo < 0 {
		lo = 0 // squares are non-negative
	}
	return Interval{Lo: lo, Hi: intervalSucc(mag * mag)}
}

// Add returns an enclosure of a + b. It is the free-function form of
// [Interval.Add].
func Add(a, b Interval) Interval { return a.Add(b) }

// Sub returns an enclosure of a - b. It is the free-function form of
// [Interval.Sub].
func Sub(a, b Interval) Interval { return a.Sub(b) }

// Mul returns an enclosure of a * b. It is the free-function form of
// [Interval.Mul].
func Mul(a, b Interval) Interval { return a.Mul(b) }

// Div returns an enclosure of a / b. It is the free-function form of
// [Interval.Div].
func Div(a, b Interval) Interval { return a.Div(b) }

// Intersect returns the intersection of a and b. It is the free-function form
// of [Interval.Intersect].
func Intersect(a, b Interval) Interval { return a.Intersect(b) }

// Hull returns the convex hull of a and b. It is the free-function form of
// [Interval.Hull].
func Hull(a, b Interval) Interval { return a.Hull(b) }
