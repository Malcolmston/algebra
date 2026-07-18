package interval

import "math"

// intervalElemUlps is the number of ULPs by which transcendental function
// results are inflated outward. The Go math library's elementary functions are
// accurate to well under this many ULPs on all supported platforms, so the
// inflation guarantees a rigorous enclosure.
const intervalElemUlps = 2

// Sqrt returns an enclosure of the square root of a. The negative part of a is
// discarded (sqrt is undefined there); if a is entirely negative the result is
// [Empty].
func (a Interval) Sqrt() Interval {
	if a.IsEmpty() || a.Hi < 0 {
		return Empty()
	}
	lo := 0.0
	if a.Lo > 0 {
		lo = a.Lo
	}
	return Interval{
		Lo: intervalDownN(math.Sqrt(lo), intervalElemUlps),
		Hi: intervalUpN(math.Sqrt(a.Hi), intervalElemUlps),
	}
}

// Exp returns an enclosure of e raised to the power of a. Because exp is
// monotone increasing the bounds come directly from the endpoints.
func (a Interval) Exp() Interval {
	if a.IsEmpty() {
		return Empty()
	}
	return Interval{
		Lo: intervalDownN(math.Exp(a.Lo), intervalElemUlps),
		Hi: intervalUpN(math.Exp(a.Hi), intervalElemUlps),
	}
}

// Log returns an enclosure of the natural logarithm of a. Only the positive
// part of a is in the domain: if Hi <= 0 the result is [Empty], and if Lo <= 0
// the lower bound is -Inf.
func (a Interval) Log() Interval {
	if a.IsEmpty() || a.Hi <= 0 {
		return Empty()
	}
	lo := math.Inf(-1)
	if a.Lo > 0 {
		lo = intervalDownN(math.Log(a.Lo), intervalElemUlps)
	}
	return Interval{Lo: lo, Hi: intervalUpN(math.Log(a.Hi), intervalElemUlps)}
}

// Sinh returns an enclosure of the hyperbolic sine of a. sinh is monotone
// increasing.
func (a Interval) Sinh() Interval {
	if a.IsEmpty() {
		return Empty()
	}
	return Interval{
		Lo: intervalDownN(math.Sinh(a.Lo), intervalElemUlps),
		Hi: intervalUpN(math.Sinh(a.Hi), intervalElemUlps),
	}
}

// Cosh returns an enclosure of the hyperbolic cosine of a. cosh is even and
// convex with minimum 1 at 0, so the lower bound comes from the mignitude and
// the upper bound from the magnitude.
func (a Interval) Cosh() Interval {
	if a.IsEmpty() {
		return Empty()
	}
	return Interval{
		Lo: intervalDownN(math.Cosh(a.Mig()), intervalElemUlps),
		Hi: intervalUpN(math.Cosh(a.Mag()), intervalElemUlps),
	}
}

// Tanh returns an enclosure of the hyperbolic tangent of a. tanh is monotone
// increasing.
func (a Interval) Tanh() Interval {
	if a.IsEmpty() {
		return Empty()
	}
	return Interval{
		Lo: intervalDownN(math.Tanh(a.Lo), intervalElemUlps),
		Hi: intervalUpN(math.Tanh(a.Hi), intervalElemUlps),
	}
}

// Atan returns an enclosure of the arctangent of a. atan is monotone
// increasing.
func (a Interval) Atan() Interval {
	if a.IsEmpty() {
		return Empty()
	}
	return Interval{
		Lo: intervalDownN(math.Atan(a.Lo), intervalElemUlps),
		Hi: intervalUpN(math.Atan(a.Hi), intervalElemUlps),
	}
}

// intervalTwoPi is an upper bound on 2*pi used to short-circuit periodic
// functions on wide intervals.
const intervalTwoPi = 6.283185307179587

// Cos returns an enclosure of the cosine of a. It accounts for the interior
// extrema of cosine (the points k*pi where cos is +-1) so the enclosure is
// tight, not merely [-1, 1].
func (a Interval) Cos() Interval {
	if a.IsEmpty() {
		return Empty()
	}
	if a.Width() >= intervalTwoPi {
		return Interval{Lo: -1, Hi: 1}
	}
	lo := math.Min(math.Cos(a.Lo), math.Cos(a.Hi))
	hi := math.Max(math.Cos(a.Lo), math.Cos(a.Hi))
	// cos(k*pi) = (-1)^k. Scan the integers k whose multiple of pi could lie in
	// [Lo, Hi]; a one-index margin keeps the scan sound against rounding in the
	// division. Including a spurious extremum only widens the result.
	kLo := int(math.Floor(a.Lo/math.Pi)) - 1
	kHi := int(math.Ceil(a.Hi/math.Pi)) + 1
	for k := kLo; k <= kHi; k++ {
		x := float64(k) * math.Pi
		if x < a.Lo || x > a.Hi {
			continue
		}
		if k%2 == 0 {
			hi = 1
		} else {
			lo = -1
		}
	}
	lo = math.Max(-1, intervalDownN(lo, intervalElemUlps))
	hi = math.Min(1, intervalUpN(hi, intervalElemUlps))
	return Interval{Lo: lo, Hi: hi}
}

// Sin returns an enclosure of the sine of a. As with [Interval.Cos] it accounts
// for the interior extrema at pi/2 + k*pi so the enclosure is tight.
func (a Interval) Sin() Interval {
	if a.IsEmpty() {
		return Empty()
	}
	if a.Width() >= intervalTwoPi {
		return Interval{Lo: -1, Hi: 1}
	}
	lo := math.Min(math.Sin(a.Lo), math.Sin(a.Hi))
	hi := math.Max(math.Sin(a.Lo), math.Sin(a.Hi))
	// Extrema of sin are at x = pi/2 + k*pi with sin = (-1)^k.
	kLo := int(math.Floor((a.Lo-math.Pi/2)/math.Pi)) - 1
	kHi := int(math.Ceil((a.Hi-math.Pi/2)/math.Pi)) + 1
	for k := kLo; k <= kHi; k++ {
		x := math.Pi/2 + float64(k)*math.Pi
		if x < a.Lo || x > a.Hi {
			continue
		}
		if k%2 == 0 {
			hi = 1
		} else {
			lo = -1
		}
	}
	lo = math.Max(-1, intervalDownN(lo, intervalElemUlps))
	hi = math.Min(1, intervalUpN(hi, intervalElemUlps))
	return Interval{Lo: lo, Hi: hi}
}

// Tan returns an enclosure of the tangent of a, computed as Sin(a) / Cos(a).
// Where the cosine enclosure contains 0 the tangent is unbounded and the result
// is [Entire].
func (a Interval) Tan() Interval {
	if a.IsEmpty() {
		return Empty()
	}
	return a.Sin().Div(a.Cos())
}

// IntPow returns an enclosure of a raised to the integer power n. It exploits
// the monotonicity of odd powers and the convex, even shape of positive even
// powers to produce a tight enclosure. For n == 0 the result is [1, 1]; for
// n < 0 it is the reciprocal of the positive power.
func (a Interval) IntPow(n int) Interval {
	if a.IsEmpty() {
		return Empty()
	}
	if n == 0 {
		return Point(1)
	}
	if n < 0 {
		return a.IntPow(-n).Recip()
	}
	fn := float64(n)
	if n%2 == 1 {
		// Odd power: monotone increasing.
		return Interval{
			Lo: intervalDownN(math.Pow(a.Lo, fn), intervalElemUlps),
			Hi: intervalUpN(math.Pow(a.Hi, fn), intervalElemUlps),
		}
	}
	// Even power: value is |x|^n, minimised at the mignitude.
	mig := a.Mig()
	mag := a.Mag()
	lo := intervalDownN(math.Pow(mig, fn), intervalElemUlps)
	if lo < 0 {
		lo = 0 // even powers are non-negative
	}
	return Interval{
		Lo: lo,
		Hi: intervalUpN(math.Pow(mag, fn), intervalElemUlps),
	}
}

// Pow returns an enclosure of a raised to the real power y, defined for a with
// positive lower bound via exp(y * log(a)). If a is not strictly positive the
// result is [Empty]. For integer exponents on intervals that may include
// non-positive values, use [Interval.IntPow].
func (a Interval) Pow(y float64) Interval {
	if a.IsEmpty() || a.Lo <= 0 {
		return Empty()
	}
	return a.Log().MulFloat(y).Exp()
}

// Sqrt returns an enclosure of the square root of a. It is the free-function
// form of [Interval.Sqrt].
func Sqrt(a Interval) Interval { return a.Sqrt() }

// Exp returns an enclosure of e**a. It is the free-function form of
// [Interval.Exp].
func Exp(a Interval) Interval { return a.Exp() }

// Log returns an enclosure of the natural logarithm of a. It is the
// free-function form of [Interval.Log].
func Log(a Interval) Interval { return a.Log() }

// Sin returns an enclosure of sin(a). It is the free-function form of
// [Interval.Sin].
func Sin(a Interval) Interval { return a.Sin() }

// Cos returns an enclosure of cos(a). It is the free-function form of
// [Interval.Cos].
func Cos(a Interval) Interval { return a.Cos() }

// Pow returns an enclosure of a**y for a real exponent y. It is the
// free-function form of [Interval.Pow].
func Pow(a Interval, y float64) Interval { return a.Pow(y) }
