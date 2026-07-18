// Package interval implements rigorous interval arithmetic over IEEE-754
// double precision using outward directed rounding.
//
// An [Interval] represents the set of real numbers {x : Lo <= x <= Hi}. Every
// operation returns an enclosure: an interval that is guaranteed to contain the
// true result of the corresponding real operation applied to every point of the
// inputs. This makes the package suitable for verified numerics, where a
// computed answer must come with a mathematically sound error bound rather than
// a heuristic estimate.
//
// # Directed rounding
//
// The Go standard library does not expose the hardware rounding-mode control
// registers, so this package obtains outward rounding portably with
// [math.Nextafter]. After a floating-point operation whose exact result is r,
// the computed value fl(r) differs from r by at most half an ULP for the
// correctly rounded primitives (+, -, *, /, sqrt). Nudging the lower bound one
// ULP toward -Inf and the upper bound one ULP toward +Inf therefore yields an
// interval that provably contains r. Elementary transcendental functions
// (exp, log, sin, ...) are not guaranteed correctly rounded, so their results
// are inflated by a small fixed number of ULPs to remain rigorous. The bounds
// are wider than strictly necessary but never unsound.
//
// # Special values
//
// The empty set is represented by [Empty] and reported by [Interval.IsEmpty].
// The whole real line is [Entire]. Operations propagate emptiness: any
// operation with an empty operand returns [Empty]. Unbounded intervals use
// infinities as their bounds.
//
// # Contents
//
//   - Interval arithmetic: [Interval.Add], [Interval.Sub], [Interval.Mul],
//     [Interval.Div] and scalar variants, all outward rounded.
//   - Elementary functions with enclosures: [Interval.Exp], [Interval.Log],
//     [Interval.Sqrt], [Interval.Sin], [Interval.Cos], [Interval.Tan],
//     [Interval.Sinh], [Interval.Cosh], [Interval.Tanh], [Interval.Atan],
//     [Interval.IntPow] and [Interval.Pow].
//   - Set operations and measures: [Interval.Width], [Interval.Midpoint],
//     [Interval.Radius], [Interval.Intersect], [Interval.Hull],
//     [Interval.Contains].
//   - Verified root finding: [Newton] and [Bisect] using the interval Newton
//     operator.
//   - Interval matrices: [Matrix] with [Matrix.Add], [Matrix.Mul],
//     [Matrix.MulVec] and related operations.
//
// All results are deterministic and depend only on the inputs; the package uses
// nothing outside the Go standard library.
package interval
