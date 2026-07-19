package metaheuristics

import (
	"errors"
	"math"
)

// ObjectiveFunc is a real-valued function of a vector argument. Every optimizer
// in the package searches for an argument that minimizes such a function.
type ObjectiveFunc func(x []float64) float64

// ErrDimMismatch is returned when the dimension of a supplied vector does not
// match the dimension of the problem being solved.
var ErrDimMismatch = errors.New("metaheuristics: dimension mismatch")

// ErrEmptyBounds is returned when a [Bounds] with zero dimension is used where
// a positive dimension is required.
var ErrEmptyBounds = errors.New("metaheuristics: empty bounds")

// ErrInvalidConfig is returned when an optimizer configuration is invalid, for
// example a non-positive population size or iteration count.
var ErrInvalidConfig = errors.New("metaheuristics: invalid configuration")

// Bounds describes an axis-aligned box in R^n by its per-coordinate lower and
// upper limits. The two slices must have equal length, which is the dimension.
type Bounds struct {
	Lower []float64
	Upper []float64
}

// NewBounds constructs a [Bounds] from separate lower and upper slices. It
// returns [ErrDimMismatch] if the slices differ in length and [ErrEmptyBounds]
// if they are empty.
func NewBounds(lower, upper []float64) (Bounds, error) {
	if len(lower) != len(upper) {
		return Bounds{}, ErrDimMismatch
	}
	if len(lower) == 0 {
		return Bounds{}, ErrEmptyBounds
	}
	lo := make([]float64, len(lower))
	hi := make([]float64, len(upper))
	copy(lo, lower)
	copy(hi, upper)
	return Bounds{Lower: lo, Upper: hi}, nil
}

// UniformBounds returns a [Bounds] of the given dimension in which every
// coordinate ranges over [lo, hi].
func UniformBounds(dim int, lo, hi float64) Bounds {
	l := make([]float64, dim)
	u := make([]float64, dim)
	for i := range l {
		l[i] = lo
		u[i] = hi
	}
	return Bounds{Lower: l, Upper: u}
}

// Dim returns the dimension of the box.
func (b Bounds) Dim() int { return len(b.Lower) }

// Valid reports whether the bounds are well formed: non-empty, of equal length
// and with every lower limit not exceeding its upper limit.
func (b Bounds) Valid() bool {
	if len(b.Lower) == 0 || len(b.Lower) != len(b.Upper) {
		return false
	}
	for i := range b.Lower {
		if b.Lower[i] > b.Upper[i] {
			return false
		}
	}
	return true
}

// Center returns the geometric center of the box.
func (b Bounds) Center() []float64 {
	c := make([]float64, b.Dim())
	for i := range c {
		c[i] = 0.5 * (b.Lower[i] + b.Upper[i])
	}
	return c
}

// Width returns the per-coordinate widths (Upper-Lower) of the box.
func (b Bounds) Width() []float64 {
	w := make([]float64, b.Dim())
	for i := range w {
		w[i] = b.Upper[i] - b.Lower[i]
	}
	return w
}

// Contains reports whether x lies within the closed box. It returns false if x
// has the wrong dimension.
func (b Bounds) Contains(x []float64) bool {
	if len(x) != b.Dim() {
		return false
	}
	for i := range x {
		if x[i] < b.Lower[i] || x[i] > b.Upper[i] {
			return false
		}
	}
	return true
}

// Clip returns a copy of x with each coordinate clamped into the box.
func (b Bounds) Clip(x []float64) []float64 {
	out := make([]float64, len(x))
	for i := range x {
		out[i] = Clamp(x[i], b.Lower[i], b.Upper[i])
	}
	return out
}

// ClipInPlace clamps each coordinate of x into the box, modifying x, and
// returns x for convenience.
func (b Bounds) ClipInPlace(x []float64) []float64 {
	for i := range x {
		x[i] = Clamp(x[i], b.Lower[i], b.Upper[i])
	}
	return x
}

// Reflect returns a copy of x in which coordinates outside the box are
// reflected back inside across the violated face. Repeated reflection keeps the
// result within the box for arbitrary overshoot.
func (b Bounds) Reflect(x []float64) []float64 {
	out := make([]float64, len(x))
	for i := range x {
		lo, hi := b.Lower[i], b.Upper[i]
		v := x[i]
		if hi == lo {
			out[i] = lo
			continue
		}
		span := hi - lo
		for v < lo || v > hi {
			if v < lo {
				v = lo + (lo - v)
			}
			if v > hi {
				v = hi - (v - hi)
			}
			_ = span
		}
		out[i] = v
	}
	return out
}

// Clamp returns x clamped to the closed interval [lo, hi]. If hi < lo the
// bounds are treated as swapped.
func Clamp(x, lo, hi float64) float64 {
	if hi < lo {
		lo, hi = hi, lo
	}
	if x < lo {
		return lo
	}
	if x > hi {
		return hi
	}
	return x
}

// Result records the outcome of an optimization run.
type Result struct {
	// X is the best argument found (a fresh copy owned by the caller).
	X []float64
	// F is the objective value at X.
	F float64
	// Iterations is the number of iterations actually performed.
	Iterations int
	// Evaluations is the number of objective-function evaluations performed.
	Evaluations int
	// Converged reports whether a configured convergence criterion was met.
	Converged bool
	// History holds the best objective value at the end of each iteration when
	// history recording is enabled; otherwise it is nil.
	History []float64
}

// Negate returns a new objective equal to -f(x). Minimizing the result is
// equivalent to maximizing f.
func Negate(f ObjectiveFunc) ObjectiveFunc {
	return func(x []float64) float64 { return -f(x) }
}

// Shift returns a new objective f(x - offset), translating the landscape so
// that the original optimum at o moves to o+offset.
func Shift(f ObjectiveFunc, offset []float64) ObjectiveFunc {
	off := append([]float64(nil), offset...)
	return func(x []float64) float64 {
		y := make([]float64, len(x))
		for i := range x {
			if i < len(off) {
				y[i] = x[i] - off[i]
			} else {
				y[i] = x[i]
			}
		}
		return f(y)
	}
}

// Scale returns a new objective a*f(x)+b.
func Scale(f ObjectiveFunc, a, b float64) ObjectiveFunc {
	return func(x []float64) float64 { return a*f(x) + b }
}

// Penalized returns a new objective that adds weight times the squared amount
// by which x violates the box b to the value of f. This turns a bounded
// problem into an unconstrained one suitable for optimizers that do not clip.
func Penalized(f ObjectiveFunc, b Bounds, weight float64) ObjectiveFunc {
	return func(x []float64) float64 {
		pen := 0.0
		for i := range x {
			if i >= b.Dim() {
				break
			}
			if x[i] < b.Lower[i] {
				d := b.Lower[i] - x[i]
				pen += d * d
			} else if x[i] > b.Upper[i] {
				d := x[i] - b.Upper[i]
				pen += d * d
			}
		}
		return f(x) + weight*pen
	}
}

// VecCopy returns a fresh copy of x.
func VecCopy(x []float64) []float64 {
	out := make([]float64, len(x))
	copy(out, x)
	return out
}

// VecAdd returns the element-wise sum a+b. The result has length min(len a, len b).
func VecAdd(a, b []float64) []float64 {
	n := minInt(len(a), len(b))
	out := make([]float64, n)
	for i := 0; i < n; i++ {
		out[i] = a[i] + b[i]
	}
	return out
}

// VecSub returns the element-wise difference a-b.
func VecSub(a, b []float64) []float64 {
	n := minInt(len(a), len(b))
	out := make([]float64, n)
	for i := 0; i < n; i++ {
		out[i] = a[i] - b[i]
	}
	return out
}

// VecScale returns the vector s*a.
func VecScale(a []float64, s float64) []float64 {
	out := make([]float64, len(a))
	for i := range a {
		out[i] = s * a[i]
	}
	return out
}

// VecAXPY returns a + s*b (a scaled vector accumulation).
func VecAXPY(a []float64, s float64, b []float64) []float64 {
	n := minInt(len(a), len(b))
	out := make([]float64, n)
	for i := 0; i < n; i++ {
		out[i] = a[i] + s*b[i]
	}
	return out
}

// VecDot returns the Euclidean inner product of a and b.
func VecDot(a, b []float64) float64 {
	n := minInt(len(a), len(b))
	s := 0.0
	for i := 0; i < n; i++ {
		s += a[i] * b[i]
	}
	return s
}

// VecNorm returns the Euclidean (L2) norm of x.
func VecNorm(x []float64) float64 {
	return math.Sqrt(VecDot(x, x))
}

// VecNorm1 returns the L1 (taxicab) norm of x.
func VecNorm1(x []float64) float64 {
	s := 0.0
	for _, v := range x {
		s += math.Abs(v)
	}
	return s
}

// VecNormInf returns the L-infinity (maximum absolute value) norm of x.
func VecNormInf(x []float64) float64 {
	m := 0.0
	for _, v := range x {
		if a := math.Abs(v); a > m {
			m = a
		}
	}
	return m
}

// VecDist returns the Euclidean distance between a and b.
func VecDist(a, b []float64) float64 {
	n := minInt(len(a), len(b))
	s := 0.0
	for i := 0; i < n; i++ {
		d := a[i] - b[i]
		s += d * d
	}
	return math.Sqrt(s)
}

// VecMean returns the coordinate-wise mean of the rows of xs. It returns nil if
// xs is empty.
func VecMean(xs [][]float64) []float64 {
	if len(xs) == 0 {
		return nil
	}
	d := len(xs[0])
	out := make([]float64, d)
	for _, x := range xs {
		for i := 0; i < d && i < len(x); i++ {
			out[i] += x[i]
		}
	}
	inv := 1.0 / float64(len(xs))
	for i := range out {
		out[i] *= inv
	}
	return out
}

// Lerp returns the linear interpolation (1-t)*a + t*b.
func Lerp(a, b, t float64) float64 { return a + t*(b-a) }

func minInt(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func maxInt(a, b int) int {
	if a > b {
		return a
	}
	return b
}
