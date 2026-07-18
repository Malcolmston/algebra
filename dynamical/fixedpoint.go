package dynamical

import "math"

// Stability classifies the linear stability of a fixed or periodic point of a
// one-dimensional map according to the magnitude of its multiplier (the map
// derivative there).
type Stability int

const (
	// Unstable indicates a multiplier of magnitude strictly greater than one:
	// nearby orbits move away from the point.
	Unstable Stability = iota
	// Stable indicates a multiplier of magnitude strictly less than one:
	// nearby orbits are attracted to the point.
	Stable
	// Neutral indicates a multiplier of magnitude exactly one (to within the
	// requested tolerance): linear analysis is inconclusive.
	Neutral
	// SuperStable indicates a zero multiplier, giving especially fast
	// (quadratic) local convergence.
	SuperStable
)

// String returns a human-readable name for the stability class.
func (s Stability) String() string {
	switch s {
	case Unstable:
		return "unstable"
	case Stable:
		return "stable"
	case Neutral:
		return "neutral"
	case SuperStable:
		return "superstable"
	default:
		return "unknown"
	}
}

// Multiplier evaluates the derivative df at x, i.e. the multiplier of the map
// at the point x. For a fixed point this determines its linear stability.
func Multiplier(df Map1D, x float64) float64 { return df(x) }

// ClassifyStability classifies a fixed point by its multiplier m using
// tolerance tol to decide the borderline neutral case |m| == 1. A multiplier
// of exactly zero is reported as [SuperStable].
func ClassifyStability(m, tol float64) Stability {
	a := math.Abs(m)
	if a == 0 {
		return SuperStable
	}
	if math.Abs(a-1) <= tol {
		return Neutral
	}
	if a < 1 {
		return Stable
	}
	return Unstable
}

// IsFixedPoint reports whether x is a fixed point of f, that is whether
// |f(x) - x| <= tol.
func IsFixedPoint(f Map1D, x, tol float64) bool {
	return math.Abs(f(x)-x) <= tol
}

// FixedPoint1D finds a fixed point of the map f near x0 by applying Newton's
// method to g(x) = f(x) - x, using the analytic derivative df of f. It returns
// the located point, the number of iterations performed and whether the
// iteration converged to within tol.
func FixedPoint1D(f, df Map1D, x0 float64, maxIter int, tol float64) (root float64, iters int, converged bool) {
	x := x0
	for i := 1; i <= maxIter; i++ {
		g := f(x) - x
		if math.Abs(g) <= tol {
			return x, i, true
		}
		d := df(x) - 1
		if d == 0 {
			return x, i, false
		}
		x -= g / d
	}
	return x, maxIter, math.Abs(f(x)-x) <= tol
}

// LogisticFixedPoints returns the fixed points of the logistic map with
// parameter r. There is always the origin 0; the second fixed point 1 - 1/r
// exists (and is distinct) for r != 0 and is returned whenever r is nonzero.
func LogisticFixedPoints(r float64) []float64 {
	if r == 0 {
		return []float64{0}
	}
	return []float64{0, 1 - 1/r}
}

// LogisticStability returns the stability class of each fixed point returned by
// [LogisticFixedPoints] for parameter r, using tolerance tol for the neutral
// case. The multiplier at 0 is r and at 1 - 1/r is 2 - r.
func LogisticStability(r, tol float64) []Stability {
	pts := LogisticFixedPoints(r)
	out := make([]Stability, len(pts))
	for i, x := range pts {
		out[i] = ClassifyStability(LogisticDeriv(r, x), tol)
	}
	return out
}

// TentFixedPoints returns the fixed points of the tent map with slope mu.
// The origin 0 is always fixed; for mu > 1 the point mu/(1+mu) in the right
// branch is fixed as well and is included.
func TentFixedPoints(mu float64) []float64 {
	if mu <= 1 {
		return []float64{0}
	}
	return []float64{0, mu / (1 + mu)}
}

// HenonFixedPoints returns the real fixed points of the Henon map with
// parameters a and b. They are the solutions of a*x^2 + (1-b)*x - 1 = 0 with
// y = b*x; the slice is empty when the discriminant is negative and has one or
// two entries otherwise.
func HenonFixedPoints(a, b float64) []Point2D {
	if a == 0 {
		// Degenerate: (1-b)x = 1.
		if 1-b == 0 {
			return nil
		}
		x := 1 / (1 - b)
		return []Point2D{{x, b * x}}
	}
	disc := (1-b)*(1-b) + 4*a
	if disc < 0 {
		return nil
	}
	s := math.Sqrt(disc)
	x1 := (-(1 - b) + s) / (2 * a)
	x2 := (-(1 - b) - s) / (2 * a)
	if disc == 0 {
		return []Point2D{{x1, b * x1}}
	}
	return []Point2D{{x1, b * x1}, {x2, b * x2}}
}
