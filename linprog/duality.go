package linprog

import "math"

// DualCanonical builds the dual of the canonical primal
//
//	minimize   c · x
//	subject to A x >= b
//	           x >= 0
//
// whose dual is
//
//	maximize   b · y
//	subject to A^T y <= c
//	           y >= 0.
//
// The returned [LP] has [Maximize] sense, one dual variable per primal
// constraint, and one dual constraint per primal variable. Strong duality
// guarantees the two optimal objective values coincide when both are finite.
func DualCanonical(c []float64, a [][]float64, b []float64) LP {
	at := Transpose(a)
	n := len(c)
	rel := make([]Relation, n)
	for i := range rel {
		rel[i] = LessEqual
	}
	return NewLP(Maximize, b, at, rel, c)
}

// DualityGap returns primalObj - dualObj. For a minimizing primal and its
// maximizing dual, weak duality makes this quantity nonnegative, and it is
// zero at optimality.
func DualityGap(primalObj, dualObj float64) float64 {
	return primalObj - dualObj
}

// WeakDualityHolds reports whether the primal (minimizing) and dual
// (maximizing) objective values respect weak duality, primalObj >= dualObj,
// within tol.
func WeakDualityHolds(primalObj, dualObj, tol float64) bool {
	return primalObj >= dualObj-tol
}

// ComplementarySlackness reports whether a primal point x and dual point y
// satisfy the complementary slackness conditions for the canonical pair from
// [DualCanonical], within tol. Specifically, for every constraint i either the
// dual variable y[i] is zero or the primal constraint A[i]·x = b[i] is tight,
// and for every variable j either x[j] is zero or the dual constraint
// (A^T y)[j] = c[j] is tight.
func ComplementarySlackness(c []float64, a [][]float64, b []float64, x, y []float64, tol float64) bool {
	// Primal-constraint / dual-variable pairs.
	ax := MatVec(a, x)
	for i := range b {
		if math.Abs(y[i]) > tol && math.Abs(ax[i]-b[i]) > tol {
			return false
		}
	}
	// Variable / dual-constraint pairs.
	aty := MatTVec(a, y)
	for j := range c {
		if math.Abs(x[j]) > tol && math.Abs(aty[j]-c[j]) > tol {
			return false
		}
	}
	return true
}

// LPKKT holds the four first-order optimality residuals for a canonical
// primal/dual LP pair. Every field is a nonnegative violation magnitude; all
// are zero at an optimal primal-dual pair.
type LPKKT struct {
	// PrimalFeas is the maximum violation of A x >= b and x >= 0.
	PrimalFeas float64
	// DualFeas is the maximum violation of A^T y <= c and y >= 0.
	DualFeas float64
	// Stationarity is |c·x - b·y|, the duality gap magnitude.
	Stationarity float64
	// CompSlack is the maximum complementary-slackness product magnitude.
	CompSlack float64
}

// Max returns the largest of the four residuals.
func (k LPKKT) Max() float64 {
	return math.Max(math.Max(k.PrimalFeas, k.DualFeas), math.Max(k.Stationarity, k.CompSlack))
}

// Satisfied reports whether every residual is within tol.
func (k LPKKT) Satisfied(tol float64) bool { return k.Max() <= tol }

// CheckKKTLP computes the [LPKKT] residuals for the canonical primal
//
//	minimize c·x s.t. A x >= b, x >= 0
//
// with dual variables y (for the constraints) at primal point x.
func CheckKKTLP(c []float64, a [][]float64, b []float64, x, y []float64) LPKKT {
	var k LPKKT
	ax := MatVec(a, x)
	for i := range b {
		if v := b[i] - ax[i]; v > k.PrimalFeas {
			k.PrimalFeas = v
		}
		if v := -y[i]; v > k.DualFeas {
			k.DualFeas = v
		}
	}
	for _, xj := range x {
		if v := -xj; v > k.PrimalFeas {
			k.PrimalFeas = v
		}
	}
	aty := MatTVec(a, y)
	for j := range c {
		if v := aty[j] - c[j]; v > k.DualFeas {
			k.DualFeas = v
		}
	}
	k.Stationarity = math.Abs(Dot(c, x) - Dot(b, y))
	// Complementary slackness magnitudes.
	for i := range b {
		if p := math.Abs(y[i] * (ax[i] - b[i])); p > k.CompSlack {
			k.CompSlack = p
		}
	}
	for j := range c {
		if p := math.Abs(x[j] * (aty[j] - c[j])); p > k.CompSlack {
			k.CompSlack = p
		}
	}
	return k
}
