package auctions

import (
	"errors"
	"math"
)

// coalFixed pins a coalition's excess to a value in the nucleolus LPs.
type coalFixed struct {
	s   Coalition
	val float64
}

// addSplitPayoff writes the coefficient row for x(S) = Σ_{i∈S}(p_i - q_i) into a
// coefficient slice whose first 2n entries are the split payoff variables p and
// q, scaled by sign.
func addSplitPayoff(row []float64, s Coalition, n int, sign float64) {
	for i := 0; i < n; i++ {
		if s.Contains(i) {
			row[i] += sign
			row[n+i] -= sign
		}
	}
}

// efficiencyConstraint returns the efficiency equality x(N) = v(N) over nv
// variables.
func (g CoopGame) efficiencyConstraint(nv int) lpConstraint {
	n := g.Players
	row := make([]float64, nv)
	addSplitPayoff(row, FullCoalition(n), n, 1)
	return lpConstraint{coef: row, rel: relEQ, rhs: g.GrandValue()}
}

// irConstraints returns the individual-rationality constraints x_i >= v({i})
// over nv variables (used to restrict to imputations for the nucleolus).
func (g CoopGame) irConstraints(nv int) []lpConstraint {
	n := g.Players
	out := make([]lpConstraint, 0, n)
	for i := 0; i < n; i++ {
		row := make([]float64, nv)
		row[i] = 1
		row[n+i] = -1
		out = append(out, lpConstraint{coef: row, rel: relGE, rhs: g.Value(SingletonCoalition(i))})
	}
	return out
}

// solveEpsLP minimizes ε subject to efficiency, optional individual
// rationality, the fixed equality constraints (excess(S) = val) and, for every
// active coalition, excess(S) <= ε. It returns ε*, the accompanying allocation
// and whether the program was solved.
func (g CoopGame) solveEpsLP(fixed []coalFixed, active []Coalition, imputation bool) (float64, []float64, bool) {
	n := g.Players
	nv := 2*n + 2
	epsP, epsN := 2*n, 2*n+1
	c := make([]float64, nv)
	c[epsP] = 1
	c[epsN] = -1
	cons := []lpConstraint{g.efficiencyConstraint(nv)}
	if imputation {
		cons = append(cons, g.irConstraints(nv)...)
	}
	for _, f := range fixed {
		row := make([]float64, nv)
		addSplitPayoff(row, f.s, n, -1)
		cons = append(cons, lpConstraint{coef: row, rel: relEQ, rhs: f.val - g.Value(f.s)})
	}
	for _, s := range active {
		row := make([]float64, nv)
		addSplitPayoff(row, s, n, -1)
		row[epsP] -= 1
		row[epsN] += 1
		cons = append(cons, lpConstraint{coef: row, rel: relLE, rhs: -g.Value(s)})
	}
	res := lpMinimize(nv, c, cons)
	if !res.feasible || !res.bounded {
		return 0, nil, false
	}
	x := make([]float64, n)
	for i := 0; i < n; i++ {
		x[i] = res.x[i] - res.x[n+i]
	}
	return res.x[epsP] - res.x[epsN], x, true
}

// minCoalitionExcess minimizes the excess of coalition t subject to efficiency,
// optional individual rationality, the fixed equalities and active excess <=
// epsStar. It returns the minimal excess value.
func (g CoopGame) minCoalitionExcess(t Coalition, fixed []coalFixed, active []Coalition, epsStar float64, imputation bool) (float64, bool) {
	n := g.Players
	nv := 2 * n
	c := make([]float64, nv)
	addSplitPayoff(c, t, n, -1) // objective = -x(t); minimizing it minimizes excess(t) = v(t)-x(t)
	cons := []lpConstraint{g.efficiencyConstraint(nv)}
	if imputation {
		cons = append(cons, g.irConstraints(nv)...)
	}
	for _, f := range fixed {
		row := make([]float64, nv)
		addSplitPayoff(row, f.s, n, -1)
		cons = append(cons, lpConstraint{coef: row, rel: relEQ, rhs: f.val - g.Value(f.s)})
	}
	for _, s := range active {
		row := make([]float64, nv)
		addSplitPayoff(row, s, n, -1)
		cons = append(cons, lpConstraint{coef: row, rel: relLE, rhs: epsStar - g.Value(s)})
	}
	res := lpMinimize(nv, c, cons)
	if !res.feasible || !res.bounded {
		return 0, false
	}
	xt := 0.0
	for i := 0; i < n; i++ {
		if t.Contains(i) {
			xt += res.x[i] - res.x[n+i]
		}
	}
	return g.Value(t) - xt, true
}

// coordinateRange returns the minimum and maximum of x_k over the feasible set
// defined by efficiency, optional individual rationality, the fixed equalities
// and active excess <= epsStar.
func (g CoopGame) coordinateRange(k int, fixed []coalFixed, active []Coalition, epsStar float64, imputation bool) (lo, hi float64, ok bool) {
	n := g.Players
	nv := 2 * n
	build := func(sign float64) lpResult {
		c := make([]float64, nv)
		c[k] += sign
		c[n+k] -= sign
		cons := []lpConstraint{g.efficiencyConstraint(nv)}
		if imputation {
			cons = append(cons, g.irConstraints(nv)...)
		}
		for _, f := range fixed {
			r := make([]float64, nv)
			addSplitPayoff(r, f.s, n, -1)
			cons = append(cons, lpConstraint{coef: r, rel: relEQ, rhs: f.val - g.Value(f.s)})
		}
		for _, s := range active {
			r := make([]float64, nv)
			addSplitPayoff(r, s, n, -1)
			cons = append(cons, lpConstraint{coef: r, rel: relLE, rhs: epsStar - g.Value(s)})
		}
		return lpMinimize(nv, c, cons)
	}
	rlo := build(1)
	rhi := build(-1)
	if !rlo.feasible || !rlo.bounded || !rhi.feasible || !rhi.bounded {
		return 0, 0, false
	}
	lo = rlo.x[k] - rlo.x[n+k]
	hi = rhi.x[k] - rhi.x[n+k]
	return lo, hi, true
}

// CoreIsNonEmpty reports whether the core of the game is non-empty. The core is
// non-empty exactly when the least-core value ε* is <= 0 (Bondareva-Shapley).
func (g CoopGame) CoreIsNonEmpty() bool {
	eps, _, ok := g.LeastCore()
	if !ok {
		return false
	}
	return eps <= 1e-7
}

// IsBalanced reports whether the game is balanced, which by the
// Bondareva-Shapley theorem is equivalent to having a non-empty core.
func (g CoopGame) IsBalanced() bool { return g.CoreIsNonEmpty() }

// LeastCoreValue returns ε*, the least-core value: the minimum, over all
// efficient allocations, of the maximum coalition excess. It is <= 0 iff the
// core is non-empty.
func (g CoopGame) LeastCoreValue() (float64, bool) {
	eps, _, ok := g.LeastCore()
	return eps, ok
}

// LeastCore returns the least-core value ε* together with an allocation that
// attains it: the allocation minimizing the maximum excess, subject to
// efficiency.
func (g CoopGame) LeastCore() (float64, []float64, bool) {
	return g.solveEpsLP(nil, g.properCoalitions(), false)
}

// properCoalitions returns every non-empty proper coalition of the game.
func (g CoopGame) properCoalitions() []Coalition {
	n := g.Players
	size := 1 << uint(n)
	out := make([]Coalition, 0, size-2)
	for m := 1; m < size-1; m++ {
		out = append(out, Coalition(m))
	}
	return out
}

// nucleolusLike computes the (pre)nucleolus by the Kopelowitz sequence of linear
// programs. When imputation is true it optimizes over the imputation set,
// yielding the nucleolus; when false it optimizes over efficient allocations,
// yielding the prenucleolus.
func (g CoopGame) nucleolusLike(imputation bool) ([]float64, error) {
	n := g.Players
	if n == 1 {
		return []float64{g.GrandValue()}, nil
	}
	active := g.properCoalitions()
	var fixed []coalFixed
	const tol = 1e-7
	for round := 0; round < len(active)+n+2; round++ {
		epsStar, x, ok := g.solveEpsLP(fixed, active, imputation)
		if !ok {
			return nil, errors.New("auctions: nucleolus linear program failed")
		}
		var stillActive []Coalition
		for _, s := range active {
			minE, ok := g.minCoalitionExcess(s, fixed, active, epsStar, imputation)
			if !ok {
				return nil, errors.New("auctions: nucleolus sub-program failed")
			}
			if minE >= epsStar-tol {
				fixed = append(fixed, coalFixed{s: s, val: epsStar})
			} else {
				stillActive = append(stillActive, s)
			}
		}
		active = stillActive
		unique := true
		result := make([]float64, n)
		for k := 0; k < n; k++ {
			lo, hi, ok := g.coordinateRange(k, fixed, active, epsStar, imputation)
			if !ok {
				unique = false
				break
			}
			if hi-lo > 1e-6 {
				unique = false
				break
			}
			result[k] = 0.5 * (lo + hi)
		}
		if unique {
			return result, nil
		}
		if len(active) == 0 {
			return x, nil
		}
	}
	_, x, _ := g.solveEpsLP(fixed, active, imputation)
	return x, nil
}

// Nucleolus returns the nucleolus of the game: the unique imputation that
// lexicographically minimizes the sorted (descending) vector of coalition
// excesses. It is computed by a sequence of linear programs and always lies in
// the core when the core is non-empty. For a one-player game it returns {v(N)}.
func (g CoopGame) Nucleolus() ([]float64, error) { return g.nucleolusLike(true) }

// Prenucleolus returns the prenucleolus: the lexicographic minimizer of the
// sorted excess vector over all efficient allocations (without the individual-
// rationality restriction of the nucleolus). It coincides with the nucleolus
// whenever the prenucleolus is individually rational.
func (g CoopGame) Prenucleolus() ([]float64, error) { return g.nucleolusLike(false) }

// InEpsilonCore reports whether x lies in the (strong) ε-core: efficient with
// every coalition's excess at most eps (within tolerance tol).
func (g CoopGame) InEpsilonCore(x []float64, eps, tol float64) bool {
	if !g.IsEfficient(x, tol) {
		return false
	}
	n := g.Players
	size := 1 << uint(n)
	for m := 1; m < size-1; m++ {
		if g.Excess(Coalition(m), x) > eps+tol {
			return false
		}
	}
	return true
}

// CoreConstraintViolations returns the coalitions whose core inequality
// x(S) >= v(S) is violated by more than tol at the allocation x.
func (g CoopGame) CoreConstraintViolations(x []float64, tol float64) []Coalition {
	n := g.Players
	size := 1 << uint(n)
	var out []Coalition
	for m := 1; m < size; m++ {
		if g.Excess(Coalition(m), x) > tol {
			out = append(out, Coalition(m))
		}
	}
	return out
}

// approxEqual reports whether a and b differ by at most tol.
func approxEqual(a, b, tol float64) bool { return math.Abs(a-b) <= tol }
