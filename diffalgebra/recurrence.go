package diffalgebra

import (
	"fmt"
	"math"
	"math/cmplx"
	"strings"
)

// RecurrenceTermKind classifies a basis solution of a constant-coefficient
// linear recurrence.
type RecurrenceTermKind int

const (
	// RealGeom is a term n^k lambda^n from a real characteristic root.
	RealGeom RecurrenceTermKind = iota
	// ComplexGeomCos is a term n^k rho^n cos(n theta) from a complex pair.
	ComplexGeomCos
	// ComplexGeomSin is a term n^k rho^n sin(n theta) from a complex pair.
	ComplexGeomSin
)

// RecurrenceTerm is one basis solution of a linear recurrence, of the form
// n^Power * Rho^n * trig(n*Theta) (with Theta zero for real roots, in which
// case Rho carries the sign of the root).
type RecurrenceTerm struct {
	Kind  RecurrenceTermKind
	Power int
	Rho   float64 // signed root for real terms, modulus for complex terms
	Theta float64
}

// Eval returns the value of the basis term at integer index n.
func (t RecurrenceTerm) Eval(n int) float64 {
	np := math.Pow(float64(n), float64(t.Power))
	if n == 0 && t.Power == 0 {
		np = 1
	}
	switch t.Kind {
	case RealGeom:
		return np * math.Pow(t.Rho, float64(n))
	case ComplexGeomCos:
		return np * math.Pow(t.Rho, float64(n)) * math.Cos(float64(n)*t.Theta)
	case ComplexGeomSin:
		return np * math.Pow(t.Rho, float64(n)) * math.Sin(float64(n)*t.Theta)
	}
	return 0
}

// String renders the recurrence basis term.
func (t RecurrenceTerm) String() string {
	var b strings.Builder
	if t.Power == 1 {
		b.WriteString("n*")
	} else if t.Power > 1 {
		fmt.Fprintf(&b, "n^%d*", t.Power)
	}
	switch t.Kind {
	case RealGeom:
		fmt.Fprintf(&b, "(%g)^n", t.Rho)
	case ComplexGeomCos:
		fmt.Fprintf(&b, "(%g)^n*cos(%g*n)", t.Rho, t.Theta)
	case ComplexGeomSin:
		fmt.Fprintf(&b, "(%g)^n*sin(%g*n)", t.Rho, t.Theta)
	}
	return b.String()
}

// RecurrenceSolution is the general solution of a constant-coefficient linear
// recurrence, described by its fundamental basis terms and characteristic root
// clusters.
type RecurrenceSolution struct {
	Terms []RecurrenceTerm
	Roots []RootCluster
}

// Dimension returns the number of basis terms (the order of the recurrence).
func (s RecurrenceSolution) Dimension() int { return len(s.Terms) }

// Basis returns the fundamental basis terms.
func (s RecurrenceSolution) Basis() []RecurrenceTerm { return s.Terms }

// EvalBasis returns the values of every basis term at index n.
func (s RecurrenceSolution) EvalBasis(n int) []float64 {
	out := make([]float64, len(s.Terms))
	for i, t := range s.Terms {
		out[i] = t.Eval(n)
	}
	return out
}

// Evaluate returns sum_i consts[i] * basis_i(n).
func (s RecurrenceSolution) Evaluate(consts []float64, n int) float64 {
	acc := 0.0
	for i, t := range s.Terms {
		if i < len(consts) {
			acc += consts[i] * t.Eval(n)
		}
	}
	return acc
}

// String renders the general solution as a linear combination.
func (s RecurrenceSolution) String() string {
	if len(s.Terms) == 0 {
		return "0"
	}
	var b strings.Builder
	for i, t := range s.Terms {
		if i > 0 {
			b.WriteString(" + ")
		}
		fmt.Fprintf(&b, "C%d*%s", i+1, t.String())
	}
	return b.String()
}

// SolveLinearRecurrence builds the general solution of the homogeneous
// constant-coefficient linear recurrence sum_i coeffs[i] a_{n+i} = 0, where
// coeffs[i] is the real coefficient of a_{n+i}. The characteristic roots are
// found numerically (seeded) and clustered into multiplicities using tol. It
// returns ErrDegree when the recurrence has order below one.
func SolveLinearRecurrence(coeffs []float64, seed int64, tol float64) (RecurrenceSolution, error) {
	n := len(coeffs)
	for n > 0 && coeffs[n-1] == 0 {
		n--
	}
	if n < 2 {
		return RecurrenceSolution{}, ErrDegree
	}
	cc := make([]complex128, n)
	for i := 0; i < n; i++ {
		cc[i] = complex(coeffs[i], 0)
	}
	roots := durandKerner(cc, seed)
	clusters := clusterRoots(roots, tol)
	sol := RecurrenceSolution{Roots: clusters}
	for _, cl := range clusters {
		re, im := real(cl.Value), imag(cl.Value)
		if math.Abs(im) < tol {
			for k := 0; k < cl.Mult; k++ {
				sol.Terms = append(sol.Terms, RecurrenceTerm{Kind: RealGeom, Power: k, Rho: re})
			}
		} else if im > 0 {
			rho := cmplx.Abs(cl.Value)
			theta := math.Atan2(im, re)
			for k := 0; k < cl.Mult; k++ {
				sol.Terms = append(sol.Terms, RecurrenceTerm{Kind: ComplexGeomCos, Power: k, Rho: rho, Theta: theta})
				sol.Terms = append(sol.Terms, RecurrenceTerm{Kind: ComplexGeomSin, Power: k, Rho: rho, Theta: theta})
			}
		}
	}
	return sol, nil
}

// SolveRecurrenceIVP fits the constants of the general solution to the initial
// data initial[i] = a_i for i = 0..order-1. It returns the fitted constants
// together with the general solution. It returns ErrDim when len(initial) does
// not match the order and ErrSingular when the fitting matrix is singular.
func SolveRecurrenceIVP(coeffs []float64, initial []float64, seed int64, tol float64) ([]float64, RecurrenceSolution, error) {
	sol, err := SolveLinearRecurrence(coeffs, seed, tol)
	if err != nil {
		return nil, RecurrenceSolution{}, err
	}
	n := sol.Dimension()
	if len(initial) != n {
		return nil, sol, ErrDim
	}
	m := make([][]float64, n)
	for r := 0; r < n; r++ {
		m[r] = make([]float64, n)
		for j := 0; j < n; j++ {
			m[r][j] = sol.Terms[j].Eval(r)
		}
	}
	c, err := solveLinearFloat(m, initial)
	if err != nil {
		return nil, sol, err
	}
	return c, sol, nil
}

// RecurrenceValues iterates the recurrence directly to produce the first count
// terms from the given seed values, for cross-checking a closed-form solution.
// coeffs must be normalised so that the highest-order term can be solved for;
// it returns ErrDegree if the leading coefficient is zero and ErrDim if seeds
// does not provide order-1 initial values.
func RecurrenceValues(coeffs []float64, seeds []float64, count int) ([]float64, error) {
	n := len(coeffs)
	for n > 0 && coeffs[n-1] == 0 {
		n--
	}
	if n < 2 {
		return nil, ErrDegree
	}
	order := n - 1
	if len(seeds) != order {
		return nil, ErrDim
	}
	lead := coeffs[order]
	out := make([]float64, 0, count)
	out = append(out, seeds...)
	for len(out) < count {
		idx := len(out)
		// coeffs[order]*a_{k+order} + sum_{i<order} coeffs[i]*a_{k+i} = 0
		var s float64
		for i := 0; i < order; i++ {
			s += coeffs[i] * out[idx-order+i]
		}
		out = append(out, -s/lead)
	}
	return out[:count], nil
}
