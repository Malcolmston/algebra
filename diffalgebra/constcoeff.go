package diffalgebra

import (
	"fmt"
	"math"
	"math/cmplx"
	"strings"
)

// ODETermKind classifies a basis solution of a constant-coefficient linear ODE.
type ODETermKind int

const (
	// RealExp is a term x^k e^{alpha x} coming from a real characteristic root.
	RealExp ODETermKind = iota
	// ComplexCos is a term x^k e^{alpha x} cos(beta x) from a complex pair.
	ComplexCos
	// ComplexSin is a term x^k e^{alpha x} sin(beta x) from a complex pair.
	ComplexSin
)

// ODETerm is one basis solution of a constant-coefficient linear ODE, of the
// form x^Power * e^{Alpha x} * trig(Beta x).
type ODETerm struct {
	Kind  ODETermKind
	Power int
	Alpha float64
	Beta  float64
}

// lambda returns the complex characteristic root associated with the term.
func (t ODETerm) lambda() complex128 { return complex(t.Alpha, t.Beta) }

// nthDerivValue returns the value of the d-th derivative of the term at x.
func (t ODETerm) nthDerivValue(d int, x float64) float64 {
	lam := t.lambda()
	k := t.Power
	sum := complex(0, 0)
	e := cmplx.Exp(lam * complex(x, 0))
	limit := d
	if k < limit {
		limit = k
	}
	for j := 0; j <= limit; j++ {
		term := complex(binom(d, j)*fallingFactorialFloat(k, j), 0)
		term *= powFloatC(x, k-j)
		term *= powC(lam, d-j)
		sum += term
	}
	sum *= e
	if t.Kind == ComplexSin {
		return imag(sum)
	}
	return real(sum)
}

// Eval returns the value of the basis term at x.
func (t ODETerm) Eval(x float64) float64 { return t.nthDerivValue(0, x) }

// String renders the term in readable form.
func (t ODETerm) String() string {
	var b strings.Builder
	if t.Power == 1 {
		b.WriteString("x*")
	} else if t.Power > 1 {
		fmt.Fprintf(&b, "x^%d*", t.Power)
	}
	if t.Alpha != 0 {
		fmt.Fprintf(&b, "e^(%g*x)", t.Alpha)
	}
	switch t.Kind {
	case ComplexCos:
		if t.Alpha != 0 {
			b.WriteString("*")
		}
		fmt.Fprintf(&b, "cos(%g*x)", t.Beta)
	case ComplexSin:
		if t.Alpha != 0 {
			b.WriteString("*")
		}
		fmt.Fprintf(&b, "sin(%g*x)", t.Beta)
	}
	s := b.String()
	if s == "" {
		return "1"
	}
	return strings.TrimSuffix(s, "*")
}

// ODESolution is the general homogeneous solution of a constant-coefficient
// linear ODE, described by its fundamental system of basis terms and the
// characteristic root clusters.
type ODESolution struct {
	Terms []ODETerm
	Roots []RootCluster
}

// Dimension returns the number of basis solutions (the order of the ODE).
func (s ODESolution) Dimension() int { return len(s.Terms) }

// Basis returns the fundamental system of basis terms.
func (s ODESolution) Basis() []ODETerm { return s.Terms }

// EvalBasis returns the values of every basis term at x.
func (s ODESolution) EvalBasis(x float64) []float64 {
	out := make([]float64, len(s.Terms))
	for i, t := range s.Terms {
		out[i] = t.Eval(x)
	}
	return out
}

// Evaluate returns sum_i consts[i] * basis_i(x). It panics only when consts has
// the wrong length is avoided by treating missing constants as zero.
func (s ODESolution) Evaluate(consts []float64, x float64) float64 {
	acc := 0.0
	for i, t := range s.Terms {
		if i < len(consts) {
			acc += consts[i] * t.Eval(x)
		}
	}
	return acc
}

// String renders the general solution as a linear combination C1*..+C2*.. .
func (s ODESolution) String() string {
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

// SolveLinearConstantODE builds the general solution of the homogeneous
// constant-coefficient linear ODE sum_i coeffs[i] y^(i) = 0, where coeffs[i] is
// the real coefficient of the i-th derivative. The characteristic roots are
// found numerically (seeded) and clustered into multiplicities using tol. It
// returns ErrDegree when the leading coefficient is zero or the equation is
// trivial.
func SolveLinearConstantODE(coeffs []float64, seed int64, tol float64) (ODESolution, error) {
	// trim trailing zeros
	n := len(coeffs)
	for n > 0 && coeffs[n-1] == 0 {
		n--
	}
	if n < 2 {
		return ODESolution{}, ErrDegree
	}
	cc := make([]complex128, n)
	for i := 0; i < n; i++ {
		cc[i] = complex(coeffs[i], 0)
	}
	roots := durandKerner(cc, seed)
	clusters := clusterRoots(roots, tol)
	sol := ODESolution{Roots: clusters}
	for _, cl := range clusters {
		re, im := real(cl.Value), imag(cl.Value)
		if math.Abs(im) < tol {
			for k := 0; k < cl.Mult; k++ {
				sol.Terms = append(sol.Terms, ODETerm{Kind: RealExp, Power: k, Alpha: re})
			}
		} else if im > 0 {
			for k := 0; k < cl.Mult; k++ {
				sol.Terms = append(sol.Terms, ODETerm{Kind: ComplexCos, Power: k, Alpha: re, Beta: im})
				sol.Terms = append(sol.Terms, ODETerm{Kind: ComplexSin, Power: k, Alpha: re, Beta: im})
			}
		}
	}
	return sol, nil
}

// SolveODEIVP solves the initial-value problem for the homogeneous
// constant-coefficient ODE sum_i coeffs[i] y^(i) = 0 with the initial data
// ic[d] = y^(d)(x0) for d = 0..n-1. It returns the fitted constants together
// with the general solution. It returns ErrDim when len(ic) does not match the
// order and ErrSingular when the Wronskian system is singular.
func SolveODEIVP(coeffs []float64, x0 float64, ic []float64, seed int64, tol float64) ([]float64, ODESolution, error) {
	sol, err := SolveLinearConstantODE(coeffs, seed, tol)
	if err != nil {
		return nil, ODESolution{}, err
	}
	n := sol.Dimension()
	if len(ic) != n {
		return nil, sol, ErrDim
	}
	m := make([][]float64, n)
	for d := 0; d < n; d++ {
		m[d] = make([]float64, n)
		for j := 0; j < n; j++ {
			m[d][j] = sol.Terms[j].nthDerivValue(d, x0)
		}
	}
	c, err := solveLinearFloat(m, ic)
	if err != nil {
		return nil, sol, err
	}
	return c, sol, nil
}

// CharacteristicComplexRoots returns the complex characteristic roots of the
// coefficient vector (coeffs[i] the coefficient of y^(i)).
func CharacteristicComplexRoots(coeffs []float64, seed int64) []complex128 {
	n := len(coeffs)
	for n > 0 && coeffs[n-1] == 0 {
		n--
	}
	if n < 2 {
		return nil
	}
	cc := make([]complex128, n)
	for i := 0; i < n; i++ {
		cc[i] = complex(coeffs[i], 0)
	}
	return durandKerner(cc, seed)
}

// solveLinearFloat solves the real linear system m x = b by Gaussian
// elimination with partial pivoting.
func solveLinearFloat(m [][]float64, b []float64) ([]float64, error) {
	n := len(m)
	a := make([][]float64, n)
	for i := range m {
		a[i] = make([]float64, n+1)
		copy(a[i], m[i])
		a[i][n] = b[i]
	}
	for col := 0; col < n; col++ {
		piv := col
		best := math.Abs(a[col][col])
		for r := col + 1; r < n; r++ {
			if v := math.Abs(a[r][col]); v > best {
				best = v
				piv = r
			}
		}
		if best < 1e-12 {
			return nil, ErrSingular
		}
		a[col], a[piv] = a[piv], a[col]
		for r := 0; r < n; r++ {
			if r == col {
				continue
			}
			factor := a[r][col] / a[col][col]
			for c := col; c <= n; c++ {
				a[r][c] -= factor * a[col][c]
			}
		}
	}
	x := make([]float64, n)
	for i := 0; i < n; i++ {
		x[i] = a[i][n] / a[i][i]
	}
	return x, nil
}

// --- small numeric helpers ---

func binom(n, k int) float64 {
	if k < 0 || k > n {
		return 0
	}
	res := 1.0
	for i := 0; i < k; i++ {
		res = res * float64(n-i) / float64(i+1)
	}
	return math.Round(res)
}

func fallingFactorialFloat(n, k int) float64 {
	res := 1.0
	for i := 0; i < k; i++ {
		res *= float64(n - i)
	}
	return res
}

func powFloatC(x float64, n int) complex128 {
	if n == 0 {
		return complex(1, 0)
	}
	return complex(math.Pow(x, float64(n)), 0)
}

func powC(z complex128, n int) complex128 {
	if n == 0 {
		return complex(1, 0)
	}
	if z == 0 {
		return complex(0, 0)
	}
	res := complex(1, 0)
	for i := 0; i < n; i++ {
		res *= z
	}
	return res
}
