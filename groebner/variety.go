package groebner

import (
	"errors"
	"math/cmplx"
	"math/rand"
)

// ErrNotZeroDimensional is returned by the variety solver when the ideal has
// infinitely many solutions (positive-dimensional) and cannot be solved by
// finite back-substitution.
var ErrNotZeroDimensional = errors.New("groebner: system is not zero-dimensional")

func ratToComplex(r interface{ Float64() (float64, bool) }) complex128 {
	f, _ := r.Float64()
	return complex(f, 0)
}

func cpow(z complex128, k int) complex128 {
	r := complex(1, 0)
	for i := 0; i < k; i++ {
		r *= z
	}
	return r
}

// SolveUnivariate returns the complex roots of a univariate polynomial given by
// its coefficient slice, where coeffs[k] is the coefficient of x^k. It uses the
// Durand–Kerner (Weierstrass) iteration seeded by the supplied random seed so
// results are reproducible. Trailing (high-degree) near-zero coefficients are
// dropped. A constant nonzero polynomial has no roots; the zero polynomial
// returns no roots as well.
func SolveUnivariate(coeffs []complex128, seed int64) []complex128 {
	// Strip leading (highest-degree) zeros.
	d := len(coeffs) - 1
	for d >= 0 && cmplx.Abs(coeffs[d]) < 1e-14 {
		d--
	}
	if d <= 0 {
		return nil
	}
	// Monic coefficients a[k] for k=0..d, a[d]=1.
	a := make([]complex128, d+1)
	lead := coeffs[d]
	for k := 0; k <= d; k++ {
		a[k] = coeffs[k] / lead
	}
	eval := func(z complex128) complex128 {
		res := a[d]
		for k := d - 1; k >= 0; k-- {
			res = res*z + a[k]
		}
		return res
	}
	rng := rand.New(rand.NewSource(seed))
	roots := make([]complex128, d)
	seed0 := complex(0.4, 0.9)
	for i := 0; i < d; i++ {
		roots[i] = cpow(seed0, i) * complex(1+0.01*rng.Float64(), 0.01*rng.Float64())
	}
	for iter := 0; iter < 2000; iter++ {
		maxDelta := 0.0
		for i := 0; i < d; i++ {
			denom := complex(1, 0)
			for j := 0; j < d; j++ {
				if j != i {
					denom *= roots[i] - roots[j]
				}
			}
			if cmplx.Abs(denom) < 1e-300 {
				continue
			}
			delta := eval(roots[i]) / denom
			roots[i] -= delta
			if ad := cmplx.Abs(delta); ad > maxDelta {
				maxDelta = ad
			}
		}
		if maxDelta < 1e-13 {
			break
		}
	}
	return roots
}

// SolveZeroDimensional numerically approximates all complex solutions of the
// zero-dimensional system defined by gens. It computes a lexicographic Gröbner
// basis (which for a zero-dimensional ideal has triangular shape), solves the
// univariate polynomial in the last variable, and back-substitutes recursively.
// Candidate tuples are verified against every generator within the given
// tolerance. The seed makes the underlying root finder reproducible. It returns
// ErrNotZeroDimensional if the system has infinitely many solutions and the
// unit ideal (no solutions) yields an empty slice.
func SolveZeroDimensional(gens []Poly, seed int64, tol float64) ([][]complex128, error) {
	if len(gens) == 0 {
		return nil, ErrNotZeroDimensional
	}
	n := gens[0].nvars
	id := NewIdealN(n, LexOrder(), gens...)
	if id.IsUnit() {
		return nil, nil // no solutions
	}
	if !id.IsZeroDimensional() {
		return nil, ErrNotZeroDimensional
	}
	gb := ReducedGroebnerBasis(gens, LexOrder())

	var solutions [][]complex128
	values := make([]complex128, n)
	assigned := make([]bool, n)

	var rec func(level int)
	rec = func(level int) {
		if level < 0 {
			// Verify against every generator.
			for _, g := range gens {
				if g.EvalComplexAbs(values) > tolCheck(tol) {
					return
				}
			}
			sol := make([]complex128, n)
			copy(sol, values)
			solutions = append(solutions, sol)
			return
		}
		coeffs, deg := pickUnivariate(gb, level, values, assigned)
		if coeffs == nil {
			// No univariate constraint found: treat as free -> not truly
			// zero-dimensional in this branch; abort without solutions here.
			return
		}
		if deg == 0 {
			// Nonzero constant constraint: inconsistent branch.
			return
		}
		roots := SolveUnivariate(coeffs, seed+int64(level)*7919)
		for _, root := range roots {
			values[level] = cleanupComplex(root)
			assigned[level] = true
			rec(level - 1)
		}
		assigned[level] = false
	}
	rec(n - 1)
	return dedupSolutions(solutions, tol), nil
}

func tolCheck(tol float64) float64 {
	if tol <= 0 {
		tol = 1e-6
	}
	return tol * 1e3
}

// pickUnivariate finds, among the Gröbner basis, a polynomial that after
// substituting the already-assigned variables becomes univariate in the
// variable at index level, returning its coefficient vector (indexed by the
// exponent of that variable) and its degree. It selects the constraint of
// smallest positive degree to limit branching; a nonzero constant constraint is
// reported with degree 0. It returns (nil, -1) if no such polynomial exists.
func pickUnivariate(gb []Poly, level int, values []complex128, assigned []bool) ([]complex128, int) {
	var best []complex128
	bestDeg := -1
	for _, p := range gb {
		coeffs, ok := substituteToUnivariate(p, level, values, assigned)
		if !ok {
			continue
		}
		deg := len(coeffs) - 1
		for deg >= 0 && cmplx.Abs(coeffs[deg]) < 1e-12 {
			deg--
		}
		if deg < 0 {
			continue // vanished identically
		}
		if deg == 0 {
			return coeffs[:1], 0 // inconsistent constant: prune immediately
		}
		if bestDeg == -1 || deg < bestDeg {
			bestDeg = deg
			best = coeffs[:deg+1]
		}
	}
	if bestDeg == -1 {
		return nil, -1
	}
	return best, bestDeg
}

// substituteToUnivariate substitutes the assigned variables into p and checks
// that only the variable at index level remains among the unassigned ones. On
// success it returns the coefficient vector indexed by the exponent of level.
func substituteToUnivariate(p Poly, level int, values []complex128, assigned []bool) ([]complex128, bool) {
	deg := p.DegreeIn(level)
	coeffs := make([]complex128, deg+1)
	for _, t := range p.terms {
		for v := 0; v < p.nvars; v++ {
			if v == level {
				continue
			}
			if t.Mono[v] > 0 && !assigned[v] {
				return nil, false
			}
		}
		c := ratToComplex(t.Coeff)
		for v := 0; v < p.nvars; v++ {
			if v == level || t.Mono[v] == 0 {
				continue
			}
			c *= cpow(values[v], t.Mono[v])
		}
		coeffs[t.Mono[level]] += c
	}
	return coeffs, true
}

func cleanupComplex(z complex128) complex128 {
	re, im := real(z), imag(z)
	if im < 1e-11 && im > -1e-11 {
		im = 0
	}
	if re < 1e-11 && re > -1e-11 {
		re = 0
	}
	return complex(re, im)
}

func dedupSolutions(sols [][]complex128, tol float64) [][]complex128 {
	if tol <= 0 {
		tol = 1e-6
	}
	var out [][]complex128
	for _, s := range sols {
		dup := false
		for _, e := range out {
			same := true
			for i := range s {
				if cmplx.Abs(s[i]-e[i]) > tol {
					same = false
					break
				}
			}
			if same {
				dup = true
				break
			}
		}
		if !dup {
			out = append(out, s)
		}
	}
	return out
}

// Solve numerically approximates all complex solutions of the ideal's variety,
// assuming the ideal is zero-dimensional. It is a convenience wrapper around
// SolveZeroDimensional using the ideal's generators.
func (id Ideal) Solve(seed int64, tol float64) ([][]complex128, error) {
	return SolveZeroDimensional(id.gens, seed, tol)
}

// RealSolutions filters a set of complex solutions, returning the real parts of
// those solutions whose imaginary parts are all below the tolerance imagTol.
func RealSolutions(sols [][]complex128, imagTol float64) [][]float64 {
	var out [][]float64
	for _, s := range sols {
		real := true
		for _, z := range s {
			if imag(z) > imagTol || imag(z) < -imagTol {
				real = false
				break
			}
		}
		if !real {
			continue
		}
		row := make([]float64, len(s))
		for i, z := range s {
			row[i] = realPart(z)
		}
		out = append(out, row)
	}
	return out
}

func realPart(z complex128) float64 { return real(z) }
