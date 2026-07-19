package rootfind

import (
	"math"
	"math/cmplx"
	"sort"
)

// DefaultRootTol is the default convergence tolerance for the simultaneous
// polynomial root solvers.
const DefaultRootTol = 1e-14

// initialGuesses returns n distinct starting points spread on a circle whose
// radius is derived from the coefficient bounds of the monic polynomial m. The
// classic (0.4+0.9i)^k spiral of Durand and Kerner is used, scaled by the mean
// root radius to keep the guesses near the roots.
func initialGuesses(m CPoly, n int) []complex128 {
	// Radius estimate: geometric mean-ish based on constant term.
	radius := 1.0
	if c0 := cmplx.Abs(m[0]); c0 > 0 {
		radius = math.Pow(c0, 1.0/float64(n))
	}
	if radius == 0 || math.IsInf(radius, 0) || math.IsNaN(radius) {
		radius = 1
	}
	seed := complex(0.4, 0.9)
	z := make([]complex128, n)
	cur := complex(1, 0)
	for k := 0; k < n; k++ {
		cur *= seed
		// Place on a circle with a small angular offset per index.
		ang := 2*math.Pi*float64(k)/float64(n) + 0.25
		z[k] = complex(radius*math.Cos(ang), radius*math.Sin(ang)) + 0.001*cur
	}
	return z
}

// prepareMonic trims c, requires degree >= 1, and returns a monic copy together
// with its degree.
func prepareMonic(c CPoly) (CPoly, int, error) {
	d := c.Degree()
	if d < 1 {
		return nil, 0, ErrDegreeTooLow
	}
	m, err := c.Monic()
	if err != nil {
		return nil, 0, err
	}
	return m, d, nil
}

// DurandKerner finds all n complex roots of a degree-n polynomial simultaneously
// using the Durand-Kerner (Weierstrass) iteration. Every estimate is refined by
//
//	z_k <- z_k - p(z_k) / prod_{j!=k} (z_k - z_j)
//
// The method converges from the standard spiral initialization for essentially
// all polynomials. It returns the roots, the number of iterations used, and an
// error when the iteration budget is exhausted before convergence.
func DurandKerner(c CPoly, tol float64, maxIter int) ([]complex128, int, error) {
	if tol <= 0 {
		tol = DefaultRootTol
	}
	maxIter = resolveMax(maxIter)
	m, n, err := prepareMonic(c)
	if err != nil {
		return nil, 0, err
	}
	return DurandKernerWithInit(m, initialGuesses(m, n), tol, maxIter)
}

// DurandKernerWithInit runs the Durand-Kerner iteration from a caller-supplied
// set of initial guesses, one per root. This is useful when good approximate
// roots are already known, for example when polishing the output of another
// method. The number of guesses must equal the degree of c.
func DurandKernerWithInit(c CPoly, init []complex128, tol float64, maxIter int) ([]complex128, int, error) {
	m, n, err := prepareMonic(c)
	if err != nil {
		return nil, 0, err
	}
	if len(init) != n {
		return nil, 0, ErrBadInput
	}
	if tol <= 0 {
		tol = DefaultRootTol
	}
	maxIter = resolveMax(maxIter)
	z := make([]complex128, n)
	copy(z, init)
	for it := 1; it <= maxIter; it++ {
		maxDelta := 0.0
		for k := 0; k < n; k++ {
			num := m.Eval(z[k])
			den := complex(1, 0)
			for j := 0; j < n; j++ {
				if j != k {
					den *= z[k] - z[j]
				}
			}
			if den == 0 {
				continue
			}
			delta := num / den
			z[k] -= delta
			if a := cmplx.Abs(delta); a > maxDelta {
				maxDelta = a
			}
		}
		if maxDelta <= tol {
			sortComplex(z)
			return z, it, nil
		}
	}
	sortComplex(z)
	return z, maxIter, ErrNoConvergence
}

// AberthEhrlich finds all complex roots of a polynomial simultaneously using the
// Aberth-Ehrlich iteration, a third-order-per-step refinement of Durand-Kerner
// that uses the logarithmic derivative p'/p. It is typically faster and more
// robust than Durand-Kerner, converging in far fewer iterations for high-degree
// polynomials.
func AberthEhrlich(c CPoly, tol float64, maxIter int) ([]complex128, int, error) {
	if tol <= 0 {
		tol = DefaultRootTol
	}
	maxIter = resolveMax(maxIter)
	m, n, err := prepareMonic(c)
	if err != nil {
		return nil, 0, err
	}
	dm := m.Derivative()
	z := initialGuesses(m, n)
	for it := 1; it <= maxIter; it++ {
		maxDelta := 0.0
		w := make([]complex128, n)
		for k := 0; k < n; k++ {
			pv := m.Eval(z[k])
			dv := dm.Eval(z[k])
			if dv == 0 {
				w[k] = 0
				continue
			}
			ratio := pv / dv // Newton correction
			sum := complex(0, 0)
			for j := 0; j < n; j++ {
				if j != k {
					sum += 1 / (z[k] - z[j])
				}
			}
			denom := 1 - ratio*sum
			if denom == 0 {
				w[k] = ratio
			} else {
				w[k] = ratio / denom
			}
		}
		for k := 0; k < n; k++ {
			z[k] -= w[k]
			if a := cmplx.Abs(w[k]); a > maxDelta {
				maxDelta = a
			}
		}
		if maxDelta <= tol {
			sortComplex(z)
			return z, it, nil
		}
	}
	sortComplex(z)
	return z, maxIter, ErrNoConvergence
}

// Laguerre finds a single root of the polynomial c starting from x0 using
// Laguerre's method, which has cubic convergence for simple roots and is
// famously reliable: it converges to some root from almost any starting point,
// even for complex roots of a real polynomial. It operates in complex
// arithmetic throughout.
func Laguerre(c CPoly, x0 complex128, tol float64, maxIter int) (complex128, int, error) {
	if tol <= 0 {
		tol = DefaultRootTol
	}
	maxIter = resolveMax(maxIter)
	nDeg := c.Degree()
	if nDeg < 1 {
		return 0, 0, ErrDegreeTooLow
	}
	n := complex(float64(nDeg), 0)
	x := x0
	for it := 1; it <= maxIter; it++ {
		p, dp, ddp := c.EvalDeriv2(x)
		if cmplx.Abs(p) <= tol {
			return x, it, nil
		}
		if dp == 0 && ddp == 0 {
			// Nudge away from a flat spot.
			x += complex(1e-3, 1e-3)
			continue
		}
		g := dp / p
		g2 := g * g
		h := g2 - ddp/p
		disc := cmplx.Sqrt((n - 1) * (n*h - g2))
		dp1 := g + disc
		dp2 := g - disc
		denom := dp1
		if cmplx.Abs(dp2) > cmplx.Abs(dp1) {
			denom = dp2
		}
		if denom == 0 {
			x += complex(1e-3, 1e-3)
			continue
		}
		a := n / denom
		xn := x - a
		if cmplx.Abs(a) <= tol*(1+cmplx.Abs(xn)) {
			return xn, it, nil
		}
		x = xn
	}
	return x, maxIter, ErrNoConvergence
}

// LaguerreRoots finds all roots of the polynomial c by repeated application of
// [Laguerre] followed by deflation: each root found is divided out and the
// search continues on the reduced polynomial. Roots are polished on the original
// polynomial to undo deflation error. This is the classic reliable driver used
// by numerical libraries for arbitrary complex polynomials.
func LaguerreRoots(c CPoly, tol float64, maxIter int) ([]complex128, error) {
	if tol <= 0 {
		tol = DefaultRootTol
	}
	maxIter = resolveMax(maxIter)
	work := c.Trim().Clone()
	d := work.Degree()
	if d < 1 {
		return nil, ErrDegreeTooLow
	}
	roots := make([]complex128, 0, d)
	for work.Degree() >= 1 {
		x0 := complex(0, 0)
		r, _, err := Laguerre(work, x0, tol, maxIter)
		if err != nil {
			// Retry from a perturbed start before giving up.
			r, _, err = Laguerre(work, complex(0.5, 0.5), tol, maxIter)
			if err != nil {
				return roots, err
			}
		}
		// Polish on the original polynomial for accuracy.
		if rp, _, perr := Laguerre(c, r, tol, maxIter); perr == nil {
			r = rp
		}
		// Snap near-real roots to the real axis for real polynomials.
		roots = append(roots, r)
		work, _ = work.Deflate(r)
	}
	sortComplex(roots)
	return roots, nil
}

// PolyRoots returns all complex roots of the real polynomial p using the
// Aberth-Ehrlich method, the recommended general-purpose driver. The result is
// sorted by real part then imaginary part. It returns ErrDegreeTooLow for
// constant polynomials.
func PolyRoots(p Poly) ([]complex128, error) {
	roots, _, err := AberthEhrlich(p.ToComplex(), DefaultRootTol, 500)
	return roots, err
}

// CPolyRoots returns all complex roots of the complex polynomial c using the
// Aberth-Ehrlich method, sorted by real then imaginary part.
func CPolyRoots(c CPoly) ([]complex128, error) {
	roots, _, err := AberthEhrlich(c, DefaultRootTol, 500)
	return roots, err
}

// RealRoots returns the real roots of p, that is the roots whose imaginary part
// is negligible relative to imagTol. Each is refined by a Newton polish on the
// real polynomial. The result is sorted in ascending order.
func RealRoots(p Poly, imagTol float64) ([]float64, error) {
	if imagTol <= 0 {
		imagTol = 1e-8
	}
	roots, err := PolyRoots(p)
	if err != nil {
		return nil, err
	}
	dp := p.Derivative()
	var out []float64
	for _, z := range roots {
		if math.Abs(imag(z)) <= imagTol*(1+math.Abs(real(z))) {
			x := real(z)
			// One Newton polish step on the real line.
			if d := dp.Eval(x); d != 0 {
				x -= p.Eval(x) / d
			}
			out = append(out, x)
		}
	}
	sort.Float64s(out)
	return out, nil
}

// sortComplex sorts a slice of complex numbers in place by real part, breaking
// ties by imaginary part, giving root lists a deterministic order.
func sortComplex(z []complex128) {
	sort.Slice(z, func(i, j int) bool {
		if real(z[i]) != real(z[j]) {
			return real(z[i]) < real(z[j])
		}
		return imag(z[i]) < imag(z[j])
	})
}

// SortComplex returns a copy of z sorted by real part then imaginary part.
func SortComplex(z []complex128) []complex128 {
	out := make([]complex128, len(z))
	copy(out, z)
	sortComplex(out)
	return out
}

// MaxResidual returns the largest modulus |c(z)| over the given candidate roots,
// a scale-free measure of how well a computed root set satisfies the polynomial.
func MaxResidual(c CPoly, roots []complex128) float64 {
	m := 0.0
	for _, z := range roots {
		if r := cmplx.Abs(c.Eval(z)); r > m {
			m = r
		}
	}
	return m
}
