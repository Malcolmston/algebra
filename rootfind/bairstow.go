package rootfind

import (
	"math"
	"math/cmplx"
)

// QuadraticRoots returns the two roots of the quadratic a*x^2 + b*x + c using a
// numerically stable formula that avoids catastrophic cancellation. The roots
// are returned as complex numbers; real roots have zero imaginary part.
func QuadraticRoots(a, b, c float64) (complex128, complex128) {
	if a == 0 {
		if b == 0 {
			return cmplx.NaN(), cmplx.NaN()
		}
		return complex(-c/b, 0), complex(-c/b, 0)
	}
	disc := b*b - 4*a*c
	if disc >= 0 {
		s := math.Sqrt(disc)
		// Stable: pick the sign of b to avoid cancellation.
		var q float64
		if b >= 0 {
			q = -0.5 * (b + s)
		} else {
			q = -0.5 * (b - s)
		}
		r1 := q / a
		if q == 0 {
			return complex(r1, 0), complex(0, 0)
		}
		r2 := c / q
		return complex(r1, 0), complex(r2, 0)
	}
	s := math.Sqrt(-disc)
	re := -b / (2 * a)
	im := s / (2 * a)
	return complex(re, im), complex(re, -im)
}

// Bairstow extracts a real quadratic factor x^2 + u*x + v from the real
// polynomial p, starting from initial guesses u0, v0, using Bairstow's method.
// It returns the coefficients u and v of the converged factor, from which the
// corresponding pair of (possibly complex-conjugate) roots follows via
// [QuadraticRoots]. Bairstow's method finds complex roots of a real polynomial
// using only real arithmetic.
func Bairstow(p Poly, u0, v0, tol float64, maxIter int) (u, v float64, iters int, err error) {
	if tol <= 0 {
		tol = 1e-14
	}
	maxIter = resolveMax(maxIter)
	n := p.Degree()
	if n < 2 {
		return 0, 0, 0, ErrDegreeTooLow
	}
	a := make([]float64, n+1)
	copy(a, p[:n+1])
	u, v = u0, v0
	b := make([]float64, n+1)
	cc := make([]float64, n+1)
	for it := 1; it <= maxIter; it++ {
		// Synthetic division of a by x^2 + u x + v (coefficients ascending,
		// processed from high to low degree).
		b[n] = a[n]
		if n-1 >= 0 {
			b[n-1] = a[n-1] - u*b[n]
		}
		for i := n - 2; i >= 0; i-- {
			b[i] = a[i] - u*b[i+1] - v*b[i+2]
		}
		cc[n] = b[n]
		if n-1 >= 0 {
			cc[n-1] = b[n-1] - u*cc[n]
		}
		for i := n - 2; i >= 1; i-- {
			cc[i] = b[i] - u*cc[i+1] - v*cc[i+2]
		}
		// The remainder is b[1]*x + b[0]; solve the 2x2 Newton system.
		c1 := cc[1]
		c2 := cc[2]
		c3 := cc[3]
		if n == 2 {
			c2 = cc[2]
			c3 = 0
		}
		det := c2*c2 - c1*c3
		if det == 0 {
			u += 1.0
			v += 1.0
			continue
		}
		// The synthetic-division recurrences above are written for the divisor
		// x^2 + u*x + v, which corresponds to the textbook (r, s) form with
		// r = -u and s = -v. The Newton corrections (dr, ds) below are therefore
		// subtracted from u and v.
		dr := (-b[1]*c2 + b[0]*c3) / det
		ds := (-b[0]*c2 + b[1]*c1) / det
		du := -dr
		dv := -ds
		u += du
		v += dv
		if math.Abs(du) <= tol*(1+math.Abs(u)) && math.Abs(dv) <= tol*(1+math.Abs(v)) {
			return u, v, it, nil
		}
	}
	return u, v, maxIter, ErrNoConvergence
}

// BairstowRoots finds all roots of the real polynomial p by repeatedly applying
// [Bairstow] to peel off quadratic factors and deflating, dropping to a linear
// solve when a single degree remains. It returns every root (real roots have
// zero imaginary part) using only real arithmetic internally. Roots are polished
// against the original polynomial.
func BairstowRoots(p Poly, tol float64, maxIter int) ([]complex128, error) {
	if tol <= 0 {
		tol = 1e-13
	}
	maxIter = resolveMax(maxIter)
	work := p.Trim().Clone()
	d := work.Degree()
	if d < 1 {
		return nil, ErrDegreeTooLow
	}
	roots := make([]complex128, 0, d)
	// Seed guesses from coefficient ratios.
	u, v := 0.0, 0.0
	for work.Degree() > 2 {
		n := work.Degree()
		if u == 0 && v == 0 {
			u = work[n-1] / work[n]
			v = work[n-2] / work[n]
		}
		uu, vv, _, err := Bairstow(work, u, v, tol, maxIter)
		if err != nil {
			// Try a fresh random-ish seed.
			uu, vv, _, err = Bairstow(work, 1.0, 1.0, tol, maxIter*2)
			if err != nil {
				return roots, err
			}
		}
		r1, r2 := QuadraticRoots(1, uu, vv)
		roots = append(roots, r1, r2)
		// Deflate by x^2 + uu x + vv.
		work = deflateQuadratic(work, uu, vv)
		u, v = uu, vv
	}
	switch work.Degree() {
	case 2:
		r1, r2 := QuadraticRoots(work[2], work[1], work[0])
		roots = append(roots, r1, r2)
	case 1:
		roots = append(roots, complex(-work[0]/work[1], 0))
	}
	// Polish real-approx roots with Newton on original; leave complex as-is.
	cp := p.ToComplex()
	dcp := cp.Derivative()
	for i, r := range roots {
		for k := 0; k < 3; k++ {
			pv := cp.Eval(r)
			dv := dcp.Eval(r)
			if dv == 0 {
				break
			}
			step := pv / dv
			r -= step
			if cmplx.Abs(step) <= tol*(1+cmplx.Abs(r)) {
				break
			}
		}
		roots[i] = r
	}
	sortComplex(roots)
	return roots, nil
}

// deflateQuadratic divides the real polynomial p by the monic quadratic
// x^2 + u*x + v exactly (assuming it is a factor), returning the quotient.
func deflateQuadratic(p Poly, u, v float64) Poly {
	n := p.Degree()
	if n < 2 {
		return Poly{}
	}
	b := make([]float64, n+1)
	b[n] = p[n]
	b[n-1] = p[n-1] - u*b[n]
	for i := n - 2; i >= 0; i-- {
		b[i] = p[i] - u*b[i+1] - v*b[i+2]
	}
	// Quotient coefficients are b[2..n] shifted down by 2.
	q := make(Poly, n-1)
	for i := 2; i <= n; i++ {
		q[i-2] = b[i]
	}
	return q.Trim().Clone()
}
