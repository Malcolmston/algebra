package rootfind

import (
	"math"
	"math/cmplx"
	"sort"
	"strconv"
	"strings"
)

// At is an alias for Eval, evaluating p(x).
func (p Poly) At(x float64) float64 { return p.Eval(x) }

// At is an alias for Eval, evaluating c(x).
func (c CPoly) At(x complex128) complex128 { return c.Eval(x) }

// IsConstant reports whether p has degree 0 or is the zero polynomial.
func (p Poly) IsConstant() bool { return p.Degree() <= 0 }

// IsMonic reports whether the leading coefficient of p equals 1 within tol.
func (p Poly) IsMonic(tol float64) bool {
	d := p.Degree()
	return d >= 0 && math.Abs(p[d]-1) <= tol
}

// NumTerms returns the number of nonzero coefficients of p.
func (p Poly) NumTerms() int {
	n := 0
	for _, c := range p {
		if c != 0 {
			n++
		}
	}
	return n
}

// InfNorm returns the largest absolute coefficient of p (the sup norm of its
// coefficient vector).
func (p Poly) InfNorm() float64 { return polyMaxAbs(p) }

// L1Norm returns the sum of absolute values of the coefficients of p.
func (p Poly) L1Norm() float64 {
	s := 0.0
	for _, c := range p {
		s += math.Abs(c)
	}
	return s
}

// L2Norm returns the Euclidean norm of the coefficient vector of p.
func (p Poly) L2Norm() float64 {
	s := 0.0
	for _, c := range p {
		s += c * c
	}
	return math.Sqrt(s)
}

// Pow returns p raised to the nonnegative integer power k by repeated squaring.
// Pow(0) is the constant polynomial 1.
func (p Poly) Pow(k int) Poly {
	if k <= 0 {
		return Poly{1}
	}
	result := Poly{1}
	base := p.Clone()
	for k > 0 {
		if k&1 == 1 {
			result = result.Mul(base)
		}
		k >>= 1
		if k > 0 {
			base = base.Mul(base)
		}
	}
	return result
}

// SumOfRoots returns the sum of all roots of p counted with multiplicity, which
// by Vieta's formulas equals -a_{n-1}/a_n. It returns 0 for constant
// polynomials.
func (p Poly) SumOfRoots() float64 {
	d := p.Degree()
	if d < 1 {
		return 0
	}
	return -p[d-1] / p[d]
}

// ProductOfRoots returns the product of all roots of p counted with
// multiplicity, which by Vieta's formulas equals (-1)^n * a_0/a_n.
func (p Poly) ProductOfRoots() float64 {
	d := p.Degree()
	if d < 1 {
		return 0
	}
	prod := p[0] / p[d]
	if d%2 == 1 {
		prod = -prod
	}
	return prod
}

// IsMonic reports whether the leading coefficient of c equals 1 within tol.
func (c CPoly) IsMonic(tol float64) bool {
	d := c.Degree()
	return d >= 0 && cmplx.Abs(c[d]-1) <= tol
}

// NumTerms returns the number of nonzero coefficients of c.
func (c CPoly) NumTerms() int {
	n := 0
	for _, v := range c {
		if v != 0 {
			n++
		}
	}
	return n
}

// Pow returns c raised to the nonnegative integer power k by repeated squaring.
func (c CPoly) Pow(k int) CPoly {
	if k <= 0 {
		return CPoly{1}
	}
	result := CPoly{1}
	base := c.Clone()
	for k > 0 {
		if k&1 == 1 {
			result = result.Mul(base)
		}
		k >>= 1
		if k > 0 {
			base = base.Mul(base)
		}
	}
	return result
}

// String renders c in descending-power notation with parenthesized complex
// coefficients, for example "(1+2i)x^2 + (3+0i)".
func (c CPoly) String() string {
	d := c.Degree()
	if d < 0 {
		return "0"
	}
	var b strings.Builder
	first := true
	for i := d; i >= 0; i-- {
		v := c[i]
		if v == 0 {
			continue
		}
		if !first {
			b.WriteString(" + ")
		}
		first = false
		b.WriteByte('(')
		b.WriteString(strconv.FormatFloat(real(v), 'g', -1, 64))
		if imag(v) >= 0 {
			b.WriteByte('+')
		}
		b.WriteString(strconv.FormatFloat(imag(v), 'g', -1, 64))
		b.WriteString("i)")
		switch i {
		case 0:
		case 1:
			b.WriteString("x")
		default:
			b.WriteString("x^")
			b.WriteString(strconv.Itoa(i))
		}
	}
	if first {
		return "0"
	}
	return b.String()
}

// Horner evaluates the polynomial given by ascending-order real coefficients at
// x using Horner's method, without constructing a Poly.
func Horner(coeffs []float64, x float64) float64 {
	if len(coeffs) == 0 {
		return 0
	}
	y := coeffs[len(coeffs)-1]
	for i := len(coeffs) - 2; i >= 0; i-- {
		y = y*x + coeffs[i]
	}
	return y
}

// HornerComplex evaluates the polynomial given by ascending-order complex
// coefficients at x using Horner's method.
func HornerComplex(coeffs []complex128, x complex128) complex128 {
	if len(coeffs) == 0 {
		return 0
	}
	y := coeffs[len(coeffs)-1]
	for i := len(coeffs) - 2; i >= 0; i-- {
		y = y*x + coeffs[i]
	}
	return y
}

// NewtonComplex finds a root of the complex polynomial c from the start x0 using
// Newton's method in complex arithmetic. It is the complex analogue of [Newton]
// and converges quadratically near a simple root.
func NewtonComplex(c CPoly, x0 complex128, tol float64, maxIter int) (complex128, int, error) {
	if tol <= 0 {
		tol = DefaultRootTol
	}
	maxIter = resolveMax(maxIter)
	x := x0
	for i := 1; i <= maxIter; i++ {
		v, d := c.EvalDeriv(x)
		if cmplx.Abs(v) <= tol {
			return x, i, nil
		}
		if d == 0 {
			return x, i, ErrZeroDerivative
		}
		step := v / d
		x -= step
		if cmplx.Abs(step) <= tol*(1+cmplx.Abs(x)) {
			return x, i, nil
		}
	}
	return x, maxIter, ErrNoConvergence
}

// HalleyComplex finds a root of the complex polynomial c from x0 using Halley's
// cubically convergent method in complex arithmetic.
func HalleyComplex(c CPoly, x0 complex128, tol float64, maxIter int) (complex128, int, error) {
	if tol <= 0 {
		tol = DefaultRootTol
	}
	maxIter = resolveMax(maxIter)
	x := x0
	for i := 1; i <= maxIter; i++ {
		v, d1, d2 := c.EvalDeriv2(x)
		if cmplx.Abs(v) <= tol {
			return x, i, nil
		}
		denom := 2*d1*d1 - v*d2
		if denom == 0 {
			return x, i, ErrZeroDerivative
		}
		step := 2 * v * d1 / denom
		x -= step
		if cmplx.Abs(step) <= tol*(1+cmplx.Abs(x)) {
			return x, i, nil
		}
	}
	return x, maxIter, ErrNoConvergence
}

// PolishComplexRoots refines each approximate root in roots with a few Newton
// steps on the complex polynomial c, improving accuracy after a global solve or
// deflation. The input slice is not modified; a refined copy is returned.
func PolishComplexRoots(c CPoly, roots []complex128, tol float64, steps int) []complex128 {
	if steps <= 0 {
		steps = 5
	}
	out := make([]complex128, len(roots))
	for i, r := range roots {
		x := r
		for k := 0; k < steps; k++ {
			v, d := c.EvalDeriv(x)
			if d == 0 {
				break
			}
			step := v / d
			x -= step
			if cmplx.Abs(step) <= tol*(1+cmplx.Abs(x)) {
				break
			}
		}
		out[i] = x
	}
	return out
}

// CountComplexRoots returns the number of roots of c counted with multiplicity,
// which by the fundamental theorem of algebra equals its degree.
func CountComplexRoots(c CPoly) int {
	d := c.Degree()
	if d < 0 {
		return 0
	}
	return d
}

// SeparateRoots partitions a list of complex roots into those that are
// effectively real (|imag| <= imagTol) and those that are genuinely complex.
// The real parts of the real roots are returned sorted in ascending order.
func SeparateRoots(roots []complex128, imagTol float64) (reals []float64, complexes []complex128) {
	if imagTol <= 0 {
		imagTol = 1e-9
	}
	for _, z := range roots {
		if math.Abs(imag(z)) <= imagTol*(1+math.Abs(real(z))) {
			reals = append(reals, real(z))
		} else {
			complexes = append(complexes, z)
		}
	}
	sort.Float64s(reals)
	sortComplex(complexes)
	return reals, complexes
}

// SolveLinear returns the root of the linear equation a*x + b = 0. It returns
// ErrDegreeTooLow when a is zero.
func SolveLinear(a, b float64) (float64, error) {
	if a == 0 {
		return 0, ErrDegreeTooLow
	}
	return -b / a, nil
}

// SolveQuadraticReal returns the real roots of a*x^2 + b*x + c in ascending
// order: two values when the discriminant is positive, one (repeated) when it is
// zero, and none when it is negative.
func SolveQuadraticReal(a, b, c float64) []float64 {
	if a == 0 {
		if b == 0 {
			return nil
		}
		return []float64{-c / b}
	}
	disc := b*b - 4*a*c
	if disc < 0 {
		return nil
	}
	if disc == 0 {
		return []float64{-b / (2 * a)}
	}
	r1, r2 := QuadraticRoots(a, b, c)
	out := []float64{real(r1), real(r2)}
	sort.Float64s(out)
	return out
}

// DiscriminantQuadratic returns the discriminant b^2 - 4ac of a*x^2 + b*x + c.
// Its sign classifies the roots: positive means two real roots, zero a repeated
// real root, negative a complex-conjugate pair.
func DiscriminantQuadratic(a, b, c float64) float64 {
	return b*b - 4*a*c
}

// DiscriminantCubic returns the discriminant of the cubic a*x^3+b*x^2+c*x+d,
// 18abcd - 4b^3d + b^2c^2 - 4ac^3 - 27a^2d^2. A positive value means three
// distinct real roots, zero a repeated root, and negative one real and two
// complex-conjugate roots.
func DiscriminantCubic(a, b, c, d float64) float64 {
	return 18*a*b*c*d - 4*b*b*b*d + b*b*c*c - 4*a*c*c*c - 27*a*a*d*d
}

// SolveCubic returns the three roots of a*x^3 + b*x^2 + c*x + d as complex
// numbers, using Cardano's method with complex cube roots so the formula is
// numerically valid in all cases. Real roots have (near) zero imaginary part.
// It returns ErrDegreeTooLow when a is zero.
func SolveCubic(a, b, c, d float64) ([]complex128, error) {
	if a == 0 {
		return nil, ErrDegreeTooLow
	}
	// Depressed cubic t^3 + p t + q via x = t - b/(3a).
	bb := b / a
	cc := c / a
	dd := d / a
	p := cc - bb*bb/3
	q := 2*bb*bb*bb/27 - bb*cc/3 + dd
	shift := complex(-bb/3, 0)
	pc := complex(p, 0)
	qc := complex(q, 0)
	disc := cmplx.Sqrt(qc*qc/4 + pc*pc*pc/27)
	u3 := -qc/2 + disc
	if u3 == 0 {
		u3 = -qc/2 - disc
	}
	u := cmplx.Pow(u3, complex(1.0/3.0, 0))
	omega := complex(-0.5, math.Sqrt(3)/2)
	roots := make([]complex128, 3)
	uk := u
	for k := 0; k < 3; k++ {
		var t complex128
		if uk == 0 {
			t = 0
		} else {
			t = uk - pc/(3*uk)
		}
		roots[k] = t + shift
		uk *= omega
	}
	sortComplex(roots)
	return roots, nil
}

// SolveCubicReal returns the real roots of the cubic a*x^3+b*x^2+c*x+d in
// ascending order, extracted from [SolveCubic] by keeping roots whose imaginary
// part is negligible.
func SolveCubicReal(a, b, c, d float64) ([]float64, error) {
	roots, err := SolveCubic(a, b, c, d)
	if err != nil {
		return nil, err
	}
	var out []float64
	for _, z := range roots {
		if math.Abs(imag(z)) <= 1e-9*(1+math.Abs(real(z))) {
			out = append(out, real(z))
		}
	}
	sort.Float64s(out)
	return out, nil
}

// TotalRealRoots returns the number of real roots of p counted with
// multiplicity, summing the multiplicities reported by
// [RealRootsWithMultiplicity].
func TotalRealRoots(p Poly, tol float64) int {
	n := 0
	for _, rm := range RealRootsWithMultiplicity(p, tol) {
		n += rm.Multiplicity
	}
	return n
}

// AllRootsReal reports whether every root of p is real, i.e. the count of
// distinct real roots accounts for the full degree once multiplicity is
// included.
func AllRootsReal(p Poly, tol float64) bool {
	return TotalRealRoots(p, tol) == p.Degree()
}
