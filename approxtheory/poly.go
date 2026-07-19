package approxtheory

import (
	"errors"
	"fmt"
	"math"
	"strings"
)

// ErrSingular is returned when a linear system encountered during a fit or
// solve is singular (or numerically so) and cannot be solved.
var ErrSingular = errors.New("approxtheory: singular linear system")

// ErrDimensionMismatch is returned when input slices that must have equal
// length do not.
var ErrDimensionMismatch = errors.New("approxtheory: dimension mismatch")

// ErrEmptyInput is returned when a routine that requires at least one point or
// coefficient is given none.
var ErrEmptyInput = errors.New("approxtheory: empty input")

// Polyval evaluates the polynomial with the given monomial coefficients
// (ascending order, coeffs[i] multiplies x**i) at x using Horner's method.
func Polyval(coeffs []float64, x float64) float64 {
	if len(coeffs) == 0 {
		return 0
	}
	acc := coeffs[len(coeffs)-1]
	for i := len(coeffs) - 2; i >= 0; i-- {
		acc = acc*x + coeffs[i]
	}
	return acc
}

// Horner is an alias for Polyval; it evaluates a monomial polynomial with
// Horner's method.
func Horner(coeffs []float64, x float64) float64 {
	return Polyval(coeffs, x)
}

// PolyvalSlice evaluates a monomial polynomial at every point in xs and
// returns the resulting slice.
func PolyvalSlice(coeffs []float64, xs []float64) []float64 {
	out := make([]float64, len(xs))
	for i, x := range xs {
		out[i] = Polyval(coeffs, x)
	}
	return out
}

// PolyDegree returns the degree of a monomial polynomial, ignoring trailing
// (high-order) coefficients that are exactly zero. The zero polynomial has
// degree 0.
func PolyDegree(coeffs []float64) int {
	d := 0
	for i := len(coeffs) - 1; i >= 0; i-- {
		if coeffs[i] != 0 {
			d = i
			break
		}
	}
	return d
}

// PolyLeadingCoeff returns the coefficient of the highest-order nonzero term,
// or 0 for the zero polynomial.
func PolyLeadingCoeff(coeffs []float64) float64 {
	for i := len(coeffs) - 1; i >= 0; i-- {
		if coeffs[i] != 0 {
			return coeffs[i]
		}
	}
	return 0
}

// PolyIsZero reports whether every coefficient is exactly zero.
func PolyIsZero(coeffs []float64) bool {
	for _, c := range coeffs {
		if c != 0 {
			return false
		}
	}
	return true
}

// PolyAdd returns the sum of two monomial polynomials. The result has the
// length of the longer input.
func PolyAdd(a, b []float64) []float64 {
	n := len(a)
	if len(b) > n {
		n = len(b)
	}
	out := make([]float64, n)
	for i := range a {
		out[i] += a[i]
	}
	for i := range b {
		out[i] += b[i]
	}
	return out
}

// PolySub returns the difference a-b of two monomial polynomials.
func PolySub(a, b []float64) []float64 {
	n := len(a)
	if len(b) > n {
		n = len(b)
	}
	out := make([]float64, n)
	for i := range a {
		out[i] += a[i]
	}
	for i := range b {
		out[i] -= b[i]
	}
	return out
}

// PolyScale returns the polynomial a scaled by the constant s.
func PolyScale(a []float64, s float64) []float64 {
	out := make([]float64, len(a))
	for i := range a {
		out[i] = a[i] * s
	}
	return out
}

// PolyMul returns the product of two monomial polynomials by convolution.
func PolyMul(a, b []float64) []float64 {
	if len(a) == 0 || len(b) == 0 {
		return nil
	}
	out := make([]float64, len(a)+len(b)-1)
	for i, ai := range a {
		if ai == 0 {
			continue
		}
		for j, bj := range b {
			out[i+j] += ai * bj
		}
	}
	return out
}

// PolyDeriv returns the derivative of a monomial polynomial.
func PolyDeriv(coeffs []float64) []float64 {
	if len(coeffs) <= 1 {
		return []float64{0}
	}
	out := make([]float64, len(coeffs)-1)
	for i := 1; i < len(coeffs); i++ {
		out[i-1] = coeffs[i] * float64(i)
	}
	return out
}

// PolyInt returns an antiderivative of a monomial polynomial whose constant
// term equals the supplied value.
func PolyInt(coeffs []float64, constant float64) []float64 {
	out := make([]float64, len(coeffs)+1)
	out[0] = constant
	for i := range coeffs {
		out[i+1] = coeffs[i] / float64(i+1)
	}
	return out
}

// PolyFromRoots returns the monic monomial polynomial whose roots are exactly
// the supplied values, i.e. product over i of (x - roots[i]).
func PolyFromRoots(roots []float64) []float64 {
	out := []float64{1}
	for _, r := range roots {
		out = PolyMul(out, []float64{-r, 1})
	}
	return out
}

// PolyNormalize returns a copy of coeffs with trailing zero coefficients
// removed; the zero polynomial is returned as a single zero coefficient.
func PolyNormalize(coeffs []float64) []float64 {
	d := len(coeffs)
	for d > 1 && coeffs[d-1] == 0 {
		d--
	}
	out := make([]float64, d)
	copy(out, coeffs)
	return out
}

// PolyEqual reports whether two monomial polynomials are equal within the
// absolute tolerance tol, ignoring trailing zero coefficients.
func PolyEqual(a, b []float64, tol float64) bool {
	n := len(a)
	if len(b) > n {
		n = len(b)
	}
	for i := 0; i < n; i++ {
		var av, bv float64
		if i < len(a) {
			av = a[i]
		}
		if i < len(b) {
			bv = b[i]
		}
		if math.Abs(av-bv) > tol {
			return false
		}
	}
	return true
}

// PolyString renders a monomial polynomial as a human readable expression in
// the variable x, e.g. "1 + 2*x - 3*x^2".
func PolyString(coeffs []float64) string {
	var b strings.Builder
	first := true
	for i, c := range coeffs {
		if c == 0 {
			continue
		}
		mag := c
		sign := "+"
		if c < 0 {
			sign = "-"
			mag = -c
		}
		if first {
			if sign == "-" {
				b.WriteString("-")
			}
			first = false
		} else {
			b.WriteString(" ")
			b.WriteString(sign)
			b.WriteString(" ")
		}
		switch i {
		case 0:
			fmt.Fprintf(&b, "%g", mag)
		case 1:
			fmt.Fprintf(&b, "%g*x", mag)
		default:
			fmt.Fprintf(&b, "%g*x^%d", mag, i)
		}
	}
	if first {
		return "0"
	}
	return b.String()
}

// Linspace returns n equally spaced points spanning the closed interval
// [a, b]. For n == 1 it returns the midpoint.
func Linspace(a, b float64, n int) []float64 {
	if n <= 0 {
		return nil
	}
	if n == 1 {
		return []float64{(a + b) / 2}
	}
	out := make([]float64, n)
	step := (b - a) / float64(n-1)
	for i := 0; i < n; i++ {
		out[i] = a + float64(i)*step
	}
	out[n-1] = b
	return out
}

// solveLinear solves the dense linear system A x = b in place using Gaussian
// elimination with partial pivoting. A is an n-by-n matrix stored row major as
// a slice of rows; b has length n. It returns the solution vector or
// ErrSingular.
func solveLinear(A [][]float64, b []float64) ([]float64, error) {
	n := len(A)
	if n == 0 {
		return nil, ErrEmptyInput
	}
	// Work on copies so callers keep their data.
	m := make([][]float64, n)
	for i := range A {
		if len(A[i]) != n {
			return nil, ErrDimensionMismatch
		}
		row := make([]float64, n+1)
		copy(row, A[i])
		row[n] = b[i]
		m[i] = row
	}
	for col := 0; col < n; col++ {
		// Partial pivot.
		piv := col
		best := math.Abs(m[col][col])
		for r := col + 1; r < n; r++ {
			if v := math.Abs(m[r][col]); v > best {
				best = v
				piv = r
			}
		}
		if best == 0 {
			return nil, ErrSingular
		}
		m[col], m[piv] = m[piv], m[col]
		pivot := m[col][col]
		for r := 0; r < n; r++ {
			if r == col {
				continue
			}
			factor := m[r][col] / pivot
			if factor == 0 {
				continue
			}
			for c := col; c <= n; c++ {
				m[r][c] -= factor * m[col][c]
			}
		}
	}
	x := make([]float64, n)
	for i := 0; i < n; i++ {
		x[i] = m[i][n] / m[i][i]
	}
	return x, nil
}
