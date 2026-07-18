// Package powerseries implements truncated formal power series over float64
// coefficients together with the classical generating-function toolkit.
//
// A Series holds a finite slice of coefficients c[0], c[1], … , c[n-1] and is
// interpreted as the polynomial
//
//	c[0] + c[1]·x + c[2]·x² + … + c[n-1]·x^(n-1)
//
// with every coefficient of higher order taken to be exactly zero. The length
// of the slice is the working precision of the series. Arithmetic that would
// generate higher-order terms (multiplication, composition, the transcendental
// functions) truncates its result back to a finite precision so that the
// package works entirely with fixed-size representations.
//
// The package provides
//
//   - construction and inspection: New, FromSlice, Zero, One, Constant,
//     Monomial, Ident, FromFunc, together with Len, Coeff, Coeffs, Order,
//     Truncate, Extend, Evaluate and Equal;
//   - ring operations: Add, Sub, Neg, Scale, Mul, Hadamard, Shift and Pow;
//   - calculus: Derivative, Integral and IntegralConst;
//   - analytic operations: Inverse, Div, Compose, Reversion, Exp, Log, Sqrt,
//     PowReal, Sin, Cos, Tan, Sinh, Cosh and Atan;
//   - ordinary and exponential generating functions and the conversions
//     between them: OGFtoEGF, EGFtoOGF and friends;
//   - a library of named generating functions (geometric, exponential,
//     binomial, Catalan, Fibonacci, Bernoulli, Bell, derangement, harmonic,
//     Motzkin and integer-partition series) and the matching integer
//     sequences;
//   - the Lagrange–Bürmann inversion formula: LagrangeInversion and
//     LagrangeInversionApply.
//
// Every routine is deterministic and depends only on the Go standard library.
package powerseries

import (
	"math"
	"strconv"
	"strings"
)

// Series is a truncated formal power series with float64 coefficients. The
// zero value is not usable; construct series with New, Zero, One and the other
// constructors. A Series is treated as immutable by every method: operations
// return freshly allocated series and never mutate their receiver or
// arguments.
type Series struct {
	coeffs []float64
}

// New returns the series whose coefficients are the supplied values in order of
// increasing degree, so New(1, 2, 3) is 1 + 2x + 3x². Called with no arguments
// it returns the zero series of length one.
func New(coeffs ...float64) Series {
	if len(coeffs) == 0 {
		return Series{coeffs: []float64{0}}
	}
	out := make([]float64, len(coeffs))
	copy(out, coeffs)
	return Series{coeffs: out}
}

// FromSlice returns a series that copies the coefficients of the supplied
// slice. An empty or nil slice yields the zero series of length one.
func FromSlice(coeffs []float64) Series {
	if len(coeffs) == 0 {
		return Series{coeffs: []float64{0}}
	}
	out := make([]float64, len(coeffs))
	copy(out, coeffs)
	return Series{coeffs: out}
}

// Zero returns the zero series carrying n coefficients (precision n). It panics
// if n is not positive.
func Zero(n int) Series {
	if n <= 0 {
		panic("powerseries: precision must be positive")
	}
	return Series{coeffs: make([]float64, n)}
}

// One returns the constant series 1 carrying n coefficients. It panics if n is
// not positive.
func One(n int) Series {
	return Constant(1, n)
}

// Constant returns the constant series c carrying n coefficients. It panics if
// n is not positive.
func Constant(c float64, n int) Series {
	s := Zero(n)
	s.coeffs[0] = c
	return s
}

// Monomial returns the series c·x^degree carrying n coefficients. Terms whose
// degree is at least n are truncated away. It panics if n is not positive or if
// degree is negative.
func Monomial(c float64, degree, n int) Series {
	if degree < 0 {
		panic("powerseries: monomial degree must be non-negative")
	}
	s := Zero(n)
	if degree < n {
		s.coeffs[degree] = c
	}
	return s
}

// Ident returns the identity series x carrying n coefficients. It is a shorthand
// for Monomial(1, 1, n) and is the natural argument for Compose and Reversion.
func Ident(n int) Series {
	return Monomial(1, 1, n)
}

// FromFunc returns the length-n series whose degree-i coefficient is f(i).
func FromFunc(f func(i int) float64, n int) Series {
	s := Zero(n)
	for i := 0; i < n; i++ {
		s.coeffs[i] = f(i)
	}
	return s
}

// Len reports the number of stored coefficients, that is the working precision
// of the series.
func (s Series) Len() int { return len(s.coeffs) }

// Coeff returns the coefficient of x^i. Indices outside the stored range return
// zero, matching the convention that higher-order and negative-order terms are
// exactly zero.
func (s Series) Coeff(i int) float64 {
	if i < 0 || i >= len(s.coeffs) {
		return 0
	}
	return s.coeffs[i]
}

// Coeffs returns a fresh copy of the coefficient slice in order of increasing
// degree.
func (s Series) Coeffs() []float64 {
	out := make([]float64, len(s.coeffs))
	copy(out, s.coeffs)
	return out
}

// Clone returns an independent copy of the series.
func (s Series) Clone() Series {
	return FromSlice(s.coeffs)
}

// Order returns the valuation of the series, that is the index of its lowest
// non-zero coefficient. Coefficients whose magnitude does not exceed tol are
// treated as zero. The zero series (no coefficient above tol) reports the
// series length.
func (s Series) Order(tol float64) int {
	for i, c := range s.coeffs {
		if math.Abs(c) > tol {
			return i
		}
	}
	return len(s.coeffs)
}

// IsZero reports whether every coefficient has magnitude at most tol.
func (s Series) IsZero(tol float64) bool {
	return s.Order(tol) == len(s.coeffs)
}

// Truncate returns the series keeping only the coefficients of degree below n.
// If n is at least the current length the series is returned unchanged (as a
// copy). It panics if n is not positive.
func (s Series) Truncate(n int) Series {
	if n <= 0 {
		panic("powerseries: precision must be positive")
	}
	if n >= len(s.coeffs) {
		return s.Clone()
	}
	return FromSlice(s.coeffs[:n])
}

// Extend returns the series padded with zero coefficients up to length n. If n
// is at most the current length the series is returned unchanged (as a copy).
func (s Series) Extend(n int) Series {
	if n <= len(s.coeffs) {
		return s.Clone()
	}
	out := make([]float64, n)
	copy(out, s.coeffs)
	return Series{coeffs: out}
}

// Equal reports whether s and other agree, coefficient by coefficient, to within
// an absolute tolerance of tol. Series of different length are compared as
// though the shorter were padded with zeros.
func (s Series) Equal(other Series, tol float64) bool {
	n := len(s.coeffs)
	if len(other.coeffs) > n {
		n = len(other.coeffs)
	}
	for i := 0; i < n; i++ {
		if math.Abs(s.Coeff(i)-other.Coeff(i)) > tol {
			return false
		}
	}
	return true
}

// Evaluate returns the value of the truncated series at x using Horner's rule.
// Because the series is truncated the result approximates the underlying
// function only for x well inside the radius of convergence.
func (s Series) Evaluate(x float64) float64 {
	var acc float64
	for i := len(s.coeffs) - 1; i >= 0; i-- {
		acc = acc*x + s.coeffs[i]
	}
	return acc
}

// String renders the series in the conventional a + b·x + c·x² notation, ending
// with an O(x^n) remainder term recording the precision.
func (s Series) String() string {
	var b strings.Builder
	first := true
	for i, c := range s.coeffs {
		if c == 0 {
			continue
		}
		if !first {
			b.WriteString(" + ")
		}
		first = false
		b.WriteString(strconv.FormatFloat(c, 'g', -1, 64))
		switch i {
		case 0:
		case 1:
			b.WriteString("·x")
		default:
			b.WriteString("·x^")
			b.WriteString(strconv.Itoa(i))
		}
	}
	if first {
		b.WriteString("0")
	}
	b.WriteString(" + O(x^")
	b.WriteString(strconv.Itoa(len(s.coeffs)))
	b.WriteString(")")
	return b.String()
}

// powerseriesMaxLen returns the larger of two lengths.
func powerseriesMaxLen(a, b int) int {
	if a > b {
		return a
	}
	return b
}
