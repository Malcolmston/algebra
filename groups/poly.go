package groups

import (
	"fmt"
	"math"
	"strings"
)

// Poly is a univariate polynomial over the real (rational) field, stored in
// coefficient slice form from lowest to highest degree: Poly{c0, c1, c2}
// represents c0 + c1·x + c2·x². The zero polynomial is represented by an empty
// (or all-zero) slice. Polynomials form a Euclidean domain, so the Euclidean
// algorithm yields a greatest common divisor via [PolyGCD].
type Poly []float64

// groupsPolyEps is the magnitude below which a coefficient is treated as zero
// when trimming and comparing polynomials.
const groupsPolyEps = 1e-9

// PolyTrim returns a copy of a with trailing (highest-degree) coefficients that
// are effectively zero removed, giving the canonical representation. The zero
// polynomial trims to an empty slice.
func PolyTrim(a Poly) Poly {
	n := len(a)
	for n > 0 && math.Abs(a[n-1]) <= groupsPolyEps {
		n--
	}
	return append(Poly(nil), a[:n]...)
}

// Degree returns the degree of a: the index of its highest non-zero
// coefficient. The zero polynomial has degree -1 by convention.
func (a Poly) Degree() int {
	t := PolyTrim(a)
	return len(t) - 1
}

// IsZero reports whether a is the zero polynomial (every coefficient within
// tolerance of zero).
func (a Poly) IsZero() bool {
	return a.Degree() < 0
}

// LeadingCoeff returns the coefficient of the highest-degree term of a, or 0
// for the zero polynomial.
func (a Poly) LeadingCoeff() float64 {
	t := PolyTrim(a)
	if len(t) == 0 {
		return 0
	}
	return t[len(t)-1]
}

// PolyAdd returns the sum a + b.
func PolyAdd(a, b Poly) Poly {
	n := len(a)
	if len(b) > n {
		n = len(b)
	}
	out := make(Poly, n)
	for i := 0; i < n; i++ {
		var av, bv float64
		if i < len(a) {
			av = a[i]
		}
		if i < len(b) {
			bv = b[i]
		}
		out[i] = av + bv
	}
	return PolyTrim(out)
}

// PolySub returns the difference a - b.
func PolySub(a, b Poly) Poly {
	return PolyAdd(a, PolyScale(b, -1))
}

// PolyScale returns the polynomial a with every coefficient multiplied by the
// scalar c.
func PolyScale(a Poly, c float64) Poly {
	out := make(Poly, len(a))
	for i, v := range a {
		out[i] = v * c
	}
	return PolyTrim(out)
}

// PolyMul returns the product a·b via the convolution of coefficient slices.
func PolyMul(a, b Poly) Poly {
	a, b = PolyTrim(a), PolyTrim(b)
	if len(a) == 0 || len(b) == 0 {
		return Poly{}
	}
	out := make(Poly, len(a)+len(b)-1)
	for i, av := range a {
		for j, bv := range b {
			out[i+j] += av * bv
		}
	}
	return PolyTrim(out)
}

// PolyDivMod returns the quotient q and remainder r of polynomial long
// division a = q·b + r with deg(r) < deg(b). The divisor b must be non-zero; it
// panics otherwise.
func PolyDivMod(a, b Poly) (q, r Poly) {
	b = PolyTrim(b)
	if len(b) == 0 {
		panic("groups: PolyDivMod division by zero polynomial")
	}
	rem := PolyTrim(a)
	db := len(b) - 1
	lcB := b[db]
	quo := make(Poly, 0)
	for len(rem)-1 >= db && len(rem) > 0 {
		dr := len(rem) - 1
		coeff := rem[dr] / lcB
		shift := dr - db
		// Ensure quotient slice is large enough.
		for len(quo) <= shift {
			quo = append(quo, 0)
		}
		quo[shift] = coeff
		// Subtract coeff·x^shift·b from rem.
		for i := 0; i <= db; i++ {
			rem[shift+i] -= coeff * b[i]
		}
		rem = PolyTrim(rem)
	}
	return PolyTrim(quo), PolyTrim(rem)
}

// PolyMod returns the remainder of a divided by b (the r from [PolyDivMod]).
func PolyMod(a, b Poly) Poly {
	_, r := PolyDivMod(a, b)
	return r
}

// PolyMonic returns a scaled to be monic (leading coefficient 1). The zero
// polynomial is returned unchanged.
func PolyMonic(a Poly) Poly {
	lc := a.LeadingCoeff()
	if lc == 0 {
		return Poly{}
	}
	return PolyScale(a, 1/lc)
}

// PolyGCD returns the monic greatest common divisor of a and b via the
// Euclidean algorithm over polynomials. PolyGCD of two zero polynomials is the
// zero polynomial.
func PolyGCD(a, b Poly) Poly {
	a, b = PolyTrim(a), PolyTrim(b)
	for len(b) > 0 {
		a, b = b, PolyMod(a, b)
	}
	if len(a) == 0 {
		return Poly{}
	}
	return PolyMonic(a)
}

// PolyEval returns the value a(x) evaluated at the point x using Horner's
// method.
func PolyEval(a Poly, x float64) float64 {
	result := 0.0
	for i := len(a) - 1; i >= 0; i-- {
		result = result*x + a[i]
	}
	return result
}

// PolyDerivative returns the formal derivative a'(x) of the polynomial a.
func PolyDerivative(a Poly) Poly {
	if len(a) <= 1 {
		return Poly{}
	}
	out := make(Poly, len(a)-1)
	for i := 1; i < len(a); i++ {
		out[i-1] = a[i] * float64(i)
	}
	return PolyTrim(out)
}

// PolyEqual reports whether a and b are equal as polynomials, comparing
// coefficients within tolerance after trimming.
func PolyEqual(a, b Poly) bool {
	ta, tb := PolyTrim(a), PolyTrim(b)
	if len(ta) != len(tb) {
		return false
	}
	for i := range ta {
		if math.Abs(ta[i]-tb[i]) > groupsPolyEps {
			return false
		}
	}
	return true
}

// String renders a in conventional descending-degree notation, e.g.
// "2x^2 + 3x - 1". The zero polynomial renders as "0".
func (a Poly) String() string {
	t := PolyTrim(a)
	if len(t) == 0 {
		return "0"
	}
	var parts []string
	for i := len(t) - 1; i >= 0; i-- {
		c := t[i]
		if math.Abs(c) <= groupsPolyEps {
			continue
		}
		mag := math.Abs(c)
		var term string
		switch {
		case i == 0:
			term = fmt.Sprintf("%g", mag)
		case i == 1 && math.Abs(mag-1) <= groupsPolyEps:
			term = "x"
		case i == 1:
			term = fmt.Sprintf("%gx", mag)
		case math.Abs(mag-1) <= groupsPolyEps:
			term = fmt.Sprintf("x^%d", i)
		default:
			term = fmt.Sprintf("%gx^%d", mag, i)
		}
		if len(parts) == 0 {
			if c < 0 {
				term = "-" + term
			}
		} else {
			if c < 0 {
				term = "- " + term
			} else {
				term = "+ " + term
			}
		}
		parts = append(parts, term)
	}
	if len(parts) == 0 {
		return "0"
	}
	return strings.Join(parts, " ")
}
