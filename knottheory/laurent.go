package knottheory

import (
	"fmt"
	"sort"
	"strings"
)

// Laurent is a Laurent polynomial in a single variable with integer
// coefficients, that is a finite sum of terms c*X^e where e may be any integer
// (positive, negative or zero). The zero polynomial is represented by an empty
// coefficient slice. A Laurent value is always kept trimmed so that neither the
// leading nor the trailing coefficient is zero.
type Laurent struct {
	min    int   // exponent of coeffs[0]
	coeffs []int // coeffs[i] is the coefficient of X^(min+i)
}

// LaurentTerm is a single monomial c*X^e of a Laurent polynomial.
type LaurentTerm struct {
	Exp   int
	Coeff int
}

// trim removes leading and trailing zero coefficients and returns a canonical
// Laurent value. It is used internally after every arithmetic operation.
func trim(min int, coeffs []int) Laurent {
	lo := 0
	for lo < len(coeffs) && coeffs[lo] == 0 {
		lo++
	}
	if lo == len(coeffs) {
		return Laurent{}
	}
	hi := len(coeffs) - 1
	for hi >= 0 && coeffs[hi] == 0 {
		hi--
	}
	out := make([]int, hi-lo+1)
	copy(out, coeffs[lo:hi+1])
	return Laurent{min: min + lo, coeffs: out}
}

// NewLaurent builds a Laurent polynomial whose coefficient of X^(min+i) is
// coeffs[i]. The slice is copied and the result is trimmed.
func NewLaurent(min int, coeffs []int) Laurent {
	cp := make([]int, len(coeffs))
	copy(cp, coeffs)
	return trim(min, cp)
}

// Monomial returns the single-term Laurent polynomial coeff*X^exp.
func Monomial(coeff, exp int) Laurent {
	if coeff == 0 {
		return Laurent{}
	}
	return Laurent{min: exp, coeffs: []int{coeff}}
}

// LaurentConst returns the constant Laurent polynomial c (that is c*X^0).
func LaurentConst(c int) Laurent { return Monomial(c, 0) }

// ZeroLaurent returns the zero Laurent polynomial.
func ZeroLaurent() Laurent { return Laurent{} }

// OneLaurent returns the constant Laurent polynomial 1.
func OneLaurent() Laurent { return Monomial(1, 0) }

// IsZero reports whether L is the zero polynomial.
func (L Laurent) IsZero() bool { return len(L.coeffs) == 0 }

// Coeff returns the coefficient of X^exp in L.
func (L Laurent) Coeff(exp int) int {
	i := exp - L.min
	if i < 0 || i >= len(L.coeffs) {
		return 0
	}
	return L.coeffs[i]
}

// MinDegree returns the smallest exponent that occurs with a non-zero
// coefficient. It returns 0 for the zero polynomial.
func (L Laurent) MinDegree() int {
	if L.IsZero() {
		return 0
	}
	return L.min
}

// MaxDegree returns the largest exponent that occurs with a non-zero
// coefficient. It returns 0 for the zero polynomial.
func (L Laurent) MaxDegree() int {
	if L.IsZero() {
		return 0
	}
	return L.min + len(L.coeffs) - 1
}

// SpanWidth returns MaxDegree-MinDegree, the breadth of the polynomial. The
// span of the Jones polynomial is a lower bound for the crossing number.
func (L Laurent) SpanWidth() int {
	if L.IsZero() {
		return 0
	}
	return len(L.coeffs) - 1
}

// NumTerms returns the number of non-zero terms of L.
func (L Laurent) NumTerms() int {
	n := 0
	for _, c := range L.coeffs {
		if c != 0 {
			n++
		}
	}
	return n
}

// LeadingCoeff returns the coefficient of the highest-degree term, or 0 for the
// zero polynomial.
func (L Laurent) LeadingCoeff() int {
	if L.IsZero() {
		return 0
	}
	return L.coeffs[len(L.coeffs)-1]
}

// TrailingCoeff returns the coefficient of the lowest-degree term, or 0 for the
// zero polynomial.
func (L Laurent) TrailingCoeff() int {
	if L.IsZero() {
		return 0
	}
	return L.coeffs[0]
}

// Terms returns the non-zero terms of L sorted by increasing exponent.
func (L Laurent) Terms() []LaurentTerm {
	out := make([]LaurentTerm, 0, len(L.coeffs))
	for i, c := range L.coeffs {
		if c != 0 {
			out = append(out, LaurentTerm{Exp: L.min + i, Coeff: c})
		}
	}
	return out
}

// Clone returns an independent copy of L.
func (L Laurent) Clone() Laurent {
	cp := make([]int, len(L.coeffs))
	copy(cp, L.coeffs)
	return Laurent{min: L.min, coeffs: cp}
}

// Neg returns -L.
func (L Laurent) Neg() Laurent {
	out := make([]int, len(L.coeffs))
	for i, c := range L.coeffs {
		out[i] = -c
	}
	return Laurent{min: L.min, coeffs: out}
}

// Add returns L+other.
func (L Laurent) Add(other Laurent) Laurent {
	if L.IsZero() {
		return other.Clone()
	}
	if other.IsZero() {
		return L.Clone()
	}
	lo := L.min
	if other.min < lo {
		lo = other.min
	}
	hiL := L.min + len(L.coeffs) - 1
	hiO := other.min + len(other.coeffs) - 1
	hi := hiL
	if hiO > hi {
		hi = hiO
	}
	out := make([]int, hi-lo+1)
	for i, c := range L.coeffs {
		out[L.min+i-lo] += c
	}
	for i, c := range other.coeffs {
		out[other.min+i-lo] += c
	}
	return trim(lo, out)
}

// Sub returns L-other.
func (L Laurent) Sub(other Laurent) Laurent { return L.Add(other.Neg()) }

// Scale returns k*L.
func (L Laurent) Scale(k int) Laurent {
	if k == 0 {
		return Laurent{}
	}
	out := make([]int, len(L.coeffs))
	for i, c := range L.coeffs {
		out[i] = c * k
	}
	return Laurent{min: L.min, coeffs: out}
}

// ShiftExp returns X^d * L, that is L with every exponent shifted by d.
func (L Laurent) ShiftExp(d int) Laurent {
	if L.IsZero() {
		return Laurent{}
	}
	return Laurent{min: L.min + d, coeffs: append([]int(nil), L.coeffs...)}
}

// Mul returns the product L*other.
func (L Laurent) Mul(other Laurent) Laurent {
	if L.IsZero() || other.IsZero() {
		return Laurent{}
	}
	out := make([]int, len(L.coeffs)+len(other.coeffs)-1)
	for i, a := range L.coeffs {
		if a == 0 {
			continue
		}
		for j, b := range other.coeffs {
			out[i+j] += a * b
		}
	}
	return trim(L.min+other.min, out)
}

// Pow returns L raised to the non-negative power n. Pow panics if n is
// negative because a general Laurent polynomial has no polynomial inverse; use
// a monomial and ShiftExp for negative powers of a variable.
func (L Laurent) Pow(n int) Laurent {
	if n < 0 {
		panic("knottheory: Laurent.Pow requires a non-negative exponent")
	}
	result := OneLaurent()
	base := L.Clone()
	for n > 0 {
		if n&1 == 1 {
			result = result.Mul(base)
		}
		base = base.Mul(base)
		n >>= 1
	}
	return result
}

// Equal reports whether L and other are the same polynomial.
func (L Laurent) Equal(other Laurent) bool {
	if L.min != other.min && !(L.IsZero() && other.IsZero()) {
		if len(L.coeffs) != len(other.coeffs) {
			return false
		}
	}
	if len(L.coeffs) != len(other.coeffs) {
		return false
	}
	if L.IsZero() {
		return true
	}
	if L.min != other.min {
		return false
	}
	for i := range L.coeffs {
		if L.coeffs[i] != other.coeffs[i] {
			return false
		}
	}
	return true
}

// Reverse returns L(X^{-1}), the polynomial obtained by replacing X with its
// inverse. A palindromic polynomial (many knot invariants are) is fixed by
// Reverse up to a shift.
func (L Laurent) Reverse() Laurent {
	if L.IsZero() {
		return Laurent{}
	}
	n := len(L.coeffs)
	out := make([]int, n)
	for i, c := range L.coeffs {
		out[n-1-i] = c
	}
	return Laurent{min: -(L.min + n - 1), coeffs: out}
}

// IsPalindromic reports whether L(X) equals L(X^{-1}) up to multiplication by a
// power of X. The Alexander polynomial of a knot is always palindromic.
func (L Laurent) IsPalindromic() bool {
	if L.IsZero() {
		return true
	}
	n := len(L.coeffs)
	for i := 0; i < n; i++ {
		if L.coeffs[i] != L.coeffs[n-1-i] {
			return false
		}
	}
	return true
}

// Content returns the greatest common divisor of the coefficients of L (the
// integer content), and 0 for the zero polynomial.
func (L Laurent) Content() int {
	g := 0
	for _, c := range L.coeffs {
		g = gcdInt(g, c)
	}
	return g
}

// Derivative returns the formal derivative dL/dX, itself a Laurent polynomial.
func (L Laurent) Derivative() Laurent {
	if L.IsZero() {
		return Laurent{}
	}
	out := make([]int, len(L.coeffs))
	for i, c := range L.coeffs {
		e := L.min + i
		out[i] = c * e
	}
	return trim(L.min-1, out)
}

// Eval evaluates L at the real number x. Negative exponents use 1/x; Eval
// panics only if x is zero and a negative exponent is present.
func (L Laurent) Eval(x float64) float64 {
	var sum float64
	for i, c := range L.coeffs {
		if c == 0 {
			continue
		}
		sum += float64(c) * ipow(x, L.min+i)
	}
	return sum
}

// EvalUnit evaluates L at an integer root of unity argument, meaning it returns
// the integer value L(x) for x equal to +1 or -1, where all powers stay
// integral. EvalUnit panics for any other x.
func (L Laurent) EvalUnit(x int) int {
	if x != 1 && x != -1 {
		panic("knottheory: EvalUnit only supports x = 1 or x = -1")
	}
	sum := 0
	for i, c := range L.coeffs {
		if c == 0 {
			continue
		}
		e := L.min + i
		term := c
		if x == -1 && ((e%2+2)%2) == 1 {
			term = -term
		}
		sum += term
	}
	return sum
}

// SubstitutePow returns L with the variable X replaced by X^k, mapping every
// term c*X^e to c*X^{e*k}. k must be non-zero.
func (L Laurent) SubstitutePow(k int) Laurent {
	if k == 0 {
		panic("knottheory: SubstitutePow requires a non-zero exponent")
	}
	if L.IsZero() {
		return Laurent{}
	}
	terms := L.Terms()
	res := ZeroLaurent()
	for _, t := range terms {
		res = res.Add(Monomial(t.Coeff, t.Exp*k))
	}
	return res
}

// String renders L in the variable X using conventional notation, for example
// "X^-2 - X + 3". The zero polynomial renders as "0".
func (L Laurent) String() string { return L.StringVar("X") }

// StringVar renders L using the supplied variable name.
func (L Laurent) StringVar(v string) string {
	if L.IsZero() {
		return "0"
	}
	terms := L.Terms()
	sort.Slice(terms, func(i, j int) bool { return terms[i].Exp > terms[j].Exp })
	var b strings.Builder
	for idx, t := range terms {
		c := t.Coeff
		sign := "+"
		if c < 0 {
			sign = "-"
			c = -c
		}
		if idx == 0 {
			if sign == "-" {
				b.WriteString("-")
			}
		} else {
			b.WriteString(" ")
			b.WriteString(sign)
			b.WriteString(" ")
		}
		mono := monomialString(v, t.Exp)
		if mono == "" {
			fmt.Fprintf(&b, "%d", c)
		} else if c == 1 {
			b.WriteString(mono)
		} else {
			fmt.Fprintf(&b, "%d%s", c, mono)
		}
	}
	return b.String()
}

// monomialString renders v^e, returning "" for e==0, v for e==1 and v^e
// otherwise.
func monomialString(v string, e int) string {
	switch e {
	case 0:
		return ""
	case 1:
		return v
	default:
		return fmt.Sprintf("%s^%d", v, e)
	}
}

// ipow raises x to the (possibly negative) integer power e.
func ipow(x float64, e int) float64 {
	if e == 0 {
		return 1
	}
	neg := e < 0
	if neg {
		e = -e
	}
	r := 1.0
	for i := 0; i < e; i++ {
		r *= x
	}
	if neg {
		return 1 / r
	}
	return r
}

// gcdInt returns the non-negative greatest common divisor of a and b.
func gcdInt(a, b int) int {
	if a < 0 {
		a = -a
	}
	if b < 0 {
		b = -b
	}
	for b != 0 {
		a, b = b, a%b
	}
	return a
}

// LaurentSum returns the sum of all the supplied Laurent polynomials.
func LaurentSum(ps ...Laurent) Laurent {
	acc := ZeroLaurent()
	for _, p := range ps {
		acc = acc.Add(p)
	}
	return acc
}

// LaurentProduct returns the product of all the supplied Laurent polynomials.
func LaurentProduct(ps ...Laurent) Laurent {
	acc := OneLaurent()
	for _, p := range ps {
		acc = acc.Mul(p)
	}
	return acc
}

// DivExact returns the quotient L/other when the division is exact, and reports
// ok=false otherwise. Division is performed in the Laurent ring: the result is
// the unique Laurent polynomial q with q*other == L when one exists.
func (L Laurent) DivExact(other Laurent) (q Laurent, ok bool) {
	if other.IsZero() {
		return Laurent{}, false
	}
	if L.IsZero() {
		return Laurent{}, true
	}
	// Shift both to non-negative exponents to perform ordinary polynomial
	// long division, then shift the quotient back.
	shift := L.min - other.min
	a := append([]int(nil), L.coeffs...)
	bcoef := append([]int(nil), other.coeffs...)
	if len(a) < len(bcoef) {
		return Laurent{}, false
	}
	quo := make([]int, len(a)-len(bcoef)+1)
	rem := append([]int(nil), a...)
	bLead := bcoef[len(bcoef)-1]
	for i := len(quo) - 1; i >= 0; i-- {
		hi := i + len(bcoef) - 1
		if rem[hi] == 0 {
			continue
		}
		if rem[hi]%bLead != 0 {
			return Laurent{}, false
		}
		factor := rem[hi] / bLead
		quo[i] = factor
		for j := 0; j < len(bcoef); j++ {
			rem[i+j] -= factor * bcoef[j]
		}
	}
	for _, r := range rem {
		if r != 0 {
			return Laurent{}, false
		}
	}
	return trim(shift, quo), true
}
