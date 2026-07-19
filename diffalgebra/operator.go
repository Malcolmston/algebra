package diffalgebra

import (
	"math/big"
	"strings"
)

// Operator is a linear differential operator sum_i a_i(x) D^i with
// rational-function coefficients, an element of the non-commutative Ore ring
// Q(x)[D] with the multiplication rule D*a = a*D + a'. Coefficient i multiplies
// D^i; the slice is kept trimmed of trailing zero coefficients.
type Operator struct {
	c []RatFunc
}

// normalizeOp trims trailing zero coefficients.
func normalizeOp(c []RatFunc) Operator {
	n := len(c)
	for n > 0 && c[n-1].IsZero() {
		n--
	}
	return Operator{c: c[:n]}
}

// NewOperator builds an operator from coefficients in ascending order of D
// (coeffs[0] is the zeroth-order term).
func NewOperator(coeffs ...RatFunc) Operator {
	c := make([]RatFunc, len(coeffs))
	copy(c, coeffs)
	return normalizeOp(c)
}

// OperatorFromPolys builds an operator whose coefficients are the given
// polynomials, in ascending order of D.
func OperatorFromPolys(coeffs ...Poly) Operator {
	c := make([]RatFunc, len(coeffs))
	for i, p := range coeffs {
		c[i] = RatFuncFromPoly(p)
	}
	return normalizeOp(c)
}

// ZeroOperator returns the zero operator.
func ZeroOperator() Operator { return Operator{} }

// IdentityOperator returns the identity operator (multiplication by 1).
func IdentityOperator() Operator { return NewOperator(OneRatFunc()) }

// DOperator returns the derivation operator D itself.
func DOperator() Operator { return NewOperator(ZeroRatFunc(), OneRatFunc()) }

// ConstOperator returns the zeroth-order operator that multiplies by the
// rational function a.
func ConstOperator(a RatFunc) Operator { return NewOperator(a) }

// Order returns the order (highest power of D) of the operator; the zero
// operator has order -1.
func (o Operator) Order() int { return len(o.c) - 1 }

// IsZero reports whether o is the zero operator.
func (o Operator) IsZero() bool { return len(o.c) == 0 }

// Coeff returns the coefficient of D^i, or zero when out of range.
func (o Operator) Coeff(i int) RatFunc {
	if i < 0 || i >= len(o.c) {
		return ZeroRatFunc()
	}
	return o.c[i]
}

// Coeffs returns a copy of the coefficient slice in ascending order of D.
func (o Operator) Coeffs() []RatFunc {
	out := make([]RatFunc, len(o.c))
	copy(out, o.c)
	return out
}

// LeadingCoeff returns the coefficient of the highest-order term.
func (o Operator) LeadingCoeff() RatFunc {
	if o.IsZero() {
		return ZeroRatFunc()
	}
	return o.c[len(o.c)-1]
}

// Equal reports whether o and p are equal operators.
func (o Operator) Equal(p Operator) bool {
	if len(o.c) != len(p.c) {
		return false
	}
	for i := range o.c {
		if !o.c[i].Equal(p.c[i]) {
			return false
		}
	}
	return true
}

// Neg returns -o.
func (o Operator) Neg() Operator {
	c := make([]RatFunc, len(o.c))
	for i := range o.c {
		c[i] = o.c[i].Neg()
	}
	return normalizeOp(c)
}

// Add returns o+p.
func (o Operator) Add(p Operator) Operator {
	n := len(o.c)
	if len(p.c) > n {
		n = len(p.c)
	}
	c := make([]RatFunc, n)
	for i := 0; i < n; i++ {
		c[i] = o.Coeff(i).Add(p.Coeff(i))
	}
	return normalizeOp(c)
}

// Sub returns o-p.
func (o Operator) Sub(p Operator) Operator { return o.Add(p.Neg()) }

// ScalarMul multiplies o on the left by the rational function a (a*o).
func (o Operator) ScalarMul(a RatFunc) Operator {
	c := make([]RatFunc, len(o.c))
	for i := range o.c {
		c[i] = a.Mul(o.c[i])
	}
	return normalizeOp(c)
}

// Mul returns the composition o∘p in the Ore ring, using D*a = a*D + a'.
func (o Operator) Mul(p Operator) Operator {
	if o.IsZero() || p.IsZero() {
		return ZeroOperator()
	}
	acc := ZeroOperator()
	// o = sum_i a_i D^i ; multiply each a_i D^i by p on the right.
	for i := len(o.c) - 1; i >= 0; i-- {
		if o.c[i].IsZero() {
			continue
		}
		term := p
		// apply D^i on the left of p: D∘q handled by dApply.
		for k := 0; k < i; k++ {
			term = dApply(term)
		}
		term = term.ScalarMul(o.c[i])
		acc = acc.Add(term)
	}
	return acc
}

// dApply returns D∘q for an operator q, i.e. the operator whose action is the
// derivative of q applied to a function. Using D∘(sum b_j D^j) = sum (b_j' D^j
// + b_j D^{j+1}).
func dApply(q Operator) Operator {
	c := make([]RatFunc, len(q.c)+1)
	for i := range c {
		c[i] = ZeroRatFunc()
	}
	for j := 0; j < len(q.c); j++ {
		c[j] = c[j].Add(q.c[j].Derivative())
		c[j+1] = c[j+1].Add(q.c[j])
	}
	return normalizeOp(c)
}

// Pow raises the operator to the non-negative integer power n (composition).
func (o Operator) Pow(n int) Operator {
	result := IdentityOperator()
	base := o
	for n > 0 {
		if n&1 == 1 {
			result = result.Mul(base)
		}
		base = base.Mul(base)
		n >>= 1
	}
	return result
}

// ApplyRatFunc applies the operator to a rational function y, returning
// sum_i a_i y^(i).
func (o Operator) ApplyRatFunc(y RatFunc) RatFunc {
	acc := ZeroRatFunc()
	deriv := y
	for i := 0; i < len(o.c); i++ {
		acc = acc.Add(o.c[i].Mul(deriv))
		deriv = deriv.Derivative()
	}
	return acc
}

// ApplyPoly applies the operator to a polynomial y.
func (o Operator) ApplyPoly(y Poly) RatFunc { return o.ApplyRatFunc(RatFuncFromPoly(y)) }

// SymbolPoly returns the symbol of the operator: the polynomial in a formal
// variable obtained by replacing D^i with the i-th power and freezing the
// coefficients at the rational point x. It is used to read off the principal
// part; here it returns the leading coefficient's numerator scaled form as a
// Poly in D with the coefficients evaluated is not well defined over Q(x), so
// instead SymbolPoly returns the vector of coefficients as a Poly when the
// operator has polynomial coefficients.
func (o Operator) SymbolPoly() (Poly, bool) {
	c := make([]*big.Rat, len(o.c))
	for i := range o.c {
		if !o.c[i].IsPolynomial() || !o.c[i].Num().IsConstant() {
			return ZeroPoly(), false
		}
		c[i] = o.c[i].Num().ConstantTerm()
	}
	return normalizePoly(c), true
}

// IndicialPolynomial returns the indicial polynomial at the ordinary point x=0
// for an operator with polynomial coefficients, computed from the lowest-order
// behaviour. It returns the polynomial whose roots are the indicial exponents
// and reports false when the operator does not have the required form.
func (o Operator) IndicialPolynomial() (Poly, bool) {
	if o.IsZero() {
		return ZeroPoly(), false
	}
	n := o.Order()
	// indicial polynomial coefficient of s: sum over terms x^k D^i acting on
	// x^s contributes falling factorial. We build it for the lowest power of x.
	// For each coefficient a_i(x) = sum c_{i,k} x^k, term c_{i,k} x^k D^i x^s =
	// c_{i,k} (s)(s-1)...(s-i+1) x^{s-i+k}. The indicial polynomial collects the
	// coefficient of the minimal exponent of x.
	minExp := 1 << 30
	type contrib struct {
		order int
		coeff *big.Rat
	}
	var contribs []contrib
	for i := 0; i <= n; i++ {
		if !o.c[i].IsPolynomial() {
			return ZeroPoly(), false
		}
		p := o.c[i].Num().ScalarMul(ratInv(o.c[i].Den().LeadingCoeff()))
		for k := 0; k <= p.Degree(); k++ {
			ck := p.Coeff(k)
			if ratZero(ck) {
				continue
			}
			exp := k - i
			if exp < minExp {
				minExp = exp
			}
		}
	}
	for i := 0; i <= n; i++ {
		p := o.c[i].Num().ScalarMul(ratInv(o.c[i].Den().LeadingCoeff()))
		for k := 0; k <= p.Degree(); k++ {
			ck := p.Coeff(k)
			if ratZero(ck) {
				continue
			}
			if k-i == minExp {
				contribs = append(contribs, contrib{order: i, coeff: ck})
			}
		}
	}
	// Build sum c * fallingFactorial(s, order) as a polynomial in s.
	acc := ZeroPoly()
	for _, ct := range contribs {
		acc = acc.Add(fallingFactorial(ct.order).ScalarMul(ct.coeff))
	}
	if acc.IsZero() {
		return ZeroPoly(), false
	}
	return acc.Monic(), true
}

// fallingFactorial returns the polynomial s(s-1)...(s-k+1) in the variable s.
func fallingFactorial(k int) Poly {
	acc := OnePoly()
	for j := 0; j < k; j++ {
		acc = acc.Mul(NewPoly(ratInt(int64(-j)), ratInt(1)))
	}
	return acc
}

// Adjoint returns the formal adjoint operator L^* defined by L^*(y) = sum_i
// (-1)^i (a_i y)^(i), obtained by the standard alternating-sign transpose.
func (o Operator) Adjoint() Operator {
	acc := ZeroOperator()
	for i := 0; i < len(o.c); i++ {
		// term: (-1)^i D^i ∘ a_i
		term := ConstOperator(o.c[i])
		for k := 0; k < i; k++ {
			term = dApply(term)
		}
		if i%2 == 1 {
			term = term.Neg()
		}
		acc = acc.Add(term)
	}
	return acc
}

// String renders the operator in descending order of D.
func (o Operator) String() string {
	if o.IsZero() {
		return "0"
	}
	var b strings.Builder
	first := true
	for i := len(o.c) - 1; i >= 0; i-- {
		if o.c[i].IsZero() {
			continue
		}
		if !first {
			b.WriteString(" + ")
		}
		first = false
		b.WriteString("(")
		b.WriteString(o.c[i].String())
		b.WriteString(")")
		if i == 1 {
			b.WriteString("*D")
		} else if i > 1 {
			b.WriteString("*D^")
			b.WriteString(itoa(i))
		}
	}
	return b.String()
}
