package groebner

import (
	"strconv"
	"strings"
)

// Monomial represents a power product x_1^e_1 * ... * x_n^e_n as its vector of
// non-negative integer exponents. The length of the slice is the number of
// variables in the ambient polynomial ring.
type Monomial []int

// NewMonomial returns a copy of the given exponent vector as a Monomial. The
// input is not retained.
func NewMonomial(exps ...int) Monomial {
	m := make(Monomial, len(exps))
	copy(m, exps)
	return m
}

// ZeroMonomial returns the constant monomial 1 in n variables, i.e. the vector
// of n zero exponents.
func ZeroMonomial(n int) Monomial {
	return make(Monomial, n)
}

// VarMonomial returns the monomial x_i (the i-th variable to the first power)
// in a ring with n variables. It panics only if i is out of range.
func VarMonomial(n, i int) Monomial {
	m := make(Monomial, n)
	m[i] = 1
	return m
}

// Nvars returns the number of variables (the length of the exponent vector).
func (m Monomial) Nvars() int { return len(m) }

// Clone returns an independent copy of the monomial.
func (m Monomial) Clone() Monomial {
	c := make(Monomial, len(m))
	copy(c, m)
	return c
}

// Degree returns the total degree of the monomial, the sum of its exponents.
func (m Monomial) Degree() int {
	s := 0
	for _, e := range m {
		s += e
	}
	return s
}

// TotalDegree is an alias for Degree, the sum of all exponents.
func (m Monomial) TotalDegree() int { return m.Degree() }

// Exp returns the exponent of the i-th variable.
func (m Monomial) Exp(i int) int { return m[i] }

// IsConstant reports whether the monomial equals 1, i.e. every exponent is zero.
func (m Monomial) IsConstant() bool {
	for _, e := range m {
		if e != 0 {
			return false
		}
	}
	return true
}

// Equal reports whether two monomials have identical exponent vectors.
func (m Monomial) Equal(o Monomial) bool {
	if len(m) != len(o) {
		return false
	}
	for i := range m {
		if m[i] != o[i] {
			return false
		}
	}
	return true
}

// Mul returns the product of two monomials, obtained by adding exponents
// componentwise.
func (m Monomial) Mul(o Monomial) Monomial {
	n := len(m)
	if len(o) > n {
		n = len(o)
	}
	r := make(Monomial, n)
	for i := 0; i < n; i++ {
		if i < len(m) {
			r[i] += m[i]
		}
		if i < len(o) {
			r[i] += o[i]
		}
	}
	return r
}

// Divides reports whether m divides o, i.e. every exponent of m is at most the
// corresponding exponent of o.
func (m Monomial) Divides(o Monomial) bool {
	if len(m) != len(o) {
		return false
	}
	for i := range m {
		if m[i] > o[i] {
			return false
		}
	}
	return true
}

// Div returns the quotient o/m when m divides o. The boolean result is false
// (and the returned monomial nil) when the division is not exact.
func (m Monomial) Div(o Monomial) (Monomial, bool) {
	if !m.Divides(o) {
		return nil, false
	}
	r := make(Monomial, len(o))
	for i := range o {
		r[i] = o[i] - m[i]
	}
	return r, true
}

// LCM returns the least common multiple of two monomials, taking the
// componentwise maximum of the exponents.
func (m Monomial) LCM(o Monomial) Monomial {
	n := len(m)
	if len(o) > n {
		n = len(o)
	}
	r := make(Monomial, n)
	for i := 0; i < n; i++ {
		a, b := 0, 0
		if i < len(m) {
			a = m[i]
		}
		if i < len(o) {
			b = o[i]
		}
		if a > b {
			r[i] = a
		} else {
			r[i] = b
		}
	}
	return r
}

// GCD returns the greatest common divisor of two monomials, taking the
// componentwise minimum of the exponents.
func (m Monomial) GCD(o Monomial) Monomial {
	n := len(m)
	if len(o) < n {
		n = len(o)
	}
	r := make(Monomial, n)
	for i := 0; i < n; i++ {
		if m[i] < o[i] {
			r[i] = m[i]
		} else {
			r[i] = o[i]
		}
	}
	return r
}

// Coprime reports whether two monomials have disjoint support, i.e. their GCD
// is 1. This is the condition used by Buchberger's first (coprime) criterion.
func (m Monomial) Coprime(o Monomial) bool {
	n := len(m)
	if len(o) < n {
		n = len(o)
	}
	for i := 0; i < n; i++ {
		if m[i] != 0 && o[i] != 0 {
			return false
		}
	}
	return true
}

// Support returns the sorted indices of variables that appear with a positive
// exponent in the monomial.
func (m Monomial) Support() []int {
	var s []int
	for i, e := range m {
		if e > 0 {
			s = append(s, i)
		}
	}
	return s
}

// Max returns the componentwise maximum of two monomials (identical to LCM);
// it is provided as a descriptive name for exponent-vector joins.
func (m Monomial) Max(o Monomial) Monomial { return m.LCM(o) }

// Pow returns the monomial raised to a non-negative integer power, multiplying
// every exponent by k.
func (m Monomial) Pow(k int) Monomial {
	r := make(Monomial, len(m))
	for i := range m {
		r[i] = m[i] * k
	}
	return r
}

// String renders the monomial using the default variable names x1, x2, ....
func (m Monomial) String() string {
	return m.Format(defaultVarNames(len(m)))
}

// Format renders the monomial using the supplied variable names. The constant
// monomial is rendered as "1".
func (m Monomial) Format(vars []string) string {
	if m.IsConstant() {
		return "1"
	}
	var b strings.Builder
	first := true
	for i, e := range m {
		if e == 0 {
			continue
		}
		if !first {
			b.WriteString("*")
		}
		first = false
		name := "x" + strconv.Itoa(i+1)
		if i < len(vars) && vars[i] != "" {
			name = vars[i]
		}
		b.WriteString(name)
		if e != 1 {
			b.WriteString("^")
			b.WriteString(strconv.Itoa(e))
		}
	}
	return b.String()
}

func defaultVarNames(n int) []string {
	v := make([]string, n)
	for i := range v {
		v[i] = "x" + strconv.Itoa(i+1)
	}
	return v
}
