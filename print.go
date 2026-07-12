package algebra

import (
	"sort"
	"strconv"
	"strings"
)

// String renders a symbol as its name.
func (s *Symbol) String() string { return s.Name }

// String renders an integer in base 10.
func (i *Integer) String() string { return i.Val.String() }

// String renders a rational as numerator/denominator.
func (r *Rational) String() string {
	return r.Val.Num().String() + "/" + r.Val.Denom().String()
}

// String renders a float using the shortest round-tripping representation.
func (f *Float) String() string {
	return strconv.FormatFloat(f.Val, 'g', -1, 64)
}

// String renders a named constant.
func (c *Constant) String() string { return c.Name }

// String renders a function application, e.g. sin(x).
func (f *fn) String() string { return f.name + "(" + f.arg.String() + ")" }

// String renders an unevaluated integral, e.g. Integral(f(x), x).
func (n *integral) String() string {
	return "Integral(" + n.arg.String() + ", " + n.v.String() + ")"
}

// String renders a power with parenthesization of compound bases and
// exponents.
func (p *power) String() string {
	return wrap(p.base) + "^" + wrap(p.exp)
}

// String renders a product using * as the separator, factoring out a leading
// -1 as a unary minus and parenthesizing embedded sums.
func (p *product) String() string {
	fs := p.factors
	prefix := ""
	if len(fs) > 1 && isMinusOne(fs[0]) {
		prefix = "-"
		fs = fs[1:]
	}
	parts := make([]string, 0, len(fs))
	for _, f := range fs {
		s := f.String()
		if _, ok := f.(*sum); ok {
			s = "(" + s + ")"
		}
		parts = append(parts, s)
	}
	return prefix + strings.Join(parts, "*")
}

// String renders a sum in descending polynomial-degree order using + and -
// so that, e.g., x^2 - 2*x + 1 prints naturally.
func (s *sum) String() string {
	terms := append([]Expr(nil), s.args...)
	// Collect all variables so equal-degree terms can be ordered by a stable
	// descending exponent vector (giving textbook order like a^2 + 2*a*b + b^2).
	varSet := map[string]bool{}
	exps := make([]map[string]int, len(terms))
	for i, t := range terms {
		exps[i] = termExps(t)
		for v := range exps[i] {
			varSet[v] = true
		}
	}
	vars := make([]string, 0, len(varSet))
	for v := range varSet {
		vars = append(vars, v)
	}
	sort.Strings(vars)
	idx := map[Expr]int{}
	for i, t := range terms {
		idx[t] = i
	}
	sort.SliceStable(terms, func(i, j int) bool {
		di, dj := polyDegree(terms[i]), polyDegree(terms[j])
		if di != dj {
			return di > dj
		}
		ei, ej := exps[idx[terms[i]]], exps[idx[terms[j]]]
		for _, v := range vars {
			if ei[v] != ej[v] {
				return ei[v] > ej[v]
			}
		}
		return compareExpr(terms[i], terms[j]) < 0
	})
	var b strings.Builder
	for i, t := range terms {
		neg, mag := splitSign(t)
		ms := mag.String()
		if _, ok := mag.(*sum); ok {
			ms = "(" + ms + ")"
		}
		switch {
		case i == 0 && neg:
			b.WriteString("-" + ms)
		case i == 0:
			b.WriteString(ms)
		case neg:
			b.WriteString(" - " + ms)
		default:
			b.WriteString(" + " + ms)
		}
	}
	return b.String()
}

// wrap parenthesizes e when it appears as the base or exponent of a power.
func wrap(e Expr) string {
	if needParen(e) {
		return "(" + e.String() + ")"
	}
	return e.String()
}

func needParen(e Expr) bool {
	switch e.(type) {
	case *sum, *product, *power, *Rational:
		return true
	}
	if isNum(e) && numSign(e) < 0 {
		return true
	}
	return false
}

// splitSign returns whether the term carries a negative leading coefficient
// together with its positive-magnitude form, used by sum printing.
func splitSign(t Expr) (bool, Expr) {
	if isNum(t) {
		if numSign(t) < 0 {
			return true, numNeg(t)
		}
		return false, t
	}
	if p, ok := t.(*product); ok && isNum(p.factors[0]) && numSign(p.factors[0]) < 0 {
		return true, mulNumTerm(numNeg(p.factors[0]), productFrom(p.factors[1:]))
	}
	return false, t
}

// termExps returns the per-variable exponent map of a monomial term, used to
// order same-degree terms for printing.
func termExps(e Expr) map[string]int {
	m := map[string]int{}
	for _, f := range factorsOf(e) {
		switch x := f.(type) {
		case *Symbol:
			m[x.Name]++
		case *power:
			if b, ok := x.base.(*Symbol); ok {
				if n, ok := x.exp.(*Integer); ok && n.Val.IsInt64() {
					m[b.Name] += int(n.Val.Int64())
				}
			}
		}
	}
	return m
}

// polyDegree estimates the total polynomial degree of e, used only to order
// terms for printing.
func polyDegree(e Expr) int {
	switch x := e.(type) {
	case *Symbol:
		return 1
	case *power:
		if n, ok := x.exp.(*Integer); ok && n.Val.Sign() > 0 {
			return int(n.Val.Int64()) * polyDegree(x.base)
		}
		return 0
	case *product:
		d := 0
		for _, f := range x.factors {
			d += polyDegree(f)
		}
		return d
	case *sum:
		d := 0
		for _, a := range x.args {
			if pd := polyDegree(a); pd > d {
				d = pd
			}
		}
		return d
	case *fn:
		return 1
	}
	return 0
}
