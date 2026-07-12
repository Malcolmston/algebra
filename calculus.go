package algebra

import "math/big"

// Diff returns the symbolic derivative of e with respect to the symbol v.
// It implements the sum, product, power, quotient and chain rules and the
// derivatives of sin, cos, tan, exp, log and sqrt. The result is returned in
// canonical (simplified) form. v must be a [Symbol].
func Diff(e, v Expr) Expr {
	s, ok := v.(*Symbol)
	if !ok {
		return Int(0)
	}
	return diff(e, s.Name)
}

func diff(e Expr, name string) Expr {
	switch x := e.(type) {
	case *Integer, *Rational, *Float, *Constant:
		return Int(0)
	case *Symbol:
		if x.Name == name {
			return Int(1)
		}
		return Int(0)
	case *sum:
		terms := make([]Expr, len(x.args))
		for i, a := range x.args {
			terms[i] = diff(a, name)
		}
		return Add(terms...)
	case *product:
		// Product rule: d(f1*f2*...*fn) = sum_i (df_i * prod_{j!=i} f_j).
		terms := make([]Expr, 0, len(x.factors))
		for i := range x.factors {
			parts := []Expr{diff(x.factors[i], name)}
			for j, f := range x.factors {
				if j != i {
					parts = append(parts, f)
				}
			}
			terms = append(terms, Mul(parts...))
		}
		return Add(terms...)
	case *power:
		return diffPow(x, name)
	case *fn:
		return diffFn(x, name)
	}
	return Int(0)
}

func diffPow(p *power, name string) Expr {
	b, x := p.base, p.exp
	db := diff(b, name)
	// Constant exponent: d(u^n) = n*u^(n-1)*u'.
	if !containsSym(x, name) {
		return Mul(x, Pow(b, Add(x, Int(-1))), db)
	}
	dx := diff(x, name)
	// Constant base: d(a^u) = a^u * ln(a) * u'.
	if !containsSym(b, name) {
		return Mul(Pow(b, x), Log(b), dx)
	}
	// General: d(u^v) = u^v * (v'*ln(u) + v*u'/u).
	return Mul(Pow(b, x), Add(Mul(dx, Log(b)), Mul(x, db, Pow(b, Int(-1)))))
}

func diffFn(f *fn, name string) Expr {
	u := f.arg
	du := diff(u, name)
	switch f.name {
	case "sin":
		return Mul(Cos(u), du)
	case "cos":
		return Mul(Int(-1), Sin(u), du)
	case "tan":
		// d/dx tan(u) = (1 + tan^2(u)) * u'.
		return Mul(Add(Int(1), Pow(Tan(u), Int(2))), du)
	case "exp":
		return Mul(Exp(u), du)
	case "log":
		return Mul(Pow(u, Int(-1)), du)
	case "sqrt":
		// d/dx sqrt(u) = u' / (2*sqrt(u)).
		return Mul(Rat(1, 2), Pow(Sqrt(u), Int(-1)), du)
	}
	return Int(0)
}

// Integrate returns a symbolic antiderivative of e with respect to the symbol
// v (the constant of integration is omitted).
//
// Coverage: constants; powers x^n for any integer n (including 1/x -> log x);
// sums (term by term); constant multiples (constants are pulled out); the
// elementary functions exp, sin and cos; and, via a linear substitution, any
// of the above whose argument or base is linear in v (a*v+b). Integrands
// outside this set are returned as an unevaluated [Integral] node.
func Integrate(e, v Expr) Expr {
	s, ok := v.(*Symbol)
	if !ok {
		return newIntegral(e, v)
	}
	return integrate(e, s.Name, v)
}

func integrate(e Expr, name string, v Expr) Expr {
	// Constant with respect to v.
	if !containsSym(e, name) {
		return Mul(e, v)
	}
	switch x := e.(type) {
	case *Symbol:
		// integral of x dx = x^2/2.
		return Mul(Rat(1, 2), Pow(v, Int(2)))
	case *sum:
		terms := make([]Expr, len(x.args))
		for i, a := range x.args {
			terms[i] = integrate(a, name, v)
		}
		return Add(terms...)
	case *product:
		// Pull out factors that do not depend on v.
		var consts, rest []Expr
		for _, f := range x.factors {
			if containsSym(f, name) {
				rest = append(rest, f)
			} else {
				consts = append(consts, f)
			}
		}
		if len(consts) > 0 && len(rest) > 0 {
			return Mul(Mul(consts...), integrate(Mul(rest...), name, v))
		}
	case *power:
		if r := integratePow(x, name, v); r != nil {
			return r
		}
	case *fn:
		if r := integrateFn(x, name, v); r != nil {
			return r
		}
	}
	return newIntegral(e, v)
}

func integratePow(p *power, name string, v Expr) Expr {
	n, ok := p.exp.(*Integer)
	if !ok || containsSym(p.exp, name) {
		return nil
	}
	a, _, ok := linearCoeffs(p.base, name, v)
	if !ok {
		return nil
	}
	if n.Val.Cmp(big.NewInt(-1)) == 0 {
		// integral of (a*v+b)^-1 dv = log(a*v+b)/a.
		return Mul(Log(p.base), Pow(a, Int(-1)))
	}
	np1 := new(big.Int).Add(n.Val, big.NewInt(1)) // n+1
	// integral of (a*v+b)^n dv = (a*v+b)^(n+1) / (a*(n+1)).
	return Mul(Pow(p.base, newInteger(np1)), Pow(Mul(a, newInteger(np1)), Int(-1)))
}

func integrateFn(f *fn, name string, v Expr) Expr {
	a, _, ok := linearCoeffs(f.arg, name, v)
	if !ok {
		return nil
	}
	inv := Pow(a, Int(-1))
	switch f.name {
	case "exp":
		return Mul(Exp(f.arg), inv)
	case "sin":
		return Mul(Int(-1), Cos(f.arg), inv)
	case "cos":
		return Mul(Sin(f.arg), inv)
	}
	return nil
}

// linearCoeffs reports whether e is a linear expression a*v+b in the symbol
// named name and, if so, returns the constant coefficients a and b. It works
// by differentiating (a must be constant) and evaluating at v=0 (giving b),
// then verifying the reconstruction equals e.
func linearCoeffs(e Expr, name string, v Expr) (a, b Expr, ok bool) {
	da := Simplify(diff(e, name))
	if containsSym(da, name) {
		return nil, nil, false
	}
	b = Simplify(subst(e, name, Int(0)))
	recon := Simplify(Expand(Add(Mul(da, v), b)))
	if !Simplify(Expand(e)).Equal(recon) {
		return nil, nil, false
	}
	return da, b, true
}
