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
	case *fn2:
		return diffFn2(x, name)
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
	case "sec":
		return Mul(Sec(u), Tan(u), du)
	case "csc":
		return Mul(Int(-1), Csc(u), Cot(u), du)
	case "cot":
		return Mul(Int(-1), Pow(Csc(u), Int(2)), du)
	case "asin":
		return Mul(du, Pow(Sqrt(Add(Int(1), neg(Pow(u, Int(2))))), Int(-1)))
	case "acos":
		return Mul(Int(-1), du, Pow(Sqrt(Add(Int(1), neg(Pow(u, Int(2))))), Int(-1)))
	case "atan":
		return Mul(du, Pow(Add(Int(1), Pow(u, Int(2))), Int(-1)))
	case "acot":
		return Mul(Int(-1), du, Pow(Add(Int(1), Pow(u, Int(2))), Int(-1)))
	case "asec":
		return Mul(du, Pow(Mul(Abs(u), Sqrt(Add(Pow(u, Int(2)), Int(-1)))), Int(-1)))
	case "acsc":
		return Mul(Int(-1), du, Pow(Mul(Abs(u), Sqrt(Add(Pow(u, Int(2)), Int(-1)))), Int(-1)))
	case "sinh":
		return Mul(Cosh(u), du)
	case "cosh":
		return Mul(Sinh(u), du)
	case "tanh":
		return Mul(Pow(Sech(u), Int(2)), du)
	case "coth":
		return Mul(Int(-1), Pow(Csch(u), Int(2)), du)
	case "sech":
		return Mul(Int(-1), Sech(u), Tanh(u), du)
	case "csch":
		return Mul(Int(-1), Csch(u), Coth(u), du)
	case "asinh":
		return Mul(du, Pow(Sqrt(Add(Pow(u, Int(2)), Int(1))), Int(-1)))
	case "acosh":
		return Mul(du, Pow(Sqrt(Add(Pow(u, Int(2)), Int(-1))), Int(-1)))
	case "atanh":
		return Mul(du, Pow(Add(Int(1), neg(Pow(u, Int(2)))), Int(-1)))
	case "exp":
		return Mul(Exp(u), du)
	case "log":
		return Mul(Pow(u, Int(-1)), du)
	case "sqrt":
		// d/dx sqrt(u) = u' / (2*sqrt(u)).
		return Mul(Rat(1, 2), Pow(Sqrt(u), Int(-1)), du)
	case "abs":
		// d/dx |u| = sign(u) * u'.
		return Mul(Sign(u), du)
	case "sign", "floor", "ceil":
		// Locally constant almost everywhere.
		return Int(0)
	case "gamma":
		return Mul(Gamma(u), newFn("digamma", u), du)
	case "factorial":
		return Mul(Gamma(Add(u, Int(1))), newFn("digamma", Add(u, Int(1))), du)
	case "erf":
		// d/dx erf(u) = 2/sqrt(pi) * exp(-u^2) * u'.
		return Mul(Int(2), Pow(Sqrt(Pi), Int(-1)), Exp(neg(Pow(u, Int(2)))), du)
	case "erfc":
		return Mul(Int(-2), Pow(Sqrt(Pi), Int(-1)), Exp(neg(Pow(u, Int(2)))), du)
	}
	return Int(0)
}

// diffFn2 differentiates the two-argument functions atan2 and beta.
func diffFn2(f *fn2, name string) Expr {
	switch f.name {
	case "atan2":
		// d/dt atan2(y, x) = (x*y' - y*x') / (x^2 + y^2).
		y, x := f.arg1, f.arg2
		dy, dx := diff(y, name), diff(x, name)
		num := Add(Mul(x, dy), neg(Mul(y, dx)))
		den := Add(Pow(x, Int(2)), Pow(y, Int(2)))
		return Mul(num, Pow(den, Int(-1)))
	case "beta":
		// dB/dt = B(a,b)*((psi(a)-psi(a+b))*a' + (psi(b)-psi(a+b))*b').
		a, b := f.arg1, f.arg2
		da, db := diff(a, name), diff(b, name)
		psiab := newFn("digamma", Add(a, b))
		termA := Mul(Add(newFn("digamma", a), neg(psiab)), da)
		termB := Mul(Add(newFn("digamma", b), neg(psiab)), db)
		return Mul(Beta(a, b), Add(termA, termB))
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
		if x.name == "sqrt" {
			// Rewrite sqrt(u) as u^(1/2) so the power rule integrates it, e.g.
			// ∫sqrt(x) dx = (2/3)*x^(3/2).
			return integrate(Pow(x.arg, Rat(1, 2)), name, v)
		}
		if r := integrateFn(x, name, v); r != nil {
			return r
		}
	}
	// Extended strategies for integrands the structural cases did not resolve.
	if r := integrateInvSqrtQuadratic(e, name, v); r != nil {
		return r
	}
	if r := integrateRational(e, name, v); r != nil {
		return r
	}
	if r := integrateByParts(e, name, v); r != nil {
		return r
	}
	return newIntegral(e, v)
}

func integratePow(p *power, name string, v Expr) Expr {
	// The exponent must be numeric (integer or rational) and free of the
	// integration variable. The power rule then applies for any such exponent,
	// so this handles fractional powers like x^(1/2) as well as integer ones.
	if containsSym(p.exp, name) {
		return nil
	}
	n, ok := toRat(p.exp)
	if !ok {
		return nil
	}
	a, _, ok := linearCoeffs(p.base, name, v)
	if !ok {
		return nil
	}
	if n.Cmp(big.NewRat(-1, 1)) == 0 {
		// integral of (a*v+b)^-1 dv = log(a*v+b)/a.
		return Mul(Log(p.base), Pow(a, Int(-1)))
	}
	np1 := new(big.Rat).Add(n, big.NewRat(1, 1)) // n+1
	// integral of (a*v+b)^n dv = (a*v+b)^(n+1) / (a*(n+1)).
	return Mul(Pow(p.base, newRational(np1)), Pow(Mul(a, newRational(np1)), Int(-1)))
}

func integrateFn(f *fn, name string, v Expr) Expr {
	a, _, ok := linearCoeffs(f.arg, name, v)
	if !ok {
		return nil
	}
	inv := Pow(a, Int(-1))
	u := f.arg
	switch f.name {
	case "exp":
		return Mul(Exp(u), inv)
	case "sin":
		return Mul(Int(-1), Cos(u), inv)
	case "cos":
		return Mul(Sin(u), inv)
	case "tan":
		// ∫tan = -log(cos).
		return Mul(Int(-1), Log(Cos(u)), inv)
	case "cot":
		// ∫cot = log(sin).
		return Mul(Log(Sin(u)), inv)
	case "sec":
		// ∫sec = log(sec + tan).
		return Mul(Log(Add(Sec(u), Tan(u))), inv)
	case "csc":
		// ∫csc = -log(csc + cot).
		return Mul(Int(-1), Log(Add(Csc(u), Cot(u))), inv)
	case "sinh":
		return Mul(Cosh(u), inv)
	case "cosh":
		return Mul(Sinh(u), inv)
	case "tanh":
		// ∫tanh = log(cosh).
		return Mul(Log(Cosh(u)), inv)
	case "coth":
		// ∫coth = log(sinh).
		return Mul(Log(Sinh(u)), inv)
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
