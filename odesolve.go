package algebra

import (
	"fmt"
	"math/big"
	"sort"
)

// ODESolution is the result of symbolically solving an ordinary differential
// equation. It mirrors the information returned by a computer-algebra dsolve.
type ODESolution struct {
	// General is the general solution. For explicit solutions it is the
	// right-hand side f such that y = f(x) (a formula that does not itself
	// mention the unknown-function symbol y). For implicit solutions it is a
	// relation expression R that is understood to satisfy R == 0, and it may
	// mention both x and y. Which form applies is indicated by Kind and can be
	// tested by whether General references the y symbol.
	General Expr
	// Kind names the solution method that succeeded, for example "separable",
	// "linear", "bernoulli", "exact", "homogeneous" or one of the
	// second-order constant-coefficient cases. Implicit first-order results use
	// the suffix "-implicit".
	Kind string
	// Constants holds the arbitrary-constant symbols introduced in General, in
	// order: [C1] for first-order equations and [C1, C2] for second-order ones.
	Constants []Expr
}

// ODEError reports that no implemented method could solve the differential
// equation, or that its inputs were malformed.
type ODEError struct{ Reason string }

// Error returns the human-readable reason the ODE could not be solved.
func (e *ODEError) Error() string { return "algebra: " + e.Reason }

// odeC1 returns the first arbitrary-constant symbol C1.
func odeC1() Expr { return Sym("C1") }

// odeC2 returns the second arbitrary-constant symbol C2.
func odeC2() Expr { return Sym("C2") }

// odeSymName returns the name of e when it is a symbol.
func odeSymName(e Expr) (string, bool) {
	if s, ok := e.(*Symbol); ok {
		return s.Name, true
	}
	return "", false
}

// --- exponential post-processing -------------------------------------------

// odeMatchLogTerm reports whether the additive term t has the form c*log(u),
// returning u and the coefficient c (the product of the remaining factors).
func odeMatchLogTerm(t Expr) (u, c Expr, ok bool) {
	var logArg Expr
	var coeff []Expr
	n := 0
	for _, f := range factorsOf(t) {
		if fnn, is := f.(*fn); is && fnn.name == "log" {
			logArg = fnn.arg
			n++
			continue
		}
		coeff = append(coeff, f)
	}
	if n != 1 {
		return nil, nil, false
	}
	return logArg, Mul(coeff...), true
}

// odeExp returns exp(e) with logarithmic terms folded into powers, using the
// identity exp(c*log(u)) = u^c. This keeps integrating factors such as
// exp(Integrate(P, x)) in a canonical rational form (for example
// exp(-log(x)) becomes x^(-1)) so that later cancellation succeeds.
func odeExp(e Expr) Expr {
	e = Simplify(e)
	var facs []Expr
	var rest []Expr
	for _, t := range termsOf(e) {
		if u, c, ok := odeMatchLogTerm(t); ok {
			facs = append(facs, Pow(u, c))
		} else {
			rest = append(rest, t)
		}
	}
	if len(rest) > 0 {
		facs = append(facs, Exp(Add(rest...)))
	}
	if len(facs) == 0 {
		return Int(1)
	}
	return Simplify(Mul(facs...))
}

// odeCombineExp merges exponential factors within every product node, folding
// exp(a)*exp(b) into exp(a+b) so that manufactured factors such as
// exp(x)*exp(-x) collapse to 1. It is a best-effort normaliser used when
// checking that a residual vanishes.
func odeCombineExp(e Expr) Expr {
	switch x := e.(type) {
	case *sum:
		args := make([]Expr, len(x.args))
		for i, a := range x.args {
			args[i] = odeCombineExp(a)
		}
		return Add(args...)
	case *product:
		var expArgs []Expr
		var others []Expr
		for _, f := range x.factors {
			cf := odeCombineExp(f)
			if fnn, ok := cf.(*fn); ok && fnn.name == "exp" {
				expArgs = append(expArgs, fnn.arg)
			} else {
				others = append(others, cf)
			}
		}
		if len(expArgs) > 0 {
			others = append(others, Exp(Add(expArgs...)))
		}
		return Mul(others...)
	case *power:
		return Pow(odeCombineExp(x.base), odeCombineExp(x.exp))
	case *fn:
		return applyFn(x.name, odeCombineExp(x.arg))
	case *fn2:
		return applyFn2(x.name, odeCombineExp(x.arg1), odeCombineExp(x.arg2))
	}
	return e
}

// odeSimplify canonicalises e for presentation and for zero-testing, expanding
// products, combining exponentials and re-simplifying to a fixed point.
func odeSimplify(e Expr) Expr {
	prev := e
	for i := 0; i < 4; i++ {
		next := odeCombineExp(Simplify(Expand(prev)))
		next = Simplify(next)
		if next.Equal(prev) {
			return next
		}
		prev = next
	}
	return prev
}

// --- first-order solver ----------------------------------------------------

// SolveODE1 solves the first-order ordinary differential equation
// dy/dx = rhs, where rhs is f(x, y). The symbols x (independent variable) and
// y (unknown function) are supplied as [Symbol] values.
//
// It attempts, in this deterministic order, the classical closed-form methods,
// returning the first that applies:
//
//  1. separable      rhs factors as g(x)*h(y);
//  2. linear         y' + P(x)y = Q(x);
//  3. Bernoulli      y' + P(x)y = Q(x)*y^n with numeric n not 0 or 1;
//  4. exact          M(x,y) + N(x,y)y' = 0 with dM/dy == dN/dx;
//  5. homogeneous    rhs depends only on the ratio y/x.
//
// An arbitrary constant C1 is introduced. Solutions are returned explicitly
// (y = f(x)) when the relation can be inverted and implicitly (R(x,y) == 0)
// otherwise. If no method applies an *[ODEError] is returned.
func SolveODE1(rhs, x, y Expr) (*ODESolution, error) {
	xn, ok := odeSymName(x)
	if !ok {
		return nil, &ODEError{Reason: "x must be a symbol"}
	}
	yn, ok := odeSymName(y)
	if !ok {
		return nil, &ODEError{Reason: "y must be a symbol"}
	}
	rhs = Simplify(rhs)
	if sol, ok := odeTrySeparable(rhs, x, y, xn, yn); ok {
		return sol, nil
	}
	if sol, ok := odeTryLinear(rhs, x, y, xn, yn); ok {
		return sol, nil
	}
	if sol, ok := odeTryBernoulli(rhs, x, y, xn, yn); ok {
		return sol, nil
	}
	if sol, ok := odeTryExact(rhs, x, y, xn, yn); ok {
		return sol, nil
	}
	if sol, ok := odeTryHomogeneous(rhs, x, y, xn, yn); ok {
		return sol, nil
	}
	return nil, &ODEError{Reason: "no first-order method applies to dy/dx = " + rhs.String()}
}

// odeSeparable splits rhs into a product g(x)*h(y). It reports false when any
// single factor depends on both x and y (so that no structural separation
// exists).
func odeSeparable(rhs Expr, xname, yname string) (g, h Expr, ok bool) {
	gg := Expr(Int(1))
	hh := Expr(Int(1))
	for _, f := range factorsOf(Simplify(rhs)) {
		cx := containsSym(f, xname)
		cy := containsSym(f, yname)
		if cx && cy {
			return nil, nil, false
		}
		if cy {
			hh = Mul(hh, f)
		} else {
			gg = Mul(gg, f)
		}
	}
	return Simplify(gg), Simplify(hh), true
}

// odeMatchLogY reports whether e equals c*log(y) for the unknown-function
// symbol named yname, returning the coefficient c.
func odeMatchLogY(e Expr, yname string) (Expr, bool) {
	var coeff []Expr
	logs := 0
	for _, f := range factorsOf(e) {
		if fnn, ok := f.(*fn); ok && fnn.name == "log" {
			if s, ok := fnn.arg.(*Symbol); ok && s.Name == yname {
				logs++
				continue
			}
			return nil, false
		}
		if containsSym(f, yname) {
			return nil, false
		}
		coeff = append(coeff, f)
	}
	if logs != 1 {
		return nil, false
	}
	return Mul(coeff...), true
}

// odeTrySeparable applies the separable-equation method.
func odeTrySeparable(rhs, x, y Expr, xname, yname string) (*ODESolution, bool) {
	g, h, ok := odeSeparable(rhs, xname, yname)
	if !ok {
		return nil, false
	}
	ih := Simplify(Integrate(Pow(h, Int(-1)), y))
	ig := Simplify(Integrate(g, x))
	rhsEq := Add(ig, odeC1())
	// Prefer an explicit solution when the y-side integral inverts.
	if sols, err := Solve(Simplify(Add(ih, neg(rhsEq))), y); err == nil && len(sols) >= 1 {
		return &ODESolution{General: odeSimplify(sols[0]), Kind: "separable", Constants: []Expr{odeC1()}}, true
	}
	if k, ok := odeMatchLogY(ih, yname); ok {
		yexpr := Exp(Mul(rhsEq, Pow(k, Int(-1))))
		return &ODESolution{General: odeSimplify(yexpr), Kind: "separable", Constants: []Expr{odeC1()}}, true
	}
	lhs := Simplify(Add(ih, neg(rhsEq)))
	return &ODESolution{General: lhs, Kind: "separable-implicit", Constants: []Expr{odeC1()}}, true
}

// odeTryLinear applies the first-order linear method y' + P(x)y = Q(x).
func odeTryLinear(rhs, x, y Expr, xname, yname string) (*ODESolution, bool) {
	a := Simplify(Diff(rhs, y))
	if containsSym(a, yname) {
		return nil, false
	}
	b := Simplify(Subs(rhs, y, Int(0)))
	// Confirm rhs == a*y + b so the equation is genuinely linear in y.
	if !Simplify(Expand(Add(rhs, neg(Add(Mul(a, y), b))))).Equal(Int(0)) {
		return nil, false
	}
	p := Simplify(neg(a))
	q := b
	mu := odeExp(Integrate(p, x))
	intPart := Simplify(Integrate(Simplify(Mul(mu, q)), x))
	yexpr := Mul(Add(intPart, odeC1()), Pow(mu, Int(-1)))
	return &ODESolution{General: odeSimplify(yexpr), Kind: "linear", Constants: []Expr{odeC1()}}, true
}

// odeTryBernoulli applies the Bernoulli method y' + P(x)y = Q(x)*y^n.
func odeTryBernoulli(rhs, x, y Expr, xname, yname string) (*ODESolution, bool) {
	coeffByExp := map[string]Expr{}
	ratByKey := map[string]*big.Rat{}
	var order []string
	for _, t := range termsOf(Simplify(rhs)) {
		yexp := big.NewRat(0, 1)
		var coeff []Expr
		pure := true
		for _, f := range factorsOf(t) {
			if s, ok := f.(*Symbol); ok && s.Name == yname {
				yexp.Add(yexp, big.NewRat(1, 1))
				continue
			}
			if p, ok := f.(*power); ok {
				if bs, ok := p.base.(*Symbol); ok && bs.Name == yname {
					if r, ok := toRat(p.exp); ok {
						yexp.Add(yexp, r)
						continue
					}
					pure = false
					break
				}
			}
			if containsSym(f, yname) {
				pure = false
				break
			}
			coeff = append(coeff, f)
		}
		if !pure {
			return nil, false
		}
		key := yexp.RatString()
		if _, ex := ratByKey[key]; !ex {
			ratByKey[key] = new(big.Rat).Set(yexp)
			order = append(order, key)
			coeffByExp[key] = Mul(coeff...)
		} else {
			coeffByExp[key] = Add(coeffByExp[key], Mul(coeff...))
		}
	}
	if len(order) != 2 {
		return nil, false
	}
	one := big.NewRat(1, 1)
	var nRat *big.Rat
	oneKey, nKey := "", ""
	for _, k := range order {
		if ratByKey[k].Cmp(one) == 0 {
			oneKey = k
		} else {
			nKey = k
			nRat = ratByKey[k]
		}
	}
	if oneKey == "" || nKey == "" {
		return nil, false
	}
	if nRat.Sign() == 0 || nRat.Cmp(one) == 0 {
		return nil, false
	}
	// rhs = -P*y + Q*y^n, so the y^1 coefficient is -P.
	p := Simplify(neg(coeffByExp[oneKey]))
	q := coeffByExp[nKey]
	n := newRational(new(big.Rat).Set(nRat))
	oneMinusN := Simplify(Add(Int(1), neg(n)))
	pw := Simplify(Mul(oneMinusN, p))
	qw := Simplify(Mul(oneMinusN, q))
	mu := odeExp(Integrate(pw, x))
	intPart := Simplify(Integrate(Simplify(Mul(mu, qw)), x))
	w := Mul(Add(intPart, odeC1()), Pow(mu, Int(-1)))
	expo := Simplify(Pow(oneMinusN, Int(-1)))
	yexpr := Pow(odeSimplify(w), expo)
	return &ODESolution{General: odeSimplify(yexpr), Kind: "bernoulli", Constants: []Expr{odeC1()}}, true
}

// odeTryExact applies the exact-equation method. Writing dy/dx = rhs as
// M(x,y) + N(x,y)y' = 0 with M = -numerator(rhs) and N = denominator(rhs), it
// checks dM/dy == dN/dx and, when exact, reconstructs the potential F with
// F(x,y) = C1.
func odeTryExact(rhs, x, y Expr, xname, yname string) (*ODESolution, bool) {
	f := Simplify(rhs)
	num, den := numDenom(f)
	m := Simplify(neg(num))
	nn := Simplify(den)
	dMdy := Simplify(Diff(m, y))
	dNdx := Simplify(Diff(nn, x))
	if !Simplify(Add(dMdy, neg(dNdx))).Equal(Int(0)) {
		return nil, false
	}
	// The relation must genuinely couple x and y; otherwise this is an already
	// handled degenerate case.
	if !containsSym(m, yname) && !containsSym(nn, xname) {
		return nil, false
	}
	fx := Integrate(m, x)
	if _, unresolved := fx.(*integral); unresolved {
		return nil, false
	}
	gy := Simplify(Add(nn, neg(Diff(fx, y))))
	if containsSym(gy, xname) {
		return nil, false
	}
	gInt := Integrate(gy, y)
	if _, unresolved := gInt.(*integral); unresolved {
		return nil, false
	}
	potential := Simplify(Add(fx, gInt))
	general := Simplify(Add(potential, neg(odeC1())))
	return &ODESolution{General: general, Kind: "exact", Constants: []Expr{odeC1()}}, true
}

// odeTryHomogeneous applies the homogeneous method: when rhs depends only on
// v = y/x, the substitution y = v*x reduces the equation to a separable one in
// v and x, x*dv/dx = F(v) - v, whose solution is left implicit in y/x.
//
// Degree-zero homogeneity is detected exactly with Euler's identity
// x*f_x + y*f_y == 0, and the reduced right-hand side F(v) is recovered as
// f(1, v). This avoids relying on the simplifier to cancel the common factor
// of x that a direct substitution y = v*x would leave behind.
func odeTryHomogeneous(rhs, x, y Expr, xname, yname string) (*ODESolution, bool) {
	if !containsSym(rhs, xname) || !containsSym(rhs, yname) {
		return nil, false
	}
	euler := Simplify(Add(Mul(x, Diff(rhs, x)), Mul(y, Diff(rhs, y))))
	if !euler.Equal(Int(0)) {
		return nil, false
	}
	v := Sym("ode_v")
	fv := Simplify(Subs(Subs(rhs, x, Int(1)), y, v)) // F(v) = f(1, v)
	if containsSym(fv, xname) {
		return nil, false
	}
	denom := Simplify(Add(fv, neg(v)))
	yOverX := Mul(y, Pow(x, Int(-1)))
	if denom.Equal(Int(0)) {
		// F(v) == v gives v' == 0, so y/x is constant.
		lhs := Simplify(Add(yOverX, neg(odeC1())))
		return &ODESolution{General: lhs, Kind: "homogeneous", Constants: []Expr{odeC1()}}, true
	}
	iv := Simplify(Integrate(Pow(denom, Int(-1)), v))
	lhs := Simplify(Add(iv, neg(Log(x)), neg(odeC1())))
	lhs = Simplify(Subs(lhs, v, yOverX))
	return &ODESolution{General: lhs, Kind: "homogeneous", Constants: []Expr{odeC1()}}, true
}

// VerifyODE1 checks a first-order solution by substitution: for an explicit
// solution y = sol.General it forms the residual y' - rhs (with y replaced by
// the solution) and reports whether it simplifies exactly to zero. It is a
// best-effort exact check and returns false for implicit solutions (whose
// General still mentions y) or when the residual cannot be reduced to zero.
func VerifyODE1(sol *ODESolution, rhs, x, y Expr) bool {
	if sol == nil {
		return false
	}
	yn, ok := odeSymName(y)
	if !ok {
		return false
	}
	if _, ok := odeSymName(x); !ok {
		return false
	}
	if containsSym(sol.General, yn) {
		return false
	}
	f := sol.General
	dfdx := Diff(f, x)
	sub := Subs(rhs, y, f)
	res := odeSimplify(Add(dfdx, neg(sub)))
	return res.Equal(Int(0))
}

// --- second-order constant-coefficient solver ------------------------------

// SolveODE2Const solves the second-order linear constant-coefficient equation
//
//	a*y'' + b*y' + c*y = g(x)
//
// where a, b and c are numeric [Expr] and x is the independent-variable
// [Symbol]. The homogeneous part is built from the characteristic roots of
// a*r^2 + b*r + c, covering the distinct-real, repeated-real and
// complex-conjugate cases with arbitrary constants C1 and C2. A particular
// solution is found by undetermined coefficients when g is a sum of terms of
// the form (polynomial)*exp(alpha*x)*{1, cos(omega*x), sin(omega*x)} with
// rational alpha and omega; resonance with the homogeneous solution is handled
// by the usual multiplication by a power of x. If g falls outside that family
// an *[ODEError] is returned.
func SolveODE2Const(a, b, c, g, x, y Expr) (*ODESolution, error) {
	xn, ok := odeSymName(x)
	if !ok {
		return nil, &ODEError{Reason: "x must be a symbol"}
	}
	if !isNum(a) || !isNum(b) || !isNum(c) {
		return nil, &ODEError{Reason: "coefficients a, b, c must be numeric"}
	}
	if numSign(a) == 0 {
		return nil, &ODEError{Reason: "leading coefficient a must be nonzero"}
	}
	yh, kind := odeHomo2(a, b, c, x)
	yp := Expr(Int(0))
	gS := Simplify(g)
	if !isZero(gS) {
		part, err := odeParticularAll(a, b, c, gS, x, xn)
		if err != nil {
			return nil, err
		}
		yp = part
	}
	general := odeSimplify(Add(yh, yp))
	return &ODESolution{General: general, Kind: kind, Constants: []Expr{odeC1(), odeC2()}}, nil
}

// odeHomo2 builds the homogeneous solution of a*y” + b*y' + c*y = 0 from the
// discriminant, returning the solution expression and a Kind naming the case.
func odeHomo2(a, b, c, x Expr) (Expr, string) {
	d := Simplify(Add(Pow(b, Int(2)), Mul(Int(-4), a, c)))
	twoa := Mul(Int(2), a)
	switch numSign(d) {
	case 1:
		s := sqrtExpr(d)
		r1 := Simplify(Mul(Add(neg(b), s), Pow(twoa, Int(-1))))
		r2 := Simplify(Mul(Add(neg(b), neg(s)), Pow(twoa, Int(-1))))
		yh := Add(Mul(odeC1(), Exp(Mul(r1, x))), Mul(odeC2(), Exp(Mul(r2, x))))
		return yh, "const-coeff-distinct-real"
	case 0:
		r := Simplify(Mul(neg(b), Pow(twoa, Int(-1))))
		yh := Mul(Add(odeC1(), Mul(odeC2(), x)), Exp(Mul(r, x)))
		return yh, "const-coeff-repeated"
	default:
		alpha := Simplify(Mul(neg(b), Pow(twoa, Int(-1))))
		beta := Simplify(Mul(sqrtExpr(neg(d)), Pow(twoa, Int(-1))))
		yh := Mul(Exp(Mul(alpha, x)),
			Add(Mul(odeC1(), Cos(Mul(beta, x))), Mul(odeC2(), Sin(Mul(beta, x)))))
		return yh, "const-coeff-complex"
	}
}

// odeClassifyTerm decomposes a single forcing term into its polynomial degree
// in x, exponential rate alpha, trig frequency omega and trig kind. It reports
// false when the term is outside the undetermined-coefficients family.
func odeClassifyTerm(t, x Expr, xname string) (deg int, alpha, omega Expr, trig string, ok bool) {
	alpha = Int(0)
	omega = Int(0)
	trig = ""
	var expArg, trigArg Expr
	for _, f := range factorsOf(Simplify(t)) {
		if isNum(f) {
			continue
		}
		switch node := f.(type) {
		case *Symbol:
			if node.Name == xname {
				deg++
				continue
			}
			return 0, nil, nil, "", false
		case *power:
			if bs, ok := node.base.(*Symbol); ok && bs.Name == xname {
				if n, ok := node.exp.(*Integer); ok && n.Val.Sign() >= 0 && n.Val.IsInt64() {
					deg += int(n.Val.Int64())
					continue
				}
			}
			return 0, nil, nil, "", false
		case *fn:
			switch node.name {
			case "exp":
				if expArg != nil {
					return 0, nil, nil, "", false
				}
				expArg = node.arg
			case "cos":
				if trig != "" {
					return 0, nil, nil, "", false
				}
				trig = "cos"
				trigArg = node.arg
			case "sin":
				if trig != "" {
					return 0, nil, nil, "", false
				}
				trig = "sin"
				trigArg = node.arg
			default:
				return 0, nil, nil, "", false
			}
		default:
			return 0, nil, nil, "", false
		}
	}
	if expArg != nil {
		al := Simplify(Diff(expArg, x))
		if containsSym(al, xname) {
			return 0, nil, nil, "", false
		}
		if !Simplify(Add(expArg, neg(Mul(al, x)))).Equal(Int(0)) {
			return 0, nil, nil, "", false
		}
		alpha = al
	}
	if trigArg != nil {
		w := Simplify(Diff(trigArg, x))
		if containsSym(w, xname) {
			return 0, nil, nil, "", false
		}
		if !Simplify(Add(trigArg, neg(Mul(w, x)))).Equal(Int(0)) {
			return 0, nil, nil, "", false
		}
		omega = w
	}
	return deg, alpha, omega, trig, true
}

// odeParticularAll finds a particular solution for the whole forcing gS by
// grouping its additive terms by (alpha, omega) family and superposing the
// undetermined-coefficients solution of each group.
func odeParticularAll(a, b, c, gS, x Expr, xname string) (Expr, error) {
	type grp struct {
		alpha, omega Expr
		hasTrig      bool
		deg          int
		terms        []Expr
	}
	var groups []*grp
	idx := map[string]int{}
	for _, t := range termsOf(gS) {
		if isZero(t) {
			continue
		}
		d, al, om, _, ok := odeClassifyTerm(t, x, xname)
		if !ok {
			return nil, &ODEError{Reason: "forcing term outside undetermined-coefficients family: " + t.String()}
		}
		key := al.String() + "||" + om.String()
		gi, exists := idx[key]
		if !exists {
			groups = append(groups, &grp{alpha: al, omega: om, hasTrig: !isZero(om), deg: d})
			idx[key] = len(groups) - 1
			gi = len(groups) - 1
		}
		gr := groups[gi]
		if d > gr.deg {
			gr.deg = d
		}
		gr.terms = append(gr.terms, t)
	}
	yp := Expr(Int(0))
	for _, gr := range groups {
		part, ok := odeUndetermined(a, b, c, Add(gr.terms...), x, xname, gr.alpha, gr.omega, gr.deg, gr.hasTrig)
		if !ok {
			return nil, &ODEError{Reason: "cannot determine coefficients for forcing group"}
		}
		yp = Add(yp, part)
	}
	return Simplify(yp), nil
}

// odeRealRootMult returns the multiplicity (0, 1 or 2) of the real number
// alpha as a root of a*r^2 + b*r + c.
func odeRealRootMult(a, b, c, alpha Expr) int {
	val := Simplify(Add(Mul(a, Pow(alpha, Int(2))), Mul(b, alpha), c))
	if !val.Equal(Int(0)) {
		return 0
	}
	d := Simplify(Add(Pow(b, Int(2)), Mul(Int(-4), a, c)))
	if d.Equal(Int(0)) {
		return 2
	}
	return 1
}

// odeComplexResonance reports whether alpha+omega*i is a root of
// a*r^2 + b*r + c (with omega nonzero), i.e. the forcing resonates with a
// complex-conjugate homogeneous mode.
func odeComplexResonance(a, b, c, alpha, omega Expr) bool {
	p := Simplify(Add(Mul(a, Add(Pow(alpha, Int(2)), neg(Pow(omega, Int(2))))), Mul(b, alpha), c))
	q := Simplify(Mul(omega, Add(Mul(Int(2), a, alpha), b)))
	return p.Equal(Int(0)) && q.Equal(Int(0))
}

// odeUndetermined solves for the particular solution of one forcing group via
// undetermined coefficients, building a template of basis functions (bumped by
// x^s for resonance), applying the operator a*D^2 + b*D + c and solving the
// resulting linear system exactly over the rationals.
func odeUndetermined(a, b, c, gGroup, x Expr, xname string, alpha, omega Expr, deg int, hasTrig bool) (Expr, bool) {
	s := 0
	if !hasTrig {
		s = odeRealRootMult(a, b, c, alpha)
	} else if odeComplexResonance(a, b, c, alpha, omega) {
		s = 1
	}
	var phis []Expr
	for k := 0; k <= deg; k++ {
		base := Mul(Pow(x, Int(int64(s+k))), Exp(Mul(alpha, x)))
		if hasTrig {
			phis = append(phis, Simplify(Mul(base, Cos(Mul(omega, x)))))
			phis = append(phis, Simplify(Mul(base, Sin(Mul(omega, x)))))
		} else {
			phis = append(phis, Simplify(base))
		}
	}
	gc, ok := odeCoords(gGroup, xname)
	if !ok {
		return nil, false
	}
	coeffMaps := make([]map[string]*big.Rat, len(phis))
	seen := map[string]bool{}
	var keys []string
	addKeys := func(mp map[string]*big.Rat) {
		for k := range mp {
			if !seen[k] {
				seen[k] = true
				keys = append(keys, k)
			}
		}
	}
	for i, phi := range phis {
		li := Simplify(Add(Mul(a, Diff(Diff(phi, x), x)), Mul(b, Diff(phi, x)), Mul(c, phi)))
		cm, ok := odeCoords(li, xname)
		if !ok {
			return nil, false
		}
		coeffMaps[i] = cm
		addKeys(cm)
	}
	addKeys(gc)
	sort.Strings(keys)
	mat := make([][]*big.Rat, len(keys))
	rhs := make([]*big.Rat, len(keys))
	for ki, k := range keys {
		row := make([]*big.Rat, len(phis))
		for i := range phis {
			if val, ok := coeffMaps[i][k]; ok {
				row[i] = val
			} else {
				row[i] = big.NewRat(0, 1)
			}
		}
		mat[ki] = row
		if val, ok := gc[k]; ok {
			rhs[ki] = val
		} else {
			rhs[ki] = big.NewRat(0, 1)
		}
	}
	sol, ok := odeSolveRows(mat, rhs, len(phis))
	if !ok {
		return nil, false
	}
	terms := make([]Expr, 0, len(phis))
	for i := range phis {
		terms = append(terms, Mul(newRational(sol[i]), phis[i]))
	}
	return Simplify(Add(terms...)), true
}

// odeCoords expresses e as a linear combination of test functions
// x^k*exp(alpha*x)*{1,cos,sin} and returns the rational coefficient of each,
// keyed by a canonical string. It reports false when a term is outside that
// family.
func odeCoords(e Expr, xname string) (map[string]*big.Rat, bool) {
	e = Simplify(Expand(e))
	out := map[string]*big.Rat{}
	for _, t := range termsOf(e) {
		if isZero(t) {
			continue
		}
		key, coeff, ok := odeTermKey(t, xname)
		if !ok {
			return nil, false
		}
		if prev, exists := out[key]; exists {
			prev.Add(prev, coeff)
		} else {
			out[key] = new(big.Rat).Set(coeff)
		}
	}
	return out, true
}

// odeTermKey returns the canonical test-function key and rational coefficient
// of a single term, or false when the term is not of the recognised shape.
func odeTermKey(t Expr, xname string) (string, *big.Rat, bool) {
	coeff := big.NewRat(1, 1)
	xdeg := 0
	var expArg, trigArg Expr
	trig := ""
	for _, f := range factorsOf(t) {
		if isNum(f) {
			r, ok := toRat(f)
			if !ok {
				return "", nil, false
			}
			coeff.Mul(coeff, r)
			continue
		}
		switch node := f.(type) {
		case *Symbol:
			if node.Name == xname {
				xdeg++
				continue
			}
			return "", nil, false
		case *power:
			if bs, ok := node.base.(*Symbol); ok && bs.Name == xname {
				if n, ok := node.exp.(*Integer); ok && n.Val.Sign() >= 0 && n.Val.IsInt64() {
					xdeg += int(n.Val.Int64())
					continue
				}
			}
			return "", nil, false
		case *fn:
			switch node.name {
			case "exp":
				if expArg != nil {
					return "", nil, false
				}
				expArg = node.arg
			case "cos":
				if trig != "" {
					return "", nil, false
				}
				trig = "cos"
				trigArg = node.arg
			case "sin":
				if trig != "" {
					return "", nil, false
				}
				trig = "sin"
				trigArg = node.arg
			default:
				return "", nil, false
			}
		default:
			return "", nil, false
		}
	}
	eaStr := ""
	if expArg != nil {
		eaStr = expArg.String()
	}
	taStr := ""
	if trigArg != nil {
		taStr = trigArg.String()
	}
	key := fmt.Sprintf("x%d|e%s|%s%s", xdeg, eaStr, trig, taStr)
	return key, coeff, true
}

// odeSolveRows solves the (possibly rectangular but consistent) rational
// linear system mat*A = rhs for the n unknowns A by Gauss-Jordan elimination,
// setting any free variable to zero. It reports false on an inconsistent
// system.
func odeSolveRows(mat [][]*big.Rat, rhs []*big.Rat, n int) ([]*big.Rat, bool) {
	m := len(mat)
	aug := make([][]*big.Rat, m)
	for i := 0; i < m; i++ {
		aug[i] = make([]*big.Rat, n+1)
		for j := 0; j < n; j++ {
			if mat[i][j] != nil {
				aug[i][j] = new(big.Rat).Set(mat[i][j])
			} else {
				aug[i][j] = big.NewRat(0, 1)
			}
		}
		if rhs[i] != nil {
			aug[i][n] = new(big.Rat).Set(rhs[i])
		} else {
			aug[i][n] = big.NewRat(0, 1)
		}
	}
	var pivCols []int
	row := 0
	for col := 0; col < n && row < m; col++ {
		sel := -1
		for r := row; r < m; r++ {
			if aug[r][col].Sign() != 0 {
				sel = r
				break
			}
		}
		if sel < 0 {
			continue
		}
		aug[row], aug[sel] = aug[sel], aug[row]
		inv := new(big.Rat).Inv(aug[row][col])
		for j := col; j <= n; j++ {
			aug[row][j].Mul(aug[row][j], inv)
		}
		for r := 0; r < m; r++ {
			if r == row || aug[r][col].Sign() == 0 {
				continue
			}
			factor := new(big.Rat).Set(aug[r][col])
			for j := col; j <= n; j++ {
				aug[r][j].Sub(aug[r][j], new(big.Rat).Mul(factor, aug[row][j]))
			}
		}
		pivCols = append(pivCols, col)
		row++
	}
	for r := 0; r < m; r++ {
		allZero := true
		for j := 0; j < n; j++ {
			if aug[r][j].Sign() != 0 {
				allZero = false
				break
			}
		}
		if allZero && aug[r][n].Sign() != 0 {
			return nil, false
		}
	}
	sol := make([]*big.Rat, n)
	for j := range sol {
		sol[j] = big.NewRat(0, 1)
	}
	for i, col := range pivCols {
		sol[col] = new(big.Rat).Set(aug[i][n])
	}
	return sol, true
}
