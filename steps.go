package algebra

import (
	"math"
	"math/big"
	"strconv"
	"strings"
)

// This file adds a fully worked, step-by-step solution facility to the CAS.
// Because it lives in package algebra it can traverse the real, unexported
// expression tree (*sum, *product, *power, *fn, ...) and reuse the existing
// LaTeX renderer (see [LaTeX]) so that every intermediate expression is
// rendered exactly as it would be anywhere else in the package.
//
// The design has two data types, [Step] and [Solution], and a family of
// generator functions (DifferentiateSteps, IntegrateSteps, SolveQuadraticSteps,
// ...) each of which returns a *[Solution] whose Steps slice records every
// intermediate transformation with a plain-language explanation. No step is
// ever skipped: leaves of the tree emit their own base-case steps and internal
// nodes emit a combining step, so the reader sees all of the work.
//
// Correctness is guaranteed by delegating the actual mathematics to the
// existing, tested engine ([Diff], [Integrate], [Simplify], [Expand], [Solve],
// [Factor], [PartialFractions], [SolveSystem], [Limit], [Series]); the
// generators only *narrate* that work, they never invent results.

// Step is a single, clearly-defined move in a worked solution. It records what
// was done (Explanation, in plain language), the name of the rule or identity
// applied (Rule), and the mathematics immediately before and after the move
// (Before and After, as ordinary [Expr] values). A step whose Before and After
// are structurally equal represents a "state" line — the starting expression or
// a restatement — rather than a transformation.
type Step struct {
	// Explanation is a plain-language description of exactly what was done and
	// why, written to stand on its own to a reader.
	Explanation string
	// Rule is the short name of the rule, law or identity applied, such as
	// "power rule", "sum rule", "quadratic formula" or "l'Hopital's rule".
	Rule string
	// Before is the expression as it stood before this step was applied.
	Before Expr
	// After is the expression that results from applying this step. For a pure
	// "state" line After equals Before.
	After Expr
}

// mkStep is the internal constructor for a [Step].
func mkStep(rule, explanation string, before, after Expr) Step {
	return Step{Rule: rule, Explanation: explanation, Before: before, After: after}
}

// String renders the step as a single human-readable line of the form
// "[rule] explanation: before -> after" (the arrow and after are omitted when
// the step does not change the expression).
func (s Step) String() string {
	var b strings.Builder
	if s.Rule != "" {
		b.WriteString("[" + s.Rule + "] ")
	}
	b.WriteString(s.Explanation)
	b.WriteString("\n    ")
	b.WriteString(s.Before.String())
	if s.After != nil && !s.Before.Equal(s.After) {
		b.WriteString("  ->  " + s.After.String())
	}
	return b.String()
}

// Solution is a complete worked solution: a Title, the ordered list of every
// [Step] taken (with no intermediate work skipped), and the final Result.
type Solution struct {
	// Title is a short description of the problem being solved.
	Title string
	// Steps is the ordered list of every step taken, from first to last.
	Steps []Step
	// Result is the final answer of the solution.
	Result Expr
}

// String renders the whole solution as numbered, human-readable text: the
// title, then each step in order with its explanation and its before/after
// mathematics, then the final result. Every step is shown.
func (sol *Solution) String() string {
	var b strings.Builder
	if sol.Title != "" {
		b.WriteString(sol.Title + "\n")
		b.WriteString(strings.Repeat("=", len([]rune(sol.Title))) + "\n\n")
	}
	for i, s := range sol.Steps {
		b.WriteString(strconv.Itoa(i+1) + ". ")
		if s.Rule != "" {
			b.WriteString("[" + s.Rule + "] ")
		}
		b.WriteString(s.Explanation + "\n")
		b.WriteString("    " + s.Before.String())
		if s.After != nil && !s.Before.Equal(s.After) {
			b.WriteString("  ->  " + s.After.String())
		}
		b.WriteString("\n\n")
	}
	if sol.Result != nil {
		b.WriteString("Result: " + sol.Result.String() + "\n")
	}
	return b.String()
}

// --- small local helpers ---------------------------------------------------

// itoa is a short alias for strconv.Itoa used throughout the generators.
func itoa(n int) string { return strconv.Itoa(n) }

// ordinal returns an English ordinal such as "0th", "1st", "2nd", "3rd",
// "4th" for the small non-negative integers used as derivative orders.
func ordinal(n int) string {
	switch n {
	case 0:
		return "0th"
	case 1:
		return "1st"
	case 2:
		return "2nd"
	case 3:
		return "3rd"
	default:
		return itoa(n) + "th"
	}
}

// isSymbolNamed reports whether e is exactly the symbol named name.
func isSymbolNamed(e Expr, name string) bool {
	s, ok := e.(*Symbol)
	return ok && s.Name == name
}

// errSolution builds a one-step Solution reporting that a problem could not be
// handled, so a generator never returns nil.
func errSolution(title, rule, msg string, e Expr) *Solution {
	return &Solution{
		Title:  title,
		Steps:  []Step{mkStep(rule, msg, e, e)},
		Result: Simplify(e),
	}
}

// ratPolyExpr rebuilds an [Expr] from rational coefficients (index i is the
// coefficient of v^i).
func ratPolyExpr(c []*big.Rat, v Expr) Expr {
	terms := make([]Expr, 0, len(c))
	for i, r := range c {
		if r.Sign() == 0 {
			continue
		}
		terms = append(terms, Mul(newRational(new(big.Rat).Set(r)), Pow(v, Int(int64(i)))))
	}
	return Add(terms...)
}

// polyDegreeIn returns the coefficients of e as a polynomial in name and its
// true degree (highest index with a non-zero coefficient), reporting ok only
// when e is polynomial in name.
func polyDegreeIn(e Expr, name string) (cs []Expr, deg int, ok bool) {
	cs, ok = polyCoeffs(e, name)
	if !ok {
		return nil, 0, false
	}
	deg = len(cs) - 1
	for deg > 0 && isZero(cs[deg]) {
		deg--
	}
	return cs, deg, true
}

// =========================================================================
// Differentiation
// =========================================================================

// DifferentiateSteps returns a [Solution] that differentiates e with respect to
// the symbol v, showing every rule applied at every node of the expression
// tree: the constant, variable, sum, product, constant-multiple, power,
// exponential and chain rules, and the derivative of each elementary function.
// Leaf nodes emit their base-case step and each internal node emits a combining
// step, so no intermediate work is skipped. v must be a [Symbol].
func DifferentiateSteps(e, v Expr) *Solution {
	s, ok := v.(*Symbol)
	if !ok {
		return errSolution("Differentiation", "error",
			"The variable of differentiation must be a symbol.", e)
	}
	res, steps := diffSteps(e, s.Name)
	return &Solution{
		Title:  "Differentiate " + e.String() + " with respect to " + s.Name,
		Steps:  steps,
		Result: res,
	}
}

// diffSteps differentiates e with respect to name, returning the (simplified)
// derivative together with the ordered steps that produced it.
func diffSteps(e Expr, name string) (Expr, []Step) {
	res := Simplify(diff(e, name))
	switch x := e.(type) {
	case *Symbol:
		if x.Name == name {
			return res, []Step{mkStep("variable rule",
				"The derivative of "+name+" with respect to "+name+" is 1.", e, res)}
		}
		return res, []Step{mkStep("constant rule",
			"The symbol "+x.Name+" does not depend on "+name+", so it is a constant and its derivative is 0.", e, res)}
	case *Integer, *Rational, *Float, *Constant:
		return res, []Step{mkStep("constant rule",
			"The derivative of the constant "+e.String()+" is 0.", e, res)}
	case *sum:
		var steps []Step
		for _, a := range x.args {
			_, ss := diffSteps(a, name)
			steps = append(steps, ss...)
		}
		steps = append(steps, mkStep("sum rule",
			"Sum rule: differentiate each term of the sum separately and add the derivatives.", e, res))
		return res, steps
	case *product:
		return res, diffProductSteps(x, name, res)
	case *power:
		return res, diffPowSteps(x, name, res)
	case *fn:
		return res, diffFnSteps(x, name, res)
	case *fn2:
		return res, []Step{mkStep("derivative",
			"Differentiate the two-argument function "+x.name+".", e, res)}
	}
	return res, []Step{mkStep("derivative", "Differentiate "+e.String()+".", e, res)}
}

// diffProductSteps differentiates a product node, using the constant-multiple
// rule when only one factor depends on name and the full product rule
// otherwise.
func diffProductSteps(p *product, name string, res Expr) []Step {
	var consts, rest []Expr
	for _, f := range p.factors {
		if containsSym(f, name) {
			rest = append(rest, f)
		} else {
			consts = append(consts, f)
		}
	}
	var steps []Step
	if len(rest) == 1 && len(consts) >= 1 {
		_, ss := diffSteps(rest[0], name)
		steps = append(steps, ss...)
		steps = append(steps, mkStep("constant multiple rule",
			"Constant multiple rule: the constant factor "+Mul(consts...).String()+
				" is pulled out and multiplies the derivative of the remaining factor.", p, res))
		return steps
	}
	for _, f := range p.factors {
		if containsSym(f, name) {
			_, ss := diffSteps(f, name)
			steps = append(steps, ss...)
		}
	}
	steps = append(steps, mkStep("product rule",
		"Product rule: the derivative of a product is the sum, over each factor, of that factor's derivative times all the other factors.", p, res))
	return steps
}

// diffPowSteps differentiates a power node, distinguishing the power rule, the
// power rule combined with the chain rule, the exponential rule and the general
// (logarithmic-differentiation) rule.
func diffPowSteps(p *power, name string, res Expr) []Step {
	if !containsSym(p, name) {
		return []Step{mkStep("constant rule",
			"The power "+p.String()+" does not depend on "+name+", so its derivative is 0.", p, res)}
	}
	b, x := p.base, p.exp
	switch {
	case !containsSym(x, name):
		var steps []Step
		if !isSymbolNamed(b, name) {
			_, ss := diffSteps(b, name)
			steps = append(steps, ss...)
			steps = append(steps, mkStep("power rule with chain rule",
				"Power rule with the chain rule: d(u^n) = n*u^(n-1)*u' with n = "+x.String()+
					" and u = "+b.String()+".", p, res))
			return steps
		}
		return []Step{mkStep("power rule",
			"Power rule: d("+name+"^n) = n*"+name+"^(n-1) with n = "+x.String()+".", p, res)}
	case !containsSym(b, name):
		var steps []Step
		if !isSymbolNamed(x, name) {
			_, ss := diffSteps(x, name)
			steps = append(steps, ss...)
		}
		steps = append(steps, mkStep("exponential rule",
			"Exponential rule: d(a^u) = a^u*ln(a)*u' with base a = "+b.String()+".", p, res))
		return steps
	default:
		return []Step{mkStep("general power rule",
			"Logarithmic differentiation of a variable base and exponent: d(u^v) = u^v*(v'*ln(u) + v*u'/u).", p, res)}
	}
}

// fnDerivDesc maps a function name to a short textual statement of its
// derivative, used in the explanation of a differentiation step.
var fnDerivDesc = map[string]string{
	"sin":  "d/dx sin(x) = cos(x)",
	"cos":  "d/dx cos(x) = -sin(x)",
	"tan":  "d/dx tan(x) = 1 + tan^2(x)",
	"exp":  "d/dx exp(x) = exp(x)",
	"log":  "d/dx log(x) = 1/x",
	"sqrt": "d/dx sqrt(x) = 1/(2*sqrt(x))",
	"sinh": "d/dx sinh(x) = cosh(x)",
	"cosh": "d/dx cosh(x) = sinh(x)",
	"asin": "d/dx asin(x) = 1/sqrt(1 - x^2)",
	"atan": "d/dx atan(x) = 1/(1 + x^2)",
}

// diffFnSteps differentiates a single-argument function application, inserting a
// chain-rule step when the argument is not simply the variable.
func diffFnSteps(f *fn, name string, res Expr) []Step {
	if !containsSym(f, name) {
		return []Step{mkStep("constant rule",
			"The expression "+f.String()+" does not depend on "+name+", so its derivative is 0.", f, res)}
	}
	desc := fnDerivDesc[f.name]
	if desc == "" {
		desc = "the derivative of " + f.name
	}
	if isSymbolNamed(f.arg, name) {
		return []Step{mkStep(f.name+" rule",
			"Apply the derivative of "+f.name+": "+desc+".", f, res)}
	}
	var steps []Step
	_, ss := diffSteps(f.arg, name)
	steps = append(steps, ss...)
	steps = append(steps, mkStep("chain rule",
		"Chain rule: differentiate the outer function ("+desc+") and multiply by the derivative of the inner expression "+f.arg.String()+".", f, res))
	return steps
}

// =========================================================================
// Integration
// =========================================================================

// IntegrateSteps returns a [Solution] that integrates e with respect to the
// symbol v, naming the technique used at each stage: the sum rule, the
// constant-multiple rule, the power and logarithm rules, the standard
// antiderivatives, integration by parts, partial fractions and substitution.
// The final antiderivative is produced by [Integrate], so it is always correct;
// the steps narrate how it is reached. v must be a [Symbol].
func IntegrateSteps(e, v Expr) *Solution {
	s, ok := v.(*Symbol)
	if !ok {
		return errSolution("Integration", "error",
			"The variable of integration must be a symbol.", e)
	}
	res, steps := integrateSteps(e, s.Name, v)
	return &Solution{
		Title:  "Integrate " + e.String() + " with respect to " + s.Name,
		Steps:  steps,
		Result: res,
	}
}

// integrateSteps integrates e with respect to name, delegating the actual
// antiderivative to [Integrate] and emitting a step per structural level and
// per leaf technique.
func integrateSteps(e Expr, name string, v Expr) (Expr, []Step) {
	res := Integrate(e, v)
	if !containsSym(e, name) {
		return res, []Step{mkStep("constant rule",
			"The integrand does not depend on "+name+", so it integrates to the constant times "+name+".", e, res)}
	}
	switch x := e.(type) {
	case *sum:
		var steps []Step
		for _, a := range x.args {
			_, ss := integrateSteps(a, name, v)
			steps = append(steps, ss...)
		}
		steps = append(steps, mkStep("sum rule",
			"Sum rule for integrals: integrate each term separately and add the antiderivatives.", e, res))
		return res, steps
	case *product:
		var consts, rest []Expr
		for _, f := range x.factors {
			if containsSym(f, name) {
				rest = append(rest, f)
			} else {
				consts = append(consts, f)
			}
		}
		if len(consts) > 0 && len(rest) > 0 {
			var steps []Step
			_, ss := integrateSteps(Mul(rest...), name, v)
			steps = append(steps, ss...)
			steps = append(steps, mkStep("constant multiple rule",
				"Pull the constant factor "+Mul(consts...).String()+" outside the integral.", e, res))
			return res, steps
		}
	case *Symbol:
		return res, []Step{mkStep("power rule",
			"Power rule for integrals: the integral of "+name+" is "+name+"^2/2.", e, res)}
	case *power:
		if bs, ok := x.base.(*Symbol); ok && bs.Name == name {
			if n, ok := x.exp.(*Integer); ok {
				if n.Val.Cmp(big.NewInt(-1)) == 0 {
					return res, []Step{mkStep("logarithm rule",
						"Because the exponent is -1, the integral of "+name+"^(-1) is log("+name+").", e, res)}
				}
				return res, []Step{mkStep("power rule",
					"Power rule for integrals: the integral of "+name+"^n is "+name+"^(n+1)/(n+1) with n = "+x.exp.String()+".", e, res)}
			}
		}
	case *fn:
		if _, isInt := res.(*integral); isInt {
			return res, []Step{mkStep("unresolved",
				"No elementary antiderivative was found; the integral is left unevaluated.", e, res)}
		}
		return res, []Step{mkStep("standard form",
			"Apply the standard antiderivative of "+x.name+" (with a linear-substitution adjustment when the argument is a*"+name+"+b).", e, res)}
	}
	if _, isInt := res.(*integral); isInt {
		return res, []Step{mkStep("unresolved",
			"No elementary antiderivative was found; the integral is left unevaluated.", e, res)}
	}
	rule, expl := integrationTechnique(e, name)
	return res, []Step{mkStep(rule, expl, e, res)}
}

// integrationTechnique classifies the leaf integrand e so its resolving step can
// be labelled with the technique used.
func integrationTechnique(e Expr, name string) (rule, expl string) {
	_, den := numDenom(e)
	if containsSym(den, name) {
		return "partial fractions",
			"The integrand is a rational function; decompose it into partial fractions and integrate each term."
	}
	for _, f := range factorsOf(e) {
		if fc, ok := f.(*fn); ok {
			switch fc.name {
			case "exp", "sin", "cos":
				return "integration by parts",
					"Integrate by parts, differentiating the polynomial factor and integrating the transcendental factor."
			}
		}
	}
	return "substitution",
		"Use a substitution to reduce the integrand to a standard antiderivative."
}

// =========================================================================
// Simplification and expansion
// =========================================================================

// SimplifySteps returns a [Solution] that simplifies e, showing each place in
// the expression tree where an identity is applied (combining like terms,
// folding numeric constants, combining repeated powers, evaluating a function
// at a known value) as its own step. The final result is produced by
// [Simplify].
func SimplifySteps(e Expr) *Solution {
	res, steps := simplifySteps(e)
	if len(steps) == 0 {
		steps = append(steps, mkStep("already simplified",
			"The expression is already in simplest canonical form.", e, res))
	}
	return &Solution{
		Title:  "Simplify " + e.String(),
		Steps:  steps,
		Result: res,
	}
}

// simplifySteps simplifies e bottom-up, emitting a step at every node whose
// canonicalization changes the expression. Raw (non-canonicalizing) node
// constructors are used to capture the "before" form so the change is visible.
func simplifySteps(e Expr) (Expr, []Step) {
	switch x := e.(type) {
	case *sum:
		var steps []Step
		simp := make([]Expr, len(x.args))
		for i, a := range x.args {
			r, ss := simplifySteps(a)
			simp[i] = r
			steps = append(steps, ss...)
		}
		before := newSum(simp)
		after := Simplify(before)
		if before.String() != after.String() {
			steps = append(steps, mkStep("combine like terms",
				"Combine like terms and fold numeric constants in the sum.", before, after))
		}
		return after, steps
	case *product:
		var steps []Step
		simp := make([]Expr, len(x.factors))
		for i, f := range x.factors {
			r, ss := simplifySteps(f)
			simp[i] = r
			steps = append(steps, ss...)
		}
		before := newProduct(simp)
		after := Simplify(before)
		if before.String() != after.String() {
			steps = append(steps, mkStep("combine factors",
				"Multiply numeric factors and combine repeated bases into a single power.", before, after))
		}
		return after, steps
	case *power:
		bR, bs := simplifySteps(x.base)
		eR, es := simplifySteps(x.exp)
		steps := append(append([]Step{}, bs...), es...)
		before := newPower(bR, eR)
		after := Simplify(before)
		if before.String() != after.String() {
			steps = append(steps, mkStep("simplify power",
				"Apply the power identities to "+before.String()+".", before, after))
		}
		return after, steps
	case *fn:
		aR, steps := simplifySteps(x.arg)
		before := newFn(x.name, aR)
		after := Simplify(before)
		if before.String() != after.String() {
			steps = append(steps, mkStep("evaluate function",
				"Evaluate "+x.name+" using its known identities.", before, after))
		}
		return after, steps
	}
	return Simplify(e), nil
}

// ExpandSteps returns a [Solution] that expands e, showing each distribution of
// a product over a sum and each expansion of an integer power of a sum as its
// own step. The final result is produced by [Expand].
func ExpandSteps(e Expr) *Solution {
	res, steps := expandSteps(e)
	if len(steps) == 0 {
		steps = append(steps, mkStep("already expanded",
			"There is nothing to distribute; the expression is already expanded.", e, res))
	}
	return &Solution{
		Title:  "Expand " + e.String(),
		Steps:  steps,
		Result: res,
	}
}

// expandSteps expands e bottom-up, emitting a step wherever a product or power
// distributes.
func expandSteps(e Expr) (Expr, []Step) {
	switch x := e.(type) {
	case *sum:
		var steps []Step
		simp := make([]Expr, len(x.args))
		for i, a := range x.args {
			r, ss := expandSteps(a)
			simp[i] = r
			steps = append(steps, ss...)
		}
		return Add(simp...), steps
	case *product:
		var steps []Step
		simp := make([]Expr, len(x.factors))
		for i, f := range x.factors {
			r, ss := expandSteps(f)
			simp[i] = r
			steps = append(steps, ss...)
		}
		before := newProduct(simp)
		after := Expand(before)
		if before.String() != after.String() {
			steps = append(steps, mkStep("distributive law",
				"Distribute the product "+before.String()+" over the sums it multiplies.", before, after))
		}
		return after, steps
	case *power:
		bR, steps := expandSteps(x.base)
		before := newPower(bR, x.exp)
		after := Expand(before)
		if before.String() != after.String() {
			steps = append(steps, mkStep("binomial expansion",
				"Expand the power "+before.String()+" by repeated multiplication.", before, after))
		}
		return after, steps
	}
	return Expand(e), nil
}

// =========================================================================
// Solving polynomial equations
// =========================================================================

// SolveLinearSteps returns a [Solution] that solves the linear equation e == 0
// for the symbol v, showing the isolate-the-variable-term and
// divide-by-the-coefficient steps explicitly. v must be a [Symbol]; if e is not
// linear in v a single explanatory step is returned instead.
func SolveLinearSteps(e, v Expr) *Solution {
	s, ok := v.(*Symbol)
	title := "Solve the linear equation " + e.String() + " = 0"
	if !ok {
		return errSolution(title, "error", "The variable must be a symbol.", e)
	}
	name := s.Name
	cs, deg, ok := polyDegreeIn(e, name)
	if !ok || deg != 1 {
		return &Solution{Title: title,
			Steps:  []Step{mkStep("not applicable", "The equation is not linear in "+name+".", e, Simplify(e))},
			Result: Simplify(e)}
	}
	a, b := cs[1], cs[0]
	root := Simplify(Mul(neg(b), Pow(a, Int(-1))))
	steps := []Step{
		mkStep("given", "Set the expression equal to zero and solve for "+name+".", e, Int(0)),
		mkStep("collect", "Write in the form a*"+name+" + b = 0 with a = "+a.String()+" and b = "+b.String()+".", e, Collect(e, v)),
		mkStep("isolate variable term", "Subtract the constant b = "+b.String()+" from both sides so that a*"+name+" = "+neg(b).String()+".", Mul(a, v), neg(b)),
		mkStep("divide by coefficient", "Divide both sides by a = "+a.String()+" to solve for "+name+".", v, root),
	}
	return &Solution{Title: title, Steps: steps, Result: root}
}

// SolveQuadraticSteps returns a [Solution] that solves the quadratic equation
// e == 0 for the symbol v: it identifies the coefficients, states the quadratic
// formula, computes the discriminant, and reports both roots (a repeated root
// once, complex conjugate roots when the discriminant is negative). v must be a
// [Symbol]; if e is not quadratic in v a single explanatory step is returned.
func SolveQuadraticSteps(e, v Expr) *Solution {
	s, ok := v.(*Symbol)
	title := "Solve the quadratic equation " + e.String() + " = 0"
	if !ok {
		return errSolution(title, "error", "The variable must be a symbol.", e)
	}
	name := s.Name
	cs, deg, ok := polyDegreeIn(e, name)
	if !ok || deg != 2 {
		return &Solution{Title: title,
			Steps:  []Step{mkStep("not applicable", "The equation is not quadratic in "+name+".", e, Simplify(e))},
			Result: Simplify(e)}
	}
	a, b, c := cs[2], cs[1], cs[0]
	steps := []Step{
		mkStep("given", "Set the expression equal to zero and solve for "+name+".", e, Int(0)),
		mkStep("identify coefficients", "Identify a = "+a.String()+", b = "+b.String()+", c = "+c.String()+".", e, Collect(e, v)),
	}
	// General symbolic quadratic formula, so the rendered derivation shows the
	// familiar (-b +/- sqrt(b^2 - 4ac)) / (2a) with a genuine sqrt and fraction.
	symA, symB, symC := Sym("a"), Sym("b"), Sym("c")
	disc := Add(Pow(symB, Int(2)), Mul(Int(-4), symA, symC))
	denom := Mul(Int(2), symA)
	rootP := Mul(Add(neg(symB), Sqrt(disc)), Pow(denom, Int(-1)))
	rootM := Mul(Add(neg(symB), neg(Sqrt(disc))), Pow(denom, Int(-1)))
	steps = append(steps, mkStep("quadratic formula",
		"Apply the quadratic formula x = (-b +/- sqrt(b^2 - 4ac)) / (2a).", rootP, rootM))
	discNum := Simplify(Add(Pow(b, Int(2)), Mul(Int(-4), a, c)))
	steps = append(steps, mkStep("discriminant",
		"Compute the discriminant D = b^2 - 4ac = "+discNum.String()+".",
		Add(Pow(b, Int(2)), Mul(Int(-4), a, c)), discNum))
	if isNum(discNum) {
		if numSign(discNum) < 0 {
			steps = append(steps, mkStep("negative discriminant",
				"The discriminant is negative, so the two roots are a complex-conjugate pair.", discNum, discNum))
		} else if isZero(discNum) {
			steps = append(steps, mkStep("zero discriminant",
				"The discriminant is zero, so there is a single repeated root.", discNum, discNum))
		}
	}
	roots := solveQuadratic(a, b, c, v)
	for i, r := range roots {
		steps = append(steps, mkStep("solution",
			"Solution "+itoa(i+1)+": substitute the coefficients into the quadratic formula.", v, r))
	}
	return &Solution{Title: title, Steps: steps, Result: roots[0]}
}

// SolveCubicSteps returns a [Solution] that solves the cubic equation e == 0 for
// the symbol v: it searches for a rational root (rational-root theorem),
// divides it out by synthetic division to leave a quadratic factor, and solves
// that quadratic with the quadratic formula. When no rational root exists the
// roots are reported from the numeric solver. v must be a [Symbol].
func SolveCubicSteps(e, v Expr) *Solution {
	s, ok := v.(*Symbol)
	title := "Solve the cubic equation " + e.String() + " = 0"
	if !ok {
		return errSolution(title, "error", "The variable must be a symbol.", e)
	}
	name := s.Name
	cs, deg, ok := polyDegreeIn(e, name)
	if !ok || deg != 3 {
		return &Solution{Title: title,
			Steps:  []Step{mkStep("not applicable", "The equation is not cubic in "+name+".", e, Simplify(e))},
			Result: Simplify(e)}
	}
	steps := []Step{mkStep("given",
		"Set the cubic equal to zero and look for a first root.", e, Int(0))}
	roots, solveErr := Solve(e, v)
	if rc, ok := ratCoeffs(cs); ok && ratDegree(rc) == 3 {
		rc = trimRat(rc)
		rr := rationalRoots(rc)
		if len(rr) > 0 {
			r := rr[0]
			rExpr := newRational(new(big.Rat).Set(r))
			steps = append(steps, mkStep("rational root theorem",
				"Testing the divisors of the constant term over the divisors of the leading coefficient, "+name+" = "+rExpr.String()+" is a root (substituting gives 0).",
				Subs(e, v, rExpr), Int(0)))
			quo := deflate(rc, r)
			quoExpr := ratPolyExpr(quo, v)
			steps = append(steps, mkStep("synthetic division",
				"Divide out the factor ("+name+" - ("+rExpr.String()+")); the quotient is the quadratic "+quoExpr.String()+".",
				e, Mul(Add(v, neg(rExpr)), quoExpr)))
			if ratDegree(quo) == 2 {
				qa := newRational(new(big.Rat).Set(quo[2]))
				qb := newRational(new(big.Rat).Set(quo[1]))
				qc := newRational(new(big.Rat).Set(quo[0]))
				steps = append(steps, mkStep("quadratic formula",
					"Solve the remaining quadratic factor with the quadratic formula.", quoExpr, Int(0)))
				for i, qr := range solveQuadratic(qa, qb, qc, v) {
					steps = append(steps, mkStep("solution",
						"Root from the quadratic factor ("+ordinal(i+1)+").", v, qr))
				}
			}
		} else {
			steps = append(steps, mkStep("no rational root",
				"There is no rational root, so the roots are found numerically (Durand-Kerner iteration).", e, e))
		}
	}
	result := Int(0)
	if solveErr == nil && len(roots) > 0 {
		result = roots[0]
		var parts []string
		for _, r := range roots {
			parts = append(parts, r.String())
		}
		steps = append(steps, mkStep("all solutions",
			"The complete solution set is { "+strings.Join(parts, ", ")+" }.", e, Int(0)))
	}
	return &Solution{Title: title, Steps: steps, Result: result}
}

// =========================================================================
// Factoring, completing the square, partial fractions
// =========================================================================

// FactorSteps returns a [Solution] that factors the univariate quadratic e in
// the symbol v: it computes the discriminant, finds the two roots and writes e
// as a*(v - r1)*(v - r2). For non-quadratic input it reports the factored form
// from [Factor] in a single step. v must be a [Symbol].
func FactorSteps(e, v Expr) *Solution {
	s, ok := v.(*Symbol)
	title := "Factor " + e.String()
	if !ok {
		return errSolution(title, "error", "The variable must be a symbol.", e)
	}
	name := s.Name
	factored := Factor(e, v)
	steps := []Step{mkStep("given", "Factor the polynomial "+e.String()+".", e, e)}
	cs, deg, ok := polyDegreeIn(e, name)
	if ok && deg == 2 {
		a, b, c := cs[2], cs[1], cs[0]
		disc := Simplify(Add(Pow(b, Int(2)), Mul(Int(-4), a, c)))
		steps = append(steps, mkStep("discriminant",
			"Compute the discriminant D = b^2 - 4ac with a = "+a.String()+", b = "+b.String()+", c = "+c.String()+".",
			Add(Pow(b, Int(2)), Mul(Int(-4), a, c)), disc))
		for i, r := range solveQuadratic(a, b, c, v) {
			steps = append(steps, mkStep("root",
				"Root "+itoa(i+1)+" of the quadratic, giving the factor ("+name+" - ("+r.String()+")).", v, r))
		}
		steps = append(steps, mkStep("factored form",
			"Write the quadratic as a*("+name+" - r1)*("+name+" - r2).", e, factored))
	} else {
		steps = append(steps, mkStep("factored form", "The factored form.", e, factored))
	}
	return &Solution{Title: title, Steps: steps, Result: factored}
}

// CompleteSquareSteps returns a [Solution] that rewrites the quadratic e in the
// symbol v in completed-square form a*(v + h)^2 + k, showing the factor-out,
// half-the-linear-coefficient and add-and-subtract steps. v must be a [Symbol].
func CompleteSquareSteps(e, v Expr) *Solution {
	s, ok := v.(*Symbol)
	title := "Complete the square for " + e.String()
	if !ok {
		return errSolution(title, "error", "The variable must be a symbol.", e)
	}
	name := s.Name
	cs, deg, ok := polyDegreeIn(e, name)
	if !ok || deg != 2 {
		return &Solution{Title: title,
			Steps:  []Step{mkStep("not applicable", "Completing the square requires a quadratic in "+name+".", e, Simplify(e))},
			Result: Simplify(e)}
	}
	a, b, c := cs[2], cs[1], cs[0]
	h := Simplify(Mul(b, Pow(Mul(Int(2), a), Int(-1)))) // b/(2a)
	k := Simplify(Add(c, neg(Mul(a, Pow(h, Int(2))))))  // c - a*h^2
	result := Simplify(Add(Mul(a, Pow(Add(v, h), Int(2))), k))
	factored := Add(Mul(a, Add(Pow(v, Int(2)), Mul(Simplify(Mul(b, Pow(a, Int(-1)))), v))), c)
	steps := []Step{
		mkStep("given", "Start with the quadratic a*"+name+"^2 + b*"+name+" + c, where a = "+a.String()+", b = "+b.String()+", c = "+c.String()+".", e, e),
		mkStep("factor leading coefficient", "Factor a = "+a.String()+" out of the "+name+"^2 and "+name+" terms.", e, factored),
		mkStep("half the linear coefficient", "Take half of b/a to get h = b/(2a) = "+h.String()+"; its square is "+Simplify(Pow(h, Int(2))).String()+".", v, h),
		mkStep("complete the square", "Add and subtract a*h^2 to form the perfect square a*("+name+" + h)^2 plus the constant k = c - a*h^2 = "+k.String()+".", e, result),
	}
	return &Solution{Title: title, Steps: steps, Result: result}
}

// PartialFractionSteps returns a [Solution] that decomposes the rational
// expression e in the symbol v into partial fractions: it factors the
// denominator, extracts any polynomial part, and reports each residue term.
// The final decomposition is produced by [ApartExpr]. v must be a [Symbol].
func PartialFractionSteps(e, v Expr) *Solution {
	_, ok := v.(*Symbol)
	title := "Partial fraction decomposition of " + e.String()
	if !ok {
		return errSolution(title, "error", "The variable must be a symbol.", e)
	}
	res, err := ApartExpr(e, v)
	if err != nil {
		return &Solution{Title: title,
			Steps:  []Step{mkStep("not applicable", "The decomposition could not be carried out: "+err.Error()+".", e, Simplify(e))},
			Result: Simplify(e)}
	}
	numE, denE := numDenom(Simplify(e))
	steps := []Step{mkStep("given", "Decompose the rational function "+e.String()+".", e, e)}
	num, e2 := PolyFrom(numE, v)
	den, e3 := PolyFrom(denE, v)
	if e2 == nil && e3 == nil {
		steps = append(steps, mkStep("factor denominator",
			"Factor the denominator "+denE.String()+" into irreducible factors.", denE, Factor(denE, v)))
		poly, terms, perr := PartialFractions(num, den)
		if perr == nil {
			if !poly.IsZero() {
				steps = append(steps, mkStep("polynomial part",
					"The fraction is improper; polynomial division gives the polynomial part.", numE, poly.Expr()))
			}
			for _, t := range terms {
				termExpr := Mul(t.Coeff, Pow(t.Denom.Expr(), Int(int64(-t.Power))))
				denText := t.Denom.String()
				if t.Power > 1 {
					denText = "(" + denText + ")^" + itoa(t.Power)
				}
				steps = append(steps, mkStep("residue",
					"Solve the linear system for the numerator of the term over "+denText+".",
					Pow(t.Denom.Expr(), Int(int64(-t.Power))), termExpr))
			}
		}
	}
	steps = append(steps, mkStep("combine", "Collect the terms into the full partial-fraction decomposition.", e, res))
	return &Solution{Title: title, Steps: steps, Result: res}
}

// =========================================================================
// Linear systems (Gaussian elimination)
// =========================================================================

// rowResidual renders one augmented-matrix row as the residual expression
// (sum of coefficient*symbol) - rhs, i.e. the "= 0" form of that equation.
func rowResidual(row []*big.Rat, syms []Expr) Expr {
	n := len(syms)
	parts := make([]Expr, 0, n+1)
	for j := 0; j < n; j++ {
		parts = append(parts, Mul(newRational(new(big.Rat).Set(row[j])), syms[j]))
	}
	parts = append(parts, neg(newRational(new(big.Rat).Set(row[n]))))
	return Add(parts...)
}

// SolveSystemSteps returns a [Solution] that solves the square linear system
// eqs (each read as eq == 0) for the unknowns syms by Gaussian elimination,
// recording every row operation — pivot selection, row swaps, pivot
// normalization and elimination — as its own step, then reading off each
// unknown. Each entry of syms must be a [Symbol].
func SolveSystemSteps(eqs, syms []Expr) *Solution {
	n := len(syms)
	title := "Solve the linear system by Gaussian elimination"
	fail := func(msg string) *Solution {
		return &Solution{Title: title,
			Steps:  []Step{mkStep("error", msg, Int(0), Int(0))},
			Result: Int(0)}
	}
	if n == 0 || len(eqs) != n {
		return fail("The number of equations must equal the number of unknowns.")
	}
	names := make([]string, n)
	for j, sy := range syms {
		sm, ok := sy.(*Symbol)
		if !ok {
			return fail("Each unknown must be a symbol.")
		}
		names[j] = sm.Name
	}
	// Build the augmented matrix A|b for A x = b.
	m := make([][]*big.Rat, n)
	for i, eq := range eqs {
		m[i] = make([]*big.Rat, n+1)
		for j := range syms {
			cij := Simplify(Diff(eq, syms[j]))
			for _, nm := range names {
				if containsSym(cij, nm) {
					return fail("The system is not linear.")
				}
			}
			r, ok := toRat(cij)
			if !ok {
				return fail("The system must have rational coefficients.")
			}
			m[i][j] = new(big.Rat).Set(r)
		}
		constExpr := eq
		for _, sy := range syms {
			constExpr = Subs(constExpr, sy, Int(0))
		}
		r, ok := toRat(Simplify(constExpr))
		if !ok {
			return fail("The system must have rational constant terms.")
		}
		m[i][n] = new(big.Rat).Neg(r)
	}
	steps := []Step{mkStep("setup",
		"Write each equation in the residual form (linear combination) = 0 and form the augmented matrix.",
		rowResidual(m[0], syms), Int(0))}
	// Gauss-Jordan elimination to reduced row echelon form.
	for col := 0; col < n; col++ {
		piv := -1
		for r := col; r < n; r++ {
			if m[r][col].Sign() != 0 {
				piv = r
				break
			}
		}
		if piv < 0 {
			return fail("The system is singular or underdetermined.")
		}
		if piv != col {
			m[col], m[piv] = m[piv], m[col]
			steps = append(steps, mkStep("row swap",
				"Swap row "+itoa(col+1)+" with row "+itoa(piv+1)+" to bring a non-zero pivot for "+names[col]+" into position.",
				rowResidual(m[col], syms), rowResidual(m[col], syms)))
		}
		if m[col][col].Cmp(big.NewRat(1, 1)) != 0 {
			inv := new(big.Rat).Inv(m[col][col])
			old := rowResidual(m[col], syms)
			for j := col; j <= n; j++ {
				m[col][j].Mul(m[col][j], inv)
			}
			steps = append(steps, mkStep("normalize pivot",
				"Divide row "+itoa(col+1)+" by its pivot so the coefficient of "+names[col]+" becomes 1.",
				old, rowResidual(m[col], syms)))
		}
		for r := 0; r < n; r++ {
			if r == col || m[r][col].Sign() == 0 {
				continue
			}
			factor := new(big.Rat).Set(m[r][col])
			old := rowResidual(m[r], syms)
			for j := col; j <= n; j++ {
				m[r][j].Sub(m[r][j], new(big.Rat).Mul(factor, m[col][j]))
			}
			steps = append(steps, mkStep("eliminate",
				"Eliminate "+names[col]+" from row "+itoa(r+1)+" by subtracting "+factor.RatString()+" times row "+itoa(col+1)+".",
				old, rowResidual(m[r], syms)))
		}
	}
	sol := make([]Expr, n)
	for i := 0; i < n; i++ {
		sol[i] = newRational(new(big.Rat).Set(m[i][n]))
	}
	for j := 0; j < n; j++ {
		steps = append(steps, mkStep("solution",
			"Read off the value of "+names[j]+" from the reduced system.", syms[j], sol[j]))
	}
	return &Solution{Title: title, Steps: steps, Result: sol[0]}
}

// =========================================================================
// Limits and series
// =========================================================================

// LimitSteps returns a [Solution] that evaluates the limit of e as the symbol v
// approaches to. It shows the direct-substitution attempt, detects the 0/0
// indeterminate form and applies l'Hopital's rule, and handles limits at
// infinity by degree comparison. The final value is produced by [Limit]. v must
// be a [Symbol].
func LimitSteps(e, v, to Expr) *Solution {
	s, ok := v.(*Symbol)
	title := "Compute the limit of " + e.String() + " as " + v.String() + " -> " + to.String()
	if !ok {
		return errSolution(title, "error", "The variable must be a symbol.", e)
	}
	name := s.Name
	res := Limit(e, v, to)
	limExpr := newLimit(e, name, to)
	steps := []Step{mkStep("given",
		"Evaluate the limit as "+name+" approaches "+to.String()+".", limExpr, limExpr)}
	if isInfinite(to) {
		steps = append(steps, mkStep("limit at infinity",
			"For a limit at infinity, compare the degrees of the numerator and denominator (the ratio of leading coefficients when the degrees are equal).",
			e, res))
		return &Solution{Title: title, Steps: steps, Result: res}
	}
	num, den := numDenom(e)
	if toVal, err := Evalf(to); err == nil {
		nv, okn := evalAtSafe(num, name, toVal)
		dv, okd := evalAtSafe(den, name, toVal)
		if okd && math.Abs(dv) > 1e-12 {
			sub := Simplify(Subs(e, v, to))
			steps = append(steps, mkStep("direct substitution",
				"The denominator is non-zero at "+name+" = "+to.String()+", so substitute directly.", e, sub))
			return &Solution{Title: title, Steps: steps, Result: res}
		}
		if okd && okn && math.Abs(dv) < 1e-9 && math.Abs(nv) < 1e-9 {
			steps = append(steps, mkStep("indeterminate form",
				"Direct substitution gives 0/0, an indeterminate form, so l'Hopital's rule applies.", e, e))
			np := Simplify(diff(num, name))
			dp := Simplify(diff(den, name))
			ratio := Mul(np, Pow(dp, Int(-1)))
			steps = append(steps, mkStep("l'Hopital's rule",
				"Differentiate the numerator and the denominator separately and take the limit of the ratio.",
				Mul(num, Pow(den, Int(-1))), ratio))
			steps = append(steps, mkStep("evaluate",
				"Evaluate the resulting limit.", ratio, res))
			return &Solution{Title: title, Steps: steps, Result: res}
		}
	}
	steps = append(steps, mkStep("evaluate", "Evaluate the limit.", e, res))
	return &Solution{Title: title, Steps: steps, Result: res}
}

// SeriesSteps returns a [Solution] that builds the Taylor expansion of e about
// the symbol v = at to n terms (degrees 0 through n-1), showing, for each k,
// the k-th derivative, its value at the expansion point, the division by k! and
// the resulting term. The final polynomial is produced by [Series]. v must be a
// [Symbol] and n must be positive.
func SeriesSteps(e, v, at Expr, n int) *Solution {
	s, ok := v.(*Symbol)
	title := "Taylor series of " + e.String() + " about " + v.String() + " = " + at.String() + " to " + itoa(n) + " terms"
	if !ok || n <= 0 {
		return errSolution(title, "error", "The variable must be a symbol and n must be positive.", e)
	}
	name := s.Name
	var steps []Step
	deriv := e
	fact := big.NewInt(1)
	for k := 0; k < n; k++ {
		if k > 0 {
			fact.Mul(fact, big.NewInt(int64(k)))
		}
		coef := Simplify(Subs(deriv, v, at))
		c := Simplify(Mul(coef, Pow(newInteger(new(big.Int).Set(fact)), Int(-1))))
		var basis Expr
		if k == 0 {
			basis = Int(1)
		} else {
			basis = Pow(Add(v, neg(at)), Int(int64(k)))
		}
		term := Simplify(Mul(c, basis))
		steps = append(steps, mkStep("taylor term",
			"Term "+itoa(k)+": the "+ordinal(k)+" derivative is "+deriv.String()+
				"; evaluated at "+name+" = "+at.String()+" it is "+coef.String()+
				", divided by "+itoa(k)+"! = "+fact.String()+".",
			deriv, term))
		deriv = Simplify(diff(deriv, name))
	}
	res := Series(e, v, at, n)
	steps = append(steps, mkStep("sum terms",
		"Add all of the terms to form the truncated Taylor polynomial.", e, res))
	return &Solution{Title: title, Steps: steps, Result: res}
}
