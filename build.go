package algebra

import "math/big"

// Add returns the canonical sum of its arguments. It flattens nested sums,
// drops zero terms, folds numeric constants, combines like terms (so that
// x + x becomes 2*x) and sorts the result. Add() with no arguments is 0.
func Add(args ...Expr) Expr {
	flat := flattenSum(args)
	var order []string
	rests := map[string]Expr{}
	coeffs := map[string]Expr{}
	for _, t := range flat {
		if isZero(t) {
			continue
		}
		c, r := splitCoeff(t)
		key := r.String()
		if _, seen := coeffs[key]; seen {
			coeffs[key] = numAdd(coeffs[key], c)
		} else {
			order = append(order, key)
			rests[key] = r
			coeffs[key] = c
		}
	}
	var out []Expr
	for _, key := range order {
		c, r := coeffs[key], rests[key]
		if isZero(c) {
			continue
		}
		switch {
		case isOne(r):
			out = append(out, c)
		case isOne(c):
			out = append(out, r)
		default:
			out = append(out, mulNumTerm(c, r))
		}
	}
	if len(out) == 0 {
		return Int(0)
	}
	if len(out) == 1 {
		return out[0]
	}
	sortExprs(out)
	return newSum(out)
}

// Mul returns the canonical product of its arguments. It flattens nested
// products, folds numeric constants, applies x*0 and x*1, combines repeated
// bases (so that x*x becomes x^2) and sorts the result. Mul() with no
// arguments is 1.
func Mul(args ...Expr) Expr {
	flat := flattenProduct(args)
	coeff := Expr(Int(1))
	var order []string
	bases := map[string]Expr{}
	exps := map[string]Expr{}
	for _, f := range flat {
		if isZero(f) {
			return Int(0)
		}
		if isNum(f) {
			coeff = numMul(coeff, f)
			continue
		}
		b, e := splitPow(f)
		key := b.String()
		if _, seen := exps[key]; seen {
			exps[key] = Add(exps[key], e)
		} else {
			order = append(order, key)
			bases[key] = b
			exps[key] = e
		}
	}
	var out []Expr
	for _, key := range order {
		p := Pow(bases[key], exps[key])
		if isOne(p) {
			continue
		}
		if isNum(p) {
			coeff = numMul(coeff, p)
			continue
		}
		out = append(out, p)
	}
	if isZero(coeff) {
		return Int(0)
	}
	if len(out) == 0 {
		return coeff
	}
	sortExprs(out)
	if !isOne(coeff) {
		out = append([]Expr{coeff}, out...)
	}
	if len(out) == 1 {
		return out[0]
	}
	return newProduct(out)
}

// Pow returns base raised to exp in canonical form, applying x^0, x^1, 1^x,
// 0^x, folding numeric powers with integer exponents and collapsing nested
// powers with integer outer exponents.
func Pow(base, exp Expr) Expr {
	if isZero(exp) {
		return Int(1)
	}
	if isOne(exp) {
		return base
	}
	if isOne(base) {
		return Int(1)
	}
	if isZero(base) {
		if isNum(exp) && numSign(exp) > 0 {
			return Int(0)
		}
	}
	if isNum(base) && isInteger(exp) {
		if v := numPowInt(base, exp.(*Integer).Val); v != nil {
			return v
		}
	}
	if p, ok := base.(*power); ok && isInteger(exp) {
		return Pow(p.base, Mul(p.exp, exp))
	}
	return newPower(base, exp)
}

// numPowInt raises a numeric base to an integer power exactly (or via float64
// for Float bases). It returns nil to signal "leave symbolic".
func numPowInt(base Expr, n *big.Int) Expr {
	if f, ok := base.(*Float); ok {
		fn := new(big.Float).SetInt(n)
		e, _ := fn.Float64()
		return newFloat(powFloat(f.Val, e))
	}
	r, ok := toRat(base)
	if !ok {
		return nil
	}
	neg := n.Sign() < 0
	abs := new(big.Int).Abs(n)
	if !abs.IsInt64() {
		return nil // absurdly large exponent; leave symbolic
	}
	num := new(big.Int).Exp(r.Num(), abs, nil)
	den := new(big.Int).Exp(r.Denom(), abs, nil)
	res := new(big.Rat).SetFrac(num, den)
	if neg {
		if res.Sign() == 0 {
			return nil
		}
		res.Inv(res)
	}
	return newRational(res)
}

func powFloat(b, e float64) float64 { return mathPow(b, e) }

// flattenSum recursively expands nested sum nodes into a flat slice.
func flattenSum(args []Expr) []Expr {
	var out []Expr
	for _, a := range args {
		if s, ok := a.(*sum); ok {
			out = append(out, flattenSum(s.args)...)
		} else {
			out = append(out, a)
		}
	}
	return out
}

// flattenProduct recursively expands nested product nodes into a flat slice.
func flattenProduct(args []Expr) []Expr {
	var out []Expr
	for _, a := range args {
		if p, ok := a.(*product); ok {
			out = append(out, flattenProduct(p.factors)...)
		} else {
			out = append(out, a)
		}
	}
	return out
}

// splitCoeff separates a term into its numeric coefficient and the remaining
// non-numeric part. A pure number yields (n, 1); a product with a leading
// numeric factor yields (n, rest); anything else yields (1, term).
func splitCoeff(t Expr) (Expr, Expr) {
	if isNum(t) {
		return t, Int(1)
	}
	if p, ok := t.(*product); ok && isNum(p.factors[0]) {
		return p.factors[0], productFrom(p.factors[1:])
	}
	return Int(1), t
}

// splitPow separates a factor into its base and exponent. A non-power factor
// has exponent 1.
func splitPow(f Expr) (Expr, Expr) {
	if p, ok := f.(*power); ok {
		return p.base, p.exp
	}
	return f, Int(1)
}

// productFrom rebuilds a product from an already-canonical factor slice.
func productFrom(factors []Expr) Expr {
	if len(factors) == 0 {
		return Int(1)
	}
	if len(factors) == 1 {
		return factors[0]
	}
	cp := append([]Expr(nil), factors...)
	return newProduct(cp)
}

// mulNumTerm multiplies a numeric coefficient c with a non-numeric remainder r
// without triggering full re-canonicalization.
func mulNumTerm(c, r Expr) Expr {
	if isOne(c) {
		return r
	}
	if isZero(c) {
		return Int(0)
	}
	factors := []Expr{c}
	if p, ok := r.(*product); ok {
		factors = append(factors, p.factors...)
	} else if !isOne(r) {
		factors = append(factors, r)
	}
	if len(factors) == 1 {
		return factors[0]
	}
	return newProduct(factors)
}

// --- elementary function constructors --------------------------------------

// Sin returns sin(x), folding sin(0) to 0.
func Sin(x Expr) Expr {
	if isZero(x) {
		return Int(0)
	}
	return newFn("sin", x)
}

// Cos returns cos(x), folding cos(0) to 1.
func Cos(x Expr) Expr {
	if isZero(x) {
		return Int(1)
	}
	return newFn("cos", x)
}

// Tan returns tan(x), folding tan(0) to 0.
func Tan(x Expr) Expr {
	if isZero(x) {
		return Int(0)
	}
	return newFn("tan", x)
}

// Exp returns e^x, folding exp(0) to 1.
func Exp(x Expr) Expr {
	if isZero(x) {
		return Int(1)
	}
	return newFn("exp", x)
}

// Log returns the natural logarithm log(x), folding log(1) to 0 and log(E)
// to 1.
func Log(x Expr) Expr {
	if isOne(x) {
		return Int(0)
	}
	if c, ok := x.(*Constant); ok && c.Name == "E" {
		return Int(1)
	}
	return newFn("log", x)
}

// Sqrt returns the square root of x. Non-negative integers have their largest
// square factor pulled out of the radical, so sqrt(9) folds to 3 and sqrt(8)
// to 2*sqrt(2).
func Sqrt(x Expr) Expr {
	if isZero(x) {
		return Int(0)
	}
	if isOne(x) {
		return Int(1)
	}
	if n, ok := x.(*Integer); ok && n.Val.Sign() >= 0 {
		out, rad := extractSquare(n.Val)
		if rad.Cmp(one) == 0 {
			return newInteger(out)
		}
		if out.Cmp(one) == 0 {
			return newFn("sqrt", newInteger(rad))
		}
		return Mul(newInteger(out), newFn("sqrt", newInteger(rad)))
	}
	return newFn("sqrt", x)
}

var one = big.NewInt(1)

// applyFn constructs the elementary function named name applied to arg, going
// through the folding constructors so numeric identities apply.
func applyFn(name string, arg Expr) Expr {
	switch name {
	case "sin":
		return Sin(arg)
	case "cos":
		return Cos(arg)
	case "tan":
		return Tan(arg)
	case "exp":
		return Exp(arg)
	case "log":
		return Log(arg)
	case "sqrt":
		return Sqrt(arg)
	}
	return newFn(name, arg)
}
