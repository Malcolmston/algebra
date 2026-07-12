package algebra

import (
	"errors"
	"fmt"
	"math"
)

// Simplify returns a canonicalized copy of e. Because the [Add], [Mul] and
// [Pow] constructors already canonicalize, Simplify simply rebuilds the tree
// bottom-up through them: this folds numeric arithmetic, applies the algebraic
// identities, combines like terms in sums, combines repeated powers in
// products and evaluates function identities such as sin(0)=0.
func Simplify(e Expr) Expr {
	switch x := e.(type) {
	case *sum:
		args := make([]Expr, len(x.args))
		for i, a := range x.args {
			args[i] = Simplify(a)
		}
		return Add(args...)
	case *product:
		fs := make([]Expr, len(x.factors))
		for i, f := range x.factors {
			fs[i] = Simplify(f)
		}
		return Mul(fs...)
	case *power:
		return Pow(Simplify(x.base), Simplify(x.exp))
	case *fn:
		return applyFn(x.name, Simplify(x.arg))
	case *integral:
		return newIntegral(Simplify(x.arg), x.v)
	}
	return e
}

// Expand distributes products over sums and expands non-negative integer
// powers of sums (binomial/multinomial expansion), then re-simplifies.
func Expand(e Expr) Expr {
	switch x := e.(type) {
	case *sum:
		args := make([]Expr, len(x.args))
		for i, a := range x.args {
			args[i] = Expand(a)
		}
		return Add(args...)
	case *product:
		acc := Expr(Int(1))
		for _, f := range x.factors {
			acc = distribute(acc, Expand(f))
		}
		return acc
	case *power:
		base := Expand(x.base)
		if n, ok := x.exp.(*Integer); ok && n.Val.Sign() > 0 && n.Val.IsInt64() {
			acc := Expr(Int(1))
			for k := int64(0); k < n.Val.Int64(); k++ {
				acc = distribute(acc, base)
			}
			return acc
		}
		return Pow(base, Expand(x.exp))
	case *fn:
		return applyFn(x.name, Expand(x.arg))
	}
	return e
}

// distribute multiplies two expanded expressions, cross-multiplying their
// terms and summing the results so like terms combine.
func distribute(a, b Expr) Expr {
	ta, tb := termsOf(a), termsOf(b)
	out := make([]Expr, 0, len(ta)*len(tb))
	for _, x := range ta {
		for _, y := range tb {
			out = append(out, Mul(x, y))
		}
	}
	return Add(out...)
}

// termsOf returns the additive terms of e (the arguments of a sum, or e).
func termsOf(e Expr) []Expr {
	if s, ok := e.(*sum); ok {
		return s.args
	}
	return []Expr{e}
}

// factorsOf returns the multiplicative factors of e (the factors of a product,
// or e).
func factorsOf(e Expr) []Expr {
	if p, ok := e.(*product); ok {
		return p.factors
	}
	return []Expr{e}
}

// Subs replaces every occurrence of the symbol sym with val and rebuilds the
// expression through the canonicalizing constructors. sym must be a [Symbol].
func Subs(e, sym, val Expr) Expr {
	s, ok := sym.(*Symbol)
	if !ok {
		return e
	}
	return subst(e, s.Name, val)
}

func subst(e Expr, name string, val Expr) Expr {
	switch x := e.(type) {
	case *Symbol:
		if x.Name == name {
			return val
		}
		return e
	case *sum:
		args := make([]Expr, len(x.args))
		for i, a := range x.args {
			args[i] = subst(a, name, val)
		}
		return Add(args...)
	case *product:
		fs := make([]Expr, len(x.factors))
		for i, f := range x.factors {
			fs[i] = subst(f, name, val)
		}
		return Mul(fs...)
	case *power:
		return Pow(subst(x.base, name, val), subst(x.exp, name, val))
	case *fn:
		return applyFn(x.name, subst(x.arg, name, val))
	case *integral:
		return newIntegral(subst(x.arg, name, val), x.v)
	}
	return e
}

// containsSym reports whether e references the symbol named name.
func containsSym(e Expr, name string) bool {
	switch x := e.(type) {
	case *Symbol:
		return x.Name == name
	case *sum:
		for _, a := range x.args {
			if containsSym(a, name) {
				return true
			}
		}
	case *product:
		for _, f := range x.factors {
			if containsSym(f, name) {
				return true
			}
		}
	case *power:
		return containsSym(x.base, name) || containsSym(x.exp, name)
	case *fn:
		return containsSym(x.arg, name)
	case *integral:
		return containsSym(x.arg, name)
	}
	return false
}

// Eval numerically evaluates e to a float64 using env to look up symbol
// values. It returns an error if e references an unbound symbol or an
// operation it cannot evaluate numerically.
func Eval(e Expr, env map[string]float64) (float64, error) {
	switch x := e.(type) {
	case *Integer, *Rational, *Float, *Constant:
		return toFloat(e), nil
	case *Symbol:
		if env != nil {
			if v, ok := env[x.Name]; ok {
				return v, nil
			}
		}
		return 0, fmt.Errorf("algebra: unbound symbol %q", x.Name)
	case *sum:
		total := 0.0
		for _, a := range x.args {
			v, err := Eval(a, env)
			if err != nil {
				return 0, err
			}
			total += v
		}
		return total, nil
	case *product:
		total := 1.0
		for _, f := range x.factors {
			v, err := Eval(f, env)
			if err != nil {
				return 0, err
			}
			total *= v
		}
		return total, nil
	case *power:
		b, err := Eval(x.base, env)
		if err != nil {
			return 0, err
		}
		p, err := Eval(x.exp, env)
		if err != nil {
			return 0, err
		}
		return math.Pow(b, p), nil
	case *fn:
		v, err := Eval(x.arg, env)
		if err != nil {
			return 0, err
		}
		switch x.name {
		case "sin":
			return math.Sin(v), nil
		case "cos":
			return math.Cos(v), nil
		case "tan":
			return math.Tan(v), nil
		case "exp":
			return math.Exp(v), nil
		case "log":
			return math.Log(v), nil
		case "sqrt":
			return math.Sqrt(v), nil
		}
	}
	return 0, errors.New("algebra: cannot evaluate " + e.String())
}

// Evalf numerically evaluates a fully numeric expression (one with no free
// symbols). It is Eval with an empty environment.
func Evalf(e Expr) (float64, error) { return Eval(e, nil) }

// mathPow is a thin wrapper so build.go can raise floats without importing
// math directly.
func mathPow(b, e float64) float64 { return math.Pow(b, e) }
