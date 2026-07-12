package algebra

import (
	"math"
	"math/big"
)

// This file adds the special functions Abs, Sign, Floor, Ceil, Factorial,
// Gamma, Beta, Erf and Erfc. Each is a symbolic node that folds exact numeric
// inputs and evaluates numerically through [Eval]/[Evalf].

// Abs returns the absolute value |x|. Numeric arguments fold exactly, |I| folds
// to 1, and a numeric coefficient is pulled out as |c|*Abs(rest). For a complex
// number a+b*I the modulus sqrt(a^2+b^2) is returned.
func Abs(x Expr) Expr {
	if isNum(x) {
		return numAbs(x)
	}
	if isImagUnit(x) {
		return Int(1)
	}
	if re, im, ok := splitComplex(x); ok && !isZero(im) {
		return Simplify(Sqrt(Add(Pow(re, Int(2)), Pow(im, Int(2)))))
	}
	if p, ok := x.(*product); ok && isNum(p.factors[0]) {
		c := numAbs(p.factors[0])
		rest := productFrom(p.factors[1:])
		if isOne(c) {
			return Abs(rest)
		}
		return Mul(c, Abs(rest))
	}
	return newFn("abs", x)
}

// Sign returns the sign of x: -1, 0 or +1 for numeric inputs. A leading numeric
// factor contributes its sign, so Sign(-3*x) folds to -Sign(x).
func Sign(x Expr) Expr {
	if isNum(x) {
		return Int(int64(numSign(x)))
	}
	if p, ok := x.(*product); ok && isNum(p.factors[0]) {
		s := numSign(p.factors[0])
		rest := productFrom(p.factors[1:])
		if s < 0 {
			return neg(Sign(rest))
		}
		return Sign(rest)
	}
	return newFn("sign", x)
}

// Floor returns the greatest integer not exceeding x for numeric arguments.
func Floor(x Expr) Expr {
	switch n := x.(type) {
	case *Integer:
		return x
	case *Rational:
		return newInteger(ratFloor(n.Val))
	case *Float:
		return newFloat(math.Floor(n.Val))
	}
	return newFn("floor", x)
}

// Ceil returns the least integer not less than x for numeric arguments.
func Ceil(x Expr) Expr {
	switch n := x.(type) {
	case *Integer:
		return x
	case *Rational:
		return newInteger(new(big.Int).Neg(ratFloor(new(big.Rat).Neg(n.Val))))
	case *Float:
		return newFloat(math.Ceil(n.Val))
	}
	return newFn("ceil", x)
}

// Factorial returns x! for a non-negative integer x (folded exactly). For other
// arguments an unevaluated node is returned; it evaluates numerically as
// Gamma(x+1).
func Factorial(x Expr) Expr {
	if n, ok := x.(*Integer); ok && n.Val.Sign() >= 0 {
		return newInteger(bigFactorial(n.Val))
	}
	return newFn("factorial", x)
}

// Gamma returns the gamma function Γ(x). Positive integers fold to a factorial
// (Γ(n) = (n-1)!) and positive half-integers fold to a rational multiple of
// sqrt(pi).
func Gamma(x Expr) Expr {
	if n, ok := x.(*Integer); ok && n.Val.Sign() > 0 {
		return newInteger(bigFactorial(new(big.Int).Sub(n.Val, big.NewInt(1))))
	}
	if r, ok := x.(*Rational); ok && r.Val.Denom().Cmp(big.NewInt(2)) == 0 && r.Val.Sign() > 0 {
		// Γ(p/2) via the recurrence Γ(t) = (t-1)Γ(t-1) down to Γ(1/2)=sqrt(pi).
		acc := Expr(Sqrt(Pi))
		cur := new(big.Rat).Set(r.Val)
		for cur.Cmp(big.NewRat(1, 2)) > 0 {
			cur.Sub(cur, big.NewRat(1, 1))
			acc = Mul(newRational(new(big.Rat).Set(cur)), acc)
		}
		return Simplify(acc)
	}
	return newFn("gamma", x)
}

// Beta returns the beta function B(a, b) = Γ(a)Γ(b)/Γ(a+b) as a symbolic node
// that folds to a [Float] when both arguments are numeric.
func Beta(a, b Expr) Expr {
	if isNum(a) && isNum(b) {
		if v, ok := evalFn2("beta", toFloat(a), toFloat(b)); ok {
			return Flt(v)
		}
	}
	return newFn2("beta", a, b)
}

// Erf returns the error function erf(x), folding erf(0) to 0 and numeric
// [Float] arguments to their value.
func Erf(x Expr) Expr {
	if isZero(x) {
		return Int(0)
	}
	if f, ok := x.(*Float); ok {
		return newFloat(math.Erf(f.Val))
	}
	return newFn("erf", x)
}

// Erfc returns the complementary error function erfc(x) = 1 - erf(x), folding
// erfc(0) to 1 and numeric [Float] arguments to their value.
func Erfc(x Expr) Expr {
	if isZero(x) {
		return Int(1)
	}
	if f, ok := x.(*Float); ok {
		return newFloat(math.Erfc(f.Val))
	}
	return newFn("erfc", x)
}

// --- numeric helpers -------------------------------------------------------

func numAbs(e Expr) Expr {
	switch n := e.(type) {
	case *Integer:
		return newInteger(new(big.Int).Abs(n.Val))
	case *Rational:
		return newRational(new(big.Rat).Abs(new(big.Rat).Set(n.Val)))
	case *Float:
		return newFloat(math.Abs(n.Val))
	}
	return e
}

// ratFloor returns the floor of r as a big.Int. big.Int.Div is floor division
// because a big.Rat always keeps a positive denominator.
func ratFloor(r *big.Rat) *big.Int {
	return new(big.Int).Div(r.Num(), r.Denom())
}

// bigFactorial returns n! for n >= 0.
func bigFactorial(n *big.Int) *big.Int {
	res := big.NewInt(1)
	i := big.NewInt(2)
	for i.Cmp(n) <= 0 {
		res.Mul(res, i)
		i.Add(i, big.NewInt(1))
	}
	return res
}
