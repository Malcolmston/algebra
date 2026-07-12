package algebra

import (
	"math"
	"math/big"
)

// Expr is the interface implemented by every node in an expression tree.
//
// Expr values are immutable and value-based. In addition to the reflective
// helpers below (String, Equal), Expr carries fluent builder methods so that
// expressions can be assembled left to right, e.g. x.Pow(Int(2)).Add(Int(1)).
type Expr interface {
	// String renders the expression as a readable infix string with correct
	// operator precedence and parenthesization.
	String() string
	// Equal reports whether other is structurally identical to the receiver.
	Equal(other Expr) bool
	// Add returns the canonical sum of the receiver and the arguments.
	Add(others ...Expr) Expr
	// Mul returns the canonical product of the receiver and the arguments.
	Mul(others ...Expr) Expr
	// Pow returns the receiver raised to exp.
	Pow(exp Expr) Expr
	// Diff returns the derivative of the receiver with respect to v.
	Diff(v Expr) Expr
	// Simplify returns a canonicalized, simplified copy of the receiver.
	Simplify() Expr
	// Expand distributes products over sums and expands integer powers.
	Expand() Expr
	// Subs substitutes every occurrence of sym with val.
	Subs(sym, val Expr) Expr

	implExpr()
}

// builders supplies the fluent Expr methods to every node type. The self field
// holds the wrapping Expr value so the promoted methods can forward to the
// package-level functions.
type builders struct{ self Expr }

func (b builders) implExpr()               {}
func (b builders) Add(o ...Expr) Expr      { return Add(append([]Expr{b.self}, o...)...) }
func (b builders) Mul(o ...Expr) Expr      { return Mul(append([]Expr{b.self}, o...)...) }
func (b builders) Pow(exp Expr) Expr       { return Pow(b.self, exp) }
func (b builders) Diff(v Expr) Expr        { return Diff(b.self, v) }
func (b builders) Simplify() Expr          { return Simplify(b.self) }
func (b builders) Expand() Expr            { return Expand(b.self) }
func (b builders) Subs(sym, val Expr) Expr { return Subs(b.self, sym, val) }

// Symbol is a named variable such as x or y.
type Symbol struct {
	builders
	Name string
}

// Integer is an arbitrary-precision integer.
type Integer struct {
	builders
	Val *big.Int
}

// Rational is an exact fraction whose reduced denominator is greater than one.
// Fractions that reduce to whole numbers are represented as [Integer] instead.
type Rational struct {
	builders
	Val *big.Rat
}

// Float is an inexact floating-point literal.
type Float struct {
	builders
	Val float64
}

// Constant is a named mathematical constant such as pi or E.
type Constant struct {
	builders
	Name string
	num  float64
}

// sum is an n-ary addition node. Its arguments are kept in canonical order.
type sum struct {
	builders
	args []Expr
}

// product is an n-ary multiplication node. Its factors are kept in canonical
// order with any numeric coefficient first.
type product struct {
	builders
	factors []Expr
}

// power is a base raised to an exponent.
type power struct {
	builders
	base, exp Expr
}

// fn is an application of a single-argument elementary function.
type fn struct {
	builders
	name string
	arg  Expr
}

// fn2 is an application of a two-argument function such as atan2 or beta.
type fn2 struct {
	builders
	name       string
	arg1, arg2 Expr
}

// integral is an unevaluated integral, returned by [Integrate] for integrands
// it cannot handle in closed form.
type integral struct {
	builders
	arg Expr
	v   Expr
}

// --- constructors that set the self back-reference -------------------------

func newInteger(v *big.Int) Expr { e := &Integer{Val: v}; e.self = e; return e }

func newRational(r *big.Rat) Expr {
	if r.IsInt() {
		return newInteger(new(big.Int).Set(r.Num()))
	}
	e := &Rational{Val: r}
	e.self = e
	return e
}

func newFloat(f float64) Expr { e := &Float{Val: f}; e.self = e; return e }
func newConst(n string, v float64) Expr {
	e := &Constant{Name: n, num: v}
	e.self = e
	return e
}
func newSymbol(n string) Expr       { e := &Symbol{Name: n}; e.self = e; return e }
func newSum(a []Expr) Expr          { e := &sum{args: a}; e.self = e; return e }
func newProduct(a []Expr) Expr      { e := &product{factors: a}; e.self = e; return e }
func newPower(b, x Expr) Expr       { e := &power{base: b, exp: x}; e.self = e; return e }
func newFn(n string, arg Expr) Expr { e := &fn{name: n, arg: arg}; e.self = e; return e }
func newFn2(n string, a1, a2 Expr) Expr {
	e := &fn2{name: n, arg1: a1, arg2: a2}
	e.self = e
	return e
}
func newIntegral(a, v Expr) Expr { e := &integral{arg: a, v: v}; e.self = e; return e }

// Sym returns the symbol (named variable) with the given name.
func Sym(name string) Expr { return newSymbol(name) }

// Int returns the integer expression for v.
func Int(v int64) Expr { return newInteger(big.NewInt(v)) }

// IntBig returns the integer expression for the given big.Int; the value is
// copied so later mutation of v does not affect the expression.
func IntBig(v *big.Int) Expr { return newInteger(new(big.Int).Set(v)) }

// Rat returns the exact rational a/b in lowest terms. It panics if b is zero.
func Rat(a, b int64) Expr {
	if b == 0 {
		panic("algebra: Rat with zero denominator")
	}
	return newRational(big.NewRat(a, b))
}

// Flt returns the inexact float literal f.
func Flt(f float64) Expr { return newFloat(f) }

// Pi is the mathematical constant pi.
var Pi = newConst("pi", math.Pi)

// E is Euler's number, the base of the natural logarithm.
var E = newConst("E", math.E)

// I is the imaginary unit, satisfying I^2 == -1. It participates in ordinary
// [Add], [Mul] and [Pow] arithmetic; integer powers of I fold (I^2 -> -1) and
// clean Euler identities such as exp(I*Pi) -> -1 are recognised by [Simplify].
var I = newConst("I", math.NaN())

// Inf is positive infinity (∞), used as a target for [Limit] and produced by
// evaluations that diverge. NegInf is negative infinity (-∞).
var (
	Inf    = newConst("oo", math.Inf(1))
	NegInf = newConst("-oo", math.Inf(-1))
)

// isImagUnit reports whether e is the imaginary unit I.
func isImagUnit(e Expr) bool {
	c, ok := e.(*Constant)
	return ok && c.Name == "I"
}

// isInfinite reports whether e is one of the infinity constants.
func isInfinite(e Expr) bool {
	c, ok := e.(*Constant)
	return ok && (c.Name == "oo" || c.Name == "-oo")
}

// --- numeric helpers -------------------------------------------------------

func isNum(e Expr) bool {
	switch e.(type) {
	case *Integer, *Rational, *Float:
		return true
	}
	return false
}

func isInteger(e Expr) bool { _, ok := e.(*Integer); return ok }

func toRat(e Expr) (*big.Rat, bool) {
	switch n := e.(type) {
	case *Integer:
		return new(big.Rat).SetInt(n.Val), true
	case *Rational:
		return n.Val, true
	}
	return nil, false
}

func toFloat(e Expr) float64 {
	switch n := e.(type) {
	case *Integer:
		f, _ := new(big.Float).SetInt(n.Val).Float64()
		return f
	case *Rational:
		f, _ := n.Val.Float64()
		return f
	case *Float:
		return n.Val
	case *Constant:
		return n.num
	}
	return math.NaN()
}

// numSign returns -1, 0 or +1 for a numeric expression.
func numSign(e Expr) int {
	switch n := e.(type) {
	case *Integer:
		return n.Val.Sign()
	case *Rational:
		return n.Val.Sign()
	case *Float:
		if n.Val < 0 {
			return -1
		} else if n.Val > 0 {
			return 1
		}
		return 0
	}
	return 0
}

func isZero(e Expr) bool {
	switch n := e.(type) {
	case *Integer:
		return n.Val.Sign() == 0
	case *Float:
		return n.Val == 0
	}
	return false
}

func isOne(e Expr) bool {
	switch n := e.(type) {
	case *Integer:
		return n.Val.Cmp(big.NewInt(1)) == 0
	case *Float:
		return n.Val == 1
	}
	return false
}

func isMinusOne(e Expr) bool {
	switch n := e.(type) {
	case *Integer:
		return n.Val.Cmp(big.NewInt(-1)) == 0
	case *Float:
		return n.Val == -1
	}
	return false
}

// numAdd returns the numeric sum of two numeric expressions.
func numAdd(a, b Expr) Expr {
	if _, ok := a.(*Float); ok {
		return newFloat(toFloat(a) + toFloat(b))
	}
	if _, ok := b.(*Float); ok {
		return newFloat(toFloat(a) + toFloat(b))
	}
	ra, _ := toRat(a)
	rb, _ := toRat(b)
	return newRational(new(big.Rat).Add(ra, rb))
}

// numMul returns the numeric product of two numeric expressions.
func numMul(a, b Expr) Expr {
	if _, ok := a.(*Float); ok {
		return newFloat(toFloat(a) * toFloat(b))
	}
	if _, ok := b.(*Float); ok {
		return newFloat(toFloat(a) * toFloat(b))
	}
	ra, _ := toRat(a)
	rb, _ := toRat(b)
	return newRational(new(big.Rat).Mul(ra, rb))
}

// numNeg returns the numeric negation of a numeric expression.
func numNeg(a Expr) Expr {
	switch n := a.(type) {
	case *Integer:
		return newInteger(new(big.Int).Neg(n.Val))
	case *Rational:
		return newRational(new(big.Rat).Neg(n.Val))
	case *Float:
		return newFloat(-n.Val)
	}
	return a
}
