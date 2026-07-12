package algebra

import "math/big"

// This file completes the trigonometric surface: the reciprocal functions
// Sec, Csc and Cot, the inverse functions Asin, Acos, Atan, Acot, Asec and
// Acsc, the two-argument Atan2, exact values at the standard rational
// multiples of Pi, and the trigonometric identities applied by Simplify.

// --- reciprocal trig functions ---------------------------------------------

// Sec returns sec(x) = 1/cos(x), folding sec(0) to 1 and returning the exact
// reciprocal of cos at the standard angles where cos is non-zero.
func Sec(x Expr) Expr {
	if isZero(x) {
		return Int(1)
	}
	if _, c, ok := exactSinCosOf(x); ok && !isZero(c) {
		return Simplify(Pow(c, Int(-1)))
	}
	return newFn("sec", x)
}

// Csc returns csc(x) = 1/sin(x), returning the exact reciprocal of sin at the
// standard angles where sin is non-zero.
func Csc(x Expr) Expr {
	if s, _, ok := exactSinCosOf(x); ok && !isZero(s) {
		return Simplify(Pow(s, Int(-1)))
	}
	return newFn("csc", x)
}

// Cot returns cot(x) = cos(x)/sin(x), returning exact values at the standard
// angles where sin is non-zero.
func Cot(x Expr) Expr {
	if s, c, ok := exactSinCosOf(x); ok && !isZero(s) {
		return Simplify(Mul(c, Pow(s, Int(-1))))
	}
	return newFn("cot", x)
}

// --- inverse trig functions ------------------------------------------------

// Asin returns the inverse sine arcsin(x), folding the exact values at
// 0, ±1/2 and ±1.
func Asin(x Expr) Expr {
	switch {
	case isZero(x):
		return Int(0)
	case x.Equal(Int(1)):
		return piMul(big.NewRat(1, 2))
	case x.Equal(Int(-1)):
		return piMul(big.NewRat(-1, 2))
	case x.Equal(Rat(1, 2)):
		return piMul(big.NewRat(1, 6))
	case x.Equal(Rat(-1, 2)):
		return piMul(big.NewRat(-1, 6))
	}
	return newFn("asin", x)
}

// Acos returns the inverse cosine arccos(x), folding the exact values at
// 0, 1/2, -1 and 1.
func Acos(x Expr) Expr {
	switch {
	case isZero(x):
		return piMul(big.NewRat(1, 2))
	case x.Equal(Int(1)):
		return Int(0)
	case x.Equal(Int(-1)):
		return Pi
	case x.Equal(Rat(1, 2)):
		return piMul(big.NewRat(1, 3))
	}
	return newFn("acos", x)
}

// Atan returns the inverse tangent arctan(x), folding the exact values at
// 0 and ±1.
func Atan(x Expr) Expr {
	switch {
	case isZero(x):
		return Int(0)
	case x.Equal(Int(1)):
		return piMul(big.NewRat(1, 4))
	case x.Equal(Int(-1)):
		return piMul(big.NewRat(-1, 4))
	}
	return newFn("atan", x)
}

// Acot returns the inverse cotangent arccot(x), folding acot(0) to Pi/2 and
// acot(1) to Pi/4.
func Acot(x Expr) Expr {
	switch {
	case isZero(x):
		return piMul(big.NewRat(1, 2))
	case x.Equal(Int(1)):
		return piMul(big.NewRat(1, 4))
	}
	return newFn("acot", x)
}

// Asec returns the inverse secant arcsec(x), folding asec(1) to 0 and
// asec(-1) to Pi.
func Asec(x Expr) Expr {
	switch {
	case x.Equal(Int(1)):
		return Int(0)
	case x.Equal(Int(-1)):
		return Pi
	}
	return newFn("asec", x)
}

// Acsc returns the inverse cosecant arccsc(x), folding acsc(1) to Pi/2 and
// acsc(-1) to -Pi/2.
func Acsc(x Expr) Expr {
	switch {
	case x.Equal(Int(1)):
		return piMul(big.NewRat(1, 2))
	case x.Equal(Int(-1)):
		return piMul(big.NewRat(-1, 2))
	}
	return newFn("acsc", x)
}

// Atan2 returns the two-argument arctangent atan2(y, x), the angle of the point
// (x, y) measured from the positive x-axis in the range (-Pi, Pi]. Numeric
// arguments are folded to a [Float]; otherwise an unevaluated node is kept.
func Atan2(y, x Expr) Expr {
	if isNum(y) && isNum(x) {
		if v, ok := evalFn2("atan2", toFloat(y), toFloat(x)); ok {
			return Flt(v)
		}
	}
	if isZero(y) && isNum(x) && numSign(x) > 0 {
		return Int(0)
	}
	return newFn2("atan2", y, x)
}

// --- exact special-angle values --------------------------------------------

// piMul returns the expression r*Pi in canonical form.
func piMul(r *big.Rat) Expr { return Mul(newRational(new(big.Rat).Set(r)), Pi) }

// piCoeff reports whether e is a rational multiple of Pi and, if so, returns
// that rational coefficient.
func piCoeff(e Expr) (*big.Rat, bool) {
	switch x := e.(type) {
	case *Constant:
		if x.Name == "pi" {
			return big.NewRat(1, 1), true
		}
	case *product:
		coeff := big.NewRat(1, 1)
		piCount := 0
		for _, f := range x.factors {
			if c, ok := f.(*Constant); ok && c.Name == "pi" {
				piCount++
				continue
			}
			r, ok := toRat(f)
			if !ok {
				return nil, false
			}
			coeff.Mul(coeff, r)
		}
		if piCount == 1 {
			return coeff, true
		}
	}
	return nil, false
}

// normMod2 reduces the rational c into the half-open interval [0, 2), matching
// the 2-period of sine and cosine measured in units of Pi.
func normMod2(c *big.Rat) *big.Rat {
	// floor(c/2); big.Int.Div is floor division because Denom is positive.
	half := new(big.Rat).Quo(c, big.NewRat(2, 1))
	fl := new(big.Int).Div(half.Num(), half.Denom())
	twoFloor := new(big.Rat).SetInt(new(big.Int).Mul(fl, big.NewInt(2)))
	return new(big.Rat).Sub(c, twoFloor)
}

var (
	sqrt3half = Mul(Rat(1, 2), Sqrt(Int(3)))
	sqrt2half = Mul(Rat(1, 2), Sqrt(Int(2)))
	oneHalf   = Rat(1, 2)
)

// sinTable6 and cosTable6 hold sin and cos at k*Pi/6 for k in [0,12).
var (
	sinTable6 = []Expr{
		Int(0), oneHalf, sqrt3half, Int(1), sqrt3half, oneHalf,
		Int(0), neg(oneHalf), neg(sqrt3half), Int(-1), neg(sqrt3half), neg(oneHalf),
	}
	cosTable6 = []Expr{
		Int(1), sqrt3half, oneHalf, Int(0), neg(oneHalf), neg(sqrt3half),
		Int(-1), neg(sqrt3half), neg(oneHalf), Int(0), oneHalf, sqrt3half,
	}
	// sinTable4 and cosTable4 hold sin and cos at k*Pi/4 for k in [0,8).
	sinTable4 = []Expr{
		Int(0), sqrt2half, Int(1), sqrt2half, Int(0), neg(sqrt2half), Int(-1), neg(sqrt2half),
	}
	cosTable4 = []Expr{
		Int(1), sqrt2half, Int(0), neg(sqrt2half), Int(-1), neg(sqrt2half), Int(0), sqrt2half,
	}
)

var (
	sqrt3full  = Sqrt(Int(3))
	sqrt3third = Mul(Rat(1, 3), Sqrt(Int(3)))
	// tanTable6 holds tan at k*Pi/6 (nil where tan is undefined).
	tanTable6 = []Expr{
		Int(0), sqrt3third, sqrt3full, nil, neg(sqrt3full), neg(sqrt3third),
		Int(0), sqrt3third, sqrt3full, nil, neg(sqrt3full), neg(sqrt3third),
	}
	// tanTable4 holds tan at k*Pi/4 (nil where tan is undefined).
	tanTable4 = []Expr{Int(0), Int(1), nil, Int(-1), Int(0), Int(1), nil, Int(-1)}
)

// exactTanOf returns the exact tangent of x at a standard angle, reporting
// ok=false where x is not such an angle or tan is undefined there.
func exactTanOf(x Expr) (Expr, bool) {
	c, ok := piCoeff(x)
	if !ok {
		return nil, false
	}
	r := normMod2(c)
	if six := new(big.Rat).Mul(r, big.NewRat(6, 1)); six.IsInt() {
		k := new(big.Int).Mod(six.Num(), big.NewInt(12)).Int64()
		if t := tanTable6[k]; t != nil {
			return t, true
		}
		return nil, false
	}
	if four := new(big.Rat).Mul(r, big.NewRat(4, 1)); four.IsInt() {
		k := new(big.Int).Mod(four.Num(), big.NewInt(8)).Int64()
		if t := tanTable4[k]; t != nil {
			return t, true
		}
		return nil, false
	}
	return nil, false
}

// exactSinCosOf returns the exact sine and cosine of x when x is a rational
// multiple of Pi landing on a standard angle (a multiple of Pi/6 or Pi/4).
func exactSinCosOf(x Expr) (sin, cos Expr, ok bool) {
	c, ok := piCoeff(x)
	if !ok {
		return nil, nil, false
	}
	return exactSinCos(c)
}

func exactSinCos(c *big.Rat) (sin, cos Expr, ok bool) {
	r := normMod2(c)
	if six := new(big.Rat).Mul(r, big.NewRat(6, 1)); six.IsInt() {
		k := new(big.Int).Mod(six.Num(), big.NewInt(12)).Int64()
		return sinTable6[k], cosTable6[k], true
	}
	if four := new(big.Rat).Mul(r, big.NewRat(4, 1)); four.IsInt() {
		k := new(big.Int).Mod(four.Num(), big.NewInt(8)).Int64()
		return sinTable4[k], cosTable4[k], true
	}
	return nil, nil, false
}

// --- trigonometric identities applied by Simplify --------------------------

// simplifyTrigSum applies the Pythagorean identity to a canonical sum,
// collapsing c*sin(u)^2 + c*cos(u)^2 to c and c*cos(u)^2 - c*sin(u)^2 to
// c*cos(2u). Only exact coefficient matches are reduced, so the transformation
// is always valid.
func simplifyTrigSum(e Expr) Expr {
	s, ok := e.(*sum)
	if !ok {
		return e
	}
	sinCoeff := map[string]Expr{}
	cosCoeff := map[string]Expr{}
	uOf := map[string]Expr{}
	var others []Expr
	for _, t := range s.args {
		c, core := splitCoeff(t)
		if u, ok2 := matchSquaredTrig(core, "sin"); ok2 {
			k := u.String()
			if prev, seen := sinCoeff[k]; seen {
				sinCoeff[k] = Add(prev, c)
			} else {
				sinCoeff[k] = c
				uOf[k] = u
			}
			continue
		}
		if u, ok2 := matchSquaredTrig(core, "cos"); ok2 {
			k := u.String()
			if prev, seen := cosCoeff[k]; seen {
				cosCoeff[k] = Add(prev, c)
			} else {
				cosCoeff[k] = c
				uOf[k] = u
			}
			continue
		}
		others = append(others, t)
	}
	used := map[string]bool{}
	var result []Expr
	for k, sc := range sinCoeff {
		cc, has := cosCoeff[k]
		if !has {
			continue
		}
		switch {
		case sc.Equal(cc):
			// c*sin^2 + c*cos^2 = c.
			result = append(result, sc)
			used[k] = true
		case sc.Equal(neg(cc)):
			// c*cos^2 - c*sin^2 = c*cos(2u).
			result = append(result, Mul(cc, Cos(Mul(Int(2), uOf[k]))))
			used[k] = true
		}
	}
	for k, sc := range sinCoeff {
		if !used[k] {
			result = append(result, Mul(sc, Pow(Sin(uOf[k]), Int(2))))
		}
	}
	for k, cc := range cosCoeff {
		if !used[k] {
			result = append(result, Mul(cc, Pow(Cos(uOf[k]), Int(2))))
		}
	}
	result = append(result, others...)
	return Add(result...)
}

// matchSquaredTrig reports whether core is name(u)^2 and returns u.
func matchSquaredTrig(core Expr, name string) (Expr, bool) {
	p, ok := core.(*power)
	if !ok {
		return nil, false
	}
	n, ok := p.exp.(*Integer)
	if !ok || n.Val.Cmp(big.NewInt(2)) != 0 {
		return nil, false
	}
	f, ok := p.base.(*fn)
	if !ok || f.name != name {
		return nil, false
	}
	return f.arg, true
}

// simplifyTrigProduct rewrites products that admit a cleaner trigonometric
// form: sin(u)*cos(u) -> sin(2u)/2 and sin(u)/cos(u) -> tan(u). Both rewrites
// are exact identities.
func simplifyTrigProduct(e Expr) Expr {
	p, ok := e.(*product)
	if !ok {
		return e
	}
	// Index factors by base name for single-power trig factors.
	var rest []Expr
	sinArg := map[string]Expr{}
	cosPos := map[string]Expr{} // cos(u)^+1
	cosNeg := map[string]Expr{} // cos(u)^-1
	for _, f := range p.factors {
		if fnn, ok := f.(*fn); ok && fnn.name == "sin" {
			sinArg[fnn.arg.String()] = fnn.arg
			rest = append(rest, f)
			continue
		}
		if fnn, ok := f.(*fn); ok && fnn.name == "cos" {
			cosPos[fnn.arg.String()] = fnn.arg
			rest = append(rest, f)
			continue
		}
		if pw, ok := f.(*power); ok {
			if fnn, ok := pw.base.(*fn); ok && fnn.name == "cos" {
				if n, ok := pw.exp.(*Integer); ok && n.Val.Cmp(big.NewInt(-1)) == 0 {
					cosNeg[fnn.arg.String()] = fnn.arg
					rest = append(rest, f)
					continue
				}
			}
		}
		rest = append(rest, f)
	}
	// sin(u)/cos(u) -> tan(u)
	for k, u := range sinArg {
		if _, ok := cosNeg[k]; ok {
			out := make([]Expr, 0, len(rest))
			for _, f := range rest {
				if isTrigFactor(f, "sin", u) || isTrigInvCos(f, u) {
					continue
				}
				out = append(out, f)
			}
			out = append(out, Tan(u))
			return Mul(out...)
		}
	}
	// sin(u)*cos(u) -> sin(2u)/2
	for k, u := range sinArg {
		if _, ok := cosPos[k]; ok {
			out := make([]Expr, 0, len(rest))
			removedSin, removedCos := false, false
			for _, f := range rest {
				if !removedSin && isTrigFactor(f, "sin", u) {
					removedSin = true
					continue
				}
				if !removedCos && isTrigFactor(f, "cos", u) {
					removedCos = true
					continue
				}
				out = append(out, f)
			}
			out = append(out, Rat(1, 2), Sin(Mul(Int(2), u)))
			return Mul(out...)
		}
	}
	return e
}

// isTrigFactor reports whether f is name(u).
func isTrigFactor(f Expr, name string, u Expr) bool {
	fnn, ok := f.(*fn)
	return ok && fnn.name == name && fnn.arg.Equal(u)
}

// isTrigInvCos reports whether f is cos(u)^(-1).
func isTrigInvCos(f Expr, u Expr) bool {
	pw, ok := f.(*power)
	if !ok {
		return false
	}
	n, ok := pw.exp.(*Integer)
	if !ok || n.Val.Cmp(big.NewInt(-1)) != 0 {
		return false
	}
	fnn, ok := pw.base.(*fn)
	return ok && fnn.name == "cos" && fnn.arg.Equal(u)
}
