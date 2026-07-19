package diffalgebra

import (
	"fmt"
	"math"
	"math/big"
)

// ReduceToNormalForm converts the second-order linear ODE
// a2 y” + a1 y' + a0 y = 0 into the reduced normal form z” = r z, returning r
// as a rational function. The substitution is y = z * exp(-1/2 integral(a1/a2)).
// It returns ErrDivByZero when a2 is the zero polynomial.
func ReduceToNormalForm(a2, a1, a0 Poly) (RatFunc, error) {
	if a2.IsZero() {
		return ZeroRatFunc(), ErrDivByZero
	}
	p, err := NewRatFunc(a1, a2)
	if err != nil {
		return ZeroRatFunc(), err
	}
	q, err := NewRatFunc(a0, a2)
	if err != nil {
		return ZeroRatFunc(), err
	}
	// r = p^2/4 + p'/2 - q
	quarter := big.NewRat(1, 4)
	half := big.NewRat(1, 2)
	r := p.Mul(p).ScalarMul(quarter).Add(p.Derivative().ScalarMul(half)).Sub(q)
	return r, nil
}

// KovacicPole describes a rational pole c of r together with its order and the
// leading Laurent coefficient b (the coefficient of 1/(x-c)^order).
type KovacicPole struct {
	C     *big.Rat
	Order int
	B     *big.Rat
}

// KovacicResult is the outcome of the Kovacic algorithm for y” = r y. When
// Found is true a Liouvillian solution y = P(x) exp(integral Omega) was
// constructed (Case 1), with ExpLogs giving integral(Omega) as a sum of
// logarithms.
type KovacicResult struct {
	Found   bool
	Case    int
	R       RatFunc
	Omega   RatFunc
	P       Poly
	ExpLogs []LogTerm
}

// SolutionString renders the constructed solution y = P * exp(integral Omega).
func (k KovacicResult) SolutionString() string {
	if !k.Found {
		return "no Liouvillian solution found"
	}
	pref := ""
	if !(k.P.IsConstant() && k.P.LeadingCoeff().Cmp(ratInt(1)) == 0) {
		pref = "(" + k.P.String() + ")*"
	}
	arg := ""
	for i, lt := range k.ExpLogs {
		if i > 0 {
			arg += " + "
		}
		arg += lt.String()
	}
	if arg == "" {
		return pref + "1"
	}
	return pref + "exp(" + arg + ")"
}

// EvalFloat evaluates the constructed solution at x, using |x-c| inside the
// logarithms; it returns 0 when no solution was found.
func (k KovacicResult) EvalFloat(x float64) float64 {
	if !k.Found {
		return 0
	}
	expArg := 0.0
	for _, lt := range k.ExpLogs {
		expArg += RatToFloat(lt.Coeff) * math.Log(math.Abs(lt.Arg.EvalFloat(x)))
	}
	return k.P.EvalFloat(x) * math.Exp(expArg)
}

// kovacicPoles returns the rational poles of r with orders and leading Laurent
// coefficients, and reports whether the denominator splits completely over Q
// into factors of order at most two.
func kovacicPoles(r RatFunc) ([]KovacicPole, bool) {
	t := r.Den()
	if t.IsConstant() {
		return nil, true
	}
	roots := t.RationalRoots()
	var poles []KovacicPole
	deg := 0
	for _, c := range roots {
		// multiplicity of c in t
		fac := NewPoly(ratNeg(c), ratInt(1)) // (x - c)
		mult := 0
		cur := t
		for {
			q, rem, _ := cur.DivMod(fac)
			if !rem.IsZero() {
				break
			}
			cur = q
			mult++
		}
		if mult == 0 {
			continue
		}
		deg += mult
		// leading Laurent coefficient b = ((x-c)^order * r)(c)
		facPow := fac.Pow(mult)
		rr := r.Mul(RatFuncFromPoly(facPow))
		b, err := rr.EvalRat(c)
		if err != nil {
			return nil, false
		}
		poles = append(poles, KovacicPole{C: cloneRat(c), Order: mult, B: b})
	}
	if deg != t.Degree() {
		return nil, false // not fully split over Q
	}
	for _, p := range poles {
		if p.Order > 2 {
			return nil, false
		}
	}
	return poles, true
}

// alphaPair returns the two Kovacic exponents alpha^+ and alpha^- for a pole or
// for infinity given the leading coefficient b and pole order. It reports
// success only when the exponents are rational.
func alphaPair(order int, b *big.Rat) (aplus, aminus *big.Rat, ok bool) {
	if order == 1 {
		return ratInt(1), ratInt(1), true
	}
	// order 2: alpha^± = 1/2 ± 1/2 sqrt(1+4b)
	disc := ratAdd(ratInt(1), ratMul(ratInt(4), b))
	s, sq := ratSqrt(disc)
	if !sq {
		return nil, nil, false
	}
	half := big.NewRat(1, 2)
	aplus = ratAdd(half, ratMul(half, s))
	aminus = ratSub(half, ratMul(half, s))
	return aplus, aminus, true
}

// Kovacic runs Case 1 of the Kovacic algorithm on y” = r y, attempting to
// construct a Liouvillian solution y = P exp(integral omega). It handles
// rational poles of order at most two and order at infinity at least two. When
// no such solution is found it returns a result with Found == false and a nil
// error; structural obstructions (non-rational poles, unsupported orders)
// likewise yield Found == false.
func Kovacic(r RatFunc) (KovacicResult, error) {
	res := KovacicResult{R: r, Case: 1}
	if r.IsZero() {
		return res, nil
	}
	poles, ok := kovacicPoles(r)
	if !ok {
		return res, nil
	}
	s := r.Num()
	t := r.Den()
	oInf := t.Degree() - s.Degree()
	if oInf < 2 {
		return res, nil // unsupported order at infinity for this Case 1 variant
	}
	// exponents at infinity
	var aInfPlus, aInfMinus *big.Rat
	if oInf > 2 {
		aInfPlus, aInfMinus = ratInt(0), ratInt(1)
	} else { // oInf == 2, b = lc(s)/lc(t)
		bInf := ratDiv(s.LeadingCoeff(), t.LeadingCoeff())
		ap, am, sok := alphaPair(2, bInf)
		if !sok {
			return res, nil
		}
		aInfPlus, aInfMinus = ap, am
	}
	// exponents at each pole
	type poleAlpha struct {
		c            *big.Rat
		aplus, minus *big.Rat
	}
	pa := make([]poleAlpha, len(poles))
	for i, p := range poles {
		ap, am, sok := alphaPair(p.Order, p.B)
		if !sok {
			return res, nil
		}
		pa[i] = poleAlpha{c: p.C, aplus: ap, minus: am}
	}
	np := len(poles)
	if np > 16 {
		return res, nil
	}
	// enumerate sign of infinity (2) and each pole (2^np)
	for infSign := 0; infSign < 2; infSign++ {
		aInf := aInfPlus
		if infSign == 1 {
			aInf = aInfMinus
		}
		for mask := 0; mask < (1 << np); mask++ {
			// compute d = aInf - sum alpha_c
			d := cloneRat(aInf)
			for i := 0; i < np; i++ {
				a := pa[i].aplus
				if mask&(1<<i) != 0 {
					a = pa[i].minus
				}
				d = ratSub(d, a)
			}
			dInt, isInt := ratIsInteger(d)
			if !isInt || dInt < 0 {
				continue
			}
			// build omega = sum alpha_c/(x-c)
			omega := ZeroRatFunc()
			var logs []LogTerm
			for i := 0; i < np; i++ {
				a := pa[i].aplus
				if mask&(1<<i) != 0 {
					a = pa[i].minus
				}
				if ratZero(a) {
					continue
				}
				fac := NewPoly(ratNeg(pa[i].c), ratInt(1)) // x - c
				omega = omega.Add(mustRat(ConstPoly(a), fac))
				logs = append(logs, LogTerm{Coeff: cloneRat(a), Arg: fac.Monic()})
			}
			P, found := solveKovacicP(omega, r, dInt)
			if found {
				res.Found = true
				res.Omega = omega
				res.P = P
				res.ExpLogs = logs
				return res, nil
			}
		}
	}
	return res, nil
}

// solveKovacicP searches for a monic polynomial P of degree d satisfying
// P” + 2 omega P' + (omega' + omega^2 - r) P = 0, returning it when found.
func solveKovacicP(omega, r RatFunc, d int) (Poly, bool) {
	coeff := omega.Derivative().Add(omega.Mul(omega)).Sub(r) // omega'+omega^2-r
	// L(x^i) as a RatFunc for i = 0..d
	lval := func(i int) RatFunc {
		mon := Monomial(ratInt(1), i)
		second := RatFuncFromPoly(mon.Derivative().Derivative())
		first := omega.ScalarMul(ratInt(2)).Mul(RatFuncFromPoly(mon.Derivative()))
		zeroth := coeff.Mul(RatFuncFromPoly(mon))
		return second.Add(first).Add(zeroth)
	}
	if d == 0 {
		return OnePoly(), lval(0).IsZero()
	}
	L := make([]RatFunc, d+1)
	lcd := OnePoly()
	for i := 0; i <= d; i++ {
		L[i] = lval(i)
		lcd = lcd.Mul(L[i].Den())
	}
	polys := make([]Poly, d+1)
	maxDeg := 0
	for i := 0; i <= d; i++ {
		pr := L[i].Mul(RatFuncFromPoly(lcd))
		polys[i] = pr.Num() // exact polynomial
		if polys[i].Degree() > maxDeg {
			maxDeg = polys[i].Degree()
		}
	}
	// sum_{i=0}^{d-1} p_i polys[i] = -polys[d]
	A := make([][]*big.Rat, maxDeg+1)
	rhs := make([]*big.Rat, maxDeg+1)
	for k := 0; k <= maxDeg; k++ {
		A[k] = make([]*big.Rat, d)
		for j := 0; j < d; j++ {
			A[k][j] = polys[j].Coeff(k)
		}
		rhs[k] = ratNeg(polys[d].Coeff(k))
	}
	sol, ok := solveRatSystem(A, rhs)
	if !ok {
		return ZeroPoly(), false
	}
	coeffs := make([]*big.Rat, d+1)
	for j := 0; j < d; j++ {
		coeffs[j] = sol[j]
	}
	coeffs[d] = ratInt(1)
	return NewPoly(coeffs...), true
}

// KovacicSolveSecondOrder is a convenience wrapper that reduces the general
// second-order equation a2 y” + a1 y' + a0 y = 0 to normal form and runs the
// Kovacic algorithm, returning both the reduced r and the result.
func KovacicSolveSecondOrder(a2, a1, a0 Poly) (RatFunc, KovacicResult, error) {
	r, err := ReduceToNormalForm(a2, a1, a0)
	if err != nil {
		return ZeroRatFunc(), KovacicResult{}, err
	}
	kr, err := Kovacic(r)
	return r, kr, err
}

// String renders a pole for diagnostics.
func (p KovacicPole) String() string {
	return fmt.Sprintf("pole c=%s order=%d b=%s", p.C.RatString(), p.Order, p.B.RatString())
}
