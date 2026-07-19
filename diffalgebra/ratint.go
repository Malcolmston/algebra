package diffalgebra

import (
	"fmt"
	"math"
	"math/big"
	"strings"
)

// LogTerm is a summand Coeff * log(Arg) of the logarithmic part of an integral,
// where Coeff is a rational residue and Arg is a monic polynomial.
type LogTerm struct {
	Coeff *big.Rat
	Arg   Poly
}

// String renders the log term as "c*log(arg)".
func (t LogTerm) String() string {
	return fmt.Sprintf("%s*log(%s)", t.Coeff.RatString(), t.Arg.String())
}

// RationalIntegral is the result of integrating a rational function: an
// elementary rational part plus a sum of logarithm terms. When some residue is
// not rational, AllResiduesRational is false and ResidueResultant holds the
// Rothstein-Trager resultant whose roots are the residues.
type RationalIntegral struct {
	Rational            RatFunc
	Logs                []LogTerm
	ResidueResultant    Poly
	AllResiduesRational bool
}

// String renders the integral as its rational part followed by its log terms.
func (r RationalIntegral) String() string {
	var parts []string
	if !r.Rational.IsZero() {
		parts = append(parts, r.Rational.String())
	}
	for _, lt := range r.Logs {
		parts = append(parts, lt.String())
	}
	if len(parts) == 0 {
		return "0"
	}
	return strings.Join(parts, " + ")
}

// EvalFloat evaluates the integral numerically at x, using log|Arg| for the
// logarithmic terms. It is intended for cross-checking against numerical
// quadrature.
func (r RationalIntegral) EvalFloat(x float64) float64 {
	acc := r.Rational.EvalFloat(x)
	for _, lt := range r.Logs {
		acc += RatToFloat(lt.Coeff) * math.Log(math.Abs(lt.Arg.EvalFloat(x)))
	}
	return acc
}

// primePowerParts CRT-splits num/den into parts A_i / F_i^{e_i}, one per
// square-free factor F_i^{e_i} of den, with deg A_i < deg F_i^{e_i}.
type primePart struct {
	factor Poly
	power  int
	numer  Poly
}

func primePowerParts(num, den Poly) []primePart {
	sqf := den.SquareFreeFactorization()
	var parts []primePart
	for _, sf := range sqf {
		parts = append(parts, primePart{factor: sf.Factor, power: sf.Mult})
	}
	for i := range parts {
		mod := parts[i].factor.Pow(parts[i].power)
		Ni, _, _ := den.DivMod(mod)
		inv := invMod(Ni, mod)
		ai := num.Mul(inv)
		_, r, _ := ai.DivMod(mod)
		parts[i].numer = r
	}
	return parts
}

// hermiteFactor reduces A / F^n (F monic square-free, n >= 1) to a rational
// part plus a remaining numerator over F (the logarithmic remainder A/F^1).
func hermiteFactor(A, F Poly, n int) (rational RatFunc, logNum Poly) {
	rational = ZeroRatFunc()
	Bm := A
	_, S, T := F.ExtendedGCD(F.Derivative()) // S*F + T*F' = 1 (F square-free)
	for m := n; m >= 2; m-- {
		invm := ratInv(ratInt(int64(m - 1)))
		W := Bm.Mul(T)
		fpow := F.Pow(m - 1)
		rational = rational.Sub(mustRat(W.ScalarMul(invm), fpow))
		Bm = Bm.Mul(S).Add(W.Derivative().ScalarMul(invm))
	}
	q, r, _ := Bm.DivMod(F)
	rational = rational.Add(RatFuncFromPoly(q.Integral()))
	return rational, r
}

// hermiteReduceProper reduces a proper rational function to a rational part
// plus a remaining proper integrand whose denominator is square-free.
func hermiteReduceProper(f RatFunc) (rational RatFunc, remaining RatFunc) {
	rational = ZeroRatFunc()
	remaining = ZeroRatFunc()
	for _, part := range primePowerParts(f.num, f.den) {
		rat, logNum := hermiteFactor(part.numer, part.factor, part.power)
		rational = rational.Add(rat)
		if !logNum.IsZero() {
			remaining = remaining.Add(mustRat(logNum, part.factor))
		}
	}
	return rational, remaining
}

// HermiteReduce performs Hermite reduction of the rational function f, returning
// the elementary rational part of its integral together with a remaining
// integrand whose denominator is square-free (so that its integral is a sum of
// logarithms). It always succeeds.
func HermiteReduce(f RatFunc) (rational RatFunc, remaining RatFunc) {
	q, proper := f.PolynomialPart()
	rational = RatFuncFromPoly(q.Integral())
	if !proper.IsZero() {
		rat, rem := hermiteReduceProper(proper)
		rational = rational.Add(rat)
		remaining = rem
	} else {
		remaining = ZeroRatFunc()
	}
	return rational, remaining
}

// RothsteinTragerResultant returns the Rothstein-Trager resultant
// R(z) = Res_x(D, C - z D') for the square-free proper integrand C/D, whose
// rational roots are the residues of the logarithmic part.
func RothsteinTragerResultant(f RatFunc) Poly {
	C := f.num
	D := f.den
	Dp := D.Derivative()
	N := D.Degree()
	if N < 1 {
		return ZeroPoly()
	}
	xs := make([]*big.Rat, N+1)
	ys := make([]*big.Rat, N+1)
	for j := 0; j <= N; j++ {
		z0 := ratInt(int64(j))
		A := C.Sub(Dp.ScalarMul(z0))
		xs[j] = z0
		ys[j] = resultant(D.Clone(), A)
	}
	return lagrangeInterp(xs, ys)
}

// lagrangeInterp returns the interpolating polynomial through the points
// (xs[i], ys[i]); the xs must be distinct.
func lagrangeInterp(xs, ys []*big.Rat) Poly {
	n := len(xs)
	acc := ZeroPoly()
	for i := 0; i < n; i++ {
		basis := OnePoly()
		denom := ratInt(1)
		for j := 0; j < n; j++ {
			if j == i {
				continue
			}
			basis = basis.Mul(NewPoly(ratNeg(xs[j]), ratInt(1)))
			denom = ratMul(denom, ratSub(xs[i], xs[j]))
		}
		coeff := ratDiv(ys[i], denom)
		acc = acc.Add(basis.ScalarMul(coeff))
	}
	return acc
}

// logPart computes the logarithmic part of the integral of the square-free
// proper integrand f = C/D via the Rothstein-Trager method. It returns the
// resultant, the rational log terms found, and whether every residue was
// rational.
func logPart(f RatFunc) (Poly, []LogTerm, bool) {
	if f.IsZero() {
		return ZeroPoly(), nil, true
	}
	C := f.num
	D := f.den
	Dp := D.Derivative()
	res := RothsteinTragerResultant(f)
	roots := res.RationalRoots()
	var logs []LogTerm
	degCovered := 0
	for _, c := range roots {
		arg := C.Sub(Dp.ScalarMul(c)).GCD(D)
		if arg.Degree() < 1 {
			continue
		}
		arg = arg.Monic()
		logs = append(logs, LogTerm{Coeff: cloneRat(c), Arg: arg})
		degCovered += arg.Degree()
	}
	allRational := degCovered == D.Degree()
	return res, logs, allRational
}

// IntegrateRational returns the elementary integral of the rational function f,
// as a rational part plus a sum of logarithms with rational coefficients. Every
// rational function has an elementary integral; when a residue happens to be
// irrational the corresponding logarithm is omitted from Logs and
// AllResiduesRational is reported false, with ResidueResultant describing the
// remaining residues.
func IntegrateRational(f RatFunc) (RationalIntegral, error) {
	rational, remaining := HermiteReduce(f)
	out := RationalIntegral{Rational: rational, AllResiduesRational: true}
	if remaining.IsZero() {
		return out, nil
	}
	res, logs, allRat := logPart(remaining)
	out.Logs = logs
	out.ResidueResultant = res
	out.AllResiduesRational = allRat
	return out, nil
}
