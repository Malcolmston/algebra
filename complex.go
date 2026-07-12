package algebra

import (
	"errors"
	"fmt"
	"math"
	"math/big"
	"math/cmplx"
)

// This file implements complex-number support: integer powers of the imaginary
// unit I, clean Euler identities for exp, the component functions Conjugate,
// Re, Im and Arg, and a complex128-backed numeric evaluator [Evalc].
//
// Free symbols are treated as real by Re, Im, Conjugate and Arg; this matches
// the common small-CAS convention and keeps the decomposition of a+b*I exact.

// reduceIPow returns I^n for an integer n, using I^4 = 1.
func reduceIPow(n *big.Int) Expr {
	m := new(big.Int).Mod(n, big.NewInt(4)) // Euclidean: result in [0,4)
	switch m.Int64() {
	case 0:
		return Int(1)
	case 1:
		return I
	case 2:
		return Int(-1)
	default: // 3
		return neg(I)
	}
}

// eulerExp recognises exp(c*I*Pi) for a rational c whose double is an integer
// (that is, integer or half-integer multiples of Pi), returning the exact
// value via exp(c*I*Pi) = I^(2c). It returns nil when no clean value applies.
func eulerExp(x Expr) Expr {
	c, ok := iPiCoeff(x)
	if !ok {
		return nil
	}
	two := new(big.Rat).Mul(c, big.NewRat(2, 1))
	if !two.IsInt() {
		return nil
	}
	return reduceIPow(two.Num())
}

// iPiCoeff reports whether e is a rational multiple of I*Pi and returns the
// rational coefficient.
func iPiCoeff(e Expr) (*big.Rat, bool) {
	coeff := big.NewRat(1, 1)
	iCount, piCount := 0, 0
	for _, f := range factorsOf(e) {
		if isImagUnit(f) {
			iCount++
			continue
		}
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
	if iCount == 1 && piCount == 1 {
		return coeff, true
	}
	return nil, false
}

// containsI reports whether e references the imaginary unit anywhere.
func containsI(e Expr) bool {
	switch x := e.(type) {
	case *Constant:
		return x.Name == "I"
	case *sum:
		for _, a := range x.args {
			if containsI(a) {
				return true
			}
		}
	case *product:
		for _, f := range x.factors {
			if containsI(f) {
				return true
			}
		}
	case *power:
		return containsI(x.base) || containsI(x.exp)
	case *fn:
		return containsI(x.arg)
	case *fn2:
		return containsI(x.arg1) || containsI(x.arg2)
	}
	return false
}

// splitComplex decomposes e into its real and imaginary parts re + im*I,
// assuming any free symbol is real. It reports ok=false when the imaginary unit
// appears in a way that cannot be reduced to that linear form (for example
// inside a function argument).
func splitComplex(e Expr) (re, im Expr, ok bool) {
	ex := Simplify(Expand(e))
	var reTerms, imTerms []Expr
	for _, t := range termsOf(ex) {
		var others []Expr
		iCount := 0
		for _, f := range factorsOf(t) {
			if isImagUnit(f) {
				iCount++
				continue
			}
			others = append(others, f)
		}
		rest := Mul(others...)
		if containsI(rest) {
			return nil, nil, false
		}
		switch iCount {
		case 0:
			reTerms = append(reTerms, rest)
		case 1:
			imTerms = append(imTerms, rest)
		default:
			return nil, nil, false
		}
	}
	return Add(reTerms...), Add(imTerms...), true
}

// Conjugate returns the complex conjugate of e. For a+b*I it returns a-b*I;
// free symbols are treated as real.
func Conjugate(e Expr) Expr {
	if isImagUnit(e) {
		return neg(I)
	}
	re, im, ok := splitComplex(e)
	if !ok {
		return newFn("conjugate", e)
	}
	if isZero(im) {
		return re
	}
	return Add(re, Mul(neg(im), I))
}

// Re returns the real part of e, treating free symbols as real.
func Re(e Expr) Expr {
	re, _, ok := splitComplex(e)
	if !ok {
		return newFn("re", e)
	}
	return re
}

// Im returns the imaginary part of e (the real coefficient of I), treating free
// symbols as real.
func Im(e Expr) Expr {
	_, im, ok := splitComplex(e)
	if !ok {
		return newFn("im", e)
	}
	return im
}

// Arg returns the argument (phase) of e in the range (-Pi, Pi]. Numeric inputs
// fold to an exact angle where possible and otherwise to a [Float].
func Arg(e Expr) Expr {
	re, im, ok := splitComplex(e)
	if !ok {
		return newFn("arg", e)
	}
	if isZero(im) {
		if isNum(re) {
			if numSign(re) >= 0 {
				return Int(0)
			}
			return Pi
		}
		return newFn("arg", e)
	}
	return Atan2(im, re)
}

// Evalc numerically evaluates a symbol-free expression to a complex128,
// interpreting I as the imaginary unit. Use it when an expression may take a
// complex value (for example sqrt(-1) or exp(I*Pi/3)).
func Evalc(e Expr) (complex128, error) { return evalc(e, nil) }

// EvalComplex is Evalc with an environment binding symbols to complex values.
func EvalComplex(e Expr, env map[string]complex128) (complex128, error) {
	return evalc(e, env)
}

func evalc(e Expr, env map[string]complex128) (complex128, error) {
	switch x := e.(type) {
	case *Integer, *Rational, *Float:
		return complex(toFloat(e), 0), nil
	case *Constant:
		if x.Name == "I" {
			return 1i, nil
		}
		return complex(x.num, 0), nil
	case *Symbol:
		if env != nil {
			if v, ok := env[x.Name]; ok {
				return v, nil
			}
		}
		return 0, fmt.Errorf("algebra: unbound symbol %q", x.Name)
	case *sum:
		total := complex(0, 0)
		for _, a := range x.args {
			v, err := evalc(a, env)
			if err != nil {
				return 0, err
			}
			total += v
		}
		return total, nil
	case *product:
		total := complex(1, 0)
		for _, f := range x.factors {
			v, err := evalc(f, env)
			if err != nil {
				return 0, err
			}
			total *= v
		}
		return total, nil
	case *power:
		b, err := evalc(x.base, env)
		if err != nil {
			return 0, err
		}
		p, err := evalc(x.exp, env)
		if err != nil {
			return 0, err
		}
		return cmplx.Pow(b, p), nil
	case *fn:
		v, err := evalc(x.arg, env)
		if err != nil {
			return 0, err
		}
		if r, ok := evalcFn1(x.name, v); ok {
			return r, nil
		}
	case *fn2:
		a, err := evalc(x.arg1, env)
		if err != nil {
			return 0, err
		}
		b, err := evalc(x.arg2, env)
		if err != nil {
			return 0, err
		}
		if r, ok := evalcFn2(x.name, a, b); ok {
			return r, nil
		}
	}
	return 0, errors.New("algebra: cannot evaluate " + e.String())
}

// evalcFn1 evaluates a single-argument function on a complex argument. Names
// without a native complex implementation fall back to the real evaluator when
// the argument is (numerically) real.
func evalcFn1(name string, v complex128) (complex128, bool) {
	switch name {
	case "sin":
		return cmplx.Sin(v), true
	case "cos":
		return cmplx.Cos(v), true
	case "tan":
		return cmplx.Tan(v), true
	case "sec":
		return 1 / cmplx.Cos(v), true
	case "csc":
		return 1 / cmplx.Sin(v), true
	case "cot":
		return cmplx.Cos(v) / cmplx.Sin(v), true
	case "asin":
		return cmplx.Asin(v), true
	case "acos":
		return cmplx.Acos(v), true
	case "atan":
		return cmplx.Atan(v), true
	case "sinh":
		return cmplx.Sinh(v), true
	case "cosh":
		return cmplx.Cosh(v), true
	case "tanh":
		return cmplx.Tanh(v), true
	case "asinh":
		return cmplx.Asinh(v), true
	case "acosh":
		return cmplx.Acosh(v), true
	case "atanh":
		return cmplx.Atanh(v), true
	case "exp":
		return cmplx.Exp(v), true
	case "log":
		return cmplx.Log(v), true
	case "sqrt":
		return cmplx.Sqrt(v), true
	case "abs":
		return complex(cmplx.Abs(v), 0), true
	case "conjugate":
		return cmplx.Conj(v), true
	case "re":
		return complex(real(v), 0), true
	case "im":
		return complex(imag(v), 0), true
	case "arg":
		return complex(cmplx.Phase(v), 0), true
	}
	// Real-only functions: fall back when the argument is real.
	if imag(v) == 0 {
		if r, ok := evalFn1(name, real(v)); ok {
			return complex(r, 0), true
		}
	}
	return 0, false
}

func evalcFn2(name string, a, b complex128) (complex128, bool) {
	if imag(a) == 0 && imag(b) == 0 {
		if r, ok := evalFn2(name, real(a), real(b)); ok {
			return complex(r, 0), true
		}
	}
	return 0, false
}

// evalfComplexFallback tries a complex evaluation and returns its real part
// when the imaginary component is negligible.
func evalfComplexFallback(e Expr) (float64, error) {
	z, err := evalc(e, nil)
	if err != nil {
		return 0, err
	}
	if math.Abs(imag(z)) > 1e-12*(1+math.Abs(real(z))) {
		return 0, fmt.Errorf("algebra: %s evaluates to a non-real value; use Evalc", e)
	}
	return real(z), nil
}
