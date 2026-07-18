package autodiff

import "math"

// This file provides forward-mode overloads of the elementary functions for
// [Dual] numbers. Each returns f(x) with the derivative slot carrying f'(x)
// times the incoming seed, following the chain rule.

// Sin returns sin(x); its derivative is cos(x).
func Sin(x Dual) Dual {
	return autodiffDualApply(x, math.Sin(x.Val), math.Cos(x.Val))
}

// Cos returns cos(x); its derivative is -sin(x).
func Cos(x Dual) Dual {
	return autodiffDualApply(x, math.Cos(x.Val), -math.Sin(x.Val))
}

// Tan returns tan(x); its derivative is sec²(x) = 1/cos²(x).
func Tan(x Dual) Dual {
	c := math.Cos(x.Val)
	return autodiffDualApply(x, math.Tan(x.Val), 1/(c*c))
}

// Cot returns the cotangent cot(x) = cos(x)/sin(x); its derivative is
// -csc²(x) = -1/sin²(x).
func Cot(x Dual) Dual {
	s := math.Sin(x.Val)
	return autodiffDualApply(x, math.Cos(x.Val)/s, -1/(s*s))
}

// Sec returns the secant sec(x) = 1/cos(x); its derivative is
// sec(x)·tan(x).
func Sec(x Dual) Dual {
	sec := 1 / math.Cos(x.Val)
	return autodiffDualApply(x, sec, sec*math.Tan(x.Val))
}

// Csc returns the cosecant csc(x) = 1/sin(x); its derivative is
// -csc(x)·cot(x).
func Csc(x Dual) Dual {
	csc := 1 / math.Sin(x.Val)
	return autodiffDualApply(x, csc, -csc*(math.Cos(x.Val)/math.Sin(x.Val)))
}

// Asin returns arcsin(x); its derivative is 1/√(1-x²).
func Asin(x Dual) Dual {
	return autodiffDualApply(x, math.Asin(x.Val), 1/math.Sqrt(1-x.Val*x.Val))
}

// Acos returns arccos(x); its derivative is -1/√(1-x²).
func Acos(x Dual) Dual {
	return autodiffDualApply(x, math.Acos(x.Val), -1/math.Sqrt(1-x.Val*x.Val))
}

// Atan returns arctan(x); its derivative is 1/(1+x²).
func Atan(x Dual) Dual {
	return autodiffDualApply(x, math.Atan(x.Val), 1/(1+x.Val*x.Val))
}

// Atan2 returns the two-argument arctangent atan2(y, x). Its differential is
// (x·dy - y·dx)/(x²+y²), giving the correct gradient of the polar angle.
func Atan2(y, x Dual) Dual {
	den := x.Val*x.Val + y.Val*y.Val
	return Dual{
		Val: math.Atan2(y.Val, x.Val),
		Der: (x.Val*y.Der - y.Val*x.Der) / den,
	}
}

// Sinh returns sinh(x); its derivative is cosh(x).
func Sinh(x Dual) Dual {
	return autodiffDualApply(x, math.Sinh(x.Val), math.Cosh(x.Val))
}

// Cosh returns cosh(x); its derivative is sinh(x).
func Cosh(x Dual) Dual {
	return autodiffDualApply(x, math.Cosh(x.Val), math.Sinh(x.Val))
}

// Tanh returns tanh(x); its derivative is sech²(x) = 1-tanh²(x).
func Tanh(x Dual) Dual {
	t := math.Tanh(x.Val)
	return autodiffDualApply(x, t, 1-t*t)
}

// Coth returns the hyperbolic cotangent coth(x) = cosh(x)/sinh(x); its
// derivative is -csch²(x) = -1/sinh²(x).
func Coth(x Dual) Dual {
	s := math.Sinh(x.Val)
	return autodiffDualApply(x, math.Cosh(x.Val)/s, -1/(s*s))
}

// Sech returns the hyperbolic secant sech(x) = 1/cosh(x); its derivative is
// -sech(x)·tanh(x).
func Sech(x Dual) Dual {
	sech := 1 / math.Cosh(x.Val)
	return autodiffDualApply(x, sech, -sech*math.Tanh(x.Val))
}

// Csch returns the hyperbolic cosecant csch(x) = 1/sinh(x); its derivative is
// -csch(x)·coth(x).
func Csch(x Dual) Dual {
	csch := 1 / math.Sinh(x.Val)
	return autodiffDualApply(x, csch, -csch*(math.Cosh(x.Val)/math.Sinh(x.Val)))
}

// Asinh returns the inverse hyperbolic sine arsinh(x); its derivative is
// 1/√(x²+1).
func Asinh(x Dual) Dual {
	return autodiffDualApply(x, math.Asinh(x.Val), 1/math.Sqrt(x.Val*x.Val+1))
}

// Acosh returns the inverse hyperbolic cosine arcosh(x); its derivative is
// 1/√(x²-1) and is defined for x > 1.
func Acosh(x Dual) Dual {
	return autodiffDualApply(x, math.Acosh(x.Val), 1/math.Sqrt(x.Val*x.Val-1))
}

// Atanh returns the inverse hyperbolic tangent artanh(x); its derivative is
// 1/(1-x²) and is defined for |x| < 1.
func Atanh(x Dual) Dual {
	return autodiffDualApply(x, math.Atanh(x.Val), 1/(1-x.Val*x.Val))
}

// Exp returns e^x; its derivative is e^x.
func Exp(x Dual) Dual {
	e := math.Exp(x.Val)
	return autodiffDualApply(x, e, e)
}

// Exp2 returns 2^x; its derivative is ln(2)·2^x.
func Exp2(x Dual) Dual {
	e := math.Exp2(x.Val)
	return autodiffDualApply(x, e, math.Ln2*e)
}

// Expm1 returns e^x - 1 accurately for small x; its derivative is e^x.
func Expm1(x Dual) Dual {
	return autodiffDualApply(x, math.Expm1(x.Val), math.Exp(x.Val))
}

// Log returns the natural logarithm ln(x); its derivative is 1/x.
func Log(x Dual) Dual {
	return autodiffDualApply(x, math.Log(x.Val), 1/x.Val)
}

// Log2 returns the base-2 logarithm log₂(x); its derivative is 1/(x·ln2).
func Log2(x Dual) Dual {
	return autodiffDualApply(x, math.Log2(x.Val), 1/(x.Val*math.Ln2))
}

// Log10 returns the base-10 logarithm log₁₀(x); its derivative is
// 1/(x·ln10).
func Log10(x Dual) Dual {
	return autodiffDualApply(x, math.Log10(x.Val), 1/(x.Val*math.Ln10))
}

// Log1p returns ln(1+x) accurately for small x; its derivative is 1/(1+x).
func Log1p(x Dual) Dual {
	return autodiffDualApply(x, math.Log1p(x.Val), 1/(1+x.Val))
}

// Sqrt returns the square root √x; its derivative is 1/(2√x).
func Sqrt(x Dual) Dual {
	r := math.Sqrt(x.Val)
	return autodiffDualApply(x, r, 0.5/r)
}

// Cbrt returns the cube root ∛x; its derivative is 1/(3·x^(2/3)).
func Cbrt(x Dual) Dual {
	r := math.Cbrt(x.Val)
	return autodiffDualApply(x, r, 1/(3*r*r))
}

// PowReal returns x^p for a constant real exponent p; its derivative is
// p·x^(p-1).
func PowReal(x Dual, p float64) Dual {
	return autodiffDualApply(x, math.Pow(x.Val, p), p*math.Pow(x.Val, p-1))
}

// Pow returns x^y where both the base and exponent are dual numbers. Using
// x^y = exp(y·ln x), its differential is
// x^y·(y'·ln x + y·x'/x). The base must be positive.
func Pow(x, y Dual) Dual {
	v := math.Pow(x.Val, y.Val)
	return Dual{
		Val: v,
		Der: v * (y.Der*math.Log(x.Val) + y.Val*x.Der/x.Val),
	}
}

// Powf returns a^x for a constant real base a and dual exponent x. Its
// derivative is ln(a)·a^x.
func Powf(a float64, x Dual) Dual {
	v := math.Pow(a, x.Val)
	return autodiffDualApply(x, v, math.Log(a)*v)
}

// Abs returns |x|; its derivative is the sign of x and is undefined at zero,
// where this implementation returns a zero derivative.
func Abs(x Dual) Dual {
	switch {
	case x.Val > 0:
		return autodiffDualApply(x, x.Val, 1)
	case x.Val < 0:
		return autodiffDualApply(x, -x.Val, -1)
	default:
		return Dual{Val: 0, Der: 0}
	}
}

// Hypot returns √(x²+y²) computed without spurious overflow; its differential
// is (x·dx + y·dy)/hypot(x,y).
func Hypot(x, y Dual) Dual {
	h := math.Hypot(x.Val, y.Val)
	return Dual{
		Val: h,
		Der: (x.Val*x.Der + y.Val*y.Der) / h,
	}
}

// Erf returns the error function erf(x); its derivative is
// (2/√π)·e^(-x²).
func Erf(x Dual) Dual {
	d := 2 / math.SqrtPi * math.Exp(-x.Val*x.Val)
	return autodiffDualApply(x, math.Erf(x.Val), d)
}

// Erfc returns the complementary error function erfc(x) = 1-erf(x); its
// derivative is -(2/√π)·e^(-x²).
func Erfc(x Dual) Dual {
	d := -2 / math.SqrtPi * math.Exp(-x.Val*x.Val)
	return autodiffDualApply(x, math.Erfc(x.Val), d)
}

// Sigmoid returns the logistic function σ(x) = 1/(1+e^(-x)); its derivative is
// σ(x)·(1-σ(x)).
func Sigmoid(x Dual) Dual {
	s := 1 / (1 + math.Exp(-x.Val))
	return autodiffDualApply(x, s, s*(1-s))
}

// Softplus returns the smooth rectifier ln(1+e^x); its derivative is the
// logistic function σ(x).
func Softplus(x Dual) Dual {
	s := 1 / (1 + math.Exp(-x.Val))
	return autodiffDualApply(x, softplusValue(x.Val), s)
}

// softplusValue evaluates ln(1+e^x) in a numerically stable way.
func softplusValue(x float64) float64 {
	if x > 0 {
		return x + math.Log1p(math.Exp(-x))
	}
	return math.Log1p(math.Exp(x))
}

// Relu returns the rectified linear unit max(0, x); its derivative is 1 for
// positive x and 0 otherwise, with the sub-gradient at zero taken as 0.
func Relu(x Dual) Dual {
	if x.Val > 0 {
		return autodiffDualApply(x, x.Val, 1)
	}
	return Dual{Val: 0, Der: 0}
}

// Max returns the larger of x and y as a dual number, propagating the
// derivative of the selected argument. Ties select x.
func Max(x, y Dual) Dual {
	if y.Val > x.Val {
		return y
	}
	return x
}

// Min returns the smaller of x and y as a dual number, propagating the
// derivative of the selected argument. Ties select x.
func Min(x, y Dual) Dual {
	if y.Val < x.Val {
		return y
	}
	return x
}
