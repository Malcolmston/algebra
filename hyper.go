package algebra

// This file adds the hyperbolic functions Sinh, Cosh, Tanh, Coth, Sech and
// Csch together with the inverse hyperbolic functions Asinh, Acosh and Atanh.

// Sinh returns the hyperbolic sine sinh(x), folding sinh(0) to 0.
func Sinh(x Expr) Expr {
	if isZero(x) {
		return Int(0)
	}
	if neg, mag := splitSign(x); neg {
		return neg1Mul(Sinh(mag)) // sinh is odd
	}
	return newFn("sinh", x)
}

// Cosh returns the hyperbolic cosine cosh(x), folding cosh(0) to 1.
func Cosh(x Expr) Expr {
	if isZero(x) {
		return Int(1)
	}
	if neg, mag := splitSign(x); neg {
		return Cosh(mag) // cosh is even
	}
	return newFn("cosh", x)
}

// Tanh returns the hyperbolic tangent tanh(x), folding tanh(0) to 0.
func Tanh(x Expr) Expr {
	if isZero(x) {
		return Int(0)
	}
	if neg, mag := splitSign(x); neg {
		return neg1Mul(Tanh(mag)) // tanh is odd
	}
	return newFn("tanh", x)
}

// Coth returns the hyperbolic cotangent coth(x) = cosh(x)/sinh(x).
func Coth(x Expr) Expr {
	if neg, mag := splitSign(x); neg {
		return neg1Mul(Coth(mag)) // coth is odd
	}
	return newFn("coth", x)
}

// Sech returns the hyperbolic secant sech(x) = 1/cosh(x), folding sech(0) to 1.
func Sech(x Expr) Expr {
	if isZero(x) {
		return Int(1)
	}
	if neg, mag := splitSign(x); neg {
		return Sech(mag) // sech is even
	}
	return newFn("sech", x)
}

// Csch returns the hyperbolic cosecant csch(x) = 1/sinh(x).
func Csch(x Expr) Expr {
	if neg, mag := splitSign(x); neg {
		return neg1Mul(Csch(mag)) // csch is odd
	}
	return newFn("csch", x)
}

// Asinh returns the inverse hyperbolic sine arsinh(x), folding asinh(0) to 0.
func Asinh(x Expr) Expr {
	if isZero(x) {
		return Int(0)
	}
	return newFn("asinh", x)
}

// Acosh returns the inverse hyperbolic cosine arcosh(x), folding acosh(1) to 0.
func Acosh(x Expr) Expr {
	if isOne(x) {
		return Int(0)
	}
	return newFn("acosh", x)
}

// Atanh returns the inverse hyperbolic tangent artanh(x), folding atanh(0)
// to 0.
func Atanh(x Expr) Expr {
	if isZero(x) {
		return Int(0)
	}
	return newFn("atanh", x)
}
