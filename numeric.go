package algebra

import "math"

// evalFn1 numerically evaluates the single-argument function name at v. It
// reports ok=false for names it does not recognise.
func evalFn1(name string, v float64) (float64, bool) {
	switch name {
	case "sin":
		return math.Sin(v), true
	case "cos":
		return math.Cos(v), true
	case "tan":
		return math.Tan(v), true
	case "sec":
		return 1 / math.Cos(v), true
	case "csc":
		return 1 / math.Sin(v), true
	case "cot":
		return math.Cos(v) / math.Sin(v), true
	case "asin":
		return math.Asin(v), true
	case "acos":
		return math.Acos(v), true
	case "atan":
		return math.Atan(v), true
	case "acot":
		// Branch with derivative -1/(1+v^2), continuous for all real v.
		return math.Atan2(1, v), true
	case "asec":
		return math.Acos(1 / v), true
	case "acsc":
		return math.Asin(1 / v), true
	case "sinh":
		return math.Sinh(v), true
	case "cosh":
		return math.Cosh(v), true
	case "tanh":
		return math.Tanh(v), true
	case "coth":
		return math.Cosh(v) / math.Sinh(v), true
	case "sech":
		return 1 / math.Cosh(v), true
	case "csch":
		return 1 / math.Sinh(v), true
	case "asinh":
		return math.Asinh(v), true
	case "acosh":
		return math.Acosh(v), true
	case "atanh":
		return math.Atanh(v), true
	case "exp":
		return math.Exp(v), true
	case "log":
		return math.Log(v), true
	case "sqrt":
		return math.Sqrt(v), true
	case "abs":
		return math.Abs(v), true
	case "sign":
		if v > 0 {
			return 1, true
		} else if v < 0 {
			return -1, true
		}
		return 0, true
	case "floor":
		return math.Floor(v), true
	case "ceil":
		return math.Ceil(v), true
	case "factorial":
		return math.Gamma(v + 1), true
	case "gamma":
		return math.Gamma(v), true
	case "digamma":
		return digamma(v), true
	case "erf":
		return math.Erf(v), true
	case "erfc":
		return math.Erfc(v), true
	case "re":
		return v, true
	case "im":
		return 0, true
	case "conjugate":
		return v, true
	case "arg":
		if v >= 0 {
			return 0, true
		}
		return math.Pi, true
	}
	return 0, false
}

// evalFn2 numerically evaluates the two-argument function name at (a, b).
func evalFn2(name string, a, b float64) (float64, bool) {
	switch name {
	case "atan2":
		return math.Atan2(a, b), true
	case "beta":
		la, sa := math.Lgamma(a)
		lb, sb := math.Lgamma(b)
		lab, sab := math.Lgamma(a + b)
		return float64(sa*sb*sab) * math.Exp(la+lb-lab), true
	}
	return 0, false
}

// digamma returns the digamma function psi(x) = d/dx log(Gamma(x)) for real x
// that is not a non-positive integer. It uses the recurrence
// psi(x) = psi(x+1) - 1/x to push the argument above 6 and then the standard
// asymptotic expansion.
func digamma(x float64) float64 {
	result := 0.0
	if x <= 0 && x == math.Trunc(x) {
		return math.NaN()
	}
	// Reflection for negative arguments keeps the recurrence well conditioned.
	if x < 0 {
		// psi(1-x) - psi(x) = pi*cot(pi*x).
		return digamma(1-x) - math.Pi/math.Tan(math.Pi*x)
	}
	for x < 6 {
		result -= 1 / x
		x++
	}
	f := 1 / (x * x)
	// psi(x) ~ ln(x) - 1/(2x) - 1/(12x^2) + 1/(120x^4) - 1/(252x^6) ...
	result += math.Log(x) - 1/(2*x) -
		f*(1.0/12-f*(1.0/120-f*(1.0/252)))
	return result
}
