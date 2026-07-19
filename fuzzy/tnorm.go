package fuzzy

import "math"

// TNorm is a triangular norm, a binary operator on [0, 1] used as the fuzzy
// logical conjunction (AND) and as fuzzy set intersection. A t-norm is
// commutative, associative, monotone and has 1 as identity.
type TNorm func(a, b float64) float64

// TConorm is a triangular conorm (s-norm), a binary operator on [0, 1] used as
// the fuzzy logical disjunction (OR) and as fuzzy set union. A t-conorm is
// commutative, associative, monotone and has 0 as identity.
type TConorm func(a, b float64) float64

// TNormMin is the minimum (Godel) t-norm min(a, b), the largest t-norm.
func TNormMin(a, b float64) float64 { return math.Min(a, b) }

// TNormProduct is the algebraic product t-norm a*b.
func TNormProduct(a, b float64) float64 { return a * b }

// TNormLukasiewicz is the Lukasiewicz t-norm max(0, a+b-1).
func TNormLukasiewicz(a, b float64) float64 { return math.Max(0, a+b-1) }

// TNormDrastic is the drastic product t-norm, the smallest t-norm: it equals a
// when b is 1, b when a is 1 and 0 otherwise.
func TNormDrastic(a, b float64) float64 {
	if a == 1 {
		return b
	}
	if b == 1 {
		return a
	}
	return 0
}

// TNormEinstein is the Einstein product t-norm a*b / (2 - (a+b-a*b)).
func TNormEinstein(a, b float64) float64 {
	return (a * b) / (2 - (a + b - a*b))
}

// TNormNilpotentMin is the nilpotent minimum t-norm: min(a, b) when a+b > 1 and
// 0 otherwise.
func TNormNilpotentMin(a, b float64) float64 {
	if a+b > 1 {
		return math.Min(a, b)
	}
	return 0
}

// TNormHamacher returns the Hamacher product t-norm with parameter gamma >= 0:
// a*b / (gamma + (1-gamma)(a+b-a*b)). gamma = 1 recovers the algebraic product
// and gamma = 2 the Einstein product.
func TNormHamacher(gamma float64) TNorm {
	return func(a, b float64) float64 {
		if a == 0 && b == 0 {
			return 0
		}
		den := gamma + (1-gamma)*(a+b-a*b)
		if den == 0 {
			return 0
		}
		return (a * b) / den
	}
}

// TConormMax is the maximum (Godel) t-conorm max(a, b), the smallest t-conorm.
func TConormMax(a, b float64) float64 { return math.Max(a, b) }

// TConormProbabilistic is the probabilistic sum t-conorm a + b - a*b, dual to
// the algebraic product.
func TConormProbabilistic(a, b float64) float64 { return a + b - a*b }

// TConormLukasiewicz is the Lukasiewicz (bounded sum) t-conorm min(1, a+b).
func TConormLukasiewicz(a, b float64) float64 { return math.Min(1, a+b) }

// TConormDrastic is the drastic sum t-conorm, the largest t-conorm: it equals a
// when b is 0, b when a is 0 and 1 otherwise.
func TConormDrastic(a, b float64) float64 {
	if a == 0 {
		return b
	}
	if b == 0 {
		return a
	}
	return 1
}

// TConormEinstein is the Einstein sum t-conorm (a+b) / (1 + a*b).
func TConormEinstein(a, b float64) float64 {
	return (a + b) / (1 + a*b)
}

// TConormNilpotentMax is the nilpotent maximum t-conorm: max(a, b) when a+b < 1
// and 1 otherwise.
func TConormNilpotentMax(a, b float64) float64 {
	if a+b < 1 {
		return math.Max(a, b)
	}
	return 1
}

// TConormHamacher returns the Hamacher sum t-conorm with parameter gamma >= 0:
// (a + b + (gamma-1)*a*b) / (1 + gamma*a*b). gamma = 1 recovers the
// probabilistic sum and gamma = 2 the Einstein sum.
func TConormHamacher(gamma float64) TConorm {
	return func(a, b float64) float64 {
		den := 1 + gamma*a*b
		if den == 0 {
			return 1
		}
		return (a + b + (gamma-1)*a*b) / den
	}
}

// DualTConorm returns the t-conorm dual to the t-norm t under the standard
// complement, s(a, b) = 1 - t(1-a, 1-b).
func DualTConorm(t TNorm) TConorm {
	return func(a, b float64) float64 { return 1 - t(1-a, 1-b) }
}

// DualTNorm returns the t-norm dual to the t-conorm s under the standard
// complement, t(a, b) = 1 - s(1-a, 1-b).
func DualTNorm(s TConorm) TNorm {
	return func(a, b float64) float64 { return 1 - s(1-a, 1-b) }
}

// TNormN aggregates the values vs with the t-norm t, folding left and starting
// from the identity 1. An empty slice yields 1.
func TNormN(t TNorm, vs ...float64) float64 {
	acc := 1.0
	for _, v := range vs {
		acc = t(acc, v)
	}
	return acc
}

// TConormN aggregates the values vs with the t-conorm s, folding left and
// starting from the identity 0. An empty slice yields 0.
func TConormN(s TConorm, vs ...float64) float64 {
	acc := 0.0
	for _, v := range vs {
		acc = s(acc, v)
	}
	return acc
}

// ComplementStandard is the standard (Zadeh) fuzzy complement 1 - a.
func ComplementStandard(a float64) float64 { return 1 - a }

// ComplementSugeno returns the Sugeno class complement (1-a)/(1+lambda*a) with
// parameter lambda > -1. lambda = 0 recovers the standard complement.
func ComplementSugeno(lambda float64) func(float64) float64 {
	return func(a float64) float64 {
		return (1 - a) / (1 + lambda*a)
	}
}

// ComplementYager returns the Yager class complement (1 - a^w)^(1/w) with
// parameter w > 0. w = 1 recovers the standard complement.
func ComplementYager(w float64) func(float64) float64 {
	return func(a float64) float64 {
		return math.Pow(1-math.Pow(a, w), 1/w)
	}
}
