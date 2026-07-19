package fuzzy

import "math"

// MF is a membership function mapping a point of the universe of discourse to
// a membership grade. Values produced by the constructors in this package are
// always clamped to the closed interval [0, 1].
type MF func(x float64) float64

// clamp01 restricts v to the closed interval [0, 1].
func clamp01(v float64) float64 {
	if v < 0 {
		return 0
	}
	if v > 1 {
		return 1
	}
	return v
}

// Clamp01 returns v restricted to the closed interval [0, 1]. It is the
// canonical way membership grades are kept within range.
func Clamp01(v float64) float64 { return clamp01(v) }

// TriangularAt evaluates the triangular membership function with feet a and c
// and peak b at the point x. It returns 0 outside (a, c), rises linearly from a
// to b and falls linearly from b to c.
func TriangularAt(x, a, b, c float64) float64 {
	switch {
	case x == b:
		return 1
	case x <= a || x >= c:
		return 0
	case x < b:
		return (x - a) / (b - a)
	default:
		return (c - x) / (c - b)
	}
}

// Triangular returns the triangular membership function with feet a and c and
// peak b. It requires a <= b <= c.
func Triangular(a, b, c float64) MF {
	return func(x float64) float64 { return clamp01(TriangularAt(x, a, b, c)) }
}

// TrapezoidalAt evaluates the trapezoidal membership function with feet a and d
// and shoulders b and c at the point x. It rises from a to b, is 1 on [b, c]
// and falls from c to d.
func TrapezoidalAt(x, a, b, c, d float64) float64 {
	switch {
	case x >= b && x <= c:
		return 1
	case x <= a || x >= d:
		return 0
	case x < b:
		return (x - a) / (b - a)
	default:
		return (d - x) / (d - c)
	}
}

// Trapezoidal returns the trapezoidal membership function with feet a and d and
// shoulders b and c. It requires a <= b <= c <= d.
func Trapezoidal(a, b, c, d float64) MF {
	return func(x float64) float64 { return clamp01(TrapezoidalAt(x, a, b, c, d)) }
}

// GaussianAt evaluates the Gaussian membership function with center mean and
// width sigma at the point x, exp(-(x-mean)^2 / (2 sigma^2)).
func GaussianAt(x, mean, sigma float64) float64 {
	if sigma == 0 {
		if x == mean {
			return 1
		}
		return 0
	}
	d := x - mean
	return math.Exp(-(d * d) / (2 * sigma * sigma))
}

// Gaussian returns the Gaussian membership function with center mean and width
// sigma. sigma must be positive for a proper bell shape.
func Gaussian(mean, sigma float64) MF {
	return func(x float64) float64 { return clamp01(GaussianAt(x, mean, sigma)) }
}

// Gaussian2At evaluates the two-sided Gaussian membership function that uses
// (mean1, sigma1) for the left flank and (mean2, sigma2) for the right flank,
// with a flat top of 1 between the two means.
func Gaussian2At(x, mean1, sigma1, mean2, sigma2 float64) float64 {
	left := 1.0
	if x < mean1 {
		left = GaussianAt(x, mean1, sigma1)
	}
	right := 1.0
	if x > mean2 {
		right = GaussianAt(x, mean2, sigma2)
	}
	return left * right
}

// Gaussian2 returns the two-sided Gaussian membership function combining a left
// flank (mean1, sigma1) and a right flank (mean2, sigma2).
func Gaussian2(mean1, sigma1, mean2, sigma2 float64) MF {
	return func(x float64) float64 { return clamp01(Gaussian2At(x, mean1, sigma1, mean2, sigma2)) }
}

// SigmoidAt evaluates the sigmoidal membership function 1/(1+exp(-a(x-c))) at
// the point x. Positive a gives a left-to-right rising curve, negative a a
// falling one, and larger |a| a steeper transition around c.
func SigmoidAt(x, a, c float64) float64 {
	return 1 / (1 + math.Exp(-a*(x-c)))
}

// Sigmoid returns the sigmoidal membership function with slope a and inflection
// point c.
func Sigmoid(a, c float64) MF {
	return func(x float64) float64 { return clamp01(SigmoidAt(x, a, c)) }
}

// DiffSigmoidAt evaluates the difference of two sigmoids, clamped to [0, 1],
// giving an asymmetric bump: sigmoid(a1, c1) minus sigmoid(a2, c2).
func DiffSigmoidAt(x, a1, c1, a2, c2 float64) float64 {
	return clamp01(SigmoidAt(x, a1, c1) - SigmoidAt(x, a2, c2))
}

// DiffSigmoid returns the difference-of-sigmoids membership function.
func DiffSigmoid(a1, c1, a2, c2 float64) MF {
	return func(x float64) float64 { return DiffSigmoidAt(x, a1, c1, a2, c2) }
}

// ProdSigmoidAt evaluates the product of two sigmoids sigmoid(a1, c1) times
// sigmoid(a2, c2) at the point x, producing a smooth bump.
func ProdSigmoidAt(x, a1, c1, a2, c2 float64) float64 {
	return SigmoidAt(x, a1, c1) * SigmoidAt(x, a2, c2)
}

// ProdSigmoid returns the product-of-sigmoids membership function.
func ProdSigmoid(a1, c1, a2, c2 float64) MF {
	return func(x float64) float64 { return clamp01(ProdSigmoidAt(x, a1, c1, a2, c2)) }
}

// BellAt evaluates the generalized bell membership function
// 1/(1+|(x-c)/a|^(2b)) at the point x. a controls the width, b the slope of the
// flanks and c the center.
func BellAt(x, a, b, c float64) float64 {
	if a == 0 {
		if x == c {
			return 1
		}
		return 0
	}
	t := math.Abs((x - c) / a)
	return 1 / (1 + math.Pow(t, 2*b))
}

// Bell returns the generalized bell membership function with width a, slope b
// and center c.
func Bell(a, b, c float64) MF {
	return func(x float64) float64 { return clamp01(BellAt(x, a, b, c)) }
}

// SShapeAt evaluates the spline-based S-shaped membership function with feet a
// and b at the point x. It rises smoothly from 0 at a to 1 at b.
func SShapeAt(x, a, b float64) float64 {
	switch {
	case x <= a:
		return 0
	case x >= b:
		return 1
	case a == b:
		return 1
	}
	mid := (a + b) / 2
	if x <= mid {
		t := (x - a) / (b - a)
		return 2 * t * t
	}
	t := (x - b) / (b - a)
	return 1 - 2*t*t
}

// SShape returns the spline-based S-shaped membership function with feet a and
// b. It requires a <= b.
func SShape(a, b float64) MF {
	return func(x float64) float64 { return clamp01(SShapeAt(x, a, b)) }
}

// ZShapeAt evaluates the spline-based Z-shaped membership function with
// shoulders a and b at the point x. It falls smoothly from 1 at a to 0 at b and
// equals 1 - SShapeAt.
func ZShapeAt(x, a, b float64) float64 {
	return 1 - SShapeAt(x, a, b)
}

// ZShape returns the spline-based Z-shaped membership function with shoulders a
// and b. It requires a <= b.
func ZShape(a, b float64) MF {
	return func(x float64) float64 { return clamp01(ZShapeAt(x, a, b)) }
}

// PiShapeAt evaluates the Pi-shaped membership function that rises from a to b,
// stays at 1 on [b, c] and falls from c to d, at the point x.
func PiShapeAt(x, a, b, c, d float64) float64 {
	return math.Min(SShapeAt(x, a, b), ZShapeAt(x, c, d))
}

// PiShape returns the Pi-shaped membership function with parameters a <= b <= c
// <= d built from an S-shaped rising edge and a Z-shaped falling edge.
func PiShape(a, b, c, d float64) MF {
	return func(x float64) float64 { return clamp01(PiShapeAt(x, a, b, c, d)) }
}

// RampUpAt evaluates the linear rising ramp that is 0 at or below a, 1 at or
// above b and interpolates linearly between, at the point x.
func RampUpAt(x, a, b float64) float64 {
	switch {
	case x <= a:
		return 0
	case x >= b:
		return 1
	case a == b:
		return 1
	default:
		return (x - a) / (b - a)
	}
}

// RampUp returns the linear rising ramp membership function from a to b.
func RampUp(a, b float64) MF {
	return func(x float64) float64 { return clamp01(RampUpAt(x, a, b)) }
}

// RampDownAt evaluates the linear falling ramp that is 1 at or below a, 0 at or
// above b and interpolates linearly between, at the point x.
func RampDownAt(x, a, b float64) float64 {
	return 1 - RampUpAt(x, a, b)
}

// RampDown returns the linear falling ramp membership function from a to b.
func RampDown(a, b float64) MF {
	return func(x float64) float64 { return clamp01(RampDownAt(x, a, b)) }
}

// LeftShoulder returns a membership function that is 1 up to a, falls linearly
// to 0 at b and remains 0 thereafter. It is an alias for RampDown suited to
// left shoulder linguistic terms.
func LeftShoulder(a, b float64) MF { return RampDown(a, b) }

// RightShoulder returns a membership function that is 0 up to a, rises linearly
// to 1 at b and remains 1 thereafter. It is an alias for RampUp suited to right
// shoulder linguistic terms.
func RightShoulder(a, b float64) MF { return RampUp(a, b) }

// RectangularAt evaluates the crisp rectangular (interval) membership function
// that is 1 on the closed interval [a, b] and 0 elsewhere, at the point x.
func RectangularAt(x, a, b float64) float64 {
	if x >= a && x <= b {
		return 1
	}
	return 0
}

// Rectangular returns the crisp rectangular membership function equal to 1 on
// [a, b] and 0 outside.
func Rectangular(a, b float64) MF {
	return func(x float64) float64 { return RectangularAt(x, a, b) }
}

// SingletonAt returns 1 when x equals c and 0 otherwise, the membership of a
// fuzzy singleton (used to fuzzify a crisp input).
func SingletonAt(x, c float64) float64 {
	if x == c {
		return 1
	}
	return 0
}

// Singleton returns the fuzzy singleton membership function concentrated at c.
func Singleton(c float64) MF {
	return func(x float64) float64 { return SingletonAt(x, c) }
}

// Constant returns a membership function equal to the constant grade v
// (clamped to [0, 1]) everywhere.
func Constant(v float64) MF {
	c := clamp01(v)
	return func(x float64) float64 { return c }
}
