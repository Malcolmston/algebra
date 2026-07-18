package complexanalysis

import (
	"math"
	"math/cmplx"
)

// Sqrt returns the principal square root of z, the branch whose result lies in
// the right half-plane (or on the non-negative imaginary axis).
func Sqrt(z complex128) complex128 { return cmplx.Sqrt(z) }

// Cbrt returns the principal cube root of z, exp(Log(z)/3), which is real and
// positive for real positive z.
func Cbrt(z complex128) complex128 {
	if z == 0 {
		return 0
	}
	return cmplx.Exp(cmplx.Log(z) / 3)
}

// Exp returns e**z.
func Exp(z complex128) complex128 { return cmplx.Exp(z) }

// Log returns the principal branch of the natural logarithm of z, with
// imaginary part (the argument) in (-pi, pi].
func Log(z complex128) complex128 { return cmplx.Log(z) }

// LogBranch returns the value of the natural logarithm of z on branch k, that
// is Log(z) + 2*pi*i*k. Branch 0 is the principal branch.
func LogBranch(z complex128, k int) complex128 {
	return cmplx.Log(z) + complex(0, 2*math.Pi*float64(k))
}

// Pow returns z raised to the power w using the principal branch, exp(w*Log(z)).
// Pow(0, 0) is defined to be 1.
func Pow(z, w complex128) complex128 {
	if z == 0 {
		if w == 0 {
			return 1
		}
		return 0
	}
	return cmplx.Pow(z, w)
}

// NthRoots returns all n distinct n-th roots of z, ordered by increasing
// argument starting from the principal root. It returns nil for n <= 0.
func NthRoots(z complex128, n int) []complex128 {
	if n <= 0 {
		return nil
	}
	roots := make([]complex128, n)
	r := math.Pow(cmplx.Abs(z), 1/float64(n))
	theta := cmplx.Phase(z)
	for k := 0; k < n; k++ {
		angle := (theta + 2*math.Pi*float64(k)) / float64(n)
		roots[k] = cmplx.Rect(r, angle)
	}
	return roots
}

// Reciprocal returns 1/z.
func Reciprocal(z complex128) complex128 { return 1 / z }

// Sign returns z/|z|, the unit complex number with the same argument as z, or 0
// when z is 0.
func Sign(z complex128) complex128 {
	a := cmplx.Abs(z)
	if a == 0 {
		return 0
	}
	return z / complex(a, 0)
}

// Sin returns the sine of z.
func Sin(z complex128) complex128 { return cmplx.Sin(z) }

// Cos returns the cosine of z.
func Cos(z complex128) complex128 { return cmplx.Cos(z) }

// Tan returns the tangent of z.
func Tan(z complex128) complex128 { return cmplx.Tan(z) }

// Cot returns the cotangent of z, cos(z)/sin(z).
func Cot(z complex128) complex128 { return cmplx.Cos(z) / cmplx.Sin(z) }

// Sec returns the secant of z, 1/cos(z).
func Sec(z complex128) complex128 { return 1 / cmplx.Cos(z) }

// Csc returns the cosecant of z, 1/sin(z).
func Csc(z complex128) complex128 { return 1 / cmplx.Sin(z) }

// Sinh returns the hyperbolic sine of z.
func Sinh(z complex128) complex128 { return cmplx.Sinh(z) }

// Cosh returns the hyperbolic cosine of z.
func Cosh(z complex128) complex128 { return cmplx.Cosh(z) }

// Tanh returns the hyperbolic tangent of z.
func Tanh(z complex128) complex128 { return cmplx.Tanh(z) }

// Coth returns the hyperbolic cotangent of z, cosh(z)/sinh(z).
func Coth(z complex128) complex128 { return cmplx.Cosh(z) / cmplx.Sinh(z) }

// Asin returns the principal inverse sine of z.
func Asin(z complex128) complex128 { return cmplx.Asin(z) }

// Acos returns the principal inverse cosine of z.
func Acos(z complex128) complex128 { return cmplx.Acos(z) }

// Atan returns the principal inverse tangent of z.
func Atan(z complex128) complex128 { return cmplx.Atan(z) }

// Asinh returns the principal inverse hyperbolic sine of z.
func Asinh(z complex128) complex128 { return cmplx.Asinh(z) }

// Acosh returns the principal inverse hyperbolic cosine of z.
func Acosh(z complex128) complex128 { return cmplx.Acosh(z) }

// Atanh returns the principal inverse hyperbolic tangent of z.
func Atanh(z complex128) complex128 { return cmplx.Atanh(z) }
