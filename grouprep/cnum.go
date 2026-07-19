package grouprep

import (
	"math"
	"math/cmplx"
)

// defaultTol is a reasonable tolerance for floating-point comparisons on the
// small groups this package targets.
const defaultTol = 1e-9

// Cplx builds a complex128 from its real and imaginary parts. It is a thin
// convenience wrapper around the built-in complex.
func Cplx(re, im float64) complex128 {
	return complex(re, im)
}

// Real returns the real part of z.
func Real(z complex128) float64 {
	return real(z)
}

// Imag returns the imaginary part of z.
func Imag(z complex128) float64 {
	return imag(z)
}

// Conj returns the complex conjugate of z.
func Conj(z complex128) complex128 {
	return cmplx.Conj(z)
}

// AbsC returns the modulus |z| of the complex number z.
func AbsC(z complex128) float64 {
	return cmplx.Abs(z)
}

// ApproxEqualC reports whether a and b are within tol of each other in
// modulus, i.e. |a-b| <= tol.
func ApproxEqualC(a, b complex128, tol float64) bool {
	return cmplx.Abs(a-b) <= tol
}

// ApproxZeroC reports whether |z| <= tol.
func ApproxZeroC(z complex128, tol float64) bool {
	return cmplx.Abs(z) <= tol
}

// IsRealC reports whether z is real to within tol, i.e. |Imag(z)| <= tol.
func IsRealC(z complex128, tol float64) bool {
	return math.Abs(imag(z)) <= tol
}

// IsIntegerC reports whether z is a (real) integer to within tol.
func IsIntegerC(z complex128, tol float64) bool {
	if math.Abs(imag(z)) > tol {
		return false
	}
	return math.Abs(real(z)-math.Round(real(z))) <= tol
}

// RoundC rounds both the real and imaginary parts of z to the given number of
// decimal places. It is useful for presenting otherwise noisy floating-point
// results.
func RoundC(z complex128, decimals int) complex128 {
	scale := math.Pow(10, float64(decimals))
	re := math.Round(real(z)*scale) / scale
	im := math.Round(imag(z)*scale) / scale
	return complex(re, im)
}

// RoundToInteger returns the nearest Gaussian integer to z. It is used to snap
// character inner products, which are known to be integers, to exact values.
func RoundToInteger(z complex128) complex128 {
	return complex(math.Round(real(z)), math.Round(imag(z)))
}

// RootOfUnity returns the primitive-n-th root exp(2πi k / n), the k-th power of
// the standard primitive n-th root of unity. It panics if n <= 0.
func RootOfUnity(n, k int) complex128 {
	if n <= 0 {
		panic("grouprep: RootOfUnity requires n > 0")
	}
	theta := 2 * math.Pi * float64(k) / float64(n)
	return cmplx.Exp(complex(0, theta))
}

// PrimitiveRootOfUnity returns exp(2πi / n), the standard primitive n-th root
// of unity. It panics if n <= 0.
func PrimitiveRootOfUnity(n int) complex128 {
	return RootOfUnity(n, 1)
}
