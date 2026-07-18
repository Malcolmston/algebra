package complexanalysis

import (
	"math"
	"math/cmplx"
)

// Re returns the real part of z.
func Re(z complex128) float64 { return real(z) }

// Im returns the imaginary part of z.
func Im(z complex128) float64 { return imag(z) }

// Abs returns the modulus (absolute value) |z|.
func Abs(z complex128) float64 { return cmplx.Abs(z) }

// Arg returns the principal argument of z in the interval (-pi, pi].
func Arg(z complex128) float64 { return cmplx.Phase(z) }

// Conj returns the complex conjugate of z.
func Conj(z complex128) complex128 { return cmplx.Conj(z) }

// Rect constructs a complex number from its real and imaginary parts.
func Rect(re, im float64) complex128 { return complex(re, im) }

// Polar constructs a complex number from a modulus r and an argument theta,
// i.e. r*(cos(theta) + i*sin(theta)).
func Polar(r, theta float64) complex128 { return cmplx.Rect(r, theta) }

// ApproxEqual reports whether |a-b| <= tol.
func ApproxEqual(a, b complex128, tol float64) bool { return cmplx.Abs(a-b) <= tol }

// IsZero reports whether |z| <= tol.
func IsZero(z complex128, tol float64) bool { return cmplx.Abs(z) <= tol }

// complexanalysisFactorialFloat returns n! as a float64 for small non-negative n.
func complexanalysisFactorialFloat(n int) float64 {
	f := 1.0
	for i := 2; i <= n; i++ {
		f *= float64(i)
	}
	return f
}

// complexanalysisTwoPiI is the constant 2*pi*i.
var complexanalysisTwoPiI = complex(0, 2*math.Pi)
