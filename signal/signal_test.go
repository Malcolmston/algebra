package signal

import "math"

const tol = 1e-9

// approx reports whether a and b are within an absolute tolerance eps.
func approx(a, b, eps float64) bool {
	if math.IsNaN(a) || math.IsNaN(b) {
		return false
	}
	return math.Abs(a-b) <= eps
}

// approxSlice reports whether two slices match element-wise within eps.
func approxSlice(a, b []float64, eps float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !approx(a[i], b[i], eps) {
			return false
		}
	}
	return true
}

// cmag returns the magnitude of a complex value, for test assertions.
func cmag(c complex128) float64 { return math.Hypot(real(c), imag(c)) }
