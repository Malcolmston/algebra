package complexanalysis

import "math/cmplx"

// LaurentCoefficient returns the k-th Laurent coefficient c_k of f expanded
// about center, computed from the contour integral
// c_k = 1/(2*pi*i) * integral of f(z)/(z-center)^(k+1) around the circle of the
// given radius with n sample points. The index k may be negative; c_{-1} is the
// residue. The circle must lie in an annulus of analyticity of f.
func LaurentCoefficient(f Function, center complex128, k int, radius float64, n int) complex128 {
	g := func(z complex128) complex128 {
		return f(z) / cmplx.Pow(z-center, complex(float64(k+1), 0))
	}
	return IntegrateCircle(g, center, radius, n) / complexanalysisTwoPiI
}

// LaurentCoefficients returns the Laurent coefficients c_k of f about center for
// k running from lo to hi inclusive, as a slice indexed so that result[k-lo] is
// c_k. It returns nil when hi < lo.
func LaurentCoefficients(f Function, center complex128, lo, hi int, radius float64, n int) []complex128 {
	if hi < lo {
		return nil
	}
	out := make([]complex128, hi-lo+1)
	for k := lo; k <= hi; k++ {
		out[k-lo] = LaurentCoefficient(f, center, k, radius, n)
	}
	return out
}

// TaylorCoefficient returns the k-th Taylor coefficient of an analytic function
// f about center, a_k = f^(k)(center)/k!, computed from the same contour
// integral as LaurentCoefficient. It panics if k < 0.
func TaylorCoefficient(f Function, center complex128, k int, radius float64, n int) complex128 {
	if k < 0 {
		panic("complexanalysis: TaylorCoefficient requires k >= 0")
	}
	return LaurentCoefficient(f, center, k, radius, n)
}

// TaylorCoefficients returns the first count Taylor coefficients (orders 0
// through count-1) of an analytic function f about center. It returns nil for
// count <= 0.
func TaylorCoefficients(f Function, center complex128, count int, radius float64, n int) []complex128 {
	if count <= 0 {
		return nil
	}
	out := make([]complex128, count)
	for k := 0; k < count; k++ {
		out[k] = LaurentCoefficient(f, center, k, radius, n)
	}
	return out
}

// PowerSeriesEval evaluates the power series with the given coefficients about
// center at the point z, returning sum_j coeffs[j] * (z-center)^j via Horner's
// method.
func PowerSeriesEval(coeffs []complex128, center, z complex128) complex128 {
	var acc complex128
	w := z - center
	for j := len(coeffs) - 1; j >= 0; j-- {
		acc = acc*w + coeffs[j]
	}
	return acc
}

// AnalyticContinuation returns the value of f at the point to, reconstructed
// from its Taylor expansion about the point from. It computes the Taylor
// coefficients up to the given order on a circle of the given radius about
// from, then sums the truncated series at to. The step |to-from| should be
// smaller than the radius of convergence for the series to be accurate.
func AnalyticContinuation(f Function, from, to complex128, order int, radius float64, n int) complex128 {
	coeffs := TaylorCoefficients(f, from, order+1, radius, n)
	return PowerSeriesEval(coeffs, from, to)
}

// AnalyticContinuationPath continues f numerically along the given path of
// centers and returns the reconstructed value at the final point. Starting from
// path[0], it repeatedly forms the Taylor expansion of f about the current
// center (using a circle of the given radius and n points, to the given order)
// and evaluates it at the next path point, which becomes the new center. This
// mirrors classical analytic continuation along a chain of overlapping disks.
// It returns 0 for an empty path and f(path[0]) for a single point.
func AnalyticContinuationPath(f Function, path []complex128, order int, radius float64, n int) complex128 {
	if len(path) == 0 {
		return 0
	}
	if len(path) == 1 {
		return f(path[0])
	}
	value := f(path[0])
	for i := 0; i+1 < len(path); i++ {
		coeffs := TaylorCoefficients(f, path[i], order+1, radius, n)
		value = PowerSeriesEval(coeffs, path[i], path[i+1])
	}
	return value
}
