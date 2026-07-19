package fuzzy

import "math"

// ComplementMF returns the standard complement 1 - mf(x) of a membership
// function.
func ComplementMF(mf MF) MF {
	return func(x float64) float64 { return clamp01(1 - mf(x)) }
}

// ComplementMFWith returns the complement of mf using an arbitrary fuzzy
// complement operator c (for example ComplementSugeno or ComplementYager).
func ComplementMFWith(mf MF, c func(float64) float64) MF {
	return func(x float64) float64 { return clamp01(c(mf(x))) }
}

// UnionMF returns the union (max) of two membership functions.
func UnionMF(a, b MF) MF {
	return func(x float64) float64 { return clamp01(math.Max(a(x), b(x))) }
}

// IntersectionMF returns the intersection (min) of two membership functions.
func IntersectionMF(a, b MF) MF {
	return func(x float64) float64 { return clamp01(math.Min(a(x), b(x))) }
}

// UnionMFWith returns the union of two membership functions using the t-conorm
// s as the disjunction operator.
func UnionMFWith(a, b MF, s TConorm) MF {
	return func(x float64) float64 { return clamp01(s(a(x), b(x))) }
}

// IntersectionMFWith returns the intersection of two membership functions using
// the t-norm t as the conjunction operator.
func IntersectionMFWith(a, b MF, t TNorm) MF {
	return func(x float64) float64 { return clamp01(t(a(x), b(x))) }
}

// ScaleMF returns the membership function mf multiplied pointwise by the factor
// k and clamped to [0, 1].
func ScaleMF(mf MF, k float64) MF {
	return func(x float64) float64 { return clamp01(mf(x) * k) }
}

// ClipMF returns the membership function mf capped at the level alpha, that is
// min(mf(x), alpha). This is the alpha clipping used by Mamdani implication.
func ClipMF(mf MF, alpha float64) MF {
	return func(x float64) float64 { return clamp01(math.Min(mf(x), alpha)) }
}

// PowMF returns the membership function mf raised pointwise to the power p,
// mf(x)^p. It underlies the concentration and dilation hedges.
func PowMF(mf MF, p float64) MF {
	return func(x float64) float64 { return clamp01(math.Pow(mf(x), p)) }
}

// ConcentrateMF returns the concentration hedge mf(x)^2, sharpening a fuzzy set
// (the linguistic hedge "very").
func ConcentrateMF(mf MF) MF { return PowMF(mf, 2) }

// DilateMF returns the dilation hedge mf(x)^0.5, spreading a fuzzy set (the
// linguistic hedge "somewhat" or "more or less").
func DilateMF(mf MF) MF { return PowMF(mf, 0.5) }

// IntensifyMF returns the contrast intensification hedge that maps grades below
// 0.5 to 2*mf^2 and grades at or above 0.5 to 1 - 2*(1-mf)^2, pushing grades
// toward 0 or 1.
func IntensifyMF(mf MF) MF {
	return func(x float64) float64 {
		m := mf(x)
		if m <= 0.5 {
			return clamp01(2 * m * m)
		}
		d := 1 - m
		return clamp01(1 - 2*d*d)
	}
}

// VeryMF returns the "very" hedge, mf(x)^2.
func VeryMF(mf MF) MF { return PowMF(mf, 2) }

// ExtremelyMF returns the "extremely" hedge, mf(x)^3.
func ExtremelyMF(mf MF) MF { return PowMF(mf, 3) }

// SomewhatMF returns the "somewhat" hedge, mf(x)^0.5.
func SomewhatMF(mf MF) MF { return PowMF(mf, 0.5) }

// SlightlyMF returns the "slightly" hedge, mf(x)^(1/3).
func SlightlyMF(mf MF) MF { return PowMF(mf, 1.0/3.0) }

// MoreOrLessMF returns the "more or less" hedge, mf(x)^0.5, a synonym of
// SomewhatMF.
func MoreOrLessMF(mf MF) MF { return PowMF(mf, 0.5) }

// NotMF returns the "not" hedge, the standard complement 1 - mf(x).
func NotMF(mf MF) MF { return ComplementMF(mf) }

// SampleMF evaluates mf at each point of xs and returns the grades in a new
// slice.
func SampleMF(mf MF, xs []float64) []float64 {
	out := make([]float64, len(xs))
	for i, x := range xs {
		out[i] = clamp01(mf(x))
	}
	return out
}

// Linspace returns n evenly spaced points spanning the closed interval [a, b].
// For n <= 1 it returns a single point a.
func Linspace(a, b float64, n int) []float64 {
	if n <= 1 {
		return []float64{a}
	}
	out := make([]float64, n)
	step := (b - a) / float64(n-1)
	for i := 0; i < n; i++ {
		out[i] = a + step*float64(i)
	}
	out[n-1] = b
	return out
}
