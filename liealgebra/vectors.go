package liealgebra

import "math"

// VecDot returns the Euclidean dot product of two equal-length vectors, or
// [ErrDim] on a length mismatch.
func VecDot(a, b []float64) (float64, error) {
	if len(a) != len(b) {
		return 0, ErrDim
	}
	s := 0.0
	for i := range a {
		s += a[i] * b[i]
	}
	return s, nil
}

// VecAdd returns a+b, or [ErrDim] on a length mismatch.
func VecAdd(a, b []float64) ([]float64, error) {
	if len(a) != len(b) {
		return nil, ErrDim
	}
	out := make([]float64, len(a))
	for i := range a {
		out[i] = a[i] + b[i]
	}
	return out, nil
}

// VecSub returns a-b, or [ErrDim] on a length mismatch.
func VecSub(a, b []float64) ([]float64, error) {
	if len(a) != len(b) {
		return nil, ErrDim
	}
	out := make([]float64, len(a))
	for i := range a {
		out[i] = a[i] - b[i]
	}
	return out, nil
}

// VecScale returns the vector scaled by s.
func VecScale(a []float64, s float64) []float64 {
	out := make([]float64, len(a))
	for i := range a {
		out[i] = s * a[i]
	}
	return out
}

// VecNorm returns the Euclidean norm of a vector.
func VecNorm(a []float64) float64 {
	s := 0.0
	for _, v := range a {
		s += v * v
	}
	return math.Sqrt(s)
}

// VecNormSquared returns the squared Euclidean norm of a vector.
func VecNormSquared(a []float64) float64 {
	s := 0.0
	for _, v := range a {
		s += v * v
	}
	return s
}

// VecEqual reports whether two vectors are equal to within tolerance tol.
func VecEqual(a, b []float64, tol float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if math.Abs(a[i]-b[i]) > tol {
			return false
		}
	}
	return true
}
