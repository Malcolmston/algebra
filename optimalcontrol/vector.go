package optimalcontrol

import "math"

// VecAdd returns the element-wise sum a+b.
func VecAdd(a, b []float64) []float64 {
	out := make([]float64, len(a))
	for i := range a {
		out[i] = a[i] + b[i]
	}
	return out
}

// VecSub returns the element-wise difference a−b.
func VecSub(a, b []float64) []float64 {
	out := make([]float64, len(a))
	for i := range a {
		out[i] = a[i] - b[i]
	}
	return out
}

// VecScale returns the vector a scaled by s.
func VecScale(a []float64, s float64) []float64 {
	out := make([]float64, len(a))
	for i := range a {
		out[i] = a[i] * s
	}
	return out
}

// VecAxpy returns a + s·b (the classic "axpy" operation).
func VecAxpy(s float64, b, a []float64) []float64 {
	out := make([]float64, len(a))
	for i := range a {
		out[i] = a[i] + s*b[i]
	}
	return out
}

// VecDot returns the inner product aᵀb.
func VecDot(a, b []float64) float64 {
	var s float64
	for i := range a {
		s += a[i] * b[i]
	}
	return s
}

// VecNorm returns the Euclidean norm of a.
func VecNorm(a []float64) float64 {
	return math.Sqrt(VecDot(a, a))
}

// VecMaxAbs returns the largest absolute component of a.
func VecMaxAbs(a []float64) float64 {
	var mx float64
	for _, v := range a {
		if x := math.Abs(v); x > mx {
			mx = x
		}
	}
	return mx
}

// VecCopy returns a copy of a.
func VecCopy(a []float64) []float64 {
	out := make([]float64, len(a))
	copy(out, a)
	return out
}
