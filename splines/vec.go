package splines

import (
	"errors"
	"math"
)

// Sentinel errors shared across the package.
var (
	// ErrTooFewPoints is returned when a constructor receives fewer sample
	// points than its algorithm requires.
	ErrTooFewPoints = errors.New("splines: too few points")
	// ErrLenMismatch is returned when parallel input slices differ in length.
	ErrLenMismatch = errors.New("splines: input length mismatch")
	// ErrNotIncreasing is returned when a parameter/abscissa slice is not
	// strictly increasing.
	ErrNotIncreasing = errors.New("splines: abscissae must be strictly increasing")
	// ErrDim is returned when control points do not share a common dimension.
	ErrDim = errors.New("splines: points must share a common dimension")
	// ErrDegree is returned for an invalid curve or surface degree.
	ErrDegree = errors.New("splines: invalid degree")
	// ErrKnots is returned for a malformed knot vector.
	ErrKnots = errors.New("splines: invalid knot vector")
	// ErrWeights is returned for invalid (for example non-positive) weights.
	ErrWeights = errors.New("splines: invalid weights")
	// ErrParam is returned when a parameter falls outside the valid domain.
	ErrParam = errors.New("splines: parameter out of range")
	// ErrEmpty is returned when an operation needs a non-empty input.
	ErrEmpty = errors.New("splines: empty input")
)

// Vec is a point or vector in R^n represented by its coordinates. A Vec never
// aliases another Vec once produced by a constructor or arithmetic method: the
// package always returns freshly allocated results, leaving inputs untouched.
type Vec []float64

// NewVec returns a Vec containing a copy of the supplied coordinates.
func NewVec(coords ...float64) Vec {
	v := make(Vec, len(coords))
	copy(v, coords)
	return v
}

// ZeroVec returns the zero vector of dimension n.
func ZeroVec(n int) Vec { return make(Vec, n) }

// Dim reports the number of coordinates in v.
func (v Vec) Dim() int { return len(v) }

// Clone returns an independent copy of v.
func (v Vec) Clone() Vec {
	out := make(Vec, len(v))
	copy(out, v)
	return out
}

// Add returns the component-wise sum v+w. It panics if the dimensions differ.
func (v Vec) Add(w Vec) Vec {
	mustSameDim(v, w)
	out := make(Vec, len(v))
	for i := range v {
		out[i] = v[i] + w[i]
	}
	return out
}

// Sub returns the component-wise difference v-w. It panics if the dimensions
// differ.
func (v Vec) Sub(w Vec) Vec {
	mustSameDim(v, w)
	out := make(Vec, len(v))
	for i := range v {
		out[i] = v[i] - w[i]
	}
	return out
}

// Scale returns v scaled component-wise by s.
func (v Vec) Scale(s float64) Vec {
	out := make(Vec, len(v))
	for i := range v {
		out[i] = v[i] * s
	}
	return out
}

// Neg returns -v.
func (v Vec) Neg() Vec { return v.Scale(-1) }

// AddScaled returns v + s*w, a fused multiply-add. It panics on a dimension
// mismatch.
func (v Vec) AddScaled(s float64, w Vec) Vec {
	mustSameDim(v, w)
	out := make(Vec, len(v))
	for i := range v {
		out[i] = v[i] + s*w[i]
	}
	return out
}

// Dot returns the Euclidean inner product of v and w. It panics on a dimension
// mismatch.
func (v Vec) Dot(w Vec) float64 {
	mustSameDim(v, w)
	var s float64
	for i := range v {
		s += v[i] * w[i]
	}
	return s
}

// Norm2 returns the squared Euclidean norm of v.
func (v Vec) Norm2() float64 { return v.Dot(v) }

// Norm returns the Euclidean norm (length) of v.
func (v Vec) Norm() float64 { return math.Sqrt(v.Norm2()) }

// Dist returns the Euclidean distance between v and w.
func (v Vec) Dist(w Vec) float64 { return v.Sub(w).Norm() }

// Dist2 returns the squared Euclidean distance between v and w.
func (v Vec) Dist2(w Vec) float64 {
	mustSameDim(v, w)
	var s float64
	for i := range v {
		d := v[i] - w[i]
		s += d * d
	}
	return s
}

// Unit returns v scaled to unit length. If v is the zero vector it is returned
// unchanged (as a copy).
func (v Vec) Unit() Vec {
	n := v.Norm()
	if n == 0 {
		return v.Clone()
	}
	return v.Scale(1 / n)
}

// Lerp returns the linear interpolation (1-t)*v + t*w.
func (v Vec) Lerp(w Vec, t float64) Vec {
	mustSameDim(v, w)
	out := make(Vec, len(v))
	for i := range v {
		out[i] = v[i] + t*(w[i]-v[i])
	}
	return out
}

// Equal reports whether v and w are equal within absolute tolerance tol.
func (v Vec) Equal(w Vec, tol float64) bool {
	if len(v) != len(w) {
		return false
	}
	for i := range v {
		if math.Abs(v[i]-w[i]) > tol {
			return false
		}
	}
	return true
}

// VecLerp returns the linear interpolation (1-t)*a + t*b of two vectors.
func VecLerp(a, b Vec, t float64) Vec { return a.Lerp(b, t) }

// CloneVecs returns an independent deep copy of a slice of vectors.
func CloneVecs(pts []Vec) []Vec {
	out := make([]Vec, len(pts))
	for i := range pts {
		out[i] = pts[i].Clone()
	}
	return out
}

// VecsEqual reports whether two vector slices are element-wise equal within
// absolute tolerance tol.
func VecsEqual(a, b []Vec, tol float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !a[i].Equal(b[i], tol) {
			return false
		}
	}
	return true
}

// VecCombine returns the linear combination sum_i coeff[i]*pts[i]. All points
// must share a dimension and coeff must match pts in length.
func VecCombine(coeff []float64, pts []Vec) (Vec, error) {
	if len(coeff) != len(pts) {
		return nil, ErrLenMismatch
	}
	if len(pts) == 0 {
		return nil, ErrEmpty
	}
	d := pts[0].Dim()
	out := make(Vec, d)
	for i, p := range pts {
		if p.Dim() != d {
			return nil, ErrDim
		}
		for j := 0; j < d; j++ {
			out[j] += coeff[i] * p[j]
		}
	}
	return out, nil
}

// commonDim returns the shared dimension of pts, or an error if pts is empty or
// the dimensions are not uniform.
func commonDim(pts []Vec) (int, error) {
	if len(pts) == 0 {
		return 0, ErrEmpty
	}
	d := pts[0].Dim()
	if d == 0 {
		return 0, ErrDim
	}
	for _, p := range pts {
		if p.Dim() != d {
			return 0, ErrDim
		}
	}
	return d, nil
}

func mustSameDim(v, w Vec) {
	if len(v) != len(w) {
		panic("splines: dimension mismatch")
	}
}

// strictlyIncreasing reports whether x is strictly increasing.
func strictlyIncreasing(x []float64) bool {
	for i := 1; i < len(x); i++ {
		if x[i] <= x[i-1] {
			return false
		}
	}
	return true
}

// nondecreasing reports whether x is non-decreasing.
func nondecreasing(x []float64) bool {
	for i := 1; i < len(x); i++ {
		if x[i] < x[i-1] {
			return false
		}
	}
	return true
}

// searchInterval returns the index i such that x[i] <= q < x[i+1], clamping q
// to the ends of the [x[0], x[len-1]] domain. x must have length >= 2 and be
// strictly increasing.
func searchInterval(x []float64, q float64) int {
	n := len(x)
	if q <= x[0] {
		return 0
	}
	if q >= x[n-1] {
		return n - 2
	}
	lo, hi := 0, n-1
	for hi-lo > 1 {
		mid := (lo + hi) / 2
		if x[mid] <= q {
			lo = mid
		} else {
			hi = mid
		}
	}
	return lo
}
