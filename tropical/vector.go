package tropical

import (
	"math"
	"strings"
)

// Vector is a dense tropical vector: a slice of scalar entries together with
// the semiring under which it is interpreted.
type Vector struct {
	data []float64
	sr   Semiring
}

// NewVector returns a Vector over the given semiring holding a copy of data.
func NewVector(sr Semiring, data []float64) Vector {
	c := make([]float64, len(data))
	copy(c, data)
	return Vector{data: c, sr: sr}
}

// MinPlusVector returns a min-plus Vector holding a copy of data.
func MinPlusVector(data []float64) Vector { return NewVector(MinPlusSemiring(), data) }

// MaxPlusVector returns a max-plus Vector holding a copy of data.
func MaxPlusVector(data []float64) Vector { return NewVector(MaxPlusSemiring(), data) }

// ZeroVector returns a length-n Vector whose entries are all the tropical zero
// of the semiring.
func ZeroVector(sr Semiring, n int) Vector {
	d := make([]float64, n)
	z := sr.Zero()
	for i := range d {
		d[i] = z
	}
	return Vector{data: d, sr: sr}
}

// OneVector returns a length-n Vector whose entries are all the tropical one
// (0) of the semiring.
func OneVector(sr Semiring, n int) Vector {
	return Vector{data: make([]float64, n), sr: sr}
}

// UnitVector returns the length-n tropical unit vector: the tropical one (0) at
// index i and the tropical zero elsewhere. It panics if i is out of range.
func UnitVector(sr Semiring, n, i int) Vector {
	if i < 0 || i >= n {
		panic("tropical: UnitVector index out of range")
	}
	v := ZeroVector(sr, n)
	v.data[i] = 0
	return v
}

// Len returns the number of entries in the vector.
func (v Vector) Len() int { return len(v.data) }

// Semiring returns the semiring under which the vector is interpreted.
func (v Vector) Semiring() Semiring { return v.sr }

// At returns the i-th entry. It panics if i is out of range.
func (v Vector) At(i int) float64 { return v.data[i] }

// Set stores value at index i. It panics if i is out of range.
func (v Vector) Set(i int, value float64) { v.data[i] = value }

// Slice returns a fresh copy of the underlying entries.
func (v Vector) Slice() []float64 {
	c := make([]float64, len(v.data))
	copy(c, v.data)
	return c
}

// Clone returns an independent copy of the vector.
func (v Vector) Clone() Vector { return NewVector(v.sr, v.data) }

// Equal reports whether v and w have the same semiring, length and identical
// entries.
func (v Vector) Equal(w Vector) bool {
	if v.sr.kind != w.sr.kind || len(v.data) != len(w.data) {
		return false
	}
	for i := range v.data {
		if v.data[i] != w.data[i] {
			return false
		}
	}
	return true
}

// EqualTol reports whether v and w have the same semiring and length and every
// pair of entries agrees to within tol (with infinities required to match
// exactly and in sign).
func (v Vector) EqualTol(w Vector, tol float64) bool {
	if v.sr.kind != w.sr.kind || len(v.data) != len(w.data) {
		return false
	}
	for i := range v.data {
		if !closeScalar(v.data[i], w.data[i], tol) {
			return false
		}
	}
	return true
}

// Add returns the elementwise tropical sum of v and w. It panics if the lengths
// or semirings differ.
func (v Vector) Add(w Vector) Vector {
	v.mustMatch(w)
	out := make([]float64, len(v.data))
	for i := range out {
		out[i] = v.sr.Add(v.data[i], w.data[i])
	}
	return Vector{data: out, sr: v.sr}
}

// ScalarMul returns the vector obtained by tropically multiplying every entry
// by c.
func (v Vector) ScalarMul(c float64) Vector {
	out := make([]float64, len(v.data))
	for i := range out {
		out[i] = v.sr.Mul(v.data[i], c)
	}
	return Vector{data: out, sr: v.sr}
}

// ScalarAdd returns the vector obtained by tropically adding c to every entry.
func (v Vector) ScalarAdd(c float64) Vector {
	out := make([]float64, len(v.data))
	for i := range out {
		out[i] = v.sr.Add(v.data[i], c)
	}
	return Vector{data: out, sr: v.sr}
}

// Dot returns the tropical inner product of v and w: the tropical sum over i of
// v[i] (*) w[i]. It panics if the lengths or semirings differ.
func (v Vector) Dot(w Vector) float64 {
	v.mustMatch(w)
	r := v.sr.Zero()
	for i := range v.data {
		r = v.sr.Add(r, v.sr.Mul(v.data[i], w.data[i]))
	}
	return r
}

// Sum returns the tropical sum of all entries (the minimum for min-plus, the
// maximum for max-plus).
func (v Vector) Sum() float64 { return v.sr.Sum(v.data...) }

// Prod returns the tropical product of all entries (their ordinary sum).
func (v Vector) Prod() float64 { return v.sr.Prod(v.data...) }

// Min returns the smallest entry, ignoring sign of infinities in the usual
// numeric order. For an empty vector it returns +Inf.
func (v Vector) Min() float64 {
	m := math.Inf(1)
	for _, x := range v.data {
		if x < m {
			m = x
		}
	}
	return m
}

// Max returns the largest entry in the usual numeric order. For an empty vector
// it returns -Inf.
func (v Vector) Max() float64 {
	m := math.Inf(-1)
	for _, x := range v.data {
		if x > m {
			m = x
		}
	}
	return m
}

// ArgBest returns the index of the entry selected by tropical addition (the
// index of the minimum for min-plus, of the maximum for max-plus) and its
// value. For an empty vector it returns (-1, tropical zero).
func (v Vector) ArgBest() (int, float64) {
	best := v.sr.Zero()
	idx := -1
	for i, x := range v.data {
		if idx == -1 || v.sr.AtLeastAsGood(x, best) && x != best {
			best = x
			idx = i
		}
	}
	return idx, best
}

// Normalize returns the vector with its tropical-product "length" shifted so
// that the best entry becomes the tropical one (0): every entry has the best
// entry subtracted. For an all-zero (tropical) vector it returns a copy.
func (v Vector) Normalize() Vector {
	shift := v.Sum()
	if v.sr.IsZero(shift) {
		return v.Clone()
	}
	out := make([]float64, len(v.data))
	for i := range out {
		out[i] = v.sr.Div(v.data[i], shift)
	}
	return Vector{data: out, sr: v.sr}
}

// IsZero reports whether every entry equals the tropical zero.
func (v Vector) IsZero() bool {
	z := v.sr.Zero()
	for _, x := range v.data {
		if x != z {
			return false
		}
	}
	return true
}

// String renders the vector as a bracketed, space-separated list.
func (v Vector) String() string {
	parts := make([]string, len(v.data))
	for i, x := range v.data {
		parts[i] = v.sr.FormatScalar(x)
	}
	return "[" + strings.Join(parts, " ") + "]"
}

func (v Vector) mustMatch(w Vector) {
	if v.sr.kind != w.sr.kind {
		panic("tropical: vector semiring mismatch")
	}
	if len(v.data) != len(w.data) {
		panic("tropical: vector length mismatch")
	}
}

func closeScalar(a, b, tol float64) bool {
	if math.IsInf(a, 0) || math.IsInf(b, 0) {
		return a == b
	}
	return math.Abs(a-b) <= tol
}
