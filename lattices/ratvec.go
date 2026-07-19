package lattices

import (
	"math/big"
	"strings"
)

// RatVec is an exact rational vector represented as a slice of *big.Rat
// coordinates. It is used where exactness matters, such as exact Gram
// determinants and dual bases.
type RatVec []*big.Rat

// NewRatVec returns an n-dimensional rational zero vector with freshly
// allocated entries.
func NewRatVec(n int) RatVec {
	v := make(RatVec, n)
	for i := range v {
		v[i] = new(big.Rat)
	}
	return v
}

// RatVecFromInts builds a RatVec from integer coordinates.
func RatVecFromInts(xs ...int64) RatVec {
	v := make(RatVec, len(xs))
	for i, x := range xs {
		v[i] = new(big.Rat).SetInt64(x)
	}
	return v
}

// RatVecFromFloats builds a RatVec from float64 coordinates, converting each
// exactly to a rational. Non-finite inputs yield a zero entry.
func RatVecFromFloats(xs ...float64) RatVec {
	v := make(RatVec, len(xs))
	for i, x := range xs {
		v[i] = ratFromFloat(x)
	}
	return v
}

// Dim returns the number of coordinates in v.
func (v RatVec) Dim() int { return len(v) }

// Clone returns an independent deep copy of v.
func (v RatVec) Clone() RatVec {
	w := make(RatVec, len(v))
	for i := range v {
		w[i] = new(big.Rat).Set(v[i])
	}
	return w
}

// Add returns the exact sum v+w. It panics if the dimensions differ.
func (v RatVec) Add(w RatVec) RatVec {
	v.mustMatch(w)
	r := make(RatVec, len(v))
	for i := range v {
		r[i] = new(big.Rat).Add(v[i], w[i])
	}
	return r
}

// Sub returns the exact difference v-w. It panics if the dimensions differ.
func (v RatVec) Sub(w RatVec) RatVec {
	v.mustMatch(w)
	r := make(RatVec, len(v))
	for i := range v {
		r[i] = new(big.Rat).Sub(v[i], w[i])
	}
	return r
}

// Scale returns the exact scalar multiple s*v.
func (v RatVec) Scale(s *big.Rat) RatVec {
	r := make(RatVec, len(v))
	for i := range v {
		r[i] = new(big.Rat).Mul(s, v[i])
	}
	return r
}

// Neg returns -v.
func (v RatVec) Neg() RatVec {
	r := make(RatVec, len(v))
	for i := range v {
		r[i] = new(big.Rat).Neg(v[i])
	}
	return r
}

// AddScaled returns the exact combination v + s*w. It panics if the dimensions
// differ.
func (v RatVec) AddScaled(w RatVec, s *big.Rat) RatVec {
	v.mustMatch(w)
	r := make(RatVec, len(v))
	for i := range v {
		t := new(big.Rat).Mul(s, w[i])
		r[i] = t.Add(v[i], t)
	}
	return r
}

// Dot returns the exact inner product <v, w>. It panics if the dimensions
// differ.
func (v RatVec) Dot(w RatVec) *big.Rat {
	v.mustMatch(w)
	s := new(big.Rat)
	t := new(big.Rat)
	for i := range v {
		s.Add(s, t.Mul(v[i], w[i]))
	}
	return s
}

// Norm2 returns the exact squared Euclidean norm <v, v>.
func (v RatVec) Norm2() *big.Rat { return v.Dot(v) }

// IsZero reports whether every coordinate of v is exactly zero.
func (v RatVec) IsZero() bool {
	for _, x := range v {
		if x.Sign() != 0 {
			return false
		}
	}
	return true
}

// Equal reports whether v and w have identical dimension and coordinates.
func (v RatVec) Equal(w RatVec) bool {
	if len(v) != len(w) {
		return false
	}
	for i := range v {
		if v[i].Cmp(w[i]) != 0 {
			return false
		}
	}
	return true
}

// Float returns the float64 approximation of v as a Vec.
func (v RatVec) Float() Vec {
	r := make(Vec, len(v))
	for i := range v {
		f, _ := v[i].Float64()
		r[i] = f
	}
	return r
}

// String renders v as a bracketed, space-separated list of rationals.
func (v RatVec) String() string {
	parts := make([]string, len(v))
	for i, x := range v {
		parts[i] = x.RatString()
	}
	return "[" + strings.Join(parts, " ") + "]"
}

func (v RatVec) mustMatch(w RatVec) {
	if len(v) != len(w) {
		panic(ErrDimMismatch)
	}
}

// ratFromFloat converts a float64 exactly to a *big.Rat, returning zero for
// non-finite inputs.
func ratFromFloat(x float64) *big.Rat {
	r := new(big.Rat).SetFloat64(x)
	if r == nil {
		return new(big.Rat)
	}
	return r
}
