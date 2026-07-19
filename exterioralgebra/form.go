package exterioralgebra

import (
	"fmt"
	"math"
	"strings"
)

// Form is an element of the exterior algebra Λ(Rⁿ): a real linear combination
// of basis blades. Each basis blade e_{i₁}∧…∧e_{i_k} (with i₁<…<i_k) is keyed
// by the bitmask whose set bits are its indices, and only nonzero coefficients
// are stored. A Form therefore represents an arbitrary, possibly mixed-grade,
// multivector or constant differential form.
//
// The zero value is not usable; construct Forms with [New] and friends.
type Form struct {
	n     int
	terms map[uint]float64
}

// Term is one basis blade of a [Form] together with its coefficient. Indices is
// the sorted index list of the blade and Mask its bitmask encoding.
type Term struct {
	Mask    uint
	Indices []int
	Coeff   float64
}

// New returns the zero Form of the exterior algebra Λ(Rⁿ). It panics if n < 0.
func New(n int) *Form {
	if n < 0 {
		panic(ErrDim)
	}
	return &Form{n: n, terms: make(map[uint]float64)}
}

// Zero is an alias for [New]: it returns the additive identity of Λ(Rⁿ).
func Zero(n int) *Form { return New(n) }

// Scalar returns the grade-0 Form with the given scalar coefficient in Λ(Rⁿ).
func Scalar(n int, c float64) *Form {
	f := New(n)
	if c != 0 {
		f.terms[0] = c
	}
	return f
}

// One returns the multiplicative identity of Λ(Rⁿ), the scalar 1.
func One(n int) *Form { return Scalar(n, 1) }

// Basis1 returns the grade-1 basis Form e_i in Λ(Rⁿ). It panics if i is out of
// range.
func Basis1(n, i int) *Form {
	if i < 0 || i >= n {
		panic(ErrIndex)
	}
	f := New(n)
	f.terms[uint(1)<<uint(i)] = 1
	return f
}

// BasisBlade returns the basis blade e_{idx[0]}∧…∧e_{idx[k-1]} in Λ(Rⁿ) with
// unit coefficient, reduced to its canonical sorted representative with the
// appropriate sign. If two indices coincide the blade is zero; the returned
// error is non-nil only when an index is out of range.
func BasisBlade(n int, idx ...int) (*Form, error) {
	return BasisBladeCoeff(n, 1, idx...)
}

// BasisBladeCoeff returns c·e_{idx[0]}∧…∧e_{idx[k-1]} in Λ(Rⁿ). It behaves like
// [BasisBlade] but scales the blade by c.
func BasisBladeCoeff(n int, c float64, idx ...int) (*Form, error) {
	for _, i := range idx {
		if i < 0 || i >= n {
			return nil, ErrIndex
		}
	}
	f := New(n)
	mask, sign, ok := IndicesToMask(n, idx...)
	if ok && c != 0 {
		f.terms[mask] = c * float64(sign)
	}
	return f, nil
}

// FromVector returns the grade-1 Form whose e_i coefficient is v[i]. The ambient
// dimension of the result equals len(v).
func FromVector(v []float64) *Form {
	f := New(len(v))
	for i, c := range v {
		if c != 0 {
			f.terms[uint(1)<<uint(i)] = c
		}
	}
	return f
}

// FromMasks builds a Form in Λ(Rⁿ) directly from a map of blade bitmasks to
// coefficients. Masks with bits at or beyond position n produce [ErrIndex].
func FromMasks(n int, terms map[uint]float64) (*Form, error) {
	full := FullMask(n)
	f := New(n)
	for m, c := range terms {
		if m & ^full != 0 {
			return nil, ErrIndex
		}
		if c != 0 {
			f.terms[m] = c
		}
	}
	return f, nil
}

// VolumeForm returns the canonical top-grade unit volume element
// e_0∧…∧e_{n-1} of Λ(Rⁿ). For n == 0 it is the scalar 1.
func VolumeForm(n int) *Form {
	f := New(n)
	f.terms[FullMask(n)] = 1
	return f
}

// Clone returns an independent deep copy of f.
func (f *Form) Clone() *Form {
	g := New(f.n)
	for m, c := range f.terms {
		g.terms[m] = c
	}
	return g
}

// Dim returns the ambient dimension n of the exterior algebra Λ(Rⁿ) that f
// belongs to.
func (f *Form) Dim() int { return f.n }

// NumTerms returns the number of nonzero basis blades in f.
func (f *Form) NumTerms() int { return len(f.terms) }

// IsZero reports whether f is the zero Form (has no nonzero terms).
func (f *Form) IsZero() bool { return len(f.terms) == 0 }

// CoeffMask returns the coefficient of the sorted basis blade encoded by mask.
// Bits of mask outside the ambient dimension are ignored.
func (f *Form) CoeffMask(mask uint) float64 { return f.terms[mask] }

// Coeff returns the coefficient of the blade e_{idx[0]}∧…∧e_{idx[k-1]},
// accounting for the sign introduced by sorting the indices. It returns 0 for a
// blade that is zero (a repeated index) or that names an out-of-range index.
func (f *Form) Coeff(idx ...int) float64 {
	mask, sign, ok := IndicesToMask(f.n, idx...)
	if !ok {
		return 0
	}
	return float64(sign) * f.terms[mask]
}

// SetMask sets the coefficient of the sorted blade encoded by mask to c,
// removing the term when c == 0. Bits of mask outside the ambient dimension are
// masked off.
func (f *Form) SetMask(mask uint, c float64) {
	mask &= FullMask(f.n)
	if c == 0 {
		delete(f.terms, mask)
		return
	}
	f.terms[mask] = c
}

// Set sets the coefficient of the blade e_{idx[0]}∧…∧e_{idx[k-1]} to c, taking
// the sorting sign into account. It returns [ErrIndex] if an index is out of
// range or repeated.
func (f *Form) Set(c float64, idx ...int) error {
	mask, sign, ok := IndicesToMask(f.n, idx...)
	if !ok {
		return ErrIndex
	}
	f.SetMask(mask, c*float64(sign))
	return nil
}

// Terms returns the nonzero blades of f as a slice of [Term], ordered by grade
// and then by mask. The returned slice is freshly allocated.
func (f *Form) Terms() []Term {
	masks := sortedMasks(f.terms)
	out := make([]Term, len(masks))
	for i, m := range masks {
		out[i] = Term{Mask: m, Indices: MaskToIndices(m), Coeff: f.terms[m]}
	}
	return out
}

// Masks returns the bitmasks of the nonzero blades of f, ordered by grade and
// then by mask.
func (f *Form) Masks() []uint { return sortedMasks(f.terms) }

// GradeProject returns the homogeneous grade-k part of f, a new Form containing
// only the blades of exactly k indices. For k outside [0,n] it returns the zero
// Form.
func (f *Form) GradeProject(k int) *Form {
	g := New(f.n)
	for m, c := range f.terms {
		if Popcount(m) == k {
			g.terms[m] = c
		}
	}
	return g
}

// ScalarPart returns the grade-0 coefficient of f as a plain float64.
func (f *Form) ScalarPart() float64 { return f.terms[0] }

// VectorPart returns the grade-1 part of f as a new Form.
func (f *Form) VectorPart() *Form { return f.GradeProject(1) }

// TopPart returns the top-grade (grade-n) part of f as a new Form.
func (f *Form) TopPart() *Form { return f.GradeProject(f.n) }

// Grades returns the sorted list of distinct grades present in f.
func (f *Form) Grades() []int {
	seen := make(map[int]bool)
	for m := range f.terms {
		seen[Popcount(m)] = true
	}
	out := make([]int, 0, len(seen))
	for g := range seen {
		out = append(out, g)
	}
	// simple insertion sort to avoid importing sort here
	for i := 1; i < len(out); i++ {
		for j := i; j > 0 && out[j-1] > out[j]; j-- {
			out[j-1], out[j] = out[j], out[j-1]
		}
	}
	return out
}

// IsHomogeneous reports whether every nonzero term of f has the same grade. The
// zero Form is considered homogeneous.
func (f *Form) IsHomogeneous() bool { return len(f.Grades()) <= 1 }

// Grade returns the common grade of a homogeneous Form. ok is false when f is
// mixed-grade; the zero Form reports grade 0 with ok true.
func (f *Form) Grade() (k int, ok bool) {
	gs := f.Grades()
	switch len(gs) {
	case 0:
		return 0, true
	case 1:
		return gs[0], true
	default:
		return 0, false
	}
}

// MaxGrade returns the largest grade present in f, or 0 for the zero Form.
func (f *Form) MaxGrade() int {
	max := 0
	for m := range f.terms {
		if g := Popcount(m); g > max {
			max = g
		}
	}
	return max
}

// MinGrade returns the smallest grade present in f, or 0 for the zero Form.
func (f *Form) MinGrade() int {
	first := true
	min := 0
	for m := range f.terms {
		g := Popcount(m)
		if first || g < min {
			min, first = g, false
		}
	}
	return min
}

// Equal reports whether f and g are exactly equal: same ambient dimension and
// identical coefficients on every blade.
func (f *Form) Equal(g *Form) bool {
	if f.n != g.n || len(f.terms) != len(g.terms) {
		return false
	}
	for m, c := range f.terms {
		if g.terms[m] != c {
			return false
		}
	}
	return true
}

// EqualTol reports whether f and g have the same ambient dimension and every
// blade coefficient agrees to within an absolute tolerance tol.
func (f *Form) EqualTol(g *Form, tol float64) bool {
	if f.n != g.n {
		return false
	}
	seen := make(map[uint]bool)
	for m, c := range f.terms {
		if math.Abs(c-g.terms[m]) > tol {
			return false
		}
		seen[m] = true
	}
	for m, c := range g.terms {
		if !seen[m] && math.Abs(c) > tol {
			return false
		}
	}
	return true
}

// ToVector returns the grade-1 part of f as a coefficient slice of length n,
// so that the i-th entry is the coefficient of e_i.
func (f *Form) ToVector() []float64 {
	v := make([]float64, f.n)
	for i := 0; i < f.n; i++ {
		v[i] = f.terms[uint(1)<<uint(i)]
	}
	return v
}

// String renders f in the usual e_{i…} basis notation, for example
// "2 + 3 e0∧e1". The zero Form renders as "0".
func (f *Form) String() string {
	if f.IsZero() {
		return "0"
	}
	var b strings.Builder
	for i, m := range sortedMasks(f.terms) {
		c := f.terms[m]
		if i > 0 {
			if c < 0 {
				b.WriteString(" - ")
				c = -c
			} else {
				b.WriteString(" + ")
			}
		} else if c < 0 {
			b.WriteString("-")
			c = -c
		}
		if m == 0 {
			fmt.Fprintf(&b, "%g", c)
			continue
		}
		fmt.Fprintf(&b, "%g ", c)
		for j, ix := range MaskToIndices(m) {
			if j > 0 {
				b.WriteString("∧")
			}
			fmt.Fprintf(&b, "e%d", ix)
		}
	}
	return b.String()
}
