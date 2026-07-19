package exterioralgebra

// HodgeStar returns the Euclidean Hodge dual ★f, computed with the positively
// oriented orthonormal frame and volume element e_0∧…∧e_{n-1}. On a grade-k
// blade e_I it returns the sign-corrected complementary blade so that
// e_I ∧ ★e_I = vol. It maps grade k to grade n−k and satisfies
// ★★ = (−1)^{k(n−k)} on grade-k Forms.
func HodgeStar(f *Form) *Form {
	full := FullMask(f.n)
	res := New(f.n)
	for m, c := range f.terms {
		j := full &^ m
		res.addTerm(j, float64(reorderSign(m, j))*c)
	}
	return res
}

// Star is the method form of [HodgeStar]: f.Star() returns ★f.
func (f *Form) Star() *Form { return HodgeStar(f) }

// InverseHodgeStar returns the inverse Euclidean Hodge dual, the unique
// operator with InverseHodgeStar(HodgeStar(f)) == f. On a grade-k Form it
// equals (−1)^{k(n−k)}★.
func InverseHodgeStar(f *Form) *Form {
	full := FullMask(f.n)
	n := f.n
	res := New(n)
	for m, c := range f.terms {
		k := Popcount(m)
		j := full &^ m
		sign := reorderSign(m, j)
		if (k*(n-k))&1 == 1 {
			sign = -sign
		}
		res.addTerm(j, float64(sign)*c)
	}
	return res
}

// HodgeSignSquared returns the scalar ε with ★★ = ε·id on grade-k Forms of
// Λ(Rⁿ) under the Euclidean metric, namely (−1)^{k(n−k)}.
func HodgeSignSquared(k, n int) int {
	if (k*(n-k))&1 == 1 {
		return -1
	}
	return 1
}

// HodgeStarMetric returns the Hodge dual of f with respect to a pseudo-Euclidean
// metric whose diagonal signature entries are each +1 or −1. The volume element
// is again e_0∧…∧e_{n-1}. It returns [ErrDim] if signature has the wrong length
// and [ErrGrade] if any entry is not ±1. With an all-+1 signature it agrees
// with [HodgeStar].
func HodgeStarMetric(f *Form, signature []int) (*Form, error) {
	if len(signature) != f.n {
		return nil, ErrDim
	}
	for _, s := range signature {
		if s != 1 && s != -1 {
			return nil, ErrGrade
		}
	}
	full := FullMask(f.n)
	res := New(f.n)
	for m, c := range f.terms {
		j := full &^ m
		prod := 1
		mm := m
		for i := 0; mm != 0; i++ {
			if mm&1 == 1 {
				prod *= signature[i]
			}
			mm >>= 1
		}
		res.addTerm(j, float64(reorderSign(m, j)*prod)*c)
	}
	return res, nil
}

// EuclideanSignature returns the all-+1 metric signature of length n.
func EuclideanSignature(n int) []int {
	s := make([]int, n)
	for i := range s {
		s[i] = 1
	}
	return s
}

// LorentzSignature returns the Lorentzian signature (−1,+1,…,+1) of length n,
// with the time index 0 carrying the −1. For n == 0 it returns an empty slice.
func LorentzSignature(n int) []int {
	s := EuclideanSignature(n)
	if n > 0 {
		s[0] = -1
	}
	return s
}

// WedgeVectors returns the exterior product of the given vectors, each a
// grade-1 Form built from a coefficient slice. All vectors must have the same
// length, which becomes the ambient dimension; it returns [ErrDim] otherwise.
func WedgeVectors(vs ...[]float64) (*Form, error) {
	if len(vs) == 0 {
		return nil, ErrDim
	}
	n := len(vs[0])
	for _, v := range vs {
		if len(v) != n {
			return nil, ErrDim
		}
	}
	acc := FromVector(vs[0])
	for _, v := range vs[1:] {
		acc = acc.Wedge(FromVector(v))
	}
	return acc, nil
}

// Determinant returns the determinant of the n×n matrix whose rows are the n
// supplied vectors of Rⁿ, computed as the top-grade coefficient of their wedge
// product. It returns [ErrDim] unless there are exactly n vectors each of
// length n.
func Determinant(rows ...[]float64) (float64, error) {
	n := len(rows)
	for _, r := range rows {
		if len(r) != n {
			return 0, ErrDim
		}
	}
	w, err := WedgeVectors(rows...)
	if err != nil {
		return 0, err
	}
	return w.terms[FullMask(n)], nil
}

// CrossProduct returns the 3-dimensional cross product a×b, computed as the
// Hodge dual of the wedge a∧b in Λ(R³). It returns [ErrDim] unless both inputs
// have length 3.
func CrossProduct(a, b []float64) ([]float64, error) {
	if len(a) != 3 || len(b) != 3 {
		return nil, ErrDim
	}
	w := FromVector(a).Wedge(FromVector(b))
	return HodgeStar(w).ToVector(), nil
}
