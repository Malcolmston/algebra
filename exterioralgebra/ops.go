package exterioralgebra

// Add returns the sum f+g. It panics if f and g have different ambient
// dimensions.
func (f *Form) Add(g *Form) *Form {
	f.requireSameDim(g)
	h := f.Clone()
	for m, c := range g.terms {
		h.addTerm(m, c)
	}
	return h
}

// Sub returns the difference f-g. It panics if f and g have different ambient
// dimensions.
func (f *Form) Sub(g *Form) *Form {
	f.requireSameDim(g)
	h := f.Clone()
	for m, c := range g.terms {
		h.addTerm(m, -c)
	}
	return h
}

// Neg returns the additive inverse -f.
func (f *Form) Neg() *Form {
	h := New(f.n)
	for m, c := range f.terms {
		h.terms[m] = -c
	}
	return h
}

// Scale returns c·f, the scalar multiple of f by c.
func (f *Form) Scale(c float64) *Form {
	h := New(f.n)
	if c == 0 {
		return h
	}
	for m, v := range f.terms {
		h.terms[m] = c * v
	}
	return h
}

// AddAll returns the sum of any number of Forms, all of which must share an
// ambient dimension. With no arguments it panics, since the dimension is
// unknown; use [New] for an explicit zero.
func AddAll(forms ...*Form) *Form {
	if len(forms) == 0 {
		panic(ErrDim)
	}
	acc := New(forms[0].n)
	for _, f := range forms {
		acc = acc.Add(f)
	}
	return acc
}

// LinearCombination returns Σ coeffs[i]·forms[i]. The slices must have equal
// length and all Forms must share an ambient dimension; it panics otherwise.
func LinearCombination(coeffs []float64, forms []*Form) *Form {
	if len(coeffs) != len(forms) || len(forms) == 0 {
		panic(ErrDim)
	}
	acc := New(forms[0].n)
	for i, f := range forms {
		acc = acc.Add(f.Scale(coeffs[i]))
	}
	return acc
}

// Wedge returns the exterior product f∧g. The product is bilinear and, on
// homogeneous grades p and q, is graded-anticommutative:
// f∧g = (−1)^{pq} g∧f. It panics if the ambient dimensions differ.
func (f *Form) Wedge(g *Form) *Form {
	f.requireSameDim(g)
	h := New(f.n)
	for ma, ca := range f.terms {
		for mb, cb := range g.terms {
			if ma&mb != 0 {
				continue // repeated index kills the blade
			}
			s := reorderSign(ma, mb)
			h.addTerm(ma|mb, float64(s)*ca*cb)
		}
	}
	return h
}

// Wedge returns the exterior product a∧b as a package-level convenience.
func Wedge(a, b *Form) *Form { return a.Wedge(b) }

// WedgeAll returns the exterior product forms[0]∧forms[1]∧… evaluated left to
// right. With no arguments it panics; with one it returns a clone.
func WedgeAll(forms ...*Form) *Form {
	if len(forms) == 0 {
		panic(ErrDim)
	}
	acc := forms[0].Clone()
	for _, f := range forms[1:] {
		acc = acc.Wedge(f)
	}
	return acc
}

// WedgePow returns the k-fold exterior power f∧…∧f (k factors). WedgePow(0) is
// the scalar 1. Because a 1-form wedged with itself vanishes, high powers are
// often zero.
func (f *Form) WedgePow(k int) *Form {
	acc := One(f.n)
	for i := 0; i < k; i++ {
		acc = acc.Wedge(f)
	}
	return acc
}

// Reverse returns the reversion f̃ of f, which reverses the order of the vector
// factors in every blade. On a grade-k blade it multiplies by
// (−1)^{k(k−1)/2}. Reversion is an anti-automorphism of the exterior algebra.
func (f *Form) Reverse() *Form {
	h := New(f.n)
	for m, c := range f.terms {
		k := Popcount(m)
		if (k*(k-1)/2)&1 == 1 {
			c = -c
		}
		h.terms[m] = c
	}
	return h
}

// GradeInvolution returns the grade involution α(f), which negates every
// odd-grade blade: on a grade-k blade it multiplies by (−1)^k.
func (f *Form) GradeInvolution() *Form {
	h := New(f.n)
	for m, c := range f.terms {
		if Popcount(m)&1 == 1 {
			c = -c
		}
		h.terms[m] = c
	}
	return h
}

// requireSameDim panics with ErrDim when g lives in a different ambient
// dimension than f.
func (f *Form) requireSameDim(g *Form) {
	if f.n != g.n {
		panic(ErrDim)
	}
}

// addTerm accumulates coefficient c onto blade mask m, pruning the entry if the
// running total lands exactly on zero.
func (f *Form) addTerm(m uint, c float64) {
	v := f.terms[m] + c
	if v == 0 {
		delete(f.terms, m)
		return
	}
	f.terms[m] = v
}
