package exterioralgebra

// Binomial returns the binomial coefficient C(n,k), the number of k-element
// subsets of an n-element set. It returns 0 when k is negative or exceeds n.
func Binomial(n, k int) int {
	if k < 0 || k > n {
		return 0
	}
	if k > n-k {
		k = n - k
	}
	res := 1
	for i := 0; i < k; i++ {
		res = res * (n - i) / (i + 1)
	}
	return res
}

// GradeDimension returns the dimension of the grade-k subspace Λ^k(Rⁿ), which
// equals C(n,k).
func GradeDimension(n, k int) int { return Binomial(n, k) }

// AlgebraDimension returns the dimension 2ⁿ of the full exterior algebra
// Λ(Rⁿ) as a real vector space.
func AlgebraDimension(n int) int {
	if n < 0 {
		return 0
	}
	return 1 << uint(n)
}

// Basis returns every grade-k basis blade of Λ(Rⁿ) as a unit-coefficient Form,
// listed in canonical (increasing mask) order. The result has GradeDimension(n,k)
// entries.
func Basis(n, k int) []*Form {
	var out []*Form
	full := FullMask(n)
	all := make(map[uint]float64)
	for m := uint(0); m <= full; m++ {
		if Popcount(m) == k {
			all[m] = 1
		}
	}
	for _, m := range sortedMasks(all) {
		f := New(n)
		f.terms[m] = 1
		out = append(out, f)
	}
	return out
}

// AllBasisBlades returns all 2ⁿ basis blades of Λ(Rⁿ) as unit Forms, ordered by
// grade and then by mask.
func AllBasisBlades(n int) []*Form {
	var out []*Form
	for k := 0; k <= n; k++ {
		out = append(out, Basis(n, k)...)
	}
	return out
}

// WedgeSign returns the graded-commutation sign (−1)^{pq} relating α∧β to β∧α
// for a grade-p and a grade-q homogeneous Form.
func WedgeSign(p, q int) int {
	if (p*q)&1 == 1 {
		return -1
	}
	return 1
}

// BivectorPart returns the grade-2 part of f.
func (f *Form) BivectorPart() *Form { return f.GradeProject(2) }

// EvenPart returns the sum of the even-grade blades of f (the even subalgebra
// component).
func (f *Form) EvenPart() *Form {
	g := New(f.n)
	for m, c := range f.terms {
		if Popcount(m)&1 == 0 {
			g.terms[m] = c
		}
	}
	return g
}

// OddPart returns the sum of the odd-grade blades of f.
func (f *Form) OddPart() *Form {
	g := New(f.n)
	for m, c := range f.terms {
		if Popcount(m)&1 == 1 {
			g.terms[m] = c
		}
	}
	return g
}

// ScalarValue is an alias for [Form.ScalarPart]: the grade-0 coefficient of f.
func (f *Form) ScalarValue() float64 { return f.ScalarPart() }

// IsScalar reports whether f has only a grade-0 part (it may be the zero Form).
func (f *Form) IsScalar() bool {
	for m := range f.terms {
		if m != 0 {
			return false
		}
	}
	return true
}

// IsVector reports whether every nonzero blade of f has grade 1.
func (f *Form) IsVector() bool {
	for m := range f.terms {
		if Popcount(m) != 1 {
			return false
		}
	}
	return true
}

// IsBlade reports whether f is a single basis blade times a scalar (at most one
// nonzero term), and hence trivially decomposable.
func (f *Form) IsBlade() bool { return len(f.terms) <= 1 }

// Coeffs returns a fresh copy of the blade→coefficient map underlying f.
func (f *Form) Coeffs() map[uint]float64 {
	out := make(map[uint]float64, len(f.terms))
	for m, c := range f.terms {
		out[m] = c
	}
	return out
}

// Vars returns the n coordinate polynomials x_0,…,x_{n-1} of arity n.
func Vars(n int) []*Poly {
	out := make([]*Poly, n)
	for i := range out {
		out[i] = Var(n, i)
	}
	return out
}

// IsConstant reports whether p has no variable dependence (degree ≤ 0).
func (p *Poly) IsConstant() bool { return p.Degree() <= 0 }

// Coefficient returns the coefficient of the monomial ∏ x_i^{exp[i]} in p. It
// returns 0 for a length mismatch or an absent monomial.
func (p *Poly) Coefficient(exp []int) float64 {
	if len(exp) != p.n {
		return 0
	}
	return p.terms[encodeExp(exp)].coeff
}

// Laplacian returns the scalar Laplacian Σ_i ∂²p/∂x_i² of the polynomial p.
func (p *Poly) Laplacian() *Poly {
	r := NewPoly(p.n)
	for i := 0; i < p.n; i++ {
		r = r.Add(p.Partial(i).Partial(i))
	}
	return r
}

// Hessian returns the n×n matrix of second partial derivatives of p, with entry
// [i][j] equal to ∂²p/∂x_i∂x_j.
func (p *Poly) Hessian() [][]*Poly {
	h := make([][]*Poly, p.n)
	for i := 0; i < p.n; i++ {
		h[i] = make([]*Poly, p.n)
		di := p.Partial(i)
		for j := 0; j < p.n; j++ {
			h[i][j] = di.Partial(j)
		}
	}
	return h
}

// PolyGradient returns the gradient of the scalar field f as its n component
// polynomials, with entry i equal to ∂f/∂x_i. It realises the first map of the
// de Rham complex, d on 0-forms.
func PolyGradient(f *Poly) []*Poly { return f.Gradient() }

// PolyDivergence returns the divergence Σ_i ∂F_i/∂x_i of the vector field whose
// components are field. It returns [ErrDim] unless there are exactly n
// components, each of arity n.
func PolyDivergence(field []*Poly) (*Poly, error) {
	n := len(field)
	for _, c := range field {
		if c.n != n {
			return nil, ErrDim
		}
	}
	if n == 0 {
		return nil, ErrDim
	}
	r := NewPoly(n)
	for i := 0; i < n; i++ {
		r = r.Add(field[i].Partial(i))
	}
	return r, nil
}

// PolyCurl returns the curl of a 3-dimensional vector field field = (P,Q,R),
// namely (∂R/∂x₁−∂Q/∂x₂, ∂P/∂x₂−∂R/∂x₀, ∂Q/∂x₀−∂P/∂x₁). It returns [ErrDim]
// unless there are exactly three components each of arity 3.
func PolyCurl(field []*Poly) ([]*Poly, error) {
	if len(field) != 3 {
		return nil, ErrDim
	}
	for _, c := range field {
		if c.n != 3 {
			return nil, ErrDim
		}
	}
	p, q, r := field[0], field[1], field[2]
	return []*Poly{
		r.Partial(1).Sub(q.Partial(2)),
		p.Partial(2).Sub(r.Partial(0)),
		q.Partial(0).Sub(p.Partial(1)),
	}, nil
}

// PolyLaplacian is a package-level alias for [Poly.Laplacian].
func PolyLaplacian(f *Poly) *Poly { return f.Laplacian() }

// Grades returns the sorted list of distinct grades present in the polynomial
// form w.
func (w *PForm) Grades() []int {
	seen := make(map[int]bool)
	for m := range w.terms {
		seen[Popcount(m)] = true
	}
	out := make([]int, 0, len(seen))
	for g := range seen {
		out = append(out, g)
	}
	for i := 1; i < len(out); i++ {
		for j := i; j > 0 && out[j-1] > out[j]; j-- {
			out[j-1], out[j] = out[j], out[j-1]
		}
	}
	return out
}

// IsHomogeneous reports whether every nonzero blade of w has the same grade.
func (w *PForm) IsHomogeneous() bool { return len(w.Grades()) <= 1 }

// MaxGrade returns the largest grade present in w, or 0 for the zero form.
func (w *PForm) MaxGrade() int {
	max := 0
	for m := range w.terms {
		if g := Popcount(m); g > max {
			max = g
		}
	}
	return max
}

// IsClosed reports whether w is a closed form, that is dw = 0.
func (w *PForm) IsClosed() bool { return w.ExteriorDerivative().IsZero() }

// InteriorConst returns the interior product ι_v w of the polynomial form w with
// the constant vector v of length n, lowering grade by one. It panics if len(v)
// differs from the ambient dimension.
func (w *PForm) InteriorConst(v []float64) *PForm {
	if len(v) != w.n {
		panic(ErrDim)
	}
	res := NewPForm(w.n)
	for i := 0; i < w.n; i++ {
		if v[i] == 0 {
			continue
		}
		for m, p := range w.terms {
			if nm, s, ok := interiorBladeVec(i, m); ok {
				res.addPoly(nm, p.Scale(v[i]*float64(s)))
			}
		}
	}
	return res
}
