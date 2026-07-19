package exterioralgebra

import "strings"

// PForm is a differential form on Rⁿ whose coefficients are polynomials: a
// finite sum Σ_I P_I(x) dx^I where each P_I is a [Poly] in the n coordinates
// and dx^I ranges over basis blades. Because polynomial partial derivatives are
// exact, PForm supports an exterior derivative that satisfies d² = 0 exactly, a
// pullback that commutes with d, and an exact Hodge–de Rham Laplacian.
//
// The zero value is not usable; construct PForms with [NewPForm], [PConst],
// [Dx] and [PBasisBlade].
type PForm struct {
	n     int
	terms map[uint]*Poly
}

// NewPForm returns the zero differential form on Rⁿ. It panics if n < 0.
func NewPForm(n int) *PForm {
	if n < 0 {
		panic(ErrDim)
	}
	return &PForm{n: n, terms: make(map[uint]*Poly)}
}

// PConst returns the grade-0 differential form whose single coefficient is the
// polynomial p. The ambient dimension equals the arity of p.
func PConst(p *Poly) *PForm {
	w := NewPForm(p.n)
	if !p.IsZero() {
		w.terms[0] = p.Clone()
	}
	return w
}

// Dx returns the constant basis 1-form dx^i on Rⁿ. It panics if i is out of
// range.
func Dx(n, i int) *PForm {
	if i < 0 || i >= n {
		panic(ErrIndex)
	}
	w := NewPForm(n)
	w.terms[uint(1)<<uint(i)] = ConstPoly(n, 1)
	return w
}

// PBasisBlade returns the differential form coeff·dx^{idx[0]}∧… on Rⁿ, reducing
// the blade to canonical sorted order with the appropriate sign. The polynomial
// coeff must have arity n. It returns [ErrIndex] for an out-of-range index and
// [ErrDim] on an arity mismatch.
func PBasisBlade(n int, coeff *Poly, idx ...int) (*PForm, error) {
	if coeff.n != n {
		return nil, ErrDim
	}
	for _, i := range idx {
		if i < 0 || i >= n {
			return nil, ErrIndex
		}
	}
	w := NewPForm(n)
	mask, sign, ok := IndicesToMask(n, idx...)
	if ok && !coeff.IsZero() {
		w.terms[mask] = coeff.Scale(float64(sign))
	}
	return w, nil
}

// Dim returns the ambient dimension n of the form.
func (w *PForm) Dim() int { return w.n }

// NumTerms returns the number of nonzero blades in w.
func (w *PForm) NumTerms() int { return len(w.terms) }

// IsZero reports whether w has no nonzero coefficient polynomials.
func (w *PForm) IsZero() bool { return len(w.terms) == 0 }

// Clone returns an independent deep copy of w.
func (w *PForm) Clone() *PForm {
	r := NewPForm(w.n)
	for m, p := range w.terms {
		r.terms[m] = p.Clone()
	}
	return r
}

// CoeffMask returns the polynomial coefficient of the sorted blade encoded by
// mask, or the zero polynomial if that blade is absent.
func (w *PForm) CoeffMask(mask uint) *Poly {
	if p, ok := w.terms[mask]; ok {
		return p.Clone()
	}
	return NewPoly(w.n)
}

// Masks returns the bitmasks of the nonzero blades of w, ordered by grade and
// then by mask.
func (w *PForm) Masks() []uint {
	tmp := make(map[uint]float64, len(w.terms))
	for m := range w.terms {
		tmp[m] = 1
	}
	return sortedMasks(tmp)
}

// addPoly accumulates polynomial p onto blade mask m, pruning the entry when
// its coefficient polynomial becomes zero.
func (w *PForm) addPoly(m uint, p *Poly) {
	if p.IsZero() {
		return
	}
	if cur, ok := w.terms[m]; ok {
		s := cur.Add(p)
		if s.IsZero() {
			delete(w.terms, m)
		} else {
			w.terms[m] = s
		}
		return
	}
	w.terms[m] = p.Clone()
}

// Add returns the sum w+u. It panics if the ambient dimensions differ.
func (w *PForm) Add(u *PForm) *PForm {
	if w.n != u.n {
		panic(ErrDim)
	}
	r := w.Clone()
	for m, p := range u.terms {
		r.addPoly(m, p)
	}
	return r
}

// Sub returns the difference w−u. It panics if the ambient dimensions differ.
func (w *PForm) Sub(u *PForm) *PForm {
	if w.n != u.n {
		panic(ErrDim)
	}
	r := w.Clone()
	for m, p := range u.terms {
		r.addPoly(m, p.Neg())
	}
	return r
}

// Neg returns the additive inverse −w.
func (w *PForm) Neg() *PForm {
	r := NewPForm(w.n)
	for m, p := range w.terms {
		r.terms[m] = p.Neg()
	}
	return r
}

// ScalePoly returns the product of w with the scalar polynomial p, multiplying
// every coefficient by p. It panics if the arity of p differs from n.
func (w *PForm) ScalePoly(p *Poly) *PForm {
	if p.n != w.n {
		panic(ErrDim)
	}
	r := NewPForm(w.n)
	for m, q := range w.terms {
		r.addPoly(m, q.Mul(p))
	}
	return r
}

// Scale returns c·w, the real scalar multiple of w.
func (w *PForm) Scale(c float64) *PForm {
	r := NewPForm(w.n)
	for m, p := range w.terms {
		r.addPoly(m, p.Scale(c))
	}
	return r
}

// Wedge returns the exterior product w∧u of two polynomial forms. It panics if
// the ambient dimensions differ.
func (w *PForm) Wedge(u *PForm) *PForm {
	if w.n != u.n {
		panic(ErrDim)
	}
	r := NewPForm(w.n)
	for ma, pa := range w.terms {
		for mb, pb := range u.terms {
			if ma&mb != 0 {
				continue
			}
			s := reorderSign(ma, mb)
			r.addPoly(ma|mb, pa.Mul(pb).Scale(float64(s)))
		}
	}
	return r
}

// GradeProject returns the homogeneous grade-k part of w.
func (w *PForm) GradeProject(k int) *PForm {
	r := NewPForm(w.n)
	for m, p := range w.terms {
		if Popcount(m) == k {
			r.terms[m] = p.Clone()
		}
	}
	return r
}

// ExteriorDerivative returns the exterior derivative dw. For w = Σ_I P_I dx^I it
// is Σ_I Σ_j (∂P_I/∂x_j) dx^j∧dx^I, computed from exact polynomial partials.
// Consequently d(dw) is exactly the zero form for every w.
func (w *PForm) ExteriorDerivative() *PForm {
	r := NewPForm(w.n)
	for m, p := range w.terms {
		for j := 0; j < w.n; j++ {
			bit := uint(1) << uint(j)
			if m&bit != 0 {
				continue
			}
			dp := p.Partial(j)
			if dp.IsZero() {
				continue
			}
			s := reorderSign(bit, m)
			r.addPoly(bit|m, dp.Scale(float64(s)))
		}
	}
	return r
}

// D is a short alias for [PForm.ExteriorDerivative].
func (w *PForm) D() *PForm { return w.ExteriorDerivative() }

// HodgeStar returns the Euclidean Hodge dual ★w of the polynomial form, applying
// the blade-level Hodge map coefficient-wise. It maps grade k to grade n−k.
func (w *PForm) HodgeStar() *PForm {
	full := FullMask(w.n)
	r := NewPForm(w.n)
	for m, p := range w.terms {
		j := full &^ m
		r.addPoly(j, p.Scale(float64(reorderSign(m, j))))
	}
	return r
}

// Codifferential returns the codifferential δw with respect to the Euclidean
// metric, defined grade-by-grade on a grade-k component by
// δ = (−1)^{n(k+1)+1} ★d★. It lowers grade by one and is the formal adjoint of
// the exterior derivative.
func (w *PForm) Codifferential() *PForm {
	r := NewPForm(w.n)
	n := w.n
	for k := 1; k <= n; k++ {
		wk := w.GradeProject(k)
		if wk.IsZero() {
			continue
		}
		t := wk.HodgeStar().ExteriorDerivative().HodgeStar()
		if (n*(k+1)+1)&1 == 1 {
			t = t.Neg()
		}
		r = r.Add(t)
	}
	return r
}

// HodgeLaplacian returns the Hodge–de Rham Laplacian Δw = (dδ + δd)w with
// respect to the Euclidean metric. On a grade-0 form P it equals −Σ_i ∂²P/∂x_i²,
// the geometer's (positive-semidefinite) sign convention.
func (w *PForm) HodgeLaplacian() *PForm {
	return w.ExteriorDerivative().Codifferential().Add(w.Codifferential().ExteriorDerivative())
}

// Pullback returns the pullback φ*w of w along the polynomial map φ whose n
// component polynomials are given in phi, each of the same arity m. The result
// is a differential form on Rᵐ. Pullback is an algebra homomorphism that
// commutes with the exterior derivative: φ*(dw) = d(φ*w). It returns [ErrMap]
// unless len(phi) == n and all components share an arity.
func (w *PForm) Pullback(phi []*Poly) (*PForm, error) {
	if len(phi) != w.n {
		return nil, ErrMap
	}
	m := -1
	for _, c := range phi {
		if m == -1 {
			m = c.n
		} else if c.n != m {
			return nil, ErrMap
		}
	}
	if m == -1 {
		return nil, ErrMap
	}
	// dphi[i] is the 1-form d(phi^i) on R^m.
	dphi := make([]*PForm, w.n)
	for i, c := range phi {
		g := c.Gradient()
		f := NewPForm(m)
		for j := 0; j < m; j++ {
			if !g[j].IsZero() {
				f.terms[uint(1)<<uint(j)] = g[j]
			}
		}
		dphi[i] = f
	}
	res := NewPForm(m)
	for mask, p := range w.terms {
		coeff, err := p.Compose(phi)
		if err != nil {
			return nil, err
		}
		wedge := PConst(ConstPoly(m, 1))
		for _, i := range MaskToIndices(mask) {
			wedge = wedge.Wedge(dphi[i])
		}
		res = res.Add(wedge.ScalePoly(coeff))
	}
	return res, nil
}

// Eval evaluates every coefficient polynomial of w at the point x and returns
// the resulting constant [Form]. It panics if len(x) != n.
func (w *PForm) Eval(x []float64) *Form {
	if len(x) != w.n {
		panic(ErrDim)
	}
	f := New(w.n)
	for m, p := range w.terms {
		if v := p.Eval(x); v != 0 {
			f.terms[m] = v
		}
	}
	return f
}

// Equal reports whether w and u are the same polynomial form.
func (w *PForm) Equal(u *PForm) bool {
	if w.n != u.n || len(w.terms) != len(u.terms) {
		return false
	}
	for m, p := range w.terms {
		q, ok := u.terms[m]
		if !ok || !p.Equal(q) {
			return false
		}
	}
	return true
}

// EqualTol reports whether w and u have the same ambient dimension and every
// coefficient polynomial agrees to within the absolute tolerance tol.
func (w *PForm) EqualTol(u *PForm, tol float64) bool {
	if w.n != u.n {
		return false
	}
	zero := NewPoly(w.n)
	for m, p := range w.terms {
		q, ok := u.terms[m]
		if !ok {
			q = zero
		}
		if !p.EqualTol(q, tol) {
			return false
		}
	}
	for m, q := range u.terms {
		if _, ok := w.terms[m]; !ok && !q.EqualTol(zero, tol) {
			return false
		}
	}
	return true
}

// String renders w with each blade shown as its coefficient polynomial times
// the differential dx^I, for example "(2 x0) dx0 + (1) dx0∧dx1".
func (w *PForm) String() string {
	if w.IsZero() {
		return "0"
	}
	var b strings.Builder
	for i, m := range w.Masks() {
		if i > 0 {
			b.WriteString(" + ")
		}
		b.WriteString("(")
		b.WriteString(w.terms[m].String())
		b.WriteString(")")
		if m == 0 {
			continue
		}
		b.WriteString(" ")
		for j, ix := range MaskToIndices(m) {
			if j > 0 {
				b.WriteString("∧")
			}
			b.WriteString("dx")
			b.WriteString(itoa(ix))
		}
	}
	return b.String()
}

// itoa is a tiny helper avoiding a strconv import in this file.
func itoa(i int) string {
	if i == 0 {
		return "0"
	}
	neg := i < 0
	if neg {
		i = -i
	}
	var buf [20]byte
	pos := len(buf)
	for i > 0 {
		pos--
		buf[pos] = byte('0' + i%10)
		i /= 10
	}
	if neg {
		pos--
		buf[pos] = '-'
	}
	return string(buf[pos:])
}
