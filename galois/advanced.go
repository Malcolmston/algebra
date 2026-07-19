package galois

import (
	"errors"
	"math/big"
)

// Frobenius returns the image of a under the Frobenius automorphism x ↦ x^p,
// the canonical generator of the Galois group of GF(p^n) over GF(p).
func (a *Element) Frobenius() *Element {
	r, _ := a.Pow(a.Field.P)
	return r
}

// FrobeniusPow returns the k-fold Frobenius image a^(p^k). Negative k are
// reduced modulo n, matching the cyclic Galois group of order n.
func (a *Element) FrobeniusPow(k int) *Element {
	n := a.Field.N
	k = ((k % n) + n) % n
	exp := IntPow(a.Field.P, k)
	r, _ := a.Pow(exp)
	return r
}

// Conjugates returns the Galois conjugates of a, namely a, a^p, a^(p^2), … up
// to but not including the first repetition. Their count is the degree of a's
// minimal polynomial over GF(p).
func (a *Element) Conjugates() []*Element {
	var conj []*Element
	cur := a.Copy()
	for {
		dup := false
		for _, c := range conj {
			if c.Equal(cur) {
				dup = true
				break
			}
		}
		if dup {
			break
		}
		conj = append(conj, cur)
		cur = cur.Frobenius()
	}
	return conj
}

// Trace returns the field trace Tr(a) = a + a^p + … + a^(p^(n-1)) of a down to
// the prime field GF(p), returned as an integer in [0, p).
func (a *Element) Trace() *big.Int {
	acc := a.Field.Zero()
	cur := a.Copy()
	for i := 0; i < a.Field.N; i++ {
		acc = acc.Add(cur)
		cur = cur.Frobenius()
	}
	return acc.Poly.Coefficient(0)
}

// Norm returns the field norm N(a) = a^((p^n-1)/(p-1)) = a · a^p · … · a^(p^(n-1))
// of a down to the prime field GF(p), returned as an integer in [0, p).
func (a *Element) Norm() *big.Int {
	acc := a.Field.One()
	cur := a.Copy()
	for i := 0; i < a.Field.N; i++ {
		acc = acc.Mul(cur)
		cur = cur.Frobenius()
	}
	return acc.Poly.Coefficient(0)
}

// MinimalPoly returns the minimal polynomial of a over GF(p): the unique monic
// polynomial of least degree with a as a root, equal to the product of
// (x - c) over the distinct conjugates c of a.
func (a *Element) MinimalPoly() *Poly {
	conj := a.Conjugates()
	f := a.Field
	// Build the product of (x - c) with coefficients that are field elements.
	coeffs := []*Element{f.One()}
	for _, c := range conj {
		next := make([]*Element, len(coeffs)+1)
		for i := range next {
			next[i] = f.Zero()
		}
		for i, coef := range coeffs {
			next[i+1] = next[i+1].Add(coef)    // x · coef
			next[i] = next[i].Sub(c.Mul(coef)) // -c · coef
		}
		coeffs = next
	}
	out := make([]*big.Int, len(coeffs))
	for i, coef := range coeffs {
		out[i] = coef.Poly.Coefficient(0)
	}
	return NewPolyBig(f.P, out)
}

// Order returns the multiplicative order of the non-zero element a: the least
// positive k with a^k = 1. It returns an error when a is zero.
func (a *Element) Order() (*big.Int, error) {
	if a.IsZero() {
		return nil, errors.New("galois: zero has no multiplicative order")
	}
	e := new(big.Int).Sub(a.Field.Order(), big1)
	order := clone(e)
	for _, r := range PrimeFactors(e) {
		for new(big.Int).Mod(order, r).Sign() == 0 {
			cand := new(big.Int).Div(order, r)
			ap, _ := a.Pow(cand)
			if ap.IsOne() {
				order = cand
			} else {
				break
			}
		}
	}
	return order, nil
}

// IsPrimitive reports whether a generates the whole multiplicative group of the
// field, i.e. its order equals p^n - 1.
func (a *Element) IsPrimitive() bool {
	if a.IsZero() {
		return false
	}
	e := new(big.Int).Sub(a.Field.Order(), big1)
	for _, r := range PrimeFactors(e) {
		exp := new(big.Int).Div(e, r)
		ap, _ := a.Pow(exp)
		if ap.IsOne() {
			return false
		}
	}
	return true
}

// PrimitiveElement returns a generator of the multiplicative group of the
// field: an element of order p^n - 1. It scans elements by index and returns
// the first primitive one.
func (f *Field) PrimitiveElement() (*Element, error) {
	ord := f.Order()
	if !ord.IsInt64() {
		return nil, errors.New("galois: field too large to search for a primitive element")
	}
	q := ord.Int64()
	for i := int64(1); i < q; i++ {
		a := f.ElementFromInt(big.NewInt(i))
		if a.IsPrimitive() {
			return a, nil
		}
	}
	return nil, errors.New("galois: no primitive element found")
}

// Generator is a synonym for PrimitiveElement, returning a multiplicative
// generator of the field.
func (f *Field) Generator() (*Element, error) { return f.PrimitiveElement() }

// FrobeniusAutomorphism returns the Frobenius map x ↦ x^p of the field as a Go
// function.
func (f *Field) FrobeniusAutomorphism() func(*Element) *Element {
	return func(a *Element) *Element { return a.Frobenius() }
}

// FieldTrace returns the trace of a to the prime field. It is a package-level
// wrapper over (*Element).Trace.
func FieldTrace(a *Element) *big.Int { return a.Trace() }

// FieldNorm returns the norm of a to the prime field. It is a package-level
// wrapper over (*Element).Norm.
func FieldNorm(a *Element) *big.Int { return a.Norm() }

// MinimalPolynomial returns the minimal polynomial of a over GF(p). It is a
// package-level wrapper over (*Element).MinimalPoly.
func MinimalPolynomial(a *Element) *Poly { return a.MinimalPoly() }

// ConjugateElements returns the Galois conjugates of a. It is a package-level
// wrapper over (*Element).Conjugates.
func ConjugateElements(a *Element) []*Element { return a.Conjugates() }

// isqrtCeil returns the ceiling of the integer square root of n (n >= 0).
func isqrtCeil(n *big.Int) *big.Int {
	if n.Sign() <= 0 {
		return big.NewInt(0)
	}
	s := new(big.Int).Sqrt(n)
	if new(big.Int).Mul(s, s).Cmp(n) < 0 {
		s.Add(s, big1)
	}
	return s
}

// DiscreteLog returns the discrete logarithm x of target to the base, the least
// non-negative x with base^x = target, computed by the baby-step giant-step
// algorithm. It returns an error when the fields differ or no logarithm exists
// (target is not a power of base).
func DiscreteLog(base, target *Element) (*big.Int, error) {
	if !sameField(base, target) {
		return nil, errors.New("galois: DiscreteLog requires elements of the same field")
	}
	if base.IsZero() {
		return nil, errors.New("galois: DiscreteLog base must be non-zero")
	}
	ord, err := base.Order()
	if err != nil {
		return nil, err
	}
	m := isqrtCeil(ord)
	if !m.IsInt64() {
		return nil, errors.New("galois: field too large for baby-step giant-step")
	}
	mi := m.Int64()
	table := make(map[string]int64, mi)
	cur := base.Field.One()
	for j := int64(0); j < mi; j++ {
		key := cur.ToInt().String()
		if _, ok := table[key]; !ok {
			table[key] = j
		}
		cur = cur.Mul(base)
	}
	bm, err := base.Pow(m)
	if err != nil {
		return nil, err
	}
	factor, err := bm.Inv()
	if err != nil {
		return nil, err
	}
	gamma := target.Copy()
	for i := int64(0); i < mi; i++ {
		if j, ok := table[gamma.ToInt().String()]; ok {
			x := new(big.Int).Add(new(big.Int).Mul(big.NewInt(i), m), big.NewInt(j))
			return x.Mod(x, ord), nil
		}
		gamma = gamma.Mul(factor)
	}
	return nil, errors.New("galois: no discrete logarithm exists")
}

// BabyStepGiantStep is an alias for DiscreteLog.
func BabyStepGiantStep(base, target *Element) (*big.Int, error) {
	return DiscreteLog(base, target)
}

// Subfields returns the degrees d dividing n for which GF(p^d) is a subfield of
// the field GF(p^n), in increasing order. GF(p) (d=1) and the field itself
// (d=n) are always included.
func (f *Field) Subfields() []int {
	return divisorsInt64(f.N)
}

// SubfieldDegrees returns the degrees of the subfields of GF(p^n): the positive
// divisors of n.
func SubfieldDegrees(n int) []int {
	return divisorsInt64(n)
}

// IsSubfieldDegree reports whether GF(p^d) embeds in GF(p^n), i.e. d divides n.
func IsSubfieldDegree(n, d int) bool {
	return d > 0 && n%d == 0
}

// SubfieldOrders returns the orders p^d of all subfields of the field, in
// increasing order.
func (f *Field) SubfieldOrders() []*big.Int {
	var out []*big.Int
	for _, d := range divisorsInt64(f.N) {
		out = append(out, IntPow(f.P, d))
	}
	return out
}
