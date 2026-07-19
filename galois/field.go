package galois

import (
	"errors"
	"fmt"
	"math/big"
	"math/rand"
)

// Field represents the finite field GF(p^n). It fixes the prime characteristic
// P, the extension degree N, and a monic irreducible modulus polynomial Mod of
// degree N over GF(P). Field elements are residues of GF(P)[x] modulo Mod.
type Field struct {
	P   *big.Int
	N   int
	Mod *Poly
}

// Element is a member of a Field, represented by its unique residue polynomial
// of degree less than the field degree.
type Element struct {
	Field *Field
	Poly  *Poly
}

// NewPrimeField returns the prime field GF(p). It returns an error when p is
// not prime.
func NewPrimeField(p *big.Int) (*Field, error) {
	if !IsPrimeInt(p) {
		return nil, errors.New("galois: NewPrimeField requires a prime")
	}
	return &Field{P: clone(p), N: 1, Mod: PolyX(p)}, nil
}

// NewField returns GF(p^n) defined by the given monic modulus polynomial, which
// must be irreducible of degree n >= 1 over GF(p). It returns an error
// otherwise.
func NewField(p *big.Int, mod *Poly) (*Field, error) {
	if !IsPrimeInt(p) {
		return nil, errors.New("galois: NewField requires a prime characteristic")
	}
	if mod.P.Cmp(p) != 0 {
		return nil, errors.New("galois: modulus is defined over a different field")
	}
	if mod.Degree() < 1 {
		return nil, errors.New("galois: modulus must have positive degree")
	}
	if !IsIrreducible(mod) {
		return nil, errors.New("galois: modulus polynomial is not irreducible")
	}
	m := mod.Monic()
	return &Field{P: clone(p), N: m.Degree(), Mod: m}, nil
}

// NewFieldFromDegree returns GF(p^n) using an automatically discovered monic
// irreducible modulus of degree n.
func NewFieldFromDegree(p *big.Int, n int) (*Field, error) {
	if n < 1 {
		return nil, errors.New("galois: degree must be at least 1")
	}
	if n == 1 {
		return NewPrimeField(p)
	}
	mod, err := FindIrreducible(p, n)
	if err != nil {
		return nil, err
	}
	return &Field{P: clone(p), N: n, Mod: mod}, nil
}

// Order returns the number of elements in the field, p^n.
func (f *Field) Order() *big.Int { return IntPow(f.P, f.N) }

// Characteristic returns the prime characteristic p of the field.
func (f *Field) Characteristic() *big.Int { return clone(f.P) }

// Degree returns the extension degree n, so the field is GF(p^n).
func (f *Field) Degree() int { return f.N }

// Modulus returns a copy of the field's defining irreducible polynomial.
func (f *Field) Modulus() *Poly { return f.Mod.Copy() }

// reduce returns p reduced modulo the field's defining polynomial.
func (f *Field) reduce(p *Poly) *Poly {
	_, r, err := p.DivMod(f.Mod)
	if err != nil {
		return PolyZero(f.P)
	}
	return r
}

// FromPoly returns the field element represented by the polynomial p reduced
// modulo the field modulus.
func (f *Field) FromPoly(p *Poly) *Element {
	return &Element{Field: f, Poly: f.reduce(p)}
}

// Element returns the field element for the given little-endian int64
// coefficients over GF(p), reduced into the field.
func (f *Field) Element(coeffs ...int64) *Element {
	return f.FromPoly(NewPoly(f.P, coeffs...))
}

// ElementFromCoeffs returns the field element for the given little-endian
// big.Int coefficients, reduced into the field.
func (f *Field) ElementFromCoeffs(coeffs []*big.Int) *Element {
	return f.FromPoly(NewPolyBig(f.P, coeffs))
}

// ElementFromInt maps an integer index k in [0, p^n) to a field element by
// writing k in base p as the element's coefficient vector. It is the inverse of
// (*Element).ToInt.
func (f *Field) ElementFromInt(k *big.Int) *Element {
	kk := new(big.Int).Mod(k, f.Order())
	coeffs := make([]*big.Int, f.N)
	for i := 0; i < f.N; i++ {
		q, r := new(big.Int).DivMod(kk, f.P, new(big.Int))
		coeffs[i] = r
		kk = q
	}
	return f.ElementFromCoeffs(coeffs)
}

// Zero returns the additive identity of the field.
func (f *Field) Zero() *Element { return &Element{Field: f, Poly: PolyZero(f.P)} }

// One returns the multiplicative identity of the field.
func (f *Field) One() *Element { return &Element{Field: f, Poly: PolyOne(f.P)} }

// PrimitiveElementX returns the residue class of x, the canonical generator of
// GF(p^n) as GF(p)[x]/(mod). It is a field generator but not necessarily a
// generator of the multiplicative group.
func (f *Field) PrimitiveElementX() *Element {
	if f.N == 1 {
		// x reduces to 0 here; use 1 as the trivial generator of GF(p).
		return f.One()
	}
	return f.FromPoly(PolyX(f.P))
}

// sameField reports whether two elements live in the same field instance shape.
func sameField(a, b *Element) bool {
	return a.Field.P.Cmp(b.Field.P) == 0 && a.Field.Mod.Equal(b.Field.Mod)
}

// Add returns a + b in the field.
func (a *Element) Add(b *Element) *Element {
	return &Element{Field: a.Field, Poly: a.Field.reduce(a.Poly.Add(b.Poly))}
}

// Sub returns a - b in the field.
func (a *Element) Sub(b *Element) *Element {
	return &Element{Field: a.Field, Poly: a.Field.reduce(a.Poly.Sub(b.Poly))}
}

// Neg returns the additive inverse of a.
func (a *Element) Neg() *Element {
	return &Element{Field: a.Field, Poly: a.Field.reduce(a.Poly.Neg())}
}

// Mul returns a * b in the field.
func (a *Element) Mul(b *Element) *Element {
	return &Element{Field: a.Field, Poly: a.Field.reduce(a.Poly.Mul(b.Poly))}
}

// Inv returns the multiplicative inverse of a, or an error when a is zero.
func (a *Element) Inv() (*Element, error) {
	if a.IsZero() {
		return nil, errors.New("galois: zero has no multiplicative inverse")
	}
	g, s, _ := a.Poly.ExtendedGcd(a.Field.Mod)
	if g.Degree() != 0 {
		return nil, errors.New("galois: element is not invertible")
	}
	inv, err := InvMod(g.Coefficient(0), a.Field.P)
	if err != nil {
		return nil, err
	}
	return &Element{Field: a.Field, Poly: a.Field.reduce(s.ScalarMul(inv))}, nil
}

// Div returns a / b in the field, or an error when b is zero.
func (a *Element) Div(b *Element) (*Element, error) {
	binv, err := b.Inv()
	if err != nil {
		return nil, err
	}
	return a.Mul(binv), nil
}

// Pow returns a raised to the integer power e. Negative exponents invert a
// first and require a to be non-zero.
func (a *Element) Pow(e *big.Int) (*Element, error) {
	if e.Sign() == 0 {
		return a.Field.One(), nil
	}
	base := a
	ee := clone(e)
	if ee.Sign() < 0 {
		inv, err := a.Inv()
		if err != nil {
			return nil, err
		}
		base = inv
		ee.Neg(ee)
	}
	result := a.Field.One()
	for ee.Sign() > 0 {
		if ee.Bit(0) == 1 {
			result = result.Mul(base)
		}
		ee.Rsh(ee, 1)
		if ee.Sign() > 0 {
			base = base.Mul(base)
		}
	}
	return result, nil
}

// Equal reports whether a and b are the same element of the same field.
func (a *Element) Equal(b *Element) bool {
	return sameField(a, b) && a.Poly.Equal(b.Poly)
}

// IsZero reports whether a is the additive identity.
func (a *Element) IsZero() bool { return a.Poly.IsZero() }

// IsOne reports whether a is the multiplicative identity.
func (a *Element) IsOne() bool { return a.Poly.IsOne() }

// Copy returns an independent copy of the element.
func (a *Element) Copy() *Element {
	return &Element{Field: a.Field, Poly: a.Poly.Copy()}
}

// Coeffs returns the little-endian coefficient vector of the element, padded to
// the field degree n.
func (a *Element) Coeffs() []*big.Int {
	c := make([]*big.Int, a.Field.N)
	for i := 0; i < a.Field.N; i++ {
		c[i] = a.Poly.Coefficient(i)
	}
	return c
}

// ToInt returns the integer index in [0, p^n) obtained by reading the element's
// coefficient vector as base-p digits. It is the inverse of ElementFromInt.
func (a *Element) ToInt() *big.Int {
	acc := big.NewInt(0)
	for i := a.Field.N - 1; i >= 0; i-- {
		acc.Mul(acc, a.Field.P)
		acc.Add(acc, a.Poly.Coefficient(i))
	}
	return acc
}

// String renders the element via its residue polynomial.
func (a *Element) String() string {
	return a.Poly.String()
}

// Field-level arithmetic wrappers ------------------------------------------------

// Add returns a + b, a convenience wrapper over Element.Add.
func (f *Field) Add(a, b *Element) *Element { return a.Add(b) }

// Sub returns a - b.
func (f *Field) Sub(a, b *Element) *Element { return a.Sub(b) }

// Mul returns a * b.
func (f *Field) Mul(a, b *Element) *Element { return a.Mul(b) }

// Neg returns -a.
func (f *Field) Neg(a *Element) *Element { return a.Neg() }

// Inv returns the inverse of a.
func (f *Field) Inv(a *Element) (*Element, error) { return a.Inv() }

// Div returns a / b.
func (f *Field) Div(a, b *Element) (*Element, error) { return a.Div(b) }

// Pow returns a^e.
func (f *Field) Pow(a *Element, e *big.Int) (*Element, error) { return a.Pow(e) }

// Elements enumerates every element of the field in index order. It is intended
// for small fields; the slice has p^n entries.
func (f *Field) Elements() []*Element {
	ord := f.Order()
	if !ord.IsInt64() {
		return nil
	}
	n := ord.Int64()
	out := make([]*Element, 0, n)
	for i := int64(0); i < n; i++ {
		out = append(out, f.ElementFromInt(big.NewInt(i)))
	}
	return out
}

// NonzeroElements enumerates every non-zero element of the field in index
// order. It is intended for small fields.
func (f *Field) NonzeroElements() []*Element {
	all := f.Elements()
	if len(all) == 0 {
		return nil
	}
	return all[1:]
}

// RandomElement returns a uniformly random element drawn from r.
func (f *Field) RandomElement(r *rand.Rand) *Element {
	k := new(big.Int).Rand(r, f.Order())
	return f.ElementFromInt(k)
}

// String describes the field as "GF(p^n)".
func (f *Field) String() string {
	if f.N == 1 {
		return fmt.Sprintf("GF(%s)", f.P.String())
	}
	return fmt.Sprintf("GF(%s^%d)", f.P.String(), f.N)
}
