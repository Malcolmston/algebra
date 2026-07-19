package galois

import (
	"math/big"
	"sort"
)

// PolyFactor pairs an irreducible (or square-free) monic polynomial factor with
// the multiplicity at which it occurs.
type PolyFactor struct {
	Factor *Poly
	Exp    int
}

// PthRootPoly returns the p-th root of f, assuming f is a perfect p-th power in
// GF(p)[x] (its only non-zero coefficients occur at exponents divisible by p).
// Over a prime field the coefficient p-th roots equal the coefficients
// themselves by Fermat's little theorem.
func PthRootPoly(f *Poly) *Poly {
	if !f.P.IsInt64() {
		return f.Copy()
	}
	pi := int(f.P.Int64())
	deg := f.Degree()
	if deg < pi {
		return f.Copy()
	}
	newDeg := deg / pi
	coeffs := make([]*big.Int, newDeg+1)
	for i := 0; i <= newDeg; i++ {
		coeffs[i] = f.Coefficient(i * pi)
	}
	return NewPolyBig(f.P, coeffs)
}

func quoP(a, b *Poly) *Poly {
	q, _ := a.Quo(b)
	return q
}

func mergeFactors(in []PolyFactor) []PolyFactor {
	var out []PolyFactor
	for _, pf := range in {
		if pf.Factor.IsConstant() {
			continue
		}
		merged := false
		for i := range out {
			if out[i].Factor.Equal(pf.Factor) {
				out[i].Exp += pf.Exp
				merged = true
				break
			}
		}
		if !merged {
			out = append(out, PolyFactor{Factor: pf.Factor.Monic(), Exp: pf.Exp})
		}
	}
	sort.Slice(out, func(i, j int) bool {
		if out[i].Factor.Degree() != out[j].Factor.Degree() {
			return out[i].Factor.Degree() < out[j].Factor.Degree()
		}
		return out[i].Factor.ToIndex().Cmp(out[j].Factor.ToIndex()) < 0
	})
	return out
}

// ToIndex encodes a polynomial's coefficient vector as a big integer for stable
// ordering: Σ Coeff[i]·p^i.
func (a *Poly) ToIndex() *big.Int {
	acc := big.NewInt(0)
	for i := a.Degree(); i >= 0; i-- {
		acc.Mul(acc, a.P)
		acc.Add(acc, a.Coefficient(i))
	}
	return acc
}

func sqfreeRec(f *Poly, pIntVal int, out *[]PolyFactor) {
	f = f.Monic()
	fp := f.Derivative()
	if !fp.IsZero() {
		c := f.Gcd(fp)
		w := quoP(f, c)
		i := 1
		for !w.IsOne() {
			y := w.Gcd(c)
			z := quoP(w, y)
			if !z.IsOne() {
				*out = append(*out, PolyFactor{Factor: z, Exp: i})
			}
			i++
			w = y
			c = quoP(c, y)
		}
		if !c.IsConstant() && pIntVal > 0 {
			root := PthRootPoly(c)
			var sub []PolyFactor
			sqfreeRec(root, pIntVal, &sub)
			for _, pf := range sub {
				*out = append(*out, PolyFactor{Factor: pf.Factor, Exp: pf.Exp * pIntVal})
			}
		}
	} else if !f.IsConstant() && pIntVal > 0 {
		root := PthRootPoly(f)
		var sub []PolyFactor
		sqfreeRec(root, pIntVal, &sub)
		for _, pf := range sub {
			*out = append(*out, PolyFactor{Factor: pf.Factor, Exp: pf.Exp * pIntVal})
		}
	}
}

// SquareFreeFactorization returns the square-free decomposition of f: a set of
// pairwise-coprime square-free monic factors together with their multiplicities,
// whose product (each raised to its multiplicity) equals the monic part of f.
func SquareFreeFactorization(f *Poly) []PolyFactor {
	pIntVal := 0
	if f.P.IsInt64() {
		pIntVal = int(f.P.Int64())
	}
	var res []PolyFactor
	sqfreeRec(f.Monic(), pIntVal, &res)
	return mergeFactors(res)
}

// IsSquareFree reports whether f has no repeated irreducible factor.
func IsSquareFree(f *Poly) bool {
	if f.IsConstant() {
		return true
	}
	fp := f.Derivative()
	if fp.IsZero() {
		return false
	}
	return f.Gcd(fp).Degree() == 0
}

// splitWith partitions the monic square-free polynomial h using a Berlekamp
// splitting polynomial g (satisfying g^p ≡ g mod the parent), returning the
// non-constant gcd(h, g - s) over all s in GF(p).
func splitWith(h, g *Poly) []*Poly {
	p := h.P
	var pieces []*Poly
	for s := big.NewInt(0); s.Cmp(p) < 0; s.Add(s, big1) {
		d := h.Gcd(g.Sub(PolyConst(p, s)))
		if d.Degree() >= 1 {
			pieces = append(pieces, d.Monic())
		}
	}
	return pieces
}

// berlekamp returns the distinct monic irreducible factors of the monic
// square-free polynomial f using Berlekamp's deterministic algorithm.
func berlekamp(f *Poly) []*Poly {
	f = f.Monic()
	n := f.Degree()
	if n <= 1 {
		return []*Poly{f}
	}
	p := f.P
	xp, _ := PolyX(p).PowMod(p, f) // x^p mod f
	Q := NewMatModP(p, n, n)
	cur := PolyOne(p)
	for i := 0; i < n; i++ {
		for j := 0; j <= cur.Degree(); j++ {
			Q.Data[i][j] = cur.Coefficient(j)
		}
		cur = PolyMulMod(cur, xp, f)
	}
	B := Q.Sub(IdentityModP(p, n))
	basis := B.Transpose().NullSpace()
	r := len(basis)
	if r <= 1 {
		return []*Poly{f}
	}
	factors := []*Poly{f}
	for _, vec := range basis {
		g := NewPolyBig(p, vec)
		if g.Degree() < 1 {
			continue
		}
		var next []*Poly
		for _, h := range factors {
			if h.Degree() <= 1 {
				next = append(next, h)
				continue
			}
			pieces := splitWith(h, g)
			if len(pieces) <= 1 {
				next = append(next, h)
			} else {
				next = append(next, pieces...)
			}
		}
		factors = next
		if len(factors) >= r {
			break
		}
	}
	return factors
}

// Factor returns the complete factorisation of f over GF(p): the leading
// coefficient together with the distinct monic irreducible factors and their
// multiplicities. The product of the factors (each raised to its multiplicity),
// scaled by the leading coefficient, reconstructs f.
func Factor(f *Poly) (lead *big.Int, factors []PolyFactor) {
	lead = f.LeadingCoeff()
	if f.IsZero() {
		return big.NewInt(0), nil
	}
	sq := SquareFreeFactorization(f)
	var res []PolyFactor
	for _, sf := range sq {
		for _, ir := range berlekamp(sf.Factor) {
			res = append(res, PolyFactor{Factor: ir, Exp: sf.Exp})
		}
	}
	return lead, mergeFactors(res)
}

// FactorCount returns the number of irreducible factors of f counted with
// multiplicity.
func FactorCount(f *Poly) int {
	_, factors := Factor(f)
	total := 0
	for _, pf := range factors {
		total += pf.Exp
	}
	return total
}

// DistinctFactorCount returns the number of distinct irreducible factors of f.
func DistinctFactorCount(f *Poly) int {
	_, factors := Factor(f)
	return len(factors)
}
