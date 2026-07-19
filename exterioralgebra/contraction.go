package exterioralgebra

import "math"

// interiorBladeVec contracts the single basis covector e_i into the blade m.
// It returns the resulting blade mask, a sign of ±1, and ok reporting whether
// e_i actually occurs in m (a zero result otherwise).
func interiorBladeVec(i int, m uint) (uint, int, bool) {
	bit := uint(1) << uint(i)
	if m&bit == 0 {
		return 0, 0, false
	}
	p := Popcount(m & (bit - 1)) // number of earlier factors to hop over
	sign := 1
	if p&1 == 1 {
		sign = -1
	}
	return m &^ bit, sign, true
}

// InteriorProduct returns the interior product (contraction) ι_v ω of the
// vector v into the Form omega, lowering grade by one. Only the grade-1 part of
// v is used, and ι_v = Σ_i v_i ι_{e_i}. The interior product is an
// antiderivation: ι_v(α∧β) = (ι_v α)∧β + (−1)^{deg α} α∧(ι_v β). It panics if
// the ambient dimensions differ.
func InteriorProduct(v, omega *Form) *Form {
	v.requireSameDim(omega)
	res := New(omega.n)
	for i := 0; i < v.n; i++ {
		vi := v.terms[uint(1)<<uint(i)]
		if vi == 0 {
			continue
		}
		for m, c := range omega.terms {
			if nm, s, ok := interiorBladeVec(i, m); ok {
				res.addTerm(nm, vi*float64(s)*c)
			}
		}
	}
	return res
}

// Interior is the method form of [InteriorProduct]: omega.Interior(v) returns
// ι_v omega.
func (f *Form) Interior(v *Form) *Form { return InteriorProduct(v, f) }

// leftContractBlade computes e_a ⌋ e_b under the Euclidean metric. It returns
// the resulting blade mask and a sign, with ok false when a is not a subset of
// b (the contraction vanishes).
func leftContractBlade(a, b uint) (uint, int, bool) {
	if a&^b != 0 {
		return 0, 0, false
	}
	c := b &^ a
	return c, reorderSign(a, c), true
}

// rightContractBlade computes e_b ⌊ e_a under the Euclidean metric.
func rightContractBlade(b, a uint) (uint, int, bool) {
	if a&^b != 0 {
		return 0, 0, false
	}
	c := b &^ a
	return c, reorderSign(c, a), true
}

// LeftContract returns the left contraction f ⌋ g under the Euclidean metric,
// the bilinear product characterised by ⟨f⌋g, h⟩ = ⟨g, f∧h⟩. On homogeneous
// grades r and s it produces grade s−r. It panics if the dimensions differ.
func (f *Form) LeftContract(g *Form) *Form {
	f.requireSameDim(g)
	res := New(f.n)
	for ma, ca := range f.terms {
		for mb, cb := range g.terms {
			if nm, s, ok := leftContractBlade(ma, mb); ok {
				res.addTerm(nm, float64(s)*ca*cb)
			}
		}
	}
	return res
}

// RightContract returns the right contraction f ⌊ g under the Euclidean metric,
// characterised by ⟨f⌊g, h⟩ = ⟨f, h∧g⟩. On homogeneous grades r and s it
// produces grade r−s. It panics if the dimensions differ.
func (f *Form) RightContract(g *Form) *Form {
	f.requireSameDim(g)
	res := New(f.n)
	for ma, ca := range f.terms {
		for mb, cb := range g.terms {
			if nm, s, ok := rightContractBlade(ma, mb); ok {
				res.addTerm(nm, float64(s)*ca*cb)
			}
		}
	}
	return res
}

// InnerProduct returns the Euclidean inner product ⟨f, g⟩ obtained by declaring
// the basis blades orthonormal: it is the sum over blades of the products of
// matching coefficients. It panics if the ambient dimensions differ.
func InnerProduct(f, g *Form) float64 {
	f.requireSameDim(g)
	// iterate over the smaller map for efficiency
	a, b := f, g
	if len(b.terms) < len(a.terms) {
		a, b = b, a
	}
	var s float64
	for m, c := range a.terms {
		s += c * b.terms[m]
	}
	return s
}

// Dot is an alias for [InnerProduct].
func Dot(f, g *Form) float64 { return InnerProduct(f, g) }

// NormSq returns the squared Euclidean norm ⟨f, f⟩ of f.
func (f *Form) NormSq() float64 { return InnerProduct(f, f) }

// Norm returns the Euclidean norm √⟨f, f⟩ of f.
func (f *Form) Norm() float64 { return math.Sqrt(f.NormSq()) }

// Normalize returns f scaled to unit Euclidean norm, together with the original
// norm. If f is zero it returns a zero Form and norm 0.
func (f *Form) Normalize() (*Form, float64) {
	nrm := f.Norm()
	if nrm == 0 {
		return New(f.n), 0
	}
	return f.Scale(1 / nrm), nrm
}

// Angle returns the angle in radians between two nonzero homogeneous Forms of
// the same grade, defined via the Euclidean inner product. It returns NaN if
// either operand is zero.
func Angle(f, g *Form) float64 {
	nf, ng := f.Norm(), g.Norm()
	if nf == 0 || ng == 0 {
		return math.NaN()
	}
	c := InnerProduct(f, g) / (nf * ng)
	if c > 1 {
		c = 1
	} else if c < -1 {
		c = -1
	}
	return math.Acos(c)
}
