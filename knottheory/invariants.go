package knottheory

import "sort"

// This file gathers cross-cutting helpers and derived invariants that build on
// the polynomial, braid and diagram machinery of the package.

// SkeinDelta returns the loop value delta = -A^2 - A^{-2} of the Kauffman
// bracket skein relation, as a Laurent polynomial in A.
func SkeinDelta() Laurent { return Monomial(-1, 2).Add(Monomial(-1, -2)) }

// IsOne reports whether L is the constant polynomial 1.
func (L Laurent) IsOne() bool {
	return len(L.coeffs) == 1 && L.min == 0 && L.coeffs[0] == 1
}

// IsConstant reports whether L is a constant (a single term of degree zero, or
// the zero polynomial).
func (L Laurent) IsConstant() bool {
	if L.IsZero() {
		return true
	}
	return len(L.coeffs) == 1 && L.min == 0
}

// ConstantTerm returns the coefficient of X^0.
func (L Laurent) ConstantTerm() int { return L.Coeff(0) }

// MulMonomial returns L multiplied by the single term c*X^e.
func (L Laurent) MulMonomial(c, e int) Laurent { return L.Mul(Monomial(c, e)) }

// MaxAbsCoeff returns the largest absolute value among the coefficients of L,
// and 0 for the zero polynomial.
func (L Laurent) MaxAbsCoeff() int {
	m := 0
	for _, c := range L.coeffs {
		a := c
		if a < 0 {
			a = -a
		}
		if a > m {
			m = a
		}
	}
	return m
}

// CoeffSlice returns the minimum exponent and a copy of the dense coefficient
// slice of L, so that coeffs[i] is the coefficient of X^(min+i).
func (L Laurent) CoeffSlice() (min int, coeffs []int) {
	return L.min, append([]int(nil), L.coeffs...)
}

// EqualUpToShift reports whether L and other agree after shifting one of them by
// a power of X, that is whether they are equal as polynomials up to
// multiplication by a monomial X^k. This is the natural equality for the
// Alexander polynomial, which is only defined up to units.
func (L Laurent) EqualUpToShift(other Laurent) bool {
	if L.IsZero() || other.IsZero() {
		return L.IsZero() && other.IsZero()
	}
	shift := other.MinDegree() - L.MinDegree()
	return L.ShiftExp(shift).Equal(other)
}

// TranspositionCount returns the minimum number of transpositions needed to
// express the permutation, equal to size minus the number of cycles.
func (p Permutation) TranspositionCount() int { return len(p) - p.NumCycles() }

// Transpositions returns a decomposition of the permutation into transpositions.
// Reading the returned list left to right and composing with Compose reproduces
// the permutation.
func (p Permutation) Transpositions() [][2]int {
	var ts [][2]int
	for _, c := range p.Cycles() {
		for j := len(c) - 1; j >= 1; j-- {
			ts = append(ts, [2]int{c[0], c[j]})
		}
	}
	return ts
}

// Generators returns the sorted list of distinct generator indices that occur
// in the braid word (ignoring sign).
func (b Braid) Generators() []int {
	seen := map[int]bool{}
	var out []int
	for _, g := range b.word {
		a := g
		if a < 0 {
			a = -a
		}
		if !seen[a] {
			seen[a] = true
			out = append(out, a)
		}
	}
	sort.Ints(out)
	return out
}

// IsPureBraid reports whether the braid induces the identity permutation, so
// that every strand returns to its starting position.
func (b Braid) IsPureBraid() bool { return b.Permutation().IsIdentity() }

// Conjugate returns the braid g * b * g^{-1}. Conjugate braids have isotopic
// closures, so all closure invariants are preserved. It returns an error if the
// strand counts differ.
func (b Braid) Conjugate(g Braid) (Braid, error) {
	if b.strands != g.strands {
		return Braid{}, errStrandMismatch(b.strands, g.strands)
	}
	res, _ := g.Concat(b)
	res, _ = res.Concat(g.Inverse())
	return res, nil
}

// errStrandMismatch is a shared error constructor for strand-count mismatches.
func errStrandMismatch(a, c int) error {
	return &strandError{a: a, c: c}
}

type strandError struct{ a, c int }

// Error implements the error interface, reporting that two braids operate on
// different strand counts and cannot be composed or multiplied.
func (e *strandError) Error() string {
	return "knottheory: braids act on different strand counts (" +
		itoa(e.a) + " and " + itoa(e.c) + ")"
}

// itoa is a tiny base-10 integer formatter avoiding a strconv import here.
func itoa(n int) string {
	if n == 0 {
		return "0"
	}
	neg := n < 0
	if neg {
		n = -n
	}
	var buf [20]byte
	i := len(buf)
	for n > 0 {
		i--
		buf[i] = byte('0' + n%10)
		n /= 10
	}
	if neg {
		i--
		buf[i] = '-'
	}
	return string(buf[i:])
}

// ToDiagram wraps a single-component Gauss code as a one-component Diagram.
func (gc GaussCode) ToDiagram() Diagram { return Diagram{Components: []GaussCode{gc}} }

// IsKnot reports whether the diagram has exactly one component.
func (d Diagram) IsKnot() bool { return len(d.Components) == 1 }

// IsLink reports whether the diagram has more than one component.
func (d Diagram) IsLink() bool { return len(d.Components) > 1 }

// SignSequence returns the signs of the crossings of the PD diagram in order.
func (pd PDCode) SignSequence() []int {
	out := make([]int, len(pd.Crossings))
	for i, c := range pd.Crossings {
		out[i] = c.Sign
	}
	return out
}

// JonesSpan returns the breadth of the Jones polynomial (in the variable
// t^{1/2}), a lower bound of twice the crossing number for a reduced
// alternating diagram by the Kauffman-Murasugi-Thistlethwaite theorem.
func (pd PDCode) JonesSpan() int { return pd.JonesPolynomialSqrt().SpanWidth() }

// ArfInvariant returns the Arf invariant (0 or 1) of the knot, computed from the
// determinant: it is 0 when the determinant is congruent to +/-1 modulo 8 and 1
// when it is congruent to +/-3 modulo 8.
func (pd PDCode) ArfInvariant() int { return arfFromDeterminant(pd.KnotDeterminant()) }

// arfFromDeterminant maps a knot determinant to the Arf invariant using the
// congruence Delta(-1) = +/-1 (mod 8) <=> Arf 0.
func arfFromDeterminant(det int) int {
	m := ((det % 8) + 8) % 8
	if m == 1 || m == 7 {
		return 0
	}
	return 1
}

// TorusKnotArf returns the Arf invariant of the torus knot T(p, q) with coprime
// p, q, derived from its determinant.
func TorusKnotArf(p, q int) int { return arfFromDeterminant(TorusKnotDeterminant(p, q)) }

// TorusKnotBraidWord returns the generator word of the torus braid whose closure
// is T(p, q).
func TorusKnotBraidWord(p, q int) []int {
	b, err := TorusBraid(abs(p), q)
	if err != nil {
		return nil
	}
	return b.Word()
}

// IsTorusLink reports whether T(p, q) is a genuine link of more than one
// component (gcd(|p|,|q|) > 1).
func IsTorusLink(p, q int) bool { return TorusLinkComponents(p, q) > 1 }
