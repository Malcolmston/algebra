package grouprep

import (
	"errors"
	"math"
	"math/cmplx"
)

// Character is a class function on a group, stored as one complex value per
// element index: chi[i] is the value on element i. The character of a
// representation is the class function i ↦ tr ρ(i), obtained from
// [Rep.Character].
type Character []complex128

// CharacterFromRep returns the character of the representation r; it is a free
// function equivalent to r.Character().
func CharacterFromRep(r *Rep) Character { return r.Character() }

// Degree returns the degree (dimension) of the character, its value on the
// identity element, as a real number.
func (chi Character) Degree() float64 {
	return real(chi[0])
}

// Value returns the value of chi on element i.
func (chi Character) Value(i int) complex128 { return chi[i] }

// Clone returns an independent copy of chi.
func (chi Character) Clone() Character {
	return append(Character(nil), chi...)
}

// IsClassFunction reports whether chi is constant on the conjugacy classes of g
// to within tol. Every character is a class function; this checks a candidate.
func (chi Character) IsClassFunction(g *Group, tol float64) bool {
	for _, cls := range g.ConjugacyClasses() {
		v := chi[cls[0]]
		for _, e := range cls {
			if cmplx.Abs(chi[e]-v) > tol {
				return false
			}
		}
	}
	return true
}

// InnerProduct returns the Hermitian inner product of two characters of g,
//
//	〈a, b〉 = (1/|G|) Σ_{x∈G} conj(a(x)) · b(x).
//
// For irreducible characters this is 1 when the characters are equal and 0
// otherwise (the first orthogonality relation). It panics if the lengths do not
// match |G|.
func InnerProduct(g *Group, a, b Character) complex128 {
	n := g.Order()
	if len(a) != n || len(b) != n {
		panic("grouprep: InnerProduct character length mismatch")
	}
	var s complex128
	for x := 0; x < n; x++ {
		s += cmplx.Conj(a[x]) * b[x]
	}
	return s / complex(float64(n), 0)
}

// NormSquared returns the real number 〈chi, chi〉, the squared norm of chi in
// the character inner product. It equals the sum of the squares of the
// multiplicities of the irreducible constituents.
func NormSquared(g *Group, chi Character) float64 {
	return real(InnerProduct(g, chi, chi))
}

// IsIrreducible reports whether chi is an irreducible character, i.e. its norm
// squared is 1 to within tol.
func IsIrreducible(g *Group, chi Character, tol float64) bool {
	return math.Abs(NormSquared(g, chi)-1) <= tol
}

// Multiplicity returns the multiplicity of the irreducible character irr in the
// character chi, the rounded value of 〈irr, chi〉.
func Multiplicity(g *Group, chi, irr Character) int {
	ip := InnerProduct(g, irr, chi)
	return int(math.Round(real(ip)))
}

// DecomposeCharacter returns the multiplicities of chi against the given list of
// irreducible characters, so that chi = Σ mult[i]·irrs[i]. This is the
// numerical form of complete reducibility (Maschke's theorem). The multiplicity
// is the rounded inner product 〈irrs[i], chi〉.
func DecomposeCharacter(g *Group, chi Character, irrs []Character) []int {
	out := make([]int, len(irrs))
	for i, irr := range irrs {
		out[i] = Multiplicity(g, chi, irr)
	}
	return out
}

// IsCharacter reports whether chi is a (non-virtual) character with respect to
// the given complete list of irreducibles: all multiplicities are non-negative
// integers to within tol, and the reconstruction matches chi.
func IsCharacter(g *Group, chi Character, irrs []Character, tol float64) bool {
	mult := make([]float64, len(irrs))
	for i, irr := range irrs {
		ip := InnerProduct(g, irr, chi)
		if math.Abs(imag(ip)) > tol {
			return false
		}
		r := real(ip)
		if math.Abs(r-math.Round(r)) > tol || math.Round(r) < 0 {
			return false
		}
		mult[i] = math.Round(r)
	}
	recon := ZeroCharacter(g)
	for i, irr := range irrs {
		recon = AddCharacters(recon, ScaleCharacter(complex(mult[i], 0), irr))
	}
	return recon.ApproxEqual(chi, tol)
}

// ApproxEqual reports whether chi and psi agree entrywise to within tol.
func (chi Character) ApproxEqual(psi Character, tol float64) bool {
	if len(chi) != len(psi) {
		return false
	}
	for i := range chi {
		if cmplx.Abs(chi[i]-psi[i]) > tol {
			return false
		}
	}
	return true
}

// ZeroCharacter returns the identically zero class function on g.
func ZeroCharacter(g *Group) Character {
	return make(Character, g.Order())
}

// TrivialCharacter returns the character of the trivial representation of g,
// identically 1.
func TrivialCharacter(g *Group) Character {
	chi := make(Character, g.Order())
	for i := range chi {
		chi[i] = 1
	}
	return chi
}

// RegularCharacter returns the character of the regular representation of g:
// |G| on the identity and 0 elsewhere. It decomposes as Σ dᵢ·χᵢ over the
// irreducibles, each appearing with multiplicity equal to its degree.
func RegularCharacter(g *Group) Character {
	chi := make(Character, g.Order())
	chi[0] = complex(float64(g.Order()), 0)
	return chi
}

// AddCharacters returns the pointwise sum a+b, the character of the direct sum.
func AddCharacters(a, b Character) Character {
	out := make(Character, len(a))
	for i := range a {
		out[i] = a[i] + b[i]
	}
	return out
}

// SubCharacters returns the pointwise difference a-b (a virtual character).
func SubCharacters(a, b Character) Character {
	out := make(Character, len(a))
	for i := range a {
		out[i] = a[i] - b[i]
	}
	return out
}

// ScaleCharacter returns z·chi.
func ScaleCharacter(z complex128, chi Character) Character {
	out := make(Character, len(chi))
	for i := range chi {
		out[i] = z * chi[i]
	}
	return out
}

// TensorCharacters returns the pointwise product a·b, the character of the
// tensor product representation.
func TensorCharacters(a, b Character) Character {
	out := make(Character, len(a))
	for i := range a {
		out[i] = a[i] * b[i]
	}
	return out
}

// ConjugateCharacter returns the complex conjugate of chi, the character of the
// dual representation.
func ConjugateCharacter(chi Character) Character {
	out := make(Character, len(chi))
	for i := range chi {
		out[i] = cmplx.Conj(chi[i])
	}
	return out
}

// ClassValues reduces chi to one value per conjugacy class of g, taking the
// value on each class representative (the smallest element). It errors if chi
// is not a class function to within tol.
func ClassValues(g *Group, chi Character, tol float64) ([]complex128, error) {
	if !chi.IsClassFunction(g, tol) {
		return nil, errors.New("grouprep: value is not a class function")
	}
	classes := g.ConjugacyClasses()
	out := make([]complex128, len(classes))
	for i, cls := range classes {
		out[i] = chi[cls[0]]
	}
	return out, nil
}

// CharacterDimension returns the dimension of the representation affording chi,
// the rounded real part of chi on the identity.
func CharacterDimension(chi Character) int {
	return int(math.Round(real(chi[0])))
}

// RestrictCharacter restricts a character of a group G to a subgroup H given by
// an embedding emb, where emb[h] is the index in G of the element h of H. The
// result is the class function on H with value chi(emb[h]).
func RestrictCharacter(chi Character, emb []int) Character {
	out := make(Character, len(emb))
	for h, gi := range emb {
		out[h] = chi[gi]
	}
	return out
}

// InduceCharacter induces a character psi of a subgroup H up to the group G,
// using the embedding emb (emb[h] the index in G of the H-element h). It applies
// the Frobenius formula
//
//	Ind(psi)(y) = (1/|H|) Σ_{x∈G} psi°(x⁻¹ y x),
//
// where psi° extends psi by zero off H. It returns an error if emb is not a
// valid injection of a subgroup of G.
func InduceCharacter(gG *Group, gH *Group, emb []int, psi Character) (Character, error) {
	if len(emb) != gH.Order() || len(psi) != gH.Order() {
		return nil, errors.New("grouprep: InduceCharacter size mismatch")
	}
	gToH := make([]int, gG.Order())
	for i := range gToH {
		gToH[i] = -1
	}
	sub := make([]int, 0, len(emb))
	for h, gi := range emb {
		if gi < 0 || gi >= gG.Order() || gToH[gi] != -1 {
			return nil, errors.New("grouprep: embedding is not injective")
		}
		gToH[gi] = h
		sub = append(sub, gi)
	}
	if !gG.IsSubgroup(sub) {
		return nil, errors.New("grouprep: embedding image is not a subgroup")
	}
	nG := gG.Order()
	hOrd := complex(float64(gH.Order()), 0)
	out := make(Character, nG)
	for y := 0; y < nG; y++ {
		var s complex128
		for x := 0; x < nG; x++ {
			conj := gG.Mul(gG.Mul(gG.Inverse(x), y), x) // x⁻¹ y x
			if h := gToH[conj]; h != -1 {
				s += psi[h]
			}
		}
		out[y] = s / hOrd
	}
	return out, nil
}
