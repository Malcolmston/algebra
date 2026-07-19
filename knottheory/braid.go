package knottheory

import (
	"fmt"
	"strings"
)

// Braid is an element of the Artin braid group B_n on n strands, represented by
// a word in the standard generators. The word is a slice of non-zero integers:
// the value +i denotes the generator sigma_i (the i-th strand crossing over the
// (i+1)-th) and -i denotes its inverse sigma_i^{-1}, for 1 <= i <= n-1.
type Braid struct {
	strands int
	word    []int
}

// NewBraid returns the braid on the given number of strands whose word is the
// supplied sequence of generators. It returns an error if the strand count is
// less than one or if any generator index is zero or out of range.
func NewBraid(strands int, word ...int) (Braid, error) {
	if strands < 1 {
		return Braid{}, fmt.Errorf("knottheory: a braid needs at least one strand, got %d", strands)
	}
	for _, g := range word {
		if g == 0 {
			return Braid{}, fmt.Errorf("knottheory: 0 is not a valid braid generator")
		}
		a := g
		if a < 0 {
			a = -a
		}
		if a > strands-1 {
			return Braid{}, fmt.Errorf("knottheory: generator %d out of range for %d strands", g, strands)
		}
	}
	cp := make([]int, len(word))
	copy(cp, word)
	return Braid{strands: strands, word: cp}, nil
}

// MustBraid is like NewBraid but panics on error. It is convenient for
// constructing constant braids in tests and examples.
func MustBraid(strands int, word ...int) Braid {
	b, err := NewBraid(strands, word...)
	if err != nil {
		panic(err)
	}
	return b
}

// IdentityBraid returns the trivial (empty) braid on n strands.
func IdentityBraid(n int) Braid { return Braid{strands: n} }

// Generator returns the braid consisting of the single generator sigma_i on n
// strands.
func Generator(n, i int) (Braid, error) { return NewBraid(n, i) }

// GeneratorInv returns the braid consisting of the single inverse generator
// sigma_i^{-1} on n strands.
func GeneratorInv(n, i int) (Braid, error) { return NewBraid(n, -i) }

// IsValidBraidWord reports whether word is a legal word for B_n, that is every
// entry is non-zero and has absolute value at most n-1.
func IsValidBraidWord(strands int, word []int) bool {
	if strands < 1 {
		return false
	}
	for _, g := range word {
		if g == 0 {
			return false
		}
		a := g
		if a < 0 {
			a = -a
		}
		if a > strands-1 {
			return false
		}
	}
	return true
}

// BraidGroupRank returns the number of Artin generators of B_n, namely n-1.
func BraidGroupRank(n int) int {
	if n < 1 {
		return 0
	}
	return n - 1
}

// Strands returns the number of strands of the braid.
func (b Braid) Strands() int { return b.strands }

// Word returns a copy of the generator word of the braid.
func (b Braid) Word() []int {
	cp := make([]int, len(b.word))
	copy(cp, b.word)
	return cp
}

// Length returns the number of letters (generators) in the braid word, which
// equals the number of crossings in the associated braid diagram.
func (b Braid) Length() int { return len(b.word) }

// IsTrivial reports whether the braid word is empty.
func (b Braid) IsTrivial() bool { return len(b.word) == 0 }

// ExponentSum returns the algebraic sum of the generator signs, that is the
// number of positive generators minus the number of negative ones. It equals
// the writhe of the braid closure and is a braid invariant.
func (b Braid) ExponentSum() int {
	s := 0
	for _, g := range b.word {
		if g > 0 {
			s++
		} else {
			s--
		}
	}
	return s
}

// Writhe returns the writhe of the closure of the braid, which equals the
// exponent sum.
func (b Braid) Writhe() int { return b.ExponentSum() }

// ClosureCrossingNumber returns the number of crossings in the standard closed
// braid diagram, which is the braid word length.
func (b Braid) ClosureCrossingNumber() int { return len(b.word) }

// IsPositive reports whether every generator in the word is positive.
func (b Braid) IsPositive() bool {
	for _, g := range b.word {
		if g < 0 {
			return false
		}
	}
	return true
}

// IsNegative reports whether every generator in the word is negative.
func (b Braid) IsNegative() bool {
	for _, g := range b.word {
		if g > 0 {
			return false
		}
	}
	return len(b.word) > 0
}

// Inverse returns the inverse braid, whose word is the reversed word with every
// generator sign flipped.
func (b Braid) Inverse() Braid {
	w := make([]int, len(b.word))
	for i, g := range b.word {
		w[len(b.word)-1-i] = -g
	}
	return Braid{strands: b.strands, word: w}
}

// Reverse returns the braid whose word is read right to left (the reverse or
// "review" of the braid), keeping every generator sign.
func (b Braid) Reverse() Braid {
	w := make([]int, len(b.word))
	for i, g := range b.word {
		w[len(b.word)-1-i] = g
	}
	return Braid{strands: b.strands, word: w}
}

// Mirror returns the braid obtained by flipping every crossing, that is by
// negating every generator. Its closure is the mirror image of the original
// link.
func (b Braid) Mirror() Braid {
	w := make([]int, len(b.word))
	for i, g := range b.word {
		w[i] = -g
	}
	return Braid{strands: b.strands, word: w}
}

// Concat returns the product b*other formed by concatenating the two braid
// words. It returns an error if the strand counts differ.
func (b Braid) Concat(other Braid) (Braid, error) {
	if b.strands != other.strands {
		return Braid{}, fmt.Errorf("knottheory: cannot concatenate braids on %d and %d strands", b.strands, other.strands)
	}
	w := make([]int, 0, len(b.word)+len(other.word))
	w = append(w, b.word...)
	w = append(w, other.word...)
	return Braid{strands: b.strands, word: w}, nil
}

// Power returns the braid raised to the integer power k. Negative powers use the
// inverse braid and k==0 gives the identity braid.
func (b Braid) Power(k int) Braid {
	base := b
	if k < 0 {
		base = b.Inverse()
		k = -k
	}
	res := IdentityBraid(b.strands)
	for i := 0; i < k; i++ {
		res, _ = res.Concat(base)
	}
	return res
}

// FreeReduce returns the braid obtained by repeatedly cancelling adjacent
// inverse pairs sigma_i sigma_i^{-1}. The result is freely reduced but is not a
// solution of the braid word problem; distinct freely reduced words can still
// represent the same braid.
func (b Braid) FreeReduce() Braid {
	stack := make([]int, 0, len(b.word))
	for _, g := range b.word {
		if n := len(stack); n > 0 && stack[n-1] == -g {
			stack = stack[:n-1]
		} else {
			stack = append(stack, g)
		}
	}
	w := make([]int, len(stack))
	copy(w, stack)
	return Braid{strands: b.strands, word: w}
}

// IsFreelyReduced reports whether the braid word contains no adjacent inverse
// pair.
func (b Braid) IsFreelyReduced() bool {
	for i := 0; i+1 < len(b.word); i++ {
		if b.word[i] == -b.word[i+1] {
			return false
		}
	}
	return true
}

// Permutation returns the underlying permutation of the braid: strand starting
// at position i ends at position Permutation()[i]. Generator sigma_i acts as the
// transposition (i, i+1) on positions numbered from 0.
func (b Braid) Permutation() Permutation {
	p := IdentityPermutation(b.strands)
	// Reading the braid from top to bottom, each generator swaps the two
	// strands currently occupying positions i and i+1.
	for _, g := range b.word {
		i := g
		if i < 0 {
			i = -i
		}
		i-- // to zero-based lower position
		p[i], p[i+1] = p[i+1], p[i]
	}
	return p
}

// NumComponents returns the number of link components of the braid closure,
// which equals the number of cycles of the underlying permutation.
func (b Braid) NumComponents() int { return b.Permutation().NumCycles() }

// ClosureIsKnot reports whether the closure of the braid is a knot (a single
// component).
func (b Braid) ClosureIsKnot() bool { return b.NumComponents() == 1 }

// String renders the braid word using sigma notation, for example
// "s1 s2^-1 s1".
func (b Braid) String() string {
	if len(b.word) == 0 {
		return fmt.Sprintf("1 (B_%d)", b.strands)
	}
	parts := make([]string, len(b.word))
	for i, g := range b.word {
		if g > 0 {
			parts[i] = fmt.Sprintf("s%d", g)
		} else {
			parts[i] = fmt.Sprintf("s%d^-1", -g)
		}
	}
	return strings.Join(parts, " ")
}

// TorusBraid returns the braid (sigma_1 sigma_2 ... sigma_{p-1})^q on p strands
// whose closure is the torus link T(p, q). It returns an error if p < 1.
func TorusBraid(p, q int) (Braid, error) {
	if p < 1 {
		return Braid{}, fmt.Errorf("knottheory: torus braid needs p >= 1, got %d", p)
	}
	base := make([]int, 0, p-1)
	for i := 1; i <= p-1; i++ {
		base = append(base, i)
	}
	word := make([]int, 0, (p-1)*abs(q))
	rep := q
	neg := false
	if rep < 0 {
		neg = true
		rep = -rep
	}
	for r := 0; r < rep; r++ {
		if neg {
			for i := len(base) - 1; i >= 0; i-- {
				word = append(word, -base[i])
			}
		} else {
			word = append(word, base...)
		}
	}
	return NewBraid(p, word...)
}

// FullTwist returns the full twist braid (sigma_1 ... sigma_{n-1})^n on n
// strands, the generator of the centre of B_n.
func FullTwist(n int) Braid {
	b, _ := TorusBraid(n, n)
	return b
}

// GarsideHalfTwist returns the Garside half-twist braid Delta_n on n strands,
// the positive half twist whose square is the full twist. It is the product
// (sigma_1)(sigma_2 sigma_1)(sigma_3 sigma_2 sigma_1)...
func GarsideHalfTwist(n int) Braid {
	word := make([]int, 0, n*(n-1)/2)
	for i := 1; i <= n-1; i++ {
		for j := i; j >= 1; j-- {
			word = append(word, j)
		}
	}
	b, _ := NewBraid(n, word...)
	return b
}

// abs returns the absolute value of x.
func abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

// reducedBurauGen returns the (n-1)x(n-1) reduced Burau matrix of the generator
// sigma_g (g>0) or its inverse (g<0) on n strands. Entries are Laurent
// polynomials in the variable t.
func reducedBurauGen(n, g int) *LaurentMatrix {
	m := n - 1
	M := IdentityLaurentMatrix(m)
	i := g
	inv := false
	if i < 0 {
		inv = true
		i = -i
	}
	negT := Monomial(-1, 1)     // -t
	negTinv := Monomial(-1, -1) // -t^{-1}
	t := Monomial(1, 1)         // t
	tinv := Monomial(1, -1)     // t^{-1}
	one := OneLaurent()         //  1
	set := func(r, c int, v Laurent) { M.Set(r, c, v) }
	switch {
	case m == 1:
		// B_2: the reduced Burau representation is one-dimensional and
		// sigma_1 acts as multiplication by -t.
		if !inv {
			set(0, 0, negT)
		} else {
			set(0, 0, negTinv)
		}
	case i == 1:
		// block [[-t,0],[1,1]] on coordinates 0,1 (1-based 1,2)
		if !inv {
			set(0, 0, negT)
			set(1, 0, one)
			set(1, 1, one)
			set(0, 1, ZeroLaurent())
		} else {
			set(0, 0, negTinv)
			set(1, 0, tinv)
			set(1, 1, one)
			set(0, 1, ZeroLaurent())
		}
	case i == n-1:
		// block on coordinates m-2, m-1 (1-based n-2, n-1)
		a := m - 2
		b := m - 1
		if !inv {
			set(a, a, one)
			set(a, b, t)
			set(b, a, ZeroLaurent())
			set(b, b, negT)
		} else {
			set(a, a, one)
			set(a, b, one)
			set(b, a, ZeroLaurent())
			set(b, b, negTinv)
		}
	default:
		// middle generator: 3x3 block on coordinates i-2, i-1, i (1-based i-1,i,i+1)
		a := i - 2
		c := i - 1
		d := i
		if !inv {
			set(a, a, one)
			set(a, c, t)
			set(a, d, ZeroLaurent())
			set(c, a, ZeroLaurent())
			set(c, c, negT)
			set(c, d, ZeroLaurent())
			set(d, a, ZeroLaurent())
			set(d, c, one)
			set(d, d, one)
		} else {
			set(a, a, one)
			set(a, c, one)
			set(a, d, ZeroLaurent())
			set(c, a, ZeroLaurent())
			set(c, c, negTinv)
			set(c, d, ZeroLaurent())
			set(d, a, ZeroLaurent())
			set(d, c, tinv)
			set(d, d, one)
		}
	}
	return M
}

// ReducedBurau returns the reduced Burau matrix of the braid, an (n-1)x(n-1)
// matrix of Laurent polynomials in t obtained as the ordered product of the
// generator matrices. For the identity braid it is the identity matrix.
func (b Braid) ReducedBurau() *LaurentMatrix {
	m := b.strands - 1
	if m <= 0 {
		return NewLaurentMatrix(0, 0)
	}
	acc := IdentityLaurentMatrix(m)
	for _, g := range b.word {
		acc = acc.Mul(reducedBurauGen(b.strands, g))
	}
	return acc
}

// AlexanderPolynomial returns the Alexander polynomial of the braid closure,
// computed from the reduced Burau matrix via the Burau formula
//
//	det(ReducedBurau - I) = ± t^k (1 + t + ... + t^{n-1}) / (1 - t) * Delta(t),
//
// implemented as det(ReducedBurau - I) divided by 1 + t + ... + t^{n-1} and then
// normalised to the symmetric representative with Delta(1) = 1. It returns an
// error if the closure is not a knot or if the division is not exact.
func (b Braid) AlexanderPolynomial() (Laurent, error) {
	if !b.ClosureIsKnot() {
		return Laurent{}, fmt.Errorf("knottheory: Burau Alexander polynomial requires a knot (single component), closure has %d components", b.NumComponents())
	}
	if b.strands == 1 {
		// The unknot.
		return OneLaurent(), nil
	}
	m := b.strands - 1
	red := b.ReducedBurau()
	diff := red.Sub(IdentityLaurentMatrix(m))
	det := diff.Determinant()
	// Divide by 1 + t + ... + t^{n-1}.
	denomCoeffs := make([]int, b.strands)
	for i := range denomCoeffs {
		denomCoeffs[i] = 1
	}
	denom := NewLaurent(0, denomCoeffs)
	q, ok := det.DivExact(denom)
	if !ok {
		return Laurent{}, fmt.Errorf("knottheory: Burau determinant not divisible by the cyclotomic denominator")
	}
	return normalizeAlexander(q), nil
}

// normalizeAlexander returns the canonical representative of an Alexander
// polynomial: it is shifted to be symmetric about exponent zero when possible
// and its overall sign is chosen so that the value at t=1 is positive.
func normalizeAlexander(L Laurent) Laurent {
	if L.IsZero() {
		return L
	}
	lo := L.MinDegree()
	hi := L.MaxDegree()
	shift := -(lo + hi) / 2
	L = L.ShiftExp(shift)
	if L.EvalUnit(1) < 0 {
		L = L.Neg()
	}
	return L
}
