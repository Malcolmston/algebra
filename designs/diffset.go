package designs

import (
	"errors"
	"sort"
)

// DifferenceSet is a subset D of the cyclic group Z_n. It is a (n,k,lambda)
// difference set when |D|=k and every non-zero element of Z_n arises as a
// difference d-d' (mod n) of ordered pairs from D in exactly lambda ways.
type DifferenceSet struct {
	N        int
	Elements []int
}

// NewDifferenceSet builds a subset of Z_n from the given residues, reducing
// them modulo n, removing duplicates and sorting. It reports an error when n<=0.
func NewDifferenceSet(n int, elems []int) (*DifferenceSet, error) {
	if n <= 0 {
		return nil, errors.New("designs: group order must be positive")
	}
	seen := make(map[int]bool)
	var out []int
	for _, e := range elems {
		e = ((e % n) + n) % n
		if !seen[e] {
			seen[e] = true
			out = append(out, e)
		}
	}
	sort.Ints(out)
	return &DifferenceSet{N: n, Elements: out}, nil
}

// Order returns the order n of the ambient cyclic group.
func (d *DifferenceSet) Order() int { return d.N }

// Size returns the number of elements k of the set.
func (d *DifferenceSet) Size() int { return len(d.Elements) }

// DifferenceCounts returns a slice indexed by g in [0,n) giving the number of
// ordered pairs (x,y) with x,y in D and x-y == g (mod n). Index 0 always holds
// k (the trivial differences).
func (d *DifferenceSet) DifferenceCounts() []int {
	counts := make([]int, d.N)
	for _, x := range d.Elements {
		for _, y := range d.Elements {
			g := ((x-y)%d.N + d.N) % d.N
			counts[g]++
		}
	}
	return counts
}

// IsDifferenceSet reports whether D is a (n,k,lambda) difference set: every
// non-zero group element occurs equally often as a difference. A set of size 0,
// 1 or n does not qualify.
func (d *DifferenceSet) IsDifferenceSet() bool {
	_, ok := d.Lambda()
	return ok
}

// Lambda returns the common difference multiplicity lambda together with true
// when D is a difference set, or (0,false) otherwise.
func (d *DifferenceSet) Lambda() (int, bool) {
	k := len(d.Elements)
	if k < 2 || k >= d.N {
		return 0, false
	}
	counts := d.DifferenceCounts()
	lambda := counts[1]
	for g := 1; g < d.N; g++ {
		if counts[g] != lambda {
			return 0, false
		}
	}
	if lambda < 1 {
		return 0, false
	}
	return lambda, true
}

// Parameters returns the (n,k,lambda) parameters of the difference set. It
// reports an error when D is not a difference set.
func (d *DifferenceSet) Parameters() (n, k, lambda int, err error) {
	l, ok := d.Lambda()
	if !ok {
		return 0, 0, 0, errors.New("designs: not a difference set")
	}
	return d.N, len(d.Elements), l, nil
}

// IsPlanar reports whether D is a planar (lambda=1) difference set, which
// develops into a projective plane.
func (d *DifferenceSet) IsPlanar() bool {
	l, ok := d.Lambda()
	return ok && l == 1
}

// Translate returns the shifted set D+s (mod n).
func (d *DifferenceSet) Translate(s int) *DifferenceSet {
	elems := make([]int, len(d.Elements))
	for i, x := range d.Elements {
		elems[i] = ((x+s)%d.N + d.N) % d.N
	}
	nd, _ := NewDifferenceSet(d.N, elems)
	return nd
}

// Multiply returns the dilated set t*D (mod n).
func (d *DifferenceSet) Multiply(t int) *DifferenceSet {
	elems := make([]int, len(d.Elements))
	for i, x := range d.Elements {
		elems[i] = ((x*t)%d.N + d.N) % d.N
	}
	nd, _ := NewDifferenceSet(d.N, elems)
	return nd
}

// Complement returns the complement of D within Z_n.
func (d *DifferenceSet) Complement() *DifferenceSet {
	in := make([]bool, d.N)
	for _, x := range d.Elements {
		in[x] = true
	}
	var elems []int
	for x := 0; x < d.N; x++ {
		if !in[x] {
			elems = append(elems, x)
		}
	}
	nd, _ := NewDifferenceSet(d.N, elems)
	return nd
}

// setEqual reports whether the sorted element lists of two subsets of the same
// Z_n coincide.
func (d *DifferenceSet) setEqual(o *DifferenceSet) bool {
	if len(d.Elements) != len(o.Elements) {
		return false
	}
	for i := range d.Elements {
		if d.Elements[i] != o.Elements[i] {
			return false
		}
	}
	return true
}

// IsMultiplier reports whether t is a numerical multiplier of D, meaning t is
// coprime to n and t*D is a translate D+s of D for some s.
func (d *DifferenceSet) IsMultiplier(t int) bool {
	if Gcd(((t%d.N)+d.N)%d.N, d.N) != 1 {
		return false
	}
	td := d.Multiply(t)
	for s := 0; s < d.N; s++ {
		if td.setEqual(d.Translate(s)) {
			return true
		}
	}
	return false
}

// Multipliers returns all numerical multipliers t in [1,n) of the difference
// set, i.e. the units t of Z_n for which t*D is a translate of D.
func (d *DifferenceSet) Multipliers() []int {
	var out []int
	for t := 1; t < d.N; t++ {
		if d.IsMultiplier(t) {
			out = append(out, t)
		}
	}
	return out
}

// Develop returns the symmetric design obtained by developing D: the point set
// is Z_n and the blocks are the n translates D+g for g in Z_n. When D is a
// (n,k,lambda) difference set this is a symmetric 2-(n,k,lambda) design.
func (d *DifferenceSet) Develop() *Design {
	blocks := make([][]int, d.N)
	for g := 0; g < d.N; g++ {
		b := make([]int, len(d.Elements))
		for i, x := range d.Elements {
			b[i] = (x + g) % d.N
		}
		blocks[g] = b
	}
	nd, _ := NewDesign(d.N, blocks)
	return nd
}

// IsDifferenceSetSlice reports whether the given residues form a (n,k,lambda)
// difference set in Z_n.
func IsDifferenceSetSlice(n int, elems []int) bool {
	d, err := NewDifferenceSet(n, elems)
	if err != nil {
		return false
	}
	return d.IsDifferenceSet()
}

// QuadraticResidueDifferenceSet returns the Paley difference set of quadratic
// residues in Z_p for a prime p congruent to 3 modulo 4, a
// (p, (p-1)/2, (p-3)/4) difference set. It reports an error when p is not such
// a prime.
func QuadraticResidueDifferenceSet(p int) (*DifferenceSet, error) {
	if !IsPrime(p) || p%4 != 3 {
		return nil, errors.New("designs: require a prime p = 3 (mod 4)")
	}
	return NewDifferenceSet(p, QuadraticResidues(p))
}

// SingerDifferenceSet returns the Singer planar difference set of order q, a
// (q^2+q+1, q+1, 1) difference set constructed from the trace-zero hyperplane
// of GF(q^3) over GF(q). Developing it yields the projective plane PG(2,q). It
// reports an error when q is not a prime power.
func SingerDifferenceSet(q int) (*DifferenceSet, error) {
	if !IsPrimePower(q) {
		return nil, errors.New("designs: q must be a prime power")
	}
	f, err := NewGaloisField(q * q * q)
	if err != nil {
		return nil, err
	}
	g, err := f.PrimitiveElement()
	if err != nil {
		return nil, err
	}
	n := q*q + q + 1
	trace := func(x int) int {
		return f.Add(f.Add(x, f.Pow(x, q)), f.Pow(x, q*q))
	}
	var elems []int
	cur := 1 // g^0
	for i := 0; i < n; i++ {
		if trace(cur) == 0 {
			elems = append(elems, i)
		}
		cur = f.Mul(cur, g)
	}
	return NewDifferenceSet(n, elems)
}
