package grouprep

import (
	"errors"
	"fmt"
	"sort"
)

// Group is a finite group given by an explicit Cayley (multiplication) table.
// Elements are the integers 0,...,n-1; element 0 is always the identity. The
// table satisfies table[i][j] = index of the product i·j. Construct groups with
// the family constructors ([CyclicGroup], [DihedralGroup], ...) rather than by
// hand.
type Group struct {
	name  string
	names []string
	table [][]int
	inv   []int
}

// newGroupFromTable validates a Cayley table and builds a Group. The identity
// is required to be element 0. It returns an error if the data do not form a
// group.
func newGroupFromTable(name string, names []string, table [][]int) (*Group, error) {
	n := len(table)
	if n == 0 {
		return nil, errors.New("grouprep: empty group")
	}
	for i := range table {
		if len(table[i]) != n {
			return nil, errors.New("grouprep: non-square Cayley table")
		}
		for _, v := range table[i] {
			if v < 0 || v >= n {
				return nil, errors.New("grouprep: Cayley entry out of range")
			}
		}
	}
	// Identity check: element 0 must be a two-sided identity.
	for i := 0; i < n; i++ {
		if table[0][i] != i || table[i][0] != i {
			return nil, errors.New("grouprep: element 0 is not the identity")
		}
	}
	// Inverses.
	inv := make([]int, n)
	for i := range inv {
		inv[i] = -1
	}
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if table[i][j] == 0 && table[j][i] == 0 {
				inv[i] = j
			}
		}
		if inv[i] == -1 {
			return nil, fmt.Errorf("grouprep: element %d has no inverse", i)
		}
	}
	// Latin-square property (each row and column a permutation).
	for i := 0; i < n; i++ {
		seenR := make([]bool, n)
		seenC := make([]bool, n)
		for j := 0; j < n; j++ {
			if seenR[table[i][j]] || seenC[table[j][i]] {
				return nil, errors.New("grouprep: Cayley table is not a Latin square")
			}
			seenR[table[i][j]] = true
			seenC[table[j][i]] = true
		}
	}
	g := &Group{name: name, names: names, table: table, inv: inv}
	return g, nil
}

// Order returns the number of elements of g.
func (g *Group) Order() int { return len(g.table) }

// Name returns the display name of g.
func (g *Group) Name() string { return g.name }

// Identity returns the index of the identity element, which is always 0.
func (g *Group) Identity() int { return 0 }

// ElementName returns a human-readable label for element i.
func (g *Group) ElementName(i int) string {
	if g.names != nil && i >= 0 && i < len(g.names) {
		return g.names[i]
	}
	return fmt.Sprintf("g%d", i)
}

// Elements returns the slice {0,1,...,Order()-1}.
func (g *Group) Elements() []int {
	out := make([]int, g.Order())
	for i := range out {
		out[i] = i
	}
	return out
}

// Mul returns the product i·j.
func (g *Group) Mul(i, j int) int { return g.table[i][j] }

// Inverse returns the inverse of element i.
func (g *Group) Inverse(i int) int { return g.inv[i] }

// Conjugate returns the conjugate g x g⁻¹ (with a the conjugator and x the
// element), i.e. a·x·a⁻¹.
func (g *Group) Conjugate(a, x int) int {
	return g.Mul(g.Mul(a, x), g.inv[a])
}

// Commutator returns the group commutator i j i⁻¹ j⁻¹.
func (g *Group) Commutator(i, j int) int {
	return g.Mul(g.Mul(i, j), g.Mul(g.inv[i], g.inv[j]))
}

// Pow returns element i raised to the integer power k (negative k uses the
// inverse).
func (g *Group) Pow(i, k int) int {
	if k < 0 {
		i = g.inv[i]
		k = -k
	}
	result := 0 // identity
	for k > 0 {
		if k&1 == 1 {
			result = g.Mul(result, i)
		}
		i = g.Mul(i, i)
		k >>= 1
	}
	return result
}

// ElementOrder returns the order of element i, the least positive k with
// i^k = identity.
func (g *Group) ElementOrder(i int) int {
	order := 1
	x := i
	for x != 0 {
		x = g.Mul(x, i)
		order++
	}
	return order
}

// Exponent returns the exponent of g, the least common multiple of all element
// orders.
func (g *Group) Exponent() int {
	e := 1
	for i := 0; i < g.Order(); i++ {
		e = lcmInt(e, g.ElementOrder(i))
	}
	return e
}

// CayleyTable returns a copy of the multiplication table of g.
func (g *Group) CayleyTable() [][]int {
	out := make([][]int, g.Order())
	for i := range g.table {
		out[i] = append([]int(nil), g.table[i]...)
	}
	return out
}

// IsAbelian reports whether g is commutative.
func (g *Group) IsAbelian() bool {
	n := g.Order()
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if g.Mul(i, j) != g.Mul(j, i) {
				return false
			}
		}
	}
	return true
}

// Center returns the sorted indices of the elements commuting with every
// element of g.
func (g *Group) Center() []int {
	n := g.Order()
	var out []int
	for i := 0; i < n; i++ {
		central := true
		for j := 0; j < n; j++ {
			if g.Mul(i, j) != g.Mul(j, i) {
				central = false
				break
			}
		}
		if central {
			out = append(out, i)
		}
	}
	return out
}

// Centralizer returns the sorted indices of the elements commuting with x.
func (g *Group) Centralizer(x int) []int {
	var out []int
	for i := 0; i < g.Order(); i++ {
		if g.Mul(i, x) == g.Mul(x, i) {
			out = append(out, i)
		}
	}
	return out
}

// ConjugacyClasses returns the conjugacy classes of g, each a sorted slice of
// element indices. The classes are ordered by their smallest element, so the
// identity class {0} comes first.
func (g *Group) ConjugacyClasses() [][]int {
	n := g.Order()
	assigned := make([]bool, n)
	var classes [][]int
	for x := 0; x < n; x++ {
		if assigned[x] {
			continue
		}
		set := map[int]bool{}
		for a := 0; a < n; a++ {
			set[g.Conjugate(a, x)] = true
		}
		cls := make([]int, 0, len(set))
		for e := range set {
			cls = append(cls, e)
			assigned[e] = true
		}
		sort.Ints(cls)
		classes = append(classes, cls)
	}
	sort.Slice(classes, func(i, j int) bool { return classes[i][0] < classes[j][0] })
	return classes
}

// NumConjugacyClasses returns the number of conjugacy classes of g, which
// equals the number of complex irreducible representations.
func (g *Group) NumConjugacyClasses() int {
	return len(g.ConjugacyClasses())
}

// ClassOf returns the conjugacy class containing element x.
func (g *Group) ClassOf(x int) []int {
	set := map[int]bool{}
	for a := 0; a < g.Order(); a++ {
		set[g.Conjugate(a, x)] = true
	}
	out := make([]int, 0, len(set))
	for e := range set {
		out = append(out, e)
	}
	sort.Ints(out)
	return out
}

// GeneratedBy returns the sorted set of elements of the subgroup generated by
// gens (the closure under multiplication and inverses).
func (g *Group) GeneratedBy(gens []int) []int {
	inSet := map[int]bool{0: true}
	frontier := []int{0}
	for _, gen := range gens {
		if !inSet[gen] {
			inSet[gen] = true
			frontier = append(frontier, gen)
		}
	}
	for len(frontier) > 0 {
		x := frontier[len(frontier)-1]
		frontier = frontier[:len(frontier)-1]
		for _, gen := range gens {
			for _, y := range []int{g.Mul(x, gen), g.Mul(gen, x)} {
				if !inSet[y] {
					inSet[y] = true
					frontier = append(frontier, y)
				}
			}
		}
	}
	out := make([]int, 0, len(inSet))
	for e := range inSet {
		out = append(out, e)
	}
	sort.Ints(out)
	return out
}

// IsSubgroup reports whether the given set of element indices is closed under
// the group operation and inverses (and hence forms a subgroup).
func (g *Group) IsSubgroup(elems []int) bool {
	set := map[int]bool{}
	for _, e := range elems {
		if e < 0 || e >= g.Order() {
			return false
		}
		set[e] = true
	}
	if !set[0] {
		return false
	}
	for a := range set {
		if !set[g.inv[a]] {
			return false
		}
		for b := range set {
			if !set[g.Mul(a, b)] {
				return false
			}
		}
	}
	return true
}

// IsValid re-checks that g really is a group. It returns nil for the groups the
// constructors produce.
func (g *Group) IsValid() error {
	_, err := newGroupFromTable(g.name, g.names, g.table)
	return err
}

// TrivialGroup returns the group with one element.
func TrivialGroup() *Group {
	g, _ := newGroupFromTable("C1", []string{"e"}, [][]int{{0}})
	return g
}

// CyclicGroup returns the cyclic group Z/nZ of order n written additively:
// element i is i mod n and i·j = (i+j) mod n. It panics if n < 1.
func CyclicGroup(n int) *Group {
	if n < 1 {
		panic("grouprep: CyclicGroup requires n >= 1")
	}
	table := make([][]int, n)
	names := make([]string, n)
	for i := 0; i < n; i++ {
		table[i] = make([]int, n)
		for j := 0; j < n; j++ {
			table[i][j] = (i + j) % n
		}
		if i == 0 {
			names[i] = "e"
		} else {
			names[i] = fmt.Sprintf("r^%d", i)
		}
	}
	g, _ := newGroupFromTable(fmt.Sprintf("C%d", n), names, table)
	return g
}

// KleinFourGroup returns the Klein four-group V = C2×C2 of order 4.
func KleinFourGroup() *Group {
	return DirectProduct(CyclicGroup(2), CyclicGroup(2))
}

// DihedralGroup returns the dihedral group D_n of order 2n, the symmetry group
// of a regular n-gon. Element a·n+b denotes s^a r^b with r a rotation of order
// n, s a reflection, and s r s = r⁻¹. It panics if n < 1.
func DihedralGroup(n int) *Group {
	if n < 1 {
		panic("grouprep: DihedralGroup requires n >= 1")
	}
	size := 2 * n
	idx := func(a, b int) int { return a*n + ((b%n)+n)%n }
	table := make([][]int, size)
	names := make([]string, size)
	for a := 0; a < 2; a++ {
		for b := 0; b < n; b++ {
			i := idx(a, b)
			table[i] = make([]int, size)
			for a2 := 0; a2 < 2; a2++ {
				for b2 := 0; b2 < n; b2++ {
					j := idx(a2, b2)
					na := (a + a2) % 2
					var nb int
					if a2 == 0 {
						nb = b + b2
					} else {
						nb = -b + b2
					}
					table[i][j] = idx(na, nb)
				}
			}
			switch {
			case a == 0 && b == 0:
				names[i] = "e"
			case a == 0:
				names[i] = fmt.Sprintf("r^%d", b)
			case b == 0:
				names[i] = "s"
			default:
				names[i] = fmt.Sprintf("s*r^%d", b)
			}
		}
	}
	g, err := newGroupFromTable(fmt.Sprintf("D%d", n), names, table)
	if err != nil {
		panic(err)
	}
	return g
}

// symmetricData builds the symmetric group S_n and returns it together with the
// list of permutations in the same element order (lexicographic). Element 0 is
// the identity permutation.
func symmetricData(n int) (*Group, []Perm) {
	if n < 1 {
		panic("grouprep: SymmetricGroup requires n >= 1")
	}
	perms := AllPermutations(n)
	// Ensure the identity is first.
	sort.Slice(perms, func(i, j int) bool {
		for k := 0; k < n; k++ {
			if perms[i][k] != perms[j][k] {
				return perms[i][k] < perms[j][k]
			}
		}
		return false
	})
	index := map[string]int{}
	key := func(p Perm) string {
		b := make([]byte, len(p))
		for i, v := range p {
			b[i] = byte(v)
		}
		return string(b)
	}
	for i, p := range perms {
		index[key(p)] = i
	}
	size := len(perms)
	table := make([][]int, size)
	names := make([]string, size)
	for i := range perms {
		table[i] = make([]int, size)
		for j := range perms {
			table[i][j] = index[key(perms[i].Compose(perms[j]))]
		}
		if i == 0 {
			names[i] = "e"
		} else {
			names[i] = permName(perms[i])
		}
	}
	g, err := newGroupFromTable(fmt.Sprintf("S%d", n), names, table)
	if err != nil {
		panic(err)
	}
	return g, perms
}

// SymmetricGroup returns the symmetric group S_n of all permutations of n
// points, of order n!. It panics if n < 1. Practical use is limited to small n.
func SymmetricGroup(n int) *Group {
	g, _ := symmetricData(n)
	return g
}

// AlternatingGroup returns the alternating group A_n of even permutations of n
// points, of order n!/2 for n >= 2. It panics if n < 1.
func AlternatingGroup(n int) *Group {
	if n < 1 {
		panic("grouprep: AlternatingGroup requires n >= 1")
	}
	all := AllPermutations(n)
	var even []Perm
	for _, p := range all {
		if p.Sign() == 1 {
			even = append(even, p)
		}
	}
	sort.Slice(even, func(i, j int) bool {
		for k := 0; k < n; k++ {
			if even[i][k] != even[j][k] {
				return even[i][k] < even[j][k]
			}
		}
		return false
	})
	index := map[string]int{}
	key := func(p Perm) string {
		b := make([]byte, len(p))
		for i, v := range p {
			b[i] = byte(v)
		}
		return string(b)
	}
	for i, p := range even {
		index[key(p)] = i
	}
	size := len(even)
	table := make([][]int, size)
	names := make([]string, size)
	for i := range even {
		table[i] = make([]int, size)
		for j := range even {
			table[i][j] = index[key(even[i].Compose(even[j]))]
		}
		if i == 0 {
			names[i] = "e"
		} else {
			names[i] = permName(even[i])
		}
	}
	g, err := newGroupFromTable(fmt.Sprintf("A%d", n), names, table)
	if err != nil {
		panic(err)
	}
	return g
}

// QuaternionGroup returns the quaternion group Q8 = {±1, ±i, ±j, ±k} of order
// 8. Elements are indexed 0:1, 1:-1, 2:i, 3:-i, 4:j, 5:-j, 6:k, 7:-k.
func QuaternionGroup() *Group {
	// units: 0=1, 1=i, 2=j, 3=k; element index = unit*2 + (sign<0 ? 1 : 0).
	unitMul := func(u1, u2 int) (int, int) { // returns (signFactor, resultUnit)
		if u1 == 0 {
			return 1, u2
		}
		if u2 == 0 {
			return 1, u1
		}
		if u1 == u2 {
			return -1, 0
		}
		// i,j,k cyclic: i*j=k, j*k=i, k*i=j (positive)
		switch [2]int{u1, u2} {
		case [2]int{1, 2}:
			return 1, 3
		case [2]int{2, 3}:
			return 1, 1
		case [2]int{3, 1}:
			return 1, 2
		case [2]int{2, 1}:
			return -1, 3
		case [2]int{3, 2}:
			return -1, 1
		case [2]int{1, 3}:
			return -1, 2
		}
		panic("grouprep: bad quaternion units")
	}
	decode := func(idx int) (unit, sign int) {
		unit = idx / 2
		if idx%2 == 0 {
			sign = 1
		} else {
			sign = -1
		}
		return
	}
	encode := func(unit, sign int) int {
		if sign < 0 {
			return unit*2 + 1
		}
		return unit * 2
	}
	names := []string{"1", "-1", "i", "-i", "j", "-j", "k", "-k"}
	table := make([][]int, 8)
	for i := 0; i < 8; i++ {
		table[i] = make([]int, 8)
		u1, s1 := decode(i)
		for j := 0; j < 8; j++ {
			u2, s2 := decode(j)
			sf, ru := unitMul(u1, u2)
			table[i][j] = encode(ru, s1*s2*sf)
		}
	}
	g, err := newGroupFromTable("Q8", names, table)
	if err != nil {
		panic(err)
	}
	return g
}

// DirectProduct returns the external direct product g×h. Element a·|h|+b
// denotes the pair (a, b) with componentwise multiplication.
func DirectProduct(g, h *Group) *Group {
	ng, nh := g.Order(), h.Order()
	size := ng * nh
	idx := func(a, b int) int { return a*nh + b }
	table := make([][]int, size)
	names := make([]string, size)
	for a := 0; a < ng; a++ {
		for b := 0; b < nh; b++ {
			i := idx(a, b)
			table[i] = make([]int, size)
			for a2 := 0; a2 < ng; a2++ {
				for b2 := 0; b2 < nh; b2++ {
					table[i][idx(a2, b2)] = idx(g.Mul(a, a2), h.Mul(b, b2))
				}
			}
			names[i] = fmt.Sprintf("(%s,%s)", g.ElementName(a), h.ElementName(b))
		}
	}
	res, err := newGroupFromTable(fmt.Sprintf("%s×%s", g.name, h.name), names, table)
	if err != nil {
		panic(err)
	}
	return res
}

// permName renders a non-identity permutation in cycle notation.
func permName(p Perm) string {
	s := ""
	for _, c := range p.CycleDecomposition() {
		if len(c) == 1 {
			continue
		}
		s += "("
		for i, v := range c {
			if i > 0 {
				s += " "
			}
			s += fmt.Sprintf("%d", v)
		}
		s += ")"
	}
	if s == "" {
		return "e"
	}
	return s
}
