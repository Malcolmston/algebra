package grouprep

import (
	"fmt"
	"math"
	"math/cmplx"
	"sort"
	"strings"
)

// CharacterTable holds a complete set of irreducible characters of a group,
// tabulated against its conjugacy classes. The number of irreducibles equals
// the number of conjugacy classes. Rows are irreducible characters and columns
// are classes.
type CharacterTable struct {
	group      *Group
	classes    [][]int
	classSizes []int
	reps       []int
	orders     []int
	irr        []Character    // full class functions, one per irreducible
	values     [][]complex128 // values[i][c] = irr i on class c
	dims       []int
}

// BuildCharacterTable assembles a character table for g from a list of
// irreducible characters (each a full class function over all |G| elements).
// The characters are sorted by degree (ascending), then lexicographically, so
// the trivial character comes first. The rows are reduced to one value per
// conjugacy class.
func BuildCharacterTable(g *Group, irrs []Character) *CharacterTable {
	classes := g.ConjugacyClasses()
	sizes := make([]int, len(classes))
	reps := make([]int, len(classes))
	orders := make([]int, len(classes))
	for i, cls := range classes {
		sizes[i] = len(cls)
		reps[i] = cls[0]
		orders[i] = g.ElementOrder(cls[0])
	}
	rows := make([]Character, len(irrs))
	copy(rows, irrs)
	sort.SliceStable(rows, func(a, b int) bool {
		da, db := real(rows[a][0]), real(rows[b][0])
		if math.Abs(da-db) > 1e-9 {
			return da < db
		}
		for c := 0; c < len(rows[a]); c++ {
			ra, rb := real(rows[a][c]), real(rows[b][c])
			if math.Abs(ra-rb) > 1e-9 {
				return ra > rb
			}
		}
		return false
	})
	values := make([][]complex128, len(rows))
	dims := make([]int, len(rows))
	for i, chi := range rows {
		values[i] = make([]complex128, len(classes))
		for c, cls := range classes {
			values[i][c] = chi[cls[0]]
		}
		dims[i] = int(math.Round(real(chi[0])))
	}
	return &CharacterTable{
		group: g, classes: classes, classSizes: sizes, reps: reps,
		orders: orders, irr: rows, values: values, dims: dims,
	}
}

// Group returns the group of the table.
func (ct *CharacterTable) Group() *Group { return ct.group }

// NumIrreducibles returns the number of irreducible characters (rows).
func (ct *CharacterTable) NumIrreducibles() int { return len(ct.irr) }

// NumClasses returns the number of conjugacy classes (columns). It equals
// NumIrreducibles for a complete table.
func (ct *CharacterTable) NumClasses() int { return len(ct.classes) }

// Dimensions returns the degrees of the irreducible characters, in row order.
func (ct *CharacterTable) Dimensions() []int {
	return append([]int(nil), ct.dims...)
}

// ClassSizes returns the number of elements in each conjugacy class, in column
// order.
func (ct *CharacterTable) ClassSizes() []int {
	return append([]int(nil), ct.classSizes...)
}

// ClassRepresentatives returns the smallest element index of each conjugacy
// class, in column order.
func (ct *CharacterTable) ClassRepresentatives() []int {
	return append([]int(nil), ct.reps...)
}

// ClassElementOrders returns the order of a representative element of each
// conjugacy class, in column order.
func (ct *CharacterTable) ClassElementOrders() []int {
	return append([]int(nil), ct.orders...)
}

// Irreducible returns the i-th irreducible character as a full class function
// over the group.
func (ct *CharacterTable) Irreducible(i int) Character {
	return ct.irr[i].Clone()
}

// Irreducibles returns all irreducible characters as full class functions.
func (ct *CharacterTable) Irreducibles() []Character {
	out := make([]Character, len(ct.irr))
	for i := range ct.irr {
		out[i] = ct.irr[i].Clone()
	}
	return out
}

// Value returns the value of irreducible i on conjugacy class c.
func (ct *CharacterTable) Value(i, c int) complex128 { return ct.values[i][c] }

// RowOrthogonality reports whether the rows satisfy the first orthogonality
// relation: 〈χᵢ, χⱼ〉 = δᵢⱼ to within tol, using the weighted class sum
// (1/|G|) Σ_c |class c|·conj(χᵢ)·χⱼ.
func (ct *CharacterTable) RowOrthogonality(tol float64) bool {
	n := ct.NumIrreducibles()
	order := float64(ct.group.Order())
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			var s complex128
			for c := range ct.classes {
				s += complex(float64(ct.classSizes[c]), 0) *
					cmplx.Conj(ct.values[i][c]) * ct.values[j][c]
			}
			s /= complex(order, 0)
			want := complex128(0)
			if i == j {
				want = 1
			}
			if cmplx.Abs(s-want) > tol {
				return false
			}
		}
	}
	return true
}

// ColumnOrthogonality reports whether the columns satisfy the second
// orthogonality relation: for classes c and d with representatives of
// centralizer order |G|/|class|,
//
//	Σᵢ conj(χᵢ(c))·χᵢ(d) = δ_{cd}·|G|/|class c|,
//
// to within tol.
func (ct *CharacterTable) ColumnOrthogonality(tol float64) bool {
	n := ct.NumIrreducibles()
	order := float64(ct.group.Order())
	for c := range ct.classes {
		for d := range ct.classes {
			var s complex128
			for i := 0; i < n; i++ {
				s += cmplx.Conj(ct.values[i][c]) * ct.values[i][d]
			}
			want := complex128(0)
			if c == d {
				want = complex(order/float64(ct.classSizes[c]), 0)
			}
			if cmplx.Abs(s-want) > tol {
				return false
			}
		}
	}
	return true
}

// SumOfSquaresOfDims returns Σ dᵢ², which for a complete table equals the group
// order |G|.
func (ct *CharacterTable) SumOfSquaresOfDims() int {
	s := 0
	for _, d := range ct.dims {
		s += d * d
	}
	return s
}

// Decompose returns the multiplicities of an arbitrary character chi against
// the irreducibles of the table, so that chi = Σ mult[i]·χᵢ.
func (ct *CharacterTable) Decompose(chi Character) []int {
	return DecomposeCharacter(ct.group, chi, ct.irr)
}

// IsComplete reports whether the table has as many irreducibles as classes and
// satisfies Σ dᵢ² = |G|.
func (ct *CharacterTable) IsComplete() bool {
	return ct.NumIrreducibles() == ct.NumClasses() &&
		ct.SumOfSquaresOfDims() == ct.group.Order()
}

// String renders the character table as an aligned text grid: one row per
// irreducible, one column per class, values rounded to two decimals.
func (ct *CharacterTable) String() string {
	var b strings.Builder
	fmt.Fprintf(&b, "Character table of %s\n", ct.group.Name())
	b.WriteString("class size:")
	for _, s := range ct.classSizes {
		fmt.Fprintf(&b, "%8d", s)
	}
	b.WriteString("\n")
	b.WriteString("elem order:")
	for _, o := range ct.orders {
		fmt.Fprintf(&b, "%8d", o)
	}
	b.WriteString("\n")
	for i := range ct.irr {
		fmt.Fprintf(&b, "chi_%d (%d):", i, ct.dims[i])
		for c := range ct.classes {
			z := RoundC(ct.values[i][c], 2)
			if math.Abs(imag(z)) < 1e-9 {
				fmt.Fprintf(&b, "%8g", real(z))
			} else {
				fmt.Fprintf(&b, "  %4.1f%+.1fi", real(z), imag(z))
			}
		}
		b.WriteString("\n")
	}
	return b.String()
}

// CharacterTableCyclic returns the character table of the cyclic group C_n. Its
// n irreducible characters send the generator to each n-th root of unity. It
// panics if n < 1.
func CharacterTableCyclic(n int) *CharacterTable {
	g := CyclicGroup(n)
	irrs := make([]Character, n)
	for k := 0; k < n; k++ {
		chi := make(Character, n)
		for i := 0; i < n; i++ {
			chi[i] = RootOfUnity(n, (k*i)%n)
		}
		irrs[k] = chi
	}
	return BuildCharacterTable(g, irrs)
}

// CharacterTableKleinFour returns the character table of the Klein four-group.
func CharacterTableKleinFour() *CharacterTable {
	g := KleinFourGroup()
	// C2×C2: elements (a,b), a,b∈{0,1}. Four 1-dim reps χ_{s,t}(a,b)=(-1)^{sa+tb}.
	irrs := make([]Character, 4)
	k := 0
	for s := 0; s < 2; s++ {
		for t := 0; t < 2; t++ {
			chi := make(Character, 4)
			idx := 0
			for a := 0; a < 2; a++ {
				for b := 0; b < 2; b++ {
					sign := 1.0
					if (s*a+t*b)%2 == 1 {
						sign = -1
					}
					chi[idx] = complex(sign, 0)
					idx++
				}
			}
			irrs[k] = chi
			k++
		}
	}
	return BuildCharacterTable(g, irrs)
}

// CharacterTableDihedral returns the character table of the dihedral group D_n.
// It is assembled from the one-dimensional characters and the two-dimensional
// characters χ_h(r^j) = 2cos(2πhj/n), χ_h(reflection) = 0. It panics if n < 1.
func CharacterTableDihedral(n int) *CharacterTable {
	g := DihedralGroup(n)
	idx := func(a, b int) int { return a*n + ((b%n)+n)%n }
	var irrs []Character

	one := func(f func(a, b int) float64) Character {
		chi := make(Character, 2*n)
		for a := 0; a < 2; a++ {
			for b := 0; b < n; b++ {
				chi[idx(a, b)] = complex(f(a, b), 0)
			}
		}
		return chi
	}
	// Trivial and the r↦1, s↦-1 rep exist for every n.
	irrs = append(irrs, one(func(a, b int) float64 { return 1 }))
	irrs = append(irrs, one(func(a, b int) float64 {
		if a == 1 {
			return -1
		}
		return 1
	}))
	if n%2 == 0 {
		irrs = append(irrs, one(func(a, b int) float64 {
			if b%2 == 1 {
				return -1
			}
			return 1
		}))
		irrs = append(irrs, one(func(a, b int) float64 {
			if (a+b)%2 == 1 {
				return -1
			}
			return 1
		}))
	}
	// Two-dimensional characters.
	var kmax int
	if n%2 == 0 {
		kmax = n/2 - 1
	} else {
		kmax = (n - 1) / 2
	}
	for h := 1; h <= kmax; h++ {
		chi := make(Character, 2*n)
		for b := 0; b < n; b++ {
			val := 2 * math.Cos(2*math.Pi*float64(h)*float64(b)/float64(n))
			chi[idx(0, b)] = complex(val, 0)
			chi[idx(1, b)] = 0
		}
		irrs = append(irrs, chi)
	}
	return BuildCharacterTable(g, irrs)
}

// CharacterTableQuaternion returns the character table of the quaternion group
// Q8: four one-dimensional characters and one two-dimensional character with
// values (2, -2, 0, 0, 0) on the classes {1}, {-1}, {±i}, {±j}, {±k}.
func CharacterTableQuaternion() *CharacterTable {
	g := QuaternionGroup()
	// Element indices: 0:1,1:-1,2:i,3:-i,4:j,5:-j,6:k,7:-k.
	// 1-dim characters take value on {1,-1,i,-i,j,-j,k,-k}.
	rows := [][]float64{
		{1, 1, 1, 1, 1, 1, 1, 1},     // trivial
		{1, 1, 1, 1, -1, -1, -1, -1}, // i↦1, j↦-1
		{1, 1, -1, -1, 1, 1, -1, -1}, // i↦-1, j↦1
		{1, 1, -1, -1, -1, -1, 1, 1}, // i↦-1, j↦-1
	}
	var irrs []Character
	for _, r := range rows {
		chi := make(Character, 8)
		for i, v := range r {
			chi[i] = complex(v, 0)
		}
		irrs = append(irrs, chi)
	}
	// Two-dimensional character: 2 on 1, -2 on -1, 0 on the rest.
	two := make(Character, 8)
	two[0] = 2
	two[1] = -2
	irrs = append(irrs, two)
	return BuildCharacterTable(g, irrs)
}

// CharacterTableSymmetric returns the character table of the symmetric group
// S_n for 1 <= n <= 4. It panics for n outside this range (the general case is
// governed by the Murnaghan-Nakayama rule and is not implemented here).
func CharacterTableSymmetric(n int) *CharacterTable {
	g, perms := symmetricData(n)
	switch n {
	case 1:
		return BuildCharacterTable(g, []Character{TrivialCharacter(g)})
	case 2:
		return BuildCharacterTable(g, c2WrapS2(g, perms))
	}
	// For n = 3, 4 assign values by cycle type.
	byType := symmetricCharacterRows(n)
	var irrs []Character
	for _, row := range byType {
		chi := make(Character, len(perms))
		for i, p := range perms {
			key := cycleTypeKey(p.CycleType())
			chi[i] = complex(row[key], 0)
		}
		irrs = append(irrs, chi)
	}
	return BuildCharacterTable(g, irrs)
}

// symmetricCharacterRows returns, for n = 3 or 4, the irreducible character
// values keyed by cycle-type signature. Values are the classical character
// tables of S3 and S4.
func symmetricCharacterRows(n int) []map[string]float64 {
	switch n {
	case 3:
		// classes by cycle type: [1,1,1] identity, [2,1] transposition, [3] 3-cycle.
		e := cycleTypeKey([]int{1, 1, 1})
		t := cycleTypeKey([]int{2, 1})
		c := cycleTypeKey([]int{3})
		return []map[string]float64{
			{e: 1, t: 1, c: 1},  // trivial
			{e: 1, t: -1, c: 1}, // sign
			{e: 2, t: 0, c: -1}, // standard (2-dim)
		}
	case 4:
		e := cycleTypeKey([]int{1, 1, 1, 1})
		t := cycleTypeKey([]int{2, 1, 1})
		d := cycleTypeKey([]int{2, 2})
		c3 := cycleTypeKey([]int{3, 1})
		c4 := cycleTypeKey([]int{4})
		return []map[string]float64{
			{e: 1, t: 1, d: 1, c3: 1, c4: 1},   // trivial
			{e: 1, t: -1, d: 1, c3: 1, c4: -1}, // sign
			{e: 2, t: 0, d: 2, c3: -1, c4: 0},  // 2-dim
			{e: 3, t: 1, d: -1, c3: 0, c4: -1}, // standard (3-dim)
			{e: 3, t: -1, d: -1, c3: 0, c4: 1}, // standard ⊗ sign
		}
	}
	panic("grouprep: symmetricCharacterRows only for n in {3,4}")
}

// cycleTypeKey renders a cycle-type partition as a stable string key.
func cycleTypeKey(t []int) string {
	parts := append([]int(nil), t...)
	sort.Sort(sort.Reverse(sort.IntSlice(parts)))
	var b strings.Builder
	for i, v := range parts {
		if i > 0 {
			b.WriteByte('.')
		}
		fmt.Fprintf(&b, "%d", v)
	}
	return b.String()
}

// c2WrapS2 builds the S2 table directly (S2 ≅ C2) using the permutation order.
func c2WrapS2(g *Group, perms []Perm) []Character {
	triv := make(Character, len(perms))
	sgn := make(Character, len(perms))
	for i, p := range perms {
		triv[i] = 1
		sgn[i] = complex(float64(p.Sign()), 0)
	}
	return []Character{triv, sgn}
}
