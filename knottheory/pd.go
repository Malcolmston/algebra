package knottheory

import "fmt"

// PDCrossing is one crossing of a planar-diagram (PD) code. The four edge
// labels A, B, C, D are listed counter-clockwise starting from the incoming
// under-strand: A is the incoming under-arc, C the outgoing under-arc, and B, D
// the over-arc. Sign is the writhe contribution of the crossing (+1 or -1).
type PDCrossing struct {
	A, B, C, D int
	Sign       int
}

// PDCode is a knot or link diagram in planar-diagram notation together with a
// count of any crossingless unknotted circles. Edges are labelled by positive
// integers, each of which appears exactly twice across all crossings.
type PDCode struct {
	Crossings []PDCrossing
	Circles   int
}

// NumCrossings returns the number of crossings of the diagram.
func (pd PDCode) NumCrossings() int { return len(pd.Crossings) }

// Writhe returns the writhe of the diagram, the sum of the crossing signs.
func (pd PDCode) Writhe() int {
	w := 0
	for _, c := range pd.Crossings {
		w += c.Sign
	}
	return w
}

// Validate checks that every edge label occurs exactly twice and that every
// sign is +1 or -1. A diagram with only circles and no crossings is valid.
func (pd PDCode) Validate() error {
	count := map[int]int{}
	for _, c := range pd.Crossings {
		if c.Sign != 1 && c.Sign != -1 {
			return fmt.Errorf("knottheory: crossing has invalid sign %d", c.Sign)
		}
		for _, e := range []int{c.A, c.B, c.C, c.D} {
			if e <= 0 {
				return fmt.Errorf("knottheory: edge labels must be positive, got %d", e)
			}
			count[e]++
		}
	}
	for e, n := range count {
		if n != 2 {
			return fmt.Errorf("knottheory: edge %d appears %d times, expected 2", e, n)
		}
	}
	return nil
}

// Mirror returns the PD code of the mirror image, obtained by swapping the
// over- and under-arcs (A,B,C,D -> B,C,D,A) and negating every sign.
func (pd PDCode) Mirror() PDCode {
	cs := make([]PDCrossing, len(pd.Crossings))
	for i, c := range pd.Crossings {
		cs[i] = PDCrossing{A: c.B, B: c.C, C: c.D, D: c.A, Sign: -c.Sign}
	}
	return PDCode{Crossings: cs, Circles: pd.Circles}
}

// dsu is a tiny union-find structure used by the Kauffman bracket state sum.
type dsu struct{ parent []int }

func newDSU(n int) *dsu {
	p := make([]int, n)
	for i := range p {
		p[i] = i
	}
	return &dsu{parent: p}
}

func (d *dsu) find(x int) int {
	for d.parent[x] != x {
		d.parent[x] = d.parent[d.parent[x]]
		x = d.parent[x]
	}
	return x
}

func (d *dsu) union(a, b int) {
	ra, rb := d.find(a), d.find(b)
	if ra != rb {
		d.parent[ra] = rb
	}
}

// stateLoops returns the number of loops in the smoothed diagram for a given
// state. The state bit for crossing i is 0 for the A-smoothing and 1 for the
// B-smoothing.
func (pd PDCode) stateLoops(state uint) int {
	n := len(pd.Crossings)
	if n == 0 {
		return pd.Circles
	}
	d := newDSU(4 * n)
	// Edge connections: ports carrying the same edge label are joined.
	edgePorts := map[int][]int{}
	for i, c := range pd.Crossings {
		edges := [4]int{c.A, c.B, c.C, c.D}
		for k := 0; k < 4; k++ {
			p := 4*i + k
			edgePorts[edges[k]] = append(edgePorts[edges[k]], p)
		}
	}
	for _, ports := range edgePorts {
		for i := 1; i < len(ports); i++ {
			d.union(ports[0], ports[i])
		}
	}
	// Smoothing connections at each crossing.
	for i := 0; i < n; i++ {
		base := 4 * i
		if state&(1<<uint(i)) == 0 {
			// A-smoothing: joins the over-arc to the arc swept
			// counter-clockwise, i.e. B-C and D-A.
			d.union(base+1, base+2)
			d.union(base+3, base+0)
		} else {
			// B-smoothing: A-B and C-D.
			d.union(base+0, base+1)
			d.union(base+2, base+3)
		}
	}
	roots := map[int]bool{}
	for p := 0; p < 4*n; p++ {
		roots[d.find(p)] = true
	}
	return len(roots) + pd.Circles
}

// KauffmanBracket returns the Kauffman bracket polynomial of the diagram as a
// Laurent polynomial in the variable A. It is computed by the state sum
//
//	<D> = sum over states A^(a-b) * (-A^2 - A^-2)^(loops-1),
//
// where a and b count A- and B-smoothings and loops is the number of circles in
// the smoothed state. The bracket of the unknot is 1.
func (pd PDCode) KauffmanBracket() Laurent {
	n := len(pd.Crossings)
	delta := Monomial(-1, 2).Add(Monomial(-1, -2)) // -A^2 - A^-2
	result := ZeroLaurent()
	for state := uint(0); state < (uint(1) << uint(n)); state++ {
		a := 0
		for i := 0; i < n; i++ {
			if state&(1<<uint(i)) == 0 {
				a++
			}
		}
		b := n - a
		loops := pd.stateLoops(state)
		term := Monomial(1, a-b)
		if loops-1 > 0 {
			term = term.Mul(delta.Pow(loops - 1))
		} else if loops-1 < 0 {
			// loops is always >= 1 for a non-empty diagram, so this cannot
			// occur; guarded for completeness.
			term = ZeroLaurent()
		}
		result = result.Add(term)
	}
	if n == 0 {
		// Only circles: <k circles> = delta^(circles-1).
		if pd.Circles <= 0 {
			return OneLaurent()
		}
		return delta.Pow(pd.Circles - 1)
	}
	return result
}

// NormalizedBracket returns the writhe-normalised Kauffman bracket
// (-A^3)^(-w) <D>, an ambient-isotopy invariant (in the variable A) from which
// the Jones polynomial is obtained by the substitution A = t^(-1/4).
func (pd PDCode) NormalizedBracket() Laurent {
	br := pd.KauffmanBracket()
	w := pd.Writhe()
	res := br.ShiftExp(-3 * w)
	if ((w%2)+2)%2 == 1 {
		res = res.Neg()
	}
	return res
}

// JonesPolynomialSqrt returns the Jones polynomial expressed in the variable
// u = t^(1/2), so that a term c*u^k denotes c*t^(k/2). Every exponent is an
// integer, which makes this form valid for links as well as knots.
func (pd PDCode) JonesPolynomialSqrt() Laurent {
	nb := pd.NormalizedBracket()
	res := ZeroLaurent()
	for _, term := range nb.Terms() {
		// A = t^(-1/4): A^e = t^(-e/4) = u^(-e/2).
		res = res.Add(Monomial(term.Coeff, -term.Exp/2))
	}
	return res
}

// JonesPolynomial returns the Jones polynomial in the variable t. It returns an
// error when the polynomial contains genuine half-integer powers of t (which
// happens for links with an even number of components); use JonesPolynomialSqrt
// in that case.
func (pd PDCode) JonesPolynomial() (Laurent, error) {
	nb := pd.NormalizedBracket()
	res := ZeroLaurent()
	for _, term := range nb.Terms() {
		if ((term.Exp%4)+4)%4 != 0 {
			return Laurent{}, fmt.Errorf("knottheory: Jones polynomial has half-integer powers of t; use JonesPolynomialSqrt")
		}
		res = res.Add(Monomial(term.Coeff, -term.Exp/4))
	}
	return res, nil
}

// arcRoots merges the two over-arcs of every crossing and returns a mapping
// from each edge label to a compact arc id, together with the total number of
// arcs. Two edges lie on the same arc when the strand runs over from one to the
// other without an intervening under-crossing.
func (pd PDCode) arcRoots() (edgeArc map[int]int, arcCount int) {
	edges := map[int]int{}
	idx := 0
	for _, c := range pd.Crossings {
		for _, e := range []int{c.A, c.B, c.C, c.D} {
			if _, ok := edges[e]; !ok {
				edges[e] = idx
				idx++
			}
		}
	}
	d := newDSU(idx)
	for _, c := range pd.Crossings {
		d.union(edges[c.B], edges[c.D]) // over-arc runs through
	}
	arcID := map[int]int{}
	edgeArc = map[int]int{}
	for e, i := range edges {
		r := d.find(i)
		id, ok := arcID[r]
		if !ok {
			id = arcCount
			arcID[r] = id
			arcCount++
		}
		edgeArc[e] = id
	}
	return edgeArc, arcCount
}

// ArcCount returns the number of arcs of the diagram (maximal over-strands
// between under-crossings).
func (pd PDCode) ArcCount() int {
	_, n := pd.arcRoots()
	return n
}

// AlexanderMatrix returns the Alexander matrix of the diagram over the Laurent
// ring in the variable t. It has one row per crossing and one column per arc;
// each crossing contributes the Fox relation among its over-arc and its two
// under-arcs.
func (pd PDCode) AlexanderMatrix() *LaurentMatrix {
	edgeArc, arcCount := pd.arcRoots()
	n := len(pd.Crossings)
	M := NewLaurentMatrix(n, arcCount)
	oneMinusT := OneLaurent().Sub(Monomial(1, 1)) // 1 - t
	t := Monomial(1, 1)
	negOne := Monomial(-1, 0)
	for r, c := range pd.Crossings {
		over := edgeArc[c.B]
		underIn := edgeArc[c.A]
		underOut := edgeArc[c.C]
		var wOver, wIn, wOut Laurent
		if c.Sign > 0 {
			wOut, wIn, wOver = negOne, t, oneMinusT
		} else {
			wIn, wOut, wOver = negOne, t, oneMinusT
		}
		M.Set(r, over, M.At(r, over).Add(wOver))
		M.Set(r, underIn, M.At(r, underIn).Add(wIn))
		M.Set(r, underOut, M.At(r, underOut).Add(wOut))
	}
	return M
}

// AlexanderPolynomial returns the Alexander polynomial of the diagram computed
// from its Alexander matrix as the determinant of the minor obtained by
// deleting the last row and last column, normalised to the symmetric
// representative with value +1 at t=1. It returns an error for a diagram with
// fewer than one crossing (use the unknot value 1 directly).
func (pd PDCode) AlexanderPolynomial() (Laurent, error) {
	n := len(pd.Crossings)
	if n == 0 {
		return OneLaurent(), nil
	}
	M := pd.AlexanderMatrix()
	if M.Cols() != n {
		return Laurent{}, fmt.Errorf("knottheory: Alexander matrix is %dx%d, expected square", M.Rows(), M.Cols())
	}
	det := minorDet(M.data, n-1, n-1)
	if det.IsZero() {
		return ZeroLaurent(), nil
	}
	return normalizeAlexander(det), nil
}

// KnotDeterminant returns the determinant of the diagram, the absolute value of
// the Alexander polynomial evaluated at t = -1. It equals the order of the first
// homology of the double branched cover.
func (pd PDCode) KnotDeterminant() int {
	alex, err := pd.AlexanderPolynomial()
	if err != nil || alex.IsZero() {
		return 0
	}
	v := alex.EvalUnit(-1)
	if v < 0 {
		return -v
	}
	return v
}

// IsPColorable reports whether the diagram admits a non-trivial Fox p-colouring,
// which for a prime p is equivalent to p dividing the determinant.
func (pd PDCode) IsPColorable(p int) bool {
	if p <= 1 {
		return false
	}
	det := pd.KnotDeterminant()
	if det == 0 {
		return true
	}
	return det%p == 0
}

// IsThreeColorable reports whether the diagram is Fox 3-colourable.
func (pd PDCode) IsThreeColorable() bool { return pd.IsPColorable(3) }

// -------- Knot and link catalogue --------

// UnknotPD returns the standard crossingless diagram of the unknot.
func UnknotPD() PDCode { return PDCode{Circles: 1} }

// UnlinkPD returns the diagram of the n-component unlink (n crossingless
// circles).
func UnlinkPD(n int) PDCode { return PDCode{Circles: n} }

// TrefoilPD returns a right-handed trefoil (3_1) diagram with writhe +3.
func TrefoilPD() PDCode {
	return PDCode{Crossings: []PDCrossing{
		{A: 1, B: 4, C: 2, D: 5, Sign: 1},
		{A: 3, B: 6, C: 4, D: 1, Sign: 1},
		{A: 5, B: 2, C: 6, D: 3, Sign: 1},
	}}
}

// LeftTrefoilPD returns the left-handed trefoil, the mirror image of TrefoilPD.
func LeftTrefoilPD() PDCode { return TrefoilPD().Mirror() }

// FigureEightPD returns the figure-eight knot (4_1) diagram, an amphichiral knot
// with writhe 0.
func FigureEightPD() PDCode {
	return PDCode{Crossings: []PDCrossing{
		{A: 4, B: 2, C: 5, D: 1, Sign: -1},
		{A: 8, B: 6, C: 1, D: 5, Sign: -1},
		{A: 6, B: 3, C: 7, D: 4, Sign: 1},
		{A: 2, B: 7, C: 3, D: 8, Sign: 1},
	}}
}

// CinquefoilPD returns the (2,5) torus knot 5_1 (Solomon's seal knot) with
// writhe +5.
func CinquefoilPD() PDCode {
	return PDCode{Crossings: []PDCrossing{
		{A: 1, B: 6, C: 2, D: 7, Sign: 1},
		{A: 3, B: 8, C: 4, D: 9, Sign: 1},
		{A: 5, B: 10, C: 6, D: 1, Sign: 1},
		{A: 7, B: 2, C: 8, D: 3, Sign: 1},
		{A: 9, B: 4, C: 10, D: 5, Sign: 1},
	}}
}

// HopfLinkPD returns a Hopf link diagram. When positive is true both crossings
// are positive (writhe +2); otherwise both are negative (writhe -2).
func HopfLinkPD(positive bool) PDCode {
	s := 1
	if !positive {
		s = -1
	}
	return PDCode{Crossings: []PDCrossing{
		{A: 1, B: 3, C: 2, D: 4, Sign: s},
		{A: 3, B: 1, C: 4, D: 2, Sign: s},
	}}
}
