package satsolver

import (
	"sort"
	"strings"
)

// PopCount returns the number of one bits (the Hamming weight) of x.
func PopCount(x uint) int {
	c := 0
	for x != 0 {
		x &= x - 1
		c++
	}
	return c
}

// BitString returns the binary representation of x in exactly numVars digits,
// most-significant bit first.
func BitString(x, numVars int) string {
	b := make([]byte, numVars)
	for i := 0; i < numVars; i++ {
		if (x>>(numVars-1-i))&1 == 1 {
			b[i] = '1'
		} else {
			b[i] = '0'
		}
	}
	return string(b)
}

// GrayEncode returns the reflected binary Gray code of n, equal to n XOR (n>>1).
func GrayEncode(n int) int { return n ^ (n >> 1) }

// GrayDecode inverts [GrayEncode].
func GrayDecode(g int) int {
	n := g
	for g >>= 1; g != 0; g >>= 1 {
		n ^= g
	}
	return n
}

// GrayCode returns the reflected binary Gray-code sequence for n bits.
func GrayCode(n int) []int {
	out := make([]int, 1<<n)
	for i := range out {
		out[i] = GrayEncode(i)
	}
	return out
}

// Implicant is a Boolean cube: a product term covering a set of minterms. Value
// holds the fixed bit pattern (dash positions are zero) and Mask marks the
// positions that are "don't care" (combined away), so the cube covers exactly
// the assignments x with (x &^ Mask) == Value.
type Implicant struct {
	// Value is the fixed part of the cube; bits at Mask positions are zero.
	Value uint
	// Mask marks eliminated (dash) positions with a set bit.
	Mask uint
	// NumVars is the number of variables in the domain.
	NumVars int
}

// Covers reports whether the cube covers the minterm m.
func (im Implicant) Covers(m int) bool {
	return uint(m)&^im.Mask == im.Value
}

// Size returns the number of minterms the cube covers, which is 2 raised to the
// number of dash positions.
func (im Implicant) Size() int {
	return 1 << PopCount(im.Mask)
}

// Minterms returns the sorted list of minterm indices covered by the cube.
func (im Implicant) Minterms() []int {
	var out []int
	for m := 0; m < (1 << im.NumVars); m++ {
		if im.Covers(m) {
			out = append(out, m)
		}
	}
	return out
}

// String renders the cube as a pattern of '0', '1' and '-' characters,
// most-significant variable first.
func (im Implicant) String() string {
	b := make([]byte, im.NumVars)
	for i := 0; i < im.NumVars; i++ {
		bit := im.NumVars - 1 - i
		switch {
		case (im.Mask>>bit)&1 == 1:
			b[i] = '-'
		case (im.Value>>bit)&1 == 1:
			b[i] = '1'
		default:
			b[i] = '0'
		}
	}
	return string(b)
}

// Term renders the cube as a product term using the given variable names, most
// significant first, e.g. "a & ~c". The all-dash cube renders as "1".
func (im Implicant) Term(names []string) string {
	var parts []string
	for i := 0; i < im.NumVars; i++ {
		bit := im.NumVars - 1 - i
		if (im.Mask>>bit)&1 == 1 {
			continue
		}
		name := names[i]
		if (im.Value>>bit)&1 == 1 {
			parts = append(parts, name)
		} else {
			parts = append(parts, "~"+name)
		}
	}
	if len(parts) == 0 {
		return "1"
	}
	return strings.Join(parts, " & ")
}

func combineImplicants(a, b Implicant) (Implicant, bool) {
	if a.Mask != b.Mask {
		return Implicant{}, false
	}
	diff := a.Value ^ b.Value
	if PopCount(diff) != 1 {
		return Implicant{}, false
	}
	return Implicant{
		Value:   a.Value &^ diff,
		Mask:    a.Mask | diff,
		NumVars: a.NumVars,
	}, true
}

func implicantKey(im Implicant) [2]uint { return [2]uint{im.Value, im.Mask} }

// PrimeImplicants returns the prime implicants of the Boolean function defined
// by the given on-set minterms and don't-care terms over numVars variables,
// computed by iterated Quine-McCluskey combination.
func PrimeImplicants(minterms, dontcares []int, numVars int) []Implicant {
	seen := map[int]bool{}
	var all []Implicant
	add := func(m int) {
		if seen[m] {
			return
		}
		seen[m] = true
		all = append(all, Implicant{Value: uint(m), Mask: 0, NumVars: numVars})
	}
	for _, m := range minterms {
		add(m)
	}
	for _, m := range dontcares {
		add(m)
	}

	primeSet := map[[2]uint]Implicant{}
	current := all
	for len(current) > 0 {
		used := make([]bool, len(current))
		nextSet := map[[2]uint]Implicant{}
		for i := 0; i < len(current); i++ {
			for j := i + 1; j < len(current); j++ {
				if c, ok := combineImplicants(current[i], current[j]); ok {
					used[i] = true
					used[j] = true
					nextSet[implicantKey(c)] = c
				}
			}
		}
		for i, im := range current {
			if !used[i] {
				primeSet[implicantKey(im)] = im
			}
		}
		current = current[:0]
		for _, im := range nextSet {
			current = append(current, im)
		}
		sort.Slice(current, func(a, b int) bool {
			if current[a].Value != current[b].Value {
				return current[a].Value < current[b].Value
			}
			return current[a].Mask < current[b].Mask
		})
	}

	out := make([]Implicant, 0, len(primeSet))
	for _, im := range primeSet {
		out = append(out, im)
	}
	sort.Slice(out, func(a, b int) bool {
		if out[a].Value != out[b].Value {
			return out[a].Value < out[b].Value
		}
		return out[a].Mask < out[b].Mask
	})
	return out
}

// QuineMcCluskey is an alias for [PrimeImplicants], returning all prime
// implicants of the function.
func QuineMcCluskey(minterms, dontcares []int, numVars int) []Implicant {
	return PrimeImplicants(minterms, dontcares, numVars)
}

// EssentialPrimeImplicants returns the prime implicants that are the sole cover
// of at least one required minterm (a minterm in the on-set, ignoring
// don't-cares).
func EssentialPrimeImplicants(minterms, dontcares []int, numVars int) []Implicant {
	primes := PrimeImplicants(minterms, dontcares, numVars)
	essentials, _ := selectEssentials(primes, minterms)
	return essentials
}

func selectEssentials(primes []Implicant, minterms []int) (essential []Implicant, covered map[int]bool) {
	covered = map[int]bool{}
	chosen := map[[2]uint]bool{}
	for _, m := range minterms {
		var coverers []Implicant
		for _, p := range primes {
			if p.Covers(m) {
				coverers = append(coverers, p)
			}
		}
		if len(coverers) == 1 {
			p := coverers[0]
			if !chosen[implicantKey(p)] {
				chosen[implicantKey(p)] = true
				essential = append(essential, p)
			}
		}
	}
	for _, p := range essential {
		for _, m := range minterms {
			if p.Covers(m) {
				covered[m] = true
			}
		}
	}
	return essential, covered
}

// MinimizeSOP returns a minimal (greedy) sum-of-products cover of the function
// as a list of implicants: all essential prime implicants plus a greedy
// selection of further primes to cover the remaining minterms.
func MinimizeSOP(minterms, dontcares []int, numVars int) []Implicant {
	primes := PrimeImplicants(minterms, dontcares, numVars)
	essentials, covered := selectEssentials(primes, minterms)
	result := append([]Implicant(nil), essentials...)
	chosen := map[[2]uint]bool{}
	for _, e := range essentials {
		chosen[implicantKey(e)] = true
	}

	remaining := map[int]bool{}
	for _, m := range minterms {
		if !covered[m] {
			remaining[m] = true
		}
	}
	for len(remaining) > 0 {
		var best Implicant
		bestCount := -1
		bestKey := [2]uint{}
		for _, p := range primes {
			if chosen[implicantKey(p)] {
				continue
			}
			cnt := 0
			for m := range remaining {
				if p.Covers(m) {
					cnt++
				}
			}
			if cnt > bestCount {
				bestCount = cnt
				best = p
				bestKey = implicantKey(p)
			}
		}
		if bestCount <= 0 {
			break
		}
		chosen[bestKey] = true
		result = append(result, best)
		for m := range remaining {
			if best.Covers(m) {
				delete(remaining, m)
			}
		}
	}
	sort.Slice(result, func(a, b int) bool {
		if result[a].Value != result[b].Value {
			return result[a].Value < result[b].Value
		}
		return result[a].Mask < result[b].Mask
	})
	return result
}

// SOPString renders a list of implicants as a sum-of-products Boolean
// expression using the given variable names, most significant first. An empty
// cover renders as "0".
func SOPString(cover []Implicant, names []string) string {
	if len(cover) == 0 {
		return "0"
	}
	parts := make([]string, len(cover))
	for i, im := range cover {
		t := im.Term(names)
		if strings.Contains(t, " & ") {
			parts[i] = "(" + t + ")"
		} else {
			parts[i] = t
		}
	}
	return strings.Join(parts, " | ")
}

// MinimizeSOPExpr returns the minimised sum-of-products cover of the function as
// a Boolean [Expr] built from the given variable names.
func MinimizeSOPExpr(minterms, dontcares []int, names []string) Expr {
	numVars := len(names)
	cover := MinimizeSOP(minterms, dontcares, numVars)
	if len(cover) == 0 {
		return False
	}
	terms := make([]Expr, 0, len(cover))
	for _, im := range cover {
		var lits []Expr
		for i := 0; i < numVars; i++ {
			bit := numVars - 1 - i
			if (im.Mask>>bit)&1 == 1 {
				continue
			}
			if (im.Value>>bit)&1 == 1 {
				lits = append(lits, Variable(names[i]))
			} else {
				lits = append(lits, Not{X: Variable(names[i])})
			}
		}
		if len(lits) == 0 {
			return True
		}
		terms = append(terms, AndAll(lits...))
	}
	return OrAll(terms...)
}

// KarnaughMap is a two-dimensional Gray-coded arrangement of the truth table of
// a Boolean function, the layout used for manual minimisation.
type KarnaughMap struct {
	// NumVars is the number of variables.
	NumVars int
	// RowVars and ColVars are the numbers of variables mapped to rows and
	// columns.
	RowVars, ColVars int
	// Cells holds the function value at each grid position: 0, 1, or 2 for a
	// don't-care.
	Cells [][]int
	// RowCodes and ColCodes give the Gray-coded axis labels.
	RowCodes, ColCodes []int
}

// NewKarnaughMap builds the Karnaugh map of the function defined by minterms and
// dontcares over numVars variables. The higher-order half of the variables
// index rows and the lower-order half index columns.
func NewKarnaughMap(minterms, dontcares []int, numVars int) KarnaughMap {
	rowVars := numVars / 2
	colVars := numVars - rowVars
	rows := 1 << rowVars
	cols := 1 << colVars
	on := map[int]bool{}
	dc := map[int]bool{}
	for _, m := range minterms {
		on[m] = true
	}
	for _, m := range dontcares {
		dc[m] = true
	}
	rowCodes := GrayCode(rowVars)
	colCodes := GrayCode(colVars)
	cells := make([][]int, rows)
	for r := 0; r < rows; r++ {
		cells[r] = make([]int, cols)
		for c := 0; c < cols; c++ {
			m := (rowCodes[r] << colVars) | colCodes[c]
			switch {
			case on[m]:
				cells[r][c] = 1
			case dc[m]:
				cells[r][c] = 2
			default:
				cells[r][c] = 0
			}
		}
	}
	return KarnaughMap{
		NumVars:  numVars,
		RowVars:  rowVars,
		ColVars:  colVars,
		Cells:    cells,
		RowCodes: rowCodes,
		ColCodes: colCodes,
	}
}

// String renders the Karnaugh map as a grid of '0', '1' and 'X' (don't-care)
// characters with Gray-coded axis labels.
func (k KarnaughMap) String() string {
	var b strings.Builder
	b.WriteString("     ")
	for _, c := range k.ColCodes {
		b.WriteString(BitString(c, k.ColVars))
		b.WriteByte(' ')
	}
	b.WriteByte('\n')
	for r, row := range k.Cells {
		b.WriteString(BitString(k.RowCodes[r], k.RowVars))
		b.WriteString("  ")
		for _, v := range row {
			switch v {
			case 1:
				b.WriteString("1")
			case 2:
				b.WriteString("X")
			default:
				b.WriteString("0")
			}
			b.WriteString(strings.Repeat(" ", k.ColVars))
		}
		b.WriteByte('\n')
	}
	return b.String()
}
