package logic

import (
	"fmt"
	"math/bits"
	"sort"
	"strings"
)

// PopCount returns the number of one bits (the Hamming weight) of x.
func PopCount(x uint) int { return bits.OnesCount(x) }

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

// GrayEncode returns the reflected binary Gray code of the non-negative integer
// n, computed as n XOR (n>>1). Consecutive values of n differ in exactly one
// bit of the result.
func GrayEncode(n int) int { return n ^ (n >> 1) }

// GrayDecode inverts [GrayEncode], recovering the integer whose Gray code is g.
func GrayDecode(g int) int {
	n := g
	for g >>= 1; g != 0; g >>= 1 {
		n ^= g
	}
	return n
}

// GrayCode returns the reflected binary Gray code sequence for n bits: the
// 2^n values GrayEncode(0), GrayEncode(1), ..., GrayEncode(2^n-1), in which
// adjacent entries (cyclically) differ in a single bit.
func GrayCode(n int) []int {
	out := make([]int, 1<<n)
	for i := range out {
		out[i] = GrayEncode(i)
	}
	return out
}

// Implicant is a product-term cube over numVars variables. Each variable
// position is either fixed to a value or eliminated (a "dash"). The zero value
// is not meaningful; construct implicants with [PrimeImplicants] or the
// minimisation routines.
type Implicant struct {
	// Value holds the fixed bit values; positions that are dashes are zero.
	Value uint32
	// Dashes marks eliminated positions with a one bit.
	Dashes uint32
	// NumVars is the number of variable positions in the cube.
	NumVars int

	covers []int
}

// Covers reports whether the cube contains the given minterm index.
func (im Implicant) Covers(minterm int) bool {
	mask := (^im.Dashes) & logicFullMask(im.NumVars)
	return (uint32(minterm) & mask) == (im.Value & mask)
}

// LiteralCount returns the number of fixed variable positions in the cube,
// i.e. the number of literals in the product term.
func (im Implicant) LiteralCount() int {
	return im.NumVars - bits.OnesCount32(im.Dashes&logicFullMask(im.NumVars))
}

// Minterms returns the original minterms that the cube covers, as recorded when
// it was produced by the minimisation routines. The slice is sorted ascending.
func (im Implicant) Minterms() []int {
	out := append([]int(nil), im.covers...)
	sort.Ints(out)
	return out
}

// String renders the cube as a pattern of numVars characters, most-significant
// position first, using '1', '0' and '-' for a dash.
func (im Implicant) String() string {
	b := make([]byte, im.NumVars)
	for i := 0; i < im.NumVars; i++ {
		pos := im.NumVars - 1 - i
		switch {
		case im.Dashes&(1<<uint(pos)) != 0:
			b[i] = '-'
		case im.Value&(1<<uint(pos)) != 0:
			b[i] = '1'
		default:
			b[i] = '0'
		}
	}
	return string(b)
}

// Expression renders the cube as a product term over the given variable names,
// where vars[0] is the most-significant position. A fixed 1 yields the plain
// variable, a fixed 0 its negation, and a dash omits the variable. The empty
// cube (all dashes) renders as the constant "T". It panics if len(vars) does
// not equal NumVars.
func (im Implicant) Expression(vars []string) string {
	if len(vars) != im.NumVars {
		panic(fmt.Sprintf("logic: Expression got %d vars, want %d", len(vars), im.NumVars))
	}
	var lits []string
	for i := 0; i < im.NumVars; i++ {
		pos := im.NumVars - 1 - i
		if im.Dashes&(1<<uint(pos)) != 0 {
			continue
		}
		if im.Value&(1<<uint(pos)) != 0 {
			lits = append(lits, vars[i])
		} else {
			lits = append(lits, "!"+vars[i])
		}
	}
	if len(lits) == 0 {
		return "T"
	}
	return strings.Join(lits, "&")
}

// logicFullMask returns a mask with the low numVars bits set.
func logicFullMask(numVars int) uint32 {
	if numVars >= 32 {
		return ^uint32(0)
	}
	return (uint32(1) << uint(numVars)) - 1
}

// logicImplKey is a canonical map key for an implicant cube.
func logicImplKey(value, dashes uint32) string {
	return fmt.Sprintf("%d:%d", value, dashes)
}

// logicSortImplicants orders cubes deterministically by dash mask then value.
func logicSortImplicants(s []Implicant) {
	sort.Slice(s, func(i, j int) bool {
		if s[i].Dashes != s[j].Dashes {
			return s[i].Dashes < s[j].Dashes
		}
		return s[i].Value < s[j].Value
	})
}

// PrimeImplicants computes the prime implicants of the Boolean function defined
// by the given minterms and don't-care terms over numVars variables, using the
// tabular Quine-McCluskey combination step. Don't-cares participate in the
// combination but the reported covers include only true minterms. The result is
// sorted deterministically.
func PrimeImplicants(minterms, dontCares []int, numVars int) []Implicant {
	mintermSet := map[int]bool{}
	for _, m := range minterms {
		mintermSet[m] = true
	}
	termSet := map[int]bool{}
	for _, m := range minterms {
		termSet[m] = true
	}
	for _, d := range dontCares {
		termSet[d] = true
	}

	var current []Implicant
	for t := range termSet {
		current = append(current, Implicant{Value: uint32(t), NumVars: numVars})
	}
	logicSortImplicants(current)

	var primes []Implicant
	seenPrime := map[string]bool{}
	for len(current) > 0 {
		used := make([]bool, len(current))
		var next []Implicant
		nextKey := map[string]bool{}
		for i := 0; i < len(current); i++ {
			for j := i + 1; j < len(current); j++ {
				a, b := current[i], current[j]
				if a.Dashes != b.Dashes {
					continue
				}
				diff := a.Value ^ b.Value
				if diff == 0 || diff&(diff-1) != 0 { // not a single-bit difference
					continue
				}
				if diff&a.Dashes != 0 {
					continue
				}
				used[i] = true
				used[j] = true
				nv := a.Value &^ diff
				nd := a.Dashes | diff
				key := logicImplKey(nv, nd)
				if !nextKey[key] {
					nextKey[key] = true
					next = append(next, Implicant{Value: nv, Dashes: nd, NumVars: numVars})
				}
			}
		}
		for i, im := range current {
			if used[i] {
				continue
			}
			key := logicImplKey(im.Value, im.Dashes)
			if !seenPrime[key] {
				seenPrime[key] = true
				primes = append(primes, im)
			}
		}
		logicSortImplicants(next)
		current = next
	}

	// Attach the true minterms each prime covers.
	uniqueMin := logicSortedKeys(mintermSet)
	for i := range primes {
		var cov []int
		for _, m := range uniqueMin {
			if primes[i].Covers(m) {
				cov = append(cov, m)
			}
		}
		primes[i].covers = cov
	}
	logicSortImplicants(primes)
	return primes
}

// logicSortedKeys returns the sorted keys of an int set.
func logicSortedKeys(set map[int]bool) []int {
	out := make([]int, 0, len(set))
	for k := range set {
		out = append(out, k)
	}
	sort.Ints(out)
	return out
}

// EssentialPrimeImplicants returns the essential prime implicants of the
// function defined by minterms and dontCares over numVars variables: those
// prime implicants that are the sole cover of some true minterm.
func EssentialPrimeImplicants(minterms, dontCares []int, numVars int) []Implicant {
	primes := PrimeImplicants(minterms, dontCares, numVars)
	mset := map[int]bool{}
	for _, m := range minterms {
		mset[m] = true
	}
	uniqueMin := logicSortedKeys(mset)
	essential := map[int]bool{}
	for _, m := range uniqueMin {
		coverer := -1
		count := 0
		for pi := range primes {
			if primes[pi].Covers(m) {
				count++
				coverer = pi
			}
		}
		if count == 1 {
			essential[coverer] = true
		}
	}
	var out []Implicant
	for pi := range primes {
		if essential[pi] {
			out = append(out, primes[pi])
		}
	}
	logicSortImplicants(out)
	return out
}

// QuineMcCluskey minimises the Boolean function given by minterms and dontCares
// over numVars variables and returns a sum-of-products cover as a set of
// implicants. It selects all essential prime implicants and then greedily adds
// prime implicants that cover the most still-uncovered minterms, breaking ties
// toward fewer literals for a deterministic result. An empty minterm list (the
// constant-false function) yields a nil cover.
func QuineMcCluskey(minterms, dontCares []int, numVars int) []Implicant {
	mset := map[int]bool{}
	for _, m := range minterms {
		mset[m] = true
	}
	if len(mset) == 0 {
		return nil
	}
	primes := PrimeImplicants(minterms, dontCares, numVars)

	remaining := map[int]bool{}
	for m := range mset {
		remaining[m] = true
	}
	chosen := map[int]bool{}
	var cover []Implicant

	// Essential prime implicants first.
	uniqueMin := logicSortedKeys(mset)
	for _, m := range uniqueMin {
		coverer := -1
		count := 0
		for pi := range primes {
			if primes[pi].Covers(m) {
				count++
				coverer = pi
			}
		}
		if count == 1 && !chosen[coverer] {
			chosen[coverer] = true
			cover = append(cover, primes[coverer])
			for _, cm := range primes[coverer].covers {
				delete(remaining, cm)
			}
		}
	}

	// Greedy cover of the rest.
	for len(remaining) > 0 {
		best := -1
		bestCount := -1
		bestLits := 1 << 30
		for pi := range primes {
			if chosen[pi] {
				continue
			}
			cnt := 0
			for _, cm := range primes[pi].covers {
				if remaining[cm] {
					cnt++
				}
			}
			if cnt == 0 {
				continue
			}
			lits := primes[pi].LiteralCount()
			if cnt > bestCount || (cnt == bestCount && lits < bestLits) {
				best = pi
				bestCount = cnt
				bestLits = lits
			}
		}
		if best == -1 {
			break // should not happen for a well-formed cover
		}
		chosen[best] = true
		cover = append(cover, primes[best])
		for _, cm := range primes[best].covers {
			delete(remaining, cm)
		}
	}

	logicSortImplicants(cover)
	return cover
}

// MinimizeSOP minimises the expression e and returns a minimal sum-of-products
// cover as a set of implicants over the variables of [Vars](e). The variable
// order is lexicographic, matching the minterm indexing used elsewhere in the
// package.
func MinimizeSOP(e Expr) []Implicant {
	vars := Vars(e)
	return QuineMcCluskey(Minterms(e), nil, len(vars))
}

// SOPString renders a cover produced by [QuineMcCluskey] or [MinimizeSOP] as a
// sum-of-products string over vars. An empty cover renders as the constant "F"
// and a cover containing the universal (all-dash) cube renders as "T".
func SOPString(cover []Implicant, vars []string) string {
	if len(cover) == 0 {
		return "F"
	}
	terms := make([]string, 0, len(cover))
	for _, im := range cover {
		t := im.Expression(vars)
		if t == "T" {
			return "T"
		}
		terms = append(terms, t)
	}
	return strings.Join(terms, " | ")
}

// MinimizeString minimises e and returns its minimal sum-of-products form as a
// string over the variables of e. A contradiction yields "F" and a tautology
// yields "T".
func MinimizeString(e Expr) string {
	return SOPString(MinimizeSOP(e), Vars(e))
}

// KarnaughMap is a Karnaugh map: a Gray-coded grid view of a Boolean function
// used to visualise adjacency and read off groupings.
type KarnaughMap struct {
	// NumVars is the number of variables.
	NumVars int
	// Vars holds optional variable labels (vars[0] most-significant). It may
	// be nil, in which case rendering uses positional headers only.
	Vars []string

	rowBits int
	colBits int
	cells   []bool
}

// NewKarnaughMap builds a Karnaugh map for numVars variables in which the given
// minterms are set to true. The higher-order half of the variables index the
// rows and the lower-order half index the columns, each in Gray-code order.
func NewKarnaughMap(minterms []int, numVars int) *KarnaughMap {
	colBits := (numVars + 1) / 2
	rowBits := numVars - colBits
	cells := make([]bool, 1<<numVars)
	for _, m := range minterms {
		if m >= 0 && m < len(cells) {
			cells[m] = true
		}
	}
	return &KarnaughMap{NumVars: numVars, rowBits: rowBits, colBits: colBits, cells: cells}
}

// Rows returns the number of rows in the map (2 raised to the row-variable
// count).
func (k *KarnaughMap) Rows() int { return 1 << k.rowBits }

// Cols returns the number of columns in the map (2 raised to the
// column-variable count).
func (k *KarnaughMap) Cols() int { return 1 << k.colBits }

// At returns the cell value at the given grid row and column. The row and
// column positions are interpreted through the Gray-code header ordering, so
// physically adjacent cells differ in a single variable.
func (k *KarnaughMap) At(row, col int) bool {
	rowVal := GrayEncode(row)
	colVal := GrayEncode(col)
	minterm := rowVal<<uint(k.colBits) | colVal
	return k.cells[minterm]
}

// Groups returns the prime implicants of the mapped function, which correspond
// exactly to the maximal rectangular groupings of adjacent true cells. The
// result is sorted deterministically.
func (k *KarnaughMap) Groups() []Implicant {
	var minterms []int
	for m, v := range k.cells {
		if v {
			minterms = append(minterms, m)
		}
	}
	return PrimeImplicants(minterms, nil, k.NumVars)
}

// String renders the map as a Gray-coded grid of 1s and 0s with binary row and
// column headers.
func (k *KarnaughMap) String() string {
	var b strings.Builder
	fmt.Fprintf(&b, "%*s", k.rowBits+2, "")
	for c := 0; c < k.Cols(); c++ {
		fmt.Fprintf(&b, " %s", BitString(GrayEncode(c), k.colBits))
	}
	b.WriteByte('\n')
	for r := 0; r < k.Rows(); r++ {
		fmt.Fprintf(&b, "%s |", BitString(GrayEncode(r), k.rowBits))
		for c := 0; c < k.Cols(); c++ {
			pad := k.colBits + 1
			fmt.Fprintf(&b, "%*s", pad, logicBit(k.At(r, c)))
		}
		b.WriteByte('\n')
	}
	return b.String()
}
