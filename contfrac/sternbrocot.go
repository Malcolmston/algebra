package contfrac

import (
	"math/big"
	"strings"
)

// Mediant returns the mediant (a+c)/(b+d) of the fractions a/b and c/d as an
// unreduced numerator and denominator. The mediant of two fractions lies
// strictly between them when they are distinct.
func Mediant(a, b, c, d int64) (int64, int64) {
	return a + c, b + d
}

// MediantFrac returns the mediant of two fractions as a [Frac] (not reduced,
// mirroring the Stern-Brocot construction).
func MediantFrac(f, g Frac) Frac {
	return Frac{f.Num + g.Num, f.Den + g.Den}
}

// MediantRat returns the mediant of two *big.Rat values.
func MediantRat(f, g *big.Rat) *big.Rat {
	num := new(big.Int).Add(f.Num(), g.Num())
	den := new(big.Int).Add(f.Denom(), g.Denom())
	return new(big.Rat).SetFrac(num, den)
}

// SternBrocotPath returns the path from the root 1/1 of the Stern-Brocot tree
// to the positive reduced fraction p/q, encoded as a string of 'L' (left) and
// 'R' (right) moves. The empty string denotes the root itself. It panics if
// p/q is not positive.
func SternBrocotPath(p, q int64) string {
	if p <= 0 || q <= 0 {
		panic("contfrac: SternBrocotPath requires a positive fraction")
	}
	p, q = ReduceFraction(p, q)
	var b strings.Builder
	for p != q {
		if p < q {
			b.WriteByte('L')
			q -= p
		} else {
			b.WriteByte('R')
			p -= q
		}
	}
	return b.String()
}

// SternBrocotFromPath returns the fraction reached by following the given
// L/R path from the root of the Stern-Brocot tree. Invalid characters are
// ignored. The empty path yields 1/1.
func SternBrocotFromPath(path string) Frac {
	// Boundaries: left = 0/1, right = 1/0 (infinity), current = mediant = 1/1.
	ln, ld := int64(0), int64(1)
	rn, rd := int64(1), int64(0)
	cn, cd := int64(1), int64(1)
	for i := 0; i < len(path); i++ {
		switch path[i] {
		case 'L', 'l':
			rn, rd = cn, cd
		case 'R', 'r':
			ln, ld = cn, cd
		default:
			continue
		}
		cn, cd = ln+rn, ld+rd
	}
	return Frac{cn, cd}
}

// SternBrocotParent returns the parent of p/q in the Stern-Brocot tree. The
// parent of the root 1/1 is itself.
func SternBrocotParent(p, q int64) Frac {
	path := SternBrocotPath(p, q)
	if path == "" {
		return Frac{1, 1}
	}
	return SternBrocotFromPath(path[:len(path)-1])
}

// SternBrocotLeftChild returns the left child of p/q in the Stern-Brocot tree.
func SternBrocotLeftChild(p, q int64) Frac {
	return SternBrocotFromPath(SternBrocotPath(p, q) + "L")
}

// SternBrocotRightChild returns the right child of p/q in the Stern-Brocot tree.
func SternBrocotRightChild(p, q int64) Frac {
	return SternBrocotFromPath(SternBrocotPath(p, q) + "R")
}

// SternBrocotChildren returns both children (left, right) of p/q.
func SternBrocotChildren(p, q int64) (left, right Frac) {
	return SternBrocotLeftChild(p, q), SternBrocotRightChild(p, q)
}

// SternBrocotDepth returns the depth of p/q in the Stern-Brocot tree (the root
// 1/1 has depth 0).
func SternBrocotDepth(p, q int64) int {
	return len(SternBrocotPath(p, q))
}

// SternBrocotAncestors returns the fractions on the path from the root down to
// (but not including) p/q, in order from the root.
func SternBrocotAncestors(p, q int64) []Frac {
	path := SternBrocotPath(p, q)
	out := make([]Frac, 0, len(path))
	for i := 0; i < len(path); i++ {
		out = append(out, SternBrocotFromPath(path[:i]))
	}
	return out
}

// CFToPath converts a continued fraction of a positive value into the
// corresponding Stern-Brocot L/R path. The moves alternate R, L, R, ... with
// run lengths given by the partial quotients, the final term being decremented
// by one (so [1; 2] -> "RL"). Both representations of the value (with the last
// term reduced or not) map to the same path.
func CFToPath(cf CF) string {
	if len(cf) == 0 {
		return ""
	}
	terms := cf.Clone()
	terms[len(terms)-1]-- // the node sits one mediant short of the infinite path
	var b strings.Builder
	letter := byte('R')
	for _, t := range terms {
		for j := int64(0); j < t; j++ {
			b.WriteByte(letter)
		}
		if letter == 'R' {
			letter = 'L'
		} else {
			letter = 'R'
		}
	}
	return b.String()
}

// PathToCF converts a Stern-Brocot L/R path into the continued fraction of the
// fraction it names. It is the inverse of [CFToPath]; the empty path yields the
// continued fraction [1].
func PathToCF(path string) CF {
	if path == "" {
		return CF{1}
	}
	// Run-length encode the path.
	type run struct {
		ch  byte
		cnt int64
	}
	var runs []run
	for i := 0; i < len(path); i++ {
		ch := path[i]
		if ch == 'l' {
			ch = 'L'
		} else if ch == 'r' {
			ch = 'R'
		}
		if len(runs) > 0 && runs[len(runs)-1].ch == ch {
			runs[len(runs)-1].cnt++
		} else {
			runs = append(runs, run{ch, 1})
		}
	}
	var cf CF
	expected := byte('R')
	idx := 0
	for idx < len(runs) {
		if runs[idx].ch == expected {
			cf = append(cf, runs[idx].cnt)
			idx++
		} else {
			cf = append(cf, 0)
		}
		if expected == 'R' {
			expected = 'L'
		} else {
			expected = 'R'
		}
	}
	cf[len(cf)-1]++ // undo the decrement applied in CFToPath
	return cf
}
