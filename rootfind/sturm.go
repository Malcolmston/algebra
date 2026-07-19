package rootfind

import (
	"math"
	"sort"
)

// SturmSequence returns the canonical Sturm chain of the polynomial p:
//
//	p0 = p, p1 = p', p_{k+1} = -(p_{k-1} mod p_k)
//
// continuing until a constant is reached. The sign-variation counts of this
// sequence, via [SturmVariations], give the exact number of distinct real roots
// of p in any interval by Sturm's theorem. The zero polynomial yields an empty
// sequence.
func SturmSequence(p Poly) []Poly {
	t := p.Trim().Clone()
	if t.Degree() < 0 {
		return nil
	}
	if t.Degree() == 0 {
		return []Poly{t}
	}
	seq := []Poly{t, t.Derivative()}
	for {
		last := seq[len(seq)-1]
		if last.Degree() <= 0 {
			break
		}
		prev := seq[len(seq)-2]
		_, r, err := prev.DivMod(last)
		if err != nil {
			break
		}
		neg := r.Neg()
		tol := 1e-12 * (1 + polyMaxAbs(neg))
		neg = polyChop(neg, tol)
		if neg.Degree() < 0 {
			break
		}
		seq = append(seq, neg)
	}
	return seq
}

// SturmVariations returns the number of sign changes in the Sturm sequence seq
// evaluated at x, ignoring terms that evaluate to zero. This is the function
// V(x) of Sturm's theorem.
func SturmVariations(seq []Poly, x float64) int {
	prev := 0
	count := 0
	for _, q := range seq {
		v := q.Eval(x)
		s := 0
		if v > 0 {
			s = 1
		} else if v < 0 {
			s = -1
		}
		if s == 0 {
			continue
		}
		if prev != 0 && s != prev {
			count++
		}
		prev = s
	}
	return count
}

// SturmVariationsAtNegInf returns the sign-variation count of the Sturm sequence
// as x -> -infinity, determined from the leading coefficients and degrees. It is
// used to count roots on unbounded left intervals.
func SturmVariationsAtNegInf(seq []Poly) int {
	prev := 0
	count := 0
	for _, q := range seq {
		d := q.Degree()
		if d < 0 {
			continue
		}
		s := 0
		lc := q[d]
		// sign of lc * (-inf)^d = lc * (-1)^d
		if d%2 == 0 {
			if lc > 0 {
				s = 1
			} else {
				s = -1
			}
		} else {
			if lc > 0 {
				s = -1
			} else {
				s = 1
			}
		}
		if prev != 0 && s != prev {
			count++
		}
		prev = s
	}
	return count
}

// SturmVariationsAtPosInf returns the sign-variation count of the Sturm sequence
// as x -> +infinity, determined from the signs of the leading coefficients.
func SturmVariationsAtPosInf(seq []Poly) int {
	prev := 0
	count := 0
	for _, q := range seq {
		d := q.Degree()
		if d < 0 {
			continue
		}
		s := 1
		if q[d] < 0 {
			s = -1
		}
		if prev != 0 && s != prev {
			count++
		}
		prev = s
	}
	return count
}

// SturmCountRoots returns the number of distinct real roots of the polynomial
// underlying the Sturm sequence seq that lie in the half-open interval (a, b],
// computed as V(a) - V(b). The endpoints must satisfy a < b.
func SturmCountRoots(seq []Poly, a, b float64) int {
	return SturmVariations(seq, a) - SturmVariations(seq, b)
}

// CountRealRoots returns the total number of distinct real roots of p (each
// counted once regardless of multiplicity) using Sturm's theorem over the whole
// real line.
func CountRealRoots(p Poly) int {
	seq := SturmSequence(p)
	if len(seq) == 0 {
		return 0
	}
	return SturmVariationsAtNegInf(seq) - SturmVariationsAtPosInf(seq)
}

// CountRealRootsInInterval returns the number of distinct real roots of p in the
// half-open interval (a, b].
func CountRealRootsInInterval(p Poly, a, b float64) int {
	return SturmCountRoots(SturmSequence(p), a, b)
}

// IsolateRoots returns a set of disjoint intervals, each containing exactly one
// distinct real root of p, covering all real roots. Bisection driven by Sturm
// root counts subdivides an interval known to contain every real root until each
// piece is isolating. The returned intervals are sorted in increasing order.
func IsolateRoots(p Poly) [][2]float64 {
	b := LagrangeBound(p)
	if b <= 0 {
		b = 1
	}
	// Widen slightly so all roots are strictly inside.
	b *= 1.0001
	return IsolateRootsInterval(p, -b, b)
}

// IsolateRootsInterval returns disjoint isolating intervals for the distinct
// real roots of p that lie in (a, b], using recursive Sturm bisection. Each
// returned interval contains exactly one root.
func IsolateRootsInterval(p Poly, a, b float64) [][2]float64 {
	seq := SturmSequence(p)
	if len(seq) == 0 {
		return nil
	}
	var out [][2]float64
	var rec func(lo, hi float64, depth int)
	rec = func(lo, hi float64, depth int) {
		n := SturmCountRoots(seq, lo, hi)
		if n == 0 {
			return
		}
		if n == 1 {
			out = append(out, [2]float64{lo, hi})
			return
		}
		if depth > 200 || hi-lo < 1e-14 {
			out = append(out, [2]float64{lo, hi})
			return
		}
		mid := 0.5 * (lo + hi)
		rec(lo, mid, depth+1)
		rec(mid, hi, depth+1)
	}
	rec(a, b, 0)
	sort.Slice(out, func(i, j int) bool { return out[i][0] < out[j][0] })
	return out
}

// SturmRefine refines an isolating interval [a, b] known to contain exactly one
// distinct real root of p down to width tol, using Sturm-count bisection. Unlike
// sign-based bisection it works for roots of even multiplicity, where p does not
// change sign. It returns the interval midpoint.
func SturmRefine(seq []Poly, a, b, tol float64) float64 {
	if tol <= 0 {
		tol = DefaultTol
	}
	for i := 0; i < 200 && b-a > tol; i++ {
		m := 0.5 * (a + b)
		if SturmCountRoots(seq, a, m) >= 1 {
			b = m
		} else {
			a = m
		}
	}
	return 0.5 * (a + b)
}

// SturmRealRoots returns all distinct real roots of p refined to tolerance tol.
// It isolates each root with [IsolateRoots] and refines it with Sturm-count
// bisection, so it is correct even in the presence of repeated (including
// even-multiplicity) roots. The roots are returned in increasing order.
func SturmRealRoots(p Poly, tol float64) []float64 {
	if tol <= 0 {
		tol = 1e-12
	}
	seq := SturmSequence(p)
	if len(seq) == 0 {
		return nil
	}
	intervals := IsolateRoots(p)
	roots := make([]float64, 0, len(intervals))
	for _, iv := range intervals {
		r := SturmRefine(seq, iv[0], iv[1], tol)
		// A Newton polish sharpens the estimate for simple roots.
		if v, d := p.EvalDeriv(r); d != 0 && math.Abs(v/d) < iv[1]-iv[0]+tol {
			r -= v / d
		}
		roots = append(roots, r)
	}
	sort.Float64s(roots)
	return roots
}
