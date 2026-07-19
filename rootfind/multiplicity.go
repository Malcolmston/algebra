package rootfind

import (
	"math"
	"math/cmplx"
	"sort"
)

// SquareFree returns the squarefree part of p, that is p divided by
// gcd(p, p'). The result has the same roots as p but each with multiplicity one.
// It returns ErrZeroPolynomial for the zero polynomial and a constant for a
// nonzero constant input.
func SquareFree(p Poly) (Poly, error) {
	if p.Degree() < 0 {
		return nil, ErrZeroPolynomial
	}
	if p.Degree() == 0 {
		return Poly{1}, nil
	}
	g := p.GCD(p.Derivative())
	if g.Degree() <= 0 {
		m, _ := p.Monic()
		return m, nil
	}
	q, _, err := p.DivMod(g)
	if err != nil {
		return nil, err
	}
	m, err := q.Monic()
	if err != nil {
		return q, nil
	}
	return m, nil
}

// SquareFreeFactor is one factor of a squarefree factorization: the polynomial
// Factor collects exactly the roots that occur with the given Multiplicity.
type SquareFreeFactor struct {
	// Factor is a monic squarefree polynomial whose roots all have the same
	// multiplicity in the original polynomial.
	Factor Poly
	// Multiplicity is the common multiplicity of the roots of Factor.
	Multiplicity int
}

// SquareFreeFactorization returns the squarefree factorization of p using Yun's
// algorithm: it decomposes p (up to its leading constant) as a product
// prod_i a_i(x)^i where each a_i is monic and squarefree and the a_i are
// pairwise coprime. Only factors of positive degree are returned, sorted by
// increasing multiplicity. This groups the roots of p by their multiplicity
// exactly, using only gcd computations.
func SquareFreeFactorization(p Poly) ([]SquareFreeFactor, error) {
	if p.Degree() < 0 {
		return nil, ErrZeroPolynomial
	}
	if p.Degree() == 0 {
		return nil, nil
	}
	pm, _ := p.Monic()
	a0 := pm.GCD(pm.Derivative())
	b, _, err := pm.DivMod(a0)
	if err != nil {
		return nil, err
	}
	c, _, err := pm.Derivative().DivMod(a0)
	if err != nil {
		return nil, err
	}
	d := c.Sub(b.Derivative())
	var out []SquareFreeFactor
	i := 1
	for b.Degree() > 0 {
		a := b.GCD(d)
		if a.Degree() > 0 {
			m, _ := a.Monic()
			out = append(out, SquareFreeFactor{Factor: m, Multiplicity: i})
		}
		nb, _, err := b.DivMod(a)
		if err != nil {
			break
		}
		nc, _, err := d.DivMod(a)
		if err != nil {
			break
		}
		b = nb
		d = nc.Sub(b.Derivative())
		i++
		if i > p.Degree()+1 {
			break
		}
	}
	return out, nil
}

// Multiplicity returns the multiplicity of the real value r as a root of p,
// i.e. the largest k with (x-r)^k dividing p. It uses the derivative test:
// r has multiplicity k when p, p', ... , p^(k-1) all vanish at r but p^(k) does
// not. Values are compared against a scaled tolerance tol. A return of 0 means r
// is not a root.
func Multiplicity(p Poly, r, tol float64) int {
	if tol <= 0 {
		tol = 1e-8
	}
	d := p.Degree()
	if d < 0 {
		return 0
	}
	scale := 1 + polyMaxAbs(p) + math.Abs(r)
	cur := p.Trim().Clone()
	k := 0
	for cur.Degree() >= 0 {
		v := cur.Eval(r)
		if math.Abs(v) > tol*scale {
			break
		}
		k++
		if cur.Degree() == 0 {
			break
		}
		cur = cur.Derivative()
	}
	return k
}

// MultiplicityComplex returns the multiplicity of the complex value r as a root
// of the complex polynomial c, via the complex derivative test with tolerance
// tol.
func MultiplicityComplex(c CPoly, r complex128, tol float64) int {
	if tol <= 0 {
		tol = 1e-8
	}
	d := c.Degree()
	if d < 0 {
		return 0
	}
	scale := 1 + cpolyMaxAbs(c) + cmplx.Abs(r)
	cur := c.Trim().Clone()
	k := 0
	for cur.Degree() >= 0 {
		v := cur.Eval(r)
		if cmplx.Abs(v) > tol*scale {
			break
		}
		k++
		if cur.Degree() == 0 {
			break
		}
		cur = cur.Derivative()
	}
	return k
}

// DistinctRealRoots returns the distinct real roots of p, each listed once,
// refined to tolerance tol. It is [SturmRealRoots] under a descriptive name and
// discards multiplicity information.
func DistinctRealRoots(p Poly, tol float64) []float64 {
	return SturmRealRoots(p, tol)
}

// RootMultiplicity pairs a real root value with its multiplicity in p.
type RootMultiplicity struct {
	// Root is the location of the root.
	Root float64
	// Multiplicity is how many times Root occurs as a root of p.
	Multiplicity int
}

// RealRootsWithMultiplicity returns the distinct real roots of p together with
// their multiplicities. Distinct roots are located by [SturmRealRoots] and each
// multiplicity is measured by the derivative test in [Multiplicity]. The result
// is sorted by root location.
func RealRootsWithMultiplicity(p Poly, tol float64) []RootMultiplicity {
	roots := SturmRealRoots(p, tol)
	mtol := 1e-6
	out := make([]RootMultiplicity, 0, len(roots))
	for _, r := range roots {
		m := Multiplicity(p, r, mtol)
		if m == 0 {
			m = 1
		}
		out = append(out, RootMultiplicity{Root: r, Multiplicity: m})
	}
	sort.Slice(out, func(i, j int) bool { return out[i].Root < out[j].Root })
	return out
}

// DeflateRoots divides p successively by (x - r) for each r in roots, returning
// the final quotient. It is the composition of repeated [Poly.DeflateReal] and
// removes the listed roots from p.
func DeflateRoots(p Poly, roots []float64) Poly {
	q := p.Trim().Clone()
	for _, r := range roots {
		if q.Degree() < 1 {
			break
		}
		q, _ = q.DeflateReal(r)
	}
	return q
}

// GroupComplexRoots clusters a list of complex roots so that roots within
// distance tol of one another are merged into a single representative (their
// centroid) carrying a multiplicity equal to the cluster size. This recovers
// multiplicity from the output of a numerical global solver, where a root of
// multiplicity m appears as a tight cluster of m nearby approximate roots.
func GroupComplexRoots(roots []complex128, tol float64) []ComplexRootMultiplicity {
	if tol <= 0 {
		tol = 1e-6
	}
	used := make([]bool, len(roots))
	var out []ComplexRootMultiplicity
	for i := range roots {
		if used[i] {
			continue
		}
		sum := roots[i]
		count := 1
		used[i] = true
		for j := i + 1; j < len(roots); j++ {
			if used[j] {
				continue
			}
			if cmplx.Abs(roots[j]-roots[i]) <= tol {
				sum += roots[j]
				count++
				used[j] = true
			}
		}
		out = append(out, ComplexRootMultiplicity{
			Root:         sum / complex(float64(count), 0),
			Multiplicity: count,
		})
	}
	sortComplexMult(out)
	return out
}

// ComplexRootMultiplicity pairs a complex root with its multiplicity.
type ComplexRootMultiplicity struct {
	// Root is the location of the root in the complex plane.
	Root complex128
	// Multiplicity is how many times Root occurs.
	Multiplicity int
}

// sortComplexMult orders complex-root/multiplicity pairs deterministically.
func sortComplexMult(s []ComplexRootMultiplicity) {
	sort.Slice(s, func(i, j int) bool {
		if real(s[i].Root) != real(s[j].Root) {
			return real(s[i].Root) < real(s[j].Root)
		}
		return imag(s[i].Root) < imag(s[j].Root)
	})
}
