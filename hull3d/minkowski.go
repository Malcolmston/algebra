package hull3d

import "errors"

// MinkowskiSumPoints returns the Minkowski sum of two point sets: every pairwise
// sum a+b for a in as and b in bs.
func MinkowskiSumPoints(as, bs []Vec3) []Vec3 {
	out := make([]Vec3, 0, len(as)*len(bs))
	for _, a := range as {
		for _, b := range bs {
			out = append(out, a.Add(b))
		}
	}
	return out
}

// MinkowskiDiffPoints returns the pointwise differences a-b for a in as and b in
// bs. Its convex hull is the Minkowski difference A ⊕ (-B), the configuration
// space obstacle whose containing the origin signals overlap of A and B.
func MinkowskiDiffPoints(as, bs []Vec3) []Vec3 {
	out := make([]Vec3, 0, len(as)*len(bs))
	for _, a := range as {
		for _, b := range bs {
			out = append(out, a.Sub(b))
		}
	}
	return out
}

// MinkowskiSum returns the convex polytope that is the Minkowski sum of two
// convex polytopes: the convex hull of all pairwise vertex sums. It returns an
// error if the result is degenerate (lower-dimensional).
func MinkowskiSum(a, b *Polytope) (*Polytope, error) {
	if a == nil || b == nil {
		return nil, errors.New("hull3d: nil polytope")
	}
	return ConvexHull(MinkowskiSumPoints(a.Vertices, b.Vertices))
}

// MinkowskiDifference returns the convex polytope A ⊕ (-B), the convex hull of
// all pairwise vertex differences. The two input bodies overlap exactly when
// this polytope contains the origin. It returns an error if the result is
// degenerate.
func MinkowskiDifference(a, b *Polytope) (*Polytope, error) {
	if a == nil || b == nil {
		return nil, errors.New("hull3d: nil polytope")
	}
	return ConvexHull(MinkowskiDiffPoints(a.Vertices, b.Vertices))
}

// PolytopesOverlap reports whether two convex polytopes overlap, by testing
// whether their Minkowski difference contains the origin.
func PolytopesOverlap(a, b *Polytope, eps float64) (bool, error) {
	diff, err := MinkowskiDifference(a, b)
	if err != nil {
		return false, err
	}
	return diff.Contains(Zero(), eps), nil
}

// Reflect returns a copy of the polytope reflected through the origin (each
// vertex negated), with face windings reversed so it stays outward-oriented.
func (p *Polytope) Reflect() *Polytope { return p.Scaled(-1) }

// ScaledSum returns the convex hull of s*a + t*b over vertices, a weighted
// Minkowski combination useful for morphing between convex bodies. It returns an
// error if the result is degenerate.
func ScaledSum(a, b *Polytope, s, t float64) (*Polytope, error) {
	if a == nil || b == nil {
		return nil, errors.New("hull3d: nil polytope")
	}
	pts := make([]Vec3, 0, len(a.Vertices)*len(b.Vertices))
	for _, va := range a.Vertices {
		for _, vb := range b.Vertices {
			pts = append(pts, va.Scale(s).Add(vb.Scale(t)))
		}
	}
	return ConvexHull(pts)
}
