package projectivegeom

import "math/rand"

// RandomPoint returns a pseudo-random projective point whose homogeneous
// coordinates are drawn uniformly from [-1, 1). The result is deterministic for
// a given source, so callers control reproducibility by seeding r.
func RandomPoint(r *rand.Rand) Point {
	return Point{Vec3{2*r.Float64() - 1, 2*r.Float64() - 1, 2*r.Float64() - 1}}
}

// RandomFinitePoint returns a pseudo-random finite point [x y 1] with x, y drawn
// uniformly from [-1, 1).
func RandomFinitePoint(r *rand.Rand) Point {
	return PointFromAffine(2*r.Float64()-1, 2*r.Float64()-1)
}

// RandomLine returns a pseudo-random projective line with coordinates drawn
// uniformly from [-1, 1).
func RandomLine(r *rand.Rand) Line {
	return Line{Vec3{2*r.Float64() - 1, 2*r.Float64() - 1, 2*r.Float64() - 1}}
}

// RandomSPoint returns a pseudo-random point of RP^3 with homogeneous
// coordinates drawn uniformly from [-1, 1).
func RandomSPoint(r *rand.Rand) SPoint {
	return SPoint{Vec4{2*r.Float64() - 1, 2*r.Float64() - 1, 2*r.Float64() - 1, 2*r.Float64() - 1}}
}

// RandomHomography returns a pseudo-random homography whose matrix entries are
// drawn uniformly from [-1, 1). With probability one the matrix is invertible;
// the second result is false in the rare singular case.
func RandomHomography(r *rand.Rand) (Homography, bool) {
	var m Mat3
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			m[i][j] = 2*r.Float64() - 1
		}
	}
	h := Homography{m}
	return h, h.IsInvertible()
}

// PerturbPoint returns p with independent uniform noise of amplitude eps added
// to each homogeneous coordinate, useful for robustness testing.
func PerturbPoint(r *rand.Rand, p Point, eps float64) Point {
	return Point{Vec3{
		p.V.X + eps*(2*r.Float64()-1),
		p.V.Y + eps*(2*r.Float64()-1),
		p.V.Z + eps*(2*r.Float64()-1),
	}}
}
