package hull3d

import (
	"math"
	"math/rand"
)

// BoxPolytope returns the axis-aligned box centred at center with the given
// half-extents as a closed triangulated [Polytope].
func BoxPolytope(center, half Vec3) *Polytope {
	var vs []Vec3
	for _, sx := range []float64{-1, 1} {
		for _, sy := range []float64{-1, 1} {
			for _, sz := range []float64{-1, 1} {
				vs = append(vs, center.Add(Vec3{sx * half.X, sy * half.Y, sz * half.Z}))
			}
		}
	}
	p, _ := ConvexHull(vs)
	return p
}

// UnitCube returns the axis-aligned cube of side 1 centred at the origin.
func UnitCube() *Polytope { return BoxPolytope(Zero(), Splat(0.5)) }

// TetrahedronPolytope returns the tetrahedron with the four given corners as a
// closed [Polytope]. It returns an error if the corners are coplanar.
func TetrahedronPolytope(a, b, c, d Vec3) (*Polytope, error) {
	return ConvexHull([]Vec3{a, b, c, d})
}

// RegularTetrahedron returns a regular tetrahedron with the given circumradius
// centred at the origin.
func RegularTetrahedron(r float64) *Polytope {
	s := r / math.Sqrt(3)
	vs := []Vec3{
		{s, s, s}, {s, -s, -s}, {-s, s, -s}, {-s, -s, s},
	}
	p, _ := ConvexHull(vs)
	return p
}

// Octahedron returns the regular octahedron with the given circumradius centred
// at the origin (vertices at ±r along each axis).
func Octahedron(r float64) *Polytope {
	vs := []Vec3{
		{r, 0, 0}, {-r, 0, 0}, {0, r, 0}, {0, -r, 0}, {0, 0, r}, {0, 0, -r},
	}
	p, _ := ConvexHull(vs)
	return p
}

// Icosahedron returns the regular icosahedron with the given circumradius
// centred at the origin.
func Icosahedron(r float64) *Polytope {
	phi := (1 + math.Sqrt(5)) / 2
	raw := []Vec3{
		{-1, phi, 0}, {1, phi, 0}, {-1, -phi, 0}, {1, -phi, 0},
		{0, -1, phi}, {0, 1, phi}, {0, -1, -phi}, {0, 1, -phi},
		{phi, 0, -1}, {phi, 0, 1}, {-phi, 0, -1}, {-phi, 0, 1},
	}
	scale := r / math.Sqrt(1+phi*phi)
	vs := make([]Vec3, len(raw))
	for i, v := range raw {
		vs[i] = v.Scale(scale)
	}
	p, _ := ConvexHull(vs)
	return p
}

// Dodecahedron returns the regular dodecahedron with the given circumradius
// centred at the origin.
func Dodecahedron(r float64) *Polytope {
	phi := (1 + math.Sqrt(5)) / 2
	inv := 1 / phi
	var raw []Vec3
	for _, sx := range []float64{-1, 1} {
		for _, sy := range []float64{-1, 1} {
			for _, sz := range []float64{-1, 1} {
				raw = append(raw, Vec3{sx, sy, sz})
			}
		}
	}
	for _, s1 := range []float64{-1, 1} {
		for _, s2 := range []float64{-1, 1} {
			raw = append(raw, Vec3{0, s1 * inv, s2 * phi})
			raw = append(raw, Vec3{s1 * inv, s2 * phi, 0})
			raw = append(raw, Vec3{s1 * phi, 0, s2 * inv})
		}
	}
	scale := r / math.Sqrt(3)
	vs := make([]Vec3, len(raw))
	for i, v := range raw {
		vs[i] = v.Scale(scale)
	}
	p, _ := ConvexHull(vs)
	return p
}

// SamplePointsInBall returns n points sampled uniformly inside the ball of the
// given radius centred at center, using the supplied seed for reproducibility.
func SamplePointsInBall(n int, center Vec3, radius float64, seed int64) []Vec3 {
	r := rand.New(rand.NewSource(seed))
	out := make([]Vec3, 0, n)
	for len(out) < n {
		p := Vec3{r.Float64()*2 - 1, r.Float64()*2 - 1, r.Float64()*2 - 1}
		if p.LengthSq() <= 1 {
			out = append(out, center.Add(p.Scale(radius)))
		}
	}
	return out
}

// SamplePointsOnSphere returns n points sampled uniformly on the sphere of the
// given radius centred at center, using the supplied seed.
func SamplePointsOnSphere(n int, center Vec3, radius float64, seed int64) []Vec3 {
	r := rand.New(rand.NewSource(seed))
	out := make([]Vec3, 0, n)
	for len(out) < n {
		p := Vec3{r.NormFloat64(), r.NormFloat64(), r.NormFloat64()}
		if u, err := p.Normalize(); err == nil {
			out = append(out, center.Add(u.Scale(radius)))
		}
	}
	return out
}

// SphereVolume returns the volume of a ball of the given radius.
func SphereVolume(radius float64) float64 { return 4.0 / 3.0 * math.Pi * radius * radius * radius }

// SphereSurfaceArea returns the surface area of a sphere of the given radius.
func SphereSurfaceArea(radius float64) float64 { return 4 * math.Pi * radius * radius }
