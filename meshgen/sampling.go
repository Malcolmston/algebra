package meshgen

import "math/rand"

// RandomPoints returns n points drawn uniformly from the axis-aligned rectangle
// with the given corners, using a deterministic generator seeded by seed. The
// same seed always yields the same points.
func RandomPoints(n int, seed int64, min, max Vec2) []Vec2 {
	r := rand.New(rand.NewSource(seed))
	out := make([]Vec2, n)
	for i := range out {
		out[i] = Vec2{
			min.X + r.Float64()*(max.X-min.X),
			min.Y + r.Float64()*(max.Y-min.Y),
		}
	}
	return out
}

// JitteredGridPoints returns nx by ny points laid out on a regular lattice and
// perturbed by up to jitter times the cell size in each direction, using a
// deterministic generator seeded by seed.
func JitteredGridPoints(nx, ny int, x0, y0, dx, dy, jitter float64, seed int64) []Vec2 {
	r := rand.New(rand.NewSource(seed))
	out := make([]Vec2, 0, nx*ny)
	for j := 0; j < ny; j++ {
		for i := 0; i < nx; i++ {
			jx := (r.Float64()*2 - 1) * jitter * dx
			jy := (r.Float64()*2 - 1) * jitter * dy
			out = append(out, Vec2{x0 + float64(i)*dx + jx, y0 + float64(j)*dy + jy})
		}
	}
	return out
}

// RandomPointsInTriangle returns n points drawn uniformly from the triangle
// (a, b, c), using a deterministic generator seeded by seed.
func RandomPointsInTriangle(n int, seed int64, a, b, c Vec2) []Vec2 {
	r := rand.New(rand.NewSource(seed))
	out := make([]Vec2, n)
	for i := range out {
		u := r.Float64()
		v := r.Float64()
		if u+v > 1 {
			u, v = 1-u, 1-v
		}
		out[i] = a.Add(b.Sub(a).Scale(u)).Add(c.Sub(a).Scale(v))
	}
	return out
}
