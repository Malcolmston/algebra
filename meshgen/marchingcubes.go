package meshgen

// Triangle3 is a triangle in space, described by three vertex positions.
type Triangle3 struct {
	A, B, C Vec3
}

// Area returns the area of the triangle.
func (t Triangle3) Area() float64 { return TriangleArea3(t.A, t.B, t.C) }

// Normal returns the unit normal of the triangle following the right-hand rule.
func (t Triangle3) Normal() Vec3 { return TriangleNormal3(t.A, t.B, t.C) }

// Centroid returns the centroid of the triangle.
func (t Triangle3) Centroid() Vec3 { return t.A.Add(t.B).Add(t.C).Div(3) }

// ScalarGrid3 is a scalar field sampled on a regular Nx by Ny by Nz lattice.
// Node (i, j, k) sits at (X0+i*Dx, Y0+j*Dy, Z0+k*Dz) with value stored at
// Values[i + Nx*(j + Ny*k)].
type ScalarGrid3 struct {
	Nx, Ny, Nz             int
	X0, Y0, Z0, Dx, Dy, Dz float64
	Values                 []float64
}

// NewScalarGrid3 returns a zero-initialised 3-D grid. It panics if any
// dimension is less than two or any spacing is not positive.
func NewScalarGrid3(nx, ny, nz int, x0, y0, z0, dx, dy, dz float64) *ScalarGrid3 {
	if nx < 2 || ny < 2 || nz < 2 {
		panic("meshgen: grid needs at least 2x2x2 nodes")
	}
	if dx <= 0 || dy <= 0 || dz <= 0 {
		panic("meshgen: grid spacing must be positive")
	}
	return &ScalarGrid3{
		Nx: nx, Ny: ny, Nz: nz,
		X0: x0, Y0: y0, Z0: z0, Dx: dx, Dy: dy, Dz: dz,
		Values: make([]float64, nx*ny*nz),
	}
}

// SampleGrid3 builds a 3-D scalar grid by evaluating f at every node.
func SampleGrid3(nx, ny, nz int, x0, y0, z0, dx, dy, dz float64, f func(x, y, z float64) float64) *ScalarGrid3 {
	g := NewScalarGrid3(nx, ny, nz, x0, y0, z0, dx, dy, dz)
	for k := 0; k < nz; k++ {
		for j := 0; j < ny; j++ {
			for i := 0; i < nx; i++ {
				g.Set(i, j, k, f(
					g.X0+float64(i)*g.Dx,
					g.Y0+float64(j)*g.Dy,
					g.Z0+float64(k)*g.Dz,
				))
			}
		}
	}
	return g
}

// At returns the value at node (i, j, k).
func (g *ScalarGrid3) At(i, j, k int) float64 {
	return g.Values[i+g.Nx*(j+g.Ny*k)]
}

// Set assigns v to node (i, j, k).
func (g *ScalarGrid3) Set(i, j, k int, v float64) {
	g.Values[i+g.Nx*(j+g.Ny*k)] = v
}

// Node returns the position of node (i, j, k).
func (g *ScalarGrid3) Node(i, j, k int) Vec3 {
	return Vec3{
		g.X0 + float64(i)*g.Dx,
		g.Y0 + float64(j)*g.Dy,
		g.Z0 + float64(k)*g.Dz,
	}
}

// cubeOffsets are the eight corner offsets of a grid cell in (i,j,k) space.
var cubeOffsets = [8][3]int{
	{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
	{0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1},
}

// cubeTetra decomposes a cube into six tetrahedra sharing the main diagonal
// between corners 0 and 6, giving a watertight surface.
var cubeTetra = [6][4]int{
	{0, 1, 2, 6},
	{0, 2, 3, 6},
	{0, 3, 7, 6},
	{0, 7, 4, 6},
	{0, 4, 5, 6},
	{0, 5, 1, 6},
}

// MarchingCubes extracts the iso-surface of the grid at the given isovalue as a
// triangle soup. Each grid cell is decomposed into six tetrahedra and processed
// with marching tetrahedra, which yields a consistent, watertight surface with
// no lookup-table ambiguity.
func MarchingCubes(g *ScalarGrid3, iso float64) []Triangle3 {
	var out []Triangle3
	for k := 0; k < g.Nz-1; k++ {
		for j := 0; j < g.Ny-1; j++ {
			for i := 0; i < g.Nx-1; i++ {
				var cp [8]Vec3
				var cv [8]float64
				for c := 0; c < 8; c++ {
					o := cubeOffsets[c]
					cp[c] = g.Node(i+o[0], j+o[1], k+o[2])
					cv[c] = g.At(i+o[0], j+o[1], k+o[2])
				}
				for _, tet := range cubeTetra {
					var p [4]Vec3
					var v [4]float64
					for m := 0; m < 4; m++ {
						p[m] = cp[tet[m]]
						v[m] = cv[tet[m]]
					}
					out = append(out, MarchingTetrahedron(p, v, iso)...)
				}
			}
		}
	}
	return out
}

// MarchingTetrahedron returns the iso-surface triangles for a single
// tetrahedron whose corners are p with scalar values v. Corners with value
// greater than or equal to iso are treated as inside.
func MarchingTetrahedron(p [4]Vec3, v [4]float64, iso float64) []Triangle3 {
	var in, out []int
	for i := 0; i < 4; i++ {
		if v[i] >= iso {
			in = append(in, i)
		} else {
			out = append(out, i)
		}
	}
	interp := func(a, b int) Vec3 { return isoInterp3(p[a], p[b], v[a], v[b], iso) }
	switch len(in) {
	case 1:
		a := in[0]
		return []Triangle3{{interp(a, out[0]), interp(a, out[1]), interp(a, out[2])}}
	case 3:
		a := out[0]
		return []Triangle3{{interp(a, in[0]), interp(a, in[1]), interp(a, in[2])}}
	case 2:
		i0, i1 := in[0], in[1]
		o0, o1 := out[0], out[1]
		q0 := interp(i0, o0)
		q1 := interp(i0, o1)
		q2 := interp(i1, o1)
		q3 := interp(i1, o0)
		return []Triangle3{{q0, q1, q2}, {q0, q2, q3}}
	default:
		return nil
	}
}

// isoInterp3 returns the point on segment a-b where the linearly interpolated
// field equals iso, falling back to the midpoint for equal endpoint values.
func isoInterp3(a, b Vec3, va, vb, iso float64) Vec3 {
	den := vb - va
	if den == 0 {
		return a.Midpoint(b)
	}
	t := (iso - va) / den
	return a.Lerp(b, t)
}

// SurfaceArea returns the total area of the given triangles.
func SurfaceArea(tris []Triangle3) float64 {
	var s float64
	for _, t := range tris {
		s += t.Area()
	}
	return s
}
