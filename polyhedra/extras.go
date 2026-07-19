package polyhedra

import "math"

// SphereVolume returns the volume 4/3 π r³ of a ball of radius r.
func SphereVolume(r float64) float64 { return 4.0 / 3.0 * math.Pi * r * r * r }

// SphereSurfaceArea returns the surface area 4 π r² of a sphere of radius r.
func SphereSurfaceArea(r float64) float64 { return 4 * math.Pi * r * r }

// DihedralFromSchlafli returns the dihedral angle in radians of the regular
// polyhedron with Schläfli symbol {p, q}, namely 2 arcsin(cos(π/q)/sin(π/p)).
func DihedralFromSchlafli(p, q int) float64 {
	return 2 * math.Asin(math.Cos(math.Pi/float64(q))/math.Sin(math.Pi/float64(p)))
}

// RegularPolygonVertices returns the vertices of a regular n-gon of circumradius
// r lying in the plane z = height, starting at angle offset (radians) and going
// counter-clockwise. It returns nil for n < 3.
func RegularPolygonVertices(n int, r, height, offset float64) []Vec3 {
	if n < 3 {
		return nil
	}
	out := make([]Vec3, n)
	for i := 0; i < n; i++ {
		t := offset + 2*math.Pi*float64(i)/float64(n)
		out[i] = Vec3{r * math.Cos(t), r * math.Sin(t), height}
	}
	return out
}

// Prism returns a right uniform n-gonal prism whose regular n-gon base has edge
// length a and whose lateral height is h, centered on the z-axis. It returns nil
// for n < 3.
func Prism(n int, a, h float64) *Polyhedron {
	if n < 3 {
		return nil
	}
	r := RegularPolygonCircumradius(n, a)
	bottom := RegularPolygonVertices(n, r, -h/2, 0)
	top := RegularPolygonVertices(n, r, h/2, 0)
	verts := append(append([]Vec3{}, bottom...), top...)
	var faces [][]int
	// bottom face
	bl := make([]int, n)
	tl := make([]int, n)
	for i := 0; i < n; i++ {
		bl[i] = i
		tl[i] = n + i
	}
	faces = append(faces, bl, tl)
	for i := 0; i < n; i++ {
		j := (i + 1) % n
		faces = append(faces, []int{i, j, n + j, n + i})
	}
	p := &Polyhedron{Verts: verts, Faces: faces}
	return p.OrientOutward()
}

// PrismVolume returns the volume of a right prism with a regular n-gon base of
// edge a and height h.
func PrismVolume(n int, a, h float64) float64 { return RegularPolygonArea(n, a) * h }

// PrismSurfaceArea returns the surface area of a right prism with a regular
// n-gon base of edge a and height h (two bases plus n rectangles).
func PrismSurfaceArea(n int, a, h float64) float64 {
	return 2*RegularPolygonArea(n, a) + float64(n)*a*h
}

// UniformAntiprism returns the uniform n-gonal antiprism with edge length a: two
// parallel regular n-gons, the top rotated by π/n relative to the bottom, joined
// by 2n equilateral triangles, centered on the z-axis. It returns nil for n < 3.
func UniformAntiprism(n int, a float64) *Polyhedron {
	if n < 3 {
		return nil
	}
	r := RegularPolygonCircumradius(n, a)
	c := math.Cos(math.Pi / (2 * float64(n)))
	h := a * math.Sqrt(1-1/(4*c*c))
	bottom := RegularPolygonVertices(n, r, -h/2, 0)
	top := RegularPolygonVertices(n, r, h/2, math.Pi/float64(n))
	verts := append(append([]Vec3{}, bottom...), top...)
	bl := make([]int, n)
	tl := make([]int, n)
	for i := 0; i < n; i++ {
		bl[i] = i
		tl[i] = n + i
	}
	faces := [][]int{bl, tl}
	for i := 0; i < n; i++ {
		j := (i + 1) % n
		// Triangle pointing up (bottom edge to a top vertex) and down.
		faces = append(faces, []int{i, j, n + i})
		faces = append(faces, []int{j, n + j, n + i})
	}
	p := &Polyhedron{Verts: verts, Faces: faces}
	return p.OrientOutward()
}

// Pyramid returns a right pyramid over a regular n-gon base of edge a with apex
// at height h above the base center, centered on the z-axis. It returns nil for
// n < 3.
func Pyramid(n int, a, h float64) *Polyhedron {
	if n < 3 {
		return nil
	}
	r := RegularPolygonCircumradius(n, a)
	base := RegularPolygonVertices(n, r, 0, 0)
	verts := append([]Vec3{}, base...)
	apex := len(verts)
	verts = append(verts, Vec3{0, 0, h})
	bl := make([]int, n)
	for i := 0; i < n; i++ {
		bl[i] = i
	}
	faces := [][]int{bl}
	for i := 0; i < n; i++ {
		j := (i + 1) % n
		faces = append(faces, []int{i, j, apex})
	}
	p := &Polyhedron{Verts: verts, Faces: faces}
	return p.OrientOutward()
}

// PyramidVolume returns the volume of a pyramid with a regular n-gon base of
// edge a and apex height h, namely one third base area times height.
func PyramidVolume(n int, a, h float64) float64 { return RegularPolygonArea(n, a) * h / 3 }

// Bipyramid returns the n-gonal bipyramid: two right pyramids joined base to
// base, with apexes at heights ±h and a regular n-gon of edge a as the shared
// equator, centered on the z-axis. It returns nil for n < 3.
func Bipyramid(n int, a, h float64) *Polyhedron {
	if n < 3 {
		return nil
	}
	r := RegularPolygonCircumradius(n, a)
	ring := RegularPolygonVertices(n, r, 0, 0)
	verts := append([]Vec3{}, ring...)
	topApex := len(verts)
	verts = append(verts, Vec3{0, 0, h})
	botApex := len(verts)
	verts = append(verts, Vec3{0, 0, -h})
	var faces [][]int
	for i := 0; i < n; i++ {
		j := (i + 1) % n
		faces = append(faces, []int{i, j, topApex})
		faces = append(faces, []int{i, j, botApex})
	}
	p := &Polyhedron{Verts: verts, Faces: faces}
	return p.OrientOutward()
}

// FaceAdjacency returns, for each face index, the sorted indices of the faces
// that share at least one edge with it.
func FaceAdjacency(p *Polyhedron) [][]int {
	edgeFaces := make(map[Edge][]int)
	for fi, f := range p.Faces {
		m := len(f)
		for i := 0; i < m; i++ {
			e := MakeEdge(f[i], f[(i+1)%m])
			edgeFaces[e] = append(edgeFaces[e], fi)
		}
	}
	adjSet := make([]map[int]struct{}, len(p.Faces))
	for i := range adjSet {
		adjSet[i] = make(map[int]struct{})
	}
	for _, fs := range edgeFaces {
		for i := 0; i < len(fs); i++ {
			for j := i + 1; j < len(fs); j++ {
				adjSet[fs[i]][fs[j]] = struct{}{}
				adjSet[fs[j]][fs[i]] = struct{}{}
			}
		}
	}
	out := make([][]int, len(p.Faces))
	for i, s := range adjSet {
		lst := make([]int, 0, len(s))
		for k := range s {
			lst = append(lst, k)
		}
		sortInts(lst)
		out[i] = lst
	}
	return out
}

// sortInts sorts a slice of ints in place (small helper to avoid importing sort
// here; insertion sort is ample for face-degree-sized slices).
func sortInts(s []int) {
	for i := 1; i < len(s); i++ {
		for j := i; j > 0 && s[j-1] > s[j]; j-- {
			s[j-1], s[j] = s[j], s[j-1]
		}
	}
}

// DihedralAngleAtEdge returns the interior dihedral angle in radians between the
// two faces sharing edge e and true, or 0 and false when e is not shared by
// exactly two non-degenerate faces. The mesh should be outward-oriented (as
// produced by the solid constructors, ConvexHull, or OrientOutward).
func DihedralAngleAtEdge(p *Polyhedron, e Edge) (float64, bool) {
	var incident []int
	for fi, f := range p.Faces {
		m := len(f)
		found := false
		for i := 0; i < m; i++ {
			if MakeEdge(f[i], f[(i+1)%m]) == e {
				found = true
				break
			}
		}
		if found {
			incident = append(incident, fi)
		}
	}
	if len(incident) != 2 {
		return 0, false
	}
	n1, ok1 := p.FaceNormal(incident[0])
	n2, ok2 := p.FaceNormal(incident[1])
	if !ok1 || !ok2 {
		return 0, false
	}
	d := n1.Dot(n2)
	if d > 1 {
		d = 1
	} else if d < -1 {
		d = -1
	}
	// Angle between outward normals is π minus the interior dihedral angle.
	return math.Pi - math.Acos(d), true
}

// IsConvex reports whether the closed mesh is convex within tolerance tol: every
// vertex must lie on the inner side of (or within tol of) every face's outward
// supporting plane. The mesh must be outward-oriented.
func IsConvex(p *Polyhedron, tol float64) bool {
	for fi := range p.Faces {
		n, ok := p.FaceNormal(fi)
		if !ok {
			return false
		}
		base := p.Verts[p.Faces[fi][0]]
		for _, v := range p.Verts {
			if v.Sub(base).Dot(n) > tol {
				return false
			}
		}
	}
	return true
}

// MeanEdgeLength returns the average edge length of the mesh (0 when it has no
// edges).
func MeanEdgeLength(p *Polyhedron) float64 {
	ls := p.EdgeLengths()
	if len(ls) == 0 {
		return 0
	}
	var s float64
	for _, l := range ls {
		s += l
	}
	return s / float64(len(ls))
}

// Sphericity returns the isoperimetric sphericity of the mesh, the ratio of the
// surface area of a sphere of equal volume to the mesh's own surface area. It
// lies in (0, 1], reaching 1 only for a true sphere.
func Sphericity(p *Polyhedron) float64 {
	a := p.SurfaceArea()
	if a <= 0 {
		return 0
	}
	v := p.Volume()
	return math.Cbrt(math.Pi) * math.Cbrt(6*v*6*v) / a
}
