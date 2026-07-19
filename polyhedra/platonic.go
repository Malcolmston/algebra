package polyhedra

import "math"

// Phi is the golden ratio (1 + sqrt 5) / 2, which governs the coordinates of the
// dodecahedron and icosahedron.
var Phi = (1 + math.Sqrt(5)) / 2

// PlatonicSolid captures the combinatorial and metric data of one of the five
// regular convex polyhedra at a given edge length.
type PlatonicSolid struct {
	Name        string  // common name
	Schlafli    [2]int  // Schläfli symbol {p, q}
	NumV        int     // number of vertices
	NumE        int     // number of edges
	NumF        int     // number of faces
	FaceSides   int     // p: sides per face
	VertexFaces int     // q: faces meeting at each vertex
	Edge        float64 // edge length a
}

// SurfaceArea returns the total surface area of the solid.
func (s PlatonicSolid) SurfaceArea() float64 {
	return float64(s.NumF) * RegularPolygonArea(s.FaceSides, s.Edge)
}

// EulerCharacteristic returns V - E + F, which is 2 for every Platonic solid.
func (s PlatonicSolid) EulerCharacteristic() int { return s.NumV - s.NumE + s.NumF }

// DihedralAngle returns the interior dihedral angle in radians between two
// adjacent faces, computed from the Schläfli symbol {p, q} as
// 2*arcsin(cos(π/q)/sin(π/p)).
func (s PlatonicSolid) DihedralAngle() float64 {
	p := float64(s.FaceSides)
	q := float64(s.VertexFaces)
	return 2 * math.Asin(math.Cos(math.Pi/q)/math.Sin(math.Pi/p))
}

// Circumradius returns the radius of the circumscribed sphere (through all
// vertices).
func (s PlatonicSolid) Circumradius() float64 {
	a := s.Edge
	switch s.Name {
	case "tetrahedron":
		return TetrahedronCircumradius(a)
	case "cube":
		return CubeCircumradius(a)
	case "octahedron":
		return OctahedronCircumradius(a)
	case "dodecahedron":
		return DodecahedronCircumradius(a)
	case "icosahedron":
		return IcosahedronCircumradius(a)
	}
	return 0
}

// Midradius returns the radius of the midsphere (tangent to every edge at its
// midpoint).
func (s PlatonicSolid) Midradius() float64 {
	a := s.Edge
	switch s.Name {
	case "tetrahedron":
		return a / (2 * math.Sqrt2)
	case "cube":
		return a * math.Sqrt2 / 2
	case "octahedron":
		return a / 2
	case "dodecahedron":
		return a * Phi * Phi / 2
	case "icosahedron":
		return a * Phi / 2
	}
	return 0
}

// Inradius returns the radius of the inscribed sphere (tangent to every face).
func (s PlatonicSolid) Inradius() float64 {
	a := s.Edge
	switch s.Name {
	case "tetrahedron":
		return a / (2 * math.Sqrt(6))
	case "cube":
		return a / 2
	case "octahedron":
		return a / math.Sqrt(6)
	case "dodecahedron":
		return a / 2 * math.Sqrt((25+11*math.Sqrt(5))/10)
	case "icosahedron":
		return a * Phi * Phi / (2 * math.Sqrt(3))
	}
	return 0
}

// Volume returns the enclosed volume of the solid.
func (s PlatonicSolid) Volume() float64 {
	a := s.Edge
	switch s.Name {
	case "tetrahedron":
		return a * a * a / (6 * math.Sqrt2)
	case "cube":
		return a * a * a
	case "octahedron":
		return math.Sqrt2 / 3 * a * a * a
	case "dodecahedron":
		return (15 + 7*math.Sqrt(5)) / 4 * a * a * a
	case "icosahedron":
		return 5 * (3 + math.Sqrt(5)) / 12 * a * a * a
	}
	return 0
}

// Mesh materialises the solid as a [Polyhedron] with exact vertex coordinates
// and outward-oriented polygonal faces at the solid's edge length, centered at
// the origin.
func (s PlatonicSolid) Mesh() *Polyhedron {
	switch s.Name {
	case "tetrahedron":
		return tetrahedronMesh(s.Edge)
	case "cube":
		return cubeMesh(s.Edge)
	case "octahedron":
		return octahedronMesh(s.Edge)
	case "dodecahedron":
		return dodecahedronMesh(s.Edge)
	case "icosahedron":
		return icosahedronMesh(s.Edge)
	}
	return nil
}

// NewTetrahedron returns the regular tetrahedron with the given edge length.
func NewTetrahedron(edge float64) PlatonicSolid {
	return PlatonicSolid{"tetrahedron", [2]int{3, 3}, 4, 6, 4, 3, 3, edge}
}

// NewCube returns the cube (regular hexahedron) with the given edge length.
func NewCube(edge float64) PlatonicSolid {
	return PlatonicSolid{"cube", [2]int{4, 3}, 8, 12, 6, 4, 3, edge}
}

// NewOctahedron returns the regular octahedron with the given edge length.
func NewOctahedron(edge float64) PlatonicSolid {
	return PlatonicSolid{"octahedron", [2]int{3, 4}, 6, 12, 8, 3, 4, edge}
}

// NewDodecahedron returns the regular dodecahedron with the given edge length.
func NewDodecahedron(edge float64) PlatonicSolid {
	return PlatonicSolid{"dodecahedron", [2]int{5, 3}, 20, 30, 12, 5, 3, edge}
}

// NewIcosahedron returns the regular icosahedron with the given edge length.
func NewIcosahedron(edge float64) PlatonicSolid {
	return PlatonicSolid{"icosahedron", [2]int{3, 5}, 12, 30, 20, 3, 5, edge}
}

// PlatonicSolids returns all five Platonic solids at the given edge length, in
// order of increasing face count.
func PlatonicSolids(edge float64) []PlatonicSolid {
	return []PlatonicSolid{
		NewTetrahedron(edge),
		NewCube(edge),
		NewOctahedron(edge),
		NewDodecahedron(edge),
		NewIcosahedron(edge),
	}
}

// ---- Standalone analytic property functions ----

// TetrahedronVolume returns the volume of a regular tetrahedron of edge a.
func TetrahedronVolume(a float64) float64 { return a * a * a / (6 * math.Sqrt2) }

// TetrahedronSurfaceArea returns the surface area of a regular tetrahedron of edge a.
func TetrahedronSurfaceArea(a float64) float64 { return math.Sqrt(3) * a * a }

// TetrahedronCircumradius returns the circumradius of a regular tetrahedron of edge a.
func TetrahedronCircumradius(a float64) float64 { return a * math.Sqrt(6) / 4 }

// TetrahedronInradius returns the inradius of a regular tetrahedron of edge a.
func TetrahedronInradius(a float64) float64 { return a / (2 * math.Sqrt(6)) }

// TetrahedronMidradius returns the midradius of a regular tetrahedron of edge a.
func TetrahedronMidradius(a float64) float64 { return a / (2 * math.Sqrt2) }

// TetrahedronDihedralAngle returns the dihedral angle arccos(1/3) of a regular tetrahedron.
func TetrahedronDihedralAngle() float64 { return math.Acos(1.0 / 3.0) }

// CubeVolume returns the volume of a cube of edge a.
func CubeVolume(a float64) float64 { return a * a * a }

// CubeSurfaceArea returns the surface area of a cube of edge a.
func CubeSurfaceArea(a float64) float64 { return 6 * a * a }

// CubeCircumradius returns the circumradius of a cube of edge a.
func CubeCircumradius(a float64) float64 { return a * math.Sqrt(3) / 2 }

// CubeInradius returns the inradius of a cube of edge a.
func CubeInradius(a float64) float64 { return a / 2 }

// CubeMidradius returns the midradius of a cube of edge a.
func CubeMidradius(a float64) float64 { return a * math.Sqrt2 / 2 }

// CubeDihedralAngle returns the dihedral angle (π/2) of a cube.
func CubeDihedralAngle() float64 { return math.Pi / 2 }

// OctahedronVolume returns the volume of a regular octahedron of edge a.
func OctahedronVolume(a float64) float64 { return math.Sqrt2 / 3 * a * a * a }

// OctahedronSurfaceArea returns the surface area of a regular octahedron of edge a.
func OctahedronSurfaceArea(a float64) float64 { return 2 * math.Sqrt(3) * a * a }

// OctahedronCircumradius returns the circumradius of a regular octahedron of edge a.
func OctahedronCircumradius(a float64) float64 { return a * math.Sqrt2 / 2 }

// OctahedronInradius returns the inradius of a regular octahedron of edge a.
func OctahedronInradius(a float64) float64 { return a / math.Sqrt(6) }

// OctahedronMidradius returns the midradius of a regular octahedron of edge a.
func OctahedronMidradius(a float64) float64 { return a / 2 }

// OctahedronDihedralAngle returns the dihedral angle arccos(-1/3) of a regular octahedron.
func OctahedronDihedralAngle() float64 { return math.Acos(-1.0 / 3.0) }

// DodecahedronVolume returns the volume of a regular dodecahedron of edge a.
func DodecahedronVolume(a float64) float64 { return (15 + 7*math.Sqrt(5)) / 4 * a * a * a }

// DodecahedronSurfaceArea returns the surface area of a regular dodecahedron of edge a.
func DodecahedronSurfaceArea(a float64) float64 { return 3 * math.Sqrt(25+10*math.Sqrt(5)) * a * a }

// DodecahedronCircumradius returns the circumradius of a regular dodecahedron of edge a.
func DodecahedronCircumradius(a float64) float64 { return a * math.Sqrt(3) / 4 * (1 + math.Sqrt(5)) }

// DodecahedronInradius returns the inradius of a regular dodecahedron of edge a.
func DodecahedronInradius(a float64) float64 { return a / 2 * math.Sqrt((25+11*math.Sqrt(5))/10) }

// DodecahedronMidradius returns the midradius of a regular dodecahedron of edge a.
func DodecahedronMidradius(a float64) float64 { return a * (3 + math.Sqrt(5)) / 4 }

// DodecahedronDihedralAngle returns the dihedral angle arccos(-1/sqrt 5) of a regular dodecahedron.
func DodecahedronDihedralAngle() float64 { return math.Acos(-1 / math.Sqrt(5)) }

// IcosahedronVolume returns the volume of a regular icosahedron of edge a.
func IcosahedronVolume(a float64) float64 { return 5 * (3 + math.Sqrt(5)) / 12 * a * a * a }

// IcosahedronSurfaceArea returns the surface area of a regular icosahedron of edge a.
func IcosahedronSurfaceArea(a float64) float64 { return 5 * math.Sqrt(3) * a * a }

// IcosahedronCircumradius returns the circumradius of a regular icosahedron of edge a.
func IcosahedronCircumradius(a float64) float64 { return a / 4 * math.Sqrt(10+2*math.Sqrt(5)) }

// IcosahedronInradius returns the inradius of a regular icosahedron of edge a.
func IcosahedronInradius(a float64) float64 { return a * Phi * Phi / (2 * math.Sqrt(3)) }

// IcosahedronMidradius returns the midradius of a regular icosahedron of edge a.
func IcosahedronMidradius(a float64) float64 { return a * Phi / 2 }

// IcosahedronDihedralAngle returns the dihedral angle arccos(-sqrt 5 / 3) of a regular icosahedron.
func IcosahedronDihedralAngle() float64 { return math.Acos(-math.Sqrt(5) / 3) }

// ---- Mesh builders ----

// scaleMesh scales every vertex of the mesh built at a reference edge so the
// resulting edge length equals target, then centers it at the origin.
func rescaleToEdge(verts []Vec3, faces [][]int, refEdge, target float64) *Polyhedron {
	f := target / refEdge
	vs := make([]Vec3, len(verts))
	for i, v := range verts {
		vs[i] = v.Scale(f)
	}
	p := &Polyhedron{Verts: vs, Faces: faces}
	return p.OrientOutward()
}

func tetrahedronMesh(a float64) *Polyhedron {
	verts := []Vec3{
		{1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1},
	}
	faces := [][]int{
		{0, 1, 2}, {0, 3, 1}, {0, 2, 3}, {1, 3, 2},
	}
	return rescaleToEdge(verts, faces, 2*math.Sqrt2, a)
}

func cubeMesh(a float64) *Polyhedron {
	verts := []Vec3{
		{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
		{-1, -1, 1}, {1, -1, 1}, {1, 1, 1}, {-1, 1, 1},
	}
	faces := [][]int{
		{0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4},
		{1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7},
	}
	return rescaleToEdge(verts, faces, 2, a)
}

func octahedronMesh(a float64) *Polyhedron {
	verts := []Vec3{
		{1, 0, 0}, {-1, 0, 0}, {0, 1, 0},
		{0, -1, 0}, {0, 0, 1}, {0, 0, -1},
	}
	faces := [][]int{
		{0, 2, 4}, {2, 1, 4}, {1, 3, 4}, {3, 0, 4},
		{2, 0, 5}, {1, 2, 5}, {3, 1, 5}, {0, 3, 5},
	}
	return rescaleToEdge(verts, faces, math.Sqrt2, a)
}

func dodecahedronMesh(a float64) *Polyhedron {
	ip := 1 / Phi
	verts := []Vec3{
		{1, 1, 1}, {1, 1, -1}, {1, -1, 1}, {1, -1, -1},
		{-1, 1, 1}, {-1, 1, -1}, {-1, -1, 1}, {-1, -1, -1},
		{0, ip, Phi}, {0, ip, -Phi}, {0, -ip, Phi}, {0, -ip, -Phi},
		{ip, Phi, 0}, {ip, -Phi, 0}, {-ip, Phi, 0}, {-ip, -Phi, 0},
		{Phi, 0, ip}, {Phi, 0, -ip}, {-Phi, 0, ip}, {-Phi, 0, -ip},
	}
	// Derive faces from the convex hull with coplanar merging so the pentagon
	// loops are correct by construction.
	tri, _ := ConvexHull(verts)
	merged := MergeCoplanarFaces(tri, 1e-6)
	// merged.Verts are the hull vertices in compacted order (all 20 are used).
	return rescaleMeshToEdge(merged, a)
}

func icosahedronMesh(a float64) *Polyhedron {
	verts := []Vec3{
		{0, 1, Phi}, {0, 1, -Phi}, {0, -1, Phi}, {0, -1, -Phi},
		{1, Phi, 0}, {1, -Phi, 0}, {-1, Phi, 0}, {-1, -Phi, 0},
		{Phi, 0, 1}, {Phi, 0, -1}, {-Phi, 0, 1}, {-Phi, 0, -1},
	}
	tri, _ := ConvexHull(verts)
	return rescaleMeshToEdge(tri, a)
}

// rescaleMeshToEdge scales an already-built mesh so its (assumed uniform) edge
// length equals target.
func rescaleMeshToEdge(p *Polyhedron, target float64) *Polyhedron {
	lengths := p.EdgeLengths()
	if len(lengths) == 0 {
		return p
	}
	ref := lengths[0]
	if ref == 0 {
		return p
	}
	return p.ScaledBy(target / ref).OrientOutward()
}
