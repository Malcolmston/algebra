# polyhedra — provides tools for three-dimensional convex geometry with a focus on polyhedra

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/polyhedra.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/polyhedra)

Package polyhedra provides tools for three-dimensional convex geometry with a
focus on polyhedra: the five Platonic solids, the thirteen Archimedean solids,
general triangulated and polygonal meshes, and 3-D convex hulls.

The package is organised around a small number of value types. `Vec3` is the
core three-component vector, used interchangeably as a point or a
displacement. `Polyhedron` is a boundary representation (a set of vertices
together with polygonal faces given as ordered index loops) on top of which
the package computes surface area, enclosed volume, centroid, face normals and
areas, edge enumeration, the Euler characteristic, and the (polar/canonical)
dual polyhedron. `PlatonicSolid` and `ArchimedeanSolid` capture the
combinatorial and metric data of the regular and semiregular convex solids and
can materialise them as `Polyhedron` meshes with exact vertex coordinates.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/polyhedra"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func CubeCircumradius(a float64) float64` | CubeCircumradius returns the circumradius of a cube of edge a. |
| `func CubeDihedralAngle() float64` | CubeDihedralAngle returns the dihedral angle (π/2) of a cube. |
| `func CubeInradius(a float64) float64` | CubeInradius returns the inradius of a cube of edge a. |
| `func CubeMidradius(a float64) float64` | CubeMidradius returns the midradius of a cube of edge a. |
| `func CubeSurfaceArea(a float64) float64` | CubeSurfaceArea returns the surface area of a cube of edge a. |
| `func CubeVolume(a float64) float64` | CubeVolume returns the volume of a cube of edge a. |
| `func DihedralAngleAtEdge(p *Polyhedron, e Edge) (float64, bool)` | DihedralAngleAtEdge returns the interior dihedral angle in radians between the two faces sharing edge e and true, or 0 and false when e is not shared… |
| `func DihedralFromSchlafli(p, q int) float64` | DihedralFromSchlafli returns the dihedral angle in radians of the regular polyhedron with Schläfli symbol {p, q}, namely 2… |
| `func DodecahedronCircumradius(a float64) float64` | DodecahedronCircumradius returns the circumradius of a regular dodecahedron of edge a. |
| `func DodecahedronDihedralAngle() float64` | DodecahedronDihedralAngle returns the dihedral angle arccos(-1/sqrt 5) of a regular dodecahedron. |
| `func DodecahedronInradius(a float64) float64` | DodecahedronInradius returns the inradius of a regular dodecahedron of edge a. |
| `func DodecahedronMidradius(a float64) float64` | DodecahedronMidradius returns the midradius of a regular dodecahedron of edge a. |
| `func DodecahedronSurfaceArea(a float64) float64` | DodecahedronSurfaceArea returns the surface area of a regular dodecahedron of edge a. |
| `func DodecahedronVolume(a float64) float64` | DodecahedronVolume returns the volume of a regular dodecahedron of edge a. |
| `func FaceAdjacency(p *Polyhedron) [][]int` | FaceAdjacency returns, for each face index, the sorted indices of the faces that share at least one edge with it. |
| `func IcosahedronCircumradius(a float64) float64` | IcosahedronCircumradius returns the circumradius of a regular icosahedron of edge a. |
| `func IcosahedronDihedralAngle() float64` | IcosahedronDihedralAngle returns the dihedral angle arccos(-sqrt 5 / 3) of a regular icosahedron. |
| `func IcosahedronInradius(a float64) float64` | IcosahedronInradius returns the inradius of a regular icosahedron of edge a. |
| `func IcosahedronMidradius(a float64) float64` | IcosahedronMidradius returns the midradius of a regular icosahedron of edge a. |
| `func IcosahedronSurfaceArea(a float64) float64` | IcosahedronSurfaceArea returns the surface area of a regular icosahedron of edge a. |
| `func IcosahedronVolume(a float64) float64` | IcosahedronVolume returns the volume of a regular icosahedron of edge a. |
| `func IsConvex(p *Polyhedron, tol float64) bool` | IsConvex reports whether the closed mesh is convex within tolerance tol: every vertex must lie on the inner side of (or within tol of) every face's… |
| `func IsPlanar(loop []Vec3, tol float64) bool` | IsPlanar reports whether every vertex of the loop lies within tol of the plane through the first vertex with the polygon's Newell normal. |
| `func MeanEdgeLength(p *Polyhedron) float64` | MeanEdgeLength returns the average edge length of the mesh (0 when it has no edges). |
| `func OctahedronCircumradius(a float64) float64` | OctahedronCircumradius returns the circumradius of a regular octahedron of edge a. |
| `func OctahedronDihedralAngle() float64` | OctahedronDihedralAngle returns the dihedral angle arccos(-1/3) of a regular octahedron. |
| `func OctahedronInradius(a float64) float64` | OctahedronInradius returns the inradius of a regular octahedron of edge a. |
| `func OctahedronMidradius(a float64) float64` | OctahedronMidradius returns the midradius of a regular octahedron of edge a. |
| `func OctahedronSurfaceArea(a float64) float64` | OctahedronSurfaceArea returns the surface area of a regular octahedron of edge a. |
| `func OctahedronVolume(a float64) float64` | OctahedronVolume returns the volume of a regular octahedron of edge a. |
| `func PolygonArea(loop []Vec3) float64` | PolygonArea returns the area of the planar polygon given by the ordered loop of vertices, using Newell's normal so the result is exact for planar… |
| `func PolygonPerimeter(loop []Vec3) float64` | PolygonPerimeter returns the total edge length of the closed polygon loop. |
| `func PrismSurfaceArea(n int, a, h float64) float64` | PrismSurfaceArea returns the surface area of a right prism with a regular n-gon base of edge a and height h (two bases plus n rectangles). |
| `func PrismVolume(n int, a, h float64) float64` | PrismVolume returns the volume of a right prism with a regular n-gon base of edge a and height h. |
| `func PyramidVolume(n int, a, h float64) float64` | PyramidVolume returns the volume of a pyramid with a regular n-gon base of edge a and apex height h, namely one third base area times height. |
| `func RegularPolygonArea(n int, a float64) float64` | RegularPolygonArea returns the area of a regular n-gon with side length a. |
| `func RegularPolygonCircumradius(n int, a float64) float64` | RegularPolygonCircumradius returns the circumradius of a regular n-gon with side length a (the distance from center to a vertex). |
| `func RegularPolygonInradius(n int, a float64) float64` | RegularPolygonInradius returns the inradius (apothem) of a regular n-gon with side length a (the distance from center to an edge midpoint). |
| `func RegularPolygonInteriorAngle(n int) float64` | RegularPolygonInteriorAngle returns the interior angle in radians of a regular n-gon. |
| `func SignedTetraVolume(a, b, c, d Vec3) float64` | SignedTetraVolume returns the signed volume of the tetrahedron with vertices a, b, c and d, equal to one sixth of (b-a) · ((c-a) × (d-a)). |
| `func SphereSurfaceArea(r float64) float64` | SphereSurfaceArea returns the surface area 4 π r² of a sphere of radius r. |
| `func SphereVolume(r float64) float64` | SphereVolume returns the volume 4/3 π r³ of a ball of radius r. |
| `func Sphericity(p *Polyhedron) float64` | Sphericity returns the isoperimetric sphericity of the mesh, the ratio of the surface area of a sphere of equal volume to the mesh's own surface area. |
| `func TetraVolume(a, b, c, d Vec3) float64` | TetraVolume returns the (non-negative) volume of the tetrahedron a, b, c, d. |
| `func TetrahedronCircumradius(a float64) float64` | TetrahedronCircumradius returns the circumradius of a regular tetrahedron of edge a. |
| `func TetrahedronDihedralAngle() float64` | TetrahedronDihedralAngle returns the dihedral angle arccos(1/3) of a regular tetrahedron. |
| `func TetrahedronInradius(a float64) float64` | TetrahedronInradius returns the inradius of a regular tetrahedron of edge a. |
| `func TetrahedronMidradius(a float64) float64` | TetrahedronMidradius returns the midradius of a regular tetrahedron of edge a. |
| `func TetrahedronSurfaceArea(a float64) float64` | TetrahedronSurfaceArea returns the surface area of a regular tetrahedron of edge a. |
| `func TetrahedronVolume(a float64) float64` | TetrahedronVolume returns the volume of a regular tetrahedron of edge a. |
| `func TriangleArea(a, b, c Vec3) float64` | TriangleArea returns the area of the triangle with vertices a, b and c. |

### Types

| Type | What it is |
| --- | --- |
| `ArchimedeanSolid` | ArchimedeanSolid captures the combinatorial and metric data of one of the thirteen semiregular convex polyhedra (uniform, with regular polygon faces… |
| `Edge` | Edge is an undirected edge between two vertex indices. |
| `FaceType` | FaceType records that a solid has Count regular faces each with Sides sides. |
| `PlatonicSolid` | PlatonicSolid captures the combinatorial and metric data of one of the five regular convex polyhedra at a given edge length. |
| `Polyhedron` | Polyhedron is a boundary representation of a polyhedron: a list of vertex positions together with faces, each face being an ordered loop of vertex… |
| `Vec3` | Vec3 is a vector in three-dimensional Euclidean space. |

<details>
<summary><code>ArchimedeanSolid</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ArchimedeanSolids() []ArchimedeanSolid` | ArchimedeanSolids returns all thirteen Archimedean solids in a canonical order (increasing face count, then vertices). |
| `func NewCuboctahedron() ArchimedeanSolid` | NewCuboctahedron returns the cuboctahedron (vertex figure 3.4.3.4). |
| `func NewIcosidodecahedron() ArchimedeanSolid` | NewIcosidodecahedron returns the icosidodecahedron (vertex figure 3.5.3.5). |
| `func NewRhombicosidodecahedron() ArchimedeanSolid` | NewRhombicosidodecahedron returns the rhombicosidodecahedron (vertex figure 3.4.5.4). |
| `func NewRhombicuboctahedron() ArchimedeanSolid` | NewRhombicuboctahedron returns the rhombicuboctahedron (vertex figure 3.4.4.4). |
| `func NewSnubCube() ArchimedeanSolid` | NewSnubCube returns the snub cube (vertex figure 3.3.3.3.4). |
| `func NewSnubDodecahedron() ArchimedeanSolid` | NewSnubDodecahedron returns the snub dodecahedron (vertex figure 3.3.3.3.5). |
| `func NewTruncatedCube() ArchimedeanSolid` | NewTruncatedCube returns the truncated cube (vertex figure 3.8.8). |
| `func NewTruncatedCuboctahedron() ArchimedeanSolid` | NewTruncatedCuboctahedron returns the truncated cuboctahedron, also called the great rhombicuboctahedron (vertex figure 4.6.8). |
| `func NewTruncatedDodecahedron() ArchimedeanSolid` | NewTruncatedDodecahedron returns the truncated dodecahedron (vertex figure 3.10.10). |
| `func NewTruncatedIcosahedron() ArchimedeanSolid` | NewTruncatedIcosahedron returns the truncated icosahedron, the classic football/soccer-ball shape (vertex figure 5.6.6). |
| `func NewTruncatedIcosidodecahedron() ArchimedeanSolid` | NewTruncatedIcosidodecahedron returns the truncated icosidodecahedron, also called the great rhombicosidodecahedron (vertex figure 4.6.10). |
| `func NewTruncatedOctahedron() ArchimedeanSolid` | NewTruncatedOctahedron returns the truncated octahedron (vertex figure 4.6.6). |
| `func NewTruncatedTetrahedron() ArchimedeanSolid` | NewTruncatedTetrahedron returns the truncated tetrahedron (vertex figure 3.6.6). |
| `func (s ArchimedeanSolid) Circumradius(a float64) float64` | Circumradius returns the circumradius (center-to-vertex distance) at edge length a. |
| `func (s ArchimedeanSolid) EulerCharacteristic() int` | EulerCharacteristic returns V - E + F, which is 2 for every Archimedean solid. |
| `func (s ArchimedeanSolid) Midradius(a float64) float64` | Midradius returns the midradius at edge length a, derived from the circumradius and edge via ρ = sqrt(R² - a²/4). |
| `func (s ArchimedeanSolid) SurfaceArea(a float64) float64` | SurfaceArea returns the surface area of the solid at edge length a, summed exactly over its regular polygonal faces. |
| `func (s ArchimedeanSolid) TotalFaceCount() int` | TotalFaceCount returns the number of faces implied by the Faces breakdown; it equals NumF for a consistent record and is useful as a self-check. |
| `func (s ArchimedeanSolid) Volume(a float64) float64` | Volume returns the volume of the solid at edge length a. |

</details>

<details>
<summary><code>Edge</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func MakeEdge(i, j int) Edge` | MakeEdge returns the canonical Edge connecting vertex indices i and j, with endpoints sorted so that A <= B. |

</details>

<details>
<summary><code>PlatonicSolid</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewCube(edge float64) PlatonicSolid` | NewCube returns the cube (regular hexahedron) with the given edge length. |
| `func NewDodecahedron(edge float64) PlatonicSolid` | NewDodecahedron returns the regular dodecahedron with the given edge length. |
| `func NewIcosahedron(edge float64) PlatonicSolid` | NewIcosahedron returns the regular icosahedron with the given edge length. |
| `func NewOctahedron(edge float64) PlatonicSolid` | NewOctahedron returns the regular octahedron with the given edge length. |
| `func NewTetrahedron(edge float64) PlatonicSolid` | NewTetrahedron returns the regular tetrahedron with the given edge length. |
| `func PlatonicSolids(edge float64) []PlatonicSolid` | PlatonicSolids returns all five Platonic solids at the given edge length, in order of increasing face count. |
| `func (s PlatonicSolid) Circumradius() float64` | Circumradius returns the radius of the circumscribed sphere (through all vertices). |
| `func (s PlatonicSolid) DihedralAngle() float64` | DihedralAngle returns the interior dihedral angle in radians between two adjacent faces, computed from the Schläfli symbol {p, q} as… |
| `func (s PlatonicSolid) EulerCharacteristic() int` | EulerCharacteristic returns V - E + F, which is 2 for every Platonic solid. |
| `func (s PlatonicSolid) Inradius() float64` | Inradius returns the radius of the inscribed sphere (tangent to every face). |
| `func (s PlatonicSolid) Mesh() *Polyhedron` | Mesh materialises the solid as a `Polyhedron` with exact vertex coordinates and outward-oriented polygonal faces at the solid's edge length, centered… |
| `func (s PlatonicSolid) Midradius() float64` | Midradius returns the radius of the midsphere (tangent to every edge at its midpoint). |
| `func (s PlatonicSolid) SurfaceArea() float64` | SurfaceArea returns the total surface area of the solid. |
| `func (s PlatonicSolid) Volume() float64` | Volume returns the enclosed volume of the solid. |

</details>

<details>
<summary><code>Polyhedron</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Bipyramid(n int, a, h float64) *Polyhedron` | Bipyramid returns the n-gonal bipyramid: two right pyramids joined base to base, with apexes at heights ±h and a regular n-gon of edge a as the… |
| `func ConvexHull(points []Vec3) (*Polyhedron, error)` | ConvexHull returns the convex hull of the given points as a triangulated `Polyhedron` whose vertices are exactly the hull vertices (interior and… |
| `func ConvexHullPolygons(points []Vec3) (*Polyhedron, error)` | ConvexHullPolygons returns the convex hull of the given points with coplanar triangles merged into maximal polygonal faces (so a box yields six… |
| `func GiftWrapHull(points []Vec3) (*Polyhedron, error)` | GiftWrapHull returns the convex hull of the given points computed with a gift-wrapping (pivoting) algorithm. |
| `func MergeCoplanarFaces(p *Polyhedron, tol float64) *Polyhedron` | MergeCoplanarFaces returns a copy of p in which adjacent triangular (or polygonal) faces that share an edge and lie in a common plane within angular… |
| `func NewPolyhedron(verts []Vec3, faces [][]int) *Polyhedron` | NewPolyhedron returns a Polyhedron with the given vertices and faces. |
| `func Prism(n int, a, h float64) *Polyhedron` | Prism returns a right uniform n-gonal prism whose regular n-gon base has edge length a and whose lateral height is h, centered on the z-axis. |
| `func Pyramid(n int, a, h float64) *Polyhedron` | Pyramid returns a right pyramid over a regular n-gon base of edge a with apex at height h above the base center, centered on the z-axis. |
| `func UniformAntiprism(n int, a float64) *Polyhedron` | UniformAntiprism returns the uniform n-gonal antiprism with edge length a: two parallel regular n-gons, the top rotated by π/n relative to the… |
| `func (p *Polyhedron) BoundingBox() (min, max Vec3, ok bool)` | BoundingBox returns the axis-aligned minimum and maximum corners of the vertex set and true, or two zero vectors and false when there are no vertices. |
| `func (p *Polyhedron) Centered() *Polyhedron` | Centered returns a copy of the polyhedron translated so that its solid centroid lies at the origin. |
| `func (p *Polyhedron) Centroid() (Vec3, bool)` | Centroid returns the centroid of the solid region enclosed by the mesh (the center of mass of a uniform-density body) and true, or the vertex… |
| `func (p *Polyhedron) Circumradius() float64` | Circumradius returns the maximum distance from the centroid c to any vertex. |
| `func (p *Polyhedron) Clone() *Polyhedron` | Clone returns a deep copy of the polyhedron. |
| `func (p *Polyhedron) Dual() (*Polyhedron, error)` | Dual returns the polar/canonical dual polyhedron. |
| `func (p *Polyhedron) EdgeLengths() []float64` | EdgeLengths returns the length of every edge in EdgeList order. |
| `func (p *Polyhedron) EdgeList() []Edge` | EdgeList returns all distinct undirected edges of the polyhedron, sorted lexicographically by endpoint indices. |
| `func (p *Polyhedron) EulerCharacteristic() int` | EulerCharacteristic returns V - E + F. For a simply connected closed polyhedron (topologically a sphere) it equals 2. |
| `func (p *Polyhedron) FaceArea(i int) float64` | FaceArea returns the area of face i. |
| `func (p *Polyhedron) FaceAreas() []float64` | FaceAreas returns the area of every face. |
| `func (p *Polyhedron) FaceCentroid(i int) (Vec3, bool)` | FaceCentroid returns the area-weighted centroid of face i and true, or the zero vector and false when face i is empty. |
| `func (p *Polyhedron) FaceCentroids() []Vec3` | FaceCentroids returns the centroid of every face (the zero vector for any empty face). |
| `func (p *Polyhedron) FaceLoop(i int) []Vec3` | FaceLoop returns the vertex positions of face i as an ordered slice. |
| `func (p *Polyhedron) FaceNormal(i int) (Vec3, bool)` | FaceNormal returns the outward unit normal of face i (via Newell's method) and true, or the zero vector and false when the face is degenerate. |
| `func (p *Polyhedron) FaceNormals() []Vec3` | FaceNormals returns the unit normal of every face. |
| `func (p *Polyhedron) FaceSides() []int` | FaceSides returns the number of sides (vertices) of every face. |
| `func (p *Polyhedron) FacesAroundVertex(v int) []int` | FacesAroundVertex returns the indices of the faces incident to vertex v, in ascending face-index order. |
| `func (p *Polyhedron) Inradius() float64` | Inradius returns the minimum distance from the centroid to any face plane, i.e. |
| `func (p *Polyhedron) IsClosed() bool` | IsClosed reports whether the mesh is closed, i.e. |
| `func (p *Polyhedron) Midradius() float64` | Midradius returns the mean distance from the centroid to the edge midpoints (the midradius of a canonical solid, where all such distances coincide). |
| `func (p *Polyhedron) MinVertexRadius() float64` | MinVertexRadius returns the minimum distance from the centroid to any vertex. |
| `func (p *Polyhedron) NumEdges() int` | NumEdges returns the number of distinct undirected edges across all faces. |
| `func (p *Polyhedron) NumFaces() int` | NumFaces returns the number of faces. |
| `func (p *Polyhedron) NumVertices() int` | NumVertices returns the number of vertices. |
| `func (p *Polyhedron) OrientOutward() *Polyhedron` | OrientOutward returns a copy of the mesh in which every face is oriented so that its normal points away from the mesh centroid. |
| `func (p *Polyhedron) ScaledBy(s float64) *Polyhedron` | ScaledBy returns a copy of the polyhedron with every vertex scaled about the origin by the factor s. |
| `func (p *Polyhedron) SignedVolume() float64` | SignedVolume returns the signed enclosed volume computed from the divergence theorem by fanning each face into triangles from its first vertex and… |
| `func (p *Polyhedron) SurfaceArea() float64` | SurfaceArea returns the total surface area of the polyhedron, the sum of its face areas. |
| `func (p *Polyhedron) Translate(d Vec3) *Polyhedron` | Translate returns a copy of the polyhedron with every vertex shifted by d. |
| `func (p *Polyhedron) Triangulate() *Polyhedron` | Triangulate returns a new polyhedron whose faces are all triangles, obtained by fanning every face from its first vertex. |
| `func (p *Polyhedron) Validate() error` | Validate reports the first structural problem with the polyhedron, or nil when it is well formed: every face must have at least three vertices and… |
| `func (p *Polyhedron) VertexCentroid() (Vec3, bool)` | VertexCentroid returns the arithmetic mean of the vertex positions and true, or the zero vector and false when the polyhedron has no vertices. |
| `func (p *Polyhedron) VertexDegrees() []int` | VertexDegrees returns, for each vertex index, the number of edges incident to it. |
| `func (p *Polyhedron) Volume() float64` | Volume returns the (non-negative) enclosed volume of the polyhedron. |

</details>

<details>
<summary><code>Vec3</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Centroid(ps []Vec3) (Vec3, bool)` | Centroid returns the arithmetic mean of the points in ps and true, or the zero vector and false when ps is empty. |
| `func NewVec3(x, y, z float64) Vec3` | NewVec3 returns the vector with the given components. |
| `func NewellNormal(loop []Vec3) Vec3` | NewellNormal returns the (unnormalised) normal of the planar or near-planar polygon given by the ordered loop of vertices, computed with Newell's… |
| `func PolygonCentroid(loop []Vec3) (Vec3, bool)` | PolygonCentroid returns the area-weighted centroid of the planar polygon given by the ordered loop, and true. |
| `func PolygonUnitNormal(loop []Vec3) (Vec3, bool)` | PolygonUnitNormal returns the unit normal of the polygon loop (following the right-hand rule around the vertex order) and true, or the zero vector… |
| `func RegularPolygonVertices(n int, r, height, offset float64) []Vec3` | RegularPolygonVertices returns the vertices of a regular n-gon of circumradius r lying in the plane z = height, starting at angle offset (radians)… |
| `func Sum(ps []Vec3) Vec3` | Sum returns the vector sum of all points in ps (the zero vector when empty). |
| `func TriangleCentroid(a, b, c Vec3) Vec3` | TriangleCentroid returns the centroid (arithmetic mean of vertices) of triangle (a, b, c). |
| `func TriangleNormal(a, b, c Vec3) Vec3` | TriangleNormal returns the (unnormalised) normal (b-a) × (c-a) of the triangle (a, b, c). |
| `func TriangleUnitNormal(a, b, c Vec3) (Vec3, bool)` | TriangleUnitNormal returns the unit normal of triangle (a, b, c) and true, or the zero vector and false when the triangle is degenerate. |
| `func XAxis() Vec3` | XAxis returns the unit vector (1, 0, 0). |
| `func YAxis() Vec3` | YAxis returns the unit vector (0, 1, 0). |
| `func ZAxis() Vec3` | ZAxis returns the unit vector (0, 0, 1). |
| `func ZeroVec3() Vec3` | ZeroVec3 returns the zero vector. |
| `func (v Vec3) Abs() Vec3` | Abs returns the component-wise absolute value of v. |
| `func (v Vec3) Add(w Vec3) Vec3` | Add returns the component-wise sum v + w. |
| `func (v Vec3) Angle(w Vec3) float64` | Angle returns the unsigned angle in radians between v and w, in [0, π]. |
| `func (v Vec3) Cross(w Vec3) Vec3` | Cross returns the cross product v × w, a vector orthogonal to both operands whose length equals the area of the parallelogram they span. |
| `func (v Vec3) Distance(w Vec3) float64` | Distance returns the Euclidean distance between points v and w. |
| `func (v Vec3) DistanceSq(w Vec3) float64` | DistanceSq returns the squared Euclidean distance between points v and w. |
| `func (v Vec3) Div(s float64) Vec3` | Div returns v with each component divided by the scalar s. |
| `func (v Vec3) Dot(w Vec3) float64` | Dot returns the Euclidean dot (inner) product v · w. |
| `func (v Vec3) Equal(w Vec3, tol float64) bool` | Equal reports whether v and w are equal within absolute tolerance tol in each component. |
| `func (v Vec3) Hadamard(w Vec3) Vec3` | Hadamard returns the component-wise (Hadamard) product of v and w. |
| `func (v Vec3) IsZero() bool` | IsZero reports whether every component of v is within Eps of zero. |
| `func (v Vec3) Len() float64` | Len returns the Euclidean length (2-norm) of v. |
| `func (v Vec3) LenSq() float64` | LenSq returns the squared Euclidean length of v. |
| `func (v Vec3) Lerp(w Vec3, t float64) Vec3` | Lerp returns the linear interpolation (1-t)*v + t*w. |
| `func (v Vec3) MaxAbsComp() float64` | MaxAbsComp returns the largest absolute value among the components of v. |
| `func (v Vec3) MaxComp(w Vec3) Vec3` | MaxComp returns the component-wise maximum of v and w. |
| `func (v Vec3) Midpoint(w Vec3) Vec3` | Midpoint returns the point halfway between v and w. |
| `func (v Vec3) MinComp(w Vec3) Vec3` | MinComp returns the component-wise minimum of v and w. |
| `func (v Vec3) Neg() Vec3` | Neg returns the additive inverse -v. |
| `func (v Vec3) Normalize() (Vec3, bool)` | Normalize returns the unit vector in the direction of v and true, or the zero vector and false when v has length below Eps. |
| `func (v Vec3) Project(w Vec3) Vec3` | Project returns the vector projection of v onto w. |
| `func (v Vec3) Reflect(n Vec3) Vec3` | Reflect returns v reflected across the plane through the origin with unit normal n, i.e. |
| `func (v Vec3) Reject(w Vec3) Vec3` | Reject returns the component of v orthogonal to w, that is v minus its projection onto w. |
| `func (v Vec3) Scale(s float64) Vec3` | Scale returns the vector v scaled by the scalar s. |
| `func (v Vec3) String() string` | String renders v in a compact, deterministic textual form. |
| `func (v Vec3) Sub(w Vec3) Vec3` | Sub returns the component-wise difference v - w. |
| `func (v Vec3) Triple(w, u Vec3) float64` | Triple returns the scalar triple product v · (w × u), equal to the signed volume of the parallelepiped spanned by the three vectors. |
| `func (v Vec3) Unit() Vec3` | Unit returns the unit vector in the direction of v, or the zero vector if v is (numerically) zero. |

</details>

### Constants

`Eps`

### Variables

`ErrDegenerate`, `ErrEmptyMesh`, `ErrNotEnoughPoints`, `Phi`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/polyhedra).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
