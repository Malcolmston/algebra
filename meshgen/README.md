# meshgen — implements two- and three-dimensional mesh generation and iso-contouring using only

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/meshgen.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/meshgen)

Package meshgen implements two- and three-dimensional mesh generation and
iso-contouring using only the Go standard library.

The package is self-contained and deterministic: identical inputs (including
any caller-supplied random seed) always produce identical outputs. All
coordinates are float64.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/meshgen"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func Barycentric(p, a, b, c Vec2) (u, v, w float64, ok bool)` | Barycentric returns the barycentric coordinates (u, v, w) of p with respect to the triangle (a, b, c), so that p = u*a + v*b + w*c and u+v+w = 1. |
| `func Circumradius(a, b, c Vec2) (float64, error)` | Circumradius returns the radius of the circle passing through a, b and c. |
| `func Collinear(a, b, c Vec2, eps float64) bool` | Collinear reports whether a, b and c are collinear within tolerance eps applied to twice the signed area. |
| `func DistancePointSegment(p, a, b Vec2) float64` | DistancePointSegment returns the distance from p to the segment a-b. |
| `func InCircle(a, b, c, d Vec2) float64` | InCircle returns a positive value when d lies inside the circle through a, b, c (assumed counterclockwise), negative when outside, and zero when the… |
| `func InCircleCCW(a, b, c, d Vec2) float64` | InCircleCCW returns the in-circle determinant with the first three points re-ordered so they are counterclockwise, so that a positive result always… |
| `func Inradius(a, b, c Vec2) float64` | Inradius returns the radius of the inscribed circle of the triangle (a,b,c). |
| `func IsCCW(a, b, c Vec2) bool` | IsCCW reports whether the triangle (a, b, c) is oriented counterclockwise. |
| `func IsCW(a, b, c Vec2) bool` | IsCW reports whether the triangle (a, b, c) is oriented clockwise. |
| `func IsSliver(a, b, c Vec2, minAngleDeg float64) bool` | IsSliver reports whether the triangle is a sliver: it has positive area but a minimum angle below minAngleDeg degrees. |
| `func IsWellShaped(a, b, c Vec2, minAngleDeg float64) bool` | IsWellShaped reports whether the triangle's minimum angle is at least minAngleDeg degrees. |
| `func MaxFloat(xs []float64) float64` | MaxFloat returns the maximum of the slice, or -Inf for an empty slice. |
| `func MeanFloat(xs []float64) float64` | MeanFloat returns the arithmetic mean of the slice, or 0 for an empty slice. |
| `func Median(xs []float64) float64` | Median returns the median of the slice, or 0 for an empty slice. |
| `func MinFloat(xs []float64) float64` | MinFloat returns the minimum of the slice, or +Inf for an empty slice. |
| `func Orient2D(a, b, c Vec2) float64` | Orient2D returns twice the signed area of the triangle (a, b, c). |
| `func Orient3D(a, b, c, d Vec3) float64` | Orient3D returns six times the signed volume of the tetrahedron (a,b,c,d). |
| `func Orientation(a, b, c Vec2, eps float64) int` | Orientation classifies the turn a -> b -> c as +1 (counterclockwise), -1 (clockwise) or 0 (collinear) using the tolerance eps on the signed area. |
| `func PointInPolygon(p Vec2, poly []Vec2) bool` | PointInPolygon reports whether p lies inside the simple polygon using the even-odd (ray casting) rule. |
| `func PointInTriangle(p, a, b, c Vec2, eps float64) bool` | PointInTriangle reports whether p lies inside or on the boundary of the triangle (a, b, c), using tolerance eps. |
| `func PolygonArea(poly []Vec2) float64` | PolygonArea returns the non-negative area of the polygon. |
| `func PolygonIsCCW(poly []Vec2) bool` | PolygonIsCCW reports whether the polygon is oriented counterclockwise. |
| `func PolygonPerimeter(poly []Vec2) float64` | PolygonPerimeter returns the total edge length of the closed polygon. |
| `func PolygonSignedArea(poly []Vec2) float64` | PolygonSignedArea returns the signed area of the polygon given by its ordered vertices; positive for counterclockwise orientation. |
| `func SegmentsProperlyIntersect(a, b, c, d Vec2) bool` | SegmentsProperlyIntersect reports whether the open segments a-b and c-d cross at a single interior point. |
| `func SignedArea2(a, b, c Vec2) float64` | SignedArea2 returns the signed area of the triangle (a, b, c); positive for counterclockwise orientation. |
| `func StdDev(xs []float64) float64` | StdDev returns the population standard deviation of the slice, or 0 for a slice with fewer than two elements. |
| `func SurfaceArea(tris []Triangle3) float64` | SurfaceArea returns the total area of the given triangles. |
| `func TetraVolume(a, b, c, d Vec3) float64` | TetraVolume returns the non-negative volume of the tetrahedron (a,b,c,d). |
| `func TotalSegmentLength(segs []Segment) float64` | TotalSegmentLength returns the sum of the lengths of the given segments. |
| `func TriAngles(a, b, c Vec2) [3]float64` | TriAngles returns the three interior angles of the triangle at vertices a, b and c respectively, in radians. |
| `func TriAspectRatio(a, b, c Vec2) float64` | TriAspectRatio returns the aspect ratio (longest edge over twice the inradius) of the triangle; 1 for an equilateral triangle up to scaling. |
| `func TriEdgeRatio(a, b, c Vec2) float64` | TriEdgeRatio returns the ratio of the longest to the shortest edge length. |
| `func TriMaxAngle(a, b, c Vec2) float64` | TriMaxAngle returns the largest interior angle, in radians, of the triangle. |
| `func TriMinAngle(a, b, c Vec2) float64` | TriMinAngle returns the smallest interior angle, in radians, of the triangle with vertices a, b, c. |
| `func TriRadiusRatio(a, b, c Vec2) float64` | TriRadiusRatio returns the radius ratio 2*inradius/circumradius, which lies in [0,1] and equals 1 for an equilateral triangle. |
| `func TriShapeRegularity(a, b, c Vec2) float64` | TriShapeRegularity returns 4*sqrt(3)*area / (sum of squared edge lengths), a scale-invariant quality in [0,1] that is 1 for an equilateral triangle. |
| `func TriangleArea(a, b, c Vec2) float64` | TriangleArea returns the non-negative area of the triangle (a, b, c). |
| `func TriangleArea3(a, b, c Vec3) float64` | TriangleArea3 returns the area of the triangle (a, b, c) in space. |

### Types

| Type | What it is |
| --- | --- |
| `Edge` | Edge is an unordered pair of vertex indices with U <= V. |
| `Mesh` | Mesh is a triangular mesh: a slice of vertices and a slice of triangles that index into it. |
| `MeshStats` | MeshStats summarises area and quality statistics over an entire mesh. |
| `Quad` | Quad is an axis-aligned rectangular cell defined by its minimum and maximum corners. |
| `QuadNode` | QuadNode is a node of a Quadtree. |
| `Quadtree` | Quadtree is an adaptively refined spatial subdivision of a rectangular region. |
| `ScalarGrid2` | ScalarGrid2 is a scalar field sampled on a regular rectangular lattice of Nx by Ny nodes. |
| `ScalarGrid3` | ScalarGrid3 is a scalar field sampled on a regular Nx by Ny by Nz lattice. |
| `Segment` | Segment is a straight line segment between two planar points. |
| `Tri` | Tri is a triangle described by three indices into a vertex slice. |
| `Triangle3` | Triangle3 is a triangle in space, described by three vertex positions. |
| `TriangleQuality` | TriangleQuality bundles several shape-quality measures of a single triangle. |
| `Triangulation` | Triangulation is a mutable planar triangulation supporting incremental point insertion (Bowyer-Watson), constraint recovery and refinement. |
| `Vec2` | Vec2 is a point or vector in the Euclidean plane. |
| `Vec3` | Vec3 is a point or vector in three-dimensional Euclidean space. |

<details>
<summary><code>Edge</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewEdge(a, b int) Edge` | NewEdge returns the canonical edge for the pair (a, b) with U <= V. |
| `func (e Edge) Other(v int) (int, bool)` | Other returns the endpoint of the edge that is not v, and true when v is an endpoint of the edge. |
| `func (e Edge) String() string` | String returns a human-readable representation of the edge. |

</details>

<details>
<summary><code>Mesh</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DelaunayMesh(points []Vec2) (*Mesh, error)` | DelaunayMesh is a convenience wrapper returning the Delaunay triangulation of the points directly as a Mesh. |
| `func GridMesh(nx, ny int, x0, y0, dx, dy float64) *Mesh` | GridMesh returns a triangulated regular grid of nx by ny nodes (so nx and ny must be at least two) with the given origin and spacing. |
| `func LaplacianSmooth(m *Mesh, iterations int, relaxation float64, pinBoundary bool) *Mesh` | LaplacianSmooth returns a copy of the mesh after applying iterations passes of Laplacian smoothing. |
| `func LloydSmooth(m *Mesh, iterations int, pinBoundary bool) *Mesh` | LloydSmooth returns a copy of the mesh after applying iterations passes of centroidal (Lloyd) smoothing. |
| `func NewMesh(verts []Vec2, tris []Tri) *Mesh` | NewMesh returns a mesh with the given vertices and triangles. |
| `func RefineToMinAngle(points []Vec2, minAngleDeg float64, maxSteiner int) (*Mesh, error)` | RefineToMinAngle triangulates the points and applies Delaunay refinement to the given minimum angle, returning the refined mesh. |
| `func TaubinSmooth(m *Mesh, iterations int, lambda, mu float64, pinBoundary bool) *Mesh` | TaubinSmooth returns a copy of the mesh after applying iterations Taubin lambda/mu passes, a low-pass smoother that reduces the shrinkage of repeated… |
| `func TriangulatePolygonEarClip(poly []Vec2) (*Mesh, error)` | TriangulatePolygonEarClip returns an ear-clipping triangulation of a simple (possibly non-convex) polygon given in counterclockwise or clockwise… |
| `func TriangulatePolygonFan(poly []Vec2) (*Mesh, error)` | TriangulatePolygonFan returns a triangle-fan triangulation of a convex polygon given by its ordered vertices, as a Mesh. |
| `func (m *Mesh) AngleHistogram(bins int) []int` | AngleHistogram bins the interior angles of all triangles (in degrees) into the given number of equal-width bins spanning [0, 180]. |
| `func (m *Mesh) BoundaryEdges() []Edge` | BoundaryEdges returns the edges incident to exactly one triangle, sorted. |
| `func (m *Mesh) BoundaryLength() float64` | BoundaryLength returns the total length of all boundary edges. |
| `func (m *Mesh) BoundaryVertices() []int` | BoundaryVertices returns the sorted indices of vertices lying on the mesh boundary. |
| `func (m *Mesh) Bounds() (min, max Vec2)` | Bounds returns the axis-aligned bounding box of the mesh vertices. |
| `func (m *Mesh) Clone() *Mesh` | Clone returns a deep copy of the mesh. |
| `func (m *Mesh) ConnectedComponents() int` | ConnectedComponents returns the number of connected components of the mesh graph, counting only vertices referenced by at least one triangle. |
| `func (m *Mesh) CountSlivers(minAngleDeg float64) int` | CountSlivers returns the number of triangles whose minimum angle is below minAngleDeg degrees. |
| `func (m *Mesh) EdgeLengths() []float64` | EdgeLengths returns the length of every distinct edge, in edge order. |
| `func (m *Mesh) EdgeTriangleCount() map[Edge]int` | EdgeTriangleCount maps each edge to the number of triangles incident to it. |
| `func (m *Mesh) Edges() []Edge` | Edges returns every distinct edge of the mesh, sorted. |
| `func (m *Mesh) EulerCharacteristic() int` | EulerCharacteristic returns V - E + F for the mesh, where F counts the triangular faces. |
| `func (m *Mesh) InteriorEdges() []Edge` | InteriorEdges returns the edges shared by two or more triangles, sorted. |
| `func (m *Mesh) IsBoundaryEdge(a, b int) bool` | IsBoundaryEdge reports whether the edge (a,b) lies on the mesh boundary. |
| `func (m *Mesh) IsManifold() bool` | IsManifold reports whether every edge is shared by at most two triangles. |
| `func (m *Mesh) MaxAngleDeg() float64` | MaxAngleDeg returns the largest interior angle over all triangles of the mesh, in degrees. |
| `func (m *Mesh) MinAngleDeg() float64` | MinAngleDeg returns the smallest interior angle over all triangles of the mesh, in degrees. |
| `func (m *Mesh) NumEdges() int` | NumEdges returns the number of distinct edges in the mesh. |
| `func (m *Mesh) NumTriangles() int` | NumTriangles returns the number of triangles in the mesh. |
| `func (m *Mesh) NumVertices() int` | NumVertices returns the number of vertices in the mesh. |
| `func (m *Mesh) OrientCCW() *Mesh` | OrientCCW returns a copy of the mesh with every triangle re-ordered to be counterclockwise. |
| `func (m *Mesh) Stats() MeshStats` | Stats computes area and quality statistics for the mesh. |
| `func (m *Mesh) TotalArea() float64` | TotalArea returns the sum of the areas of all triangles. |
| `func (m *Mesh) TriangleAdjacency() [][]int` | TriangleAdjacency returns, for each triangle, the sorted indices of triangles sharing an edge with it. |
| `func (m *Mesh) TriangleArea(i int) float64` | TriangleArea returns the area of triangle i. |
| `func (m *Mesh) TriangleAreas() []float64` | TriangleAreas returns the area of every triangle in triangle order. |
| `func (m *Mesh) TriangleCentroid(i int) Vec2` | TriangleCentroid returns the centroid of triangle i. |
| `func (m *Mesh) TriangleMinAngles() []float64` | TriangleMinAngles returns the minimum interior angle (radians) of every triangle in triangle order. |
| `func (m *Mesh) TriangleQualities() []float64` | TriangleQualities returns the radius-ratio quality of every triangle in triangle order. |
| `func (m *Mesh) TriangleQualityAt(i int) TriangleQuality` | TriangleQualityAt returns the quality of triangle i of the mesh. |
| `func (m *Mesh) TriangleVertices(i int) (a, b, c Vec2)` | TriangleVertices returns the three vertex positions of triangle i. |
| `func (m *Mesh) Validate() error` | Validate reports the first structural error found in the mesh, or nil when every triangle references distinct, in-range vertices. |
| `func (m *Mesh) VertexDegree(v int) int` | VertexDegree returns the number of edges incident to vertex v. |
| `func (m *Mesh) VertexNeighbors() [][]int` | VertexNeighbors returns, for each vertex, the sorted set of adjacent vertex indices connected by an edge. |
| `func (m *Mesh) VertexTriangles() [][]int` | VertexTriangles returns, for each vertex index, the sorted list of triangle indices incident to it. |

</details>

<details>
<summary><code>Quad</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewQuad(min, max Vec2) Quad` | NewQuad returns the quad with the given corners. |
| `func (q Quad) Area() float64` | Area returns the area of the quad. |
| `func (q Quad) Center() Vec2` | Center returns the centre point of the quad. |
| `func (q Quad) Contains(p Vec2) bool` | Contains reports whether p lies within the quad, including its boundary. |
| `func (q Quad) Corners() [4]Vec2` | Corners returns the four corners of the quad in counterclockwise order starting from the minimum corner. |
| `func (q Quad) Height() float64` | Height returns the extent of the quad along Y. |
| `func (q Quad) Subdivide() [4]Quad` | Subdivide splits the quad into its four child quadrants in the order SW, SE, NE, NW. |
| `func (q Quad) Width() float64` | Width returns the extent of the quad along X. |

</details>

<details>
<summary><code>QuadNode</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (n *QuadNode) IsLeaf() bool` | IsLeaf reports whether the node has no children. |

</details>

<details>
<summary><code>Quadtree</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BuildQuadtree(bounds Quad, maxDepth int, refine func(cell Quad, depth int) bool) *Quadtree` | BuildQuadtree builds a quadtree over bounds, subdividing a cell whenever refine reports true for it, down to at most maxDepth levels. |
| `func (t *Quadtree) Depth() int` | Depth returns the maximum depth reached by any node. |
| `func (t *Quadtree) LeafCenters() []Vec2` | LeafCenters returns the centre point of every leaf cell. |
| `func (t *Quadtree) LeafCorners() []Vec2` | LeafCorners returns the deduplicated set of all leaf-cell corner points of the quadtree, suitable as vertices for a background mesh. |
| `func (t *Quadtree) LeafCount() int` | LeafCount returns the number of leaf cells. |
| `func (t *Quadtree) Leaves() []Quad` | Leaves returns the leaf cells of the quadtree in depth-first order. |

</details>

<details>
<summary><code>ScalarGrid2</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewScalarGrid2(nx, ny int, x0, y0, dx, dy float64) *ScalarGrid2` | NewScalarGrid2 returns a zero-initialised grid of nx by ny nodes with the given origin and spacing. |
| `func SampleGrid2(nx, ny int, x0, y0, dx, dy float64, f func(x, y float64) float64) *ScalarGrid2` | SampleGrid2 builds a scalar grid by evaluating f at every node of the nx by ny lattice. |
| `func (g *ScalarGrid2) At(i, j int) float64` | At returns the value at node (i, j). |
| `func (g *ScalarGrid2) Node(i, j int) Vec2` | Node returns the position of node (i, j). |
| `func (g *ScalarGrid2) Set(i, j int, v float64)` | Set assigns v to node (i, j). |

</details>

<details>
<summary><code>ScalarGrid3</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewScalarGrid3(nx, ny, nz int, x0, y0, z0, dx, dy, dz float64) *ScalarGrid3` | NewScalarGrid3 returns a zero-initialised 3-D grid. |
| `func SampleGrid3(nx, ny, nz int, x0, y0, z0, dx, dy, dz float64, f func(x, y, z float64) float64) *ScalarGrid3` | SampleGrid3 builds a 3-D scalar grid by evaluating f at every node. |
| `func (g *ScalarGrid3) At(i, j, k int) float64` | At returns the value at node (i, j, k). |
| `func (g *ScalarGrid3) Node(i, j, k int) Vec3` | Node returns the position of node (i, j, k). |
| `func (g *ScalarGrid3) Set(i, j, k int, v float64)` | Set assigns v to node (i, j, k). |

</details>

<details>
<summary><code>Segment</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func IsoContour(nx, ny int, x0, y0, dx, dy float64, f func(x, y float64) float64, iso float64) []Segment` | IsoContour is a convenience wrapper that samples f on an nx by ny lattice and returns the marching-squares contour at the given isovalue. |
| `func MarchingSquares(g *ScalarGrid2, iso float64) []Segment` | MarchingSquares extracts the iso-contour of the grid at the given isovalue as a set of line segments. |
| `func NewSegment(a, b Vec2) Segment` | NewSegment returns the segment from a to b. |
| `func (s Segment) Direction() Vec2` | Direction returns the unit direction from A to B. |
| `func (s Segment) Length() float64` | Length returns the Euclidean length of the segment. |
| `func (s Segment) Midpoint() Vec2` | Midpoint returns the midpoint of the segment. |

</details>

<details>
<summary><code>Tri</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewTri(a, b, c int) Tri` | NewTri returns the triangle with the given vertex indices. |
| `func (t Tri) EdgesOf() [3]Edge` | EdgesOf returns the three (unordered, sorted) edges of the triangle. |
| `func (t Tri) Has(v int) bool` | Has reports whether the triangle references vertex index v. |
| `func (t Tri) Indices() [3]int` | Indices returns the three vertex indices as an array. |
| `func (t Tri) Opposite(u, v int) (int, bool)` | Opposite returns the vertex of the triangle that is neither u nor v, and a boolean reporting whether exactly one such vertex exists. |
| `func (t Tri) SharesEdge(u Tri) (Edge, bool)` | SharesEdge reports whether triangles t and u share an edge, returning that edge and true when they do. |

</details>

<details>
<summary><code>Triangle3</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func MarchingCubes(g *ScalarGrid3, iso float64) []Triangle3` | MarchingCubes extracts the iso-surface of the grid at the given isovalue as a triangle soup. |
| `func MarchingTetrahedron(p [4]Vec3, v [4]float64, iso float64) []Triangle3` | MarchingTetrahedron returns the iso-surface triangles for a single tetrahedron whose corners are p with scalar values v. |
| `func (t Triangle3) Area() float64` | Area returns the area of the triangle. |
| `func (t Triangle3) Centroid() Vec3` | Centroid returns the centroid of the triangle. |
| `func (t Triangle3) Normal() Vec3` | Normal returns the unit normal of the triangle following the right-hand rule. |

</details>

<details>
<summary><code>TriangleQuality</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func TriangleQualityOf(a, b, c Vec2) TriangleQuality` | TriangleQualityOf computes the quality measures of the triangle (a, b, c). |

</details>

<details>
<summary><code>Triangulation</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Triangulate(points []Vec2) (*Triangulation, error)` | Triangulate returns the Delaunay triangulation of the given planar points using the incremental Bowyer-Watson algorithm. |
| `func TriangulateConstrained(points []Vec2, segments [][2]int) (*Triangulation, error)` | TriangulateConstrained returns the constrained Delaunay triangulation of the points in which every segment (a pair of input point indices) is… |
| `func (t *Triangulation) InsertConstraint(a, b int) error` | InsertConstraint forces the edge between input vertices a and b to appear in the triangulation by flipping the diagonals that cross it. |
| `func (t *Triangulation) IsDelaunay(eps float64) bool` | IsDelaunay reports whether every real triangle satisfies the empty-circle property with respect to all input vertices, within tolerance eps. |
| `func (t *Triangulation) Mesh() *Mesh` | Mesh returns the triangulation as a Mesh, dropping triangles that touch the super-triangle and renumbering vertices to a compact 0-based index set. |
| `func (t *Triangulation) NumInput() int` | NumInput returns the number of original input points. |
| `func (t *Triangulation) NumTriangles() int` | NumTriangles returns the number of real (non-super) triangles. |
| `func (t *Triangulation) Points() []Vec2` | Points returns the input points (excluding internal helper vertices). |
| `func (t *Triangulation) Refine(minAngleDeg float64, maxSteiner int) int` | Refine performs Ruppert-style Delaunay refinement guarded for robustness: it repeatedly picks the worst real triangle whose minimum angle is below… |

</details>

<details>
<summary><code>Vec2</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BoundingBox2(pts []Vec2) (min, max Vec2)` | BoundingBox2 returns the axis-aligned minimum and maximum corners enclosing the given points. |
| `func CentroidVec2(pts []Vec2) Vec2` | CentroidVec2 returns the arithmetic mean of the given points. |
| `func Circumcenter(a, b, c Vec2) (Vec2, error)` | Circumcenter returns the center of the circle passing through a, b and c. |
| `func Circumcircle(a, b, c Vec2) (center Vec2, radius float64, err error)` | Circumcircle returns the center and radius of the circle through a, b and c. |
| `func ClosestPointOnSegment(p, a, b Vec2) Vec2` | ClosestPointOnSegment returns the point on the segment a-b closest to p. |
| `func ConvexHull(points []Vec2) []Vec2` | ConvexHull returns the counterclockwise convex hull of the points using Andrew's monotone chain algorithm. |
| `func GridPoints(nx, ny int, x0, y0, dx, dy float64) []Vec2` | GridPoints returns the nodes of a regular nx by ny lattice in row-major order. |
| `func Incenter(a, b, c Vec2) Vec2` | Incenter returns the incenter of the triangle (a, b, c), the center of its inscribed circle. |
| `func JitteredGridPoints(nx, ny int, x0, y0, dx, dy, jitter float64, seed int64) []Vec2` | JitteredGridPoints returns nx by ny points laid out on a regular lattice and perturbed by up to jitter times the cell size in each direction, using a… |
| `func NewVec2(x, y float64) Vec2` | NewVec2 returns the Vec2 with the given coordinates. |
| `func PolygonCentroid(poly []Vec2) Vec2` | PolygonCentroid returns the area centroid of a simple polygon. |
| `func RandomPoints(n int, seed int64, min, max Vec2) []Vec2` | RandomPoints returns n points drawn uniformly from the axis-aligned rectangle with the given corners, using a deterministic generator seeded by seed. |
| `func RandomPointsInTriangle(n int, seed int64, a, b, c Vec2) []Vec2` | RandomPointsInTriangle returns n points drawn uniformly from the triangle (a, b, c), using a deterministic generator seeded by seed. |
| `func SegmentIntersection(a, b, c, d Vec2) (Vec2, bool)` | SegmentIntersection returns the intersection point of the lines through a-b and c-d. |
| `func Vec2Zero() Vec2` | Vec2Zero is the origin (0, 0). |
| `func (p Vec2) Abs() Vec2` | Abs returns the componentwise absolute value of p. |
| `func (p Vec2) Add(q Vec2) Vec2` | Add returns the componentwise sum p + q. |
| `func (p Vec2) Angle() float64` | Angle returns the angle of p measured counterclockwise from the +X axis, in the range (-pi, pi]. |
| `func (p Vec2) AngleTo(q Vec2) float64` | AngleTo returns the signed angle from p to q in radians, in (-pi, pi]. |
| `func (p Vec2) ApproxEqual(q Vec2, eps float64) bool` | ApproxEqual reports whether p and q agree to within absolute tolerance eps. |
| `func (p Vec2) Cross(q Vec2) float64` | Cross returns the scalar 2-D cross product p.X*q.Y - p.Y*q.X. |
| `func (p Vec2) Distance(q Vec2) float64` | Distance returns the Euclidean distance between p and q. |
| `func (p Vec2) DistanceSq(q Vec2) float64` | DistanceSq returns the squared Euclidean distance between p and q. |
| `func (p Vec2) Div(s float64) Vec2` | Div returns p divided componentwise by the nonzero scalar s. |
| `func (p Vec2) Dot(q Vec2) float64` | Dot returns the dot product p . |
| `func (p Vec2) Equal(q Vec2) bool` | Equal reports whether p and q are exactly equal. |
| `func (p Vec2) Hadamard(q Vec2) Vec2` | Hadamard returns the componentwise product of p and q. |
| `func (p Vec2) IsFinite() bool` | IsFinite reports whether both coordinates of p are finite. |
| `func (p Vec2) IsZero() bool` | IsZero reports whether p is exactly the zero vector. |
| `func (p Vec2) Lerp(q Vec2, t float64) Vec2` | Lerp returns the linear interpolation (1-t)*p + t*q. |
| `func (p Vec2) Max(q Vec2) Vec2` | Max returns the componentwise maximum of p and q. |
| `func (p Vec2) Midpoint(q Vec2) Vec2` | Midpoint returns the midpoint of the segment p-q. |
| `func (p Vec2) Min(q Vec2) Vec2` | Min returns the componentwise minimum of p and q. |
| `func (p Vec2) Neg() Vec2` | Neg returns the vector reflected through the origin, -p. |
| `func (p Vec2) Norm() float64` | Norm returns the Euclidean norm \|p\|. |
| `func (p Vec2) NormSq() float64` | NormSq returns the squared Euclidean norm \|p\|^2. |
| `func (p Vec2) Normalize() Vec2` | Normalize returns the unit vector in the direction of p. |
| `func (p Vec2) Perp() Vec2` | Perp returns the vector p rotated 90 degrees counterclockwise, (-p.Y, p.X). |
| `func (p Vec2) Rotate(theta float64) Vec2` | Rotate returns p rotated counterclockwise by theta radians about the origin. |
| `func (p Vec2) RotateAbout(c Vec2, theta float64) Vec2` | RotateAbout returns p rotated counterclockwise by theta radians about c. |
| `func (p Vec2) Scale(s float64) Vec2` | Scale returns p scaled by the factor s. |
| `func (p Vec2) String() string` | String returns a human-readable representation of p. |
| `func (p Vec2) Sub(q Vec2) Vec2` | Sub returns the componentwise difference p - q. |

</details>

<details>
<summary><code>Vec3</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BoundingBox3(pts []Vec3) (min, max Vec3)` | BoundingBox3 returns the axis-aligned minimum and maximum corners enclosing the given points. |
| `func CentroidVec3(pts []Vec3) Vec3` | CentroidVec3 returns the arithmetic mean of the given points. |
| `func NewVec3(x, y, z float64) Vec3` | NewVec3 returns the Vec3 with the given coordinates. |
| `func TriangleNormal3(a, b, c Vec3) Vec3` | TriangleNormal3 returns the unit normal of the triangle (a, b, c), following the right-hand rule for the given vertex order. |
| `func Vec3Zero() Vec3` | Vec3Zero is the origin (0, 0, 0). |
| `func (p Vec3) Add(q Vec3) Vec3` | Add returns the componentwise sum p + q. |
| `func (p Vec3) AngleBetween(q Vec3) float64` | AngleBetween returns the unsigned angle in radians between p and q, in the range [0, pi]. |
| `func (p Vec3) ApproxEqual(q Vec3, eps float64) bool` | ApproxEqual reports whether p and q agree to within absolute tolerance eps. |
| `func (p Vec3) Cross(q Vec3) Vec3` | Cross returns the vector cross product p x q. |
| `func (p Vec3) Distance(q Vec3) float64` | Distance returns the Euclidean distance between p and q. |
| `func (p Vec3) DistanceSq(q Vec3) float64` | DistanceSq returns the squared Euclidean distance between p and q. |
| `func (p Vec3) Div(s float64) Vec3` | Div returns p divided componentwise by the nonzero scalar s. |
| `func (p Vec3) Dot(q Vec3) float64` | Dot returns the dot product p . |
| `func (p Vec3) Equal(q Vec3) bool` | Equal reports whether p and q are exactly equal. |
| `func (p Vec3) IsFinite() bool` | IsFinite reports whether every coordinate of p is finite. |
| `func (p Vec3) Lerp(q Vec3, t float64) Vec3` | Lerp returns the linear interpolation (1-t)*p + t*q. |
| `func (p Vec3) Max(q Vec3) Vec3` | Max returns the componentwise maximum of p and q. |
| `func (p Vec3) Midpoint(q Vec3) Vec3` | Midpoint returns the midpoint of the segment p-q. |
| `func (p Vec3) Min(q Vec3) Vec3` | Min returns the componentwise minimum of p and q. |
| `func (p Vec3) Neg() Vec3` | Neg returns the vector reflected through the origin, -p. |
| `func (p Vec3) Norm() float64` | Norm returns the Euclidean norm \|p\|. |
| `func (p Vec3) NormSq() float64` | NormSq returns the squared Euclidean norm \|p\|^2. |
| `func (p Vec3) Normalize() Vec3` | Normalize returns the unit vector in the direction of p. |
| `func (p Vec3) Scale(s float64) Vec3` | Scale returns p scaled by the factor s. |
| `func (p Vec3) String() string` | String returns a human-readable representation of p. |
| `func (p Vec3) Sub(q Vec3) Vec3` | Sub returns the componentwise difference p - q. |
| `func (p Vec3) TripleProduct(q, r Vec3) float64` | TripleProduct returns the scalar triple product p . |

</details>

### Constants

`Eps`

### Variables

`ContourError`, `ErrCollinear`, `ErrConstraintFailed`, `ErrDegenerate`, `ErrNotEnoughPoints`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/meshgen).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
