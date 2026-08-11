# voronoi — implements planar computational geometry centred on the Delaunay triangulation

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/voronoi.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/voronoi)

Package voronoi implements planar computational geometry centred on the
Delaunay triangulation and its dual, the Voronoi diagram.

The package is self-contained and depends only on the Go standard library. All
coordinates are float64 and every routine is deterministic: identical inputs
yield identical outputs.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/voronoi"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AllNearestNeighbors(pts []Point) []int` | AllNearestNeighbors returns, for each point in pts, the index of its nearest other point. |
| `func Barycentric(p, a, b, c Point) (u, v, w float64, ok bool)` | Barycentric returns the barycentric coordinates (u, v, w) of p with respect to the triangle (a, b, c), so that p = u*a + v*b + w*c and u+v+w = 1. |
| `func Circumradius(a, b, c Point) (float64, error)` | Circumradius returns the radius of the circumcircle of the triangle (a, b, c). |
| `func CircumradiusSpectrum(tri *Triangulation) []float64` | CircumradiusSpectrum returns the sorted list of circumradii of all Delaunay triangles. |
| `func ClosestPair(pts []Point) (int, int, float64)` | ClosestPair returns the two indices of the closest pair of points and their distance. |
| `func Collinear(a, b, c Point, eps float64) bool` | Collinear reports whether a, b and c are collinear within tolerance eps (applied to twice the signed area). |
| `func ConvexHullArea(pts []Point) float64` | ConvexHullArea returns the area enclosed by the convex hull of the points. |
| `func ConvexHullIndices(pts []Point) []int` | ConvexHullIndices returns the indices (into pts) of the convex-hull vertices in counterclockwise order. |
| `func ConvexHullPerimeter(pts []Point) float64` | ConvexHullPerimeter returns the perimeter of the convex hull of the points. |
| `func Diameter(pts []Point) (dist float64, i, j int)` | Diameter returns the maximum distance between any two of the given points (the diameter of the set) together with the indices of the realizing pair. |
| `func EnclosingCircleRadius(pts []Point) float64` | EnclosingCircleRadius returns the radius of the minimum enclosing circle of the given points. |
| `func InCircle(a, b, c, d Point) float64` | InCircle returns a positive value when d lies inside the circle through a, b, c (assumed counterclockwise), negative when outside, and zero when the… |
| `func InCircleTest(a, b, c, d Point) bool` | InCircleTest reports whether d lies strictly inside the circumcircle of the triangle (a, b, c). |
| `func Inradius(a, b, c Point) float64` | Inradius returns the radius of the inscribed circle of the triangle (a, b, c). |
| `func IsCCW(a, b, c Point) bool` | IsCCW reports whether the triangle (a, b, c) is oriented counterclockwise. |
| `func IsCW(a, b, c Point) bool` | IsCW reports whether the triangle (a, b, c) is oriented clockwise. |
| `func IsConvexPolygon(poly []Point) bool` | IsConvexPolygon reports whether the polygon given by its vertices in order is convex. |
| `func KNearest(pts []Point, q Point, k, skip int) []int` | KNearest returns the indices of the k points in pts closest to q, ordered from nearest to farthest. |
| `func NearestNeighbor(pts []Point, q Point, skip int) (int, float64)` | NearestNeighbor returns the index of the point in pts closest to the query point q, together with the distance. |
| `func Orient2D(a, b, c Point) float64` | Orient2D returns twice the signed area of the triangle (a, b, c). |
| `func Orientation(a, b, c Point, eps float64) int` | Orientation classifies the turn a -> b -> c as +1 (counterclockwise), -1 (clockwise) or 0 (collinear) using the tolerance eps on the signed area. |
| `func PointInConvexPolygon(p Point, poly []Point, eps float64) bool` | PointInConvexPolygon reports whether p lies inside or on the convex polygon whose vertices are given in counterclockwise order, within tolerance eps. |
| `func PointInPolygon(p Point, poly []Point) bool` | PointInPolygon reports whether p lies inside the simple polygon, using the even-odd (ray casting) rule. |
| `func PointInTriangle(p, a, b, c Point, eps float64) bool` | PointInTriangle reports whether p lies inside or on the boundary of the triangle (a, b, c), using tolerance eps on the orientation sign. |
| `func PolygonArea(poly []Point) float64` | PolygonArea returns the (non-negative) area of the simple polygon. |
| `func PolygonIsCCW(poly []Point) bool` | PolygonIsCCW reports whether the polygon's vertices wind counterclockwise. |
| `func PolygonPerimeter(poly []Point) float64` | PolygonPerimeter returns the perimeter of the polygon given by its vertices in order (treated as closed). |
| `func PolygonSignedArea(poly []Point) float64` | PolygonSignedArea returns the signed area of the simple polygon given by its vertices in order (via the shoelace formula). |
| `func PolygonWindingNumber(p Point, poly []Point) int` | PolygonWindingNumber returns the winding number of the polygon around p. |
| `func SegmentsIntersect(p1, p2, q1, q2 Point) bool` | SegmentsIntersect reports whether the closed segments p1p2 and q1q2 share at least one point. |
| `func TotalEdgeLength(pts []Point, edges []Edge) float64` | TotalEdgeLength returns the sum of the lengths of the given edges over the point set. |
| `func TriangleArea(a, b, c Point) float64` | TriangleArea returns the (non-negative) area of the triangle (a, b, c). |
| `func TrianglePerimeter(a, b, c Point) float64` | TrianglePerimeter returns the perimeter of the triangle (a, b, c). |
| `func TriangleSignedArea(a, b, c Point) float64` | TriangleSignedArea returns the signed area of the triangle (a, b, c), positive for counterclockwise orientation. |

### Types

| Type | What it is |
| --- | --- |
| `AlphaShape` | AlphaShape is the boundary of an alpha complex: the family of Delaunay triangles whose circumradius does not exceed the parameter Alpha, together… |
| `Circle` | Circle is a circle defined by its centre and radius. |
| `Edge` | Edge is an undirected edge between two vertices identified by their indices into a point slice. |
| `Point` | Point is a location in the Euclidean plane. |
| `Rect` | Rect is an axis-aligned rectangle described by its minimum and maximum corners. |
| `Triangle` | Triangle references three vertices by their indices into a point slice. |
| `Triangulation` | Triangulation is a Delaunay triangulation of a planar point set. |
| `VoronoiCell` | VoronoiCell is the region of the plane assigned to a single site. |
| `VoronoiDiagram` | VoronoiDiagram is the dual of a Delaunay triangulation: it partitions the plane into cells, one per input site, where each cell is the set of points… |
| `VoronoiEdge` | VoronoiEdge is a finite edge of the Voronoi diagram, joining the circumcentres of two adjacent Delaunay triangles. |

<details>
<summary><code>AlphaShape</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AlphaShapeFromTriangulation(tri *Triangulation, alpha float64) *AlphaShape` | AlphaShapeFromTriangulation computes the alpha shape from an existing triangulation without retriangulating. |
| `func AlphaShapeOf(pts []Point, alpha float64) (*AlphaShape, error)` | AlphaShapeOf computes the alpha shape of the Delaunay triangulation for the given alpha radius. |
| `func (as *AlphaShape) Area() float64` | Area returns the total area covered by the alpha complex (the sum of its triangle areas). |
| `func (as *AlphaShape) BoundaryPoints() []int` | BoundaryPoints returns the indices of points that lie on the alpha-shape boundary, sorted. |
| `func (as *AlphaShape) Edges() []Edge` | Edges returns the boundary edges of the alpha shape. |
| `func (as *AlphaShape) NumEdges() int` | NumEdges returns the number of boundary edges of the alpha shape. |
| `func (as *AlphaShape) NumTriangles() int` | NumTriangles returns the number of triangles in the alpha complex. |
| `func (as *AlphaShape) Perimeter() float64` | Perimeter returns the total length of the alpha shape's boundary edges. |
| `func (as *AlphaShape) Triangles() []int` | Triangles returns the indices of the Delaunay triangles that make up the alpha complex (those with circumradius <= Alpha). |

</details>

<details>
<summary><code>Circle</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func CircleFromDiameter(a, b Point) Circle` | CircleFromDiameter returns the circle whose diameter is the segment ab. |
| `func CircleThrough(a, b, c Point) (Circle, error)` | CircleThrough returns the unique circle passing through the three points. |
| `func Circumcircle(a, b, c Point) (Circle, error)` | Circumcircle returns the circumscribed circle of the triangle (a, b, c). |
| `func EmptyCircleAt(sites []Point, c Point) Circle` | EmptyCircleAt returns the largest circle centred at c that contains no site in its interior; its radius is the distance from c to the nearest site. |
| `func LargestEmptyCircle(sites []Point) (Circle, error)` | LargestEmptyCircle returns the largest circle whose centre lies inside the convex hull of the sites and whose interior contains no site. |
| `func MinimumEnclosingCircle(pts []Point) Circle` | MinimumEnclosingCircle returns the smallest circle that contains all of the given points. |
| `func NewCircle(center Point, radius float64) Circle` | NewCircle returns the circle with the given centre and radius. |
| `func (c Circle) Area() float64` | Area returns the area enclosed by the circle. |
| `func (c Circle) BoundingBox() Rect` | BoundingBox returns the axis-aligned bounding box of the circle. |
| `func (c Circle) Circumference() float64` | Circumference returns the circle's circumference. |
| `func (c Circle) Contains(p Point, eps float64) bool` | Contains reports whether p lies inside or on the circle, within tolerance eps applied to the radius. |
| `func (c Circle) ContainsStrict(p Point, eps float64) bool` | ContainsStrict reports whether p lies strictly inside the circle, within tolerance eps. |
| `func (c Circle) Diameter() float64` | Diameter returns the circle's diameter. |
| `func (c Circle) Intersects(o Circle) bool` | Intersects reports whether two circles overlap or touch. |
| `func (c Circle) OnBoundary(p Point, eps float64) bool` | OnBoundary reports whether p lies on the circle boundary within tolerance eps. |
| `func (c Circle) PointAt(theta float64) Point` | PointAt returns the point on the circle at angle theta (radians) measured counterclockwise from the positive X axis. |
| `func (c Circle) Sample(n int) []Point` | Sample returns n points spaced evenly around the circle, starting at angle 0. |
| `func (c Circle) String() string` | String returns a human-readable representation of the circle. |

</details>

<details>
<summary><code>Edge</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AlphaShapeEdges(pts []Point, alpha float64) ([]Edge, error)` | AlphaShapeEdges is a convenience wrapper returning just the boundary edges of the alpha shape of pts for the given alpha. |
| `func EuclideanMST(pts []Point) []Edge` | EuclideanMST returns the edges of a Euclidean minimum spanning tree of the points. |
| `func GabrielEdges(pts []Point) []Edge` | GabrielEdges returns the edges (i, j) whose diametral circle (the circle with segment ij as diameter) contains no other point. |
| `func NearestNeighborEdges(pts []Point) []Edge` | NearestNeighborEdges returns the directed nearest-neighbour edges of the set: entry i is the edge from point i to its nearest other point. |
| `func NewEdge(a, b int) Edge` | NewEdge returns the edge joining vertex indices a and b. |
| `func RelativeNeighborhoodEdges(pts []Point) []Edge` | RelativeNeighborhoodEdges returns the edges (i, j) such that no third point k is simultaneously closer to both i and j than they are to each other. |
| `func (e Edge) Canonical() Edge` | Canonical returns the edge with its endpoints ordered so that A <= B. This gives a stable key for undirected edges. |
| `func (e Edge) Contains(v int) bool` | Contains reports whether vertex index v is an endpoint of e. |
| `func (e Edge) Equal(f Edge) bool` | Equal reports whether e and f join the same pair of vertices, ignoring direction. |
| `func (e Edge) Other(v int) int` | Other returns the endpoint of e that is not v. |
| `func (e Edge) Reversed() Edge` | Reversed returns the edge with its endpoints swapped. |

</details>

<details>
<summary><code>Point</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Centroid(pts []Point) Point` | Centroid returns the arithmetic mean of the given points. |
| `func Circumcenter(a, b, c Point) (Point, error)` | Circumcenter returns the centre of the unique circle passing through a, b and c. |
| `func ClonePoints(pts []Point) []Point` | ClonePoints returns a fresh copy of the given slice. |
| `func ConvexHull(pts []Point) []Point` | ConvexHull returns the vertices of the convex hull of the given points in counterclockwise order, using Andrew's monotone-chain algorithm. |
| `func DedupePoints(pts []Point) []Point` | DedupePoints returns the points with exact duplicates removed, preserving the order of first appearance. |
| `func GrahamHullByAngle(pts []Point) []Point` | GrahamHullByAngle returns the convex hull computed by sorting points by polar angle about the lowest point. |
| `func Incenter(a, b, c Point) Point` | Incenter returns the incentre of the triangle (a, b, c): the centre of its inscribed circle, weighted by the opposite side lengths. |
| `func NewPoint(x, y float64) Point` | NewPoint returns the Point with the given coordinates. |
| `func Origin() Point` | Origin is the point (0, 0). |
| `func PolygonCentroid(poly []Point) Point` | PolygonCentroid returns the area centroid of the simple polygon. |
| `func SortPointsXY(pts []Point) []Point` | SortPointsXY returns a copy of pts sorted lexicographically by X then Y. |
| `func TriangleCentroid(a, b, c Point) Point` | TriangleCentroid returns the centroid (barycentre) of the triangle (a, b, c). |
| `func (p Point) Add(q Point) Point` | Add returns the componentwise sum p + q. |
| `func (p Point) Angle() float64` | Angle returns the angle of the vector p measured from the positive X axis, in the range (-pi, pi]. |
| `func (p Point) AngleTo(q Point) float64` | AngleTo returns the angle of the vector q - p from the positive X axis. |
| `func (p Point) ApproxEqual(q Point, eps float64) bool` | ApproxEqual reports whether p and q are within eps (per coordinate) of each other. |
| `func (p Point) Cross(q Point) float64` | Cross returns the scalar (2D) cross product p x q = p.X*q.Y - p.Y*q.X. |
| `func (p Point) Distance(q Point) float64` | Distance returns the Euclidean distance between p and q. |
| `func (p Point) DistanceSq(q Point) float64` | DistanceSq returns the squared Euclidean distance between p and q. |
| `func (p Point) Dot(q Point) float64` | Dot returns the dot product p . |
| `func (p Point) Equal(q Point) bool` | Equal reports whether p and q are exactly equal. |
| `func (p Point) IsFinite() bool` | IsFinite reports whether both coordinates of p are finite. |
| `func (p Point) Lerp(q Point, t float64) Point` | Lerp returns the linear interpolation (1-t)*p + t*q. |
| `func (p Point) Midpoint(q Point) Point` | Midpoint returns the point halfway between p and q. |
| `func (p Point) Neg() Point` | Neg returns the point reflected through the origin, -p. |
| `func (p Point) Norm() float64` | Norm returns the Euclidean norm \|p\|. |
| `func (p Point) NormSq() float64` | NormSq returns the squared Euclidean norm \|p\|^2. |
| `func (p Point) Normalize() Point` | Normalize returns the unit vector in the direction of p. |
| `func (p Point) Perp() Point` | Perp returns p rotated by +90 degrees, i.e. |
| `func (p Point) Rotate(theta float64) Point` | Rotate returns p rotated about the origin by theta radians (counterclockwise). |
| `func (p Point) RotateAbout(c Point, theta float64) Point` | RotateAbout returns p rotated by theta radians about the centre c. |
| `func (p Point) Scale(s float64) Point` | Scale returns p scaled by the factor s. |
| `func (p Point) String() string` | String returns a human-readable representation of p. |
| `func (p Point) Sub(q Point) Point` | Sub returns the componentwise difference p - q. |

</details>

<details>
<summary><code>Rect</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BoundingBox(pts []Point) (Rect, error)` | BoundingBox returns the axis-aligned bounding box of the given points. |
| `func NewRect(a, b Point) Rect` | NewRect returns the axis-aligned rectangle spanning the two corners, normalizing so that Min holds the smaller coordinates. |
| `func PolygonBoundingBox(poly []Point) (Rect, error)` | PolygonBoundingBox returns the axis-aligned bounding box of the polygon vertices. |
| `func (r Rect) Area() float64` | Area returns the area of the rectangle. |
| `func (r Rect) Center() Point` | Center returns the centre point of the rectangle. |
| `func (r Rect) Contains(p Point, eps float64) bool` | Contains reports whether p lies inside or on the rectangle, within tolerance eps. |
| `func (r Rect) Corners() [4]Point` | Corners returns the four corners of the rectangle in counterclockwise order starting from Min. |
| `func (r Rect) Diagonal() float64` | Diagonal returns the length of the rectangle's diagonal. |
| `func (r Rect) Expand(delta float64) Rect` | Expand returns the rectangle grown outward by delta on every side. |
| `func (r Rect) Height() float64` | Height returns the vertical extent of the rectangle. |
| `func (r Rect) Perimeter() float64` | Perimeter returns the perimeter of the rectangle. |
| `func (r Rect) String() string` | String returns a human-readable representation of the rectangle. |
| `func (r Rect) Union(s Rect) Rect` | Union returns the smallest rectangle containing both r and s. |
| `func (r Rect) Width() float64` | Width returns the horizontal extent of the rectangle. |

</details>

<details>
<summary><code>Triangle</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewTriangle(a, b, c int) Triangle` | NewTriangle returns the triangle on vertex indices a, b and c. |
| `func (t Triangle) Edges() [3]Edge` | Edges returns the triangle's three directed edges in order A->B, B->C, C->A. |
| `func (t Triangle) HasEdge(e Edge) bool` | HasEdge reports whether e is one of the triangle's edges (ignoring direction). |
| `func (t Triangle) HasVertex(v int) bool` | HasVertex reports whether v is one of the triangle's vertices. |
| `func (t Triangle) OppositeVertex(e Edge) int` | OppositeVertex returns the triangle vertex not incident to edge e. |
| `func (t Triangle) SharesEdge(u Triangle) (Edge, bool)` | SharesEdge reports whether triangles t and u share an edge, and returns that shared edge. |
| `func (t Triangle) Vertices() []int` | Vertices returns the triangle's three vertex indices as a slice. |

</details>

<details>
<summary><code>Triangulation</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DelaunayTriangulation(pts []Point) (*Triangulation, error)` | DelaunayTriangulation is an alias for Triangulate provided for readability at call sites that want the fully qualified name. |
| `func Triangulate(pts []Point) (*Triangulation, error)` | Triangulate computes the Delaunay triangulation of the given points using the incremental Bowyer-Watson algorithm. |
| `func (tr *Triangulation) Adjacency() [][]int` | Adjacency returns the undirected adjacency lists of the point graph induced by the triangulation edges. |
| `func (tr *Triangulation) BoundaryEdges() []Edge` | BoundaryEdges returns the edges that belong to exactly one triangle; these form the boundary of the triangulated region (the convex hull). |
| `func (tr *Triangulation) Circumcenter(i int) (Point, error)` | Circumcenter returns the circumcentre of triangle i. |
| `func (tr *Triangulation) Circumcenters() []Point` | Circumcenters returns the circumcentres of all triangles, in triangle order. |
| `func (tr *Triangulation) Circumcircle(i int) (Circle, error)` | Circumcircle returns the circumcircle of triangle i. |
| `func (tr *Triangulation) Contains(p Point) bool` | Contains reports whether p lies within the triangulated region (its convex hull). |
| `func (tr *Triangulation) ConvexHull() []Point` | ConvexHull returns the convex hull of the triangulation's points in counterclockwise order. |
| `func (tr *Triangulation) EdgeLengths() []float64` | EdgeLengths returns the length of every unique edge, in the order returned by Edges. |
| `func (tr *Triangulation) Edges() []Edge` | Edges returns the unique undirected edges of the triangulation, each in canonical (A <= B) form, sorted. |
| `func (tr *Triangulation) EulerCharacteristic() int` | EulerCharacteristic returns V - E + F for the triangulation, counting the unbounded face. |
| `func (tr *Triangulation) IsDelaunay() bool` | IsDelaunay reports whether the triangulation satisfies the empty-circumcircle property: no point lies strictly inside any triangle's circumcircle. |
| `func (tr *Triangulation) Locate(p Point) int` | Locate returns the index of a triangle that contains the query point p, or -1 if p lies outside the triangulated region. |
| `func (tr *Triangulation) LongestEdge() (Edge, float64, bool)` | LongestEdge returns the length of the longest edge and its endpoints as an edge. |
| `func (tr *Triangulation) Neighbors() [][]int` | Neighbors returns, for each triangle index, the indices of triangles sharing an edge with it. |
| `func (tr *Triangulation) NumEdges() int` | NumEdges returns the number of unique edges in the triangulation. |
| `func (tr *Triangulation) NumPoints() int` | NumPoints returns the number of points in the triangulation. |
| `func (tr *Triangulation) NumTriangles() int` | NumTriangles returns the number of triangles in the triangulation. |
| `func (tr *Triangulation) Point(i int) Point` | Point returns the i-th point of the triangulation. |
| `func (tr *Triangulation) Points() []Point` | Points returns the (deduplicated) points the triangulation was built on. |
| `func (tr *Triangulation) TotalArea() float64` | TotalArea returns the sum of the areas of all triangles, i.e. |
| `func (tr *Triangulation) Triangle(i int) Triangle` | Triangle returns the i-th triangle. |
| `func (tr *Triangulation) TriangleArea(i int) float64` | TriangleArea returns the area of triangle i. |
| `func (tr *Triangulation) TrianglePoints(i int) (a, b, c Point)` | TrianglePoints returns the three corner points of triangle i. |
| `func (tr *Triangulation) Triangles() []Triangle` | Triangles returns the triangulation's triangles as index triples. |
| `func (tr *Triangulation) VertexDegree(v int) int` | VertexDegree returns the number of distinct neighbours of point index v in the triangulation graph. |
| `func (tr *Triangulation) VertexTriangles() [][]int` | VertexTriangles returns, for each point index, the indices of the triangles incident to it. |

</details>

<details>
<summary><code>VoronoiCell</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (c VoronoiCell) Area() float64` | Area returns the area of the cell. |
| `func (c VoronoiCell) Centroid() Point` | Centroid returns the centroid of a bounded cell. |
| `func (c VoronoiCell) NumVertices() int` | NumVertices returns the number of finite vertices of the cell. |
| `func (c VoronoiCell) Perimeter() float64` | Perimeter returns the perimeter of the finite portion of the cell boundary. |

</details>

<details>
<summary><code>VoronoiDiagram</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Voronoi(sites []Point) (*VoronoiDiagram, error)` | Voronoi builds the Voronoi diagram of the given sites by first constructing their Delaunay triangulation. |
| `func VoronoiFromTriangulation(tri *Triangulation) *VoronoiDiagram` | VoronoiFromTriangulation constructs the Voronoi diagram dual to an existing Delaunay triangulation, avoiding a second triangulation pass. |
| `func (vd *VoronoiDiagram) BoundedCells() []VoronoiCell` | BoundedCells returns only the cells that are finite (bounded) polygons. |
| `func (vd *VoronoiDiagram) Cell(s int) VoronoiCell` | Cell returns the Voronoi cell of site index s. |
| `func (vd *VoronoiDiagram) Cells() []VoronoiCell` | Cells returns all Voronoi cells, indexed by site. |
| `func (vd *VoronoiDiagram) DelaunayEdges() []Edge` | DelaunayEdges returns the edges of the dual Delaunay triangulation. |
| `func (vd *VoronoiDiagram) Edges() []VoronoiEdge` | Edges returns the finite Voronoi edges. |
| `func (vd *VoronoiDiagram) NearestSite(p Point) (int, float64)` | NearestSite returns the index of the site closest to p and the distance to it. |
| `func (vd *VoronoiDiagram) NumEdges() int` | NumEdges returns the number of finite Voronoi edges. |
| `func (vd *VoronoiDiagram) NumSites() int` | NumSites returns the number of sites. |
| `func (vd *VoronoiDiagram) NumVertices() int` | NumVertices returns the number of finite Voronoi vertices. |
| `func (vd *VoronoiDiagram) Sites() []Point` | Sites returns the diagram's sites (the deduplicated input points). |
| `func (vd *VoronoiDiagram) Triangulation() *Triangulation` | Triangulation returns the underlying Delaunay triangulation. |
| `func (vd *VoronoiDiagram) Vertices() []Point` | Vertices returns the finite Voronoi vertices (the Delaunay circumcentres), indexed by triangle. |

</details>

<details>
<summary><code>VoronoiEdge</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (e VoronoiEdge) Length() float64` | Length returns the Euclidean length of the Voronoi edge. |
| `func (e VoronoiEdge) Midpoint() Point` | Midpoint returns the midpoint of the Voronoi edge. |

</details>

### Constants

`Eps`

### Variables

`ErrCollinear`, `ErrDegenerate`, `ErrEmpty`, `ErrTooFewPoints`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/voronoi).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
