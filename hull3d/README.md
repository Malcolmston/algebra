# hull3d — implements three-dimensional computational geometry over float64 coordinates using only

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/hull3d.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/hull3d)

Package hull3d implements three-dimensional computational geometry over
float64 coordinates using only the Go standard library.

- Convex hulls of finite point clouds via an incremental / QuickHull style
construction (`ConvexHull`, `QuickHull`) with a gift-wrapping variant
(`GiftWrapHull`) provided as an independent cross-check. - Delaunay
tetrahedralization (`Delaunay`) using the incremental Bowyer–Watson
algorithm and its dual, the 3-D Voronoi diagram (`VoronoiCells`). - Half-space
intersection (`HalfSpaceIntersection`) computed as the polar dual of a convex
hull. - Minkowski sums and differences of convex polytopes (`MinkowskiSum`,
`MinkowskiDifference`). - Collision and proximity queries between convex
bodies via the Gilbert–Johnson–Keerthi distance algorithm (`GJKDistance`,
`GJKIntersect`) and Expanding Polytope Algorithm penetration depth
(`EPAPenetration`). - Metric and combinatorial queries on polytopes: volume,
surface area and centroid (`Polytope.Volume`, `Polytope.SurfaceArea`,
`Polytope.Centroid`), face / edge / vertex enumeration, point containment, and
supporting and separating planes.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/hull3d"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AreDisjoint(a, b *Polytope, eps float64) bool` | AreDisjoint reports whether two convex polytopes are disjoint (do not overlap) within tolerance eps, via their GJK distance. |
| `func BarycentricTetrahedron(q, a, b, c, d Vec3) ([4]float64, error)` | BarycentricTetrahedron returns the four barycentric coordinates of q with respect to the tetrahedron (a,b,c,d). |
| `func Collinear(a, b, c Vec3, eps float64) bool` | Collinear reports whether the three points a, b, c are collinear to within tolerance eps on the squared area of the triangle they span. |
| `func Coplanar(a, b, c, d Vec3, eps float64) bool` | Coplanar reports whether the four points a, b, c, d are coplanar to within tolerance eps on the orientation determinant. |
| `func DistancePointToPolytope(p *Polytope, q Vec3) (float64, error)` | DistancePointToPolytope returns the Euclidean distance from q to the polytope (zero if q is inside). |
| `func GJKIntersect(a, b ConvexShape) bool` | GJKIntersect reports whether the two convex shapes overlap, using the Gilbert–Johnson–Keerthi algorithm to test whether the origin lies in their… |
| `func HullSurfaceArea(points []Vec3) (float64, error)` | HullSurfaceArea computes the surface area of the convex hull of the given points. |
| `func HullVolume(points []Vec3) (float64, error)` | HullVolume computes the volume of the convex hull of the given points. |
| `func InSphere(a, b, c, d, e Vec3) float64` | InSphere returns a value whose sign tells whether point e lies inside the sphere through a, b, c, d. |
| `func IsConvexPosition(p Vec3, pts []Vec3, eps float64) bool` | IsConvexPosition reports whether p is a vertex (extreme point) of the convex hull of pts, i.e. |
| `func Orient3D(a, b, c, d Vec3) float64` | Orient3D returns a value whose sign tells the orientation of point d relative to the plane through a, b, c: positive when d lies on the positive side… |
| `func OrientSign(a, b, c, d Vec3, eps float64) int` | OrientSign returns the sign of `Orient3D` with a coplanarity tolerance eps: +1, -1, or 0 when \|Orient3D\| <= eps. |
| `func PointInTetrahedron(q, a, b, c, d Vec3, eps float64) bool` | PointInTetrahedron reports whether q lies inside (or on) the tetrahedron with the given corners, within tolerance eps on the orientation determinant. |
| `func PolytopesOverlap(a, b *Polytope, eps float64) (bool, error)` | PolytopesOverlap reports whether two convex polytopes overlap, by testing whether their Minkowski difference contains the origin. |
| `func RayPolytopeIntersection(p *Polytope, origin, dir Vec3) (tEnter, tExit float64, ok bool)` | RayPolytopeIntersection computes the entry and exit parameters of the ray origin + t*dir with the convex polytope, using the slab method against… |
| `func SegmentPolytopeIntersect(p *Polytope, s0, s1 Vec3) bool` | SegmentPolytopeIntersect reports whether the segment from s0 to s1 meets the convex polytope. |
| `func SphereSurfaceArea(radius float64) float64` | SphereSurfaceArea returns the surface area of a sphere of the given radius. |
| `func SphereVolume(radius float64) float64` | SphereVolume returns the volume of a ball of the given radius. |
| `func TetrahedronVolume(a, b, c, d Vec3) float64` | TetrahedronVolume returns the unsigned volume of the tetrahedron with the given four corners. |
| `func TriangleArea(a, b, c Vec3) float64` | TriangleArea returns the area of the triangle with the given three corners. |
| `func Triple(a, b, c Vec3) float64` | Triple returns the scalar triple product a·(b×c), the signed volume of the parallelepiped spanned by a, b and c. |

### Types

| Type | What it is |
| --- | --- |
| `BoxShape` | BoxShape is an axis-aligned box with the given centre and half-extents. |
| `CapsuleShape` | CapsuleShape is a capsule (a swept sphere): the set of points within Radius of the segment from A to B. It is the Minkowski sum of that segment and a… |
| `ConvexShape` | ConvexShape is any convex body that can report its support point: the point of the body farthest in a given direction. |
| `Edge` | Edge is an undirected edge between two vertex indices, stored with Lo <= Hi so that it can be used as a canonical map key. |
| `Face` | Face is a triangular face of a `Polytope`, given as three indices into the polytope's vertex slice. |
| `HalfSpace` | HalfSpace is a closed half-space, the feasible side of a `Plane`. |
| `MinkowskiSumShape` | MinkowskiSumShape is the Minkowski sum of two convex shapes, whose support in direction d is the sum of the two supports in d. |
| `Penetration` | Penetration describes the result of an `EPAPenetration` query for two overlapping convex shapes. |
| `Plane` | Plane represents an oriented plane in R^3 as the set of points p satisfying Normal·p = Offset. |
| `PointCloudShape` | PointCloudShape is the convex hull of a finite set of points, given only by the points themselves (no explicit hull needed for support queries). |
| `PointShape` | PointShape is a single point treated as a degenerate convex body. |
| `Polytope` | Polytope is a boundary representation of a (typically convex) polyhedron: a list of vertices together with the oriented triangular faces of its… |
| `PolytopeShape` | PolytopeShape adapts a `Polytope` to `ConvexShape` using its vertices. |
| `SegmentShape` | SegmentShape is the line segment between A and B. |
| `SphereShape` | SphereShape is a ball with the given centre and radius. |
| `Tetrahedralization` | Tetrahedralization is a decomposition of the convex hull of a point set into tetrahedra, produced by `Delaunay`. |
| `Tetrahedron` | Tetrahedron is a solid tetrahedron given as four indices into a point slice. |
| `TranslatedShape` | TranslatedShape is a convex shape shifted by Offset. |
| `Vec3` | Vec3 is a three-dimensional vector over float64, used interchangeably as a point or a displacement in R^3. |
| `VoronoiCell` | VoronoiCell is the Voronoi region of a single site: the set of points closer to Site than to any other site, represented by its bounding half-spaces… |
| `VoronoiDiagram` | VoronoiDiagram is the 3-D Voronoi diagram of a set of sites, the geometric dual of the Delaunay tetrahedralization. |

<details>
<summary><code>BoxShape</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (s BoxShape) Support(d Vec3) Vec3` | Support implements `ConvexShape`. |

</details>

<details>
<summary><code>CapsuleShape</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (s CapsuleShape) Support(d Vec3) Vec3` | Support implements `ConvexShape`. |

</details>

<details>
<summary><code>Edge</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func MakeEdge(i, j int) Edge` | MakeEdge returns the canonical `Edge` joining vertex indices i and j. |

</details>

<details>
<summary><code>HalfSpace</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func HalfSpaceFromPointNormal(p, normal Vec3) HalfSpace` | HalfSpaceFromPointNormal returns the half-space with the given boundary point and outward normal: {x : normal·(x-p) <= 0}. |
| `func NewHalfSpace(normal Vec3, offset float64) HalfSpace` | NewHalfSpace returns the half-space {x : normal·x <= offset}. |
| `func (h HalfSpace) Contains(q Vec3, eps float64) bool` | Contains reports whether q lies in the (closed, eps-relaxed) half-space. |

</details>

<details>
<summary><code>MinkowskiSumShape</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (s MinkowskiSumShape) Support(d Vec3) Vec3` | Support implements `ConvexShape`. |

</details>

<details>
<summary><code>Penetration</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func EPAPenetration(a, b ConvexShape) (Penetration, error)` | EPAPenetration computes the penetration depth and contact normal of two overlapping convex shapes with the Expanding Polytope Algorithm. |

</details>

<details>
<summary><code>Plane</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewPlane(normal Vec3, offset float64) Plane` | NewPlane returns the plane with the given normal and offset. |
| `func PlaneFromPointNormal(p, normal Vec3) Plane` | PlaneFromPointNormal returns the plane through point p with the given normal. |
| `func PlaneFromPoints(a, b, c Vec3) (Plane, error)` | PlaneFromPoints returns the plane through the three points a, b, c, oriented so that its normal is (b-a)×(c-a). |
| `func SeparatingPlane(a, b *Polytope, eps float64) (Plane, error)` | SeparatingPlane returns a plane strictly separating two disjoint convex polytopes, with a and b on the negative and positive sides respectively. |
| `func SupportingPlane(a, b, c, interior Vec3) (Plane, error)` | SupportingPlane returns the plane through the triangle (a,b,c) oriented so that its positive half-space does not contain the interior reference point… |
| `func (p Plane) Contains(q Vec3, eps float64) bool` | Contains reports whether q lies on the plane to within tolerance eps. |
| `func (p Plane) Distance(q Vec3) float64` | Distance returns the unsigned Euclidean distance from q to the plane. |
| `func (p Plane) Eval(q Vec3) float64` | Eval returns Normal·q - Offset, which is positive for q in the positive half-space, zero on the plane, and negative below. |
| `func (p Plane) Flip() Plane` | Flip returns the plane describing the same point set with the opposite orientation (negated normal and offset). |
| `func (p Plane) Normalized() (Plane, error)` | Normalized returns an equivalent plane whose normal has unit length, so that SignedDistance reports true Euclidean distance. |
| `func (p Plane) PointOnPlane() (Vec3, error)` | PointOnPlane returns some point lying on the plane (the foot of the perpendicular from the origin). |
| `func (p Plane) Project(q Vec3) (Vec3, error)` | Project returns the orthogonal projection of q onto the plane. |
| `func (p Plane) Side(q Vec3, eps float64) int` | Side returns +1 if q lies in the positive half-space beyond tolerance eps, -1 if it lies in the negative half-space, and 0 if it is within eps of the… |
| `func (p Plane) SignedDistance(q Vec3) float64` | SignedDistance returns the signed distance from q to the plane, positive on the side the normal points to. |
| `func (p Plane) String() string` | String renders the plane's equation. |

</details>

<details>
<summary><code>PointCloudShape</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (s PointCloudShape) Support(d Vec3) Vec3` | Support implements `ConvexShape`. |

</details>

<details>
<summary><code>PointShape</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (s PointShape) Support(Vec3) Vec3` | Support implements `ConvexShape`. |

</details>

<details>
<summary><code>Polytope</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BoxPolytope(center, half Vec3) *Polytope` | BoxPolytope returns the axis-aligned box centred at center with the given half-extents as a closed triangulated `Polytope`. |
| `func ConvexHull(points []Vec3) (*Polytope, error)` | ConvexHull computes the convex hull of the given points and returns it as a closed, outward-oriented triangulated `Polytope`. |
| `func ConvexHullSeeded(points []Vec3, seed int64) (*Polytope, error)` | ConvexHullSeeded behaves like `ConvexHull` but shuffles the input point order using the supplied seed before construction. |
| `func Dodecahedron(r float64) *Polytope` | Dodecahedron returns the regular dodecahedron with the given circumradius centred at the origin. |
| `func GiftWrapHull(points []Vec3) (*Polytope, error)` | GiftWrapHull computes the convex hull by enumerating supporting facets and returns it as a closed, outward-oriented triangulated `Polytope`. |
| `func HalfSpaceIntersection(halfspaces []HalfSpace, interior Vec3) (*Polytope, error)` | HalfSpaceIntersection computes the bounded convex polytope that is the intersection of the given half-spaces, using the polar-dual construction: with… |
| `func Icosahedron(r float64) *Polytope` | Icosahedron returns the regular icosahedron with the given circumradius centred at the origin. |
| `func IntersectHalfSpacesWithBox(halfspaces []HalfSpace, interior Vec3, extent float64) (*Polytope, error)` | IntersectHalfSpacesWithBox intersects the given half-spaces together with the six faces of a large axis-aligned bounding box centred at interior with… |
| `func MinkowskiDifference(a, b *Polytope) (*Polytope, error)` | MinkowskiDifference returns the convex polytope A ⊕ (-B), the convex hull of all pairwise vertex differences. |
| `func MinkowskiSum(a, b *Polytope) (*Polytope, error)` | MinkowskiSum returns the convex polytope that is the Minkowski sum of two convex polytopes: the convex hull of all pairwise vertex sums. |
| `func NewPolytope(vertices []Vec3, faces []Face) *Polytope` | NewPolytope returns a polytope with the given vertices and faces. |
| `func Octahedron(r float64) *Polytope` | Octahedron returns the regular octahedron with the given circumradius centred at the origin (vertices at ±r along each axis). |
| `func QuickHull(points []Vec3) (*Polytope, error)` | QuickHull is an alias for `ConvexHull`; the incremental construction used here selects, for each face, the farthest outside point to add next, which… |
| `func RegularTetrahedron(r float64) *Polytope` | RegularTetrahedron returns a regular tetrahedron with the given circumradius centred at the origin. |
| `func ScaledSum(a, b *Polytope, s, t float64) (*Polytope, error)` | ScaledSum returns the convex hull of s*a + t*b over vertices, a weighted Minkowski combination useful for morphing between convex bodies. |
| `func TetrahedronPolytope(a, b, c, d Vec3) (*Polytope, error)` | TetrahedronPolytope returns the tetrahedron with the four given corners as a closed `Polytope`. |
| `func UnitCube() *Polytope` | UnitCube returns the axis-aligned cube of side 1 centred at the origin. |
| `func (p *Polytope) BoundingBox() (min, max Vec3, err error)` | BoundingBox returns the axis-aligned bounding box of the vertices. |
| `func (p *Polytope) Centroid() Vec3` | Centroid returns the centroid of the solid region enclosed by the surface, computed by integrating position over the tetrahedra formed with the… |
| `func (p *Polytope) Contains(q Vec3, eps float64) bool` | Contains reports whether point q lies inside or on the polytope, treating it as the intersection of the half-spaces below each face plane. |
| `func (p *Polytope) Copy() *Polytope` | Copy returns a deep copy of the polytope. |
| `func (p *Polytope) Diameter() float64` | Diameter returns the maximum distance between any two vertices (the geometric diameter of the vertex set). |
| `func (p *Polytope) Edges() []Edge` | Edges returns the distinct undirected edges of the polytope's surface. |
| `func (p *Polytope) EulerCharacteristic() int` | EulerCharacteristic returns V - E + F. For a simple closed polyhedron (topologically a sphere) it equals 2. |
| `func (p *Polytope) FaceArea(f Face) float64` | FaceArea returns the area of triangular face f. |
| `func (p *Polytope) FaceCentroid(f Face) Vec3` | FaceCentroid returns the centroid (barycentre) of triangular face f. |
| `func (p *Polytope) FaceNormal(f Face) Vec3` | FaceNormal returns the (unnormalised) outward normal (B-A)×(C-A) of face f. |
| `func (p *Polytope) FacePlane(f Face) (Plane, error)` | FacePlane returns the plane of face f, oriented so its normal is the outward normal. |
| `func (p *Polytope) FacePlanes() []Plane` | FacePlanes returns the outward-oriented supporting plane of every face, skipping degenerate faces. |
| `func (p *Polytope) FaceVertices(f Face) (a, b, c Vec3)` | FaceVertices returns the three vertex positions of face f. |
| `func (p *Polytope) HalfSpaces() []HalfSpace` | HalfSpaces returns one outward half-space per face of the polytope, so that the polytope equals the intersection of the returned half-spaces (for a… |
| `func (p *Polytope) IsClosed() bool` | IsClosed reports whether every undirected edge is shared by exactly two faces, a necessary condition for the surface to bound a solid. |
| `func (p *Polytope) IsConvex(eps float64) bool` | IsConvex reports whether the polytope is convex and consistently outward-oriented: every vertex lies on the inner (non-positive) side of every face's… |
| `func (p *Polytope) NumEdges() int` | NumEdges returns the number of distinct undirected edges. |
| `func (p *Polytope) NumFaces() int` | NumFaces returns the number of triangular faces. |
| `func (p *Polytope) NumVertices() int` | NumVertices returns the number of vertices. |
| `func (p *Polytope) OrientOutward() *Polytope` | OrientOutward re-orders each face so its normal points away from the polytope centroid, making the surface consistently outward-oriented. |
| `func (p *Polytope) Reflect() *Polytope` | Reflect returns a copy of the polytope reflected through the origin (each vertex negated), with face windings reversed so it stays outward-oriented. |
| `func (p *Polytope) Scaled(s float64) *Polytope` | Scaled returns a copy of the polytope scaled about the origin by factor s. |
| `func (p *Polytope) SignedVolume() float64` | SignedVolume returns the signed volume enclosed by the surface via the divergence theorem. |
| `func (p *Polytope) Support(d Vec3) (Vec3, error)` | Support returns the vertex of the polytope farthest in direction d (a support point of the convex hull of its vertices). |
| `func (p *Polytope) SupportingPlaneOf(d Vec3) (Plane, error)` | SupportingPlaneOf returns the supporting plane of the polytope in direction d: the plane with outward normal d touching the polytope at its farthest… |
| `func (p *Polytope) SurfaceArea() float64` | SurfaceArea returns the total surface area of the polytope, the sum of its face areas. |
| `func (p *Polytope) Translate(t Vec3) *Polytope` | Translate returns a copy of the polytope with every vertex shifted by t. |
| `func (p *Polytope) Validate() error` | Validate reports an error if any face references an out-of-range vertex index or is degenerate (repeats a vertex). |
| `func (p *Polytope) VertexCentroid() Vec3` | VertexCentroid returns the arithmetic mean of the vertices. |
| `func (p *Polytope) Volume() float64` | Volume returns the (unsigned) volume enclosed by the surface. |
| `func (p *Polytope) Width(d Vec3) (float64, error)` | Width returns the extent of the polytope along direction d: the difference between its maximum and minimum support in d. |

</details>

<details>
<summary><code>PolytopeShape</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (s PolytopeShape) Support(d Vec3) Vec3` | Support implements `ConvexShape`. |

</details>

<details>
<summary><code>SegmentShape</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (s SegmentShape) Support(d Vec3) Vec3` | Support implements `ConvexShape`. |

</details>

<details>
<summary><code>SphereShape</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Circumsphere(a, b, c, d Vec3) (SphereShape, error)` | Circumsphere returns the circumscribed sphere of the tetrahedron with the given corners as a `SphereShape`. |
| `func (s SphereShape) Support(d Vec3) Vec3` | Support implements `ConvexShape`. |

</details>

<details>
<summary><code>Tetrahedralization</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Delaunay(points []Vec3) (*Tetrahedralization, error)` | Delaunay computes a Delaunay tetrahedralization of the given points using the incremental Bowyer–Watson algorithm. |
| `func (d *Tetrahedralization) BoundaryFaces() [][3]int` | BoundaryFaces returns the triangular faces that belong to exactly one tetrahedron; collectively they form the surface of the convex hull. |
| `func (d *Tetrahedralization) Faces() [][3]int` | Faces returns the distinct triangular faces of the tetrahedralization, each as a sorted index triple. |
| `func (d *Tetrahedralization) IsDelaunay(eps float64) bool` | IsDelaunay reports whether the tetrahedralization satisfies the empty-sphere property: no point lies strictly inside the circumsphere of any… |
| `func (d *Tetrahedralization) NumTets() int` | NumTets returns the number of tetrahedra. |
| `func (d *Tetrahedralization) TetVertices(t Tetrahedron) (a, b, c, e Vec3)` | TetVertices returns the four corner positions of tetrahedron t. |
| `func (d *Tetrahedralization) TotalVolume() float64` | TotalVolume returns the sum of the volumes of all tetrahedra, which equals the volume of the convex hull of the point set. |

</details>

<details>
<summary><code>TranslatedShape</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (s TranslatedShape) Support(d Vec3) Vec3` | Support implements `ConvexShape`. |

</details>

<details>
<summary><code>Vec3</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BoundingBox(pts []Vec3) (min, max Vec3, err error)` | BoundingBox returns the axis-aligned minimum and maximum corners enclosing the given points. |
| `func Centroid(pts []Vec3) Vec3` | Centroid returns the arithmetic mean of the given points. |
| `func ChebyshevCenterEstimate(halfspaces []HalfSpace) (center Vec3, radius float64, err error)` | ChebyshevCenterEstimate returns an estimate of the deepest interior point and its clearance for the intersection of the half-spaces, by hill-climbing… |
| `func Circumcenter(a, b, c, d Vec3) (center Vec3, radius float64, err error)` | Circumcenter returns the centre of the sphere passing through the four points and its radius. |
| `func ClosestPointInPolytope(p *Polytope, q Vec3) (Vec3, float64, error)` | ClosestPointInPolytope returns the point of the (convex) polytope closest to q and the distance to it, using GJK between the polytope and the single… |
| `func ClosestPointOnTriangle(q, a, b, c Vec3) Vec3` | ClosestPointOnTriangle returns the point on triangle (a,b,c) closest to q. |
| `func DedupPoints(pts []Vec3, eps float64) []Vec3` | DedupPoints returns the points of pts with near-duplicates (within Euclidean distance eps) removed, preserving first-seen order. |
| `func ExtremePoints(points []Vec3) ([]Vec3, error)` | ExtremePoints returns the subset of input points that are vertices of the convex hull. |
| `func FarthestPoint(pts []Vec3, d Vec3) (int, Vec3, error)` | FarthestPoint returns the index and value of the point in pts that maximises the dot product with direction d (the support point in direction d). |
| `func FeasiblePoint(halfspaces []HalfSpace) (Vec3, error)` | FeasiblePoint returns a strictly interior point of the intersection of the half-spaces, if one is found by the Chebyshev-centre heuristic of… |
| `func GJKDistance(a, b ConvexShape) (dist float64, ptA, ptB Vec3, err error)` | GJKDistance returns the Euclidean distance between two convex shapes together with the closest points on each. |
| `func HullCentroid(points []Vec3) (Vec3, error)` | HullCentroid computes the centroid of the solid convex hull of the given points, that is the centre of mass of the filled body rather than of the… |
| `func LineSegmentClosest(p0, p1, q0, q1 Vec3) (cp, cq Vec3, dist float64)` | LineSegmentClosest returns the closest points between two segments (p0,p1) and (q0,q1) and the distance between them. |
| `func MinkowskiDiffPoints(as, bs []Vec3) []Vec3` | MinkowskiDiffPoints returns the pointwise differences a-b for a in as and b in bs. |
| `func MinkowskiSumPoints(as, bs []Vec3) []Vec3` | MinkowskiSumPoints returns the Minkowski sum of two point sets: every pairwise sum a+b for a in as and b in bs. |
| `func NewVec3(x, y, z float64) Vec3` | NewVec3 returns the vector with the given components. |
| `func SamplePointsInBall(n int, center Vec3, radius float64, seed int64) []Vec3` | SamplePointsInBall returns n points sampled uniformly inside the ball of the given radius centred at center, using the supplied seed for… |
| `func SamplePointsOnSphere(n int, center Vec3, radius float64, seed int64) []Vec3` | SamplePointsOnSphere returns n points sampled uniformly on the sphere of the given radius centred at center, using the supplied seed. |
| `func SortByAxis(pts []Vec3, axis int) []Vec3` | SortByAxis returns a copy of pts sorted ascending by component axis (0, 1 or 2). |
| `func Splat(s float64) Vec3` | Splat returns the vector with every component equal to s. |
| `func SumVec3(pts []Vec3) Vec3` | SumVec3 returns the sum of the given points. |
| `func UnitX() Vec3` | UnitX returns the first standard basis vector (1, 0, 0). |
| `func UnitY() Vec3` | UnitY returns the second standard basis vector (0, 1, 0). |
| `func UnitZ() Vec3` | UnitZ returns the third standard basis vector (0, 0, 1). |
| `func Zero() Vec3` | Zero returns the zero vector (0, 0, 0). |
| `func (a Vec3) Abs() Vec3` | Abs returns the component-wise absolute value of a. |
| `func (a Vec3) Add(b Vec3) Vec3` | Add returns the component-wise sum a + b. |
| `func (a Vec3) AngleBetween(b Vec3) (float64, error)` | AngleBetween returns the unsigned angle in radians between a and b, in the range [0, pi]. |
| `func (a Vec3) AnyPerpendicular() (Vec3, error)` | AnyPerpendicular returns a non-zero vector orthogonal to a. |
| `func (a Vec3) ApproxEqual(b Vec3, eps float64) bool` | ApproxEqual reports whether a and b agree to within an absolute tolerance eps in every component. |
| `func (a Vec3) Clamp(lo, hi Vec3) Vec3` | Clamp returns a with each component clamped to the box [lo, hi]. |
| `func (a Vec3) Cross(b Vec3) Vec3` | Cross returns the vector cross product a×b. |
| `func (a Vec3) Distance(b Vec3) float64` | Distance returns the Euclidean distance between a and b. |
| `func (a Vec3) DistanceSq(b Vec3) float64` | DistanceSq returns the squared Euclidean distance between a and b. |
| `func (a Vec3) Div(b Vec3) (Vec3, error)` | Div returns the component-wise quotient a / b. |
| `func (a Vec3) Dot(b Vec3) float64` | Dot returns the Euclidean inner product a·b. |
| `func (a Vec3) Equal(b Vec3) bool` | Equal reports whether a and b are exactly equal component-wise. |
| `func (a Vec3) Get(i int) float64` | Get returns component i of a for i in {0,1,2}. |
| `func (a Vec3) IsFinite() bool` | IsFinite reports whether all components of a are finite (neither Inf nor NaN). |
| `func (a Vec3) IsZero() bool` | IsZero reports whether every component of a is exactly zero. |
| `func (a Vec3) L1() float64` | L1 returns the Manhattan (taxicab) norm \|x\|+\|y\|+\|z\|. |
| `func (a Vec3) LInf() float64` | LInf returns the Chebyshev (maximum) norm max(\|x\|,\|y\|,\|z\|). |
| `func (a Vec3) Length() float64` | Length returns the Euclidean length of a. |
| `func (a Vec3) LengthSq() float64` | LengthSq returns the squared Euclidean length of a. |
| `func (a Vec3) Lerp(b Vec3, t float64) Vec3` | Lerp returns the linear interpolation (1-t)*a + t*b. |
| `func (a Vec3) Max(b Vec3) Vec3` | Max returns the component-wise maximum of a and b. |
| `func (a Vec3) MaxAxis() int` | MaxAxis returns the index (0, 1 or 2) of the component of a with the largest absolute value. |
| `func (a Vec3) MaxComponent() float64` | MaxComponent returns the largest of the three components of a. |
| `func (a Vec3) Midpoint(b Vec3) Vec3` | Midpoint returns the arithmetic mean of a and b. |
| `func (a Vec3) Min(b Vec3) Vec3` | Min returns the component-wise minimum of a and b. |
| `func (a Vec3) MinComponent() float64` | MinComponent returns the smallest of the three components of a. |
| `func (a Vec3) Mul(b Vec3) Vec3` | Mul returns the component-wise (Hadamard) product of a and b. |
| `func (a Vec3) Neg() Vec3` | Neg returns the additive inverse -a. |
| `func (a Vec3) Norm() float64` | Norm is an alias for Length. |
| `func (a Vec3) Normalize() (Vec3, error)` | Normalize returns a unit vector in the direction of a. |
| `func (a Vec3) NormalizeOr(fallback Vec3) Vec3` | NormalizeOr returns a unit vector in the direction of a, or the fallback vector if a has effectively zero length. |
| `func (a Vec3) OrthonormalBasis() (u, v Vec3, err error)` | OrthonormalBasis returns two unit vectors u, v that together with a normalized copy of a form a right-handed orthonormal basis. |
| `func (a Vec3) ProjectOnto(b Vec3) (Vec3, error)` | ProjectOnto returns the orthogonal projection of a onto the direction b. |
| `func (a Vec3) Reflect(n Vec3) Vec3` | Reflect returns a reflected across the plane through the origin with unit normal n. |
| `func (a Vec3) RejectFrom(b Vec3) (Vec3, error)` | RejectFrom returns the component of a orthogonal to b (a minus its projection onto b). |
| `func (a Vec3) Scale(s float64) Vec3` | Scale returns the scalar multiple s*a. |
| `func (a Vec3) String() string` | String renders a in a compact parenthesised form. |
| `func (a Vec3) Sub(b Vec3) Vec3` | Sub returns the component-wise difference a - b. |
| `func (a Vec3) With(i int, v float64) Vec3` | With returns a copy of a with component i replaced by v. |

</details>

<details>
<summary><code>VoronoiCell</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (c VoronoiCell) Contains(q Vec3, eps float64) bool` | Contains reports whether q lies in the Voronoi cell of its site, i.e. |

</details>

<details>
<summary><code>VoronoiDiagram</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func VoronoiCells(sites []Vec3) (*VoronoiDiagram, error)` | VoronoiCells builds the Voronoi diagram of the given sites as the dual of the Delaunay tetrahedralization: each interior Delaunay vertex maps to a… |
| `func (v *VoronoiDiagram) Cell(i int) (VoronoiCell, error)` | Cell returns the Voronoi cell for site index i. |
| `func (v *VoronoiDiagram) NearestSite(q Vec3) (int, error)` | NearestSite returns the index of the site nearest to q by brute-force search, which is the Voronoi cell containing q. |
| `func (v *VoronoiDiagram) NumBoundedCells() int` | NumBoundedCells returns the number of bounded (interior) Voronoi cells. |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/hull3d).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
