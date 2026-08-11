# geom2d — provides primitives for 2D computational geometry

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/geom2d.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/geom2d)

Package geom2d provides primitives for 2D computational geometry.

It defines the value types `Point2` (a location in the plane) and `Vec2` (a
displacement), together with the supporting `Line` and `Segment` types. On top
of these it offers the standard toolbox of planar geometry: vector arithmetic
(add, subtract, scale, dot, cross), norms and angles, rotation and reflection,
projection, linear interpolation, point/point, point/line and point/segment
distances, closest-point queries, orientation and collinearity predicates, and
line and segment intersection.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/geom2d"
```

## Usage

This is the package's own `ExampleCircleCircleIntersection`, so it compiles and its output is
asserted on every `go test ./geom2d/`.

```go
c1 := geom2d.NewCircle(geom2d.NewPoint2(0, 0), 1)
	c2 := geom2d.NewCircle(geom2d.NewPoint2(1, 0), 1)

	pts, ok := geom2d.CircleCircleIntersection(c1, c2)
	fmt.Println(len(pts), ok)
	for _, p := range pts {
		fmt.Printf("(%.4f, %.4f)\n", p.X, p.Y)
	}

	_, ok = geom2d.CircleCircleIntersection(c1, c1)
	fmt.Println(ok)
```

```
2 true
(0.5000, 0.8660)
(0.5000, -0.8660)
false
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AngleAt(a, b, c Point2) float64` | AngleAt returns the interior angle in [0, π] at vertex b of the path a-b-c, that is the angle between the vectors b→a and b→c. |
| `func Barycentric(p, a, b, c Point2) (u, v, w float64, ok bool)` | Barycentric returns the barycentric coordinates of point p with respect to the triangle a, b, c. |
| `func CCW(a, b, c Point2, eps float64) bool` | CCW reports whether the points a, b, c are in strict counter-clockwise order, using eps to reject near-collinear triples. |
| `func Collinear(a, b, c Point2, eps float64) bool` | Collinear reports whether the points a, b, c lie on a common line, judged by twice the triangle area being at most eps in magnitude. |
| `func ConvexDiameter(pts []Point2) float64` | ConvexDiameter returns the diameter (maximum pairwise distance) of the input point set, computed with rotating calipers over the convex hull. |
| `func ConvexDiameterPoints(pts []Point2) (Point2, Point2, float64)` | ConvexDiameterPoints returns the two points of the input set that are farthest apart (the diameter of the point set) together with that distance,… |
| `func Cross(u, v Vec2) float64` | Cross returns the 2D cross product u×v. |
| `func Distance(p, q Point2) float64` | Distance returns the Euclidean distance between points p and q. |
| `func DistanceSq(p, q Point2) float64` | DistanceSq returns the squared Euclidean distance between points p and q, avoiding the square root. |
| `func Dot(u, v Vec2) float64` | Dot returns the dot product u·v. |
| `func IsConvexPolygon(pts []Point2) bool` | IsConvexPolygon reports whether the closed polygon pts is convex. |
| `func OnPolygonBoundary(p Point2, pts []Point2, eps float64) bool` | OnPolygonBoundary reports whether p lies within eps of any edge of the closed polygon pts. |
| `func OnSegment(p Point2, s Segment, eps float64) bool` | OnSegment reports whether point p, already assumed to be collinear with the segment's endpoints within eps, lies within the segment's extent. |
| `func Orientation(a, b, c Point2, eps float64) int` | Orientation returns +1 if the points a, b, c make a counter-clockwise turn, -1 if they make a clockwise turn, and 0 if they are collinear, judged… |
| `func PointInPolygon(p Point2, pts []Point2) bool` | PointInPolygon reports whether point p lies inside the closed polygon pts using the ray-casting (even-odd / crossing-number) rule. |
| `func PointInPolygonWinding(p Point2, pts []Point2) bool` | PointInPolygonWinding reports whether p lies inside the closed polygon pts using the non-zero winding rule (the winding number is non-zero). |
| `func PointLineDistance(p Point2, l Line) float64` | PointLineDistance returns the perpendicular distance from p to the infinite line l. |
| `func PointSegmentDistance(p Point2, s Segment) float64` | PointSegmentDistance returns the distance from p to the closest point on segment s. |
| `func PolygonArea(pts []Point2) float64` | PolygonArea returns the non-negative area enclosed by the polygon pts, independent of vertex winding order. |
| `func PolygonPerimeter(pts []Point2) float64` | PolygonPerimeter returns the total edge length of the closed polygon pts, including the closing edge from the last vertex back to the first. |
| `func PolygonSignedArea(pts []Point2) float64` | PolygonSignedArea returns the signed area of the closed polygon defined by pts using the shoelace (Gauss) formula. |
| `func SegmentSegmentDistance(s1, s2 Segment, eps float64) float64` | SegmentSegmentDistance returns the minimum distance between segments s1 and s2, which is zero when they intersect. |
| `func SegmentsIntersect(s1, s2 Segment, eps float64) bool` | SegmentsIntersect reports whether segments s1 and s2 share at least one point. |
| `func SignedArea2(a, b, c Point2) float64` | SignedArea2 returns twice the signed area of triangle abc, equal to the cross product (b-a)×(c-a). |
| `func SignedPointLineDistance(p Point2, l Line) float64` | SignedPointLineDistance returns the signed perpendicular distance from p to line l. |
| `func Slope(p, q Point2) (float64, bool)` | Slope returns the slope (dy/dx) of the line through p and q and whether it is finite. |
| `func TriangleArea(a, b, c Point2) float64` | TriangleArea returns the (non-negative) area of triangle abc. |
| `func WindingNumber(p Point2, pts []Point2) int` | WindingNumber returns the winding number of the closed polygon pts around point p: the signed number of times the polygon boundary travels… |

### Types

| Type | What it is |
| --- | --- |
| `BoundingBox` | BoundingBox is an axis-aligned rectangle spanning from Min (lower-left) to Max (upper-right). |
| `Circle` | Circle is the set of points at exactly Radius distance from Center. |
| `Line` | Line is the infinite straight line passing through the two distinct points A and B. The direction of the line is B-A. |
| `Point2` | Point2 is a point (location) in the Euclidean plane. |
| `Polygon` | Polygon is a planar polygon described by an ordered list of vertices. |
| `Segment` | Segment is the finite straight segment with endpoints A and B. |
| `Triangle` | Triangle is the triangle with vertices A, B and C. |
| `Vec2` | Vec2 is a vector (displacement or direction) in the Euclidean plane. |

<details>
<summary><code>BoundingBox</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BoundingBoxOf(pts ...Point2) (BoundingBox, bool)` | BoundingBoxOf returns the smallest axis-aligned box containing all the given points. |
| `func NewBoundingBox(p, q Point2) BoundingBox` | NewBoundingBox returns the axis-aligned box spanning the two opposite corners p and q. |
| `func PolygonBoundingBox(pts []Point2) (BoundingBox, bool)` | PolygonBoundingBox returns the axis-aligned bounding box of the polygon pts. |
| `func (b BoundingBox) Area() float64` | Area returns the area of the box. |
| `func (b BoundingBox) Center() Point2` | Center returns the centre point of the box. |
| `func (b BoundingBox) ContainsPoint(p Point2) bool` | ContainsPoint reports whether p lies within the closed box (boundary included). |
| `func (b BoundingBox) Corners() [4]Point2` | Corners returns the four corners of the box in counter-clockwise order starting from Min (lower-left). |
| `func (b BoundingBox) Height() float64` | Height returns the vertical extent of the box. |
| `func (b BoundingBox) Intersection(o BoundingBox) (BoundingBox, bool)` | Intersection returns the overlapping box of b and o. |
| `func (b BoundingBox) Intersects(o BoundingBox) bool` | Intersects reports whether the two boxes overlap (sharing an edge counts as intersecting). |
| `func (b BoundingBox) Perimeter() float64` | Perimeter returns the perimeter of the box. |
| `func (b BoundingBox) Union(o BoundingBox) BoundingBox` | Union returns the smallest box containing both b and o. |
| `func (b BoundingBox) Width() float64` | Width returns the horizontal extent of the box. |

</details>

<details>
<summary><code>Circle</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func CircleFrom3Points(a, b, c Point2) (Circle, bool)` | CircleFrom3Points returns the unique circle passing through three points. |
| `func NewCircle(c Point2, radius float64) Circle` | NewCircle returns the circle centred at c with the given radius. |
| `func TriangleCircumcircle(a, b, c Point2) (Circle, bool)` | TriangleCircumcircle returns the circle through points a, b and c. |
| `func TriangleIncircle(a, b, c Point2) (Circle, bool)` | TriangleIncircle returns the inscribed circle of the triangle with vertices a, b and c. |
| `func (c Circle) Area() float64` | Area returns the area enclosed by the circle. |
| `func (c Circle) BoundingBox() BoundingBox` | BoundingBox returns the axis-aligned bounding box that tightly encloses the circle. |
| `func (c Circle) Circumference() float64` | Circumference returns the perimeter (circumference) of the circle. |
| `func (c Circle) ContainsPoint(p Point2) bool` | ContainsPoint reports whether p lies inside or on the circle. |
| `func (c Circle) Diameter() float64` | Diameter returns the diameter of the circle. |
| `func (c Circle) OnBoundary(p Point2, eps float64) bool` | OnBoundary reports whether p lies within eps of the circle boundary. |
| `func (c Circle) PointAtAngle(theta float64) Point2` | PointAtAngle returns the point on the circle at the given angle in radians, measured counter-clockwise from the positive x-axis. |

</details>

<details>
<summary><code>Line</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func PerpendicularBisector(p, q Point2) Line` | PerpendicularBisector returns the line that passes through the midpoint of p and q at right angles to the segment pq. |
| `func (l Line) ClosestPoint(p Point2) Point2` | ClosestPoint returns the point on the line nearest to p. |
| `func (l Line) Contains(p Point2, eps float64) bool` | Contains reports whether p lies on the infinite line l within tolerance eps (measured as perpendicular distance). |
| `func (l Line) Direction() Vec2` | Direction returns the unit direction vector of the line (from A toward B). |
| `func (l Line) Distance(p Point2) float64` | Distance returns the perpendicular distance from p to the line. |
| `func (l Line) PointAt(t float64) Point2` | PointAt returns the point A + t*(B-A) on the line, so t=0 gives A and t=1 gives B. |
| `func (l Line) Vec() Vec2` | Vec returns the displacement vector from A to B along the line. |

</details>

<details>
<summary><code>Point2</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Centroid(pts ...Point2) Point2` | Centroid returns the arithmetic mean (centroid) of the given points. |
| `func CircleCircleIntersection(c1, c2 Circle) ([]Point2, bool)` | CircleCircleIntersection returns the intersection points of two circles. |
| `func CircleLineIntersection(c Circle, l Line) []Point2` | CircleLineIntersection returns the intersection points of circle c with the infinite line l. |
| `func CircleSegmentIntersection(c Circle, s Segment) []Point2` | CircleSegmentIntersection returns the intersection points of circle c with the finite segment s, restricted to points that actually lie on the… |
| `func ClosestPointOnLine(p Point2, l Line) Point2` | ClosestPointOnLine returns the point on the infinite line l nearest to p (the foot of the perpendicular from p). |
| `func ClosestPointOnSegment(p Point2, s Segment) Point2` | ClosestPointOnSegment returns the point on segment s nearest to p, obtained by projecting p onto the segment's line and clamping to the endpoints. |
| `func Lerp(p, q Point2, t float64) Point2` | Lerp returns the point interpolated between p and q by parameter t: (1-t)*p + t*q. |
| `func LineIntersection(l1, l2 Line, eps float64) (Point2, bool)` | LineIntersection returns the intersection point of the infinite lines l1 and l2 and reports whether it exists. |
| `func LineSegmentIntersection(l Line, s Segment, eps float64) (Point2, bool)` | LineSegmentIntersection returns the intersection point of infinite line l with segment s and reports whether it exists. |
| `func Midpoint(p, q Point2) Point2` | Midpoint returns the point halfway between p and q. |
| `func NewPoint2(x, y float64) Point2` | NewPoint2 returns the point (x, y). |
| `func PointFromPolar(r, theta float64) Point2` | PointFromPolar returns the point at distance r and angle theta (radians, measured counter-clockwise from the positive X axis) from the origin. |
| `func PolygonCentroid(pts []Point2) (Point2, bool)` | PolygonCentroid returns the area centroid (centre of mass of a uniform lamina) of the closed polygon pts. |
| `func ProjectPointOntoLine(p Point2, l Line) Point2` | ProjectPointOntoLine returns the foot of the perpendicular from p to line l. |
| `func ReflectPointAcrossLine(p Point2, l Line) Point2` | ReflectPointAcrossLine returns p reflected across the infinite line l. |
| `func RotateAround(p, c Point2, theta float64) Point2` | RotateAround returns p rotated by theta radians counter-clockwise about the pivot point c. |
| `func SegmentIntersection(s1, s2 Segment, eps float64) (Point2, bool)` | SegmentIntersection returns the unique intersection point of segments s1 and s2 and reports whether such a single point exists. |
| `func (p Point2) Add(v Vec2) Point2` | Add returns the point p displaced by the vector v. |
| `func (p Point2) ApproxEqual(q Point2, eps float64) bool` | ApproxEqual reports whether p and q are within eps (absolute) in each coordinate. |
| `func (p Point2) DistanceSqTo(q Point2) float64` | DistanceSqTo returns the squared Euclidean distance from p to q. |
| `func (p Point2) DistanceTo(q Point2) float64` | DistanceTo returns the Euclidean distance from p to q. |
| `func (p Point2) Equal(q Point2) bool` | Equal reports whether p and q are exactly equal in both coordinates. |
| `func (p Point2) Sub(q Point2) Vec2` | Sub returns the vector from q to p (that is, p-q). |
| `func (p Point2) Vec() Vec2` | Vec returns the position vector of p (the displacement from the origin to p). |

</details>

<details>
<summary><code>Polygon</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ClipPolygon(subject, clip []Point2) Polygon` | ClipPolygon clips the subject polygon against a convex clip polygon using the Sutherland-Hodgman algorithm and returns the intersection as a new… |
| `func ClipPolygonHalfPlane(subject []Point2, a, b Point2) Polygon` | ClipPolygonHalfPlane clips the subject polygon against the half-plane lying to the left of the directed line through a and b (the inside, for a… |
| `func ConvexHull(pts []Point2) Polygon` | ConvexHull returns the convex hull of the input points as a polygon wound counter-clockwise, using Andrew's monotone-chain algorithm in O(n log n)… |
| `func (poly Polygon) Area() float64` | Area returns the non-negative area enclosed by the polygon. |
| `func (poly Polygon) BoundingBox() (BoundingBox, bool)` | BoundingBox returns the axis-aligned bounding box of the polygon. |
| `func (poly Polygon) Centroid() (Point2, bool)` | Centroid returns the area centroid of the polygon. |
| `func (poly Polygon) ContainsPoint(p Point2) bool` | ContainsPoint reports whether p lies inside the polygon using the ray-casting rule. |
| `func (poly Polygon) ContainsPointWinding(p Point2) bool` | ContainsPointWinding reports whether p lies inside the polygon using the non-zero winding rule. |
| `func (poly Polygon) Edges() []Segment` | Edges returns the polygon's boundary edges as segments, including the closing edge from the last vertex back to the first. |
| `func (poly Polygon) EnsureCCW() Polygon` | EnsureCCW returns a polygon with the same shape guaranteed to be wound counter-clockwise, reversing the receiver only when necessary. |
| `func (poly Polygon) EnsureCW() Polygon` | EnsureCW returns a polygon with the same shape guaranteed to be wound clockwise, reversing the receiver only when necessary. |
| `func (poly Polygon) IsClockwise() bool` | IsClockwise reports whether the polygon vertices are ordered clockwise (negative signed area). |
| `func (poly Polygon) IsConvex() bool` | IsConvex reports whether the polygon is convex. |
| `func (poly Polygon) IsCounterClockwise() bool` | IsCounterClockwise reports whether the polygon vertices are ordered counter-clockwise (positive signed area). |
| `func (poly Polygon) NumVertices() int` | NumVertices returns the number of vertices in the polygon. |
| `func (poly Polygon) OnBoundary(p Point2, eps float64) bool` | OnBoundary reports whether p lies within eps of the polygon boundary. |
| `func (poly Polygon) Perimeter() float64` | Perimeter returns the total boundary length of the polygon. |
| `func (poly Polygon) Reverse() Polygon` | Reverse returns a new polygon with the vertex order reversed, flipping the winding direction. |
| `func (poly Polygon) RotateAbout(theta float64, center Point2) Polygon` | RotateAbout returns a copy of the polygon rotated by theta radians (positive counter-clockwise) about the fixed point center. |
| `func (poly Polygon) ScaleAbout(factor float64, center Point2) Polygon` | ScaleAbout returns a copy of the polygon scaled by factor about the fixed point center. |
| `func (poly Polygon) SignedArea() float64` | SignedArea returns the signed shoelace area of the polygon; positive for counter-clockwise winding, negative for clockwise. |
| `func (poly Polygon) Translate(v Vec2) Polygon` | Translate returns a copy of the polygon with every vertex displaced by v. |
| `func (poly Polygon) WindingNumber(p Point2) int` | WindingNumber returns the winding number of the polygon around p. |

</details>

<details>
<summary><code>Segment</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (s Segment) AsLine() Line` | AsLine returns the infinite line supporting the segment. |
| `func (s Segment) ClosestPoint(p Point2) Point2` | ClosestPoint returns the point on the segment nearest to p. |
| `func (s Segment) Contains(p Point2, eps float64) bool` | Contains reports whether p lies on segment s within tolerance eps (measured as distance to the closest point on the segment). |
| `func (s Segment) Direction() Vec2` | Direction returns the unit direction vector of the segment (from A toward B). |
| `func (s Segment) Distance(p Point2) float64` | Distance returns the distance from p to the closest point on the segment. |
| `func (s Segment) Length() float64` | Length returns the Euclidean length of the segment. |
| `func (s Segment) LengthSq() float64` | LengthSq returns the squared Euclidean length of the segment. |
| `func (s Segment) Midpoint() Point2` | Midpoint returns the midpoint of the segment. |
| `func (s Segment) PointAt(t float64) Point2` | PointAt returns the point A + t*(B-A); t is not clamped, so values outside [0, 1] lie on the segment's supporting line beyond its endpoints. |
| `func (s Segment) Vec() Vec2` | Vec returns the displacement vector from A to B along the segment. |

</details>

<details>
<summary><code>Triangle</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewTriangle(a, b, c Point2) Triangle` | NewTriangle returns the triangle with the given vertices. |
| `func (t Triangle) Area() float64` | Area returns the non-negative area of the triangle. |
| `func (t Triangle) Barycentric(p Point2) (u, v, w float64, ok bool)` | Barycentric returns the barycentric coordinates (u, v, w) of point p with respect to the triangle, satisfying p = u*A + v*B + w*C and u+v+w = 1. |
| `func (t Triangle) BoundingBox() BoundingBox` | BoundingBox returns the axis-aligned bounding box of the triangle. |
| `func (t Triangle) Centroid() Point2` | Centroid returns the triangle's centroid (the average of its vertices, which coincides with the centre of mass). |
| `func (t Triangle) Circumcenter() (Point2, bool)` | Circumcenter returns the centre of the triangle's circumscribed circle (the point equidistant from all three vertices). |
| `func (t Triangle) Circumcircle() (Circle, bool)` | Circumcircle returns the circle passing through all three vertices of the triangle. |
| `func (t Triangle) Circumradius() (float64, bool)` | Circumradius returns the radius of the triangle's circumscribed circle. |
| `func (t Triangle) ContainsPoint(p Point2) bool` | ContainsPoint reports whether p lies inside or on the boundary of the triangle, using barycentric coordinates. |
| `func (t Triangle) Incenter() (Point2, bool)` | Incenter returns the centre of the triangle's inscribed circle, the intersection of the angle bisectors. |
| `func (t Triangle) Incircle() (Circle, bool)` | Incircle returns the largest circle contained in the triangle, tangent to all three sides. |
| `func (t Triangle) Inradius() (float64, bool)` | Inradius returns the radius of the triangle's inscribed circle, equal to the area divided by the semiperimeter. |
| `func (t Triangle) Orthocenter() (Point2, bool)` | Orthocenter returns the point where the triangle's three altitudes meet. |
| `func (t Triangle) Perimeter() float64` | Perimeter returns the sum of the triangle's three side lengths. |
| `func (t Triangle) PointFromBarycentric(u, v, w float64) Point2` | PointFromBarycentric returns the Cartesian point with barycentric coordinates (u, v, w) relative to the triangle. |
| `func (t Triangle) SignedArea() float64` | SignedArea returns the signed area of the triangle; positive when A, B, C are in counter-clockwise order, negative when clockwise. |

</details>

<details>
<summary><code>Vec2</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Add(u, v Vec2) Vec2` | Add returns the vector sum u+v. |
| `func NewVec2(x, y float64) Vec2` | NewVec2 returns the vector (x, y). |
| `func Scale(u Vec2, s float64) Vec2` | Scale returns the vector u multiplied by the scalar s. |
| `func Sub(u, v Vec2) Vec2` | Sub returns the vector difference u-v. |
| `func VecFromPolar(r, theta float64) Vec2` | VecFromPolar returns the vector of magnitude r pointing at angle theta (radians, counter-clockwise from the positive X axis). |
| `func (u Vec2) Add(v Vec2) Vec2` | Add returns the vector sum u+v. |
| `func (u Vec2) Angle() float64` | Angle returns the direction of u in radians in (-π, π], measured counter-clockwise from the positive X axis. |
| `func (u Vec2) AngleBetween(v Vec2) float64` | AngleBetween returns the unsigned angle in [0, π] between u and v. |
| `func (u Vec2) ApproxEqual(v Vec2, eps float64) bool` | ApproxEqual reports whether u and v are within eps (absolute) in each component. |
| `func (v Vec2) AsPoint() Point2` | AsPoint returns v interpreted as the point at position v from the origin. |
| `func (u Vec2) Cross(v Vec2) float64` | Cross returns the 2D cross product u×v = u.X*v.Y - u.Y*v.X, the signed area of the parallelogram spanned by u and v. |
| `func (u Vec2) Dot(v Vec2) float64` | Dot returns the dot (scalar) product u·v = u.X*v.X + u.Y*v.Y. |
| `func (u Vec2) Equal(v Vec2) bool` | Equal reports whether u and v are exactly equal in both components. |
| `func (u Vec2) Hadamard(v Vec2) Vec2` | Hadamard returns the component-wise product (u.X*v.X, u.Y*v.Y). |
| `func (u Vec2) IsZero(eps float64) bool` | IsZero reports whether u is the zero vector within eps: its length is at most eps. |
| `func (u Vec2) Lerp(v Vec2, t float64) Vec2` | Lerp returns the linear interpolation (1-t)*u + t*v. |
| `func (u Vec2) Neg() Vec2` | Neg returns the vector -u. |
| `func (u Vec2) Norm() float64` | Norm returns the Euclidean length \|u\| of the vector. |
| `func (u Vec2) NormSq() float64` | NormSq returns the squared Euclidean length \|u\|², avoiding the square root. |
| `func (u Vec2) Normalize() (Vec2, float64)` | Normalize returns the unit vector in the direction of u and its original length. |
| `func (u Vec2) Perp() Vec2` | Perp returns u rotated 90° counter-clockwise: (-u.Y, u.X). |
| `func (u Vec2) PerpCW() Vec2` | PerpCW returns u rotated 90° clockwise: (u.Y, -u.X). |
| `func (u Vec2) Project(v Vec2) Vec2` | Project returns the vector projection of u onto v (the component of u parallel to v). |
| `func (u Vec2) Reflect(axis Vec2) Vec2` | Reflect returns u reflected across the line through the origin with direction axis. |
| `func (u Vec2) ReflectAcrossNormal(n Vec2) Vec2` | ReflectAcrossNormal returns u reflected as if bouncing off a surface whose unit normal is n: u - 2(u·n)n. |
| `func (u Vec2) Reject(v Vec2) Vec2` | Reject returns the component of u orthogonal to v, that is u minus the projection of u onto v. |
| `func (u Vec2) Rotate(theta float64) Vec2` | Rotate returns u rotated by theta radians counter-clockwise about the origin. |
| `func (u Vec2) ScalarProject(v Vec2) float64` | ScalarProject returns the signed length of the projection of u onto v. |
| `func (u Vec2) Scale(s float64) Vec2` | Scale returns the vector u multiplied by the scalar s. |
| `func (u Vec2) SignedAngle(v Vec2) float64` | SignedAngle returns the signed angle in (-π, π] from u to v, positive when v is counter-clockwise from u. |
| `func (u Vec2) Sub(v Vec2) Vec2` | Sub returns the vector difference u-v. |
| `func (u Vec2) Unit() Vec2` | Unit returns the unit vector in the direction of u, or the zero vector if u has zero length. |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/geom2d).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
