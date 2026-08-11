# geom3d — provides primitives for 3D computational geometry

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/geom3d.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/geom3d)

Package geom3d provides primitives for 3D computational geometry.

The core value type is `Vec3`, a three-component vector used both as a
location (point) and as a displacement (direction). On top of it the package
offers the standard toolbox of spatial geometry: vector arithmetic (add,
subtract, scale, dot, cross, triple product), norms, distances, angles,
normalization, projection/rejection and reflection; 3x3 matrix algebra
(`Mat3`) with determinant, inverse, transpose and multiplication; rotation
matrices about the coordinate axes and about an arbitrary axis (Rodrigues'
rotation formula); conversion between Euler angles and rotation matrices;
planes (`Plane`) and lines (`Line3`) with closest-point and distance queries;
ray casting (`Ray`) against planes, triangles (Moeller-Trumbore), axis-aligned
boxes (`AABB`) and spheres (`Sphere`); barycentric coordinates; and
intersection tests between boxes and spheres.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/geom3d"
```

## Usage

This is the package's own `ExampleAABBFromPoints`, so it compiles and its output is
asserted on every `go test ./geom3d/`.

```go
bx, ok := geom3d.AABBFromPoints(
		geom3d.NewVec3(0, 0, 0),
		geom3d.NewVec3(3, -1, 2),
		geom3d.NewVec3(-4, 5, 1),
	)
	fmt.Println(ok, bx.Min, bx.Max)

	_, ok = geom3d.AABBFromPoints()
	fmt.Println(ok)
```

```
true {-4 -1 0} {3 5 2}
false
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func Barycentric(p, a, b, c Vec3) (u, v, w float64, ok bool)` | Barycentric returns the barycentric coordinates (u, v, w) of point p with respect to triangle (a, b, c), so that p ≈ u*a + v*b + w*c and u+v+w = 1. |
| `func PointInTriangle(p, a, b, c Vec3) bool` | PointInTriangle reports whether point p lies inside (or on the edge of) triangle (a, b, c). |
| `func RayTriangle(r Ray, a, b, c Vec3) (t, u, v float64, ok bool)` | RayTriangle intersects a ray with triangle (a, b, c) using the Moeller-Trumbore algorithm. |
| `func ScalarTriple(a, b, c Vec3) float64` | ScalarTriple returns the scalar triple product a·(b×c), which equals the signed volume of the parallelepiped spanned by a, b and c. |

### Types

| Type | What it is |
| --- | --- |
| `AABB` | AABB is an axis-aligned bounding box defined by its minimum and maximum corners. |
| `EulerAngles` | EulerAngles holds a set of Tait-Bryan angles in radians using the aerospace convention: Roll about the x-axis, Pitch about the y-axis and Yaw about… |
| `Line3` | Line3 is the infinite straight line passing through Point in the direction Dir. |
| `Mat3` | Mat3 is a 3x3 matrix of float64 stored in row-major order: element m`r``c` lives in row r and column c. |
| `Plane` | Plane is an oriented plane in space, represented by the equation Normal·x = D. Points on the positive side of the plane (in the direction of Normal)… |
| `Ray` | Ray is a half-line starting at Origin and extending in the direction Dir. |
| `Sphere` | Sphere is a solid ball defined by its Center and Radius. |
| `Vec3` | Vec3 is a vector in three-dimensional Euclidean space. |

<details>
<summary><code>AABB</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AABBFromPoints(pts ...Vec3) (AABB, bool)` | AABBFromPoints returns the smallest axis-aligned box containing all the given points, and true. |
| `func NewAABB(a, b Vec3) AABB` | NewAABB returns the axis-aligned box with the given corners, normalizing them so that Min holds the component-wise minimum and Max the component-wise… |
| `func (bx AABB) Center() Vec3` | Center returns the geometric center of the box. |
| `func (bx AABB) ClosestPoint(p Vec3) Vec3` | ClosestPoint returns the point within the box (on its surface or interior) nearest to p, obtained by clamping p to the box on each axis. |
| `func (bx AABB) Contains(p Vec3) bool` | Contains reports whether point p lies inside or on the boundary of the box. |
| `func (bx AABB) Expand(v Vec3) AABB` | Expand returns the box grown outward by the vector v on every side; each component of v is subtracted from Min and added to Max. |
| `func (bx AABB) Extents() Vec3` | Extents returns the half-size of the box along each axis. |
| `func (bx AABB) IntersectRay(r Ray) (tmin, tmax float64, ok bool)` | IntersectRay intersects a ray with the box using the slab method. |
| `func (bx AABB) IntersectsAABB(other AABB) bool` | IntersectsAABB reports whether the box overlaps the other box, touching boundaries included. |
| `func (bx AABB) IntersectsSphere(s Sphere) bool` | IntersectsSphere reports whether the box overlaps the sphere s. |
| `func (bx AABB) Size() Vec3` | Size returns the full extent of the box along each axis (Max - Min). |
| `func (bx AABB) Union(other AABB) AABB` | Union returns the smallest axis-aligned box containing both bx and other. |

</details>

<details>
<summary><code>EulerAngles</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func MatrixToEuler(m Mat3) EulerAngles` | MatrixToEuler extracts the Tait-Bryan angles (intrinsic Z-Y-X convention) that produce the rotation matrix m, inverting `EulerToMatrix`. |

</details>

<details>
<summary><code>Line3</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func LineThrough(a, b Vec3) Line3` | LineThrough returns the line passing through the two distinct points a and b, with direction b-a. |
| `func (l Line3) At(t float64) Vec3` | At returns the point Point + t*Dir on the line. |
| `func (l Line3) ClosestPoint(p Vec3) Vec3` | ClosestPoint returns the point on the line nearest to p. |
| `func (l Line3) Distance(p Vec3) float64` | Distance returns the perpendicular distance from point p to the line. |
| `func (l Line3) Param(p Vec3) float64` | Param returns the parameter t of the point on the line closest to p, so that `Line3.At`(t) is the closest point. |

</details>

<details>
<summary><code>Mat3</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AxisAngleMatrix(axis Vec3, angle float64) Mat3` | AxisAngleMatrix returns the rotation matrix for a rotation by angle radians about the given axis, using Rodrigues' rotation formula. |
| `func EulerToMatrix(e EulerAngles) Mat3` | EulerToMatrix returns the rotation matrix for the given Tait-Bryan angles using the intrinsic Z-Y-X convention, i.e. |
| `func Mat3FromColumns(c0, c1, c2 Vec3) Mat3` | Mat3FromColumns builds a matrix whose columns are the given vectors c0, c1 and c2. |
| `func Mat3FromRows(r0, r1, r2 Vec3) Mat3` | Mat3FromRows builds a matrix whose rows are the given vectors r0, r1 and r2. |
| `func Mat3Identity() Mat3` | Mat3Identity returns the 3x3 identity matrix. |
| `func RotationX(angle float64) Mat3` | RotationX returns the matrix that rotates a vector by angle radians about the x-axis, following the right-hand rule. |
| `func RotationY(angle float64) Mat3` | RotationY returns the matrix that rotates a vector by angle radians about the y-axis, following the right-hand rule. |
| `func RotationZ(angle float64) Mat3` | RotationZ returns the matrix that rotates a vector by angle radians about the z-axis, following the right-hand rule. |
| `func (m Mat3) Add(n Mat3) Mat3` | Add returns the component-wise sum m+n. |
| `func (m Mat3) Col(j int) Vec3` | Col returns column j (0, 1 or 2) of m as a vector. |
| `func (m Mat3) Determinant() float64` | Determinant returns the determinant of m. |
| `func (m Mat3) Equal(n Mat3, eps float64) bool` | Equal reports whether m and n are equal within the absolute tolerance eps, compared element-wise. |
| `func (m Mat3) Inverse() (Mat3, bool)` | Inverse returns the matrix inverse of m and true, or the zero matrix and false if m is singular (its determinant is zero within the package… |
| `func (m Mat3) Mul(n Mat3) Mat3` | Mul returns the matrix product m*n. |
| `func (m Mat3) MulVec(v Vec3) Vec3` | MulVec returns the matrix-vector product m*v, treating v as a column vector. |
| `func (m Mat3) Row(i int) Vec3` | Row returns row i (0, 1 or 2) of m as a vector. |
| `func (m Mat3) Scale(s float64) Mat3` | Scale returns m with every element multiplied by the scalar s. |
| `func (m Mat3) String() string` | String returns a human-readable, row-by-row representation of m. |
| `func (m Mat3) Sub(n Mat3) Mat3` | Sub returns the component-wise difference m-n. |
| `func (m Mat3) Trace() float64` | Trace returns the sum of the diagonal elements of m. |
| `func (m Mat3) Transpose() Mat3` | Transpose returns the transpose of m. |

</details>

<details>
<summary><code>Plane</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func PlaneFromPointNormal(p, normal Vec3) Plane` | PlaneFromPointNormal returns the plane passing through point p with the given normal. |
| `func PlaneFromPoints(a, b, c Vec3) (Plane, bool)` | PlaneFromPoints returns the plane through the three points a, b and c and true, or the zero plane and false if the points are collinear. |
| `func (pl Plane) ClosestPoint(p Vec3) Vec3` | ClosestPoint returns the orthogonal projection of point p onto the plane, the point on the plane nearest to p. |
| `func (pl Plane) Contains(p Vec3, eps float64) bool` | Contains reports whether point p lies on the plane within the absolute tolerance eps. |
| `func (pl Plane) Distance(p Vec3) float64` | Distance returns the unsigned (absolute) distance from point p to the plane. |
| `func (pl Plane) Normalize() Plane` | Normalize returns the plane rescaled so that Normal has unit length while describing the same set of points. |
| `func (pl Plane) Side(p Vec3, eps float64) int` | Side reports which side of the plane point p lies on: +1 for the positive (normal) side, -1 for the negative side, and 0 if p lies on the plane… |
| `func (pl Plane) SignedDistance(p Vec3) float64` | SignedDistance returns the signed distance from point p to the plane. |

</details>

<details>
<summary><code>Ray</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewRay(origin, dir Vec3) Ray` | NewRay returns the ray from origin with the given direction. |
| `func (r Ray) At(t float64) Vec3` | At returns the point Origin + t*Dir on the ray. |
| `func (r Ray) IntersectPlane(pl Plane) (t float64, ok bool)` | IntersectPlane returns the parameter t at which the ray meets the plane and true, or 0 and false if the ray is parallel to the plane (or points away… |

</details>

<details>
<summary><code>Sphere</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (s Sphere) Contains(p Vec3) bool` | Contains reports whether point p lies inside or on the surface of the sphere. |
| `func (s Sphere) IntersectRay(r Ray) (tmin, tmax float64, ok bool)` | IntersectRay intersects a ray with the sphere. |
| `func (s Sphere) IntersectsAABB(bx AABB) bool` | IntersectsAABB reports whether the sphere overlaps the axis-aligned box bx, by comparing the box's closest point to the sphere's radius. |
| `func (s Sphere) IntersectsSphere(other Sphere) bool` | IntersectsSphere reports whether the sphere overlaps the other sphere, touching surfaces included. |
| `func (s Sphere) IsEmpty() bool` | IsEmpty reports whether the sphere encloses no points at all, which is the case exactly when Radius is negative. |

</details>

<details>
<summary><code>Vec3</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ClosestPointsBetweenLines(l1, l2 Line3) (p1, p2 Vec3, ok bool)` | ClosestPointsBetweenLines returns the pair of points, one on line l1 and one on line l2, that are closest to each other, together with true. |
| `func NewVec3(x, y, z float64) Vec3` | NewVec3 returns the vector with the given components. |
| `func Rodrigues(v, axis Vec3, angle float64) Vec3` | Rodrigues rotates the vector v by angle radians about the given axis using Rodrigues' rotation formula, without forming a matrix. |
| `func (v Vec3) Abs() Vec3` | Abs returns the component-wise absolute value of v. |
| `func (v Vec3) Add(w Vec3) Vec3` | Add returns the component-wise sum v+w. |
| `func (v Vec3) Angle(w Vec3) float64` | Angle returns the unsigned angle in radians between v and w, in the range [0, pi]. |
| `func (v Vec3) Cross(w Vec3) Vec3` | Cross returns the cross (vector) product v×w. |
| `func (v Vec3) Distance(w Vec3) float64` | Distance returns the Euclidean distance between the points v and w. |
| `func (v Vec3) DistanceSq(w Vec3) float64` | DistanceSq returns the squared Euclidean distance between the points v and w. |
| `func (v Vec3) Div(s float64) Vec3` | Div returns v with every component divided by the scalar s. |
| `func (v Vec3) Dot(w Vec3) float64` | Dot returns the dot (inner) product v·w. |
| `func (v Vec3) Equal(w Vec3, eps float64) bool` | Equal reports whether v and w are equal within the absolute tolerance eps, compared component-wise. |
| `func (v Vec3) Hadamard(w Vec3) Vec3` | Hadamard returns the component-wise (Hadamard) product of v and w. |
| `func (v Vec3) IsZero() bool` | IsZero reports whether every component of v is zero within the package default tolerance. |
| `func (v Vec3) Length() float64` | Length returns the Euclidean length (magnitude) of v. |
| `func (v Vec3) LengthSq() float64` | LengthSq returns the squared Euclidean length v·v. |
| `func (v Vec3) Lerp(w Vec3, t float64) Vec3` | Lerp returns the linear interpolation (1-t)*v + t*w. |
| `func (v Vec3) Max(w Vec3) Vec3` | Max returns the component-wise maximum of v and w. |
| `func (v Vec3) Min(w Vec3) Vec3` | Min returns the component-wise minimum of v and w. |
| `func (v Vec3) Neg() Vec3` | Neg returns the additive inverse -v. |
| `func (v Vec3) Normalize() (Vec3, float64)` | Normalize returns the unit vector pointing in the same direction as v, together with the original length. |
| `func (v Vec3) Project(w Vec3) Vec3` | Project returns the vector projection of v onto w, that is the component of v parallel to w. |
| `func (v Vec3) Reflect(n Vec3) Vec3` | Reflect returns v reflected across the plane whose unit normal is n. |
| `func (v Vec3) Reject(w Vec3) Vec3` | Reject returns the vector rejection of v from w, that is the component of v orthogonal to w: v minus its projection onto w. |
| `func (v Vec3) Scale(s float64) Vec3` | Scale returns v with every component multiplied by the scalar s. |
| `func (v Vec3) Sub(w Vec3) Vec3` | Sub returns the component-wise difference v-w. |
| `func (v Vec3) Unit() Vec3` | Unit returns the unit vector in the direction of v, or the zero vector if v has zero length. |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/geom3d).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
