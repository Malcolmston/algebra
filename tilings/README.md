# tilings — provides tools for planar symmetry and tilings

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/tilings.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/tilings)

Package tilings provides tools for planar symmetry and tilings: the 17
wallpaper groups and 7 frieze groups (their generators, point groups and
orbits), planar isometries (rotations, reflections, translations and glide
reflections) together with the finite point groups they build, the regular and
semiregular (Archimedean) tilings described by vertex configurations,
substitution tilings (Penrose via Robinson triangles, the chair/L-tromino
rep-tile and the pinwheel triangle), and two-dimensional lattice and orbifold
utilities.

The geometric core is the `Isometry` type, an affine map x -> M x + t whose
linear part M is orthogonal (determinant +/-1). Every planar isometry is one
of five kinds - identity, translation, rotation, reflection or glide
reflection - and `Isometry.Classify` recovers the kind together with the
invariant data (rotation centre and angle, mirror axis, glide vector). Finite
collections of origin-fixing isometries form the cyclic groups `CyclicGroup`
and dihedral groups `DihedralGroup` used as planar point groups.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/tilings"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AngleBetween(u, v Point) float64` | AngleBetween returns the unsigned angle between vectors u and v in [0, pi]. |
| `func Circumradius(n int, s float64) float64` | Circumradius returns the circumradius of a regular n-gon with side length s. |
| `func Deg2Rad(deg float64) float64` | Deg2Rad converts an angle in degrees to radians. |
| `func DualTilingName(name string) (string, bool)` | DualTilingName returns the name of the dual (Laves/Catalan) tiling of a regular tiling: triangular and hexagonal are dual to each other and the… |
| `func ExteriorAngle(n int) float64` | ExteriorAngle returns the exterior angle of a regular n-gon in radians, 2*pi/n. |
| `func FriezeIUCNames() []string` | FriezeIUCNames returns the 7 Hermann-Mauguin symbols in conventional order. |
| `func GyrationOrders(symbol string) []int` | GyrationOrders returns the orders of the cone (gyration) points listed before any '*' in the orbifold symbol; '∞' is reported as 0. |
| `func HasCrossCap(symbol string) bool` | HasCrossCap reports whether the orbifold symbol contains a cross-cap ('x'). |
| `func HasMirror(symbol string) bool` | HasMirror reports whether the orbifold symbol contains a mirror boundary '*'. |
| `func Inradius(n int, s float64) float64` | Inradius returns the inradius (apothem) of a regular n-gon with side length s. |
| `func InteriorAngle(n int) float64` | InteriorAngle returns the interior angle of a regular n-gon in radians, (n-2)*pi/n. |
| `func InteriorAngleDeg(n int) float64` | InteriorAngleDeg returns the interior angle of a regular n-gon in degrees. |
| `func IsArchimedeanConfig(v VertexConfiguration) bool` | IsArchimedeanConfig reports whether the given vertex configuration matches one of the eleven uniform tilings, comparing as a cyclic sequence up to… |
| `func IsCrystallographicRestriction(n int) bool` | IsCrystallographicRestriction reports whether n is an allowed rotation order for a two-dimensional crystallographic (lattice-preserving) symmetry,… |
| `func IsWallpaperSignature(symbol string, eps float64) bool` | IsWallpaperSignature reports whether the orbifold symbol has total cost 2 (to within eps), the signature of a two-dimensional space group. |
| `func KaleidoscopeCorners(symbol string) []int` | KaleidoscopeCorners returns the orders of the corner (kaleidoscope) points listed after a '*' in the orbifold symbol; '∞' is reported as 0. |
| `func NormalizeAngle(theta float64) float64` | NormalizeAngle returns theta reduced to the half-open interval [0, 2*pi). |
| `func NormalizeAngleSigned(theta float64) float64` | NormalizeAngleSigned returns theta reduced to the interval (-pi, pi]. |
| `func OrbifoldCost(symbol string) (float64, error)` | OrbifoldCost returns the total cost of a Conway orbifold symbol, the sum of the individual feature costs used in Conway's magic theorem. |
| `func OrbifoldEulerCharacteristic(symbol string) (float64, error)` | OrbifoldEulerCharacteristic returns the orbifold Euler characteristic 2 - OrbifoldCost(symbol). |
| `func PolygonArea(vertices []Point) float64` | PolygonArea returns the signed area of the polygon with the given vertices in order; a positive value indicates counter-clockwise orientation. |
| `func PolygonAreaFromSide(n int, s float64) float64` | PolygonAreaFromSide returns the area of a regular n-gon with side length s. |
| `func PolygonAreaRegular(n int, r float64) float64` | PolygonArea2 returns the area of a regular n-gon with circumradius r. |
| `func PolygonContains(vertices []Point, p Point) bool` | PolygonContains reports whether p lies inside the simple polygon with the given vertices, using the even-odd ray-casting rule. |
| `func Rad2Deg(rad float64) float64` | Rad2Deg converts an angle in radians to degrees. |
| `func SignedAngle(u, v Point) float64` | SignedAngle returns the signed angle from u to v in (-pi, pi], positive for a counter-clockwise turn. |
| `func WallpaperIUCNames() []string` | WallpaperIUCNames returns the 17 international symbols in conventional order. |

### Types

| Type | What it is |
| --- | --- |
| `Affine2` | Affine2 is a general affine map of the plane, x -> M x + t, with a 2x2 linear part that need not be orthogonal (so it may include scaling and shear). |
| `Chair` | Chair is one tile of the chair (L-tromino) substitution tiling, represented by the affine Frame that maps the canonical L-tromino into the plane. |
| `Classification` | Classification bundles the invariant data recovered from an isometry. |
| `FriezeGroup` | FriezeGroup describes one of the 7 frieze groups, the discrete symmetry groups of an infinitely long strip with translations along a single axis… |
| `Isometry` | Isometry is a planar isometry, the affine map x -> M x + t whose linear part M = [[M00, M01], [M10, M11]] is orthogonal (determinant +/-1). |
| `IsometryKind` | IsometryKind enumerates the five kinds of planar isometry. |
| `Lattice2D` | Lattice2D is a two-dimensional lattice generated by integer combinations of two basis vectors V1 and V2. |
| `LatticeType` | LatticeType enumerates the five Bravais lattice classes in the plane. |
| `PinwheelTriangle` | PinwheelTriangle is a right triangle similar to the (1, 2, sqrt5) prototile of the pinwheel tiling, with the right angle at A, the long leg ending at… |
| `Point` | Point is a point in the plane, also used as a two-dimensional vector. |
| `PointGroup` | PointGroup is a finite group of planar isometries that all fix a common point (taken to be the origin). |
| `Rhombus` | Rhombus is a Penrose rhombus formed from two Robinson triangles of the same kind sharing their long edge. |
| `RobinsonKind` | RobinsonKind labels the two Robinson (half-rhombus) triangle types used by the Penrose deflation. |
| `RobinsonTriangle` | RobinsonTriangle is an oriented Robinson triangle with vertices A, B and C. It is one half of a Penrose rhombus; two triangles of the same kind… |
| `UniformTiling` | UniformTiling names a regular or semiregular (Archimedean) tiling of the plane by its vertex configuration. |
| `VertexConfiguration` | VertexConfiguration is the cyclic sequence of polygon sizes meeting at a vertex of a tiling; for example {4, 8, 8} for the truncated square tiling. |
| `WallpaperGroup` | WallpaperGroup describes one of the 17 two-dimensional crystallographic (wallpaper) groups: its names, translation lattice and a finite set of coset… |

<details>
<summary><code>Affine2</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func FromIsometry(a Isometry) Affine2` | FromIsometry lifts an isometry to an affine map. |
| `func IdentityAffine() Affine2` | IdentityAffine returns the identity affine map. |
| `func ScaleAffine(s float64) Affine2` | ScaleAffine returns the uniform scaling by factor s about the origin. |
| `func TranslateAffine(dx, dy float64) Affine2` | TranslateAffine returns the translation by (dx, dy). |
| `func (a Affine2) Apply(p Point) Point` | Apply returns the image of the point p under the affine map. |
| `func (a Affine2) Compose(inner Affine2) Affine2` | Compose returns the affine map a∘inner (inner applied first). |
| `func (a Affine2) Determinant() float64` | Determinant returns the determinant of the linear part; its absolute value is the area scaling factor of the map. |

</details>

<details>
<summary><code>Chair</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ChairSubstitute(chairs []Chair) []Chair` | ChairSubstitute applies one deflation step to a set of chair tiles. |
| `func ChairTiling(n int) []Chair` | ChairTiling returns the chair tiles obtained by deflating the canonical chair n times. |
| `func NewChair() Chair` | NewChair returns the canonical chair tile (identity frame). |
| `func (c Chair) Area() float64` | Area returns the area of the chair tile (three times the square of its linear scale). |
| `func (c Chair) Polygon() []Point` | Polygon returns the six boundary vertices of the chair tile in order. |
| `func (c Chair) Substitute() []Chair` | Substitute returns the four half-scale chair tiles that subdivide c. |

</details>

<details>
<summary><code>FriezeGroup</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func FriezeGroupByIUC(iuc string) (FriezeGroup, bool)` | FriezeGroupByIUC returns the frieze group with the given Hermann-Mauguin symbol at unit period and reports whether it was found. |
| `func FriezeGroups() []FriezeGroup` | FriezeGroups returns all 7 frieze groups with unit period in the conventional order (p1, p11g, p1m1, p2, p2mg, p11m, p2mm). |
| `func (g FriezeGroup) Cosets() []Isometry` | Cosets returns the coset representatives of the frieze group modulo axial translation. |
| `func (g FriezeGroup) GroupElements(n int) []Isometry` | GroupElements returns the isometries of the group whose translation lies within n periods of the origin. |
| `func (g FriezeGroup) HasGlideReflection() bool` | HasGlideReflection reports whether the group contains a glide reflection. |
| `func (g FriezeGroup) HasReflection() bool` | HasReflection reports whether the group contains a reflection. |
| `func (g FriezeGroup) HasRotation() bool` | HasRotation reports whether the group contains a 180-degree rotation. |
| `func (g FriezeGroup) Orbit(p Point, n int, eps float64) []Point` | Orbit returns the images of the point p under GroupElements(n), deduplicated to within eps and sorted for determinism. |
| `func (g FriezeGroup) PointGroupOrder() int` | PointGroupOrder returns the order of the associated point group, equal to the number of coset representatives. |

</details>

<details>
<summary><code>Isometry</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ComposeAll(seq ...Isometry) Isometry` | ComposeAll returns the composition of the given isometries applied left to right: the first element is applied first. |
| `func GlideReflection(p Point, axis, glide float64) Isometry` | GlideReflection returns the glide reflection whose mirror line passes through point p in the direction axis (radians) and whose glide vector has… |
| `func Identity() Isometry` | Identity returns the identity isometry. |
| `func Reflection(axis float64) Isometry` | Reflection returns the reflection across the line through the origin whose direction makes angle axis (radians) with the positive x-axis. |
| `func ReflectionLine(p Point, axis float64) Isometry` | ReflectionLine returns the reflection across the line through point p whose direction makes angle axis (radians) with the positive x-axis. |
| `func Rotation(theta float64) Isometry` | Rotation returns the rotation about the origin by theta radians. |
| `func RotationAbout(center Point, theta float64) Isometry` | RotationAbout returns the rotation about the point center by theta radians. |
| `func Translation(dx, dy float64) Isometry` | Translation returns the translation by the vector (dx, dy). |
| `func TranslationVec(v Point) Isometry` | TranslationVec returns the translation by the vector v. |
| `func (a Isometry) Apply(p Point) Point` | Apply returns the image of the point p under the isometry. |
| `func (a Isometry) ApplyVec(v Point) Point` | ApplyVec returns the image of the vector v under the linear part only, ignoring the translation. |
| `func (a Isometry) ApproxEqual(b Isometry, eps float64) bool` | ApproxEqual reports whether a and b agree to within eps in every entry. |
| `func (a Isometry) Classify(eps float64) Classification` | Classify decomposes the isometry into its kind and invariant data using the tolerance eps. |
| `func (a Isometry) Compose(inner Isometry) Isometry` | Compose returns the isometry a∘inner, i.e. |
| `func (a Isometry) Determinant() float64` | Determinant returns the determinant of the linear part, either +1 or -1 for a true isometry. |
| `func (a Isometry) FixedPoint() (Point, bool)` | FixedPoint returns a point fixed by the isometry and reports whether one exists. |
| `func (a Isometry) Inverse() Isometry` | Inverse returns the inverse isometry. |
| `func (a Isometry) IsDirect() bool` | IsDirect reports whether the isometry preserves orientation (determinant +1). |
| `func (a Isometry) IsIdentity(eps float64) bool` | IsIdentity reports whether the isometry is the identity to within eps. |
| `func (a Isometry) Kind(eps float64) IsometryKind` | Kind returns just the kind of the isometry using tolerance eps. |
| `func (a Isometry) LinearAngle() float64` | LinearAngle returns atan2(M10, M00), the rotation angle of the linear part for a direct isometry and, for an opposite isometry, twice the mirror-axis… |
| `func (a Isometry) Order(max int, eps float64) int` | Order returns the smallest positive integer k such that a^k is the identity to within eps, searching up to max. |
| `func (a Isometry) RotationAngle() float64` | RotationAngle returns the rotation angle in (-pi, pi] for a direct isometry. |
| `func (a Isometry) String() string` | String returns a short human-readable description of the isometry. |
| `func (a Isometry) Then(b Isometry) Isometry` | Then returns the isometry b∘a, i.e. |
| `func (a Isometry) Translation() Point` | Translation returns the translation part (Tx, Ty) of the isometry as a point. |

</details>

<details>
<summary><code>IsometryKind</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (k IsometryKind) String() string` | String returns the name of the isometry kind. |

</details>

<details>
<summary><code>Lattice2D</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func HexagonalLattice(a float64) Lattice2D` | HexagonalLattice returns the hexagonal (triangular) lattice of spacing a, with basis vectors at 0 and 60 degrees. |
| `func NewLattice2D(v1, v2 Point) Lattice2D` | NewLattice2D returns the lattice with the given basis vectors. |
| `func ObliqueLattice(a, b, gamma float64) Lattice2D` | ObliqueLattice returns a generic oblique lattice with the given basis lengths and included angle gamma (radians). |
| `func RectangularLattice(a, b float64) Lattice2D` | RectangularLattice returns the rectangular lattice with spacings a and b. |
| `func RhombicLattice(a, h float64) Lattice2D` | RhombicLattice returns the centred-rectangular (rhombic) lattice with primitive basis (a/2, h) and (-a/2, h). |
| `func SquareLattice(a float64) Lattice2D` | SquareLattice returns the square lattice of spacing a. |
| `func (l Lattice2D) Classify(eps float64) LatticeType` | Classify returns the Bravais lattice type of the lattice using tolerance eps, judged from a reduced basis. |
| `func (l Lattice2D) Dual() Lattice2D` | Dual returns the dual (reciprocal) lattice, whose basis vectors W satisfy Wi·Vj = delta_ij. |
| `func (l Lattice2D) FromFractional(a, b float64) Point` | FromFractional returns the point a*V1 + b*V2. |
| `func (l Lattice2D) FundamentalArea() float64` | FundamentalArea returns the area of the lattice's fundamental parallelogram, the absolute value of the cross product of the basis vectors. |
| `func (l Lattice2D) GaussReduced() Lattice2D` | GaussReduced returns a Gauss/Lagrange-reduced basis for the same lattice: two shortest vectors realising the successive minima, with V1 no longer… |
| `func (l Lattice2D) Gram() (g11, g12, g22 float64)` | Gram returns the Gram matrix entries (g11, g12, g22) of the basis, namely V1·V1, V1·V2 and V2·V2. |
| `func (l Lattice2D) LatticePoints(n int) []Point` | LatticePoints returns all lattice points i*V1 + j*V2 with i, j in [-n, n], in row-major order. |
| `func (l Lattice2D) NearestPoint(p Point) Point` | NearestPoint returns the lattice point closest to p (by rounding basis coordinates); it is exact for well-conditioned bases. |
| `func (l Lattice2D) PackingRadius() float64` | PackingRadius returns the packing radius of the lattice, half the length of a shortest nonzero vector. |
| `func (l Lattice2D) Point(i, j int) Point` | Point returns the lattice point i*V1 + j*V2. |
| `func (l Lattice2D) Reduce(p Point) Point` | Reduce returns the representative of p modulo the lattice lying in the fundamental parallelogram [0,1) x [0,1) in basis coordinates. |
| `func (l Lattice2D) ShortestLength() float64` | ShortestLength returns the length of a shortest nonzero lattice vector, the square root of the minimum squared distance. |
| `func (l Lattice2D) ShortestVector() Point` | ShortestVector returns a nonzero lattice vector of minimum length. |
| `func (l Lattice2D) ToFractional(p Point) (float64, float64)` | ToFractional returns the coordinates (a, b) of p in the basis, so that p = a*V1 + b*V2. |
| `func (l Lattice2D) VoronoiRelevantVectors() []Point` | VoronoiRelevantVectors returns the lattice vectors whose bisecting lines bound the Voronoi cell of the origin (the Wigner-Seitz cell). |

</details>

<details>
<summary><code>LatticeType</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (t LatticeType) PointGroupOrder() int` | PointGroupOrder returns the order of the holohedry (full point-symmetry group) of the lattice type: 2, 4, 4, 8 or 12 for oblique, rectangular,… |
| `func (t LatticeType) String() string` | String returns the name of the lattice type. |

</details>

<details>
<summary><code>PinwheelTriangle</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func CanonicalPinwheel() PinwheelTriangle` | CanonicalPinwheel returns the prototile with the right angle at the origin, the long leg of length 2 along the x-axis and the short leg of length 1… |
| `func PinwheelSubstitute(tris []PinwheelTriangle) []PinwheelTriangle` | PinwheelSubstitute applies one deflation step to a set of pinwheel triangles. |
| `func PinwheelTiling(n int) []PinwheelTriangle` | PinwheelTiling returns the pinwheel triangles obtained by deflating the canonical prototile n times. |
| `func (t PinwheelTriangle) Area() float64` | Area returns the unsigned area of the triangle. |
| `func (t PinwheelTriangle) Centroid() Point` | Centroid returns the centroid of the triangle. |
| `func (t PinwheelTriangle) Substitute() []PinwheelTriangle` | Substitute returns the five sub-triangles, each similar to the prototile at scale 1/sqrt5, that exactly tile t. |
| `func (t PinwheelTriangle) Vertices() []Point` | Vertices returns the three vertices A, B, C. |

</details>

<details>
<summary><code>Point</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Centroid(pts ...Point) Point` | Centroid returns the arithmetic mean of the given points; the empty input yields the origin. |
| `func NewPoint(x, y float64) Point` | NewPoint returns the point (x, y). |
| `func Origin() Point` | Origin returns the point (0, 0). |
| `func Polar(r, theta float64) Point` | Polar returns the point at radius r and angle theta (radians) from the origin. |
| `func PolygonCentroid(vertices []Point) Point` | PolygonCentroid returns the area centroid of the simple polygon with the given vertices in order; for a degenerate (zero-area) polygon it returns the… |
| `func RegularPolygon(n int, c Point, r, theta0 float64) []Point` | RegularPolygon returns the vertices of a regular n-gon centred at c with circumradius r; the first vertex is placed at angle theta0 and vertices run… |
| `func (p Point) Add(q Point) Point` | Add returns the vector sum p+q. |
| `func (p Point) Angle() float64` | Angle returns the direction of p measured counter-clockwise from the positive x-axis, in radians in (-pi, pi]. |
| `func (p Point) ApproxEqual(q Point, eps float64) bool` | ApproxEqual reports whether p and q are within eps of one another. |
| `func (p Point) Cross(q Point) float64` | Cross returns the scalar cross product p.X*q.Y - p.Y*q.X. |
| `func (p Point) Distance(q Point) float64` | Distance returns the Euclidean distance between p and q. |
| `func (p Point) Distance2(q Point) float64` | Distance2 returns the squared Euclidean distance between p and q. |
| `func (p Point) Dot(q Point) float64` | Dot returns the Euclidean dot product p·q. |
| `func (p Point) IsZero(eps float64) bool` | IsZero reports whether p is within eps of the origin. |
| `func (p Point) Lerp(q Point, t float64) Point` | Lerp returns the linear interpolation (1-t)*p + t*q. |
| `func (p Point) Midpoint(q Point) Point` | Midpoint returns the midpoint of p and q. |
| `func (p Point) Neg() Point` | Neg returns the additive inverse -p. |
| `func (p Point) Norm() float64` | Norm returns the Euclidean length of p. |
| `func (p Point) Norm2() float64` | Norm2 returns the squared Euclidean length of p. |
| `func (p Point) Normalize() Point` | Normalize returns the unit vector in the direction of p; the zero vector is returned unchanged. |
| `func (p Point) Perp() Point` | Perp returns p rotated by +90 degrees, i.e. |
| `func (p Point) Rotate(theta float64) Point` | Rotate returns p rotated about the origin by theta radians. |
| `func (p Point) RotateAbout(c Point, theta float64) Point` | RotateAbout returns p rotated about the point c by theta radians. |
| `func (p Point) Scale(s float64) Point` | Scale returns p scaled by the factor s. |
| `func (p Point) Sub(q Point) Point` | Sub returns the vector difference p-q. |

</details>

<details>
<summary><code>PointGroup</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func CrystallographicPointGroups() []PointGroup` | CrystallographicPointGroups returns the ten two-dimensional crystallographic point groups (C1, C2, C3, C4, C6 and D1, D2, D3, D4, D6) in that order. |
| `func CyclicGroup(n int) PointGroup` | CyclicGroup returns the cyclic point group C_n of rotations about the origin by multiples of 2*pi/n. |
| `func DihedralGroup(n int) PointGroup` | DihedralGroup returns the dihedral point group D_n of order 2n: the n rotations of C_n together with n reflections whose axes pass through the… |
| `func (g PointGroup) Contains(a Isometry, eps float64) bool` | Contains reports whether the point group contains an isometry equal to a within eps. |
| `func (g PointGroup) HasReflection() bool` | HasReflection reports whether the point group contains any reflection. |
| `func (g PointGroup) IsClosed(eps float64) bool` | IsClosed reports whether the point group is closed under composition to within eps: the product of any two elements is again an element. |
| `func (g PointGroup) MaxRotationOrder() int` | MaxRotationOrder returns the highest rotation order in the point group, or 1 if the group contains no nontrivial rotation. |
| `func (g PointGroup) NumReflections() int` | NumReflections returns the number of reflections (opposite isometries) in the point group. |
| `func (g PointGroup) Orbit(p Point, eps float64) []Point` | Orbit returns the distinct images of the point p under every element of the point group, deduplicated to within eps. |
| `func (g PointGroup) Order() int` | Order returns the number of elements in the point group. |
| `func (g PointGroup) RotationOrders() []int` | RotationOrders returns the sorted set of rotation orders present in the point group, where a rotation by 2*pi/k contributes k. |
| `func (g PointGroup) StabilizerOrder(p Point, eps float64) int` | StabilizerOrder returns the number of group elements fixing the point p to within eps. |

</details>

<details>
<summary><code>Rhombus</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func PenroseKitesDarts(tris []RobinsonTriangle) []Rhombus` | PenroseKitesDarts is a synonym assembly returning the Robinson triangles grouped into rhombi; the P2 kite/dart tiles are locally derivable from these. |
| `func PenroseRhombi(tris []RobinsonTriangle) []Rhombus` | PenroseRhombi assembles Robinson triangles into Penrose rhombi (the P3 tiles) by pairing triangles of the same kind that share their longest edge. |

</details>

<details>
<summary><code>RobinsonKind</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (k RobinsonKind) String() string` | String returns the name of the Robinson triangle kind. |

</details>

<details>
<summary><code>RobinsonTriangle</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DeflateRobinson(tris []RobinsonTriangle) []RobinsonTriangle` | DeflateRobinson applies one step of the Penrose deflation (with inflation factor Phi) to each triangle, returning the finer set of triangles. |
| `func DeflateRobinsonN(tris []RobinsonTriangle, n int) []RobinsonTriangle` | DeflateRobinsonN applies n rounds of Penrose deflation. |
| `func PenroseP2(n int, r float64) []RobinsonTriangle` | PenroseP2 returns the Robinson triangles of a Penrose P2 (kite and dart) tiling obtained by deflating the sun seed of radius r a total of n times. |
| `func PenroseP3(n int, r float64) []RobinsonTriangle` | PenroseP3 returns the Robinson triangles of a Penrose P3 (rhombus) tiling obtained by deflating the sun seed of radius r a total of n times. |
| `func PenroseSun(r float64) []RobinsonTriangle` | PenroseSun returns the "sun" seed patch: ten Robinson triangles arranged with five-fold symmetry around the origin, each with circumradius r. |
| `func (t RobinsonTriangle) Area() float64` | Area returns the unsigned area of the triangle. |
| `func (t RobinsonTriangle) Centroid() Point` | Centroid returns the centroid of the triangle. |
| `func (t RobinsonTriangle) Vertices() []Point` | Vertices returns the three vertices A, B, C of the triangle. |

</details>

<details>
<summary><code>UniformTiling</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ArchimedeanTilings() []UniformTiling` | ArchimedeanTilings returns all eleven uniform tilings of the plane (the three regular tilings followed by the eight semiregular ones). |
| `func RegularTilings() []UniformTiling` | RegularTilings returns the three regular tilings of the plane (triangular, square, hexagonal). |
| `func SemiregularTilings() []UniformTiling` | SemiregularTilings returns the eight semiregular (Archimedean) tilings of the plane. |
| `func TilingByName(name string) (UniformTiling, bool)` | TilingByName returns the uniform tiling with the given name and reports whether it was found. |
| `func (t UniformTiling) String() string` | String returns the tiling's name and configuration. |

</details>

<details>
<summary><code>VertexConfiguration</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func EnumerateVertexTypes() []VertexConfiguration` | EnumerateVertexTypes returns every multiset of regular-polygon sizes (each at least 3) whose interior angles sum to exactly 2*pi, i.e. |
| `func (v VertexConfiguration) AngleSum() float64` | AngleSum returns the sum of the interior angles of the polygons in the configuration, in radians. |
| `func (v VertexConfiguration) AngleSumDeg() float64` | AngleSumDeg returns the interior-angle sum in degrees. |
| `func (v VertexConfiguration) Degree() int` | Degree returns the number of polygons meeting at the vertex. |
| `func (v VertexConfiguration) EqualMultiset(w VertexConfiguration) bool` | EqualMultiset reports whether two configurations use the same polygons the same number of times, ignoring order. |
| `func (v VertexConfiguration) IsPlanar(eps float64) bool` | IsPlanar reports whether the configuration's interior angles sum to exactly 2*pi (360 degrees) to within eps, the condition for a flat vertex. |
| `func (v VertexConfiguration) IsRegular(eps float64) bool` | IsRegular reports whether the configuration consists of three or more copies of a single polygon type meeting at a flat vertex (a regular tiling… |
| `func (v VertexConfiguration) Multiset() []int` | Multiset returns the polygon sizes in sorted order, discarding the cyclic arrangement. |
| `func (v VertexConfiguration) String() string` | String returns the configuration in dotted notation, e.g. |

</details>

<details>
<summary><code>WallpaperGroup</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func WallpaperGroupByIUC(iuc string) (WallpaperGroup, bool)` | WallpaperGroupByIUC returns the wallpaper group with the given international symbol and reports whether it was found. |
| `func WallpaperGroupByOrbifold(orb string) (WallpaperGroup, bool)` | WallpaperGroupByOrbifold returns the wallpaper group with the given orbifold symbol and reports whether it was found. |
| `func WallpaperGroups() []WallpaperGroup` | WallpaperGroups returns all 17 wallpaper groups in the conventional order (p1, p2, pm, pg, cm, pmm, pmg, pgg, cmm, p4, p4m, p4g, p3, p3m1, p31m, p6,… |
| `func (g WallpaperGroup) Cosets() []Isometry` | Cosets returns the coset representatives of the group modulo its translation lattice. |
| `func (g WallpaperGroup) Generators() []Isometry` | Generators returns the generating isometries used to build the group (not including the lattice translations). |
| `func (g WallpaperGroup) GroupElements(n int) []Isometry` | GroupElements returns the isometries of the group whose translation parts lie within n lattice steps of the origin: every coset representative… |
| `func (g WallpaperGroup) HasGlideReflection() bool` | HasGlideReflection reports whether the group contains a glide reflection. |
| `func (g WallpaperGroup) HasReflection() bool` | HasReflection reports whether the group contains a reflection. |
| `func (g WallpaperGroup) MaxRotationOrder() int` | MaxRotationOrder returns the highest order of rotation in the group. |
| `func (g WallpaperGroup) Orbit(p Point, n int, eps float64) []Point` | Orbit returns the images of the point p under GroupElements(n), deduplicated to within eps and sorted for determinism. |
| `func (g WallpaperGroup) PointGroupOrder() int` | PointGroupOrder returns the order of the point group, equal to the number of coset representatives. |
| `func (g WallpaperGroup) TranslationLattice() Lattice2D` | TranslationLattice returns the concrete translation lattice of the group. |

</details>

### Constants

`Epsilon`, `Phi`

### Variables

`ErrInvalidOrbifold`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/tilings).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
