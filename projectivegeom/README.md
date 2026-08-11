# projectivegeom — implements projective and inversive geometry in the real projective plane RP^2 and real

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/projectivegeom.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/projectivegeom)

Package projectivegeom implements projective and inversive geometry in the
real projective plane RP^2 and real projective space RP^3.

The plane RP^2 is modelled with homogeneous coordinates: a `Point` is a
non-zero triple [X Y Z] identified up to non-zero scale, and a `Line` is a
non-zero triple [A B C] denoting the locus A*X+B*Y+C*Z=0. Points and lines are
dual: joining two points (`Join`) and meeting two lines (`Meet`) are the same
cross-product operation, and every construction has a dual obtained by
swapping the two roles.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/projectivegeom"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AreHarmonic(a, b, c, d Point, tol float64) bool` | AreHarmonic reports whether the four collinear points form a harmonic range, i.e. |
| `func Collinear(p, q, r Point, tol float64) bool` | Collinear reports whether the three points are collinear within tol, i.e. |
| `func Concurrent(l, m, n Line, tol float64) bool` | Concurrent reports whether the three lines pass through a common point within tol (the dual of Collinear). |
| `func Concurrent4Planes(a, b, c, d SPlane, tol float64) bool` | Concurrent4Planes reports whether four planes of RP^3 pass through a common point within tol (the dual of Coplanar). |
| `func ConjugateLines(q Conic, l1, l2 Line, tol float64) bool` | ConjugateLines reports whether two lines are conjugate with respect to the non-degenerate conic, i.e. |
| `func ConjugatePoints(q Conic, p1, p2 Point, tol float64) bool` | ConjugatePoints reports whether p1 and p2 are conjugate with respect to the conic, i.e. |
| `func Coplanar(p, q, r, s SPoint, tol float64) bool` | Coplanar reports whether four points of RP^3 are coplanar within tol, i.e. |
| `func CrossRatio(a, b, c, d Point, tol float64) (float64, error)` | CrossRatio returns the projective cross ratio (a, b; c, d) of four collinear points. |
| `func CrossRatio1D(a, b, c, d float64) float64` | CrossRatio1D returns the cross ratio (a, b; c, d) of four scalars on the affine line, defined as ((c-a)(d-b)) / ((c-b)(d-a)). |
| `func CrossRatioLines(a, b, c, d Line, tol float64) (float64, error)` | CrossRatioLines returns the cross ratio of four concurrent lines, equal to the cross ratio of their four intersection points with any transversal. |
| `func CrossRatioPermuteSwapFirstPair(x float64) float64` | CrossRatioPermuteSwapFirstPair returns the cross ratio value that results from swapping the first pair, given the original value x = (a,b;c,d). |
| `func CrossRatioPermuteSwapInner(x float64) float64` | CrossRatioPermuteSwapInner returns (a,c;b,d) given x = (a,b;c,d), equal to 1 - x, one of the six values of the anharmonic group. |
| `func CrossRatioPermuteSwapPairs(x float64) float64` | CrossRatioPermuteSwapPairs returns the cross ratio (c,d;a,b) given x = (a,b;c,d); the two are equal, so this returns x unchanged. |
| `func CrossRatioSixValues(x float64) [6]float64` | CrossRatioSixValues returns the six values of the cross ratio taken by the symmetric group acting on the four points, namely {x, 1/x, 1-x, 1/(1-x),… |
| `func Incidence(p Point, l Line) float64` | Incidence returns the (scale-dependent) value p.V·l.V, which is zero exactly when p lies on l. |
| `func IncidencePointPlane(p SPoint, s SPlane) float64` | IncidencePointPlane returns the value p.V·s.V, which is zero exactly when the point lies on the plane. |
| `func LinesMeet(l, m LinePluecker, tol float64) bool` | LinesMeet reports whether two lines of RP^3 are coplanar (intersect) within tol. |
| `func OnLine(p Point, l Line, tol float64) bool` | OnLine reports whether p lies on l within tol, using coordinates scaled by their largest magnitudes so the test is independent of homogeneous scale. |
| `func OnPlane(p SPoint, s SPlane, tol float64) bool` | OnPlane reports whether the point lies on the plane within tol using scale-invariant normalization. |
| `func ParallelLines(l, m Line, tol float64) bool` | ParallelLines reports whether l and m are parallel in the affine sense, i.e. |
| `func PerpendicularLines(l, m Line, tol float64) bool` | PerpendicularLines reports whether the finite lines l and m are perpendicular within tol (using their affine normals (a, b)). |
| `func ReciprocalProduct(l, m LinePluecker) float64` | ReciprocalProduct returns the bilinear side operator of two lines, P01*Q23 + P23*Q01 - P02*Q13 - P13*Q02 + P03*Q12 + P12*Q03. |

### Types

| Type | What it is |
| --- | --- |
| `Camera` | Camera is a finite pinhole camera given by a 3x4 projection matrix P = K[R\|t], mapping scene points of RP^3 to image points of RP^2. |
| `Collineation3` | Collineation3 is a projective transformation of RP^3, an invertible 4x4 matrix acting on homogeneous point coordinates. |
| `Conic` | Conic is a plane conic represented by a symmetric 3x3 matrix M, denoting the locus of homogeneous points x with x^T M x = 0. |
| `ConicType` | ConicType classifies the affine shape of a conic. |
| `Homography` | Homography is a projective transformation (collineation) of RP^2, represented by an invertible 3x3 matrix acting on homogeneous point coordinates. |
| `Line` | Line is a line of RP^2 in homogeneous (dual) coordinates [A B C] denoting the locus A*X + B*Y + C*Z = 0, defined up to a non-zero scalar. |
| `LinePluecker` | LinePluecker holds the six Pluecker coordinates of a line in RP^3, in the order (P01, P02, P03, P12, P13, P23) where Pij = p_i*q_j - p_j*q_i for two… |
| `Mat3` | Mat3 is a 3x3 real matrix in row-major order. |
| `Mat3x4` | Mat3x4 is a 3x4 real matrix in row-major order, the shape of a pinhole camera projection matrix mapping homogeneous points of RP^3 to homogeneous… |
| `Mat4` | Mat4 is a 4x4 real matrix in row-major order. |
| `PappusConfig` | PappusConfig holds two triples of collinear points, one triple on each of two lines, as required by the theorem of Pappus. |
| `Point` | Point is a point of the real projective plane RP^2 in homogeneous coordinates [X Y Z], defined up to a non-zero scalar. |
| `SPlane` | SPlane is a plane of RP^3 in homogeneous (dual) coordinates [A B C D] denoting the locus A*X + B*Y + C*Z + D*W = 0. |
| `SPoint` | SPoint is a point of real projective space RP^3 in homogeneous coordinates [X Y Z W], defined up to a non-zero scalar. |
| `Triangle` | Triangle holds three points assumed to be in general position (no two equal, not collinear). |
| `TrianglePair` | TrianglePair holds two triangles for the theorem of Desargues. |
| `Vec3` | Vec3 is a 3-component real vector. |
| `Vec4` | Vec4 is a 4-component real vector, the carrier for homogeneous points and planes of RP^3. |

<details>
<summary><code>Camera</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func CameraFromKRt(k, r Mat3, t Vec3) Camera` | CameraFromKRt builds the camera P = K[R \| t] from calibration matrix K, rotation R and translation t. |
| `func CanonicalCamera() Camera` | CanonicalCamera returns the camera [I \| 0] whose center is the origin and whose axes coincide with the world axes. |
| `func LookAt(k Mat3, eye, target, up Vec3) (Camera, error)` | LookAt builds a camera positioned at eye, viewing target, with the given world up direction and calibration K. It constructs a right-handed camera… |
| `func NewCamera(p Mat3x4) Camera` | NewCamera wraps a 3x4 matrix as a camera. |
| `func (c Camera) Center() SPoint` | Center returns the camera center, the scene point C with P*C = 0, computed as the null vector of P via its four 3x3 minors. |
| `func (c Camera) Depth(x SPoint) float64` | Depth returns the signed depth (cheirality value) of a scene point: it is positive when the finite point lies in front of the camera for a camera… |
| `func (c Camera) PrincipalPlane() SPlane` | PrincipalPlane returns the plane of scene points that project onto the image line at infinity, namely the third row of P. |
| `func (c Camera) Project(x SPoint) Point` | Project returns the image of a scene point under the camera, as a homogeneous point of RP^2. |
| `func (c Camera) ProjectAffine(x SPoint) (u, v float64, err error)` | ProjectAffine returns the pixel coordinates of a scene point, or ErrAtInfinity when the projected point lies on the image line at infinity (the point… |

</details>

<details>
<summary><code>Collineation3</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Compose3(h, g Collineation3) Collineation3` | Compose3 returns the collineation g followed by h. |
| `func IdentityCollineation3() Collineation3` | IdentityCollineation3 returns the identity transformation of RP^3. |
| `func Translation3(dx, dy, dz float64) Collineation3` | Translation3 returns the affine translation of RP^3 by (dx, dy, dz). |
| `func (c Collineation3) Apply(p SPoint) SPoint` | Apply returns the image of the point p under the collineation. |
| `func (c Collineation3) ApplyPlane(s SPlane) (SPlane, error)` | ApplyPlane returns the image of the plane s, transformed by the inverse transpose so that incidence is preserved. |
| `func (c Collineation3) Inverse() (Collineation3, error)` | Inverse returns the inverse collineation, or ErrSingular when the matrix is not invertible. |

</details>

<details>
<summary><code>Conic</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func CircleConic(cx, cy, r float64) Conic` | CircleConic returns the conic for the circle of radius r centered at (cx, cy). |
| `func ConicFromFivePoints(pts [5]Point) (Conic, error)` | ConicFromFivePoints returns the unique conic through the five given points, computed from the 5x5 minors of their monomial matrix. |
| `func ConicFromMatrix(m Mat3) Conic` | ConicFromMatrix wraps a matrix as a conic, replacing it with its symmetric part so that only the quadratic form matters. |
| `func NewConic(a, b, c, d, e, f float64) Conic` | NewConic builds the conic with the given point-equation coefficients a, b, c, d, e, f. |
| `func UnitCircle() Conic` | UnitCircle returns the conic X^2 + Y^2 - Z^2 = 0, the unit circle. |
| `func (q Conic) Center() (Point, error)` | Center returns the center of the conic, the pole of the line at infinity. |
| `func (q Conic) Classify(tol float64) ConicType` | Classify returns the affine type of the conic. |
| `func (q Conic) Coeffs() (a, b, c, d, e, f float64)` | Coeffs returns the point-equation coefficients (a, b, c, d, e, f). |
| `func (q Conic) Det() float64` | Det returns the determinant of the conic matrix. |
| `func (q Conic) Discriminant() float64` | Discriminant returns the determinant of the upper-left 2x2 block, a*c - (b/2)^2, whose sign classifies a non-degenerate conic as an ellipse… |
| `func (q Conic) Dual() (Conic, error)` | Dual returns the dual (line) conic, whose points are the tangent lines of the original. |
| `func (q Conic) Evaluate(p Point) float64` | Evaluate returns the value x^T M x of the conic's quadratic form at p, which is zero exactly when p lies on the conic. |
| `func (q Conic) IntersectLine(l Line) []Point` | IntersectLine returns the intersection points of the line l with the conic: zero points when the line misses the conic, one (returned twice) when it… |
| `func (q Conic) IsDegenerate(tol float64) bool` | IsDegenerate reports whether the conic matrix is singular within tol, measured relative to the cube of its largest entry. |
| `func (q Conic) IsEllipse(tol float64) bool` | IsEllipse reports whether the conic classifies as an ellipse within tol. |
| `func (q Conic) IsHyperbola(tol float64) bool` | IsHyperbola reports whether the conic classifies as a hyperbola within tol. |
| `func (q Conic) IsParabola(tol float64) bool` | IsParabola reports whether the conic classifies as a parabola within tol. |
| `func (q Conic) OnConic(p Point, tol float64) bool` | OnConic reports whether p lies on the conic within tol, using scale-invariant normalization. |
| `func (q Conic) Polar(p Point) Line` | Polar returns the polar line of the point p with respect to the conic, namely the line M*p. |
| `func (q Conic) Pole(l Line) (Point, error)` | Pole returns the pole of the line l with respect to the conic, the point whose polar is l. |
| `func (q Conic) Rank(tol float64) int` | Rank returns the numerical rank (0..3) of the conic matrix, distinguishing a proper conic (rank 3) from a line pair (rank 2) or double line (rank 1). |
| `func (q Conic) TangentAt(p Point) (Line, error)` | TangentAt returns the tangent line to the conic at a point p that lies on it, or ErrNotOnConic when p is not on the conic within a mild tolerance. |
| `func (q Conic) TangentsFrom(p Point) ([2]Line, error)` | TangentsFrom returns the two tangent lines to the conic through an external point p, found as the joins of p with the two contact points where the… |
| `func (q Conic) Transform(h Homography) (Conic, error)` | Transform returns the image of the conic under the homography h. |

</details>

<details>
<summary><code>ConicType</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (t ConicType) String() string` | String returns the conic type name. |

</details>

<details>
<summary><code>Homography</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Affine(a, b, c, d, tx, ty float64) Homography` | Affine returns the affine map [[a b tx],[c d ty],[0 0 1]]. |
| `func Compose(h, g Homography) Homography` | Compose returns the homography g followed by h, whose matrix is h.M * g.M so that Compose(h, g).Apply(p) equals h.Apply(g.Apply(p)). |
| `func HomographyFromPoints(src, dst [4]Point) (Homography, error)` | HomographyFromPoints returns the unique homography mapping the four source points to the four destination points, given no three points in either set… |
| `func IdentityHomography() Homography` | IdentityHomography returns the identity transformation. |
| `func NewHomography(m Mat3) Homography` | NewHomography wraps a 3x3 matrix as a homography without checking invertibility. |
| `func RandomHomography(r *rand.Rand) (Homography, bool)` | RandomHomography returns a pseudo-random homography whose matrix entries are drawn uniformly from [-1, 1). |
| `func ReflectionX() Homography` | ReflectionX returns the Euclidean reflection across the x-axis. |
| `func ReflectionY() Homography` | ReflectionY returns the Euclidean reflection across the y-axis. |
| `func Rotation(theta float64) Homography` | Rotation returns the Euclidean rotation by theta radians about the origin. |
| `func RotationAbout(theta, cx, cy float64) Homography` | RotationAbout returns the Euclidean rotation by theta radians about the finite point (cx, cy). |
| `func Scaling(s float64) Homography` | Scaling returns the uniform affine scaling by s about the origin. |
| `func ScalingXY(sx, sy float64) Homography` | ScalingXY returns the affine scaling by sx along x and sy along y about the origin. |
| `func ShearX(k float64) Homography` | ShearX returns the affine shear x' = x + k*y. |
| `func ShearY(k float64) Homography` | ShearY returns the affine shear y' = y + k*x. |
| `func Similarity(s, theta, dx, dy float64) Homography` | Similarity returns the direct similarity that scales by s, rotates by theta and then translates by (dx, dy). |
| `func Translation(dx, dy float64) Homography` | Translation returns the affine translation by (dx, dy). |
| `func (h Homography) Apply(p Point) Point` | Apply returns the image of the point p under the homography. |
| `func (h Homography) ApplyLine(l Line) (Line, error)` | ApplyLine returns the image of the line l under the homography. |
| `func (h Homography) ApproxEqual(g Homography, tol float64) bool` | ApproxEqual reports whether h and g denote the same homography, comparing their matrices after scaling each by its largest-magnitude entry. |
| `func (h Homography) Determinant() float64` | Determinant returns the determinant of the homography matrix. |
| `func (h Homography) FixedPointsReal(tol float64) []Point` | FixedPointsReal returns the real fixed points of the homography, that is the points p with h.Apply(p) proportional to p. |
| `func (h Homography) Inverse() (Homography, error)` | Inverse returns the inverse homography, or ErrSingular when the matrix is not invertible. |
| `func (h Homography) IsAffine(tol float64) bool` | IsAffine reports whether the homography fixes the line at infinity, i.e. |
| `func (h Homography) IsInvertible() bool` | IsInvertible reports whether the homography matrix is numerically invertible. |
| `func (h Homography) IsSimilarity(tol float64) bool` | IsSimilarity reports whether the homography is affine and its linear 2x2 part is a scaled orthogonal matrix within tol (a Euclidean similarity). |
| `func (h Homography) Then(g Homography) Homography` | Then returns the homography that applies h and then g, i.e. |

</details>

<details>
<summary><code>Line</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DualOfPoint(p Point) Line` | DualOfPoint returns the line whose homogeneous coordinates equal those of the point p. |
| `func Join(p, q Point) (Line, error)` | Join returns the line through the two points p and q. |
| `func LineAtInfinity() Line` | LineAtInfinity returns the line at infinity [0 0 1]. |
| `func LineThrough(p, q Point) (Line, error)` | LineThrough is an alias for Join reading naturally as "the line through p and q". |
| `func NewLine(a, b, c float64) Line` | NewLine returns the homogeneous line [a b c] for a*x + b*y + c = 0. |
| `func PolarPoint(q Conic, p Point) Line` | PolarPoint is a synonym for Conic.Polar reading the point as the pole; it returns the polar line of p with respect to q. |
| `func RandomLine(r *rand.Rand) Line` | RandomLine returns a pseudo-random projective line with coordinates drawn uniformly from [-1, 1). |
| `func (l Line) Equal(m Line, tol float64) bool` | Equal reports whether l and m denote the same projective line. |
| `func (l Line) IsValid() bool` | IsValid reports whether the line is a non-zero homogeneous triple. |
| `func (l Line) LineAngle() float64` | LineAngle returns the direction angle of the finite line [a b c] in radians, in the range (-pi/2, pi/2], measured as atan2(-a, b) folded into that… |
| `func (l Line) Normalized() (Line, bool)` | Normalized returns l rescaled to unit Euclidean length. |

</details>

<details>
<summary><code>LinePluecker</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func LineFromPlanes(s, t SPlane) (LinePluecker, error)` | LineFromPlanes returns the Pluecker coordinates of the line where two distinct planes meet, obtained from their dual coordinates by the Hodge duality. |
| `func LineFromPoints(p, q SPoint) (LinePluecker, error)` | LineFromPoints returns the Pluecker coordinates of the line through two distinct points of RP^3. |
| `func (l LinePluecker) ContainsPoint(x SPoint, tol float64) bool` | ContainsPoint reports whether the point x lies on the line within tol, using the wedge condition x ∧ L = 0. |
| `func (l LinePluecker) Direction() Vec3` | Direction returns a direction vector of the line, taken from its (X,W),(Y,W), (Z,W) coordinates. |
| `func (l LinePluecker) Equal(m LinePluecker, tol float64) bool` | Equal reports whether two lines have proportional Pluecker coordinates within tol. |
| `func (l LinePluecker) GrassmannRelation() float64` | GrassmannRelation returns P01*P23 - P02*P13 + P03*P12, which vanishes exactly for genuine lines. |
| `func (l LinePluecker) IsLine(tol float64) bool` | IsLine reports whether the coordinates satisfy the Grassmann-Pluecker relation within tol and are non-zero, hence represent an actual line. |
| `func (l LinePluecker) MaxAbs() float64` | MaxAbs returns the largest absolute coordinate of the line. |
| `func (l LinePluecker) Moment() Vec3` | Moment returns the moment vector of the line about the origin, from its spatial coordinates. |

</details>

<details>
<summary><code>Mat3</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Diag3(a, b, c float64) Mat3` | Diag3 returns the diagonal matrix with entries a, b, c. |
| `func Identity3() Mat3` | Identity3 returns the 3x3 identity matrix. |
| `func Intrinsics(fx, fy, s, cx, cy float64) Mat3` | Intrinsics returns the upper-triangular camera calibration matrix K with focal lengths fx, fy, principal point (cx, cy) and axis skew s. |
| `func NewMat3(a, b, c, d, e, f, g, h, i float64) Mat3` | NewMat3 builds a matrix from its nine entries given in row-major order. |
| `func Zero3() Mat3` | Zero3 returns the 3x3 zero matrix. |
| `func (m Mat3) Add(n Mat3) Mat3` | Add returns the entry-wise sum m+n. |
| `func (m Mat3) Adjugate() Mat3` | Adjugate returns the classical adjugate (transpose of the cofactor matrix), so that m.Mul(m.Adjugate()) equals Det*Identity. |
| `func (m Mat3) ApproxEqual(n Mat3, tol float64) bool` | ApproxEqual reports whether m and n agree entry-wise within tol. |
| `func (m Mat3) At(i, j int) float64` | At returns the entry in row i, column j. |
| `func (m Mat3) Col(j int) Vec3` | Col returns column j as a Vec3. |
| `func (m Mat3) Det() float64` | Det returns the determinant of m. |
| `func (m Mat3) Inverse() (Mat3, bool)` | Inverse returns the inverse of m and true, or the zero matrix and false when m is (numerically) singular. |
| `func (m Mat3) IsSymmetric(tol float64) bool` | IsSymmetric reports whether m equals its transpose within tol. |
| `func (m Mat3) MaxAbs() float64` | MaxAbs returns the largest absolute entry of m. |
| `func (m Mat3) Mul(n Mat3) Mat3` | Mul returns the matrix product m*n. |
| `func (m Mat3) MulVec(v Vec3) Vec3` | MulVec returns the matrix-vector product m*v. |
| `func (m Mat3) Quad(u, v Vec3) float64` | Quad returns the quadratic form u^T m v. |
| `func (m Mat3) Row(i int) Vec3` | Row returns row i as a Vec3. |
| `func (m Mat3) Scale(s float64) Mat3` | Scale returns m multiplied by the scalar s. |
| `func (m Mat3) String() string` | String renders m on three bracketed rows. |
| `func (m Mat3) Sub(n Mat3) Mat3` | Sub returns the entry-wise difference m-n. |
| `func (m Mat3) Symmetrize() Mat3` | Symmetrize returns (m + m^T)/2, the symmetric part of m. |
| `func (m Mat3) Trace() float64` | Trace returns the sum of the diagonal entries. |
| `func (m Mat3) Transpose() Mat3` | Transpose returns the transpose of m. |

</details>

<details>
<summary><code>Mat3x4</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (m Mat3x4) MulVec4(v Vec4) Vec3` | MulVec4 returns the 3-vector product m*v. |

</details>

<details>
<summary><code>Mat4</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Identity4() Mat4` | Identity4 returns the 4x4 identity matrix. |
| `func (m Mat4) ApproxEqual(n Mat4, tol float64) bool` | ApproxEqual reports whether m and n agree entry-wise within tol. |
| `func (m Mat4) At(i, j int) float64` | At returns the entry in row i, column j. |
| `func (m Mat4) Det() float64` | Det returns the determinant of m. |
| `func (m Mat4) Inverse() (Mat4, bool)` | Inverse returns the inverse of m and true, or the zero matrix and false when m is (numerically) singular. |
| `func (m Mat4) Mul(n Mat4) Mat4` | Mul returns the matrix product m*n. |
| `func (m Mat4) MulVec(v Vec4) Vec4` | MulVec returns the matrix-vector product m*v. |
| `func (m Mat4) Quad(u, v Vec4) float64` | Quad returns the quadratic form u^T m v. |
| `func (m Mat4) Trace() float64` | Trace returns the sum of the diagonal entries. |
| `func (m Mat4) Transpose() Mat4` | Transpose returns the transpose of m. |

</details>

<details>
<summary><code>PappusConfig</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (p PappusConfig) CrossPoints() (x, y, z Point, err error)` | CrossPoints returns the three Pappus cross points |
| `func (p PappusConfig) Holds(tol float64) bool` | Holds reports whether the theorem of Pappus holds numerically for this configuration, i.e. |
| `func (p PappusConfig) PappusLine() (Line, error)` | PappusLine returns the line containing the three Pappus cross points, the Pappus line of the configuration. |

</details>

<details>
<summary><code>Point</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DualOfLine(l Line) Point` | DualOfLine returns the point whose homogeneous coordinates equal those of the line l, the dual of DualOfPoint. |
| `func HarmonicConjugate(a, b, c Point, tol float64) (Point, error)` | HarmonicConjugate returns the fourth point d that forms a harmonic range with a, b, c, i.e. |
| `func Intersect(l, m Line) (Point, error)` | Intersect is an alias for Meet reading naturally as "the intersection of l and m". |
| `func Meet(l, m Line) (Point, error)` | Meet returns the intersection point of the two lines l and m. |
| `func NewPoint(x, y, z float64) Point` | NewPoint returns the homogeneous point [x y z]. |
| `func PerturbPoint(r *rand.Rand, p Point, eps float64) Point` | PerturbPoint returns p with independent uniform noise of amplitude eps added to each homogeneous coordinate, useful for robustness testing. |
| `func PointAtInfinity(dx, dy float64) Point` | PointAtInfinity returns the point at infinity in the direction (dx, dy), namely [dx dy 0]. |
| `func PointFromAffine(x, y float64) Point` | PointFromAffine returns the finite point [x y 1]. |
| `func RandomFinitePoint(r *rand.Rand) Point` | RandomFinitePoint returns a pseudo-random finite point [x y 1] with x, y drawn uniformly from [-1, 1). |
| `func RandomPoint(r *rand.Rand) Point` | RandomPoint returns a pseudo-random projective point whose homogeneous coordinates are drawn uniformly from [-1, 1). |
| `func (p Point) Affine() (x, y float64, err error)` | Affine returns the Euclidean coordinates (x, y) = (X/Z, Y/Z), or an error when p lies on the line at infinity. |
| `func (p Point) Equal(q Point, tol float64) bool` | Equal reports whether p and q denote the same projective point, i.e. |
| `func (p Point) IsAtInfinity(tol float64) bool` | IsAtInfinity reports whether p lies on the line at infinity (Z within tol of zero after scaling by the largest component). |
| `func (p Point) IsValid() bool` | IsValid reports whether the point is a non-zero homogeneous triple. |
| `func (p Point) Normalized() (Point, bool)` | Normalized returns p rescaled so that Z = 1 when finite, or scaled to unit Euclidean length when at infinity. |

</details>

<details>
<summary><code>SPlane</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DualPlaneOfPoint(p SPoint) SPlane` | DualPlaneOfPoint returns the plane whose coordinates equal those of the point p, the dual of DualPointOfPlane. |
| `func NewSPlane(a, b, c, d float64) SPlane` | NewSPlane returns the homogeneous plane [a b c d] for a*x + b*y + c*z + d = 0. |
| `func PlaneAtInfinity() SPlane` | PlaneAtInfinity returns the plane at infinity [0 0 0 1] of RP^3. |
| `func PlaneThrough3Points(p, q, r SPoint) (SPlane, error)` | PlaneThrough3Points returns the plane through three points of RP^3. |
| `func (s SPlane) Equal(t SPlane, tol float64) bool` | Equal reports whether the two planes denote the same projective plane. |
| `func (s SPlane) IsValid() bool` | IsValid reports whether the plane is a non-zero homogeneous 4-tuple. |

</details>

<details>
<summary><code>SPoint</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DualPointOfPlane(s SPlane) SPoint` | DualPointOfPlane returns the point of the dual space RP^3 whose coordinates equal those of the plane s. |
| `func NewSPoint(x, y, z, w float64) SPoint` | NewSPoint returns the homogeneous point [x y z w]. |
| `func PointOf3Planes(a, b, c SPlane) (SPoint, error)` | PointOf3Planes returns the common point of three planes of RP^3. |
| `func RandomSPoint(r *rand.Rand) SPoint` | RandomSPoint returns a pseudo-random point of RP^3 with homogeneous coordinates drawn uniformly from [-1, 1). |
| `func SPointFromAffine(x, y, z float64) SPoint` | SPointFromAffine returns the finite point [x y z 1]. |
| `func (p SPoint) Affine() (x, y, z float64, err error)` | Affine returns the Euclidean coordinates (x, y, z) = (X/W, Y/W, Z/W), or ErrAtInfinity when the point lies on the plane at infinity. |
| `func (p SPoint) Equal(q SPoint, tol float64) bool` | Equal reports whether the two points denote the same projective point. |
| `func (p SPoint) IsAtInfinity(tol float64) bool` | IsAtInfinity reports whether the point lies on the plane at infinity within tol. |
| `func (p SPoint) IsValid() bool` | IsValid reports whether the point is a non-zero homogeneous 4-tuple. |

</details>

<details>
<summary><code>Triangle</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (t Triangle) Centroid() (Point, error)` | Centroid returns the affine centroid of the triangle's finite vertices. |
| `func (t Triangle) IsDegenerate(tol float64) bool` | IsDegenerate reports whether the three vertices are collinear (or coincident) within tol. |
| `func (t Triangle) Sides() (bc, ca, ab Line, err error)` | Sides returns the three side lines BC, CA, AB of the triangle. |

</details>

<details>
<summary><code>TrianglePair</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (tp TrianglePair) AxisPoints() (p, q, r Point, err error)` | AxisPoints returns the three intersections of corresponding sides, |
| `func (tp TrianglePair) DesarguesAxis() (Line, error)` | DesarguesAxis returns the axis of perspectivity, the line through the three AxisPoints. |
| `func (tp TrianglePair) IsPerspectiveFromLine(tol float64) bool` | IsPerspectiveFromLine reports whether the three AxisPoints are collinear within tol. |
| `func (tp TrianglePair) IsPerspectiveFromPoint(tol float64) bool` | IsPerspectiveFromPoint reports whether the three lines joining corresponding vertices are concurrent within tol. |
| `func (tp TrianglePair) PerspectiveCenter() (Point, error)` | PerspectiveCenter returns the center of perspectivity of the two triangles, the common point of the lines A A', B B', C C'. |

</details>

<details>
<summary><code>Vec3</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewVec3(x, y, z float64) Vec3` | NewVec3 returns the vector [x y z]. |
| `func (v Vec3) Add(w Vec3) Vec3` | Add returns the component-wise sum v+w. |
| `func (v Vec3) AngleBetween(w Vec3) float64` | AngleBetween returns the unsigned angle in radians between v and w treated as direction vectors in R^3, in the range [0, pi]. |
| `func (v Vec3) ApproxEqual(w Vec3, tol float64) bool` | ApproxEqual reports whether v and w agree component-wise within tol. |
| `func (v Vec3) Cross(w Vec3) Vec3` | Cross returns the vector cross product v×w. |
| `func (v Vec3) Dot(w Vec3) float64` | Dot returns the Euclidean dot product v·w. |
| `func (v Vec3) Get(i int) float64` | Get returns component i (0=X, 1=Y, 2=Z). |
| `func (v Vec3) IsZero(tol float64) bool` | IsZero reports whether every component of v is within tol of zero. |
| `func (v Vec3) Max() float64` | Max returns the largest component of v. |
| `func (v Vec3) MaxAbs() float64` | MaxAbs returns the largest absolute value among the components of v. |
| `func (v Vec3) Neg() Vec3` | Neg returns the negation -v. |
| `func (v Vec3) Norm() float64` | Norm returns the Euclidean norm \|v\|. |
| `func (v Vec3) Norm2() float64` | Norm2 returns the squared Euclidean norm v·v. |
| `func (v Vec3) Normalized() (Vec3, bool)` | Normalized returns v scaled to unit Euclidean length and true, or the zero vector and false when v is (numerically) zero. |
| `func (v Vec3) Parallel(w Vec3, tol float64) bool` | Parallel reports whether v and w are scalar multiples of one another within tol, i.e. |
| `func (v Vec3) Scale(s float64) Vec3` | Scale returns v scaled by the scalar s. |
| `func (v Vec3) String() string` | String renders v in bracketed homogeneous form. |
| `func (v Vec3) Sub(w Vec3) Vec3` | Sub returns the component-wise difference v-w. |

</details>

<details>
<summary><code>Vec4</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Cross4(a, b, c Vec4) Vec4` | Cross4 returns the generalized cross product of three 4-vectors: the unique (up to scale) vector orthogonal to a, b and c. |
| `func NewVec4(x, y, z, w float64) Vec4` | NewVec4 returns the vector [x y z w]. |
| `func (v Vec4) Add(w Vec4) Vec4` | Add returns the component-wise sum v+w. |
| `func (v Vec4) ApproxEqual(w Vec4, tol float64) bool` | ApproxEqual reports whether v and w agree component-wise within tol. |
| `func (v Vec4) Dot(w Vec4) float64` | Dot returns the Euclidean dot product v·w. |
| `func (v Vec4) Get(i int) float64` | Get returns component i (0=X, 1=Y, 2=Z, 3=W). |
| `func (v Vec4) IsZero(tol float64) bool` | IsZero reports whether every component of v is within tol of zero. |
| `func (v Vec4) MaxAbs() float64` | MaxAbs returns the largest absolute value among the components of v. |
| `func (v Vec4) Neg() Vec4` | Neg returns the negation -v. |
| `func (v Vec4) Norm() float64` | Norm returns the Euclidean norm \|v\|. |
| `func (v Vec4) Norm2() float64` | Norm2 returns the squared Euclidean norm v·v. |
| `func (v Vec4) Normalized() (Vec4, bool)` | Normalized returns v scaled to unit Euclidean length and true, or the zero vector and false when v is (numerically) zero. |
| `func (v Vec4) Scale(s float64) Vec4` | Scale returns v scaled by the scalar s. |
| `func (v Vec4) String() string` | String renders v in bracketed homogeneous form. |
| `func (v Vec4) Sub(w Vec4) Vec4` | Sub returns the component-wise difference v-w. |

</details>

### Constants

`Eps`

### Variables

`ErrAtInfinity`, `ErrConfiguration`, `ErrDegenerate`, `ErrNotCollinear`, `ErrNotOnConic`, `ErrSingular`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/projectivegeom).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
