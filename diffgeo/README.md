# diffgeo — implements classical differential geometry of curves and surfaces in three-dimensional

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/diffgeo.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/diffgeo)

Package diffgeo implements classical differential geometry of curves and
surfaces in three-dimensional Euclidean space using only the Go standard
library.

- A `Curve` is a function t ↦ (x, y, z) tracing a parametric space curve.
From it the package computes the Frenet apparatus: unit tangent, principal
normal and binormal (`FrenetFrame`), signed `Curvature` and `Torsion`, the
curvature vector, radius of curvature and `ArcLength`.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/diffgeo"
```

## Usage

This is the package's own `ExampleFrenetFrame`, so it compiles and its output is
asserted on every `go test ./diffgeo/`.

```go
f := diffgeo.FrenetFrame(diffgeo.Circle(1), 0)
	fmt.Printf("T=(%.3f,%.3f,%.3f)\n", nz(f.T.X), nz(f.T.Y), nz(f.T.Z))
	fmt.Printf("N=(%.3f,%.3f,%.3f)\n", nz(f.N.X), nz(f.N.Y), nz(f.N.Z))
	fmt.Printf("B=(%.3f,%.3f,%.3f)\n", nz(f.B.X), nz(f.B.Y), nz(f.B.Z))
```

```
T=(0.000,1.000,0.000)
N=(-1.000,0.000,0.000)
B=(0.000,0.000,1.000)
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func ArcLength(c Curve, a, b float64, n int) float64` | ArcLength returns the arc length of the curve between parameters a and b, computed as the integral of `Speed` using composite Simpson's rule with n… |
| `func AreaElement(s Surface, u, v float64) float64` | AreaElement returns the scalar surface area element √(EG − F²) = \|rᵤ×rᵥ\| at (u, v), the local factor relating parameter area to surface… |
| `func Curvature(c Curve, t float64) float64` | Curvature returns the (unsigned) curvature κ(t) = \|r'×r”\| / \|r'\|³, the reciprocal of the radius of the osculating circle. |
| `func GaussianCurvature(s Surface, u, v float64) float64` | GaussianCurvature returns the Gaussian curvature K = (LN − M²)/(EG − F²) at (u, v). |
| `func GeodesicAcceleration(s Surface, st GeodesicState) (accU, accV float64)` | GeodesicAcceleration returns the parameter-space acceleration (d²u, d²v) dictated by the geodesic equations at state st on surface s: |
| `func GeodesicCurvature(s Surface, st GeodesicState, accU, accV float64) float64` | GeodesicCurvature returns the geodesic curvature κ_g of the curve traced in the surface by the parameter-space path (u(t), v(t)) at the given state. |
| `func GeodesicLength(s Surface, path []GeodesicState) float64` | GeodesicLength returns the arc length in space of the geodesic sampled by `GeodesicPath`, measured as the total length of the polyline through the… |
| `func IsUmbilic(s Surface, u, v, tol float64) bool` | IsUmbilic reports whether (u, v) is an umbilic point, that is whether the two `PrincipalCurvatures` agree to within tol. |
| `func MeanCurvature(s Surface, u, v float64) float64` | MeanCurvature returns the mean curvature H = (EN + GL − 2FM)/(2(EG − F²)) at (u, v), the average of the two `PrincipalCurvatures`. |
| `func NormalCurvature(s Surface, u, v, du, dv float64) float64` | NormalCurvature returns the normal curvature at (u, v) in the tangent direction (du, dv), that is II(du,dv)/I(du,dv) = (L du² + 2M du dv + N… |
| `func PrincipalCurvatures(s Surface, u, v float64) (k1, k2 float64)` | PrincipalCurvatures returns the principal curvatures k1 ≥ k2 at (u, v), the extreme values of the `NormalCurvature` over all tangent directions. |
| `func PrincipalDirections(s Surface, u, v float64) (d1, d2 [2]float64)` | PrincipalDirections returns the two principal directions at (u, v) as parameter-space components (du, dv), the eigenvectors of the `ShapeOperator`… |
| `func RadiusOfCurvature(c Curve, t float64) float64` | RadiusOfCurvature returns 1/κ(t), the radius of the osculating circle. |
| `func ShapeOperator(s Surface, u, v float64) [2][2]float64` | ShapeOperator returns the matrix of the Weingarten map (shape operator) S = I⁻¹·II at (u, v), expressed in the parameter basis {rᵤ, rᵥ} and… |
| `func Speed(c Curve, t float64) float64` | Speed returns the parametric speed \|r'(t)\|, the rate at which arc length accumulates per unit parameter. |
| `func SurfaceArea(s Surface, u0, u1, v0, v1 float64, n int) float64` | SurfaceArea returns the area of the patch over the parameter rectangle [u0,u1]×[v0,v1], the double integral of the `AreaElement` evaluated by a… |
| `func Torsion(c Curve, t float64) float64` | Torsion returns the torsion τ(t) = ((r'×r”)·r”') / \|r'×r”\|², which measures how sharply the curve twists out of its osculating plane. |
| `func TotalCurvature(c Curve, a, b float64, n int) float64` | TotalCurvature returns the integral of the curvature with respect to arc length, ∫ κ(t)\|r'(t)\| dt, over the parameter interval [a, b], using… |
| `func TotalTorsion(c Curve, a, b float64, n int) float64` | TotalTorsion returns the integral of the torsion with respect to arc length, ∫ τ(t)\|r'(t)\| dt, over [a, b], using composite Simpson's rule with… |

### Types

| Type | What it is |
| --- | --- |
| `Christoffel` | Christoffel holds the Christoffel symbols of the second kind Γᵏᵢⱼ at a point, indexed as Gamma`k``i``j` with parameter indices 0 = u and 1 = v. |
| `Curve` | Curve is a parametric space curve: a function mapping a scalar parameter t to a point r(t) in three-dimensional space. |
| `FirstForm` | FirstForm is the first fundamental form (metric) of a surface at a point, I = E du² + 2F du dv + G dv², with coefficients E = rᵤ·rᵤ, F =… |
| `Frame` | Frame is an orthonormal Frenet frame at a point of a curve: the unit tangent T, principal normal N and binormal B, forming a right-handed basis with… |
| `GeodesicState` | GeodesicState is a point on a curve in a surface's parameter domain together with its parametric velocity: the parameters (U, V) and their rates of… |
| `SecondForm` | SecondForm is the second fundamental form of a surface at a point, II = L du² + 2M du dv + N dv², with coefficients L = rᵤᵤ·n, M = rᵤᵥ·n… |
| `Surface` | Surface is a parametric surface patch: a function mapping parameters (u, v) to a point r(u, v) in three-dimensional space. |
| `Vec3` | Vec3 is a vector in three-dimensional Euclidean space. |

<details>
<summary><code>Christoffel</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ChristoffelSymbols(s Surface, u, v float64) Christoffel` | ChristoffelSymbols returns the `Christoffel` symbols of the second kind at (u, v). |
| `func (c Christoffel) At(k, i, j int) float64` | At returns the symbol Γᵏᵢⱼ = Gamma`k``i``j`. |

</details>

<details>
<summary><code>Curve</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Circle(radius float64) Curve` | Circle returns the planar circle of the given radius in the xy-plane, r(t) = (radius·cos t, radius·sin t, 0). |
| `func Ellipse(a, b float64) Curve` | Ellipse returns the planar ellipse r(t) = (a·cos t, b·sin t, 0) with the given semi-axes a and b. |
| `func Helix(radius, rise float64) Curve` | Helix returns the circular helix r(t) = (radius·cos t, radius·sin t, rise·t). |
| `func Lemniscate(a float64) Curve` | Lemniscate returns the lemniscate of Bernoulli with half-width a, the planar figure-eight r(t) = a·(cos t, sin t·cos t, 0)/(1 + sin²t). |
| `func Line(p, d Vec3) Curve` | Line returns the straight-line curve r(t) = p + t·d through point p with direction d. |
| `func TorusKnot(p, q int, major, minor float64) Curve` | TorusKnot returns the (p, q) torus knot lying on the torus with center-circle radius major and tube radius minor, winding p times around the axis of… |
| `func Viviani(a float64) Curve` | Viviani returns Viviani's curve of radius a, the intersection of the sphere of radius 2a centred at the origin with the cylinder of radius a through… |

</details>

<details>
<summary><code>FirstForm</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func FirstFundamental(s Surface, u, v float64) FirstForm` | FirstFundamental returns the `FirstForm` coefficients at (u, v). |
| `func (f FirstForm) Angle(a0, a1, b0, b1 float64) float64` | Angle returns the angle in radians between two tangent vectors given by their parameter-space components, measured in the metric. |
| `func (f FirstForm) Determinant() float64` | Determinant returns EG − F², the squared area of the parallelogram spanned by rᵤ and rᵥ. |
| `func (f FirstForm) InnerProduct(a0, a1, b0, b1 float64) float64` | InnerProduct returns the metric inner product of two tangent vectors given by their parameter-space components, ⟨(a0,a1),(b0,b1)⟩ = E a0 b0 + F… |
| `func (f FirstForm) Norm(w0, w1 float64) float64` | Norm returns the length √(E w0² + 2F w0 w1 + G w1²) of the tangent vector with parameter-space components (w0, w1), measured in the metric. |

</details>

<details>
<summary><code>Frame</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func FrenetFrame(c Curve, t float64) Frame` | FrenetFrame returns the `Frame` {T, N, B} at parameter t. |

</details>

<details>
<summary><code>GeodesicState</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func GeodesicPath(s Surface, start GeodesicState, tEnd float64, steps int) []GeodesicState` | GeodesicPath integrates a geodesic on surface s from the initial `GeodesicState` start over a curve-parameter length of tEnd, using steps equal-sized… |
| `func GeodesicStep(s Surface, st GeodesicState, h float64) GeodesicState` | GeodesicStep advances one geodesic `GeodesicState` by a step h of the curve parameter using the classical fourth-order Runge-Kutta method. |
| `func ParallelTransport(s Surface, start GeodesicState, w0, w1, tEnd float64, steps int) ([]GeodesicState, [][2]float64)` | ParallelTransport transports the tangent vector with parameter-space components (w0, w1) along the geodesic starting at start on surface s, over a… |

</details>

<details>
<summary><code>SecondForm</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func SecondFundamental(s Surface, u, v float64) SecondForm` | SecondFundamental returns the `SecondForm` coefficients at (u, v). |

</details>

<details>
<summary><code>Surface</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Catenoid(c float64) Surface` | Catenoid returns the catenoid of waist radius c, the surface of revolution of a catenary: r(u, v) = (c·cosh(v/c)·cos u, c·cosh(v/c)·sin u, v). |
| `func Cylinder(radius float64) Surface` | Cylinder returns the parametrization of a circular cylinder of the given radius aligned with the z-axis: r(u,v) = (radius·cos u, radius·sin u, v). |
| `func Graph(f func(u, v float64) float64) Surface` | Graph returns the Monge patch r(u,v) = (u, v, f(u,v)) of the height function f, representing the surface z = f(x, y). |
| `func Helicoid(c float64) Surface` | Helicoid returns the helicoid of pitch c, r(u, v) = (u·cos v, u·sin v, c·v), the surface swept by a horizontal line rotating uniformly as it rises. |
| `func MonkeySaddle() Surface` | MonkeySaddle returns the monkey saddle z = x³ − 3xy² as the Monge patch r(u, v) = (u, v, u³ − 3uv²). |
| `func Paraboloid(a, b float64) Surface` | Paraboloid returns the elliptic paraboloid z = a·x² + b·y² as the Monge patch r(u, v) = (u, v, a·u² + b·v²). |
| `func Sphere(radius float64) Surface` | Sphere returns the parametrization of a sphere of the given radius by longitude u and latitude v: r(u,v) = radius·(cos v cos u, cos v sin u, sin v),… |
| `func SurfaceOfRevolution(radius, height func(v float64) float64) Surface` | SurfaceOfRevolution returns the surface swept by revolving the profile curve (radius(v), height(v)) in the xz-half-plane about the z-axis: |
| `func SurfacePlane(origin, du, dv Vec3) Surface` | SurfacePlane returns the affine plane r(u,v) = origin + u·du + v·dv spanned by the direction vectors du and dv through origin. |
| `func Torus(major, minor float64) Surface` | Torus returns the parametrization of a torus with center-circle radius major and tube radius minor: r(u,v) = ((major+minor·cos v)cos u,… |

</details>

<details>
<summary><code>Vec3</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Acceleration(c Curve, t float64) Vec3` | Acceleration returns the second derivative r”(t) by a central difference. |
| `func Binormal(c Curve, t float64) Vec3` | Binormal returns the unit binormal B(t) = T×N = (r'×r”)/\|r'×r”\|, the normal to the osculating plane. |
| `func CurvatureVector(c Curve, t float64) Vec3` | CurvatureVector returns dT/ds, the derivative of the unit tangent with respect to arc length. |
| `func Jerk(c Curve, t float64) Vec3` | Jerk returns the third derivative r”'(t), used by `Torsion`, via a five-point central difference. |
| `func NewVec3(x, y, z float64) Vec3` | NewVec3 returns the vector with the given components. |
| `func PartialU(s Surface, u, v float64) Vec3` | PartialU returns the first partial derivative rᵤ = ∂r/∂u at (u, v). |
| `func PartialUU(s Surface, u, v float64) Vec3` | PartialUU returns the second partial derivative rᵤᵤ = ∂²r/∂u² at (u, v). |
| `func PartialUV(s Surface, u, v float64) Vec3` | PartialUV returns the mixed second partial derivative rᵤᵥ = ∂²r/∂u∂v at (u, v), computed by the four-point central difference. |
| `func PartialV(s Surface, u, v float64) Vec3` | PartialV returns the first partial derivative rᵥ = ∂r/∂v at (u, v). |
| `func PartialVV(s Surface, u, v float64) Vec3` | PartialVV returns the second partial derivative rᵥᵥ = ∂²r/∂v² at (u, v). |
| `func PrincipalNormal(c Curve, t float64) Vec3` | PrincipalNormal returns the unit principal normal N(t), the direction in which the curve is turning. |
| `func SurfaceNormal(s Surface, u, v float64) Vec3` | SurfaceNormal returns the unit normal vector n = (rᵤ×rᵥ)/\|rᵤ×rᵥ\| at (u, v). |
| `func UnitTangent(c Curve, t float64) Vec3` | UnitTangent returns the unit tangent vector T(t) = r'(t)/\|r'(t)\|. |
| `func Velocity(c Curve, t float64) Vec3` | Velocity returns the first derivative r'(t), the velocity vector, computed by a central difference. |
| `func (v Vec3) Add(w Vec3) Vec3` | Add returns the component-wise sum v+w. |
| `func (v Vec3) Angle(w Vec3) float64` | Angle returns the unsigned angle in radians between v and w, in [0, π]. |
| `func (v Vec3) Cross(w Vec3) Vec3` | Cross returns the cross product v×w, a vector orthogonal to both operands. |
| `func (v Vec3) Distance(w Vec3) float64` | Distance returns the Euclidean distance between the points v and w. |
| `func (v Vec3) Dot(w Vec3) float64` | Dot returns the Euclidean inner product v·w. |
| `func (v Vec3) Equal(w Vec3, eps float64) bool` | Equal reports whether v and w agree component-wise within the tolerance eps. |
| `func (v Vec3) IsZero(eps float64) bool` | IsZero reports whether every component of v is within eps of zero. |
| `func (v Vec3) Lerp(w Vec3, t float64) Vec3` | Lerp returns the linear interpolation (1-t)·v + t·w. |
| `func (v Vec3) Neg() Vec3` | Neg returns the additive inverse -v. |
| `func (v Vec3) Norm() float64` | Norm returns the Euclidean length \|v\|. |
| `func (v Vec3) Norm2() float64` | Norm2 returns the squared Euclidean length \|v\|², avoiding the square root. |
| `func (v Vec3) Normalize() Vec3` | Normalize returns the unit vector pointing along v. |
| `func (v Vec3) Scale(s float64) Vec3` | Scale returns v with every component multiplied by the scalar s. |
| `func (v Vec3) Sub(w Vec3) Vec3` | Sub returns the component-wise difference v-w. |

</details>

### Constants

`Eps`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/diffgeo).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
