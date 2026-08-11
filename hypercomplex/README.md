# hypercomplex — implements arithmetic over several hypercomplex number systems that extend the real

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/hypercomplex.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/hypercomplex)

Package hypercomplex implements arithmetic over several hypercomplex number
systems that extend the real and complex numbers.

- Quaternion: Hamilton's quaternions, the four-dimensional normed division
algebra. Includes multiplication, conjugation, norm, inverse, the exponential
and logarithm maps, power, spherical linear interpolation (Slerp), rotation of
3-vectors, and conversions to and from axis-angle, Euler angles, and rotation
matrices. - DualQuaternion: the tensor of the dual numbers with the
quaternions, used to represent rigid-body transforms (screw motions) in a
single algebraic object. Includes composition, inversion, point transformation
and the screw decomposition into axis, angle and pitch. - Octonion: the
eight-dimensional non-associative normed division algebra, built from
quaternions by the Cayley-Dickson construction. - SplitComplex: the
split-complex (hyperbolic) numbers x + y*j with j*j = +1, the algebra
underlying Lorentz boosts in 1+1 dimensions.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/hypercomplex"
```

## Exported surface

### Types

| Type | What it is |
| --- | --- |
| `DualQuaternion` | DualQuaternion is a dual quaternion Real + Dual*epsilon, where epsilon is the dual unit with epsilon^2 = 0 and Real and Dual are ordinary… |
| `Octonion` | Octonion is an element of the eight-dimensional real algebra of octonions, stored as the coefficients E0..E7 of the basis units e0=1, e1..e7. |
| `Quaternion` | Quaternion is a Hamilton quaternion q = W + X*i + Y*j + Z*k, where i, j and k satisfy i*i = j*j = k*k = i*j*k = -1. |
| `SplitComplex` | SplitComplex is a split-complex (hyperbolic) number z = Re + Im*j, where the unit j satisfies j*j = +1 (in contrast to the ordinary imaginary unit). |
| `Vec3` | Vec3 is a three-dimensional real vector. |

<details>
<summary><code>DualQuaternion</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DualQuat(real, dual Quaternion) DualQuaternion` | DualQuat constructs a dual quaternion from its real and dual quaternion parts. |
| `func DualQuatFromRotationTranslation(rot Quaternion, t Vec3) DualQuaternion` | DualQuatFromRotationTranslation returns the unit dual quaternion for the rigid transform that first rotates by the (unit) quaternion rot and then… |
| `func DualQuatFromTranslation(t Vec3) DualQuaternion` | DualQuatFromTranslation returns the unit dual quaternion representing a pure translation by t with no rotation. |
| `func IdentityDualQuat() DualQuaternion` | IdentityDualQuat returns the identity dual quaternion (identity rotation, zero translation). |
| `func (d DualQuaternion) Add(e DualQuaternion) DualQuaternion` | Add returns the component-wise sum of the two dual quaternions. |
| `func (d DualQuaternion) Conjugate() DualQuaternion` | Conjugate returns the quaternion conjugate applied to both parts, (Real*, Dual*). |
| `func (d DualQuaternion) DualNumberConjugate() DualQuaternion` | DualNumberConjugate returns the dual-number conjugate (Real, -Dual), which negates the dual part while leaving the real part unchanged. |
| `func (d DualQuaternion) Equal(e DualQuaternion, tol float64) bool` | Equal reports whether d and e agree in both parts to within the absolute tolerance tol. |
| `func (d DualQuaternion) Inverse() DualQuaternion` | Inverse returns the multiplicative inverse of d, the dual quaternion e with d*e = e*d = 1. |
| `func (d DualQuaternion) Mul(e DualQuaternion) DualQuaternion` | Mul returns the dual quaternion product d*e. |
| `func (d DualQuaternion) Neg() DualQuaternion` | Neg returns the additive inverse -d, negating both parts. |
| `func (d DualQuaternion) Norm() float64` | Norm returns the norm of a unit dual quaternion's real part, \|Real\|. |
| `func (d DualQuaternion) Normalize() DualQuaternion` | Normalize returns d rescaled so that its real part has unit norm, dividing both parts by \|Real\|. |
| `func (d DualQuaternion) Rotation() Quaternion` | Rotation returns the rotation component of d as a unit quaternion (its normalized real part). |
| `func (d DualQuaternion) Scale(s float64) DualQuaternion` | Scale returns d with both parts multiplied by the real factor s. |
| `func (d DualQuaternion) ScrewAngle() float64` | ScrewAngle returns the rotation angle in radians of the rigid transform encoded by d, in the range [0, pi]. |
| `func (d DualQuaternion) ScrewAxis() Vec3` | ScrewAxis returns the unit rotation axis of the rigid transform encoded by d, the axis of the screw motion. |
| `func (d DualQuaternion) ScrewPitch() float64` | ScrewPitch returns the translation of the rigid transform encoded by d along its own rotation axis, the pitch of the screw motion. |
| `func (d DualQuaternion) Sub(e DualQuaternion) DualQuaternion` | Sub returns the component-wise difference of the two dual quaternions. |
| `func (d DualQuaternion) TransformPoint(p Vec3) Vec3` | TransformPoint applies the rigid transform encoded by the unit dual quaternion d to the point p, returning the rotated-then-translated point. |
| `func (d DualQuaternion) Translation() Vec3` | Translation returns the translation vector encoded by the unit dual quaternion d, recovered as 2 * Dual * Real*. |

</details>

<details>
<summary><code>Octonion</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Associator(o, p, q Octonion) Octonion` | Associator returns the associator (o*p)*q - o*(p*q), a measure of the failure of associativity. |
| `func Commutator(o, p Octonion) Octonion` | Commutator returns the commutator o*p - p*o, which vanishes exactly when o and p commute. |
| `func IdentityOct() Octonion` | IdentityOct returns the multiplicative identity octonion e0 = 1. |
| `func Oct(e0, e1, e2, e3, e4, e5, e6, e7 float64) Octonion` | Oct constructs an octonion from its eight real coefficients. |
| `func OctFromScalar(s float64) Octonion` | OctFromScalar returns the octonion whose real part is s and whose other seven coefficients are zero. |
| `func (o Octonion) Add(p Octonion) Octonion` | Add returns the component-wise sum o + p. |
| `func (o Octonion) Conj() Octonion` | Conj returns the octonion conjugate, which negates all seven imaginary coefficients and leaves the real part unchanged. |
| `func (o Octonion) Dot(p Octonion) float64` | Dot returns the Euclidean dot product of o and p treated as 8-vectors. |
| `func (o Octonion) Equal(p Octonion, tol float64) bool` | Equal reports whether o and p agree in every coefficient to within the absolute tolerance tol. |
| `func (o Octonion) Imag() Octonion` | Imag returns the imaginary part of o, the octonion with the same seven imaginary coefficients and zero real part. |
| `func (o Octonion) Inverse() Octonion` | Inverse returns the multiplicative inverse conj(o)/\|o\|^2. |
| `func (o Octonion) IsUnit(tol float64) bool` | IsUnit reports whether o has unit norm to within the absolute tolerance tol. |
| `func (o Octonion) Mul(p Octonion) Octonion` | Mul returns the octonion product o*p computed by the Cayley-Dickson formula (a,b)(c,d) = (a*c - conj(d)*b, d*a + b*conj(c)) over the quaternions. |
| `func (o Octonion) Neg() Octonion` | Neg returns the additive inverse -o. |
| `func (o Octonion) Norm() float64` | Norm returns the Euclidean norm \|o\| = sqrt(NormSq). |
| `func (o Octonion) NormSq() float64` | NormSq returns the squared norm, the sum of the squares of the eight coefficients, equal to the real part of o * conj(o). |
| `func (o Octonion) Normalize() Octonion` | Normalize returns o scaled to unit norm. |
| `func (o Octonion) Real() float64` | Real returns the real (scalar) part E0 of o. |
| `func (o Octonion) Scale(s float64) Octonion` | Scale returns o with every coefficient multiplied by the real factor s. |
| `func (o Octonion) Sub(p Octonion) Octonion` | Sub returns the component-wise difference o - p. |

</details>

<details>
<summary><code>Quaternion</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func IdentityQuat() Quaternion` | IdentityQuat returns the multiplicative identity quaternion 1 + 0i + 0j + 0k. |
| `func Quat(w, x, y, z float64) Quaternion` | Quat constructs a quaternion from its scalar part w and vector part (x, y, z). |
| `func QuatBetweenVectors(from, to Vec3) Quaternion` | QuatBetweenVectors returns the unit quaternion representing the shortest-arc rotation that takes the direction of from to the direction of to. |
| `func QuatFromAxisAngle(axis Vec3, angle float64) Quaternion` | QuatFromAxisAngle returns the unit quaternion representing a rotation of angle radians about axis. |
| `func QuatFromEuler(roll, pitch, yaw float64) Quaternion` | QuatFromEuler returns the unit quaternion for the intrinsic Tait-Bryan rotation with the given roll (about X), pitch (about Y) and yaw (about Z),… |
| `func QuatFromRotationMatrix(m [9]float64) Quaternion` | QuatFromRotationMatrix returns the unit quaternion corresponding to the row-major 3x3 rotation matrix m (nine elements m00, m01, m02, m10, ...). |
| `func QuatFromScalar(s float64) Quaternion` | QuatFromScalar returns the quaternion whose scalar part is s and whose vector part is zero. |
| `func QuatFromVector(v Vec3) Quaternion` | QuatFromVector returns the pure quaternion 0 + v (a quaternion with zero scalar part), used to embed a 3-vector into the quaternions. |
| `func Slerp(q0, q1 Quaternion, t float64) Quaternion` | Slerp returns the spherical linear interpolation between the unit quaternions q0 and q1 at parameter t in [0, 1]. |
| `func (q Quaternion) Add(r Quaternion) Quaternion` | Add returns the component-wise sum q + r. |
| `func (q Quaternion) AngleTo(r Quaternion) float64` | AngleTo returns the angle in radians of the relative rotation taking q to r, that is the rotation angle of r * q^-1, in the range [0, pi]. |
| `func (q Quaternion) Conj() Quaternion` | Conj returns the conjugate q* = W - X*i - Y*j - Z*k. |
| `func (q Quaternion) Div(r Quaternion) Quaternion` | Div returns the right quotient q * r^-1. |
| `func (q Quaternion) Dot(r Quaternion) float64` | Dot returns the Euclidean dot product of q and r treated as 4-vectors. |
| `func (q Quaternion) Equal(r Quaternion, tol float64) bool` | Equal reports whether q and r agree component-wise to within the absolute tolerance tol. |
| `func (q Quaternion) Exp() Quaternion` | Exp returns the quaternion exponential exp(q). |
| `func (q Quaternion) Inverse() Quaternion` | Inverse returns the multiplicative inverse q^-1 = q* / \|q\|^2. |
| `func (q Quaternion) IsUnit(tol float64) bool` | IsUnit reports whether q has unit norm to within the absolute tolerance tol. |
| `func (q Quaternion) Log() Quaternion` | Log returns the principal quaternion logarithm log(q). |
| `func (q Quaternion) Mul(r Quaternion) Quaternion` | Mul returns the Hamilton product q*r. |
| `func (q Quaternion) Neg() Quaternion` | Neg returns the additive inverse -q. |
| `func (q Quaternion) Norm() float64` | Norm returns the Euclidean norm \|q\| = sqrt(W^2 + X^2 + Y^2 + Z^2). |
| `func (q Quaternion) NormSq() float64` | NormSq returns the squared norm W^2 + X^2 + Y^2 + Z^2, equal to the scalar part of q * q*. |
| `func (q Quaternion) Normalize() Quaternion` | Normalize returns q scaled to unit norm. |
| `func (q Quaternion) Pow(t float64) Quaternion` | Pow returns q raised to the real power t, computed as exp(t*log(q)). |
| `func (q Quaternion) RotateVector(v Vec3) Vec3` | RotateVector returns the image of v under the rotation represented by q, computed as the vector part of q * (0,v) * q^-1. |
| `func (q Quaternion) Scalar() float64` | Scalar returns the scalar (real) part W of q. |
| `func (q Quaternion) Scale(s float64) Quaternion` | Scale returns q with every component multiplied by the real factor s. |
| `func (q Quaternion) Sub(r Quaternion) Quaternion` | Sub returns the component-wise difference q - r. |
| `func (q Quaternion) ToAxisAngle() (axis Vec3, angle float64)` | ToAxisAngle returns the rotation axis (a unit vector) and angle in radians (in the range [0, pi]) represented by q. |
| `func (q Quaternion) ToEuler() (roll, pitch, yaw float64)` | ToEuler returns the roll (about X), pitch (about Y) and yaw (about Z) angles in radians that reproduce the rotation of q under the aerospace Z-Y-X… |
| `func (q Quaternion) ToRotationMatrix() [9]float64` | ToRotationMatrix returns the row-major 3x3 rotation matrix (nine elements m00, m01, m02, m10, ...) equivalent to the rotation represented by q. |
| `func (q Quaternion) Vector() Vec3` | Vector returns the vector (imaginary) part (X, Y, Z) of q as a `Vec3`. |

</details>

<details>
<summary><code>SplitComplex</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Split(re, im float64) SplitComplex` | Split constructs a split-complex number Re + Im*j. |
| `func SplitFromLightCone(u, v float64) SplitComplex` | SplitFromLightCone reconstructs a split-complex number from its null-basis coordinates (u, v), returning ((u+v)/2) + ((u-v)/2)*j. |
| `func SplitFromModulusArgument(r, phi float64) SplitComplex` | SplitFromModulusArgument returns the time-like split-complex number with the given non-negative modulus r and rapidity phi, namely r*(cosh phi + j… |
| `func (z SplitComplex) Abs() float64` | Abs returns sqrt(\|ModulusSq\|), a non-negative magnitude that is multiplicative in absolute value. |
| `func (z SplitComplex) Add(w SplitComplex) SplitComplex` | Add returns the sum z + w. |
| `func (z SplitComplex) Argument() float64` | Argument returns the hyperbolic argument (rapidity) atanh(Im/Re) of z. |
| `func (z SplitComplex) Conj() SplitComplex` | Conj returns the conjugate Re - Im*j. |
| `func (z SplitComplex) Div(w SplitComplex) (SplitComplex, bool)` | Div returns the quotient z/w together with a boolean that is false when w is a zero divisor (its squared modulus vanishes), in which case the… |
| `func (z SplitComplex) Equal(w SplitComplex, tol float64) bool` | Equal reports whether z and w agree in both parts to within the absolute tolerance tol. |
| `func (z SplitComplex) Exp() SplitComplex` | Exp returns the split-complex exponential exp(z) = e^Re (cosh Im + j sinh Im). |
| `func (z SplitComplex) Inverse() (SplitComplex, bool)` | Inverse returns the multiplicative inverse conj(z)/ModulusSq, together with a boolean that is false when z is a zero divisor (ModulusSq = 0), in… |
| `func (z SplitComplex) IsZeroDivisor(tol float64) bool` | IsZeroDivisor reports whether z lies on the null cone Re = ±Im (with z not the zero element) to within the absolute tolerance tol, in which case z… |
| `func (z SplitComplex) LightCone() (u, v float64)` | LightCone returns the coordinates (u, v) of z in the null (diagonal) basis, with u = Re + Im and v = Re - Im. |
| `func (z SplitComplex) Log() (SplitComplex, bool)` | Log returns the split-complex logarithm ln(z) = ln\|z\| + j*atanh(Im/Re), the inverse of `SplitComplex.Exp` on the right time-like sector Re >… |
| `func (z SplitComplex) ModulusSq() float64` | ModulusSq returns the split-complex squared modulus Re^2 - Im^2, equal to the real part of z * conj(z). |
| `func (z SplitComplex) Mul(w SplitComplex) SplitComplex` | Mul returns the product z*w. |
| `func (z SplitComplex) Neg() SplitComplex` | Neg returns the additive inverse -z. |
| `func (z SplitComplex) Scale(s float64) SplitComplex` | Scale returns z with both parts multiplied by the real factor s. |
| `func (z SplitComplex) Sub(w SplitComplex) SplitComplex` | Sub returns the difference z - w. |

</details>

<details>
<summary><code>Vec3</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func V3(x, y, z float64) Vec3` | V3 constructs a `Vec3` from its three components. |
| `func (v Vec3) Add(w Vec3) Vec3` | Add returns the component-wise sum v + w. |
| `func (v Vec3) Cross(w Vec3) Vec3` | Cross returns the cross product v×w. |
| `func (v Vec3) Dot(w Vec3) float64` | Dot returns the Euclidean dot product v·w. |
| `func (v Vec3) Equal(w Vec3, tol float64) bool` | Equal reports whether v and w agree component-wise to within the absolute tolerance tol. |
| `func (v Vec3) Neg() Vec3` | Neg returns the additive inverse -v. |
| `func (v Vec3) Norm() float64` | Norm returns the Euclidean length \|v\|. |
| `func (v Vec3) Normalize() Vec3` | Normalize returns v scaled to unit length. |
| `func (v Vec3) Scale(s float64) Vec3` | Scale returns the vector v scaled by the real factor s. |
| `func (v Vec3) Sub(w Vec3) Vec3` | Sub returns the component-wise difference v - w. |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/hypercomplex).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
