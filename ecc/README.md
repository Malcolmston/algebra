# ecc — implements elliptic-curve arithmetic in short Weierstrass form y^2 = x^3 + A*x + B

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/ecc.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/ecc)

Package ecc implements elliptic-curve arithmetic in short Weierstrass form y^2
= x^3 + A*x + B, both over a prime finite field GF(p) and over the field of
rational numbers Q.

The package provides the geometric group law (point addition, doubling,
negation and scalar multiplication), curve invariants (discriminant and
j-invariant), point-counting utilities (naive enumeration and a
baby-step/giant-step point-order routine), SEC1 point compression and
decompression, and the finite-field number theory required to build higher
level protocols on top: an Elliptic-Curve Diffie-Hellman shared-secret
derivation and the ECDSA signing and verification equations over GF(p).

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/ecc"
```

## Usage

This is the package's own `ExampleCurveFp_Add`, so it compiles and its output is
asserted on every `go test ./ecc/`.

```go
c, _ := ecc.NewCurveFp(big.NewInt(2), big.NewInt(3), big.NewInt(97))
	g, _ := c.NewPoint(big.NewInt(3), big.NewInt(6))
	two := c.Double(g)
	fmt.Printf("2G = (%v, %v)\n", two.X, two.Y)
	three := c.Add(two, g)
	fmt.Printf("3G = (%v, %v)\n", three.X, three.Y)

	fmt.Println(c.Add(g, c.Negate(g)).Infinity)

	fmt.Println(c.Equal(c.Add(g, c.Identity()), g))
```

```
2G = (80, 10)
3G = (80, 87)
true
true
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func DiscriminantFp(a, b, p *big.Int) *big.Int` | DiscriminantFp returns the discriminant -16*(4*a^3 + 27*b^2) of the curve y^2 = x^3 + a*x + b modulo p, without constructing a CurveFp. |
| `func IsQuadraticResidue(a, p *big.Int) bool` | IsQuadraticResidue reports whether a is a non-zero quadratic residue modulo the odd prime p, i.e. |
| `func JInvariantFp(a, b, p *big.Int) (*big.Int, bool)` | JInvariantFp returns the j-invariant 1728 * 4*a^3 / (4*a^3 + 27*b^2) of the curve y^2 = x^3 + a*x + b modulo p, without constructing a CurveFp,… |
| `func Legendre(a, p *big.Int) int` | Legendre returns the Legendre symbol (a / p) for an odd prime p, one of the values -1, 0 or +1: 0 when p divides a, +1 when a is a non-zero quadratic… |
| `func Mod(a, p *big.Int) *big.Int` | Mod returns the representative of a modulo p in the range [0, p). |
| `func ModAdd(a, b, p *big.Int) *big.Int` | ModAdd returns (a + b) mod p in the range [0, p) for a positive modulus p. |
| `func ModDiv(a, b, p *big.Int) (*big.Int, bool)` | ModDiv returns (a / b) mod p, that is a times the modular inverse of b, along with a boolean reporting whether b is invertible modulo p. |
| `func ModExp(base, exp, p *big.Int) *big.Int` | ModExp returns (base ^ exp) mod p for a positive modulus p. |
| `func ModInverse(a, p *big.Int) (*big.Int, bool)` | ModInverse returns the multiplicative inverse of a modulo p together with a boolean reporting whether the inverse exists. |
| `func ModMul(a, b, p *big.Int) *big.Int` | ModMul returns (a * b) mod p in the range [0, p) for a positive modulus p. |
| `func ModNeg(a, p *big.Int) *big.Int` | ModNeg returns (-a) mod p in the range [0, p) for a positive modulus p. |
| `func ModSqrt(a, p *big.Int) (*big.Int, bool)` | ModSqrt returns a square root of a modulo an odd prime p together with a boolean reporting whether a is a quadratic residue. |
| `func ModSub(a, b, p *big.Int) *big.Int` | ModSub returns (a - b) mod p in the range [0, p) for a positive modulus p. |

### Types

| Type | What it is |
| --- | --- |
| `CompressedPoint` | CompressedPoint holds the SEC1-style compressed representation of an affine point: its x-coordinate together with a single parity bit selecting one… |
| `CurveFp` | CurveFp is a short Weierstrass curve y^2 = x^3 + A*x + B defined over the prime field GF(P). |
| `CurveQ` | CurveQ is a short Weierstrass curve y^2 = x^3 + A*x + B defined over the rational numbers Q. The parameters must satisfy 4*A^3 + 27*B^2 != 0, which… |
| `ECDSASignature` | ECDSASignature holds the pair (R, S) produced by the ECDSA signing equation. |
| `NamedCurve` | NamedCurve bundles a curve over GF(p) with the standard domain parameters used by cryptographic protocols: a fixed base point G that generates a… |
| `PointFp` | PointFp is an affine point on a curve over GF(p). |
| `PointQ` | PointQ is an affine point on a curve over the rational field Q. The zero of the group, the point at infinity, is represented by Infinity set to true,… |

<details>
<summary><code>CurveFp</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewCurveFp(a, b, p *big.Int) (*CurveFp, error)` | NewCurveFp constructs a curve y^2 = x^3 + a*x + b over GF(p). |
| `func (c *CurveFp) Add(p, q PointFp) PointFp` | Add returns the sum p + q under the elliptic-curve group law. |
| `func (c *CurveFp) Cofactor(n *big.Int) (*big.Int, bool)` | Cofactor returns the cofactor #E(GF(p)) / n, the index of the order-n subgroup in the full group, together with a boolean that is false when n is… |
| `func (c *CurveFp) Compress(pt PointFp) CompressedPoint` | Compress returns the compressed representation of pt, storing its x-coordinate and the parity of its y-coordinate. |
| `func (c *CurveFp) CompressBytes(pt PointFp) []byte` | CompressBytes returns the SEC1 compressed octet string of pt: a single byte 0x00 for the point at infinity, or the prefix 0x02 (even y) or 0x03 (odd… |
| `func (c *CurveFp) CountPointsNaive() *big.Int` | CountPointsNaive returns the order #E(GF(p)) of the curve group, that is the number of affine points plus one for the point at infinity, by… |
| `func (c *CurveFp) Decompress(cp CompressedPoint) (PointFp, error)` | Decompress reconstructs the affine point from a compressed representation by solving y^2 = x^3 + A*x + B for y and selecting the root matching the… |
| `func (c *CurveFp) DecompressBytes(b []byte) (PointFp, error)` | DecompressBytes parses a SEC1 compressed octet string produced by CompressBytes and returns the corresponding affine point. |
| `func (c *CurveFp) Discriminant() *big.Int` | Discriminant returns the curve discriminant -16*(4*A^3 + 27*B^2) reduced modulo P. A curve is non-singular exactly when this value is non-zero. |
| `func (c *CurveFp) Double(pt PointFp) PointFp` | Double returns 2*pt using the tangent-line group law. |
| `func (c *CurveFp) ECDHSharedSecret(d *big.Int, q PointFp) (secret *big.Int, shared PointFp, ok bool)` | ECDHSharedSecret computes the Elliptic-Curve Diffie-Hellman shared secret d*Q, where d is the local private scalar and Q is the remote party's public… |
| `func (c *CurveFp) ECDSASign(g PointFp, n, d, e, k *big.Int) (ECDSASignature, error)` | ECDSASign produces an ECDSA signature (r, s) for the integer message digest e under the private key d, using the base point G of prime subgroup order… |
| `func (c *CurveFp) ECDSAVerify(g, q PointFp, n, e *big.Int, sig ECDSASignature) bool` | ECDSAVerify checks an ECDSA signature (r, s) against the message digest e and the public key Q, using the base point G of prime subgroup order n. |
| `func (c *CurveFp) Equal(p, q PointFp) bool` | Equal reports whether two points on the curve are equal, treating all point-at-infinity values as a single element. |
| `func (c *CurveFp) GenerateKeyPair(g PointFp, n *big.Int, reader io.Reader) (priv *big.Int, pub PointFp, err error)` | GenerateKeyPair draws a private scalar uniformly from [1, n) using the supplied reader and returns it together with the public point priv*G. It is… |
| `func (c *CurveFp) HasseInterval() (lo, hi *big.Int)` | HasseInterval returns the closed interval [lo, hi] of admissible group orders #E(GF(p)) guaranteed by the Hasse theorem: lo = p + 1 -… |
| `func (c *CurveFp) Identity() PointFp` | Identity returns the point at infinity, the neutral element of the curve group. |
| `func (c *CurveFp) IsGenerator(pt PointFp, n *big.Int) bool` | IsGenerator reports whether pt has order exactly n in the curve group, i.e. |
| `func (c *CurveFp) IsOnCurve(pt PointFp) bool` | IsOnCurve reports whether pt lies on the curve. |
| `func (c *CurveFp) IsSmooth() bool` | IsSmooth reports whether the curve is non-singular over GF(P), i.e. |
| `func (c *CurveFp) JInvariant() *big.Int` | JInvariant returns the j-invariant 1728 * 4*A^3 / (4*A^3 + 27*B^2) of the curve, reduced modulo P. Curves with equal j-invariant are isomorphic over… |
| `func (c *CurveFp) Negate(pt PointFp) PointFp` | Negate returns the additive inverse -pt, the reflection of pt across the x-axis. |
| `func (c *CurveFp) NewPoint(x, y *big.Int) (PointFp, error)` | NewPoint constructs the affine point (x, y) on the curve, reducing the coordinates modulo P. It returns an error if the point does not satisfy the… |
| `func (c *CurveFp) OnCurveX(x *big.Int) bool` | OnCurveX reports whether the field element x is the abscissa of some point on the curve, i.e. |
| `func (c *CurveFp) PointOrderBSGS(pt PointFp) *big.Int` | PointOrderBSGS returns the order of the point pt using the baby-step giant-step algorithm, which runs in time proportional to sqrt(#E) rather than… |
| `func (c *CurveFp) PointOrderNaive(pt PointFp) *big.Int` | PointOrderNaive returns the order of the point pt in the curve group: the least positive n with n*pt equal to the point at infinity. |
| `func (c *CurveFp) PublicKeyFp(g PointFp, d *big.Int) PointFp` | PublicKeyFp derives the public key d*G from a private scalar d and a base point G on the curve. |
| `func (c *CurveFp) RHS(x *big.Int) *big.Int` | RHS returns the right-hand side x^3 + A*x + B of the curve equation evaluated at x, reduced modulo P. |
| `func (c *CurveFp) SameJInvariant(other *CurveFp) bool` | SameJInvariant reports whether this curve and other share the same j-invariant over the same base field, a necessary condition for the two curves to… |
| `func (c *CurveFp) ScalarMul(k *big.Int, pt PointFp) PointFp` | ScalarMul returns k*pt using the constant-structure double-and-add algorithm. |
| `func (c *CurveFp) TraceOfFrobenius() *big.Int` | TraceOfFrobenius returns the trace of Frobenius t = p + 1 - #E(GF(p)), where #E is the group order. |

</details>

<details>
<summary><code>CurveQ</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewCurveQ(a, b *big.Rat) (*CurveQ, error)` | NewCurveQ constructs a curve y^2 = x^3 + a*x + b over Q. The curve must be non-singular (its discriminant must be non-zero); otherwise a non-nil… |
| `func (c *CurveQ) Add(p, q PointQ) PointQ` | Add returns the sum p + q over Q under the elliptic-curve group law. |
| `func (c *CurveQ) Discriminant() *big.Rat` | Discriminant returns the exact rational discriminant -16*(4*A^3 + 27*B^2) of the curve over Q. A curve is non-singular exactly when this value is… |
| `func (c *CurveQ) Double(pt PointQ) PointQ` | Double returns 2*pt over Q using the tangent-line group law. |
| `func (c *CurveQ) Equal(p, q PointQ) bool` | Equal reports whether two rational points are equal, treating all point-at-infinity values as a single element. |
| `func (c *CurveQ) Identity() PointQ` | Identity returns the point at infinity, the neutral element of the curve group over Q. |
| `func (c *CurveQ) IsOnCurve(pt PointQ) bool` | IsOnCurve reports whether the rational point pt lies on the curve. |
| `func (c *CurveQ) JInvariant() *big.Rat` | JInvariant returns the exact rational j-invariant 1728 * 4*A^3 / (4*A^3 + 27*B^2) of the curve over Q. It panics only if the curve is singular, which… |
| `func (c *CurveQ) Negate(pt PointQ) PointQ` | Negate returns the additive inverse -pt, the reflection of pt across the x-axis. |
| `func (c *CurveQ) NewPoint(x, y *big.Rat) (PointQ, error)` | NewPoint constructs the affine rational point (x, y) on the curve. |
| `func (c *CurveQ) RHS(x *big.Rat) *big.Rat` | RHS returns the exact value x^3 + A*x + B of the curve equation's right-hand side evaluated at the rational x. |
| `func (c *CurveQ) ScalarMul(k *big.Int, pt PointQ) PointQ` | ScalarMul returns k*pt over Q using double-and-add. |

</details>

<details>
<summary><code>NamedCurve</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func P256() NamedCurve` | P256 returns the NIST P-256 curve (also called secp256r1 or prime256v1), y^2 = x^3 - 3*x + B over a 256-bit prime field, together with its generator,… |
| `func Secp256k1() NamedCurve` | Secp256k1 returns the standard SEC 2 curve secp256k1, y^2 = x^3 + 7 over a 256-bit prime field, together with its generator, subgroup order and… |
| `func (nc NamedCurve) GenerateKeyPair(reader io.Reader) (priv *big.Int, pub PointFp, err error)` | GenerateKeyPair draws a private key for the named curve from the supplied reader and returns it together with the corresponding public point. |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/ecc).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
