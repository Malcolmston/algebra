// Package bigfloat implements arbitrary-precision real analysis on top of
// math/big.Float.
//
// Everything in the package is computed to a caller-requested precision,
// expressed in bits of mantissa. Each routine accepts a prec argument and
// returns a freshly allocated *big.Float carrying exactly that precision,
// rounded to nearest. Internally the routines work at a slightly higher
// "guard" precision and, where the argument magnitude demands it, add extra
// bits so that argument reduction does not lose significance; the final
// result is then rounded down to the requested precision.
//
// The package is organised in families:
//
//   - Mathematical constants to N bits: [Pi], [E], [Ln2], [Ln10],
//     [EulerGamma], [Catalan], [Apery], [GoldenRatio] and many more.
//   - Elementary transcendental functions built from range-reduced power
//     series and the arithmetic-geometric mean: [Exp], [Log], [Sqrt], [Cbrt],
//     [Pow], [Sin], [Cos], [Tan], [Atan], [Asin], [Acos], [Sinh], [Cosh],
//     [Tanh] and their relatives.
//   - The gamma family: [Factorial], [Gamma], [Lgamma], [Digamma], [Beta] and
//     the Bernoulli numbers via [Bernoulli].
//   - Newton-based inverse functions such as [LambertW], [PlasticNumber] and
//     the AGM-based elliptic integrals [EllipticK] and [EllipticE].
//   - Utility layers for construction, conversion, comparison and rounding of
//     *big.Float values.
//
// Only the Go standard library is used. Functions whose mathematical domain is
// restricted (for example [Log] of a non-positive number, or [Asin] of an
// argument outside [-1, 1]) return an error rather than panicking; total
// functions return a *big.Float directly.
//
// The package is a self-contained sibling of github.com/malcolmston/algebra and
// does not depend on the parent module.
package bigfloat
