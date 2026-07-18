// Package crypto implements the mathematical primitives that underpin
// classical public-key cryptography, built exclusively on top of
// math/big.Int and the Go standard library.
//
// The package is deliberately educational and self-contained: it exposes the
// number-theoretic machinery directly rather than hiding it behind opaque,
// hardened implementations. Nothing here draws on crypto/rand, crypto/rsa or
// any other package under crypto/*; every routine is expressed in terms of
// arbitrary-precision integer arithmetic so that the algorithms can be read,
// tested against closed-form values and reused for teaching.
//
// The following families are covered:
//
//   - modular arithmetic: fast modular exponentiation (ModExp), modular
//     inverse, ModMul/ModAdd/ModSub, GCD/LCM, the extended Euclidean
//     algorithm and the Chinese Remainder Theorem (CRT, CRTPair);
//   - quadratic residues: the Jacobi and Legendre symbols, a Tonelli-Shanks
//     modular square root, multiplicative order and primitive roots;
//   - primality: the Fermat and Miller-Rabin tests (both probabilistic and a
//     deterministic witness set for moderate inputs), trial division, the
//     sieve of Eratosthenes, next/previous prime and random (safe) prime
//     generation;
//   - factoring: Pollard's rho (classic and Brent variants), full integer
//     factorization, Euler's totient and the Carmichael function;
//   - discrete logarithms: the baby-step giant-step algorithm;
//   - protocols: RSA key generation, encryption, decryption (with a CRT fast
//     path) and signatures; Diffie-Hellman key agreement; ElGamal
//     encryption; and Shamir's threshold secret sharing.
//
// Randomness. Routines that need randomness accept an explicit *math/rand.Rand
// source. This keeps the whole package deterministic and reproducible: seeding
// the source identically always yields identical keys, primes and shares. Such
// determinism is invaluable for tests and demonstrations but means the package
// is NOT suitable for producing real-world secrets. Do not use it to protect
// anything of value.
//
// Conventions. Moduli are required to be positive. Functions that cannot
// return a meaningful result for out-of-domain integer inputs (a non-positive
// modulus, a negative factorial argument and so on) panic with a message
// prefixed by "crypto:"; functions whose failure depends on runtime values the
// caller cannot cheaply check in advance (a non-invertible element, a discrete
// log that does not exist) return an error instead.
package crypto
