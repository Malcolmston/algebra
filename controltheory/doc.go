// Package controltheory provides classical linear control-systems primitives
// implemented with the Go standard library only.
//
// The package models single-input single-output (SISO) linear time-invariant
// systems in two equivalent representations:
//
//   - Transfer functions G(s) = N(s)/D(s), where N and D are real polynomials
//     stored as [Poly] values.
//   - State-space realizations (A, B, C, D) stored as [StateSpace] values.
//
// On top of these representations it offers block-diagram algebra
// (series, parallel, feedback), pole/zero extraction, time-domain
// step and impulse responses computed by numerical integration,
// controllability and observability analysis, Routh-Hurwitz stability
// testing, frequency-domain Bode and Nyquist sampling with gain and phase
// margins, and a discrete PID controller.
//
// Polynomials use the ascending-power convention: a [Poly] value p has
// p[i] as the coefficient of s^i, so Poly{2, 3, 1} represents 2 + 3s + s^2.
//
// Every routine is deterministic and depends on nothing outside the standard
// library.
package controltheory
