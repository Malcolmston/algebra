// Package codingtheory implements classical error-correcting codes, source
// codes and the finite-field arithmetic they are built on.
//
// The package is self-contained and depends only on the Go standard library.
// It is organised into the following families:
//
//   - Galois-field arithmetic: carry-less polynomial arithmetic over GF(2)
//     (GF2Poly* helpers) and a table-driven extension field GF(2^m) exposed as
//     the Field type, with the usual add/mul/div/inv/pow/log operations and
//     ready-made constructors NewGF4 … NewGF256.
//   - Block codes over GF(2): general linear codes described by generator and
//     parity-check matrices (LinearCode), Hamming and extended-Hamming codes,
//     the perfect binary Golay(23,12) and extended Golay(24,12) codes,
//     repetition and parity-check codes, and Hadamard/Walsh codes.
//   - Cyclic and polynomial codes: generator-polynomial encoding, syndrome
//     computation, BCH codes with bounded-distance decoding, and Reed-Solomon
//     codes over GF(2^m) with Berlekamp-Massey / Chien / Forney decoding.
//   - Stream codes: configurable cyclic redundancy checks (CRC of arbitrary
//     width and polynomial) and rate-1/n convolutional codes with hard and
//     soft Viterbi decoding.
//   - Modulation and framing helpers: binary and multi-symbol Gray codes,
//     block and convolutional interleavers, and LDPC parity-check utilities.
//   - Source coding: Huffman coding and a fixed-point arithmetic coder, plus
//     Hamming-weight / Hamming-distance metrics.
//
// Unless documented otherwise, code words and messages are represented as
// []int slices of bits (each entry 0 or 1) or, for codes over GF(2^m), as
// []int slices of field elements in [0, 2^m). Polynomials over GF(2) are also
// represented compactly as plain ints whose bit i is the coefficient of x^i.
//
// Functions never mutate their slice arguments unless the documentation says
// so explicitly; results are freshly allocated.
package codingtheory
