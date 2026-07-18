package ntheory

import "sort"

// This file provides small unexported helpers shared by the uint64 fast-path
// routines (factorization, primitive roots) and the modular square-root code.
// They are collected here so the feature files that were added independently can
// rely on a single definition of each.

// PrimePowerU64 pairs a prime with its exponent in a uint64 factorization. It is
// the unsigned counterpart of [PrimePower].
type PrimePowerU64 struct {
	// Prime is the prime factor.
	Prime uint64
	// Exponent is the multiplicity of Prime in the factorization (>= 1).
	Exponent uint
}

// FactorListUint64 returns the prime factorization of n as a slice of
// [PrimePowerU64] sorted by ascending prime. It is the ordered counterpart of
// [FactorizeU64]. FactorListUint64(0) and FactorListUint64(1) return nil.
func FactorListUint64(n uint64) []PrimePowerU64 {
	factors := FactorizeU64(n)
	if len(factors) == 0 {
		return nil
	}
	list := make([]PrimePowerU64, 0, len(factors))
	for p, e := range factors {
		list = append(list, PrimePowerU64{Prime: p, Exponent: uint(e)})
	}
	sort.Slice(list, func(i, j int) bool { return list[i].Prime < list[j].Prime })
	return list
}

// ntheoryGCDU64 returns the greatest common divisor of two unsigned integers
// using the iterative Euclidean algorithm. ntheoryGCDU64(0, 0) is 0.
func ntheoryGCDU64(a, b uint64) uint64 {
	for b != 0 {
		a, b = b, a%b
	}
	return a
}

// ntheoryAbsDiffU64 returns |a - b| for unsigned integers without overflow.
func ntheoryAbsDiffU64(a, b uint64) uint64 {
	if a >= b {
		return a - b
	}
	return b - a
}

// ntheorySmallerRoot returns the smaller of the two square roots x and p-x of a
// residue modulo the odd prime p. Because a modular square root always comes in
// the pair {x, p-x}, returning the smaller one makes the result deterministic.
// x is assumed to lie in [0, p).
func ntheorySmallerRoot(x, p int64) int64 {
	if other := p - x; other < x {
		return other
	}
	return x
}

// ntheoryTrialPrimes lists the small primes used to strip easy factors before
// invoking the heavier Pollard rho/Brent routine. They are kept in ascending
// order so trial division can stop as soon as p*p exceeds the remaining
// cofactor.
var ntheoryTrialPrimes = []uint64{
	2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
	73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
}
