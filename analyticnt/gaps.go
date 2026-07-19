package analyticnt

import (
	"math"
	"math/rand"
)

// PrimeGaps returns the sequence of gaps p_{n+1} − p_n between consecutive
// primes up to x. The i-th entry is the gap after the i-th prime.
func PrimeGaps(x int) []int64 {
	primes := PrimesUpTo(x)
	if len(primes) < 2 {
		return []int64{}
	}
	gaps := make([]int64, len(primes)-1)
	for i := 1; i < len(primes); i++ {
		gaps[i-1] = primes[i] - primes[i-1]
	}
	return gaps
}

// MaximalPrimeGaps returns the record ("maximal") prime gaps up to x: each entry
// is a gap strictly larger than every earlier gap, paired with the prime that
// starts it. The result is a slice of [gap, startingPrime] pairs.
func MaximalPrimeGaps(x int) [][2]int64 {
	primes := PrimesUpTo(x)
	var out [][2]int64
	var record int64
	for i := 1; i < len(primes); i++ {
		g := primes[i] - primes[i-1]
		if g > record {
			record = g
			out = append(out, [2]int64{g, primes[i-1]})
		}
	}
	return out
}

// MaxPrimeGap returns the largest gap between consecutive primes up to x, along
// with the prime that begins it.
func MaxPrimeGap(x int) (gap int64, startingPrime int64) {
	primes := PrimesUpTo(x)
	for i := 1; i < len(primes); i++ {
		g := primes[i] - primes[i-1]
		if g > gap {
			gap = g
			startingPrime = primes[i-1]
		}
	}
	return gap, startingPrime
}

// AverageGap returns the mean gap between consecutive primes up to x, which by
// the prime number theorem is asymptotic to ln x.
func AverageGap(x int) float64 {
	gaps := PrimeGaps(x)
	if len(gaps) == 0 {
		return 0
	}
	var sum int64
	for _, g := range gaps {
		sum += g
	}
	return float64(sum) / float64(len(gaps))
}

// TwinPrimesUpTo returns all twin-prime pairs (p, p+2) with p+2 ≤ x, as a slice
// of [p, p+2] pairs.
func TwinPrimesUpTo(x int) [][2]int64 {
	return primePairsUpTo(x, 2)
}

// CousinPrimesUpTo returns all cousin-prime pairs (p, p+4) with p+4 ≤ x.
func CousinPrimesUpTo(x int) [][2]int64 {
	return primePairsUpTo(x, 4)
}

// SexyPrimesUpTo returns all sexy-prime pairs (p, p+6) with p+6 ≤ x.
func SexyPrimesUpTo(x int) [][2]int64 {
	return primePairsUpTo(x, 6)
}

// primePairsUpTo returns pairs (p, p+d) of primes with p+d ≤ x.
func primePairsUpTo(x int, d int) [][2]int64 {
	if x < 2+d {
		return [][2]int64{}
	}
	s := Sieve(x)
	var out [][2]int64
	for p := 2; p+d <= x; p++ {
		if s[p] && s[p+d] {
			out = append(out, [2]int64{int64(p), int64(p + d)})
		}
	}
	return out
}

// TwinPrimeCount returns the number of twin-prime pairs (p, p+2) with p+2 ≤ x.
func TwinPrimeCount(x int) int64 {
	return int64(len(TwinPrimesUpTo(x)))
}

// PrimeTripletsUpTo returns prime triplets (p, p+2, p+6) or (p, p+4, p+6) — the
// two admissible constellations of diameter 6 — with the largest member ≤ x.
func PrimeTripletsUpTo(x int) [][3]int64 {
	if x < 8 {
		return [][3]int64{}
	}
	s := Sieve(x)
	var out [][3]int64
	for p := 2; p+6 <= x; p++ {
		if !s[p] {
			continue
		}
		if s[p+2] && s[p+6] {
			out = append(out, [3]int64{int64(p), int64(p + 2), int64(p + 6)})
		} else if s[p+4] && s[p+6] {
			out = append(out, [3]int64{int64(p), int64(p + 4), int64(p + 6)})
		}
	}
	return out
}

// MeritOfGap returns the merit of a prime gap g following prime p, defined as
// g/ln p; merits above about 30 correspond to exceptionally large gaps.
func MeritOfGap(g int64, p int64) float64 {
	return float64(g) / math.Log(float64(p))
}

// GoldbachPartitions returns the number of ways to write the even integer n as
// an unordered sum of two primes p + q with p ≤ q (the Goldbach partition
// count, Goldbach's comet). n must be even and >= 4.
func GoldbachPartitions(n int) int64 {
	if n < 4 || n%2 != 0 {
		panic("analyticnt: GoldbachPartitions requires an even n >= 4")
	}
	s := Sieve(n)
	var count int64
	for p := 2; p <= n/2; p++ {
		if s[p] && s[n-p] {
			count++
		}
	}
	return count
}

// GoldbachPair returns a pair of primes summing to the even integer n, choosing
// the pair with the smallest first prime. It returns (0,0,false) if none is
// found (which cannot happen for even n >= 4 within tested ranges).
func GoldbachPair(n int) (p, q int64, ok bool) {
	if n < 4 || n%2 != 0 {
		panic("analyticnt: GoldbachPair requires an even n >= 4")
	}
	s := Sieve(n)
	for a := 2; a <= n/2; a++ {
		if s[a] && s[n-a] {
			return int64(a), int64(n - a), true
		}
	}
	return 0, 0, false
}

// RandomPrimeBelow returns a uniformly random prime strictly less than n using a
// deterministic math/rand source seeded by seed. It returns 0 if there is no
// prime below n.
func RandomPrimeBelow(n int, seed int64) int64 {
	primes := PrimesUpTo(n - 1)
	if len(primes) == 0 {
		return 0
	}
	r := rand.New(rand.NewSource(seed))
	return primes[r.Intn(len(primes))]
}

// PrimeGapHistogram returns a map from gap size to its frequency among
// consecutive primes up to x.
func PrimeGapHistogram(x int) map[int64]int {
	gaps := PrimeGaps(x)
	h := make(map[int64]int)
	for _, g := range gaps {
		h[g]++
	}
	return h
}

// HardyLittlewoodTwinEstimate returns the Hardy–Littlewood asymptotic estimate
// 2·C₂·x/(ln x)² for the number of twin primes up to x.
func HardyLittlewoodTwinEstimate(x float64) float64 {
	if x < 3 {
		return 0
	}
	l := math.Log(x)
	return 2 * TwinPrimeConstant * x / (l * l)
}
