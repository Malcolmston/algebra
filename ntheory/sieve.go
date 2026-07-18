package ntheory

import (
	"math"
	"sort"
	"sync"
)

// This file provides a memory-lean, uint64 prime-enumeration path that
// complements the int64 [PrimesUpTo]/[PrimePi] helpers. Where a full boolean
// sieve needs O(n) memory spanning [0, n], the routines here sieve fixed-size
// windows one at a time, so their working set depends on the window width
// rather than on the largest value examined. That makes it practical to
// enumerate, count, or index primes in high ranges that a single n+1 array
// could not hold.
//
// Performance technique: a mod-30 (2, 3, 5) wheel. Every prime other than
// 2, 3 and 5 is congruent to one of the eight residues {1, 7, 11, 13, 17, 19,
// 23, 29} modulo 30, so each segment stores only 8 flags per 30 integers
// (8/30 of the naive size) and stepping across a prime's multiples skips the
// two-thirds that are divisible by 2, 3 or 5. The base primes up to sqrt(hi)
// are computed once and memoized in a package-level cache guarded by a mutex,
// then reused across every segment and every call.

// ntheorySieveBlocks is the number of mod-30 wheel blocks in one segment. Each
// block covers 30 consecutive integers with 8 flag bytes, so a segment's
// reusable window is ntheorySieveBlocks*8 == 1<<16 bytes, a cache-friendly
// fixed size independent of how large the sieved values grow.
const ntheorySieveBlocks = 1 << 13

// ntheorySieveSpan is the count of integers covered by one segment.
const ntheorySieveSpan uint64 = ntheorySieveBlocks * 30

// ntheorySieveWheelResidues lists, in ascending order, the residues modulo 30
// that are coprime to 30. Every prime greater than 5 reduces to one of these.
var ntheorySieveWheelResidues = [8]uint64{1, 7, 11, 13, 17, 19, 23, 29}

// ntheorySieveWheelGaps[i] is the distance from ntheorySieveWheelResidues[i] to
// the next integer coprime to 30 (wrapping past 30). Advancing a multiple by
// prime*gap jumps straight to the prime's next multiple that is itself coprime
// to 30, skipping the ones divisible by 2, 3 or 5.
var ntheorySieveWheelGaps = [8]uint64{6, 4, 2, 4, 2, 4, 6, 2}

// ntheorySieveWheelPos maps a residue modulo 30 to its index within
// ntheorySieveWheelResidues, or 255 when the residue is not coprime to 30.
var ntheorySieveWheelPos = func() [30]uint8 {
	var pos [30]uint8
	for i := range pos {
		pos[i] = 255
	}
	for i, r := range ntheorySieveWheelResidues {
		pos[r] = uint8(i)
	}
	return pos
}()

// ntheorySieveSmallPrimes are the primes excluded by the mod-30 wheel; they are
// emitted or counted directly whenever they fall inside the requested range.
var ntheorySieveSmallPrimes = [3]uint64{2, 3, 5}

// Package-level memoized cache of base primes. It only ever grows by full
// recomputation (a fresh slice), never by in-place mutation, so a slice handed
// back to a caller stays valid even if a later call replaces the cache.
var (
	ntheorySieveBaseMu    sync.Mutex
	ntheorySieveBaseCache []uint64
	ntheorySieveBaseLimit uint64
	ntheorySieveBaseInit  bool
)

// ntheorySieveBasePrimes returns every prime p with p <= limit, computing them
// with a plain sieve of Eratosthenes the first time a given limit is needed and
// caching the result for reuse across segments and calls. The returned slice
// must not be modified by the caller.
func ntheorySieveBasePrimes(limit uint64) []uint64 {
	ntheorySieveBaseMu.Lock()
	defer ntheorySieveBaseMu.Unlock()
	if !ntheorySieveBaseInit || limit > ntheorySieveBaseLimit {
		ntheorySieveBaseCache = ntheorySieveSimple(limit)
		ntheorySieveBaseLimit = limit
		ntheorySieveBaseInit = true
		return ntheorySieveBaseCache
	}
	if limit == ntheorySieveBaseLimit {
		return ntheorySieveBaseCache
	}
	// The cache reaches beyond limit; return the prefix of primes <= limit.
	i := sort.Search(len(ntheorySieveBaseCache), func(k int) bool {
		return ntheorySieveBaseCache[k] > limit
	})
	return ntheorySieveBaseCache[:i]
}

// ntheorySieveSimple returns all primes p with p <= limit using a direct
// boolean sieve of Eratosthenes. It backs the base-prime cache and is only ever
// called with limit around sqrt(hi), so its O(limit) memory stays small
// relative to the ranges the segmented sieve enumerates. It returns nil for
// limit < 2.
func ntheorySieveSimple(limit uint64) []uint64 {
	if limit < 2 {
		return nil
	}
	composite := make([]bool, limit+1)
	var primes []uint64
	for i := uint64(2); i <= limit; i++ {
		if composite[i] {
			continue
		}
		primes = append(primes, i)
		for j := i * i; j <= limit; j += i {
			composite[j] = true
		}
	}
	return primes
}

// ntheorySieveIsqrt returns floor(sqrt(n)) exactly, using math.Sqrt for a seed
// and integer comparisons (division rather than multiplication) to correct any
// floating-point rounding without risking overflow near the top of uint64.
func ntheorySieveIsqrt(n uint64) uint64 {
	if n < 2 {
		return n
	}
	x := uint64(math.Sqrt(float64(n)))
	for x > 1 && x > n/x {
		x--
	}
	for (x + 1) <= n/(x+1) {
		x++
	}
	return x
}

// ntheorySieveSegment resets seg and marks, as composite, every wheel candidate
// in the window [segStart, segStart+ntheorySieveSpan) that is a multiple of a
// base prime. segStart must be a multiple of 30 and seg must have length
// ntheorySieveBlocks*8. base holds the primes up to sqrt of the window's top;
// primes below 7 are skipped because the wheel already excludes their
// multiples.
func ntheorySieveSegment(seg []bool, segStart uint64, base []uint64) {
	for i := range seg {
		seg[i] = false
	}
	segEnd := segStart + ntheorySieveSpan // exclusive
	if segEnd < segStart {                // overflow at the very top of uint64
		segEnd = ^uint64(0)
	}
	for _, p := range base {
		if p < 7 {
			continue
		}
		pp := p * p
		if pp >= segEnd {
			// Larger base primes only have larger squares, so none can
			// contribute a multiple within this window.
			break
		}
		// L is the first value we might mark: the larger of the window start
		// and p*p (smaller multiples of p were crossed off by smaller primes).
		L := segStart
		if pp > L {
			L = pp
		}
		t := (L + p - 1) / p // smallest t with p*t >= L
		// Advance t to the next integer coprime to 30 so p*t is a wheel
		// candidate; at most a handful of steps.
		var tIdx uint8
		for {
			tIdx = ntheorySieveWheelPos[t%30]
			if tIdx != 255 {
				break
			}
			t++
		}
		m := p * t
		for m < segEnd {
			seg[(m-segStart)/30*8+uint64(ntheorySieveWheelPos[m%30])] = true
			step := p * ntheorySieveWheelGaps[tIdx]
			tIdx = (tIdx + 1) & 7
			next := m + step
			if next < m { // overflow guard
				break
			}
			m = next
		}
	}
}

// ntheorySieveCollect walks a sieved segment in ascending order and invokes emit
// for every wheel candidate n with lo <= n <= hi that survived as prime.
// segStart must be the window's start (a multiple of 30) and seg the window
// just processed by ntheorySieveSegment.
func ntheorySieveCollect(seg []bool, segStart, lo, hi uint64, emit func(uint64)) {
	for b := 0; b < ntheorySieveBlocks; b++ {
		blockBase := segStart + uint64(b)*30
		if blockBase > hi {
			return
		}
		row := b * 8
		for r := 0; r < 8; r++ {
			if seg[row+r] {
				continue
			}
			n := blockBase + ntheorySieveWheelResidues[r]
			// n < 7 filters out 1 (coprime to 30 but not prime); 2, 3 and 5
			// are handled separately by the small-prime path.
			if n < 7 || n < lo || n > hi {
				continue
			}
			emit(n)
		}
	}
}

// SegmentedSieve returns all primes p with lo <= p <= hi in ascending order.
//
// It uses a segmented sieve of Eratosthenes: the base primes up to sqrt(hi) are
// computed once and memoized in a package-level cache, then each fixed-size
// window of the interval is sieved into a reusable mod-30 wheel array holding
// only 8 flags per 30 integers. Because it never allocates an array spanning
// [0, hi], it handles ranges whose upper bound is far beyond what a single n+1
// boolean sieve could fit in memory. lo and hi are treated as non-negative;
// the result is deterministic. It returns nil when the interval contains no
// primes, including when hi < 2 or lo > hi.
func SegmentedSieve(lo, hi uint64) []uint64 {
	if hi < 2 || lo > hi {
		return nil
	}
	effLo := lo
	if effLo < 2 {
		effLo = 2
	}
	var primes []uint64
	for _, p := range ntheorySieveSmallPrimes {
		if p >= effLo && p <= hi {
			primes = append(primes, p)
		}
	}
	if hi < 7 {
		return primes
	}
	base := ntheorySieveBasePrimes(ntheorySieveIsqrt(hi))
	seg := make([]bool, ntheorySieveBlocks*8)
	segStart := effLo / 30 * 30
	for segStart <= hi {
		ntheorySieveSegment(seg, segStart, base)
		ntheorySieveCollect(seg, segStart, effLo, hi, func(n uint64) {
			primes = append(primes, n)
		})
		next := segStart + ntheorySieveSpan
		if next < segStart { // reached the top of uint64
			break
		}
		segStart = next
	}
	return primes
}

// PrimesInRange returns all primes in the half-open interval [lo, hi) in
// ascending order. It is a convenience wrapper over [SegmentedSieve] with an
// exclusive upper bound. It returns nil when the interval contains no primes.
func PrimesInRange(lo, hi uint64) []uint64 {
	if hi == 0 {
		return nil
	}
	return SegmentedSieve(lo, hi-1)
}

// PrimePiRange counts the primes p with lo <= p <= hi without materializing
// them, sieving one segment at a time and tallying survivors. Its memory use is
// bounded by a single window regardless of hi, making it suitable for a
// prime-counting function over huge ranges. It returns 0 when hi < 2 or
// lo > hi.
func PrimePiRange(lo, hi uint64) uint64 {
	if hi < 2 || lo > hi {
		return 0
	}
	effLo := lo
	if effLo < 2 {
		effLo = 2
	}
	var count uint64
	for _, p := range ntheorySieveSmallPrimes {
		if p >= effLo && p <= hi {
			count++
		}
	}
	if hi < 7 {
		return count
	}
	base := ntheorySieveBasePrimes(ntheorySieveIsqrt(hi))
	seg := make([]bool, ntheorySieveBlocks*8)
	segStart := effLo / 30 * 30
	for segStart <= hi {
		ntheorySieveSegment(seg, segStart, base)
		ntheorySieveCollect(seg, segStart, effLo, hi, func(uint64) {
			count++
		})
		next := segStart + ntheorySieveSpan
		if next < segStart {
			break
		}
		segStart = next
	}
	return count
}

// ntheorySieveNthUpperBound returns a value guaranteed to be at least the n-th
// prime. For n >= 6 it uses the bound p_n < n(ln n + ln ln n); the caller
// handles smaller n with an exact table.
func ntheorySieveNthUpperBound(n uint64) uint64 {
	fn := float64(n)
	ln := math.Log(fn)
	lnln := math.Log(ln)
	return uint64(fn*(ln+lnln)) + 3
}

// NthPrime returns the n-th prime, 1-indexed, so NthPrime(1) == 2,
// NthPrime(2) == 3 and so on. It estimates an upper bound with
// n(ln n + ln ln n), warms the base-prime cache up to the square root of that
// bound, then streams a segmented sieve until the n-th prime appears. It panics
// if n < 1.
func NthPrime(n uint64) uint64 {
	if n < 1 {
		panic("ntheory: NthPrime requires n >= 1")
	}
	// The estimate below needs n >= 6 to be a valid upper bound; small n use an
	// exact table.
	table := [...]uint64{2, 3, 5, 7, 11}
	if n <= uint64(len(table)) {
		return table[n-1]
	}
	bound := ntheorySieveNthUpperBound(n)
	// Warm the memoized base primes to the size this sieve will need.
	_ = ntheorySieveBasePrimes(ntheorySieveIsqrt(bound))
	s := NewPrimeSieve()
	var count uint64
	for {
		p := s.Next()
		count++
		if count == n {
			return p
		}
		if p > bound { // safety net; the bound is a proven over-estimate
			panic("ntheory: NthPrime exceeded its prime upper bound")
		}
	}
}

// PrimeSieve is a stateful generator that yields the primes in ascending order,
// one per call to [PrimeSieve.Next], without ever holding them all in memory. It
// advances through fixed-size segments lazily, reusing a single mod-30 wheel
// window and the shared base-prime cache. The zero value is not usable; obtain
// one from [NewPrimeSieve]. A PrimeSieve is not safe for concurrent use.
type PrimeSieve struct {
	seg      []bool   // reusable wheel window for the current segment
	segStart uint64   // start of the next segment to sieve (multiple of 30)
	buf      []uint64 // primes found in the current segment, awaiting emission
	bufIdx   int      // cursor into buf
	smallIdx int      // cursor into the small primes 2, 3, 5
}

// NewPrimeSieve returns a PrimeSieve positioned before the first prime, so the
// first call to Next returns 2.
func NewPrimeSieve() *PrimeSieve {
	return &PrimeSieve{
		seg: make([]bool, ntheorySieveBlocks*8),
	}
}

// Next returns the next prime in ascending order, advancing to a further
// segment only when the current one is exhausted. Successive calls yield
// 2, 3, 5, 7, 11, ... deterministically.
func (s *PrimeSieve) Next() uint64 {
	if s.smallIdx < len(ntheorySieveSmallPrimes) {
		p := ntheorySieveSmallPrimes[s.smallIdx]
		s.smallIdx++
		return p
	}
	for s.bufIdx >= len(s.buf) {
		s.fill()
	}
	p := s.buf[s.bufIdx]
	s.bufIdx++
	return p
}

// fill sieves the next segment and loads its primes into the buffer.
func (s *PrimeSieve) fill() {
	segEnd := s.segStart + ntheorySieveSpan
	if segEnd < s.segStart {
		segEnd = ^uint64(0)
	}
	base := ntheorySieveBasePrimes(ntheorySieveIsqrt(segEnd))
	ntheorySieveSegment(s.seg, s.segStart, base)
	s.buf = s.buf[:0]
	segStart := s.segStart
	ntheorySieveCollect(s.seg, segStart, 0, segEnd-1, func(n uint64) {
		s.buf = append(s.buf, n)
	})
	s.bufIdx = 0
	s.segStart = segEnd
}
