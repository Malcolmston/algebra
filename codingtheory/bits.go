package codingtheory

import "errors"

// ErrLength is returned when an input slice has the wrong length for a code.
var ErrLength = errors.New("codingtheory: wrong length for code")

// ErrBit is returned when a slice expected to contain bits holds a value other
// than 0 or 1.
var ErrBit = errors.New("codingtheory: slice entry is not a bit")

// validBits reports whether every entry of v is 0 or 1.
func validBits(v []int) bool {
	for _, x := range v {
		if x != 0 && x != 1 {
			return false
		}
	}
	return true
}

// HammingWeight returns the number of non-zero entries of v.
func HammingWeight(v []int) int {
	w := 0
	for _, x := range v {
		if x != 0 {
			w++
		}
	}
	return w
}

// HammingDistance returns the number of positions in which a and b differ. It
// panics if the slices have different lengths.
func HammingDistance(a, b []int) int {
	if len(a) != len(b) {
		panic("codingtheory: HammingDistance on unequal lengths")
	}
	d := 0
	for i := range a {
		if a[i] != b[i] {
			d++
		}
	}
	return d
}

// HammingWeightUint returns the population count of an unsigned integer's bits.
func HammingWeightUint(x uint64) int {
	c := 0
	for x != 0 {
		x &= x - 1
		c++
	}
	return c
}

// HammingDistanceUint returns the number of differing bits between a and b.
func HammingDistanceUint(a, b uint64) int { return HammingWeightUint(a ^ b) }

// XORVectors returns the entrywise exclusive-or of two equal-length bit
// vectors. It panics on a length mismatch.
func XORVectors(a, b []int) []int {
	if len(a) != len(b) {
		panic("codingtheory: XORVectors on unequal lengths")
	}
	out := make([]int, len(a))
	for i := range a {
		out[i] = a[i] ^ b[i]
	}
	return out
}

// DotGF2 returns the GF(2) inner product (parity of the pairwise ands) of two
// equal-length bit vectors. It panics on a length mismatch.
func DotGF2(a, b []int) int {
	if len(a) != len(b) {
		panic("codingtheory: DotGF2 on unequal lengths")
	}
	s := 0
	for i := range a {
		s ^= (a[i] & 1) & (b[i] & 1)
	}
	return s
}

// MatVecGF2 multiplies the GF(2) matrix m (rows of equal length) by the column
// vector v, returning m*v over GF(2). It panics on a dimension mismatch.
func MatVecGF2(m [][]int, v []int) []int {
	out := make([]int, len(m))
	for i, row := range m {
		if len(row) != len(v) {
			panic("codingtheory: MatVecGF2 dimension mismatch")
		}
		out[i] = DotGF2(row, v)
	}
	return out
}

// VecMatGF2 multiplies the row vector v by the GF(2) matrix m, returning v*m
// over GF(2). It panics on a dimension mismatch.
func VecMatGF2(v []int, m [][]int) []int {
	if len(m) == 0 {
		return nil
	}
	cols := len(m[0])
	out := make([]int, cols)
	for i, row := range m {
		if len(row) != cols {
			panic("codingtheory: VecMatGF2 ragged matrix")
		}
		if i >= len(v) {
			panic("codingtheory: VecMatGF2 dimension mismatch")
		}
		if v[i]&1 == 0 {
			continue
		}
		for j := 0; j < cols; j++ {
			out[j] ^= row[j] & 1
		}
	}
	if len(v) != len(m) {
		panic("codingtheory: VecMatGF2 dimension mismatch")
	}
	return out
}

// MatMulGF2 returns the product a*b of two GF(2) matrices. It panics on a
// dimension mismatch.
func MatMulGF2(a, b [][]int) [][]int {
	if len(a) == 0 || len(b) == 0 {
		return nil
	}
	inner := len(b)
	cols := len(b[0])
	out := make([][]int, len(a))
	for i := range a {
		if len(a[i]) != inner {
			panic("codingtheory: MatMulGF2 dimension mismatch")
		}
		out[i] = make([]int, cols)
		for k := 0; k < inner; k++ {
			if a[i][k]&1 == 0 {
				continue
			}
			for j := 0; j < cols; j++ {
				out[i][j] ^= b[k][j] & 1
			}
		}
	}
	return out
}

// TransposeGF2 returns the transpose of a rectangular GF(2) matrix.
func TransposeGF2(m [][]int) [][]int {
	if len(m) == 0 {
		return nil
	}
	rows := len(m)
	cols := len(m[0])
	out := make([][]int, cols)
	for j := 0; j < cols; j++ {
		out[j] = make([]int, rows)
		for i := 0; i < rows; i++ {
			out[j][i] = m[i][j]
		}
	}
	return out
}

// IdentityGF2 returns the n-by-n identity matrix over GF(2).
func IdentityGF2(n int) [][]int {
	out := make([][]int, n)
	for i := 0; i < n; i++ {
		out[i] = make([]int, n)
		out[i][i] = 1
	}
	return out
}

// RankGF2 returns the rank of a GF(2) matrix computed by Gaussian elimination.
// The input matrix is not modified.
func RankGF2(m [][]int) int {
	if len(m) == 0 {
		return 0
	}
	rows := len(m)
	cols := len(m[0])
	a := make([][]int, rows)
	for i := range m {
		a[i] = append([]int(nil), m[i]...)
	}
	rank := 0
	for col := 0; col < cols && rank < rows; col++ {
		pivot := -1
		for r := rank; r < rows; r++ {
			if a[r][col] != 0 {
				pivot = r
				break
			}
		}
		if pivot == -1 {
			continue
		}
		a[rank], a[pivot] = a[pivot], a[rank]
		for r := 0; r < rows; r++ {
			if r != rank && a[r][col] != 0 {
				for c := 0; c < cols; c++ {
					a[r][c] ^= a[rank][c]
				}
			}
		}
		rank++
	}
	return rank
}

// BitsToUint packs a bit slice (most-significant bit first) into a uint64. It
// panics if the slice is longer than 64 entries.
func BitsToUint(bits []int) uint64 {
	if len(bits) > 64 {
		panic("codingtheory: BitsToUint slice too long")
	}
	var x uint64
	for _, b := range bits {
		x = (x << 1) | uint64(b&1)
	}
	return x
}

// UintToBits unpacks the low n bits of x into a bit slice, most-significant bit
// first.
func UintToBits(x uint64, n int) []int {
	out := make([]int, n)
	for i := 0; i < n; i++ {
		out[n-1-i] = int((x >> uint(i)) & 1)
	}
	return out
}
