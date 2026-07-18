package logic

import "testing"

// benchMinterms is a fixed six-variable minterm set exercising the
// Quine-McCluskey engine (the heaviest routine in the package).
var benchMinterms = []int{
	0, 1, 2, 3, 5, 7, 8, 10, 13, 15,
	18, 20, 21, 23, 26, 29, 31, 34, 37, 40,
	42, 45, 47, 50, 53, 55, 58, 61, 63,
}

func BenchmarkQuineMcCluskey(b *testing.B) {
	for i := 0; i < b.N; i++ {
		_ = QuineMcCluskey(benchMinterms, nil, 6)
	}
}
