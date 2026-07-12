package ntheory

import (
	"math/big"
	"reflect"
	"testing"
)

func TestGCD(t *testing.T) {
	cases := []struct{ a, b, want int64 }{
		{48, 36, 12},
		{0, 0, 0},
		{7, 0, 7},
		{0, 5, 5},
		{-48, 36, 12},
		{270, 192, 6},
		{17, 5, 1},
	}
	for _, c := range cases {
		if got := GCD(c.a, c.b); got != c.want {
			t.Errorf("GCD(%d,%d)=%d want %d", c.a, c.b, got, c.want)
		}
	}
}

func TestLCM(t *testing.T) {
	cases := []struct{ a, b, want int64 }{
		{4, 6, 12},
		{21, 6, 42},
		{0, 5, 0},
		{-4, 6, 12},
	}
	for _, c := range cases {
		if got := LCM(c.a, c.b); got != c.want {
			t.Errorf("LCM(%d,%d)=%d want %d", c.a, c.b, got, c.want)
		}
	}
}

func TestExtendedGCD(t *testing.T) {
	cases := []struct{ a, b int64 }{{48, 36}, {240, 46}, {-17, 5}, {7, 0}}
	for _, c := range cases {
		g, x, y := ExtendedGCD(c.a, c.b)
		if g != GCD(c.a, c.b) {
			t.Errorf("ExtendedGCD(%d,%d) g=%d want %d", c.a, c.b, g, GCD(c.a, c.b))
		}
		if c.a*x+c.b*y != g {
			t.Errorf("ExtendedGCD(%d,%d): %d*%d+%d*%d != %d", c.a, c.b, c.a, x, c.b, y, g)
		}
	}
}

func TestDivisors(t *testing.T) {
	got := Divisors(360)
	want := []int64{1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 18, 20, 24, 30, 36, 40, 45, 60, 72, 90, 120, 180, 360}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("Divisors(360)=%v want %v", got, want)
	}
	if d := Divisors(1); !reflect.DeepEqual(d, []int64{1}) {
		t.Errorf("Divisors(1)=%v want [1]", d)
	}
	if d := Divisors(0); d != nil {
		t.Errorf("Divisors(0)=%v want nil", d)
	}
}

func TestSumAndCountDivisors(t *testing.T) {
	if got := SumDivisors(28); got != 56 {
		t.Errorf("SumDivisors(28)=%d want 56", got)
	}
	if got := SumDivisors(12); got != 28 {
		t.Errorf("SumDivisors(12)=%d want 28", got)
	}
	if got := CountDivisors(360); got != 24 {
		t.Errorf("CountDivisors(360)=%d want 24", got)
	}
	if got := CountDivisors(1); got != 1 {
		t.Errorf("CountDivisors(1)=%d want 1", got)
	}
}

func TestIsPerfect(t *testing.T) {
	perfect := map[int64]bool{6: true, 28: true, 496: true, 8128: true}
	for n := int64(1); n <= 8200; n++ {
		if IsPerfect(n) != perfect[n] {
			t.Errorf("IsPerfect(%d)=%v want %v", n, IsPerfect(n), perfect[n])
		}
	}
}

func TestIsPrime(t *testing.T) {
	primesTo30 := map[int64]bool{2: true, 3: true, 5: true, 7: true, 11: true, 13: true, 17: true, 19: true, 23: true, 29: true}
	for n := int64(0); n <= 30; n++ {
		if IsPrime(n) != primesTo30[n] {
			t.Errorf("IsPrime(%d)=%v want %v", n, IsPrime(n), primesTo30[n])
		}
	}
	// Carmichael numbers must be reported composite.
	for _, c := range []int64{561, 1105, 1729, 2465, 2821, 6601, 8911, 10585, 41041} {
		if IsPrime(c) {
			t.Errorf("IsPrime(%d)=true; Carmichael number must be composite", c)
		}
	}
	// Known large primes.
	for _, p := range []int64{7919, 104729, 1000003, 2147483647, 999999999989} {
		if !IsPrime(p) {
			t.Errorf("IsPrime(%d)=false want true", p)
		}
	}
	// Known large composite.
	for _, c := range []int64{7919 * 7919, 1000003 * 1000033} {
		if IsPrime(c) {
			t.Errorf("IsPrime(%d)=true want false", c)
		}
	}
	if IsPrime(-7) {
		t.Errorf("IsPrime(-7)=true want false")
	}
}

func TestIsPrimeBig(t *testing.T) {
	if !IsPrimeBig(big.NewInt(1000003)) {
		t.Errorf("IsPrimeBig(1000003)=false want true")
	}
	// A 100-digit prime.
	p, _ := new(big.Int).SetString("2074722246773485207821695222107608587480996474721117292752992589912196684750549658310084416732550077", 10)
	if !IsPrimeBig(p) {
		t.Errorf("IsPrimeBig(100-digit prime)=false want true")
	}
	if IsPrimeBig(big.NewInt(561)) {
		t.Errorf("IsPrimeBig(561)=true want false")
	}
}

func TestNextPrime(t *testing.T) {
	cases := []struct{ n, want int64 }{{0, 2}, {2, 3}, {13, 17}, {17, 19}, {100, 101}, {-5, 2}}
	for _, c := range cases {
		if got := NextPrime(c.n); got != c.want {
			t.Errorf("NextPrime(%d)=%d want %d", c.n, got, c.want)
		}
	}
}

func TestPrimesUpTo(t *testing.T) {
	got := PrimesUpTo(30)
	want := []int64{2, 3, 5, 7, 11, 13, 17, 19, 23, 29}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("PrimesUpTo(30)=%v want %v", got, want)
	}
	if PrimesUpTo(1) != nil {
		t.Errorf("PrimesUpTo(1) want nil")
	}
}

func TestPrimePi(t *testing.T) {
	cases := []struct{ n, want int64 }{{10, 4}, {100, 25}, {1000, 168}, {1, 0}}
	for _, c := range cases {
		if got := PrimePi(c.n); got != c.want {
			t.Errorf("PrimePi(%d)=%d want %d", c.n, got, c.want)
		}
	}
}

func TestFactorize(t *testing.T) {
	got := Factorize(360)
	want := map[int64]int{2: 3, 3: 2, 5: 1}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("Factorize(360)=%v want %v", got, want)
	}
	if got := Factorize(97); !reflect.DeepEqual(got, map[int64]int{97: 1}) {
		t.Errorf("Factorize(97)=%v want {97:1}", got)
	}
	if got := Factorize(1); len(got) != 0 {
		t.Errorf("Factorize(1)=%v want empty", got)
	}
}

func TestFactorList(t *testing.T) {
	got := FactorList(360)
	want := []PrimePower{{2, 3}, {3, 2}, {5, 1}}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("FactorList(360)=%v want %v", got, want)
	}
	if FactorList(1) != nil {
		t.Errorf("FactorList(1) want nil")
	}
}

func TestEulerPhi(t *testing.T) {
	cases := []struct{ n, want int64 }{{1, 1}, {10, 4}, {36, 12}, {97, 96}, {1, 1}}
	for _, c := range cases {
		if got := EulerPhi(c.n); got != c.want {
			t.Errorf("EulerPhi(%d)=%d want %d", c.n, got, c.want)
		}
	}
}

func TestMobiusMu(t *testing.T) {
	cases := []struct {
		n    int64
		want int
	}{{1, 1}, {2, -1}, {6, 1}, {4, 0}, {30, -1}, {12, 0}, {97, -1}}
	for _, c := range cases {
		if got := MobiusMu(c.n); got != c.want {
			t.Errorf("MobiusMu(%d)=%d want %d", c.n, got, c.want)
		}
	}
}

func TestRadical(t *testing.T) {
	cases := []struct{ n, want int64 }{{1, 1}, {360, 30}, {12, 6}, {97, 97}, {8, 2}}
	for _, c := range cases {
		if got := Radical(c.n); got != c.want {
			t.Errorf("Radical(%d)=%d want %d", c.n, got, c.want)
		}
	}
}

func TestModPow(t *testing.T) {
	cases := []struct{ b, e, m, want int64 }{
		{2, 10, 1000, 24},
		{3, 4, 5, 1},
		{7, 0, 13, 1},
		{2, 100, 1000000007, 976371285},
		{10, 3, 1, 0},
	}
	for _, c := range cases {
		if got := ModPow(c.b, c.e, c.m); got != c.want {
			t.Errorf("ModPow(%d,%d,%d)=%d want %d", c.b, c.e, c.m, got, c.want)
		}
	}
	// Negative exponent uses the modular inverse: 3^-1 mod 11 = 4.
	if got := ModPow(3, -1, 11); got != 4 {
		t.Errorf("ModPow(3,-1,11)=%d want 4", got)
	}
}

func TestModInverse(t *testing.T) {
	if inv, ok := ModInverse(3, 11); !ok || inv != 4 {
		t.Errorf("ModInverse(3,11)=%d,%v want 4,true", inv, ok)
	}
	if inv, ok := ModInverse(10, 17); !ok || inv != 12 {
		t.Errorf("ModInverse(10,17)=%d,%v want 12,true", inv, ok)
	}
	if _, ok := ModInverse(2, 4); ok {
		t.Errorf("ModInverse(2,4) want ok=false")
	}
	if inv, ok := ModInverse(-8, 11); !ok || inv != 4 {
		t.Errorf("ModInverse(-8,11)=%d,%v want 4,true", inv, ok)
	}
}

func TestCRT(t *testing.T) {
	// x ≡ 2 (mod 3), x ≡ 3 (mod 5), x ≡ 2 (mod 7) -> 23 mod 105.
	if x, m, ok := CRT([]int64{2, 3, 2}, []int64{3, 5, 7}); !ok || x != 23 || m != 105 {
		t.Errorf("CRT=%d,%d,%v want 23,105,true", x, m, ok)
	}
	// Non-coprime consistent: x ≡ 1 (mod 4), x ≡ 3 (mod 6) -> 9 mod 12.
	if x, m, ok := CRT([]int64{1, 3}, []int64{4, 6}); !ok || x != 9 || m != 12 {
		t.Errorf("CRT non-coprime=%d,%d,%v want 9,12,true", x, m, ok)
	}
	// Inconsistent: x ≡ 0 (mod 2), x ≡ 1 (mod 4).
	if _, _, ok := CRT([]int64{0, 1}, []int64{2, 4}); ok {
		t.Errorf("CRT inconsistent want ok=false")
	}
}

func TestLegendreSymbol(t *testing.T) {
	// Residues mod 7 are {1,2,4}; non-residues {3,5,6}.
	want := map[int64]int{0: 0, 1: 1, 2: 1, 3: -1, 4: 1, 5: -1, 6: -1}
	for a, w := range want {
		if got := LegendreSymbol(a, 7); got != w {
			t.Errorf("LegendreSymbol(%d,7)=%d want %d", a, got, w)
		}
	}
}

func TestJacobiSymbol(t *testing.T) {
	cases := []struct {
		a, n int64
		want int
	}{
		{1, 9, 1}, {2, 9, 1}, {3, 9, 0}, {5, 21, 1}, {2, 15, 1}, {1001, 9907, -1},
	}
	for _, c := range cases {
		if got := JacobiSymbol(c.a, c.n); got != c.want {
			t.Errorf("JacobiSymbol(%d,%d)=%d want %d", c.a, c.n, got, c.want)
		}
	}
	// For a prime modulus, Jacobi and Legendre agree.
	for a := int64(0); a < 11; a++ {
		if JacobiSymbol(a, 11) != LegendreSymbol(a, 11) {
			t.Errorf("Jacobi/Legendre disagree at a=%d, p=11", a)
		}
	}
}

func TestIsQuadraticResidue(t *testing.T) {
	residues := map[int64]bool{0: true, 1: true, 2: true, 3: false, 4: true, 5: false, 6: false}
	for a, w := range residues {
		if got := IsQuadraticResidue(a, 7); got != w {
			t.Errorf("IsQuadraticResidue(%d,7)=%v want %v", a, got, w)
		}
	}
	if !IsQuadraticResidue(1, 2) {
		t.Errorf("IsQuadraticResidue(1,2) want true")
	}
}

func TestDiscreteLog(t *testing.T) {
	// 3 is a primitive root mod 7; 3^x ≡ 5 gives x = 5.
	if x, ok := DiscreteLog(3, 5, 7); !ok || ModPow(3, x, 7) != 5 {
		t.Errorf("DiscreteLog(3,5,7)=%d,%v invalid", x, ok)
	}
	// 2^x ≡ 22 (mod 29) -> x = 5.
	if x, ok := DiscreteLog(2, 22, 29); !ok || ModPow(2, x, 29) != 22 {
		t.Errorf("DiscreteLog(2,22,29)=%d,%v invalid", x, ok)
	}
	// No solution: 2^x ≡ 3 (mod 7) has none since <2> = {1,2,4}.
	if _, ok := DiscreteLog(2, 3, 7); ok {
		t.Errorf("DiscreteLog(2,3,7) want ok=false")
	}
}

func TestOrder(t *testing.T) {
	// ord_7(3) = 6 (primitive root), ord_7(2) = 3.
	if o, ok := Order(3, 7); !ok || o != 6 {
		t.Errorf("Order(3,7)=%d,%v want 6,true", o, ok)
	}
	if o, ok := Order(2, 7); !ok || o != 3 {
		t.Errorf("Order(2,7)=%d,%v want 3,true", o, ok)
	}
	if _, ok := Order(2, 4); ok {
		t.Errorf("Order(2,4) want ok=false")
	}
	if o, ok := Order(1, 5); !ok || o != 1 {
		t.Errorf("Order(1,5)=%d,%v want 1,true", o, ok)
	}
}

// ---- combinatorics ----

func bi(x int64) *big.Int { return big.NewInt(x) }

func TestFactorial(t *testing.T) {
	if Factorial(0).Cmp(bi(1)) != 0 {
		t.Errorf("Factorial(0) want 1")
	}
	if Factorial(5).Cmp(bi(120)) != 0 {
		t.Errorf("Factorial(5) want 120")
	}
	want, _ := new(big.Int).SetString("2432902008176640000", 10)
	if Factorial(20).Cmp(want) != 0 {
		t.Errorf("Factorial(20)=%v want %v", Factorial(20), want)
	}
}

func TestDoubleFactorial(t *testing.T) {
	cases := []struct {
		n    int64
		want int64
	}{{0, 1}, {-1, 1}, {5, 15}, {6, 48}, {7, 105}, {8, 384}}
	for _, c := range cases {
		if DoubleFactorial(c.n).Cmp(bi(c.want)) != 0 {
			t.Errorf("DoubleFactorial(%d)=%v want %d", c.n, DoubleFactorial(c.n), c.want)
		}
	}
}

func TestBinomial(t *testing.T) {
	if Binomial(10, 3).Cmp(bi(120)) != 0 {
		t.Errorf("Binomial(10,3) want 120")
	}
	if Binomial(52, 5).Cmp(bi(2598960)) != 0 {
		t.Errorf("Binomial(52,5) want 2598960")
	}
	if Binomial(5, 6).Sign() != 0 {
		t.Errorf("Binomial(5,6) want 0")
	}
	if Binomial(5, -1).Sign() != 0 {
		t.Errorf("Binomial(5,-1) want 0")
	}
}

func TestMultinomial(t *testing.T) {
	// 10! / (2!3!5!) = 2520.
	if Multinomial(2, 3, 5).Cmp(bi(2520)) != 0 {
		t.Errorf("Multinomial(2,3,5) want 2520")
	}
	// Reduces to a binomial: (3+2 choose 2) = 10.
	if Multinomial(3, 2).Cmp(bi(10)) != 0 {
		t.Errorf("Multinomial(3,2) want 10")
	}
	if Multinomial().Cmp(bi(1)) != 0 {
		t.Errorf("Multinomial() want 1")
	}
}

func TestPermutations(t *testing.T) {
	if Permutations(10, 3).Cmp(bi(720)) != 0 {
		t.Errorf("Permutations(10,3) want 720")
	}
	if Permutations(5, 5).Cmp(bi(120)) != 0 {
		t.Errorf("Permutations(5,5) want 120")
	}
	if Permutations(5, 6).Sign() != 0 {
		t.Errorf("Permutations(5,6) want 0")
	}
}

func TestCatalanNumber(t *testing.T) {
	want := []int64{1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862}
	for n, w := range want {
		if CatalanNumber(int64(n)).Cmp(bi(w)) != 0 {
			t.Errorf("CatalanNumber(%d)=%v want %d", n, CatalanNumber(int64(n)), w)
		}
	}
}

func TestStirlingSecond(t *testing.T) {
	cases := []struct {
		n, k int64
		want int64
	}{
		{0, 0, 1}, {4, 2, 7}, {5, 3, 25}, {5, 5, 1}, {5, 0, 0}, {6, 3, 90}, {4, 5, 0},
	}
	for _, c := range cases {
		if StirlingSecond(c.n, c.k).Cmp(bi(c.want)) != 0 {
			t.Errorf("StirlingSecond(%d,%d)=%v want %d", c.n, c.k, StirlingSecond(c.n, c.k), c.want)
		}
	}
}

func TestPartition(t *testing.T) {
	want := []int64{1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42}
	for n, w := range want {
		if Partition(int64(n)).Cmp(bi(w)) != 0 {
			t.Errorf("Partition(%d)=%v want %d", n, Partition(int64(n)), w)
		}
	}
	if Partition(100).Cmp(bi(190569292)) != 0 {
		t.Errorf("Partition(100)=%v want 190569292", Partition(100))
	}
}

// ---- sequences ----

func TestFibonacci(t *testing.T) {
	want := []int64{0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55}
	for n, w := range want {
		if Fibonacci(int64(n)).Cmp(bi(w)) != 0 {
			t.Errorf("Fibonacci(%d)=%v want %d", n, Fibonacci(int64(n)), w)
		}
	}
	want100, _ := new(big.Int).SetString("354224848179261915075", 10)
	if Fibonacci(100).Cmp(want100) != 0 {
		t.Errorf("Fibonacci(100)=%v want %v", Fibonacci(100), want100)
	}
}

func TestLucas(t *testing.T) {
	want := []int64{2, 1, 3, 4, 7, 11, 18, 29, 47, 76}
	for n, w := range want {
		if Lucas(int64(n)).Cmp(bi(w)) != 0 {
			t.Errorf("Lucas(%d)=%v want %d", n, Lucas(int64(n)), w)
		}
	}
}

func TestTribonacci(t *testing.T) {
	want := []int64{0, 1, 1, 2, 4, 7, 13, 24, 44, 81, 149}
	for n, w := range want {
		if Tribonacci(int64(n)).Cmp(bi(w)) != 0 {
			t.Errorf("Tribonacci(%d)=%v want %d", n, Tribonacci(int64(n)), w)
		}
	}
}

func TestIsSquare(t *testing.T) {
	squares := map[int64]bool{0: true, 1: true, 4: true, 9: true, 16: true, 144: true, 1000000: true}
	for n := int64(0); n <= 20; n++ {
		want := squares[n]
		if IsSquare(n) != want {
			t.Errorf("IsSquare(%d)=%v want %v", n, IsSquare(n), want)
		}
	}
	if !IsSquare(1000000) {
		t.Errorf("IsSquare(1000000) want true")
	}
	if IsSquare(-4) {
		t.Errorf("IsSquare(-4) want false")
	}
}

func TestIsqrtBig(t *testing.T) {
	cases := []struct{ n, want int64 }{{0, 0}, {1, 1}, {15, 3}, {16, 4}, {17, 4}, {1000000, 1000}}
	for _, c := range cases {
		if got := IsqrtBig(bi(c.n)); got.Cmp(bi(c.want)) != 0 {
			t.Errorf("IsqrtBig(%d)=%v want %d", c.n, got, c.want)
		}
	}
}

func TestBernoulli(t *testing.T) {
	cases := []struct {
		n          int64
		num, denom int64
	}{
		{0, 1, 1}, {1, -1, 2}, {2, 1, 6}, {3, 0, 1}, {4, -1, 30}, {6, 1, 42}, {8, -1, 30}, {10, 5, 66},
	}
	for _, c := range cases {
		want := big.NewRat(c.num, c.denom)
		if got := Bernoulli(c.n); got.Cmp(want) != 0 {
			t.Errorf("Bernoulli(%d)=%v want %v", c.n, got, want)
		}
	}
}
