package combin

import (
	"math"
	"math/big"
	"testing"
)

// bigStr is a helper that compares a *big.Int result to an expected decimal
// string, failing the test with context when they differ.
func bigStr(t *testing.T, name string, got *big.Int, want string) {
	t.Helper()
	if got.String() != want {
		t.Errorf("%s = %s, want %s", name, got.String(), want)
	}
}

func TestFactorials(t *testing.T) {
	factCases := []struct {
		n    int
		want string
	}{
		{0, "1"}, {1, "1"}, {5, "120"}, {10, "3628800"},
		{20, "2432902008176640000"},
		{25, "15511210043330985984000000"},
	}
	for _, c := range factCases {
		bigStr(t, "Factorial", Factorial(c.n), c.want)
	}

	dfCases := []struct {
		n    int
		want string
	}{
		{-1, "1"}, {0, "1"}, {1, "1"}, {5, "15"}, {6, "48"}, {9, "945"}, {10, "3840"},
	}
	for _, c := range dfCases {
		bigStr(t, "DoubleFactorial", DoubleFactorial(c.n), c.want)
	}

	// Subfactorial / derangements: 1, 0, 1, 2, 9, 44, 265, 1854.
	subCases := []string{"1", "0", "1", "2", "9", "44", "265", "1854"}
	for n, want := range subCases {
		bigStr(t, "Subfactorial", Subfactorial(n), want)
		bigStr(t, "Derangement", Derangement(n), want)
	}

	bigStr(t, "Multifactorial(9,3)", Multifactorial(9, 3), "162")   // 9*6*3
	bigStr(t, "Multifactorial(10,3)", Multifactorial(10, 3), "280") // 10*7*4*1
	bigStr(t, "Superfactorial(4)", Superfactorial(4), "288")
	bigStr(t, "HyperFactorial(4)", HyperFactorial(4), "27648")
	bigStr(t, "Primorial(10)", Primorial(10), "210")
	bigStr(t, "Primorial(11)", Primorial(11), "2310")

	if got := SubfactorialFloat(4); got != 9 {
		t.Errorf("SubfactorialFloat(4) = %v, want 9", got)
	}
	if got := LogFactorial(10); math.Abs(got-math.Log(3628800)) > 1e-9 {
		t.Errorf("LogFactorial(10) = %v, want %v", got, math.Log(3628800))
	}
	if got := FactorialFloat(6); math.Abs(got-720) > 1e-6 {
		t.Errorf("FactorialFloat(6) = %v, want 720", got)
	}
}

func TestPermutationsCombinations(t *testing.T) {
	bigStr(t, "Permutations(5,2)", Permutations(5, 2), "20")
	bigStr(t, "Permutations(10,3)", Permutations(10, 3), "720")
	bigStr(t, "Permutations(5,6)", Permutations(5, 6), "0")
	bigStr(t, "PermutationsWithRepetition(3,4)", PermutationsWithRepetition(3, 4), "81")

	binCases := []struct {
		n, k int
		want string
	}{
		{5, 2, "10"}, {10, 5, "252"}, {52, 5, "2598960"},
		{6, 0, "1"}, {6, 6, "1"}, {6, 7, "0"},
		{100, 50, "100891344545564193334812497256"},
	}
	for _, c := range binCases {
		bigStr(t, "Binomial", Binomial(c.n, c.k), c.want)
		bigStr(t, "Combinations", Combinations(c.n, c.k), c.want)
	}

	bigStr(t, "CombinationsWithRepetition(5,3)", CombinationsWithRepetition(5, 3), "35")
	bigStr(t, "CentralBinomial(5)", CentralBinomial(5), "252")

	if got := BinomialFloat(52, 5); got != 2598960 {
		t.Errorf("BinomialFloat(52,5) = %v, want 2598960", got)
	}
	if got := PermutationsFloat(10, 3); got != 720 {
		t.Errorf("PermutationsFloat(10,3) = %v, want 720", got)
	}
	if got := LogBinomial(10, 5); math.Abs(got-math.Log(252)) > 1e-9 {
		t.Errorf("LogBinomial(10,5) = %v, want %v", got, math.Log(252))
	}
	if got := GammaBinomial(5, 2); math.Abs(got-10) > 1e-9 {
		t.Errorf("GammaBinomial(5,2) = %v, want 10", got)
	}
	if got := GammaBinomial(6.5, 2); math.Abs(got-17.875) > 1e-9 {
		t.Errorf("GammaBinomial(6.5,2) = %v, want 17.875", got)
	}
}

func TestMultinomial(t *testing.T) {
	bigStr(t, "Multinomial(2,1)", Multinomial(2, 1), "3")
	bigStr(t, "Multinomial(1,2,3)", Multinomial(1, 2, 3), "60")
	bigStr(t, "Multinomial(2,3,4)", Multinomial(2, 3, 4), "1260")
	// MISSISSIPPI: 11!/(1!4!4!2!) = 34650.
	bigStr(t, "Multinomial MISSISSIPPI", Multinomial(1, 4, 4, 2), "34650")

	if got := MultinomialFloat(1, 2, 3); got != 60 {
		t.Errorf("MultinomialFloat(1,2,3) = %v, want 60", got)
	}
	if got := LogMultinomial(2, 3, 4); math.Abs(got-math.Log(1260)) > 1e-9 {
		t.Errorf("LogMultinomial(2,3,4) = %v, want %v", got, math.Log(1260))
	}
}

func TestRisingFalling(t *testing.T) {
	bigStr(t, "RisingFactorial(2,3)", RisingFactorial(2, 3), "24")   // 2*3*4
	bigStr(t, "FallingFactorial(5,3)", FallingFactorial(5, 3), "60") // 5*4*3
	bigStr(t, "RisingFactorial(1,5)", RisingFactorial(1, 5), "120")  // = 5!
	bigStr(t, "RisingFactorial(3,0)", RisingFactorial(3, 0), "1")

	if got := RisingFactorialFloat(2.5, 3); math.Abs(got-2.5*3.5*4.5) > 1e-9 {
		t.Errorf("RisingFactorialFloat(2.5,3) = %v", got)
	}
	if got := FallingFactorialFloat(5.5, 2); math.Abs(got-5.5*4.5) > 1e-9 {
		t.Errorf("FallingFactorialFloat(5.5,2) = %v", got)
	}
}

func TestCatalanMotzkinSchroeder(t *testing.T) {
	catalan := []string{"1", "1", "2", "5", "14", "42", "132", "429", "1430", "4862", "16796"}
	for n, want := range catalan {
		bigStr(t, "Catalan", Catalan(n), want)
	}
	if got := CatalanFloat(5); got != 42 {
		t.Errorf("CatalanFloat(5) = %v, want 42", got)
	}

	motzkin := []string{"1", "1", "2", "4", "9", "21", "51", "127", "323"}
	for n, want := range motzkin {
		bigStr(t, "MotzkinNumber", MotzkinNumber(n), want)
	}

	schroeder := []string{"1", "2", "6", "22", "90", "394", "1806", "8558"}
	for n, want := range schroeder {
		bigStr(t, "SchroederNumber", SchroederNumber(n), want)
	}

	// Catalan's triangle row 4: 1, 4, 9, 14, 14 (sums to Catalan(5)=42).
	catTri := []string{"1", "4", "9", "14", "14"}
	for k, want := range catTri {
		bigStr(t, "CatalanTriangle(4,k)", CatalanTriangle(4, k), want)
	}
}

func TestBell(t *testing.T) {
	bell := []string{"1", "1", "2", "5", "15", "52", "203", "877", "4140", "21147", "115975"}
	for n, want := range bell {
		bigStr(t, "Bell", Bell(n), want)
	}
	// Row 3 of the Bell triangle: 5, 7, 10, 15.
	row := BellTriangleRow(3)
	wantRow := []string{"5", "7", "10", "15"}
	if len(row) != len(wantRow) {
		t.Fatalf("BellTriangleRow(3) length = %d, want %d", len(row), len(wantRow))
	}
	for i, want := range wantRow {
		bigStr(t, "BellTriangleRow(3)", row[i], want)
	}
}

func TestStirling(t *testing.T) {
	// Unsigned Stirling first kind, row 4: 0, 6, 11, 6, 1.
	first4 := []string{"0", "6", "11", "6", "1"}
	for k, want := range first4 {
		bigStr(t, "StirlingFirstUnsigned(4,k)", StirlingFirstUnsigned(4, k), want)
	}
	// Signed variants: s(4,1) = -6, s(4,2) = 11, s(4,3) = -6.
	bigStr(t, "StirlingFirst(4,1)", StirlingFirst(4, 1), "-6")
	bigStr(t, "StirlingFirst(4,2)", StirlingFirst(4, 2), "11")
	bigStr(t, "StirlingFirst(4,3)", StirlingFirst(4, 3), "-6")

	// Stirling second kind, row 4: 0, 1, 7, 6, 1.
	second4 := []string{"0", "1", "7", "6", "1"}
	for k, want := range second4 {
		bigStr(t, "StirlingSecond(4,k)", StirlingSecond(4, k), want)
	}
	bigStr(t, "StirlingSecond(5,3)", StirlingSecond(5, 3), "25")

	// Row sums of Stirling second kind equal Bell numbers.
	sum := big.NewInt(0)
	for _, v := range StirlingSecondRow(6) {
		sum.Add(sum, v)
	}
	bigStr(t, "sum StirlingSecondRow(6)", sum, "203")

	row := StirlingFirstRow(4)
	if len(row) != 5 {
		t.Fatalf("StirlingFirstRow(4) length = %d, want 5", len(row))
	}
}

func TestLah(t *testing.T) {
	// Lah numbers row 4: L(4,1)=24, L(4,2)=36, L(4,3)=12, L(4,4)=1.
	lah := map[int]string{1: "24", 2: "36", 3: "12", 4: "1"}
	for k, want := range lah {
		bigStr(t, "LahNumber(4,k)", LahNumber(4, k), want)
	}
	bigStr(t, "LahNumber(0,0)", LahNumber(0, 0), "1")
	bigStr(t, "SignedLahNumber(4,2)", SignedLahNumber(4, 2), "36")
	bigStr(t, "SignedLahNumber(3,2)", SignedLahNumber(3, 2), "-6")
}

func TestPascal(t *testing.T) {
	row := PascalRow(6)
	want := []string{"1", "6", "15", "20", "15", "6", "1"}
	if len(row) != len(want) {
		t.Fatalf("PascalRow(6) length = %d, want %d", len(row), len(want))
	}
	for i, w := range want {
		bigStr(t, "PascalRow(6)", row[i], w)
		bigStr(t, "BinomialRow(6)", BinomialRow(6)[i], w)
	}
	tri := PascalTriangle(4)
	if len(tri) != 5 {
		t.Fatalf("PascalTriangle(4) has %d rows, want 5", len(tri))
	}
	bigStr(t, "PascalTriangle[4][2]", tri[4][2], "6")
}

func TestSequences(t *testing.T) {
	fib := []string{"0", "1", "1", "2", "3", "5", "8", "13", "21", "34", "55"}
	for n, want := range fib {
		bigStr(t, "Fibonacci", Fibonacci(n), want)
	}
	lucas := []string{"2", "1", "3", "4", "7", "11", "18", "29", "47"}
	for n, want := range lucas {
		bigStr(t, "Lucas", Lucas(n), want)
	}
	trib := []string{"0", "0", "1", "1", "2", "4", "7", "13", "24", "44"}
	for n, want := range trib {
		bigStr(t, "Tribonacci", Tribonacci(n), want)
	}
	tel := []string{"1", "1", "2", "4", "10", "26", "76", "232", "764"}
	for n, want := range tel {
		bigStr(t, "TelephoneNumber", TelephoneNumber(n), want)
	}
}

func TestPartitionsCompositions(t *testing.T) {
	part := []string{"1", "1", "2", "3", "5", "7", "11", "15", "22", "30", "42"}
	for n, want := range part {
		bigStr(t, "PartitionCount", PartitionCount(n), want)
	}
	bigStr(t, "PartitionCount(50)", PartitionCount(50), "204226")
	bigStr(t, "PartitionCount(100)", PartitionCount(100), "190569292")

	bigStr(t, "PartitionsIntoKParts(5,2)", PartitionsIntoKParts(5, 2), "2")
	bigStr(t, "PartitionsIntoKParts(7,3)", PartitionsIntoKParts(7, 3), "4")

	// Sum over k of p(n,k) equals p(n).
	sum := big.NewInt(0)
	for k := 0; k <= 7; k++ {
		sum.Add(sum, PartitionsIntoKParts(7, k))
	}
	bigStr(t, "sum PartitionsIntoKParts(7,*)", sum, "15")

	bigStr(t, "Compositions(0)", Compositions(0), "1")
	bigStr(t, "Compositions(4)", Compositions(4), "8")
	bigStr(t, "WeakCompositions(3,3)", WeakCompositions(3, 3), "10")
}

func TestTriangles(t *testing.T) {
	// Narayana row 4: 1, 6, 6, 1 (sums to Catalan(4)=14).
	nar := []string{"1", "6", "6", "1"}
	for i, want := range nar {
		bigStr(t, "NarayanaNumber(4,k)", NarayanaNumber(4, i+1), want)
	}

	// Eulerian row 4: 1, 11, 11, 1.
	eul := []string{"1", "11", "11", "1"}
	for k, want := range eul {
		bigStr(t, "EulerianNumber(4,k)", EulerianNumber(4, k), want)
	}
	// Row sum equals n!.
	sum := big.NewInt(0)
	for k := 0; k < 5; k++ {
		sum.Add(sum, EulerianNumber(5, k))
	}
	bigStr(t, "sum EulerianNumber(5,*)", sum, "120")

	// Second-order Eulerian row 4: 1, 22, 58, 24.
	eul2 := []string{"1", "22", "58", "24"}
	for k, want := range eul2 {
		bigStr(t, "EulerianSecondOrder(4,k)", EulerianSecondOrder(4, k), want)
	}

	bigStr(t, "DelannoyNumber(3,3)", DelannoyNumber(3, 3), "63")
	bigStr(t, "DelannoyNumber(2,2)", DelannoyNumber(2, 2), "13")

	// Rencontres row 4: 9, 8, 6, 0, 1 (sums to 4! = 24).
	ren := []string{"9", "8", "6", "0", "1"}
	rsum := big.NewInt(0)
	for k, want := range ren {
		bigStr(t, "RencontresNumber(4,k)", RencontresNumber(4, k), want)
		rsum.Add(rsum, RencontresNumber(4, k))
	}
	bigStr(t, "sum RencontresNumber(4,*)", rsum, "24")
}

func TestHarmonicAndAsymptotics(t *testing.T) {
	if got := HarmonicNumber(4); math.Abs(got-(1+0.5+1.0/3+0.25)) > 1e-12 {
		t.Errorf("HarmonicNumber(4) = %v", got)
	}
	if got := HarmonicNumber(0); got != 0 {
		t.Errorf("HarmonicNumber(0) = %v, want 0", got)
	}
	if got := HarmonicNumberGeneralized(3, 2); math.Abs(got-(1+0.25+1.0/9)) > 1e-12 {
		t.Errorf("HarmonicNumberGeneralized(3,2) = %v", got)
	}
	// Stirling's approximation should be within 1% of the true factorial.
	exact, _ := new(big.Float).SetInt(Factorial(20)).Float64()
	approx := StirlingApproximation(20)
	if math.Abs(approx-exact)/exact > 0.01 {
		t.Errorf("StirlingApproximation(20) = %v, exact = %v (relative error too large)", approx, exact)
	}
}

// BenchmarkPartitionCount exercises the heaviest routine, the arbitrary
// precision integer-partition count driven by Euler's pentagonal recurrence.
func BenchmarkPartitionCount(b *testing.B) {
	for i := 0; i < b.N; i++ {
		_ = PartitionCount(500)
	}
}
