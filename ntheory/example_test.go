package ntheory_test

import (
	"fmt"

	"github.com/malcolmston/algebra/ntheory"
)

func ExampleGCD() {
	fmt.Println(ntheory.GCD(48, 36))
	// Output: 12
}

func ExampleLCM() {
	fmt.Println(ntheory.LCM(4, 6))
	// Output: 12
}

func ExampleExtendedGCD() {
	g, x, y := ntheory.ExtendedGCD(240, 46)
	fmt.Printf("gcd=%d, 240*%d + 46*%d = %d\n", g, x, y, 240*x+46*y)
	// Output: gcd=2, 240*-9 + 46*47 = 2
}

func ExampleDivisors() {
	fmt.Println(ntheory.Divisors(28))
	// Output: [1 2 4 7 14 28]
}

func ExampleSumDivisors() {
	fmt.Println(ntheory.SumDivisors(28))
	// Output: 56
}

func ExampleIsPerfect() {
	fmt.Println(ntheory.IsPerfect(28), ntheory.IsPerfect(12))
	// Output: true false
}

func ExampleIsPrime() {
	// 561 is a Carmichael number: composite despite passing Fermat's test.
	fmt.Println(ntheory.IsPrime(97), ntheory.IsPrime(561))
	// Output: true false
}

func ExamplePrimesUpTo() {
	fmt.Println(ntheory.PrimesUpTo(30))
	// Output: [2 3 5 7 11 13 17 19 23 29]
}

func ExamplePrimePi() {
	fmt.Println(ntheory.PrimePi(100))
	// Output: 25
}

func ExampleFactorize() {
	fmt.Println(ntheory.FactorList(360))
	// Output: [{2 3} {3 2} {5 1}]
}

func ExampleEulerPhi() {
	fmt.Println(ntheory.EulerPhi(10))
	// Output: 4
}

func ExampleMobiusMu() {
	fmt.Println(ntheory.MobiusMu(30), ntheory.MobiusMu(12))
	// Output: -1 0
}

func ExampleModPow() {
	fmt.Println(ntheory.ModPow(2, 10, 1000))
	// Output: 24
}

func ExampleModInverse() {
	inv, ok := ntheory.ModInverse(3, 11)
	fmt.Println(inv, ok)
	// Output: 4 true
}

func ExampleCRT() {
	x, m, ok := ntheory.CRT([]int64{2, 3, 2}, []int64{3, 5, 7})
	fmt.Println(x, m, ok)
	// Output: 23 105 true
}

func ExampleBinomial() {
	fmt.Println(ntheory.Binomial(10, 3))
	// Output: 120
}

func ExampleFactorial() {
	fmt.Println(ntheory.Factorial(20))
	// Output: 2432902008176640000
}

func ExampleFibonacci() {
	fmt.Println(ntheory.Fibonacci(10), ntheory.Fibonacci(100))
	// Output: 55 354224848179261915075
}

func ExampleCatalanNumber() {
	fmt.Println(ntheory.CatalanNumber(5))
	// Output: 42
}

func ExamplePartition() {
	fmt.Println(ntheory.Partition(100))
	// Output: 190569292
}

func ExampleBernoulli() {
	fmt.Println(ntheory.Bernoulli(2).RatString(), ntheory.Bernoulli(10).RatString())
	// Output: 1/6 5/66
}
