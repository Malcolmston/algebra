package padic

import (
	"fmt"
	"math/big"
)

// ExampleFromRational shows the p-adic value of a rational and its valuation.
func ExampleFromRational() {
	p := big.NewInt(5)
	x, _ := FromRational(p, big.NewInt(1), big.NewInt(2), 4)
	// 1/2 in Q_5 is a unit; its canonical integer representative mod 5^4.
	rep, _ := x.BigInt()
	fmt.Printf("val=%d rep=%d\n", x.Valuation(), rep)
	// 2 * rep == 1 mod 5^4
	fmt.Println(new(big.Int).Mod(new(big.Int).Mul(rep, big.NewInt(2)), PPow(p, 4)))
	// Output:
	// val=0 rep=313
	// 1
}

// ExamplePadic_Sqrt computes a square root of 2 in the 7-adic numbers.
func ExamplePadic_Sqrt() {
	p := big.NewInt(7)
	two := FromInt(p, 2, 4)
	root, _ := two.Sqrt()
	back := root.Square()
	fmt.Println(back.Equal(two))
	// Output: true
}

// ExampleHenselLift lifts a root of x^2 - 2 modulo 7 to higher precision.
func ExampleHenselLift() {
	p := big.NewInt(7)
	f := []*big.Int{big.NewInt(-2), big.NewInt(0), big.NewInt(1)} // x^2 - 2
	root, _ := HenselLift(f, big.NewInt(3), p, 4)
	val := PolyEval(f, root)
	val.Mod(val, PPow(p, 4))
	fmt.Printf("f(%d) == 0 mod 7^4: %v\n", root, val.Sign() == 0)
	// Output: f(2166) == 0 mod 7^4: true
}

// ExampleNewtonPolygonFromInts reads root valuations off a Newton polygon.
func ExampleNewtonPolygonFromInts() {
	// x^2 - 6 over the prime 3; 6 = 2*3, so both roots have valuation 1/2.
	np := NewtonPolygonFromInts(big.NewInt(3), []*big.Int{
		big.NewInt(-6), big.NewInt(0), big.NewInt(1),
	})
	fmt.Println(np.RootValuations())
	// Output: [1/2 1/2]
}

// ExampleTeichmullerRep computes a Teichmuller representative.
func ExampleTeichmullerRep() {
	p := big.NewInt(5)
	rep, _ := TeichmullerRep(p, big.NewInt(2), 3)
	// The representative is a 4th root of unity congruent to 2 mod 5.
	fmt.Println(rep)
	// Output: 57
}
