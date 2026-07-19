package contfrac_test

import (
	"fmt"

	"github.com/malcolmston/algebra/contfrac"
)

func ExampleFromRational() {
	cf := contfrac.FromRational(415, 93)
	fmt.Println(cf)
	// Output: [4; 2, 6, 7]
}

func ExampleBestApproximation() {
	p, q := contfrac.BestApproximation(3.141592653589793, 113)
	fmt.Printf("%d/%d\n", p, q)
	// Output: 355/113
}

func ExampleSqrtCF() {
	fmt.Println(contfrac.SqrtCF(7))
	// Output: [2; (1, 1, 1, 4)]
}

func ExamplePellFundamental() {
	x, y, _ := contfrac.PellFundamental(61)
	fmt.Printf("%s^2 - 61*%s^2 = 1\n", x, y)
	// Output: 1766319049^2 - 61*226153980^2 = 1
}

func ExampleSternBrocotPath() {
	fmt.Println(contfrac.SternBrocotPath(3, 2))
	// Output: RL
}

func ExampleFareySequence() {
	fmt.Println(contfrac.FareySequence(4))
	// Output: [0 1/4 1/3 1/2 2/3 3/4 1]
}

func ExampleEgyptianFraction() {
	dens, _ := contfrac.EgyptianFraction(5, 6)
	fmt.Println(dens)
	// Output: [2 3]
}

func ExamplePiCF() {
	fmt.Println(contfrac.PiCF(5))
	// Output: [3; 7, 15, 1, 292]
}
