// Command examples demonstrates the algebra package end to end: parsing,
// differentiation, simplification, expansion, substitution, numeric
// evaluation, integration and solving a quadratic.
//
// Run it with:
//
//	go run ./examples
package main

import (
	"fmt"

	alg "github.com/malcolmston/algebra"
)

func main() {
	x := alg.Sym("x")

	fmt.Println("== Parse, differentiate, simplify, evaluate ==")
	expr, err := alg.Parse("x^2 + 2*x + 1")
	if err != nil {
		panic(err)
	}
	fmt.Printf("f(x)          = %s\n", expr)

	d := alg.Diff(expr, x)
	fmt.Printf("f'(x)         = %s\n", d)

	fmt.Printf("f'(x) simpl.  = %s\n", alg.Simplify(d))

	val, _ := alg.Eval(expr, map[string]float64{"x": 3})
	fmt.Printf("f(3)          = %g\n", val)

	fmt.Println()
	fmt.Println("== Expand and factor ==")
	sq := alg.MustParse("(x + 1)^2")
	fmt.Printf("(x+1)^2       = %s\n", alg.Expand(sq))
	fmt.Printf("factor x^2-1  = %s\n", alg.Factor(alg.MustParse("x^2 - 1"), x))

	fmt.Println()
	fmt.Println("== Chain / product rules and elementary functions ==")
	fmt.Printf("d/dx sin(x^2) = %s\n", alg.Diff(alg.MustParse("sin(x^2)"), x))
	fmt.Printf("d/dx x*exp(x) = %s\n", alg.Diff(alg.MustParse("x*exp(x)"), x))
	fmt.Printf("d/dx log(x)   = %s\n", alg.Diff(alg.MustParse("log(x)"), x))

	fmt.Println()
	fmt.Println("== Substitution and integration ==")
	fmt.Printf("subs y->x+1   = %s\n",
		alg.Subs(alg.MustParse("y^2"), alg.Sym("y"), alg.MustParse("x + 1")).Expand())
	fmt.Printf("int (x^2) dx  = %s\n", alg.Integrate(alg.MustParse("x^2"), x))
	fmt.Printf("int cos(x) dx = %s\n", alg.Integrate(alg.MustParse("cos(x)"), x))

	fmt.Println()
	fmt.Println("== Solve a quadratic ==")
	roots, err := alg.Solve(alg.MustParse("x^2 - 5*x + 6"), x)
	if err != nil {
		panic(err)
	}
	fmt.Printf("x^2-5x+6 = 0  -> ")
	for i, r := range roots {
		if i > 0 {
			fmt.Print(", ")
		}
		fmt.Print(r)
	}
	fmt.Println()

	roots2, _ := alg.Solve(alg.MustParse("x^2 - 2"), x)
	fmt.Printf("x^2-2 = 0     -> ")
	for i, r := range roots2 {
		if i > 0 {
			fmt.Print(", ")
		}
		fmt.Print(r)
	}
	fmt.Println()
}
