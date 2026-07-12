package matrix_test

import (
	"fmt"

	"github.com/malcolmston/algebra"
	"github.com/malcolmston/algebra/matrix"
)

func ExampleIdentity() {
	fmt.Println(matrix.Identity(3))
	// Output:
	// [ 1  0  0 ]
	// [ 0  1  0 ]
	// [ 0  0  1 ]
}

func ExampleMatrix_Det() {
	m := matrix.FromInts([][]int64{{1, 2}, {3, 4}})
	d, _ := m.Det()
	fmt.Println(d)
	// Output: -2
}

func ExampleMatrix_Inverse() {
	m := matrix.FromInts([][]int64{{4, 7}, {2, 6}})
	inv, _ := m.Inverse()
	prod, _ := m.Mul(inv)
	fmt.Println(prod)
	// Output:
	// [ 1  0 ]
	// [ 0  1 ]
}

func ExampleSolve() {
	// 2x + y = 5 ; x - y = 1.
	a := matrix.FromInts([][]int64{{2, 1}, {1, -1}})
	b := matrix.VectorFromInts(5, 1)
	x, _ := matrix.Solve(a, b)
	fmt.Println(x)
	// Output: [2, 1]
}

func ExampleVector_Cross() {
	e1 := matrix.VectorFromInts(1, 0, 0)
	e2 := matrix.VectorFromInts(0, 1, 0)
	c, _ := e1.Cross(e2)
	fmt.Println(c)
	// Output: [0, 0, 1]
}

func ExampleMatrix_CharPoly() {
	m := matrix.FromInts([][]int64{{2, 0}, {0, 3}})
	p, _ := m.CharPoly()
	fmt.Println(p)
	// Output: lambda^2 - 5*lambda + 6
}

func ExampleMatrix_Eigenvalues() {
	m := matrix.FromInts([][]int64{{2, 0}, {0, 3}})
	ev, _ := m.Eigenvalues()
	for _, e := range ev {
		fmt.Println(e)
	}
	// Output:
	// 2
	// 3
}

func ExampleMatrix_symbolic() {
	// Symbolic 2x2 determinant.
	a := algebra.Sym("a")
	b := algebra.Sym("b")
	c := algebra.Sym("c")
	d := algebra.Sym("d")
	m := matrix.FromExpr([][]algebra.Expr{{a, b}, {c, d}})
	det, _ := m.Det()
	fmt.Println(det)
	// Output: a*d - b*c
}
