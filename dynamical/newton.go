package dynamical

import (
	"math"
	"math/cmplx"
)

// NewtonResult holds the outcome of a complex Newton iteration: the point Root
// reached, the number of Iterations performed, and whether the iteration
// Converged to within the requested tolerance.
type NewtonResult struct {
	Root       complex128
	Iterations int
	Converged  bool
}

// NewtonComplex applies Newton's method z -> z - f(z)/f'(z) in the complex
// plane, starting at z0, using the function f and its derivative df. It stops
// when |f(z)| <= tol (reporting convergence) or after maxIter steps.
func NewtonComplex(f, df func(complex128) complex128, z0 complex128, maxIter int, tol float64) NewtonResult {
	z := z0
	for i := 1; i <= maxIter; i++ {
		fz := f(z)
		if cmplx.Abs(fz) <= tol {
			return NewtonResult{Root: z, Iterations: i, Converged: true}
		}
		d := df(z)
		if d == 0 {
			return NewtonResult{Root: z, Iterations: i, Converged: false}
		}
		z -= fz / d
	}
	return NewtonResult{Root: z, Iterations: maxIter, Converged: cmplx.Abs(f(z)) <= tol}
}

// NewtonReal applies Newton's method x -> x - f(x)/f'(x) on the real line,
// starting at x0, using the function f and its derivative df. It returns the
// located root, the number of iterations performed, and whether |f(x)| <= tol
// was achieved.
func NewtonReal(f, df Map1D, x0 float64, maxIter int, tol float64) (root float64, iters int, converged bool) {
	x := x0
	for i := 1; i <= maxIter; i++ {
		fx := f(x)
		if math.Abs(fx) <= tol {
			return x, i, true
		}
		d := df(x)
		if d == 0 {
			return x, i, false
		}
		x -= fx / d
	}
	return x, maxIter, math.Abs(f(x)) <= tol
}

// ClassifyRoot returns the index in roots of the entry nearest to z, provided
// that nearest entry lies within tol; otherwise it returns -1. It is used to
// label which basin of attraction a Newton iterate has fallen into.
func ClassifyRoot(roots []complex128, z complex128, tol float64) int {
	best := -1
	bestD := math.Inf(1)
	for i, r := range roots {
		d := cmplx.Abs(z - r)
		if d < bestD {
			bestD = d
			best = i
		}
	}
	if best >= 0 && bestD <= tol {
		return best
	}
	return -1
}

// NewtonBasin computes the Newton basin-of-attraction grid for the complex map
// f with derivative df over the rectangle [xmin,xmax] x [ymin,ymax], sampled on
// an nx-by-ny grid. Each grid point is used as a starting value for Newton's
// method; the resulting root is matched against the supplied roots by
// [ClassifyRoot]. The returned matrix has ny rows and nx columns in row-major
// order: entry [j][i] is the index of the root that the sample at column i,
// row j converged to, or -1 if it did not converge to any listed root.
func NewtonBasin(f, df func(complex128) complex128, roots []complex128, xmin, xmax, ymin, ymax float64, nx, ny, maxIter int, tol float64) [][]int {
	grid := make([][]int, ny)
	for j := 0; j < ny; j++ {
		row := make([]int, nx)
		var y float64
		if ny == 1 {
			y = ymin
		} else {
			y = ymin + (ymax-ymin)*float64(j)/float64(ny-1)
		}
		for i := 0; i < nx; i++ {
			var x float64
			if nx == 1 {
				x = xmin
			} else {
				x = xmin + (xmax-xmin)*float64(i)/float64(nx-1)
			}
			res := NewtonComplex(f, df, complex(x, y), maxIter, tol)
			if !res.Converged {
				row[i] = -1
				continue
			}
			row[i] = ClassifyRoot(roots, res.Root, tol*1e3)
		}
		grid[j] = row
	}
	return grid
}
