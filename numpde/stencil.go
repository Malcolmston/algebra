package numpde

// Stencil holds the offset/weight pairs of a finite-difference approximation to
// a derivative on a uniform grid. Offsets are integer multiples of the grid
// spacing relative to the point being updated, and Weights are the matching
// coefficients already divided by the appropriate power of the spacing, so that
// the derivative estimate is sum_k Weights[k]*u[i+Offsets[k]].
type Stencil struct {
	Offsets []int
	Weights []float64
}

// Apply evaluates the stencil at interior index i of the slice u. The caller
// must ensure every i+Offsets[k] is a valid index; Apply performs no bounds
// checking beyond what the Go runtime provides.
func (s Stencil) Apply(u []float64, i int) float64 {
	sum := 0.0
	for k, off := range s.Offsets {
		sum += s.Weights[k] * u[i+off]
	}
	return sum
}

// CentralFirstDerivativeStencil returns the second-order accurate central
// difference stencil for d/dx on a grid of spacing dx: weights [-1,0,1]/(2dx)
// at offsets [-1,0,1].
func CentralFirstDerivativeStencil(dx float64) Stencil {
	c := 1.0 / (2 * dx)
	return Stencil{Offsets: []int{-1, 0, 1}, Weights: []float64{-c, 0, c}}
}

// ForwardFirstDerivativeStencil returns the first-order accurate forward
// difference stencil for d/dx: weights [-1,1]/dx at offsets [0,1].
func ForwardFirstDerivativeStencil(dx float64) Stencil {
	c := 1.0 / dx
	return Stencil{Offsets: []int{0, 1}, Weights: []float64{-c, c}}
}

// BackwardFirstDerivativeStencil returns the first-order accurate backward
// difference stencil for d/dx: weights [-1,1]/dx at offsets [-1,0].
func BackwardFirstDerivativeStencil(dx float64) Stencil {
	c := 1.0 / dx
	return Stencil{Offsets: []int{-1, 0}, Weights: []float64{-c, c}}
}

// SecondDerivativeStencil returns the second-order accurate central difference
// stencil for d^2/dx^2: weights [1,-2,1]/dx^2 at offsets [-1,0,1].
func SecondDerivativeStencil(dx float64) Stencil {
	c := 1.0 / (dx * dx)
	return Stencil{Offsets: []int{-1, 0, 1}, Weights: []float64{c, -2 * c, c}}
}

// FivePointLaplacianStencil returns the coefficients of the standard
// second-order five-point discrete Laplacian on a grid with spacings dx and dy.
// The returned values are the weights applied to the centre node and its four
// axis-aligned neighbours: laplacian(u)[i][j] = center*u[i][j] +
// east*u[i+1][j] + west*u[i-1][j] + north*u[i][j+1] + south*u[i][j-1].
func FivePointLaplacianStencil(dx, dy float64) (center, east, west, north, south float64) {
	ix := 1.0 / (dx * dx)
	iy := 1.0 / (dy * dy)
	center = -2*ix - 2*iy
	east, west = ix, ix
	north, south = iy, iy
	return
}

// NinePointLaplacianStencil returns the weights of the fourth-order accurate
// nine-point discrete Laplacian on a square grid of spacing h. The order of the
// returned weights is: centre, the four edge neighbours (E, W, N, S share the
// value edge), and the four diagonal neighbours (share the value corner). The
// operator is laplacian(u) = (corner*(sum of 4 diagonals) + edge*(sum of 4
// edges) + center*u_ij)/(6 h^2) already folded into the returned weights.
func NinePointLaplacianStencil(h float64) (center, edge, corner float64) {
	inv := 1.0 / (6 * h * h)
	center = -20 * inv
	edge = 4 * inv
	corner = 1 * inv
	return
}

// Laplacian1D returns the discrete second derivative of the field u on a grid
// of spacing dx using the three-point central stencil at interior points. The
// endpoints of the returned slice are set to zero, matching a homogeneous
// Dirichlet treatment; adjust them afterwards if a different boundary is
// required. It panics if len(u) < 3.
func Laplacian1D(u []float64, dx float64) []float64 {
	n := len(u)
	if n < 3 {
		panic("numpde: Laplacian1D requires len(u) >= 3")
	}
	out := make([]float64, n)
	inv := 1.0 / (dx * dx)
	for i := 1; i < n-1; i++ {
		out[i] = (u[i-1] - 2*u[i] + u[i+1]) * inv
	}
	return out
}

// Laplacian2D returns the discrete Laplacian of the field u (stored as u[i][j])
// using the five-point stencil with spacings dx and dy. Interior nodes receive
// the standard approximation; boundary nodes are left at zero. It panics if the
// grid is smaller than 3x3.
func Laplacian2D(u [][]float64, dx, dy float64) [][]float64 {
	nx := len(u)
	if nx < 3 || len(u[0]) < 3 {
		panic("numpde: Laplacian2D requires at least a 3x3 grid")
	}
	ny := len(u[0])
	out := Zeros2D(nx, ny)
	ix := 1.0 / (dx * dx)
	iy := 1.0 / (dy * dy)
	for i := 1; i < nx-1; i++ {
		for j := 1; j < ny-1; j++ {
			out[i][j] = (u[i-1][j]-2*u[i][j]+u[i+1][j])*ix +
				(u[i][j-1]-2*u[i][j]+u[i][j+1])*iy
		}
	}
	return out
}

// Gradient1D returns the discrete first derivative of u on a grid of spacing dx.
// Interior points use the second-order central difference; the endpoints use
// first-order one-sided differences so the result is defined everywhere. It
// panics if len(u) < 2.
func Gradient1D(u []float64, dx float64) []float64 {
	n := len(u)
	if n < 2 {
		panic("numpde: Gradient1D requires len(u) >= 2")
	}
	out := make([]float64, n)
	for i := 1; i < n-1; i++ {
		out[i] = (u[i+1] - u[i-1]) / (2 * dx)
	}
	out[0] = (u[1] - u[0]) / dx
	out[n-1] = (u[n-1] - u[n-2]) / dx
	return out
}
