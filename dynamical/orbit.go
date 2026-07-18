package dynamical

// Orbit returns the forward orbit x_0, x_1, ..., x_n of the one-dimensional
// map f starting from x0, that is a slice of length n+1 whose first element is
// x0 and whose k-th element is the k-th iterate of f applied to x0.
func Orbit(f Map1D, x0 float64, n int) []float64 {
	if n < 0 {
		n = 0
	}
	out := make([]float64, n+1)
	x := x0
	out[0] = x
	for i := 1; i <= n; i++ {
		x = f(x)
		out[i] = x
	}
	return out
}

// OrbitTransient iterates the map f from x0 discarding the first transient
// points, then returns the following n+1 points of the orbit. It is used to
// sample the attractor of f while suppressing initial transient behavior.
func OrbitTransient(f Map1D, x0 float64, transient, n int) []float64 {
	x := x0
	for i := 0; i < transient; i++ {
		x = f(x)
	}
	return Orbit(f, x, n)
}

// NthIterate returns the n-th iterate f^n(x0) of the map f. For n <= 0 it
// returns x0 unchanged.
func NthIterate(f Map1D, x0 float64, n int) float64 {
	x := x0
	for i := 0; i < n; i++ {
		x = f(x)
	}
	return x
}

// Orbit2D returns the forward orbit p_0, p_1, ..., p_n of the two-dimensional
// map f starting from p0, as a slice of length n+1.
func Orbit2D(f Map2D, p0 Point2D, n int) []Point2D {
	if n < 0 {
		n = 0
	}
	out := make([]Point2D, n+1)
	p := p0
	out[0] = p
	for i := 1; i <= n; i++ {
		p = f(p)
		out[i] = p
	}
	return out
}

// OrbitTransient2D iterates the map f from p0 discarding the first transient
// points, then returns the following n+1 points of the orbit.
func OrbitTransient2D(f Map2D, p0 Point2D, transient, n int) []Point2D {
	p := p0
	for i := 0; i < transient; i++ {
		p = f(p)
	}
	return Orbit2D(f, p, n)
}

// NthIterate2D returns the n-th iterate f^n(p0) of the two-dimensional map f.
func NthIterate2D(f Map2D, p0 Point2D, n int) Point2D {
	p := p0
	for i := 0; i < n; i++ {
		p = f(p)
	}
	return p
}
