package numint

import (
	"math"
	"math/rand"
	"sort"
)

// Path2 is a parametric plane curve t -> (x(t), y(t)).
type Path2 func(t float64) (x, y float64)

// Path3 is a parametric space curve t -> (x(t), y(t), z(t)).
type Path3 func(t float64) (x, y, z float64)

// Field2 is a two-dimensional vector field (x, y) -> (u, v).
type Field2 func(x, y float64) (u, v float64)

// Field3 is a three-dimensional vector field (x, y, z) -> (u, v, w).
type Field3 func(x, y, z float64) (u, v, w float64)

// --- unexported helpers (all prefixed numintmv to avoid sibling collisions) ---

// numintmvSign returns |a| with the sign of b.
func numintmvSign(a, b float64) float64 { return math.Copysign(a, b) }

// numintmvEvenUp returns n rounded up to an even value that is at least 2.
func numintmvEvenUp(n int) int {
	if n < 2 {
		return 2
	}
	if n%2 == 1 {
		return n + 1
	}
	return n
}

// numintmvMul4Up returns n rounded up to a positive multiple of four.
func numintmvMul4Up(n int) int {
	if n < 4 {
		return 4
	}
	if r := n % 4; r != 0 {
		return n + (4 - r)
	}
	return n
}

// numintmvTrapW builds nodes and weights for the composite trapezoidal rule
// with n subintervals on [a, b].
func numintmvTrapW(a, b float64, n int) (xs, ws []float64) {
	if n < 1 {
		n = 1
	}
	h := (b - a) / float64(n)
	xs = make([]float64, n+1)
	ws = make([]float64, n+1)
	for i := 0; i <= n; i++ {
		xs[i] = a + float64(i)*h
		ws[i] = h
	}
	ws[0] = h / 2
	ws[n] = h / 2
	return
}

// numintmvMidW builds nodes and weights for the composite midpoint rule with
// n subintervals on [a, b].
func numintmvMidW(a, b float64, n int) (xs, ws []float64) {
	if n < 1 {
		n = 1
	}
	h := (b - a) / float64(n)
	xs = make([]float64, n)
	ws = make([]float64, n)
	for i := 0; i < n; i++ {
		xs[i] = a + (float64(i)+0.5)*h
		ws[i] = h
	}
	return
}

// numintmvSimpW builds nodes and weights for the composite Simpson rule; n is
// rounded up to an even number of subintervals on [a, b].
func numintmvSimpW(a, b float64, n int) (xs, ws []float64) {
	n = numintmvEvenUp(n)
	h := (b - a) / float64(n)
	xs = make([]float64, n+1)
	ws = make([]float64, n+1)
	for i := 0; i <= n; i++ {
		xs[i] = a + float64(i)*h
	}
	ws[0] = h / 3
	ws[n] = h / 3
	for i := 1; i < n; i++ {
		if i%2 == 1 {
			ws[i] = 4 * h / 3
		} else {
			ws[i] = 2 * h / 3
		}
	}
	return
}

// numintmvBooleW builds nodes and weights for the composite Boole rule; n is
// rounded up to a multiple of four on [a, b].
func numintmvBooleW(a, b float64, n int) (xs, ws []float64) {
	n = numintmvMul4Up(n)
	h := (b - a) / float64(n)
	xs = make([]float64, n+1)
	ws = make([]float64, n+1)
	for i := 0; i <= n; i++ {
		xs[i] = a + float64(i)*h
	}
	c := 2 * h / 45
	for p := 0; p < n; p += 4 {
		ws[p] += 7 * c
		ws[p+1] += 32 * c
		ws[p+2] += 12 * c
		ws[p+3] += 32 * c
		ws[p+4] += 7 * c
	}
	return
}

// numintmvLegendre computes the n Gauss-Legendre nodes and weights on [-1, 1]
// using Newton iteration on the Legendre polynomial.
func numintmvLegendre(n int) (nodes, weights []float64) {
	nodes = make([]float64, n)
	weights = make([]float64, n)
	if n < 1 {
		return
	}
	if n == 1 {
		nodes[0] = 0
		weights[0] = 2
		return
	}
	m := (n + 1) / 2
	for i := 0; i < m; i++ {
		x := math.Cos(math.Pi * (float64(i) + 0.75) / (float64(n) + 0.5))
		var p0, p1 float64
		for iter := 0; iter < 100; iter++ {
			p0, p1 = 1.0, 0.0
			for j := 0; j < n; j++ {
				p2 := p1
				p1 = p0
				p0 = ((2*float64(j)+1)*x*p1 - float64(j)*p2) / float64(j+1)
			}
			dp := float64(n) * (x*p0 - p1) / (x*x - 1)
			dx := p0 / dp
			x -= dx
			if math.Abs(dx) < 1e-15 {
				break
			}
		}
		p0, p1 = 1.0, 0.0
		for j := 0; j < n; j++ {
			p2 := p1
			p1 = p0
			p0 = ((2*float64(j)+1)*x*p1 - float64(j)*p2) / float64(j+1)
		}
		dp := float64(n) * (x*p0 - p1) / (x*x - 1)
		nodes[i] = -x
		nodes[n-1-i] = x
		w := 2 / ((1 - x*x) * dp * dp)
		weights[i] = w
		weights[n-1-i] = w
	}
	return
}

// numintmvGLW maps the n Gauss-Legendre nodes and weights to the interval
// [a, b].
func numintmvGLW(a, b float64, n int) (xs, ws []float64) {
	nodes, weights := numintmvLegendre(n)
	xs = make([]float64, len(nodes))
	ws = make([]float64, len(nodes))
	c1 := 0.5 * (b - a)
	c2 := 0.5 * (a + b)
	for i := range nodes {
		xs[i] = c1*nodes[i] + c2
		ws[i] = c1 * weights[i]
	}
	return
}

// numintmvSimp1D integrates f over [a, b] with the composite Simpson rule.
func numintmvSimp1D(f Func, a, b float64, n int) float64 {
	xs, ws := numintmvSimpW(a, b, n)
	s := 0.0
	for i := range xs {
		s += ws[i] * f(xs[i])
	}
	return s
}

// numintmvTensor2 evaluates the tensor-product quadrature of f2 over the two
// node/weight sets.
func numintmvTensor2(f Func2, xs, wx, ys, wy []float64) float64 {
	total := 0.0
	for i := range xs {
		xi := xs[i]
		inner := 0.0
		for j := range ys {
			inner += wy[j] * f2eval(f, xi, ys[j])
		}
		total += wx[i] * inner
	}
	return total
}

// f2eval evaluates a Func2; kept separate so the tensor helper stays readable.
func f2eval(f Func2, x, y float64) float64 { return f(x, y) }

// numintmvTensor3 evaluates the tensor-product quadrature of f3 over the three
// node/weight sets.
func numintmvTensor3(f Func3, xs, wx, ys, wy, zs, wz []float64) float64 {
	total := 0.0
	for i := range xs {
		xi := xs[i]
		sy := 0.0
		for j := range ys {
			yj := ys[j]
			sz := 0.0
			for k := range zs {
				sz += wz[k] * f(xi, yj, zs[k])
			}
			sy += wy[j] * sz
		}
		total += wx[i] * sy
	}
	return total
}

// numintmvTqli diagonalizes a symmetric tridiagonal matrix using the implicit
// QL algorithm; d holds the diagonal (eigenvalues on return), e the
// subdiagonal, and z the first row of the accumulated eigenvector matrix.
func numintmvTqli(d, e, z []float64) {
	n := len(d)
	for i := 1; i < n; i++ {
		e[i-1] = e[i]
	}
	e[n-1] = 0
	for l := 0; l < n; l++ {
		iter := 0
		for {
			var m int
			for m = l; m < n-1; m++ {
				dd := math.Abs(d[m]) + math.Abs(d[m+1])
				if math.Abs(e[m])+dd == dd {
					break
				}
			}
			if m == l {
				break
			}
			iter++
			if iter > 60 {
				break
			}
			g := (d[l+1] - d[l]) / (2 * e[l])
			r := math.Hypot(g, 1)
			g = d[m] - d[l] + e[l]/(g+numintmvSign(r, g))
			s, c := 1.0, 1.0
			p := 0.0
			var i int
			broke := false
			for i = m - 1; i >= l; i-- {
				f := s * e[i]
				b := c * e[i]
				r = math.Hypot(f, g)
				e[i+1] = r
				if r == 0 {
					d[i+1] -= p
					e[m] = 0
					broke = true
					break
				}
				s = f / r
				c = g / r
				g = d[i+1] - p
				r = (d[i]-g)*s + 2*c*b
				p = s * r
				d[i+1] = g + p
				g = c*r - b
				fz := z[i+1]
				z[i+1] = s*z[i] + c*fz
				z[i] = c*z[i] - s*fz
			}
			if broke && r == 0 && i >= l {
				continue
			}
			d[l] -= p
			e[l] = g
			e[m] = 0
		}
	}
}

// numintmvGolub applies the Golub-Welsch algorithm to a Jacobi matrix with the
// given diagonal and subdiagonal and zeroth moment mu0, returning sorted nodes
// and weights of the associated Gaussian quadrature.
func numintmvGolub(diag, sub []float64, mu0 float64) (nodes, weights []float64) {
	n := len(diag)
	nodes = make([]float64, n)
	weights = make([]float64, n)
	if n == 0 {
		return
	}
	d := make([]float64, n)
	copy(d, diag)
	e := make([]float64, n)
	for i := 1; i < n; i++ {
		e[i] = sub[i-1]
	}
	z := make([]float64, n)
	z[0] = 1
	numintmvTqli(d, e, z)
	idx := make([]int, n)
	for i := range idx {
		idx[i] = i
	}
	sort.Slice(idx, func(a, b int) bool { return d[idx[a]] < d[idx[b]] })
	for k, i := range idx {
		nodes[k] = d[i]
		weights[k] = mu0 * z[i] * z[i]
	}
	return
}

// numintmvPathDeriv2 approximates the derivative of a plane curve at t via a
// central difference.
func numintmvPathDeriv2(p Path2, t float64) (dx, dy float64) {
	const h = 1e-6
	x1, y1 := p(t + h)
	x0, y0 := p(t - h)
	return (x1 - x0) / (2 * h), (y1 - y0) / (2 * h)
}

// numintmvPathDeriv3 approximates the derivative of a space curve at t via a
// central difference.
func numintmvPathDeriv3(p Path3, t float64) (dx, dy, dz float64) {
	const h = 1e-6
	x1, y1, z1 := p(t + h)
	x0, y0, z0 := p(t - h)
	return (x1 - x0) / (2 * h), (y1 - y0) / (2 * h), (z1 - z0) / (2 * h)
}

// numintmvPrimes returns the first k prime numbers.
func numintmvPrimes(k int) []int {
	if k < 1 {
		return nil
	}
	primes := make([]int, 0, k)
	cand := 2
	for len(primes) < k {
		isP := true
		for _, p := range primes {
			if p*p > cand {
				break
			}
			if cand%p == 0 {
				isP = false
				break
			}
		}
		if isP {
			primes = append(primes, cand)
		}
		cand++
	}
	return primes
}

// ============================================================================
// Double and triple integrals over rectangles
// ============================================================================

// DoubleTrapezoid approximates the double integral of f over the rectangle
// [ax, bx] x [ay, by] using the composite trapezoidal rule with nx and ny
// subintervals along each axis.
func DoubleTrapezoid(f Func2, ax, bx, ay, by float64, nx, ny int) float64 {
	xs, wx := numintmvTrapW(ax, bx, nx)
	ys, wy := numintmvTrapW(ay, by, ny)
	return numintmvTensor2(f, xs, wx, ys, wy)
}

// DoubleMidpoint approximates the double integral of f over the rectangle
// [ax, bx] x [ay, by] using the composite midpoint rule.
func DoubleMidpoint(f Func2, ax, bx, ay, by float64, nx, ny int) float64 {
	xs, wx := numintmvMidW(ax, bx, nx)
	ys, wy := numintmvMidW(ay, by, ny)
	return numintmvTensor2(f, xs, wx, ys, wy)
}

// DoubleSimpson approximates the double integral of f over the rectangle
// [ax, bx] x [ay, by] using the composite Simpson rule.
func DoubleSimpson(f Func2, ax, bx, ay, by float64, nx, ny int) float64 {
	xs, wx := numintmvSimpW(ax, bx, nx)
	ys, wy := numintmvSimpW(ay, by, ny)
	return numintmvTensor2(f, xs, wx, ys, wy)
}

// DoubleBoole approximates the double integral of f over the rectangle
// [ax, bx] x [ay, by] using the composite Boole rule.
func DoubleBoole(f Func2, ax, bx, ay, by float64, nx, ny int) float64 {
	xs, wx := numintmvBooleW(ax, bx, nx)
	ys, wy := numintmvBooleW(ay, by, ny)
	return numintmvTensor2(f, xs, wx, ys, wy)
}

// DoubleGaussLegendre approximates the double integral of f over the rectangle
// [ax, bx] x [ay, by] using an n-by-n tensor Gauss-Legendre rule.
func DoubleGaussLegendre(f Func2, ax, bx, ay, by float64, n int) float64 {
	xs, wx := numintmvGLW(ax, bx, n)
	ys, wy := numintmvGLW(ay, by, n)
	return numintmvTensor2(f, xs, wx, ys, wy)
}

// DoubleIntegral approximates the double integral of f over the rectangle
// [ax, bx] x [ay, by] using a fine composite Simpson rule with sensible
// default resolution.
func DoubleIntegral(f Func2, ax, bx, ay, by float64) float64 {
	return DoubleSimpson(f, ax, bx, ay, by, 64, 64)
}

// DoubleAverage returns the mean value of f over the rectangle
// [ax, bx] x [ay, by], i.e. its double integral divided by the area.
func DoubleAverage(f Func2, ax, bx, ay, by float64) float64 {
	area := (bx - ax) * (by - ay)
	if area == 0 {
		return 0
	}
	return DoubleIntegral(f, ax, bx, ay, by) / area
}

// TripleTrapezoid approximates the triple integral of f over the box
// [ax, bx] x [ay, by] x [az, bz] using the composite trapezoidal rule.
func TripleTrapezoid(f Func3, ax, bx, ay, by, az, bz float64, nx, ny, nz int) float64 {
	xs, wx := numintmvTrapW(ax, bx, nx)
	ys, wy := numintmvTrapW(ay, by, ny)
	zs, wz := numintmvTrapW(az, bz, nz)
	return numintmvTensor3(f, xs, wx, ys, wy, zs, wz)
}

// TripleMidpoint approximates the triple integral of f over the box using the
// composite midpoint rule.
func TripleMidpoint(f Func3, ax, bx, ay, by, az, bz float64, nx, ny, nz int) float64 {
	xs, wx := numintmvMidW(ax, bx, nx)
	ys, wy := numintmvMidW(ay, by, ny)
	zs, wz := numintmvMidW(az, bz, nz)
	return numintmvTensor3(f, xs, wx, ys, wy, zs, wz)
}

// TripleSimpson approximates the triple integral of f over the box using the
// composite Simpson rule.
func TripleSimpson(f Func3, ax, bx, ay, by, az, bz float64, nx, ny, nz int) float64 {
	xs, wx := numintmvSimpW(ax, bx, nx)
	ys, wy := numintmvSimpW(ay, by, ny)
	zs, wz := numintmvSimpW(az, bz, nz)
	return numintmvTensor3(f, xs, wx, ys, wy, zs, wz)
}

// TripleGaussLegendre approximates the triple integral of f over the box using
// an n-by-n-by-n tensor Gauss-Legendre rule.
func TripleGaussLegendre(f Func3, ax, bx, ay, by, az, bz float64, n int) float64 {
	xs, wx := numintmvGLW(ax, bx, n)
	ys, wy := numintmvGLW(ay, by, n)
	zs, wz := numintmvGLW(az, bz, n)
	return numintmvTensor3(f, xs, wx, ys, wy, zs, wz)
}

// TripleIntegral approximates the triple integral of f over the box using a
// composite Simpson rule with sensible default resolution.
func TripleIntegral(f Func3, ax, bx, ay, by, az, bz float64) float64 {
	return TripleSimpson(f, ax, bx, ay, by, az, bz, 24, 24, 24)
}

// TripleAverage returns the mean value of f over the box, i.e. its triple
// integral divided by the volume.
func TripleAverage(f Func3, ax, bx, ay, by, az, bz float64) float64 {
	vol := (bx - ax) * (by - ay) * (bz - az)
	if vol == 0 {
		return 0
	}
	return TripleIntegral(f, ax, bx, ay, by, az, bz) / vol
}

// ============================================================================
// Integration over general two-dimensional regions
// ============================================================================

// DoubleIntegralRegion approximates the integral of f over the region
// { (x, y) : a <= x <= b, ylo(x) <= y <= yhi(x) } using the composite Simpson
// rule in the outer variable and an ny-panel inner Simpson rule.
func DoubleIntegralRegion(f Func2, a, b float64, ylo, yhi Func, nx, ny int) float64 {
	xs, wx := numintmvSimpW(a, b, nx)
	total := 0.0
	for i := range xs {
		xi := xs[i]
		g := func(y float64) float64 { return f(xi, y) }
		total += wx[i] * numintmvSimp1D(g, ylo(xi), yhi(xi), ny)
	}
	return total
}

// DoubleIntegralRegionGL approximates the integral of f over the region
// { (x, y) : a <= x <= b, ylo(x) <= y <= yhi(x) } using an n-point
// Gauss-Legendre rule in each direction.
func DoubleIntegralRegionGL(f Func2, a, b float64, ylo, yhi Func, n int) float64 {
	xs, wx := numintmvGLW(a, b, n)
	total := 0.0
	for i := range xs {
		xi := xs[i]
		ys, wy := numintmvGLW(ylo(xi), yhi(xi), n)
		inner := 0.0
		for j := range ys {
			inner += wy[j] * f(xi, ys[j])
		}
		total += wx[i] * inner
	}
	return total
}

// AreaBetween returns the area of the region enclosed between the curves
// y = ylo(x) and y = yhi(x) for x in [a, b].
func AreaBetween(ylo, yhi Func, a, b float64, n int) float64 {
	g := func(x float64) float64 { return yhi(x) - ylo(x) }
	return numintmvSimp1D(g, a, b, n)
}

// CentroidRegion returns the centroid (xbar, ybar) of the region between the
// curves y = ylo(x) and y = yhi(x) over [a, b].
func CentroidRegion(ylo, yhi Func, a, b float64, n int) (xbar, ybar float64) {
	area := AreaBetween(ylo, yhi, a, b, n)
	if area == 0 {
		return 0, 0
	}
	mx := func(x float64) float64 { return x * (yhi(x) - ylo(x)) }
	my := func(x float64) float64 { return 0.5 * (yhi(x)*yhi(x) - ylo(x)*ylo(x)) }
	xbar = numintmvSimp1D(mx, a, b, n) / area
	ybar = numintmvSimp1D(my, a, b, n) / area
	return
}

// IntegrateTriangle approximates the integral of f over the triangle with the
// given vertices using an n-point Gauss-Legendre rule under a Duffy transform.
func IntegrateTriangle(f Func2, x1, y1, x2, y2, x3, y3 float64, n int) float64 {
	jac := math.Abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
	us, wu := numintmvGLW(0, 1, n)
	vs, wv := numintmvGLW(0, 1, n)
	total := 0.0
	for i := range us {
		u := us[i]
		for j := range vs {
			v := vs[j]
			s := u
			t := v * (1 - u)
			x := x1 + s*(x2-x1) + t*(x3-x1)
			y := y1 + s*(y2-y1) + t*(y3-y1)
			total += wu[i] * wv[j] * f(x, y) * (1 - u)
		}
	}
	return jac * total
}

// IntegrateDisk approximates the integral of f over the disk of radius r
// centered at (cx, cy), using polar coordinates with nr radial and nt angular
// subintervals.
func IntegrateDisk(f Func2, cx, cy, r float64, nr, nt int) float64 {
	rs, wr := numintmvSimpW(0, r, nr)
	ts, wt := numintmvSimpW(0, 2*math.Pi, nt)
	total := 0.0
	for i := range rs {
		rho := rs[i]
		inner := 0.0
		for j := range ts {
			th := ts[j]
			inner += wt[j] * f(cx+rho*math.Cos(th), cy+rho*math.Sin(th))
		}
		total += wr[i] * rho * inner
	}
	return total
}

// IntegrateAnnulus approximates the integral of f over the annulus centered at
// (cx, cy) with inner radius rin and outer radius rout.
func IntegrateAnnulus(f Func2, cx, cy, rin, rout float64, nr, nt int) float64 {
	rs, wr := numintmvSimpW(rin, rout, nr)
	ts, wt := numintmvSimpW(0, 2*math.Pi, nt)
	total := 0.0
	for i := range rs {
		rho := rs[i]
		inner := 0.0
		for j := range ts {
			th := ts[j]
			inner += wt[j] * f(cx+rho*math.Cos(th), cy+rho*math.Sin(th))
		}
		total += wr[i] * rho * inner
	}
	return total
}

// IntegratePolar approximates the integral of g(r, theta) * r over the polar
// rectangle r in [r0, r1], theta in [t0, t1]; the Jacobian factor r is applied
// automatically so g is the integrand expressed in polar coordinates.
func IntegratePolar(g Func2, r0, r1, t0, t1 float64, nr, nt int) float64 {
	rs, wr := numintmvSimpW(r0, r1, nr)
	ts, wt := numintmvSimpW(t0, t1, nt)
	total := 0.0
	for i := range rs {
		rho := rs[i]
		inner := 0.0
		for j := range ts {
			inner += wt[j] * g(rho, ts[j])
		}
		total += wr[i] * rho * inner
	}
	return total
}

// ============================================================================
// Integration in spherical and cylindrical coordinates
// ============================================================================

// IntegrateSphereSurface approximates the surface integral of f over the
// sphere of radius r centered at (cx, cy, cz), with nphi polar and ntheta
// azimuthal subintervals.
func IntegrateSphereSurface(f Func3, cx, cy, cz, r float64, nphi, ntheta int) float64 {
	phis, wp := numintmvSimpW(0, math.Pi, nphi)
	ths, wt := numintmvSimpW(0, 2*math.Pi, ntheta)
	total := 0.0
	for i := range phis {
		phi := phis[i]
		sp := math.Sin(phi)
		cpz := math.Cos(phi)
		inner := 0.0
		for j := range ths {
			th := ths[j]
			x := cx + r*sp*math.Cos(th)
			y := cy + r*sp*math.Sin(th)
			z := cz + r*cpz
			inner += wt[j] * f(x, y, z)
		}
		total += wp[i] * r * r * sp * inner
	}
	return total
}

// IntegrateBall approximates the volume integral of f over the solid ball of
// radius r centered at (cx, cy, cz) using spherical coordinates.
func IntegrateBall(f Func3, cx, cy, cz, r float64, nr, nphi, ntheta int) float64 {
	rs, wr := numintmvSimpW(0, r, nr)
	phis, wp := numintmvSimpW(0, math.Pi, nphi)
	ths, wt := numintmvSimpW(0, 2*math.Pi, ntheta)
	total := 0.0
	for a := range rs {
		rho := rs[a]
		sr := 0.0
		for b := range phis {
			phi := phis[b]
			sp := math.Sin(phi)
			cpz := math.Cos(phi)
			st := 0.0
			for c := range ths {
				th := ths[c]
				x := cx + rho*sp*math.Cos(th)
				y := cy + rho*sp*math.Sin(th)
				z := cz + rho*cpz
				st += wt[c] * f(x, y, z)
			}
			sr += wp[b] * sp * st
		}
		total += wr[a] * rho * rho * sr
	}
	return total
}

// IntegrateCylinder approximates the volume integral of f over the cylinder of
// radius r centered on the axis x = cx, y = cy and spanning z in [z0, z1],
// using cylindrical coordinates.
func IntegrateCylinder(f Func3, cx, cy, r, z0, z1 float64, nr, ntheta, nz int) float64 {
	rs, wr := numintmvSimpW(0, r, nr)
	ths, wt := numintmvSimpW(0, 2*math.Pi, ntheta)
	zs, wz := numintmvSimpW(z0, z1, nz)
	total := 0.0
	for a := range rs {
		rho := rs[a]
		sth := 0.0
		for b := range ths {
			th := ths[b]
			x := cx + rho*math.Cos(th)
			y := cy + rho*math.Sin(th)
			sz := 0.0
			for c := range zs {
				sz += wz[c] * f(x, y, zs[c])
			}
			sth += wt[b] * sz
		}
		total += wr[a] * rho * sth
	}
	return total
}

// ============================================================================
// Monte-Carlo integration (seeded, deterministic)
// ============================================================================

// MonteCarlo3D estimates the triple integral of f over the box
// [ax, bx] x [ay, by] x [az, bz] by averaging n uniform random samples drawn
// from the given seed.
func MonteCarlo3D(f Func3, ax, bx, ay, by, az, bz float64, n int, seed int64) float64 {
	if n < 1 {
		return 0
	}
	rng := rand.New(rand.NewSource(seed))
	vol := (bx - ax) * (by - ay) * (bz - az)
	sum := 0.0
	for i := 0; i < n; i++ {
		x := ax + (bx-ax)*rng.Float64()
		y := ay + (by-ay)*rng.Float64()
		z := az + (bz-az)*rng.Float64()
		sum += f(x, y, z)
	}
	return vol * sum / float64(n)
}

// MonteCarloWithError estimates the integral of f over [a, b] together with a
// standard-error estimate and the evaluation count.
func MonteCarloWithError(f Func, a, b float64, n int, seed int64) QuadResult {
	if n < 1 {
		return QuadResult{}
	}
	rng := rand.New(rand.NewSource(seed))
	vol := b - a
	sum, sum2 := 0.0, 0.0
	for i := 0; i < n; i++ {
		v := f(a + vol*rng.Float64())
		sum += v
		sum2 += v * v
	}
	mean := sum / float64(n)
	variance := sum2/float64(n) - mean*mean
	if variance < 0 {
		variance = 0
	}
	return QuadResult{
		Value:  vol * mean,
		ErrEst: vol * math.Sqrt(variance/float64(n)),
		Evals:  n,
	}
}

// MonteCarlo2DWithError estimates the double integral of f over the rectangle
// [ax, bx] x [ay, by] together with a standard-error estimate.
func MonteCarlo2DWithError(f Func2, ax, bx, ay, by float64, n int, seed int64) QuadResult {
	if n < 1 {
		return QuadResult{}
	}
	rng := rand.New(rand.NewSource(seed))
	area := (bx - ax) * (by - ay)
	sum, sum2 := 0.0, 0.0
	for i := 0; i < n; i++ {
		x := ax + (bx-ax)*rng.Float64()
		y := ay + (by-ay)*rng.Float64()
		v := f(x, y)
		sum += v
		sum2 += v * v
	}
	mean := sum / float64(n)
	variance := sum2/float64(n) - mean*mean
	if variance < 0 {
		variance = 0
	}
	return QuadResult{
		Value:  area * mean,
		ErrEst: area * math.Sqrt(variance/float64(n)),
		Evals:  n,
	}
}

// MonteCarloAntithetic estimates the integral of f over [a, b] using n
// antithetic sample pairs, which reduces variance for monotone integrands.
func MonteCarloAntithetic(f Func, a, b float64, n int, seed int64) float64 {
	if n < 1 {
		return 0
	}
	rng := rand.New(rand.NewSource(seed))
	vol := b - a
	sum := 0.0
	for i := 0; i < n; i++ {
		u := rng.Float64()
		sum += f(a+vol*u) + f(a+vol*(1-u))
	}
	return vol * sum / float64(2*n)
}

// MonteCarloStratified estimates the integral of f over [a, b] using stratified
// sampling with the given number of equal strata and samples per stratum.
func MonteCarloStratified(f Func, a, b float64, strata, perStratum int, seed int64) float64 {
	if strata < 1 || perStratum < 1 {
		return 0
	}
	rng := rand.New(rand.NewSource(seed))
	h := (b - a) / float64(strata)
	total := 0.0
	for k := 0; k < strata; k++ {
		lo := a + float64(k)*h
		s := 0.0
		for j := 0; j < perStratum; j++ {
			s += f(lo + h*rng.Float64())
		}
		total += h * s / float64(perStratum)
	}
	return total
}

// MonteCarloImportance estimates the integral of f using importance sampling:
// sample draws a point from a density proportional to pdf, and the estimator
// averages f(x)/pdf(x) over n samples.
func MonteCarloImportance(f Func, sample func(rng *rand.Rand) float64, pdf Func, n int, seed int64) float64 {
	if n < 1 {
		return 0
	}
	rng := rand.New(rand.NewSource(seed))
	sum := 0.0
	for i := 0; i < n; i++ {
		x := sample(rng)
		p := pdf(x)
		if p > 0 {
			sum += f(x) / p
		}
	}
	return sum / float64(n)
}

// MonteCarloDisk estimates the integral of f over the disk of radius r centered
// at (cx, cy) by averaging n samples drawn uniformly over the disk.
func MonteCarloDisk(f Func2, cx, cy, r float64, n int, seed int64) float64 {
	if n < 1 {
		return 0
	}
	rng := rand.New(rand.NewSource(seed))
	area := math.Pi * r * r
	sum := 0.0
	for i := 0; i < n; i++ {
		rho := r * math.Sqrt(rng.Float64())
		th := 2 * math.Pi * rng.Float64()
		sum += f(cx+rho*math.Cos(th), cy+rho*math.Sin(th))
	}
	return area * sum / float64(n)
}

// MonteCarloBall estimates the integral of f over the solid ball of radius r
// centered at (cx, cy, cz) by averaging n uniform samples over the ball.
func MonteCarloBall(f Func3, cx, cy, cz, r float64, n int, seed int64) float64 {
	if n < 1 {
		return 0
	}
	rng := rand.New(rand.NewSource(seed))
	vol := 4.0 / 3.0 * math.Pi * r * r * r
	sum := 0.0
	for i := 0; i < n; i++ {
		rho := r * math.Cbrt(rng.Float64())
		u := 2*rng.Float64() - 1
		phi := 2 * math.Pi * rng.Float64()
		sp := math.Sqrt(1 - u*u)
		x := cx + rho*sp*math.Cos(phi)
		y := cy + rho*sp*math.Sin(phi)
		z := cz + rho*u
		sum += f(x, y, z)
	}
	return vol * sum / float64(n)
}

// MonteCarloVolume estimates the volume of the subset of the box [lower, upper]
// for which pred reports true, by sampling n uniform random points.
func MonteCarloVolume(pred func(x []float64) bool, lower, upper []float64, n int, seed int64) float64 {
	d := len(lower)
	if d == 0 || len(upper) != d || n < 1 {
		return 0
	}
	rng := rand.New(rand.NewSource(seed))
	boxVol := 1.0
	for i := 0; i < d; i++ {
		boxVol *= upper[i] - lower[i]
	}
	point := make([]float64, d)
	hits := 0
	for i := 0; i < n; i++ {
		for k := 0; k < d; k++ {
			point[k] = lower[k] + (upper[k]-lower[k])*rng.Float64()
		}
		if pred(point) {
			hits++
		}
	}
	return boxVol * float64(hits) / float64(n)
}

// MonteCarloRegion estimates the integral of f over the region
// { (x, y) : a <= x <= b, ylo(x) <= y <= yhi(x) } using n random samples.
func MonteCarloRegion(f Func2, a, b float64, ylo, yhi Func, n int, seed int64) float64 {
	if n < 1 {
		return 0
	}
	rng := rand.New(rand.NewSource(seed))
	span := b - a
	sum := 0.0
	for i := 0; i < n; i++ {
		x := a + span*rng.Float64()
		lo := ylo(x)
		hi := yhi(x)
		y := lo + (hi-lo)*rng.Float64()
		sum += (hi - lo) * f(x, y)
	}
	return span * sum / float64(n)
}

// ============================================================================
// Quasi-Monte-Carlo integration via the Halton sequence
// ============================================================================

// VanDerCorput returns the radical-inverse (van der Corput) value of index in
// the given base, a number in [0, 1).
func VanDerCorput(index, base int) float64 {
	if base < 2 {
		base = 2
	}
	result := 0.0
	f := 1.0 / float64(base)
	i := index
	for i > 0 {
		result += f * float64(i%base)
		i /= base
		f /= float64(base)
	}
	return result
}

// HaltonPoint returns the Halton point at the given index using one base per
// coordinate.
func HaltonPoint(index int, bases []int) []float64 {
	p := make([]float64, len(bases))
	for i, b := range bases {
		p[i] = VanDerCorput(index, b)
	}
	return p
}

// HaltonSequence returns count Halton points (indices 1..count) using the given
// per-coordinate bases.
func HaltonSequence(count int, bases []int) [][]float64 {
	pts := make([][]float64, count)
	for i := 0; i < count; i++ {
		pts[i] = HaltonPoint(i+1, bases)
	}
	return pts
}

// HaltonPrimes returns the first k prime numbers, suitable as Halton bases for
// a k-dimensional sequence.
func HaltonPrimes(k int) []int {
	return numintmvPrimes(k)
}

// QuasiMonteCarlo estimates the integral of f over [a, b] using n points of the
// base-2 Halton (van der Corput) sequence.
func QuasiMonteCarlo(f Func, a, b float64, n int) float64 {
	if n < 1 {
		return 0
	}
	span := b - a
	sum := 0.0
	for i := 1; i <= n; i++ {
		sum += f(a + span*VanDerCorput(i, 2))
	}
	return span * sum / float64(n)
}

// QuasiMonteCarlo2D estimates the double integral of f over the rectangle
// [ax, bx] x [ay, by] using n points of the Halton sequence with bases 2 and 3.
func QuasiMonteCarlo2D(f Func2, ax, bx, ay, by float64, n int) float64 {
	if n < 1 {
		return 0
	}
	area := (bx - ax) * (by - ay)
	sum := 0.0
	for i := 1; i <= n; i++ {
		x := ax + (bx-ax)*VanDerCorput(i, 2)
		y := ay + (by-ay)*VanDerCorput(i, 3)
		sum += f(x, y)
	}
	return area * sum / float64(n)
}

// QuasiMonteCarloND estimates the integral of f over the box [lower, upper]
// using n Halton points whose bases are the first len(lower) primes.
func QuasiMonteCarloND(f FuncND, lower, upper []float64, n int) float64 {
	d := len(lower)
	if d == 0 || len(upper) != d || n < 1 {
		return 0
	}
	bases := numintmvPrimes(d)
	vol := 1.0
	for i := 0; i < d; i++ {
		vol *= upper[i] - lower[i]
	}
	point := make([]float64, d)
	sum := 0.0
	for i := 1; i <= n; i++ {
		for k := 0; k < d; k++ {
			point[k] = lower[k] + (upper[k]-lower[k])*VanDerCorput(i, bases[k])
		}
		sum += f(point)
	}
	return vol * sum / float64(n)
}

// ============================================================================
// Gauss-Hermite quadrature on (-inf, inf)
// ============================================================================

// GaussHermiteNodes returns the n nodes and weights of Gauss-Hermite
// quadrature for the weight exp(-x^2) on (-inf, inf).
func GaussHermiteNodes(n int) (nodes, weights []float64) {
	if n < 1 {
		return nil, nil
	}
	diag := make([]float64, n)
	sub := make([]float64, n-1)
	for k := 1; k < n; k++ {
		sub[k-1] = math.Sqrt(float64(k) / 2)
	}
	return numintmvGolub(diag, sub, math.Sqrt(math.Pi))
}

// GaussHermite approximates the integral of exp(-x^2) * f(x) over
// (-inf, inf) using n-point Gauss-Hermite quadrature.
func GaussHermite(f Func, n int) float64 {
	nodes, weights := GaussHermiteNodes(n)
	s := 0.0
	for i := range nodes {
		s += weights[i] * f(nodes[i])
	}
	return s
}

// ExpectationNormal approximates E[f(X)] for a normal random variable X with
// the given mean and standard deviation, using n-point Gauss-Hermite
// quadrature.
func ExpectationNormal(f Func, mean, stddev float64, n int) float64 {
	nodes, weights := GaussHermiteNodes(n)
	scale := math.Sqrt2 * stddev
	s := 0.0
	for i := range nodes {
		s += weights[i] * f(mean+scale*nodes[i])
	}
	return s / math.Sqrt(math.Pi)
}

// IntegrateGaussianWeighted approximates the integral of exp(-x^2/2) * f(x)
// over (-inf, inf) using n-point Gauss-Hermite quadrature.
func IntegrateGaussianWeighted(f Func, n int) float64 {
	nodes, weights := GaussHermiteNodes(n)
	s := 0.0
	for i := range nodes {
		s += weights[i] * f(math.Sqrt2*nodes[i])
	}
	return math.Sqrt2 * s
}

// ============================================================================
// Gauss-Laguerre quadrature on (0, inf)
// ============================================================================

// GaussLaguerreNodes returns the n nodes and weights of Gauss-Laguerre
// quadrature for the weight exp(-x) on (0, inf).
func GaussLaguerreNodes(n int) (nodes, weights []float64) {
	return GaussLaguerreGenNodes(n, 0)
}

// GaussLaguerreGenNodes returns the n nodes and weights of generalized
// Gauss-Laguerre quadrature for the weight x^alpha * exp(-x) on (0, inf).
func GaussLaguerreGenNodes(n int, alpha float64) (nodes, weights []float64) {
	if n < 1 {
		return nil, nil
	}
	diag := make([]float64, n)
	sub := make([]float64, n-1)
	for k := 0; k < n; k++ {
		diag[k] = 2*float64(k) + alpha + 1
	}
	for k := 1; k < n; k++ {
		sub[k-1] = math.Sqrt(float64(k) * (float64(k) + alpha))
	}
	return numintmvGolub(diag, sub, math.Gamma(alpha+1))
}

// GaussLaguerre approximates the integral of exp(-x) * f(x) over (0, inf)
// using n-point Gauss-Laguerre quadrature.
func GaussLaguerre(f Func, n int) float64 {
	nodes, weights := GaussLaguerreNodes(n)
	s := 0.0
	for i := range nodes {
		s += weights[i] * f(nodes[i])
	}
	return s
}

// GaussLaguerreGen approximates the integral of x^alpha * exp(-x) * f(x) over
// (0, inf) using n-point generalized Gauss-Laguerre quadrature.
func GaussLaguerreGen(f Func, alpha float64, n int) float64 {
	nodes, weights := GaussLaguerreGenNodes(n, alpha)
	s := 0.0
	for i := range nodes {
		s += weights[i] * f(nodes[i])
	}
	return s
}

// IntegrateExpWeighted approximates the integral of exp(-rate*x) * f(x) over
// (0, inf) using n-point Gauss-Laguerre quadrature; rate must be positive.
func IntegrateExpWeighted(f Func, rate float64, n int) float64 {
	if rate <= 0 {
		return math.NaN()
	}
	nodes, weights := GaussLaguerreNodes(n)
	s := 0.0
	for i := range nodes {
		s += weights[i] * f(nodes[i]/rate)
	}
	return s / rate
}

// LaplaceTransform approximates the Laplace transform of f evaluated at s,
// that is the integral of exp(-s*t) * f(t) over (0, inf), using n-point
// Gauss-Laguerre quadrature; s must be positive.
func LaplaceTransform(f Func, s float64, n int) float64 {
	return IntegrateExpWeighted(f, s, n)
}

// ============================================================================
// Line integrals
// ============================================================================

// ArcLength2D approximates the arc length of the plane curve p over the
// parameter interval [a, b] using n Simpson panels.
func ArcLength2D(p Path2, a, b float64, n int) float64 {
	g := func(t float64) float64 {
		dx, dy := numintmvPathDeriv2(p, t)
		return math.Hypot(dx, dy)
	}
	return numintmvSimp1D(g, a, b, n)
}

// ArcLength3D approximates the arc length of the space curve p over the
// parameter interval [a, b] using n Simpson panels.
func ArcLength3D(p Path3, a, b float64, n int) float64 {
	g := func(t float64) float64 {
		dx, dy, dz := numintmvPathDeriv3(p, t)
		return math.Sqrt(dx*dx + dy*dy + dz*dz)
	}
	return numintmvSimp1D(g, a, b, n)
}

// LineIntegralScalar2D approximates the scalar line integral of f along the
// plane curve p over [a, b], i.e. the integral of f ds.
func LineIntegralScalar2D(f Func2, p Path2, a, b float64, n int) float64 {
	g := func(t float64) float64 {
		x, y := p(t)
		dx, dy := numintmvPathDeriv2(p, t)
		return f(x, y) * math.Hypot(dx, dy)
	}
	return numintmvSimp1D(g, a, b, n)
}

// LineIntegralScalar3D approximates the scalar line integral of f along the
// space curve p over [a, b], i.e. the integral of f ds.
func LineIntegralScalar3D(f Func3, p Path3, a, b float64, n int) float64 {
	g := func(t float64) float64 {
		x, y, z := p(t)
		dx, dy, dz := numintmvPathDeriv3(p, t)
		return f(x, y, z) * math.Sqrt(dx*dx+dy*dy+dz*dz)
	}
	return numintmvSimp1D(g, a, b, n)
}

// LineIntegralVector2D approximates the line integral of the vector field along
// the plane curve p over [a, b], i.e. the work integral of F dot dr.
func LineIntegralVector2D(field Field2, p Path2, a, b float64, n int) float64 {
	g := func(t float64) float64 {
		x, y := p(t)
		u, v := field(x, y)
		dx, dy := numintmvPathDeriv2(p, t)
		return u*dx + v*dy
	}
	return numintmvSimp1D(g, a, b, n)
}

// LineIntegralVector3D approximates the line integral of the vector field along
// the space curve p over [a, b], i.e. the work integral of F dot dr.
func LineIntegralVector3D(field Field3, p Path3, a, b float64, n int) float64 {
	g := func(t float64) float64 {
		x, y, z := p(t)
		u, v, w := field(x, y, z)
		dx, dy, dz := numintmvPathDeriv3(p, t)
		return u*dx + v*dy + w*dz
	}
	return numintmvSimp1D(g, a, b, n)
}

// Circulation2D approximates the circulation of the vector field around the
// closed plane curve p, i.e. the line integral of F dot dr over [a, b].
func Circulation2D(field Field2, p Path2, a, b float64, n int) float64 {
	return LineIntegralVector2D(field, p, a, b, n)
}

// Flux2D approximates the outward flux of the vector field across the plane
// curve p over [a, b], i.e. the integral of F dot n ds for a counterclockwise
// parametrization.
func Flux2D(field Field2, p Path2, a, b float64, n int) float64 {
	g := func(t float64) float64 {
		x, y := p(t)
		u, v := field(x, y)
		dx, dy := numintmvPathDeriv2(p, t)
		return u*dy - v*dx
	}
	return numintmvSimp1D(g, a, b, n)
}

// LineIntegralPolyline2D approximates the scalar line integral of f along the
// polyline through the given vertices, using the midpoint of each segment.
func LineIntegralPolyline2D(f Func2, xs, ys []float64) float64 {
	n := len(xs)
	if n < 2 || len(ys) != n {
		return 0
	}
	total := 0.0
	for i := 0; i+1 < n; i++ {
		dx := xs[i+1] - xs[i]
		dy := ys[i+1] - ys[i]
		length := math.Hypot(dx, dy)
		mx := 0.5 * (xs[i] + xs[i+1])
		my := 0.5 * (ys[i] + ys[i+1])
		total += f(mx, my) * length
	}
	return total
}

// ArcLengthPolyline2D returns the total length of the polyline through the
// given vertices.
func ArcLengthPolyline2D(xs, ys []float64) float64 {
	n := len(xs)
	if n < 2 || len(ys) != n {
		return 0
	}
	total := 0.0
	for i := 0; i+1 < n; i++ {
		total += math.Hypot(xs[i+1]-xs[i], ys[i+1]-ys[i])
	}
	return total
}

// ============================================================================
// Improper integrals via substitution
// ============================================================================

// IntegrateInfinite approximates the integral of f over the whole real line
// using the substitution x = t/(1-t^2) with n-point Gauss-Legendre quadrature.
func IntegrateInfinite(f Func, n int) float64 {
	nodes, weights := numintmvLegendre(n)
	s := 0.0
	for i := range nodes {
		t := nodes[i]
		d := 1 - t*t
		x := t / d
		jac := (1 + t*t) / (d * d)
		s += weights[i] * f(x) * jac
	}
	return s
}

// IntegrateUpperInfinite approximates the integral of f over [a, inf) using the
// substitution x = a + t/(1-t) with n-point Gauss-Legendre quadrature.
func IntegrateUpperInfinite(f Func, a float64, n int) float64 {
	nodes, weights := numintmvLegendre(n)
	s := 0.0
	for i := range nodes {
		u := 0.5 * (nodes[i] + 1)
		w := 0.5 * weights[i]
		d := 1 - u
		x := a + u/d
		jac := 1 / (d * d)
		s += w * f(x) * jac
	}
	return s
}

// IntegrateLowerInfinite approximates the integral of f over (-inf, b] using
// the substitution x = b - u/(1-u) with n-point Gauss-Legendre quadrature.
func IntegrateLowerInfinite(f Func, b float64, n int) float64 {
	g := func(u float64) float64 { return f(b - u) }
	return IntegrateUpperInfinite(g, 0, n)
}

// IntegrateSingularEndpoint approximates the integral of f over [a, b] when f
// has an integrable singularity at the lower endpoint a, using the
// substitution x = a + t^2 to tame the singularity together with an n-point
// Gauss-Legendre rule whose nodes never touch the endpoint.
func IntegrateSingularEndpoint(f Func, a, b float64, n int) float64 {
	if b < a {
		return -IntegrateSingularEndpoint(f, b, a, n)
	}
	tb := math.Sqrt(b - a)
	ts, ws := numintmvGLW(0, tb, n)
	s := 0.0
	for i := range ts {
		t := ts[i]
		s += ws[i] * f(a+t*t) * 2 * t
	}
	return s
}
