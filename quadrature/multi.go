package quadrature

import (
	"math"
	"math/rand"
)

// Interval represents the closed one-dimensional interval [A, B].
type Interval struct {
	A, B float64
}

// Length returns B - A.
func (iv Interval) Length() float64 { return iv.B - iv.A }

// Mid returns the midpoint of the interval.
func (iv Interval) Mid() float64 { return 0.5 * (iv.A + iv.B) }

// Contains reports whether x lies within the closed interval.
func (iv Interval) Contains(x float64) bool { return x >= iv.A && x <= iv.B }

// Box represents an axis-aligned rectangular region in n dimensions given by
// its lower and upper corners, which must have equal length.
type Box struct {
	Lower, Upper []float64
}

// NewBox constructs a Box from lower and upper corner slices, copying them. It
// panics if the lengths differ.
func NewBox(lower, upper []float64) Box {
	if len(lower) != len(upper) {
		panic("quadrature: NewBox corner length mismatch")
	}
	lo := make([]float64, len(lower))
	hi := make([]float64, len(upper))
	copy(lo, lower)
	copy(hi, upper)
	return Box{Lower: lo, Upper: hi}
}

// Dim reports the dimensionality of the box.
func (b Box) Dim() int { return len(b.Lower) }

// Volume returns the product of the side lengths of the box.
func (b Box) Volume() float64 {
	v := 1.0
	for i := range b.Lower {
		v *= b.Upper[i] - b.Lower[i]
	}
	return v
}

// Contains reports whether the point x lies within the closed box.
func (b Box) Contains(x []float64) bool {
	if len(x) != len(b.Lower) {
		return false
	}
	for i := range x {
		if x[i] < b.Lower[i] || x[i] > b.Upper[i] {
			return false
		}
	}
	return true
}

// IntegrateProduct2 applies the tensor product of two one-dimensional rules to
// a bivariate function, evaluating f at every (node_i, node_j) pair weighted
// by the product of the corresponding weights.
func IntegrateProduct2(f Func2, rx, ry Rule) float64 {
	var s float64
	for i, xi := range rx.Nodes {
		wi := rx.Weights[i]
		var inner float64
		for j, yj := range ry.Nodes {
			inner += ry.Weights[j] * f(xi, yj)
		}
		s += wi * inner
	}
	return s
}

// DoubleGaussLegendre approximates the double integral of f over the rectangle
// [ax, bx] x [ay, by] using the n-point Gauss-Legendre rule in each direction.
func DoubleGaussLegendre(f Func2, ax, bx, ay, by float64, n int) float64 {
	rx := GaussLegendreRule(n).Scale(ax, bx)
	ry := GaussLegendreRule(n).Scale(ay, by)
	return IntegrateProduct2(f, rx, ry)
}

// TripleGaussLegendre approximates the triple integral of f over the box
// [ax, bx] x [ay, by] x [az, bz] using the n-point Gauss-Legendre rule in each
// direction.
func TripleGaussLegendre(f Func3, ax, bx, ay, by, az, bz float64, n int) float64 {
	rx := GaussLegendreRule(n).Scale(ax, bx)
	ry := GaussLegendreRule(n).Scale(ay, by)
	rz := GaussLegendreRule(n).Scale(az, bz)
	var s float64
	for i, xi := range rx.Nodes {
		for j, yj := range ry.Nodes {
			wij := rx.Weights[i] * ry.Weights[j]
			for k, zk := range rz.Nodes {
				s += wij * rz.Weights[k] * f(xi, yj, zk)
			}
		}
	}
	return s
}

// DoubleTrapezoid approximates the double integral of f over the rectangle
// [ax, bx] x [ay, by] using the composite trapezoidal rule with nx and ny
// subintervals.
func DoubleTrapezoid(f Func2, ax, bx, ay, by float64, nx, ny int) float64 {
	return CompositeTrapezoid(func(x float64) float64 {
		return CompositeTrapezoid(func(y float64) float64 { return f(x, y) }, ay, by, ny)
	}, ax, bx, nx)
}

// DoubleMidpoint approximates the double integral of f over the rectangle
// using the composite midpoint rule with nx and ny subintervals.
func DoubleMidpoint(f Func2, ax, bx, ay, by float64, nx, ny int) float64 {
	return CompositeMidpoint(func(x float64) float64 {
		return CompositeMidpoint(func(y float64) float64 { return f(x, y) }, ay, by, ny)
	}, ax, bx, nx)
}

// DoubleSimpson approximates the double integral of f over the rectangle using
// the composite Simpson rule with nx and ny subintervals.
func DoubleSimpson(f Func2, ax, bx, ay, by float64, nx, ny int) float64 {
	return CompositeSimpson(func(x float64) float64 {
		return CompositeSimpson(func(y float64) float64 { return f(x, y) }, ay, by, ny)
	}, ax, bx, nx)
}

// DoubleBoole approximates the double integral of f over the rectangle using
// the composite Boole rule with nx and ny subintervals.
func DoubleBoole(f Func2, ax, bx, ay, by float64, nx, ny int) float64 {
	return CompositeBoole(func(x float64) float64 {
		return CompositeBoole(func(y float64) float64 { return f(x, y) }, ay, by, ny)
	}, ax, bx, nx)
}

// TripleSimpson approximates the triple integral of f over the box using the
// composite Simpson rule with nx, ny and nz subintervals.
func TripleSimpson(f Func3, ax, bx, ay, by, az, bz float64, nx, ny, nz int) float64 {
	return CompositeSimpson(func(x float64) float64 {
		return CompositeSimpson(func(y float64) float64 {
			return CompositeSimpson(func(z float64) float64 { return f(x, y, z) }, az, bz, nz)
		}, ay, by, ny)
	}, ax, bx, nx)
}

// IntegrateBoxGaussLegendre approximates the integral of f over the box given
// by lower and upper corners using an nPerDim-point tensor Gauss-Legendre rule
// in every dimension. The cost grows as nPerDim^dim, so it is intended for low
// dimensions.
func IntegrateBoxGaussLegendre(f FuncN, lower, upper []float64, nPerDim int) float64 {
	dim := len(lower)
	nodes, weights := GaussLegendre(nPerDim)
	// scaled nodes/weights per dimension
	sn := make([][]float64, dim)
	sw := make([][]float64, dim)
	for d := 0; d < dim; d++ {
		half := 0.5 * (upper[d] - lower[d])
		mid := 0.5 * (upper[d] + lower[d])
		sn[d] = make([]float64, nPerDim)
		sw[d] = make([]float64, nPerDim)
		for i := 0; i < nPerDim; i++ {
			sn[d][i] = mid + half*nodes[i]
			sw[d][i] = half * weights[i]
		}
	}
	pt := make([]float64, dim)
	idx := make([]int, dim)
	var total float64
	for {
		w := 1.0
		for d := 0; d < dim; d++ {
			pt[d] = sn[d][idx[d]]
			w *= sw[d][idx[d]]
		}
		total += w * f(pt)
		// increment odometer
		d := dim - 1
		for d >= 0 {
			idx[d]++
			if idx[d] < nPerDim {
				break
			}
			idx[d] = 0
			d--
		}
		if d < 0 {
			break
		}
	}
	return total
}

// MonteCarloBox estimates the integral of f over the box given by lower and
// upper corners using plain Monte-Carlo sampling with the given number of
// samples and random seed. The returned Result carries the estimate together
// with the standard error of the mean scaled by the box volume.
func MonteCarloBox(f FuncN, lower, upper []float64, samples int, seed int64) Result {
	dim := len(lower)
	rng := rand.New(rand.NewSource(seed))
	vol := 1.0
	for d := 0; d < dim; d++ {
		vol *= upper[d] - lower[d]
	}
	pt := make([]float64, dim)
	var mean, m2 float64
	for i := 0; i < samples; i++ {
		for d := 0; d < dim; d++ {
			pt[d] = lower[d] + rng.Float64()*(upper[d]-lower[d])
		}
		v := f(pt)
		// Welford online mean/variance
		delta := v - mean
		mean += delta / float64(i+1)
		m2 += delta * (v - mean)
	}
	value := vol * mean
	var stderr float64
	if samples > 1 {
		variance := m2 / float64(samples-1)
		stderr = vol * math.Sqrt(variance/float64(samples))
	}
	return Result{Value: value, AbsErr: stderr, Evals: samples, Success: true}
}

// MonteCarlo2 estimates the double integral of f over [ax, bx] x [ay, by] by
// plain Monte-Carlo sampling.
func MonteCarlo2(f Func2, ax, bx, ay, by float64, samples int, seed int64) Result {
	return MonteCarloBox(func(p []float64) float64 { return f(p[0], p[1]) },
		[]float64{ax, ay}, []float64{bx, by}, samples, seed)
}

// MonteCarlo3 estimates the triple integral of f over the box
// [ax, bx] x [ay, by] x [az, bz] by plain Monte-Carlo sampling.
func MonteCarlo3(f Func3, ax, bx, ay, by, az, bz float64, samples int, seed int64) Result {
	return MonteCarloBox(func(p []float64) float64 { return f(p[0], p[1], p[2]) },
		[]float64{ax, ay, az}, []float64{bx, by, bz}, samples, seed)
}

// StratifiedMonteCarloBox estimates the integral of f over the box using
// stratified sampling: each dimension is split into strata subintervals,
// forming strata^dim cells, and one uniform sample is drawn from each cell.
// Stratification reduces variance relative to plain Monte-Carlo for smooth
// integrands.
func StratifiedMonteCarloBox(f FuncN, lower, upper []float64, strata int, seed int64) Result {
	dim := len(lower)
	if strata < 1 {
		strata = 1
	}
	rng := rand.New(rand.NewSource(seed))
	vol := 1.0
	for d := 0; d < dim; d++ {
		vol *= upper[d] - lower[d]
	}
	cells := 1
	for d := 0; d < dim; d++ {
		cells *= strata
	}
	pt := make([]float64, dim)
	idx := make([]int, dim)
	var sum float64
	for c := 0; c < cells; c++ {
		for d := 0; d < dim; d++ {
			cellW := (upper[d] - lower[d]) / float64(strata)
			pt[d] = lower[d] + (float64(idx[d])+rng.Float64())*cellW
		}
		sum += f(pt)
		d := dim - 1
		for d >= 0 {
			idx[d]++
			if idx[d] < strata {
				break
			}
			idx[d] = 0
			d--
		}
	}
	value := vol * sum / float64(cells)
	return Result{Value: value, Evals: cells, Success: true}
}

// VanDerCorput returns the i-th element (i >= 0) of the van der Corput
// low-discrepancy sequence in the given base, a number in [0, 1) formed by
// reflecting the base-b digits of i about the radix point.
func VanDerCorput(i, base int) float64 {
	if base < 2 {
		base = 2
	}
	result := 0.0
	f := 1.0 / float64(base)
	n := i
	for n > 0 {
		result += float64(n%base) * f
		n /= base
		f /= float64(base)
	}
	return result
}

// Halton returns the i-th point (i >= 0) of the Halton low-discrepancy
// sequence in dim dimensions, using the first dim prime numbers as bases.
func Halton(i, dim int) []float64 {
	p := firstPrimes(dim)
	pt := make([]float64, dim)
	for d := 0; d < dim; d++ {
		pt[d] = VanDerCorput(i, p[d])
	}
	return pt
}

// HaltonSequence returns the first n points of the dim-dimensional Halton
// sequence, skipping the degenerate index 0.
func HaltonSequence(n, dim int) [][]float64 {
	out := make([][]float64, n)
	for i := 0; i < n; i++ {
		out[i] = Halton(i+1, dim)
	}
	return out
}

// QuasiMonteCarloBox estimates the integral of f over the box using the Halton
// quasi-random sequence. For smooth integrands the deterministic Halton points
// give an error that decreases faster than plain Monte-Carlo.
func QuasiMonteCarloBox(f FuncN, lower, upper []float64, samples int) Result {
	dim := len(lower)
	vol := 1.0
	for d := 0; d < dim; d++ {
		vol *= upper[d] - lower[d]
	}
	primes := firstPrimes(dim)
	pt := make([]float64, dim)
	var sum float64
	for i := 1; i <= samples; i++ {
		for d := 0; d < dim; d++ {
			u := VanDerCorput(i, primes[d])
			pt[d] = lower[d] + u*(upper[d]-lower[d])
		}
		sum += f(pt)
	}
	value := vol * sum / float64(samples)
	return Result{Value: value, Evals: samples, Success: true}
}

// firstPrimes returns the first k prime numbers.
func firstPrimes(k int) []int {
	if k <= 0 {
		return nil
	}
	primes := make([]int, 0, k)
	candidate := 2
	for len(primes) < k {
		isPrime := true
		for _, p := range primes {
			if p*p > candidate {
				break
			}
			if candidate%p == 0 {
				isPrime = false
				break
			}
		}
		if isPrime {
			primes = append(primes, candidate)
		}
		candidate++
	}
	return primes
}
