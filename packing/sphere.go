package packing

import "math"

// BallVolume returns the volume of an n-dimensional Euclidean ball of radius r,
// pi^(n/2) r^n / Gamma(n/2 + 1). For n = 0 the "ball" is a point of volume 1.
func BallVolume(n int, r float64) float64 {
	if n < 0 {
		return math.NaN()
	}
	return UnitBallVolume(n) * math.Pow(r, float64(n))
}

// UnitBallVolume returns the volume V_n of the unit ball in R^n,
// pi^(n/2) / Gamma(n/2 + 1). V_0 = 1, V_1 = 2, V_2 = pi, V_3 = 4pi/3.
func UnitBallVolume(n int) float64 {
	if n < 0 {
		return math.NaN()
	}
	return math.Pow(math.Pi, float64(n)/2) / math.Gamma(float64(n)/2+1)
}

// BallSurface returns the surface area of the (n-1)-sphere bounding an
// n-dimensional ball of radius r, n * V_n * r^(n-1).
func BallSurface(n int, r float64) float64 {
	if n < 1 {
		return math.NaN()
	}
	return float64(n) * UnitBallVolume(n) * math.Pow(r, float64(n-1))
}

// UnitSphereSurface returns the surface area of the unit (n-1)-sphere in R^n,
// 2 pi^(n/2) / Gamma(n/2). S_1 = 2 (two points), S_2 = 2pi, S_3 = 4pi.
func UnitSphereSurface(n int) float64 {
	if n < 1 {
		return math.NaN()
	}
	return 2 * math.Pow(math.Pi, float64(n)/2) / math.Gamma(float64(n)/2)
}

// PackingRadius returns the packing radius of a lattice whose minimal squared
// vector length (minimal norm) is minNorm, namely sqrt(minNorm)/2.
func PackingRadius(minNorm float64) float64 { return math.Sqrt(minNorm) / 2 }

// MinimalDistance returns the shortest distance between distinct lattice points
// given the minimal norm (minimal squared vector length), sqrt(minNorm).
func MinimalDistance(minNorm float64) float64 { return math.Sqrt(minNorm) }

// DensityToCenterDensity converts a packing density Delta back to a center
// density delta = Delta / V_n in dimension n.
func DensityToCenterDensity(density float64, n int) float64 {
	return density / UnitBallVolume(n)
}

// CenterDensity returns the center density delta = rho^n / V of a lattice with
// packing radius rho and covolume V in dimension n.
func CenterDensity(rho, covolume float64, n int) float64 {
	return math.Pow(rho, float64(n)) / covolume
}

// PackingDensity returns the sphere-packing density Delta = delta * V_n given a
// center density delta in dimension n.
func PackingDensity(centerDensity float64, n int) float64 {
	return centerDensity * UnitBallVolume(n)
}

// Thickness returns the covering thickness Theta = V_n * R^n / V of a lattice
// with covering radius R and covolume V in dimension n. A thickness of 1 would
// mean a perfect tiling by balls; every real covering has thickness >= 1.
func Thickness(coveringRadius, covolume float64, n int) float64 {
	return UnitBallVolume(n) * math.Pow(coveringRadius, float64(n)) / covolume
}

// ----------------------------------------------------------------------------
// The integer lattice Z^n.
// ----------------------------------------------------------------------------

// ZnCovolume returns the covolume of the integer lattice Z^n, which is 1.
func ZnCovolume(n int) float64 { return 1 }

// ZnMinimalNorm returns the minimal squared vector length of Z^n, which is 1.
func ZnMinimalNorm(n int) float64 { return 1 }

// ZnPackingRadius returns the packing radius of Z^n, which is 1/2.
func ZnPackingRadius(n int) float64 { return 0.5 }

// ZnKissingNumber returns the kissing number of Z^n, 2n (the +-1 unit vectors).
func ZnKissingNumber(n int) int { return 2 * n }

// ZnCenterDensity returns the center density of Z^n, 2^(-n).
func ZnCenterDensity(n int) float64 { return math.Pow(2, -float64(n)) }

// ZnDensity returns the sphere-packing density of Z^n, V_n / 2^n.
func ZnDensity(n int) float64 { return PackingDensity(ZnCenterDensity(n), n) }

// ZnCoveringRadius returns the covering radius of Z^n, sqrt(n)/2 (the deep hole
// is the all-halves center of a unit cube).
func ZnCoveringRadius(n int) float64 { return math.Sqrt(float64(n)) / 2 }

// ZnThickness returns the covering thickness of Z^n, V_n * (sqrt(n)/2)^n.
func ZnThickness(n int) float64 { return Thickness(ZnCoveringRadius(n), 1, n) }

// ----------------------------------------------------------------------------
// The root lattice A_n (n >= 1).
// ----------------------------------------------------------------------------

// AnCovolume returns the covolume of the root lattice A_n, sqrt(n+1).
func AnCovolume(n int) float64 { return math.Sqrt(float64(n + 1)) }

// AnMinimalNorm returns the minimal squared vector length of A_n, which is 2.
func AnMinimalNorm(n int) float64 { return 2 }

// AnPackingRadius returns the packing radius of A_n, 1/sqrt(2).
func AnPackingRadius(n int) float64 { return 1 / math.Sqrt2 }

// AnKissingNumber returns the kissing number of A_n, n(n+1).
func AnKissingNumber(n int) int { return n * (n + 1) }

// AnCenterDensity returns the center density of A_n, 2^(-n/2) / sqrt(n+1).
func AnCenterDensity(n int) float64 {
	return CenterDensity(AnPackingRadius(n), AnCovolume(n), n)
}

// AnDensity returns the sphere-packing density of A_n.
func AnDensity(n int) float64 { return PackingDensity(AnCenterDensity(n), n) }

// AnCoveringRadius returns the covering radius of A_n. Its deep holes are the
// glue vectors [a]; the covering radius squared is max over a in 1..n of
// a(n+1-a)/(n+1), attained at a = floor((n+1)/2).
func AnCoveringRadius(n int) float64 {
	a := (n + 1) / 2
	r2 := float64(a) * float64(n+1-a) / float64(n+1)
	return math.Sqrt(r2)
}

// AnThickness returns the covering thickness of A_n.
func AnThickness(n int) float64 {
	return Thickness(AnCoveringRadius(n), AnCovolume(n), n)
}

// ----------------------------------------------------------------------------
// The dual lattice A_n^* (best known lattice covering in low dimensions).
// ----------------------------------------------------------------------------

// AnStarCovolume returns the covolume of the dual lattice A_n^*, 1/sqrt(n+1).
func AnStarCovolume(n int) float64 { return 1 / math.Sqrt(float64(n+1)) }

// AnStarCoveringRadius returns the covering radius of A_n^*, whose square is
// n(n+2)/(12(n+1)). A_n^* is the thinnest known lattice covering for n <= 5.
func AnStarCoveringRadius(n int) float64 {
	r2 := float64(n) * float64(n+2) / (12 * float64(n+1))
	return math.Sqrt(r2)
}

// AnStarThickness returns the covering thickness of A_n^*.
func AnStarThickness(n int) float64 {
	return Thickness(AnStarCoveringRadius(n), AnStarCovolume(n), n)
}

// ----------------------------------------------------------------------------
// The root lattice D_n (checkerboard lattice, n >= 2; densest for n = 3,4,5).
// ----------------------------------------------------------------------------

// DnCovolume returns the covolume of the checkerboard lattice D_n, which is 2.
func DnCovolume(n int) float64 { return 2 }

// DnMinimalNorm returns the minimal squared vector length of D_n, which is 2.
func DnMinimalNorm(n int) float64 { return 2 }

// DnPackingRadius returns the packing radius of D_n, 1/sqrt(2).
func DnPackingRadius(n int) float64 { return 1 / math.Sqrt2 }

// DnKissingNumber returns the kissing number of D_n, 2n(n-1).
func DnKissingNumber(n int) int { return 2 * n * (n - 1) }

// DnCenterDensity returns the center density of D_n, 2^(-(n+2)/2).
func DnCenterDensity(n int) float64 {
	return CenterDensity(DnPackingRadius(n), DnCovolume(n), n)
}

// DnDensity returns the sphere-packing density of D_n.
func DnDensity(n int) float64 { return PackingDensity(DnCenterDensity(n), n) }

// DnCoveringRadius returns the covering radius of D_n. The deep holes are the
// glue vector [1] at distance 1 and the all-halves glue vector [2] at distance
// sqrt(n)/2, so the covering radius is max(1, sqrt(n)/2).
func DnCoveringRadius(n int) float64 {
	return math.Max(1, math.Sqrt(float64(n))/2)
}

// DnThickness returns the covering thickness of D_n.
func DnThickness(n int) float64 {
	return Thickness(DnCoveringRadius(n), DnCovolume(n), n)
}

// ----------------------------------------------------------------------------
// The exceptional root lattices E6, E7, E8.
// ----------------------------------------------------------------------------

// E6Dimension returns the dimension of the lattice E6, which is 6.
func E6Dimension() int { return 6 }

// E6Covolume returns the covolume of E6, sqrt(3).
func E6Covolume() float64 { return math.Sqrt(3) }

// E6MinimalNorm returns the minimal squared vector length of E6, which is 2.
func E6MinimalNorm() float64 { return 2 }

// E6PackingRadius returns the packing radius of E6, 1/sqrt(2).
func E6PackingRadius() float64 { return 1 / math.Sqrt2 }

// E6KissingNumber returns the kissing number of E6, which is 72.
func E6KissingNumber() int { return 72 }

// E6CenterDensity returns the center density of E6, 1/(8 sqrt(3)).
func E6CenterDensity() float64 { return CenterDensity(E6PackingRadius(), E6Covolume(), 6) }

// E6Density returns the sphere-packing density of E6.
func E6Density() float64 { return PackingDensity(E6CenterDensity(), 6) }

// E6CoveringRadius returns the covering radius of E6, whose square is 4/3.
func E6CoveringRadius() float64 { return math.Sqrt(4.0 / 3.0) }

// E6Thickness returns the covering thickness of E6.
func E6Thickness() float64 { return Thickness(E6CoveringRadius(), E6Covolume(), 6) }

// E7Dimension returns the dimension of the lattice E7, which is 7.
func E7Dimension() int { return 7 }

// E7Covolume returns the covolume of E7, sqrt(2).
func E7Covolume() float64 { return math.Sqrt2 }

// E7MinimalNorm returns the minimal squared vector length of E7, which is 2.
func E7MinimalNorm() float64 { return 2 }

// E7PackingRadius returns the packing radius of E7, 1/sqrt(2).
func E7PackingRadius() float64 { return 1 / math.Sqrt2 }

// E7KissingNumber returns the kissing number of E7, which is 126.
func E7KissingNumber() int { return 126 }

// E7CenterDensity returns the center density of E7, 1/16.
func E7CenterDensity() float64 { return CenterDensity(E7PackingRadius(), E7Covolume(), 7) }

// E7Density returns the sphere-packing density of E7.
func E7Density() float64 { return PackingDensity(E7CenterDensity(), 7) }

// E7CoveringRadius returns the covering radius of E7, whose square is 3/2.
func E7CoveringRadius() float64 { return math.Sqrt(1.5) }

// E7Thickness returns the covering thickness of E7.
func E7Thickness() float64 { return Thickness(E7CoveringRadius(), E7Covolume(), 7) }

// E8Dimension returns the dimension of the lattice E8, which is 8.
func E8Dimension() int { return 8 }

// E8Covolume returns the covolume of E8, which is 1 (E8 is unimodular).
func E8Covolume() float64 { return 1 }

// E8MinimalNorm returns the minimal squared vector length of E8, which is 2.
func E8MinimalNorm() float64 { return 2 }

// E8PackingRadius returns the packing radius of E8, 1/sqrt(2).
func E8PackingRadius() float64 { return 1 / math.Sqrt2 }

// E8KissingNumber returns the kissing number of E8, which is 240. This is
// optimal: no arrangement of unit spheres in R^8 touches a central sphere more
// than 240 times.
func E8KissingNumber() int { return 240 }

// E8CenterDensity returns the center density of E8, which is 1/16.
func E8CenterDensity() float64 { return CenterDensity(E8PackingRadius(), E8Covolume(), 8) }

// E8Density returns the sphere-packing density of E8, pi^4/384. Viazovska
// proved in 2016 that this is the densest packing of R^8.
func E8Density() float64 { return PackingDensity(E8CenterDensity(), 8) }

// E8CoveringRadius returns the covering radius of E8, which is 1.
func E8CoveringRadius() float64 { return 1 }

// E8Thickness returns the covering thickness of E8, V_8 = pi^4/24.
func E8Thickness() float64 { return Thickness(E8CoveringRadius(), E8Covolume(), 8) }

// ----------------------------------------------------------------------------
// The Leech lattice Lambda_24.
// ----------------------------------------------------------------------------

// LeechDimension returns the dimension of the Leech lattice, which is 24.
func LeechDimension() int { return 24 }

// LeechCovolume returns the covolume of the Leech lattice, which is 1 (it is
// unimodular).
func LeechCovolume() float64 { return 1 }

// LeechMinimalNorm returns the minimal squared vector length of the Leech
// lattice, which is 4.
func LeechMinimalNorm() float64 { return 4 }

// LeechPackingRadius returns the packing radius of the Leech lattice, which is
// 1 (minimal distance 2).
func LeechPackingRadius() float64 { return 1 }

// LeechKissingNumber returns the kissing number of the Leech lattice, 196560.
// This is optimal in dimension 24.
func LeechKissingNumber() int { return 196560 }

// LeechCenterDensity returns the center density of the Leech lattice, which is
// 1.
func LeechCenterDensity() float64 { return CenterDensity(LeechPackingRadius(), LeechCovolume(), 24) }

// LeechDensity returns the sphere-packing density of the Leech lattice, the
// volume of the unit 24-ball pi^12/12!. Cohn, Kumar, Miller, Radchenko and
// Viazovska proved in 2016 that this is the densest packing of R^24.
func LeechDensity() float64 { return PackingDensity(LeechCenterDensity(), 24) }

// LeechCoveringRadius returns the covering radius of the Leech lattice, sqrt(2)
// (its deep holes lie at distance sqrt(2) and come in 23 types).
func LeechCoveringRadius() float64 { return math.Sqrt2 }

// LeechThickness returns the covering thickness of the Leech lattice.
func LeechThickness() float64 { return Thickness(LeechCoveringRadius(), LeechCovolume(), 24) }

// ----------------------------------------------------------------------------
// Kissing numbers.
// ----------------------------------------------------------------------------

// KissingNumberOptimal returns the exact optimal kissing number tau(n) in
// dimension n for the dimensions in which it is known (1, 2, 3, 4, 8, 24), and
// -1 otherwise. tau(1)=2, tau(2)=6, tau(3)=12, tau(4)=24, tau(8)=240,
// tau(24)=196560.
func KissingNumberOptimal(n int) int {
	switch n {
	case 1:
		return 2
	case 2:
		return 6
	case 3:
		return 12
	case 4:
		return 24
	case 8:
		return 240
	case 24:
		return 196560
	default:
		return -1
	}
}

// KissingNumberKnown reports whether the optimal kissing number is known in
// dimension n.
func KissingNumberKnown(n int) bool { return KissingNumberOptimal(n) >= 0 }

// KissingNumberLowerBound returns the best kissing number attainable by the
// standard lattice family used in this package for dimension n: Z^n gives 2n,
// A_n gives n(n+1), D_n gives 2n(n-1) for n >= 3, plus the exceptional lattices
// in dimensions 8 and 24. It is a lower bound on the true optimum.
func KissingNumberLowerBound(n int) int {
	best := ZnKissingNumber(n)
	if k := AnKissingNumber(n); k > best {
		best = k
	}
	if n >= 3 {
		if k := DnKissingNumber(n); k > best {
			best = k
		}
	}
	switch n {
	case 6:
		if E6KissingNumber() > best {
			best = E6KissingNumber()
		}
	case 7:
		if E7KissingNumber() > best {
			best = E7KissingNumber()
		}
	case 8:
		if E8KissingNumber() > best {
			best = E8KissingNumber()
		}
	case 24:
		if LeechKissingNumber() > best {
			best = LeechKissingNumber()
		}
	}
	return best
}

// ----------------------------------------------------------------------------
// Densest known / proven packing densities.
// ----------------------------------------------------------------------------

// BestLatticeDensity returns the packing density of the densest lattice among
// the standard family (Z^n, A_n, D_n, E6, E7, E8, Leech) implemented here in
// dimension n. It is a lower bound on the optimal sphere-packing density and is
// in fact optimal for n in {1, 2, 3, 8, 24}.
func BestLatticeDensity(n int) float64 {
	best := ZnDensity(n)
	if d := AnDensity(n); d > best {
		best = d
	}
	if n >= 3 {
		if d := DnDensity(n); d > best {
			best = d
		}
	}
	switch n {
	case 6:
		if d := E6Density(); d > best {
			best = d
		}
	case 7:
		if d := E7Density(); d > best {
			best = d
		}
	case 8:
		if d := E8Density(); d > best {
			best = d
		}
	case 24:
		if d := LeechDensity(); d > best {
			best = d
		}
	}
	return best
}

// BestLatticeCenterDensity returns the center density corresponding to
// [BestLatticeDensity] in dimension n.
func BestLatticeCenterDensity(n int) float64 {
	return BestLatticeDensity(n) / UnitBallVolume(n)
}

// ----------------------------------------------------------------------------
// Hermite constant and density bounds.
// ----------------------------------------------------------------------------

// HermiteConstant returns the Hermite constant gamma_n, the maximum over all
// n-dimensional lattices of the ratio (minimal norm) / covolume^(2/n), for the
// dimensions in which it is proven (1..8 and 24), and -1 otherwise. The densest
// lattice center density satisfies delta = (gamma_n/4)^(n/2).
func HermiteConstant(n int) float64 {
	switch n {
	case 1:
		return 1
	case 2:
		return math.Sqrt(4.0 / 3.0)
	case 3:
		return math.Cbrt(2)
	case 4:
		return math.Sqrt2
	case 5:
		return math.Pow(2, 3.0/5.0)
	case 6:
		return math.Pow(64.0/3.0, 1.0/6.0)
	case 7:
		return math.Pow(2, 6.0/7.0)
	case 8:
		return 2
	case 24:
		return 4
	default:
		return -1
	}
}

// HermiteConstantKnown reports whether the Hermite constant is proven in
// dimension n.
func HermiteConstantKnown(n int) bool { return HermiteConstant(n) >= 0 }

// riemannZeta returns the Riemann zeta function zeta(s) for real s > 1 using
// Euler-Maclaurin summation, accurate to roughly machine precision.
func riemannZeta(s float64) float64 {
	const N = 20
	sum := 0.0
	for k := 1; k < N; k++ {
		sum += math.Pow(float64(k), -s)
	}
	nf := float64(N)
	sum += math.Pow(nf, 1-s) / (s - 1)
	sum += 0.5 * math.Pow(nf, -s)
	// Euler-Maclaurin corrections with Bernoulli numbers B2=1/6, B4=-1/30.
	sum += (1.0 / 12.0) * s * math.Pow(nf, -s-1)
	sum += (-1.0 / 720.0) * s * (s + 1) * (s + 2) * math.Pow(nf, -s-3)
	return sum
}

// MinkowskiHlawkaBound returns the Minkowski-Hlawka lower bound on the optimal
// lattice packing density in dimension n, zeta(n) / 2^(n-1). Every dimension
// n >= 2 admits a lattice at least this dense. Returns NaN for n < 2.
func MinkowskiHlawkaBound(n int) float64 {
	if n < 2 {
		return math.NaN()
	}
	return riemannZeta(float64(n)) / math.Pow(2, float64(n-1))
}

// ----------------------------------------------------------------------------
// Aggregated lattice information.
// ----------------------------------------------------------------------------

// LatticeInfo bundles the standard packing and covering invariants of a
// lattice: its name, dimension, covolume, minimal norm (squared minimal vector
// length), packing radius, kissing number, center density, packing density,
// covering radius and covering thickness.
type LatticeInfo struct {
	Name           string
	Dimension      int
	Covolume       float64
	MinimalNorm    float64
	PackingRadius  float64
	KissingNumber  int
	CenterDensity  float64
	Density        float64
	CoveringRadius float64
	Thickness      float64
}

// ZnInfo returns the [LatticeInfo] for the integer lattice Z^n.
func ZnInfo(n int) LatticeInfo {
	return LatticeInfo{
		Name:           "Z^n",
		Dimension:      n,
		Covolume:       ZnCovolume(n),
		MinimalNorm:    ZnMinimalNorm(n),
		PackingRadius:  ZnPackingRadius(n),
		KissingNumber:  ZnKissingNumber(n),
		CenterDensity:  ZnCenterDensity(n),
		Density:        ZnDensity(n),
		CoveringRadius: ZnCoveringRadius(n),
		Thickness:      ZnThickness(n),
	}
}

// AnInfo returns the [LatticeInfo] for the root lattice A_n.
func AnInfo(n int) LatticeInfo {
	return LatticeInfo{
		Name:           "A_n",
		Dimension:      n,
		Covolume:       AnCovolume(n),
		MinimalNorm:    AnMinimalNorm(n),
		PackingRadius:  AnPackingRadius(n),
		KissingNumber:  AnKissingNumber(n),
		CenterDensity:  AnCenterDensity(n),
		Density:        AnDensity(n),
		CoveringRadius: AnCoveringRadius(n),
		Thickness:      AnThickness(n),
	}
}

// DnInfo returns the [LatticeInfo] for the checkerboard lattice D_n.
func DnInfo(n int) LatticeInfo {
	return LatticeInfo{
		Name:           "D_n",
		Dimension:      n,
		Covolume:       DnCovolume(n),
		MinimalNorm:    DnMinimalNorm(n),
		PackingRadius:  DnPackingRadius(n),
		KissingNumber:  DnKissingNumber(n),
		CenterDensity:  DnCenterDensity(n),
		Density:        DnDensity(n),
		CoveringRadius: DnCoveringRadius(n),
		Thickness:      DnThickness(n),
	}
}

// E6Info returns the [LatticeInfo] for the exceptional lattice E6.
func E6Info() LatticeInfo {
	return LatticeInfo{
		Name:           "E6",
		Dimension:      6,
		Covolume:       E6Covolume(),
		MinimalNorm:    E6MinimalNorm(),
		PackingRadius:  E6PackingRadius(),
		KissingNumber:  E6KissingNumber(),
		CenterDensity:  E6CenterDensity(),
		Density:        E6Density(),
		CoveringRadius: E6CoveringRadius(),
		Thickness:      E6Thickness(),
	}
}

// E7Info returns the [LatticeInfo] for the exceptional lattice E7.
func E7Info() LatticeInfo {
	return LatticeInfo{
		Name:           "E7",
		Dimension:      7,
		Covolume:       E7Covolume(),
		MinimalNorm:    E7MinimalNorm(),
		PackingRadius:  E7PackingRadius(),
		KissingNumber:  E7KissingNumber(),
		CenterDensity:  E7CenterDensity(),
		Density:        E7Density(),
		CoveringRadius: E7CoveringRadius(),
		Thickness:      E7Thickness(),
	}
}

// E8Info returns the [LatticeInfo] for the exceptional lattice E8.
func E8Info() LatticeInfo {
	return LatticeInfo{
		Name:           "E8",
		Dimension:      8,
		Covolume:       E8Covolume(),
		MinimalNorm:    E8MinimalNorm(),
		PackingRadius:  E8PackingRadius(),
		KissingNumber:  E8KissingNumber(),
		CenterDensity:  E8CenterDensity(),
		Density:        E8Density(),
		CoveringRadius: E8CoveringRadius(),
		Thickness:      E8Thickness(),
	}
}

// LeechInfo returns the [LatticeInfo] for the Leech lattice Lambda_24.
func LeechInfo() LatticeInfo {
	return LatticeInfo{
		Name:           "Leech",
		Dimension:      24,
		Covolume:       LeechCovolume(),
		MinimalNorm:    LeechMinimalNorm(),
		PackingRadius:  LeechPackingRadius(),
		KissingNumber:  LeechKissingNumber(),
		CenterDensity:  LeechCenterDensity(),
		Density:        LeechDensity(),
		CoveringRadius: LeechCoveringRadius(),
		Thickness:      LeechThickness(),
	}
}
