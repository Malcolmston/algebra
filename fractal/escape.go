package fractal

import (
	"math"
	"math/cmplx"
)

// EscapeResult records the outcome of iterating a quadratic map z -> z^2 + c
// until the orbit escapes a bailout radius or a maximum iteration count is
// reached. Escaped reports whether the bailout was exceeded; Iterations is the
// number of completed iterations at that moment (equal to the maximum when the
// orbit did not escape); FinalZ is the last computed orbit value.
type EscapeResult struct {
	Escaped    bool
	Iterations int
	FinalZ     complex128
}

// fractalEscape iterates z -> z^2 + c starting from z0 and returns an
// EscapeResult. The orbit is deemed to have escaped as soon as |z| > bailout.
func fractalEscape(z0, c complex128, maxIter int, bailout float64) EscapeResult {
	z := z0
	b2 := bailout * bailout
	for n := 0; n < maxIter; n++ {
		re, im := real(z), imag(z)
		if re*re+im*im > b2 {
			return EscapeResult{Escaped: true, Iterations: n, FinalZ: z}
		}
		z = z*z + c
	}
	return EscapeResult{Escaped: false, Iterations: maxIter, FinalZ: z}
}

// MandelbrotEscape computes the escape-time result for the Mandelbrot map
// z -> z^2 + c starting from z = 0. The point c is considered to have escaped
// once |z| exceeds bailout (which should be at least 2, the Mandelbrot escape
// radius). Iteration stops after maxIter steps.
func MandelbrotEscape(c complex128, maxIter int, bailout float64) EscapeResult {
	return fractalEscape(0, c, maxIter, bailout)
}

// JuliaEscape computes the escape-time result for the filled Julia set of the
// map z -> z^2 + c starting from the seed point z. The orbit escapes once |z|
// exceeds bailout (at least 2). Iteration stops after maxIter steps.
func JuliaEscape(z, c complex128, maxIter int, bailout float64) EscapeResult {
	return fractalEscape(z, c, maxIter, bailout)
}

// Smooth returns the fractional (normalized) iteration count for an escaped
// result, giving continuous values suitable for smooth coloring. It uses the
// standard renormalization mu = n + 1 - log2(log|z| / log(bailout)). For a
// result that did not escape it returns the integer iteration count unchanged.
// The estimate is most accurate when bailout is large (for example 256).
func (r EscapeResult) Smooth(bailout float64) float64 {
	if !r.Escaped {
		return float64(r.Iterations)
	}
	az := cmplx.Abs(r.FinalZ)
	if az <= 1 {
		return float64(r.Iterations)
	}
	nu := math.Log(math.Log(az)/math.Log(bailout)) / math.Log(2)
	return float64(r.Iterations) + 1 - nu
}

// MandelbrotSmooth returns the smooth (fractional) escape iteration count for
// the Mandelbrot map at c. See [EscapeResult.Smooth]. A large bailout (for
// example 256) gives the smoothest result.
func MandelbrotSmooth(c complex128, maxIter int, bailout float64) float64 {
	return MandelbrotEscape(c, maxIter, bailout).Smooth(bailout)
}

// JuliaSmooth returns the smooth (fractional) escape iteration count for the
// Julia map z -> z^2 + c seeded at z. See [EscapeResult.Smooth].
func JuliaSmooth(z, c complex128, maxIter int, bailout float64) float64 {
	return JuliaEscape(z, c, maxIter, bailout).Smooth(bailout)
}

// InMandelbrotSet reports whether c lies in the Mandelbrot set, i.e. whether
// the orbit of z -> z^2 + c starting at 0 fails to escape radius 2 within
// maxIter iterations. Because the test is truncated at maxIter, points very
// near the boundary may be misclassified as members.
func InMandelbrotSet(c complex128, maxIter int) bool {
	return !MandelbrotEscape(c, maxIter, 2).Escaped
}

// InJuliaSet reports whether the seed z lies in the filled Julia set of
// z -> z^2 + c, i.e. whether its orbit fails to escape radius 2 within maxIter
// iterations.
func InJuliaSet(z, c complex128, maxIter int) bool {
	return !JuliaEscape(z, c, maxIter, 2).Escaped
}

// InMainCardioid reports, in closed form, whether c lies inside the main
// cardioid of the Mandelbrot set — the large heart-shaped region containing all
// c for which z -> z^2 + c has an attracting fixed point. Writing
// p = |c - 1/4| = sqrt((Re c - 1/4)^2 + (Im c)^2), the membership criterion is
// Re(c) <= p - 2*p^2 + 1/4, which holds exactly on the closed cardioid.
func InMainCardioid(c complex128) bool {
	re, im := real(c), imag(c)
	dx := re - 0.25
	p := math.Hypot(dx, im)
	return re <= p-2*p*p+0.25
}

// InPeriod2Bulb reports, in closed form, whether c lies inside the period-2
// bulb of the Mandelbrot set: the disk of radius 1/4 centered at -1. Membership
// is exactly |c + 1| <= 1/4.
func InPeriod2Bulb(c complex128) bool {
	d := c + 1
	return real(d)*real(d)+imag(d)*imag(d) <= 1.0/16.0
}

// Orbit returns the first n+1 points of the orbit of z -> z^2 + c starting at
// z0, that is [z0, f(z0), f(f(z0)), ...] of length n+1. For n < 0 it returns
// nil.
func Orbit(z0, c complex128, n int) []complex128 {
	if n < 0 {
		return nil
	}
	out := make([]complex128, n+1)
	z := z0
	out[0] = z
	for i := 1; i <= n; i++ {
		z = z*z + c
		out[i] = z
	}
	return out
}

// MandelbrotGrid rasterizes the Mandelbrot set over the given viewport into a
// width×height [Grid]. Each cell holds the smooth escape iteration count (see
// [EscapeResult.Smooth]) of the complex value sampled at that pixel, so interior
// (non-escaping) points hold maxIter. Row 0 corresponds to the top of the
// viewport. See [Viewport.PixelToComplex] for the pixel-to-complex mapping.
func MandelbrotGrid(v Viewport, width, height, maxIter int, bailout float64) *Grid {
	g := NewGrid(width, height)
	for py := 0; py < height; py++ {
		for px := 0; px < width; px++ {
			c := v.PixelToComplex(px, py, width, height)
			g.Data[py*width+px] = MandelbrotEscape(c, maxIter, bailout).Smooth(bailout)
		}
	}
	return g
}

// JuliaGrid rasterizes the filled Julia set of z -> z^2 + c over the given
// viewport into a width×height [Grid]. Each cell holds the smooth escape
// iteration count of the seed sampled at that pixel. Row 0 corresponds to the
// top of the viewport.
func JuliaGrid(c complex128, v Viewport, width, height, maxIter int, bailout float64) *Grid {
	g := NewGrid(width, height)
	for py := 0; py < height; py++ {
		for px := 0; px < width; px++ {
			z := v.PixelToComplex(px, py, width, height)
			g.Data[py*width+px] = JuliaEscape(z, c, maxIter, bailout).Smooth(bailout)
		}
	}
	return g
}
