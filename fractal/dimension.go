package fractal

import "math"

// FitSlope returns the slope of the ordinary least-squares line y = m*x + b
// fitted to the paired samples xs, ys. It panics if the slices differ in length
// or contain fewer than two points, and returns 0 when every x is identical.
func FitSlope(xs, ys []float64) float64 {
	m, _ := FitLine(xs, ys)
	return m
}

// FitLine returns the slope and intercept of the ordinary least-squares line
// y = slope*x + intercept fitted to the paired samples xs, ys. It panics if the
// slices differ in length or contain fewer than two points. When every x is
// identical the slope is 0 and the intercept is the mean of ys.
func FitLine(xs, ys []float64) (slope, intercept float64) {
	if len(xs) != len(ys) {
		panic("fractal: FitLine length mismatch")
	}
	n := len(xs)
	if n < 2 {
		panic("fractal: FitLine needs at least two points")
	}
	var sx, sy, sxx, sxy float64
	for i := 0; i < n; i++ {
		sx += xs[i]
		sy += ys[i]
		sxx += xs[i] * xs[i]
		sxy += xs[i] * ys[i]
	}
	fn := float64(n)
	denom := fn*sxx - sx*sx
	if denom == 0 {
		return 0, sy / fn
	}
	slope = (fn*sxy - sx*sy) / denom
	intercept = (sy - slope*sx) / fn
	return slope, intercept
}

// BoxCount returns the number of distinct axis-aligned grid boxes of side
// boxSize that contain at least one of the given points. It panics if boxSize
// is not positive. This is the count N(epsilon) used in box-counting dimension
// estimation.
func BoxCount(points []Point2D, boxSize float64) int {
	if boxSize <= 0 {
		panic("fractal: BoxCount needs positive boxSize")
	}
	type cell struct{ ix, iy int }
	seen := make(map[cell]struct{}, len(points))
	inv := 1.0 / boxSize
	for _, p := range points {
		seen[cell{int(math.Floor(p.X * inv)), int(math.Floor(p.Y * inv))}] = struct{}{}
	}
	return len(seen)
}

// BoxCountingDimension estimates the box-counting (Minkowski–Bouligand)
// dimension of a finite point set by counting occupied boxes at each of the
// given box sizes and returning the slope of log N(epsilon) versus
// log(1/epsilon). At least two box sizes are required; the sizes must be
// positive. Smaller boxes give a better estimate for a true fractal but require
// densely sampled points.
func BoxCountingDimension(points []Point2D, boxSizes []float64) float64 {
	if len(boxSizes) < 2 {
		panic("fractal: BoxCountingDimension needs at least two box sizes")
	}
	xs := make([]float64, len(boxSizes))
	ys := make([]float64, len(boxSizes))
	for i, s := range boxSizes {
		xs[i] = math.Log(1 / s)
		ys[i] = math.Log(float64(BoxCount(points, s)))
	}
	return FitSlope(xs, ys)
}

// BoxCountingDimensionFromCounts estimates a fractal dimension directly from a
// table of (box size, occupied-box count) samples: it returns the slope of
// log(count) versus log(1/size). This is useful when counts were obtained
// externally. The slices must have equal length of at least two, and all
// entries must be positive. For an exact power law count = k*size^(-D) the
// returned value is exactly D.
func BoxCountingDimensionFromCounts(boxSizes, counts []float64) float64 {
	if len(boxSizes) != len(counts) {
		panic("fractal: BoxCountingDimensionFromCounts length mismatch")
	}
	if len(boxSizes) < 2 {
		panic("fractal: BoxCountingDimensionFromCounts needs at least two samples")
	}
	xs := make([]float64, len(boxSizes))
	ys := make([]float64, len(counts))
	for i := range boxSizes {
		xs[i] = math.Log(1 / boxSizes[i])
		ys[i] = math.Log(counts[i])
	}
	return FitSlope(xs, ys)
}

// HausdorffDimensionSelfSimilar returns the exact similarity (Hausdorff)
// dimension of a strictly self-similar set composed of n non-overlapping copies
// of itself, each scaled by the ratio scale (0 < scale < 1). The dimension is
// log(n)/log(1/scale). For example n=3, scale=1/2 gives the Sierpinski triangle
// dimension log2(3) ≈ 1.585; n=4, scale=1/3 gives the Koch curve dimension
// log(4)/log(3) ≈ 1.262. It panics if n < 1 or scale is not in (0,1).
func HausdorffDimensionSelfSimilar(n int, scale float64) float64 {
	if n < 1 {
		panic("fractal: HausdorffDimensionSelfSimilar needs n >= 1")
	}
	if scale <= 0 || scale >= 1 {
		panic("fractal: HausdorffDimensionSelfSimilar needs 0 < scale < 1")
	}
	return math.Log(float64(n)) / math.Log(1/scale)
}
