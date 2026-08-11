# fractal — provides pure–Go, standard-library-only tools for computing and exploring classic

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/fractal.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/fractal)

Package fractal provides pure–Go, standard-library-only tools for computing
and exploring classic fractals.

- Escape-time fractals: the Mandelbrot set and Julia sets, including per-point
escape iteration counts, smooth (normalized) iteration counts for continuous
coloring, closed-form interior tests (main cardioid and period-2 bulb), and
evaluation over a rectangular grid via a `Viewport`. - Fractal dimension:
box-counting dimension of a finite point set, a generic log-log least-squares
slope helper, and the exact self-similar (Hausdorff) dimension of strictly
self-similar sets. - L-systems: deterministic context-free string rewriting
(`LSystem`) plus a turtle-graphics interpreter that turns a command string
into line `Segment`s or a vertex path, with presets for the Koch curve,
Sierpinski triangle, dragon curve, and a branching plant. - Iterated function
systems: affine contraction maps (`AffineMap`), the random "chaos game"
(`IFS.ChaosGame`), and presets for the Barnsley fern and Sierpinski
triangle/carpet. - Deterministic geometric fractals: the Koch curve and
snowflake, the Sierpinski triangle by recursive subdivision, and the Cantor
set.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/fractal"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func BoxCount(points []Point2D, boxSize float64) int` | BoxCount returns the number of distinct axis-aligned grid boxes of side boxSize that contain at least one of the given points. |
| `func BoxCountingDimension(points []Point2D, boxSizes []float64) float64` | BoxCountingDimension estimates the box-counting (Minkowski–Bouligand) dimension of a finite point set by counting occupied boxes at each of the… |
| `func BoxCountingDimensionFromCounts(boxSizes, counts []float64) float64` | BoxCountingDimensionFromCounts estimates a fractal dimension directly from a table of (box size, occupied-box count) samples: it returns the slope of… |
| `func FitLine(xs, ys []float64) (slope, intercept float64)` | FitLine returns the slope and intercept of the ordinary least-squares line y = slope*x + intercept fitted to the paired samples xs, ys. |
| `func FitSlope(xs, ys []float64) float64` | FitSlope returns the slope of the ordinary least-squares line y = m*x + b fitted to the paired samples xs, ys. |
| `func HausdorffDimensionSelfSimilar(n int, scale float64) float64` | HausdorffDimensionSelfSimilar returns the exact similarity (Hausdorff) dimension of a strictly self-similar set composed of n non-overlapping copies… |
| `func InJuliaSet(z, c complex128, maxIter int) bool` | InJuliaSet reports whether the seed z lies in the filled Julia set of z -> z^2 + c, i.e. |
| `func InMainCardioid(c complex128) bool` | InMainCardioid reports, in closed form, whether c lies inside the main cardioid of the Mandelbrot set — the large heart-shaped region containing… |
| `func InMandelbrotSet(c complex128, maxIter int) bool` | InMandelbrotSet reports whether c lies in the Mandelbrot set, i.e. |
| `func InPeriod2Bulb(c complex128) bool` | InPeriod2Bulb reports, in closed form, whether c lies inside the period-2 bulb of the Mandelbrot set: the disk of radius 1/4 centered at -1. |
| `func JuliaSmooth(z, c complex128, maxIter int, bailout float64) float64` | JuliaSmooth returns the smooth (fractional) escape iteration count for the Julia map z -> z^2 + c seeded at z. |
| `func MandelbrotSmooth(c complex128, maxIter int, bailout float64) float64` | MandelbrotSmooth returns the smooth (fractional) escape iteration count for the Mandelbrot map at c. |
| `func Orbit(z0, c complex128, n int) []complex128` | Orbit returns the first n+1 points of the orbit of z -> z^2 + c starting at z0, that is [z0, f(z0), f(f(z0)), ...] of length n+1. |

### Types

| Type | What it is |
| --- | --- |
| `AffineMap` | AffineMap is a 2-D affine transformation of the form |
| `EscapeResult` | EscapeResult records the outcome of iterating a quadratic map z -> z^2 + c until the orbit escapes a bailout radius or a maximum iteration count is… |
| `Grid` | Grid is a dense Width×Height array of scalar values stored in row-major order (index = y*Width + x). |
| `IFS` | IFS is an iterated function system: a list of affine contraction Maps with an equal-length list of selection Weights. |
| `Interval` | Interval is a closed real interval [Start, End]. |
| `LSystem` | LSystem is a deterministic context-free Lindenmayer system: an Axiom string and a set of production Rules mapping a single symbol (rune) to a… |
| `Point2D` | Point2D is a point in the Euclidean plane with coordinates X and Y. |
| `Segment` | Segment is a directed line segment from (X0,Y0) to (X1,Y1), the drawing primitive produced by the turtle interpreter. |
| `Triangle` | Triangle is a triangle with vertices A, B and C. |
| `TurtleConfig` | TurtleConfig parameterizes the turtle-graphics interpreter. |
| `TurtleState` | TurtleState is the pose of the turtle: its position and its heading in degrees measured counterclockwise from the positive x-axis. |
| `Viewport` | Viewport is an axis-aligned rectangle in the complex plane, used to map integer pixel coordinates onto complex values. |

<details>
<summary><code>AffineMap</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewAffineMap(a, b, c, d, e, f float64) AffineMap` | NewAffineMap constructs an AffineMap from its six coefficients in the order a, b, c, d, e, f matching the field layout x' = a*x+b*y+e, y' = c*x+d*y+f. |
| `func (m AffineMap) Apply(p Point2D) Point2D` | Apply returns the image of p under the affine map. |

</details>

<details>
<summary><code>EscapeResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func JuliaEscape(z, c complex128, maxIter int, bailout float64) EscapeResult` | JuliaEscape computes the escape-time result for the filled Julia set of the map z -> z^2 + c starting from the seed point z. |
| `func MandelbrotEscape(c complex128, maxIter int, bailout float64) EscapeResult` | MandelbrotEscape computes the escape-time result for the Mandelbrot map z -> z^2 + c starting from z = 0. |
| `func (r EscapeResult) Smooth(bailout float64) float64` | Smooth returns the fractional (normalized) iteration count for an escaped result, giving continuous values suitable for smooth coloring. |

</details>

<details>
<summary><code>Grid</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func JuliaGrid(c complex128, v Viewport, width, height, maxIter int, bailout float64) *Grid` | JuliaGrid rasterizes the filled Julia set of z -> z^2 + c over the given viewport into a width×height `Grid`. |
| `func MandelbrotGrid(v Viewport, width, height, maxIter int, bailout float64) *Grid` | MandelbrotGrid rasterizes the Mandelbrot set over the given viewport into a width×height `Grid`. |
| `func NewGrid(width, height int) *Grid` | NewGrid allocates a zero-filled Grid with the given dimensions. |
| `func (g *Grid) At(x, y int) float64` | At returns the value at column x, row y. |
| `func (g *Grid) Clone() *Grid` | Clone returns a deep copy of the grid. |
| `func (g *Grid) CountFinite() int` | CountFinite returns the number of grid cells whose value is finite (neither NaN nor infinite). |
| `func (g *Grid) Max() float64` | Max returns the largest value in the grid, or 0 for an empty grid. |
| `func (g *Grid) Mean() float64` | Mean returns the arithmetic mean of all grid values, or 0 for an empty grid. |
| `func (g *Grid) Min() float64` | Min returns the smallest value in the grid, or 0 for an empty grid. |
| `func (g *Grid) MinMax() (min, max float64)` | MinMax returns the minimum and maximum values in the grid. |
| `func (g *Grid) Normalize() *Grid` | Normalize returns a new Grid whose values are linearly rescaled so that the original minimum maps to 0 and the original maximum maps to 1. |
| `func (g *Grid) Row(y int) []float64` | Row returns the y-th row of the grid as a slice that aliases the underlying storage. |
| `func (g *Grid) Set(x, y int, v float64)` | Set stores v at column x, row y. |

</details>

<details>
<summary><code>IFS</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BarnsleyFern() IFS` | BarnsleyFern returns the four-map iterated function system that generates Michael Barnsley's fern. |
| `func SierpinskiCarpetIFS() IFS` | SierpinskiCarpetIFS returns the eight-map iterated function system whose attractor is the Sierpinski carpet: the unit square is divided into a 3x3… |
| `func SierpinskiTriangleIFS() IFS` | SierpinskiTriangleIFS returns the three-map iterated function system whose attractor is the Sierpinski triangle with vertices (0,0), (1,0) and (0.5,… |
| `func (s IFS) ChaosGame(n int, seed int64) []Point2D` | ChaosGame runs the chaos game starting from the origin, discarding 20 transient points, and returns n attractor points. |
| `func (s IFS) ChaosGameFrom(start Point2D, n, transient int, seed int64) []Point2D` | ChaosGameFrom runs the random-iteration ("chaos game") algorithm starting from the point start: it repeatedly picks a map according to the system's… |

</details>

<details>
<summary><code>Interval</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func CantorSet(start, end float64, iterations int) []Interval` | CantorSet returns the intervals remaining after iterations steps of the middle-thirds Cantor construction applied to [start, end]. |
| `func (iv Interval) Length() float64` | Length returns End - Start. |

</details>

<details>
<summary><code>LSystem</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DragonLSystem() LSystem` | DragonLSystem returns the Heighway dragon curve L-system: axiom "FX" with rules X -> X+YF+ and Y -> -FX-Y, rendered with a turn angle of 90 degrees. |
| `func KochLSystem() LSystem` | KochLSystem returns the classic Koch curve L-system: axiom "F" with the rule F -> F+F--F+F. Rendered with a turn angle of 60 degrees it produces the… |
| `func NewLSystem(axiom string, rules map[rune]string) LSystem` | NewLSystem constructs an LSystem from an axiom and a rule map. |
| `func PlantLSystem() LSystem` | PlantLSystem returns a branching plant L-system: axiom "X" with rules X -> F+[`X`-X]-F[-FX]+X and F -> FF, rendered with a turn angle of about 25… |
| `func SierpinskiLSystem() LSystem` | SierpinskiLSystem returns an L-system whose turtle rendering draws the Sierpinski triangle. |
| `func (l LSystem) Expand(iterations int) string` | Expand applies the production rules iterations times starting from the axiom and returns the resulting string. |
| `func (l LSystem) Step(s string) string` | Step applies every production rule once to s in parallel, replacing each symbol by its rule's right-hand side (or leaving it unchanged if it has no… |

</details>

<details>
<summary><code>Point2D</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BoundingBox(pts []Point2D) (min, max Point2D)` | BoundingBox returns the axis-aligned bounding box of pts as its minimum and maximum corners. |
| `func KochCurve(p0, p1 Point2D, iterations int) []Point2D` | KochCurve returns the polyline vertices of the Koch curve constructed between the endpoints p0 and p1 after the given number of iterations. |
| `func KochSnowflake(center Point2D, radius float64, iterations int) []Point2D` | KochSnowflake returns the closed polyline of the Koch snowflake inscribed in the circle of the given radius centered at center, after the given… |
| `func TurtlePath(commands string, cfg TurtleConfig) []Point2D` | TurtlePath interprets a command string and returns the ordered vertices visited by the turtle: the start position followed by the endpoint of every… |
| `func (p Point2D) Add(q Point2D) Point2D` | Add returns the vector sum p+q. |
| `func (p Point2D) Dist(q Point2D) float64` | Dist returns the Euclidean distance between p and q. |
| `func (p Point2D) Lerp(q Point2D, t float64) Point2D` | Lerp returns the linear interpolation between p (t=0) and q (t=1). |
| `func (p Point2D) Midpoint(q Point2D) Point2D` | Midpoint returns the point halfway between p and q. |
| `func (p Point2D) Norm() float64` | Norm returns the Euclidean length of p treated as a vector from the origin. |
| `func (p Point2D) Rotate(angle float64) Point2D` | Rotate returns p rotated counterclockwise about the origin by angle radians. |
| `func (p Point2D) Scale(s float64) Point2D` | Scale returns p with both coordinates multiplied by s. |
| `func (p Point2D) Sub(q Point2D) Point2D` | Sub returns the vector difference p-q. |

</details>

<details>
<summary><code>Segment</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Turtle(commands string, cfg TurtleConfig) []Segment` | Turtle interprets a command string as turtle graphics and returns the drawn line segments. |

</details>

<details>
<summary><code>Triangle</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func SierpinskiTriangle(a, b, c Point2D, depth int) []Triangle` | SierpinskiTriangle returns the list of filled sub-triangles produced by recursively removing the central triangle from the triangle (a, b, c) to the… |

</details>

<details>
<summary><code>Viewport</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewViewport(centerRe, centerIm, radius float64) Viewport` | NewViewport returns the square viewport centered at (centerRe, centerIm) with the given half-width (radius) along each axis. |
| `func (v Viewport) Aspect() float64` | Aspect returns the width-to-height ratio of the viewport. |
| `func (v Viewport) Center() complex128` | Center returns the complex number at the center of the viewport. |
| `func (v Viewport) PixelToComplex(px, py, width, height int) complex128` | PixelToComplex maps the pixel at column px, row py of a width×height raster to the complex value it samples. |
| `func (v Viewport) SpanX() float64` | SpanX returns the width of the viewport along the real axis. |
| `func (v Viewport) SpanY() float64` | SpanY returns the height of the viewport along the imaginary axis. |
| `func (v Viewport) Zoom(factor float64) Viewport` | Zoom returns a new viewport with the same center whose span along each axis is multiplied by factor. |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/fractal).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
