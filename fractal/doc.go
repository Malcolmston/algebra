// Package fractal provides pure–Go, standard-library-only tools for computing
// and exploring classic fractals.
//
// The package is organized around several independent topics:
//
//   - Escape-time fractals: the Mandelbrot set and Julia sets, including
//     per-point escape iteration counts, smooth (normalized) iteration counts
//     for continuous coloring, closed-form interior tests (main cardioid and
//     period-2 bulb), and evaluation over a rectangular grid via a [Viewport].
//   - Fractal dimension: box-counting dimension of a finite point set, a
//     generic log-log least-squares slope helper, and the exact self-similar
//     (Hausdorff) dimension of strictly self-similar sets.
//   - L-systems: deterministic context-free string rewriting ([LSystem]) plus a
//     turtle-graphics interpreter that turns a command string into line
//     [Segment]s or a vertex path, with presets for the Koch curve, Sierpinski
//     triangle, dragon curve, and a branching plant.
//   - Iterated function systems: affine contraction maps ([AffineMap]), the
//     random "chaos game" ([IFS.ChaosGame]), and presets for the Barnsley fern
//     and Sierpinski triangle/carpet.
//   - Deterministic geometric fractals: the Koch curve and snowflake,
//     the Sierpinski triangle by recursive subdivision, and the Cantor set.
//
// All routines use only the Go standard library and are deterministic: any
// randomness (the chaos game) is driven by an explicit caller-supplied seed.
//
// Complex numbers are represented with the builtin complex128 type; planar
// points use [Point2D].
package fractal
