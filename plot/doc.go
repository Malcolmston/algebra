// Package plot is a small, Matplotlib-shaped plotting library that renders
// natively in pure Go using only the standard library (image, image/png,
// image/color) plus an in-package SVG string writer. It has no third-party
// dependencies. The only non-standard package it touches is the parent
// [github.com/malcolmston/algebra] package, and only in [Axes.PlotExpr], which
// samples a parsed symbolic expression.
//
// # Model
//
// The API mirrors Matplotlib's Figure/Axes structure. A [Figure] is a fixed
// pixel canvas holding a single [Axes]; the axes owns the data limits, labels
// and the ordered list of series. Create a figure with [New] or [NewWithSize],
// obtain its axes with [Figure.Axes], add series, then render.
//
//	fig := plot.New()
//	ax := fig.Axes()
//	ax.PlotFunc(math.Sin, 0, 2*math.Pi, 200)
//	ax.Title("sine").XLabel("x").YLabel("sin(x)").Grid(true)
//	png, _ := fig.RenderPNG()
//
// # Series types
//
// The axes supports the common Matplotlib chart types, each returning a typed
// series value whose exported fields can be tweaked after creation:
//
//   - [Axes.Plot]        — a connected poly-line ([LineSeries]).
//   - [Axes.Scatter]     — unconnected markers ([ScatterSeries]).
//   - [Axes.Bar]         — a vertical bar chart ([BarSeries]).
//   - [Axes.Hist]        — a histogram of raw samples ([HistSeries]).
//   - [Axes.Step]        — a piecewise-constant step plot ([StepSeries]).
//   - [Axes.FillBetween] — a shaded band between two curves ([FillSeries]).
//   - [Axes.PlotFunc]    — samples a Go func(float64)float64 over a range.
//   - [Axes.PlotExpr]    — samples an algebra expression string over a range.
//
// Series colors are assigned automatically from [DefaultColors] in the order
// the series are added, matching Matplotlib's property cycle.
//
// # Configuration
//
// [Axes.Title], [Axes.XLabel], [Axes.YLabel], [Axes.Grid], [Axes.Legend],
// [Axes.XLim] and [Axes.YLim] configure the axes and return the receiver so
// they can be chained. When no explicit limits are set the axes autoscales to
// the data with a small margin, and axis ticks are placed at "nice" round
// values.
//
// # Rendering
//
// Two fully native backends are provided. The raster backend
// ([Figure.RenderPNG], [Figure.SavePNG]) draws into an [image.RGBA] using
// Bresenham line rasterization and a built-in 5x7 bitmap font, then encodes the
// result with image/png. The vector backend ([Figure.RenderSVG],
// [Figure.SaveSVG]) emits an SVG document with real paths and text. Both
// backends are deterministic for a given figure.
//
// # Matplotlib export
//
// [Figure.ExportMatplotlib] emits a runnable Python script that reproduces the
// figure with the real Matplotlib library. This keeps the Go package
// dependency-free while still offering a route to genuine Matplotlib output for
// users who have Python and matplotlib available.
package plot
