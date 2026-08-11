# plot — small, Matplotlib-shaped plotting library that renders natively in pure Go using only

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/plot.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/plot)

Package plot is a small, Matplotlib-shaped plotting library that renders
natively in pure Go using only the standard library (image, image/png,
image/color) plus an in-package SVG string writer. It has no third-party
dependencies. The only non-standard package it touches is the parent
`github.com/malcolmston/algebra` package, and only in `Axes.PlotExpr`, which
samples a parsed symbolic expression.

The API mirrors Matplotlib's Figure/Axes structure. A `Figure` is a fixed
pixel canvas holding a single `Axes`; the axes owns the data limits, labels
and the ordered list of series. Create a figure with `New` or `NewWithSize`,
obtain its axes with `Figure.Axes`, add series, then render.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/plot"
```

## Usage

This is the package's own `ExampleAxes_Plot`, so it compiles and its output is
asserted on every `go test ./plot/`.

```go
fig := plot.New()
	line := fig.Axes().Plot([]float64{0, 1, 2, 3}, []float64{0, 1, 4, 9})
	line.Label = "x squared"
	fig.Axes().Legend(true)
	fmt.Println(line.Label, len(line.Xs), len(line.Ys))
```

```
x squared 4 4
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func ColorCycle(i int) color.RGBA` | ColorCycle returns the color at position i in `DefaultColors`, wrapping around modulo the palette length. |

### Types

| Type | What it is |
| --- | --- |
| `Axes` | Axes is a single set of x and y coordinate axes together with the series drawn on them, analogous to a Matplotlib Axes. |
| `BarSeries` | BarSeries is a vertical bar chart, the output of `Axes.Bar`. |
| `Figure` | Figure is the top-level drawing surface, analogous to a Matplotlib Figure. |
| `FillSeries` | FillSeries is a shaded region between two y curves sharing common x values, the output of `Axes.FillBetween`. |
| `HistSeries` | HistSeries is a histogram computed from raw samples by `Axes.Hist`. |
| `LineSeries` | LineSeries is a poly-line connecting successive (X, Y) data points, the output of `Axes.Plot`, `Axes.PlotFunc` and `Axes.PlotExpr`. |
| `ScatterSeries` | ScatterSeries is a set of unconnected point markers, the output of `Axes.Scatter`. |
| `StepSeries` | StepSeries is a piecewise-constant step plot, the output of `Axes.Step`. |

<details>
<summary><code>Axes</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (a *Axes) Bar(xs, heights []float64, width float64) *BarSeries` | Bar adds a vertical bar chart with bars centred on xs and the given heights, rising from y = 0. |
| `func (a *Axes) FillBetween(xs, y1, y2 []float64) *FillSeries` | FillBetween adds a shaded band between the two curves y1 and y2 sampled at the shared x values and returns the new series. |
| `func (a *Axes) Grid(on bool) *Axes` | Grid enables or disables the background grid lines at the tick locations. |
| `func (a *Axes) Hist(samples []float64, bins int) *HistSeries` | Hist bins the raw samples into the given number of equal-width bins spanning the sample range and adds the resulting histogram. |
| `func (a *Axes) Legend(on bool) *Axes` | Legend enables or disables the legend box, which lists every series that has a non-empty label. |
| `func (a *Axes) Plot(xs, ys []float64) *LineSeries` | Plot adds a line connecting the given points and returns the new series for further customization. |
| `func (a *Axes) PlotExpr(src, varName string, xmin, xmax float64, n int) (*LineSeries, error)` | PlotExpr parses the algebra expression src, samples it at n evenly spaced points on [xmin, xmax] by binding the given variable name to each sample,… |
| `func (a *Axes) PlotFunc(f func(float64) float64, xmin, xmax float64, n int) *LineSeries` | PlotFunc samples the Go function f at n evenly spaced points on [xmin, xmax] and adds the resulting line, returning the new series. |
| `func (a *Axes) Scatter(xs, ys []float64) *ScatterSeries` | Scatter adds unconnected point markers and returns the new series. |
| `func (a *Axes) Step(xs, ys []float64) *StepSeries` | Step adds a piecewise-constant step plot through the given points and returns the new series. |
| `func (a *Axes) Title(s string) *Axes` | Title sets the axes title, drawn centred above the plot. |
| `func (a *Axes) XLabel(s string) *Axes` | XLabel sets the x-axis label, drawn centred below the plot. |
| `func (a *Axes) XLim(lo, hi float64) *Axes` | XLim fixes the x-axis data range to [lo, hi], disabling autoscaling in x. |
| `func (a *Axes) YLabel(s string) *Axes` | YLabel sets the y-axis label, drawn rotated along the left edge. |
| `func (a *Axes) YLim(lo, hi float64) *Axes` | YLim fixes the y-axis data range to [lo, hi], disabling autoscaling in y. |

</details>

<details>
<summary><code>Figure</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func New() *Figure` | New returns a new `Figure` of the default size (640x480) with a white background and an empty axes. |
| `func NewWithSize(w, h int) *Figure` | NewWithSize returns a new `Figure` of the given pixel size with a white background and an empty axes. |
| `func (f *Figure) Axes() *Axes` | Axes returns the figure's single axes, on which all series are plotted. |
| `func (f *Figure) ExportMatplotlib() string` | ExportMatplotlib returns a runnable Python script that reproduces the figure using the real Matplotlib library. |
| `func (f *Figure) RenderPNG() ([]byte, error)` | RenderPNG renders the figure and returns the encoded PNG bytes. |
| `func (f *Figure) RenderSVG() string` | RenderSVG renders the figure and returns the SVG document as a string. |
| `func (f *Figure) SavePNG(path string) error` | SavePNG renders the figure and writes the PNG to the file at path. |
| `func (f *Figure) SaveSVG(path string) error` | SaveSVG renders the figure and writes the SVG document to the file at path. |

</details>

### Variables

`Blue`, `Orange`, `Green`, `Red`, `Purple`, `Brown`, `Pink`, `Gray`, `Olive`, `Cyan`, `Black`, `White`, `LightGray`, `DefaultColors`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/plot).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
