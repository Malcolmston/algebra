package plot

import (
	"image/color"
	"math"

	"github.com/malcolmston/algebra"
)

// Default figure geometry, in pixels.
const (
	defaultWidth  = 640
	defaultHeight = 480

	marginLeft   = 60 // room for the y tick labels and y axis label
	marginRight  = 24
	marginTop    = 36 // room for the title
	marginBottom = 48 // room for the x tick labels and x axis label
)

// Figure is the top-level drawing surface, analogous to a Matplotlib Figure.
// It has a pixel size and a single [Axes] on which series are plotted. Create
// one with [New] or [NewWithSize], add series through the axes, then render
// with [Figure.RenderPNG], [Figure.RenderSVG], [Figure.SavePNG],
// [Figure.SaveSVG] or [Figure.ExportMatplotlib].
type Figure struct {
	// Width is the figure width in pixels.
	Width int
	// Height is the figure height in pixels.
	Height int
	// Background is the fill color of the whole figure.
	Background color.RGBA

	ax *Axes
}

// Axes is a single set of x and y coordinate axes together with the series
// drawn on them, analogous to a Matplotlib Axes. Obtain the axes of a figure
// with [Figure.Axes]. All plotting methods hang off this type.
type Axes struct {
	title  string
	xlabel string
	ylabel string

	grid   bool
	legend bool

	xlimSet bool
	ylimSet bool
	xmin    float64
	xmax    float64
	ymin    float64
	ymax    float64

	series []series
}

// New returns a new [Figure] of the default size (640x480) with a white
// background and an empty axes.
func New() *Figure {
	return NewWithSize(defaultWidth, defaultHeight)
}

// NewWithSize returns a new [Figure] of the given pixel size with a white
// background and an empty axes. Non-positive dimensions are replaced by the
// defaults.
func NewWithSize(w, h int) *Figure {
	if w <= 0 {
		w = defaultWidth
	}
	if h <= 0 {
		h = defaultHeight
	}
	return &Figure{Width: w, Height: h, Background: White, ax: &Axes{}}
}

// Axes returns the figure's single axes, on which all series are plotted.
func (f *Figure) Axes() *Axes { return f.ax }

// series is the common internal interface implemented by every plottable item.
type series interface {
	// bounds returns the data-space bounding box of the series and whether it
	// contains any finite points.
	bounds() (xmin, xmax, ymin, ymax float64, ok bool)
	// styleColor returns the series color.
	styleColor() color.RGBA
	// legendLabel returns the legend label ("" if the series has none).
	legendLabel() string
}

// --- Line ---------------------------------------------------------------

// LineSeries is a poly-line connecting successive (X, Y) data points, the
// output of [Axes.Plot], [Axes.PlotFunc] and [Axes.PlotExpr]. Fields may be
// adjusted after creation to change the rendered appearance.
type LineSeries struct {
	// Xs and Ys are the point coordinates; they must have equal length.
	Xs, Ys []float64
	// Color is the line color.
	Color color.RGBA
	// Label is the legend label; an empty string omits the series from the
	// legend.
	Label string
	// Width is the line thickness in pixels (minimum 1).
	Width int
}

func (s *LineSeries) styleColor() color.RGBA { return s.Color }
func (s *LineSeries) legendLabel() string    { return s.Label }
func (s *LineSeries) bounds() (float64, float64, float64, float64, bool) {
	return xyBounds(s.Xs, s.Ys)
}

// --- Scatter ------------------------------------------------------------

// ScatterSeries is a set of unconnected point markers, the output of
// [Axes.Scatter].
type ScatterSeries struct {
	// Xs and Ys are the point coordinates; they must have equal length.
	Xs, Ys []float64
	// Color is the marker fill color.
	Color color.RGBA
	// Label is the legend label.
	Label string
	// Size is the marker radius in pixels (minimum 1).
	Size int
}

func (s *ScatterSeries) styleColor() color.RGBA { return s.Color }
func (s *ScatterSeries) legendLabel() string    { return s.Label }
func (s *ScatterSeries) bounds() (float64, float64, float64, float64, bool) {
	return xyBounds(s.Xs, s.Ys)
}

// --- Bar ----------------------------------------------------------------

// BarSeries is a vertical bar chart, the output of [Axes.Bar]. Each bar is
// centred on Xs[i] with height Heights[i] measured from Baseline.
type BarSeries struct {
	// Xs are the bar centers.
	Xs []float64
	// Heights are the bar heights measured from Baseline.
	Heights []float64
	// Width is the bar width in data units.
	Width float64
	// Baseline is the y value the bars rise from (usually 0).
	Baseline float64
	// Color is the bar fill color.
	Color color.RGBA
	// Label is the legend label.
	Label string
}

func (s *BarSeries) styleColor() color.RGBA { return s.Color }
func (s *BarSeries) legendLabel() string    { return s.Label }
func (s *BarSeries) bounds() (float64, float64, float64, float64, bool) {
	if len(s.Xs) == 0 {
		return 0, 0, 0, 0, false
	}
	xmin, xmax := math.Inf(1), math.Inf(-1)
	ymin, ymax := s.Baseline, s.Baseline
	ok := false
	for i := range s.Xs {
		x, h := s.Xs[i], s.Heights[i]
		if isBad(x) || isBad(h) {
			continue
		}
		ok = true
		l, r := x-s.Width/2, x+s.Width/2
		xmin, xmax = math.Min(xmin, l), math.Max(xmax, r)
		top := s.Baseline + h
		ymin, ymax = math.Min(ymin, top), math.Max(ymax, top)
	}
	return xmin, xmax, ymin, ymax, ok
}

// --- Histogram ----------------------------------------------------------

// HistSeries is a histogram computed from raw samples by [Axes.Hist]. Edges
// has one more element than Counts; bin i spans [Edges[i], Edges[i+1]) and has
// height Counts[i].
type HistSeries struct {
	// Edges are the bin boundaries; len(Edges) == len(Counts)+1.
	Edges []float64
	// Counts are the per-bin sample counts.
	Counts []float64
	// Color is the bar fill color.
	Color color.RGBA
	// Label is the legend label.
	Label string
}

func (s *HistSeries) styleColor() color.RGBA { return s.Color }
func (s *HistSeries) legendLabel() string    { return s.Label }
func (s *HistSeries) bounds() (float64, float64, float64, float64, bool) {
	if len(s.Edges) < 2 {
		return 0, 0, 0, 0, false
	}
	xmin, xmax := s.Edges[0], s.Edges[len(s.Edges)-1]
	ymin, ymax := 0.0, 0.0
	for _, c := range s.Counts {
		if c > ymax {
			ymax = c
		}
	}
	return xmin, xmax, ymin, ymax, true
}

// --- Step ---------------------------------------------------------------

// StepSeries is a piecewise-constant step plot, the output of [Axes.Step]. The
// value Ys[i] is held from Xs[i] up to Xs[i+1] (a "post" step).
type StepSeries struct {
	// Xs and Ys are the point coordinates; they must have equal length.
	Xs, Ys []float64
	// Color is the line color.
	Color color.RGBA
	// Label is the legend label.
	Label string
	// Width is the line thickness in pixels (minimum 1).
	Width int
}

func (s *StepSeries) styleColor() color.RGBA { return s.Color }
func (s *StepSeries) legendLabel() string    { return s.Label }
func (s *StepSeries) bounds() (float64, float64, float64, float64, bool) {
	return xyBounds(s.Xs, s.Ys)
}

// --- FillBetween --------------------------------------------------------

// FillSeries is a shaded region between two y curves sharing common x values,
// the output of [Axes.FillBetween]. The band at Xs[i] spans Y1[i] to Y2[i].
type FillSeries struct {
	// Xs are the shared x coordinates.
	Xs []float64
	// Y1 and Y2 are the lower and upper curves; both match len(Xs).
	Y1, Y2 []float64
	// Color is the fill color; it is drawn semi-transparent.
	Color color.RGBA
	// Label is the legend label.
	Label string
}

func (s *FillSeries) styleColor() color.RGBA { return s.Color }
func (s *FillSeries) legendLabel() string    { return s.Label }
func (s *FillSeries) bounds() (float64, float64, float64, float64, bool) {
	if len(s.Xs) == 0 {
		return 0, 0, 0, 0, false
	}
	xmin, xmax := math.Inf(1), math.Inf(-1)
	ymin, ymax := math.Inf(1), math.Inf(-1)
	ok := false
	for i := range s.Xs {
		x := s.Xs[i]
		if isBad(x) {
			continue
		}
		for _, y := range []float64{s.Y1[i], s.Y2[i]} {
			if isBad(y) {
				continue
			}
			ok = true
			xmin, xmax = math.Min(xmin, x), math.Max(xmax, x)
			ymin, ymax = math.Min(ymin, y), math.Max(ymax, y)
		}
	}
	return xmin, xmax, ymin, ymax, ok
}

// --- plotting methods ---------------------------------------------------

// nextColor returns the automatic color for the next series based on the
// current series count.
func (a *Axes) nextColor() color.RGBA { return ColorCycle(len(a.series)) }

// Plot adds a line connecting the given points and returns the new series for
// further customization. The x and y slices must have equal length. The series
// is assigned the next color from the default cycle.
func (a *Axes) Plot(xs, ys []float64) *LineSeries {
	s := &LineSeries{Xs: xs, Ys: ys, Color: a.nextColor(), Width: 2}
	a.series = append(a.series, s)
	return s
}

// Scatter adds unconnected point markers and returns the new series. The x and
// y slices must have equal length.
func (a *Axes) Scatter(xs, ys []float64) *ScatterSeries {
	s := &ScatterSeries{Xs: xs, Ys: ys, Color: a.nextColor(), Size: 3}
	a.series = append(a.series, s)
	return s
}

// Bar adds a vertical bar chart with bars centred on xs and the given heights,
// rising from y = 0. If width <= 0 a width of 80% of the smallest center
// spacing is chosen automatically. It returns the new series.
func (a *Axes) Bar(xs, heights []float64, width float64) *BarSeries {
	if width <= 0 {
		width = autoBarWidth(xs)
	}
	s := &BarSeries{Xs: xs, Heights: heights, Width: width, Color: a.nextColor()}
	a.series = append(a.series, s)
	return s
}

// Hist bins the raw samples into the given number of equal-width bins spanning
// the sample range and adds the resulting histogram. If bins < 1 it defaults
// to 10. It returns the new series; with no finite samples the series is empty.
func (a *Axes) Hist(samples []float64, bins int) *HistSeries {
	if bins < 1 {
		bins = 10
	}
	lo, hi := math.Inf(1), math.Inf(-1)
	n := 0
	for _, v := range samples {
		if isBad(v) {
			continue
		}
		lo, hi = math.Min(lo, v), math.Max(hi, v)
		n++
	}
	s := &HistSeries{Color: a.nextColor()}
	if n == 0 || lo == hi {
		if n > 0 {
			// All identical: single unit-wide bin.
			s.Edges = []float64{lo - 0.5, lo + 0.5}
			s.Counts = []float64{float64(n)}
		}
		a.series = append(a.series, s)
		return s
	}
	edges := make([]float64, bins+1)
	step := (hi - lo) / float64(bins)
	for i := 0; i <= bins; i++ {
		edges[i] = lo + float64(i)*step
	}
	counts := make([]float64, bins)
	for _, v := range samples {
		if isBad(v) {
			continue
		}
		idx := int((v - lo) / step)
		if idx < 0 {
			idx = 0
		}
		if idx >= bins {
			idx = bins - 1
		}
		counts[idx]++
	}
	s.Edges, s.Counts = edges, counts
	a.series = append(a.series, s)
	return s
}

// Step adds a piecewise-constant step plot through the given points and returns
// the new series. The x and y slices must have equal length.
func (a *Axes) Step(xs, ys []float64) *StepSeries {
	s := &StepSeries{Xs: xs, Ys: ys, Color: a.nextColor(), Width: 2}
	a.series = append(a.series, s)
	return s
}

// FillBetween adds a shaded band between the two curves y1 and y2 sampled at
// the shared x values and returns the new series. All three slices must have
// equal length.
func (a *Axes) FillBetween(xs, y1, y2 []float64) *FillSeries {
	c := a.nextColor()
	c.A = 0x66
	s := &FillSeries{Xs: xs, Y1: y1, Y2: y2, Color: c}
	a.series = append(a.series, s)
	return s
}

// PlotFunc samples the Go function f at n evenly spaced points on [xmin, xmax]
// and adds the resulting line, returning the new series. Points where f returns
// a non-finite value are skipped, breaking the line there. If n < 2 it defaults
// to 200.
func (a *Axes) PlotFunc(f func(float64) float64, xmin, xmax float64, n int) *LineSeries {
	xs, ys := sampleFunc(f, xmin, xmax, n)
	return a.Plot(xs, ys)
}

// PlotExpr parses the algebra expression src, samples it at n evenly spaced
// points on [xmin, xmax] by binding the given variable name to each sample, and
// adds the resulting line. It returns the series and any parse error; on a
// parse error the returned series is nil. Points that fail to evaluate to a
// finite real number are skipped. If varName is empty it defaults to "x", and
// if n < 2 it defaults to 200.
//
// This is the only method that uses the parent algebra package; it lets a
// figure be driven directly from a symbolic expression such as "sin(x)/x".
func (a *Axes) PlotExpr(src, varName string, xmin, xmax float64, n int) (*LineSeries, error) {
	if varName == "" {
		varName = "x"
	}
	e, err := algebra.Parse(src)
	if err != nil {
		return nil, err
	}
	f := func(x float64) float64 {
		v, err := algebra.Eval(e, map[string]float64{varName: x})
		if err != nil {
			return math.NaN()
		}
		return v
	}
	s := a.PlotFunc(f, xmin, xmax, n)
	s.Label = src
	return s, nil
}

// --- configuration ------------------------------------------------------

// Title sets the axes title, drawn centred above the plot. It returns the axes
// for chaining.
func (a *Axes) Title(s string) *Axes { a.title = s; return a }

// XLabel sets the x-axis label, drawn centred below the plot. It returns the
// axes for chaining.
func (a *Axes) XLabel(s string) *Axes { a.xlabel = s; return a }

// YLabel sets the y-axis label, drawn rotated along the left edge. It returns
// the axes for chaining.
func (a *Axes) YLabel(s string) *Axes { a.ylabel = s; return a }

// Grid enables or disables the background grid lines at the tick locations. It
// returns the axes for chaining.
func (a *Axes) Grid(on bool) *Axes { a.grid = on; return a }

// Legend enables or disables the legend box, which lists every series that has
// a non-empty label. It returns the axes for chaining.
func (a *Axes) Legend(on bool) *Axes { a.legend = on; return a }

// XLim fixes the x-axis data range to [lo, hi], disabling autoscaling in x. It
// returns the axes for chaining.
func (a *Axes) XLim(lo, hi float64) *Axes {
	a.xlimSet, a.xmin, a.xmax = true, lo, hi
	return a
}

// YLim fixes the y-axis data range to [lo, hi], disabling autoscaling in y. It
// returns the axes for chaining.
func (a *Axes) YLim(lo, hi float64) *Axes {
	a.ylimSet, a.ymin, a.ymax = true, lo, hi
	return a
}

// dataRange returns the final x and y ranges used for rendering, honoring any
// limits set with XLim/YLim and otherwise autoscaling from the series data
// with a small margin. When there is no data it falls back to [0,1]x[0,1].
func (a *Axes) dataRange() (xmin, xmax, ymin, ymax float64) {
	dxmin, dxmax := math.Inf(1), math.Inf(-1)
	dymin, dymax := math.Inf(1), math.Inf(-1)
	any := false
	for _, s := range a.series {
		xl, xr, yl, yr, ok := s.bounds()
		if !ok {
			continue
		}
		any = true
		dxmin, dxmax = math.Min(dxmin, xl), math.Max(dxmax, xr)
		dymin, dymax = math.Min(dymin, yl), math.Max(dymax, yr)
	}
	if !any {
		dxmin, dxmax, dymin, dymax = 0, 1, 0, 1
	}
	xmin, xmax = dxmin, dxmax
	ymin, ymax = dymin, dymax
	if a.xlimSet {
		xmin, xmax = a.xmin, a.xmax
	} else {
		xmin, xmax = padRange(dxmin, dxmax, 0.03)
	}
	if a.ylimSet {
		ymin, ymax = a.ymin, a.ymax
	} else {
		ymin, ymax = padRange(dymin, dymax, 0.05)
	}
	if xmax <= xmin {
		xmin, xmax = xmin-0.5, xmin+0.5
	}
	if ymax <= ymin {
		ymin, ymax = ymin-0.5, ymin+0.5
	}
	return xmin, xmax, ymin, ymax
}

// --- helpers ------------------------------------------------------------

// isBad reports whether v is NaN or infinite.
func isBad(v float64) bool { return math.IsNaN(v) || math.IsInf(v, 0) }

// xyBounds returns the bounding box of the finite points in xs, ys.
func xyBounds(xs, ys []float64) (float64, float64, float64, float64, bool) {
	n := len(xs)
	if len(ys) < n {
		n = len(ys)
	}
	xmin, xmax := math.Inf(1), math.Inf(-1)
	ymin, ymax := math.Inf(1), math.Inf(-1)
	ok := false
	for i := 0; i < n; i++ {
		if isBad(xs[i]) || isBad(ys[i]) {
			continue
		}
		ok = true
		xmin, xmax = math.Min(xmin, xs[i]), math.Max(xmax, xs[i])
		ymin, ymax = math.Min(ymin, ys[i]), math.Max(ymax, ys[i])
	}
	return xmin, xmax, ymin, ymax, ok
}

// sampleFunc evaluates f at n evenly spaced points across [xmin, xmax].
func sampleFunc(f func(float64) float64, xmin, xmax float64, n int) ([]float64, []float64) {
	if n < 2 {
		n = 200
	}
	xs := make([]float64, n)
	ys := make([]float64, n)
	for i := 0; i < n; i++ {
		t := float64(i) / float64(n-1)
		x := xmin + t*(xmax-xmin)
		xs[i] = x
		ys[i] = f(x)
	}
	return xs, ys
}

// autoBarWidth chooses a bar width equal to 80% of the smallest gap between
// sorted, distinct centers, or 0.8 when fewer than two distinct centers exist.
func autoBarWidth(xs []float64) float64 {
	minGap := math.Inf(1)
	for i := 0; i < len(xs); i++ {
		for j := i + 1; j < len(xs); j++ {
			g := math.Abs(xs[i] - xs[j])
			if g > 0 && g < minGap {
				minGap = g
			}
		}
	}
	if math.IsInf(minGap, 1) {
		return 0.8
	}
	return minGap * 0.8
}
