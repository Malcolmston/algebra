package plot

import (
	"bytes"
	"image/png"
	"math"
	"strings"
	"testing"
)

// pngSignature is the 8-byte magic number at the start of every PNG file.
var pngSignature = []byte{0x89, 'P', 'N', 'G', 0x0d, 0x0a, 0x1a, 0x0a}

func sampleFigure() *Figure {
	fig := NewWithSize(400, 300)
	ax := fig.Axes()
	ax.PlotFunc(math.Sin, 0, 2*math.Pi, 100).Label = "sin"
	ax.Title("Trig").XLabel("x").YLabel("y").Grid(true).Legend(true)
	return fig
}

func TestRenderPNGSignature(t *testing.T) {
	fig := sampleFigure()
	b, err := fig.RenderPNG()
	if err != nil {
		t.Fatalf("RenderPNG: %v", err)
	}
	if len(b) < 8 || !bytes.Equal(b[:8], pngSignature) {
		t.Fatalf("output does not start with PNG signature: % x", b[:min(8, len(b))])
	}
	// It must decode as a real PNG of the requested size.
	img, err := png.Decode(bytes.NewReader(b))
	if err != nil {
		t.Fatalf("png.Decode: %v", err)
	}
	if got := img.Bounds().Dx(); got != 400 {
		t.Fatalf("width = %d, want 400", got)
	}
	if got := img.Bounds().Dy(); got != 300 {
		t.Fatalf("height = %d, want 300", got)
	}
}

func TestRenderPNGDeterministic(t *testing.T) {
	a, err := sampleFigure().RenderPNG()
	if err != nil {
		t.Fatal(err)
	}
	b, err := sampleFigure().RenderPNG()
	if err != nil {
		t.Fatal(err)
	}
	if !bytes.Equal(a, b) {
		t.Fatal("PNG output is not deterministic")
	}
}

func TestRenderSVGContents(t *testing.T) {
	fig := sampleFigure()
	svg := fig.RenderSVG()
	for _, want := range []string{"<svg", "</svg>", "<path", "<text", "Trig", "sin"} {
		if !strings.Contains(svg, want) {
			t.Errorf("SVG missing %q", want)
		}
	}
}

func TestRenderSVGEscaping(t *testing.T) {
	fig := New()
	fig.Axes().Title("a < b & c")
	svg := fig.RenderSVG()
	if strings.Contains(svg, "a < b & c") {
		t.Error("SVG title was not XML-escaped")
	}
	if !strings.Contains(svg, "a &lt; b &amp; c") {
		t.Error("SVG escaped title not found")
	}
}

func TestExportMatplotlib(t *testing.T) {
	fig := NewWithSize(500, 400)
	ax := fig.Axes()
	ax.Plot([]float64{0, 1, 2}, []float64{0, 1, 4}).Label = "quad"
	ax.Scatter([]float64{0, 1}, []float64{1, 2})
	ax.Bar([]float64{0, 1, 2}, []float64{3, 1, 2}, 0)
	ax.Hist([]float64{1, 1, 2, 3, 3, 3}, 3)
	ax.Step([]float64{0, 1, 2}, []float64{1, 3, 2})
	ax.FillBetween([]float64{0, 1, 2}, []float64{0, 0, 0}, []float64{1, 2, 1})
	ax.Title("all").XLabel("x").YLabel("y").Grid(true).Legend(true).XLim(-1, 3).YLim(0, 5)

	py := fig.ExportMatplotlib()
	for _, want := range []string{
		"import matplotlib", "import matplotlib.pyplot as plt",
		"ax.plot(", "ax.scatter(", "ax.bar(", "ax.step(", "ax.fill_between(",
		"ax.set_title('all')", "ax.set_xlabel('x')", "ax.set_ylabel('y')",
		"ax.set_xlim(-1, 3)", "ax.set_ylim(0, 5)", "ax.grid(True)", "ax.legend()",
		"plt.show()",
	} {
		if !strings.Contains(py, want) {
			t.Errorf("export missing %q", want)
		}
	}
}

func TestPlotExpr(t *testing.T) {
	fig := New()
	s, err := fig.Axes().PlotExpr("x^2 + 1", "x", -2, 2, 5)
	if err != nil {
		t.Fatalf("PlotExpr: %v", err)
	}
	// Samples at x = -2,-1,0,1,2 => y = 5,2,1,2,5 (closed form x^2+1).
	wantX := []float64{-2, -1, 0, 1, 2}
	wantY := []float64{5, 2, 1, 2, 5}
	if len(s.Ys) != len(wantY) {
		t.Fatalf("got %d samples, want %d", len(s.Ys), len(wantY))
	}
	for i := range wantY {
		if math.Abs(s.Xs[i]-wantX[i]) > 1e-9 {
			t.Errorf("Xs[%d] = %v, want %v", i, s.Xs[i], wantX[i])
		}
		if math.Abs(s.Ys[i]-wantY[i]) > 1e-9 {
			t.Errorf("Ys[%d] = %v, want %v", i, s.Ys[i], wantY[i])
		}
	}
	if s.Label != "x^2 + 1" {
		t.Errorf("label = %q", s.Label)
	}
}

func TestPlotExprParseError(t *testing.T) {
	_, err := New().Axes().PlotExpr("x +* ", "x", 0, 1, 3)
	if err == nil {
		t.Fatal("expected parse error, got nil")
	}
}

func TestPlotFuncSampling(t *testing.T) {
	xs, ys := sampleFunc(func(x float64) float64 { return 2 * x }, 0, 10, 11)
	if len(xs) != 11 {
		t.Fatalf("len = %d, want 11", len(xs))
	}
	for i := range xs {
		if math.Abs(xs[i]-float64(i)) > 1e-9 {
			t.Errorf("xs[%d] = %v", i, xs[i])
		}
		if math.Abs(ys[i]-2*float64(i)) > 1e-9 {
			t.Errorf("ys[%d] = %v, want %v", i, ys[i], 2*float64(i))
		}
	}
}

func TestNiceNum(t *testing.T) {
	cases := []struct {
		x     float64
		round bool
		want  float64
	}{
		{10, false, 10},
		{9.9, false, 10},
		{5.1, false, 10},
		{0.9, false, 1},
		{100, true, 100},
		{1.4, true, 1},
		{2.9, true, 2},
		{6, true, 5},
	}
	for _, c := range cases {
		if got := niceNum(c.x, c.round); math.Abs(got-c.want) > 1e-9 {
			t.Errorf("niceNum(%v, %v) = %v, want %v", c.x, c.round, got, c.want)
		}
	}
}

func TestTicksFor(t *testing.T) {
	ticks, step := ticksFor(0, 10, 5)
	if math.Abs(step-2) > 1e-9 {
		t.Fatalf("step = %v, want 2", step)
	}
	want := []float64{0, 2, 4, 6, 8, 10}
	if len(ticks) != len(want) {
		t.Fatalf("got %d ticks %v, want %v", len(ticks), ticks, want)
	}
	for i := range want {
		if math.Abs(ticks[i]-want[i]) > 1e-9 {
			t.Errorf("tick[%d] = %v, want %v", i, ticks[i], want[i])
		}
	}
	// All ticks lie within the requested interval.
	for _, tk := range ticks {
		if tk < -1e-9 || tk > 10+1e-9 {
			t.Errorf("tick %v out of range", tk)
		}
	}
}

func TestFormatTick(t *testing.T) {
	cases := []struct {
		v, step float64
		want    string
	}{
		{0, 1, "0"},
		{5, 1, "5"},
		{0.5, 0.1, "0.5"},
		{-2, 1, "-2"},
	}
	for _, c := range cases {
		if got := formatTick(c.v, c.step); got != c.want {
			t.Errorf("formatTick(%v,%v) = %q, want %q", c.v, c.step, got, c.want)
		}
	}
}

func TestHistCounts(t *testing.T) {
	samples := []float64{0, 0, 1, 1, 1, 2}
	h := New().Axes().Hist(samples, 3)
	if len(h.Counts) != 3 {
		t.Fatalf("got %d bins, want 3", len(h.Counts))
	}
	total := 0.0
	for _, c := range h.Counts {
		total += c
	}
	if total != 6 {
		t.Fatalf("total count = %v, want 6", total)
	}
	// Last bin includes the maximum sample.
	if h.Counts[len(h.Counts)-1] < 1 {
		t.Errorf("max sample not counted, counts = %v", h.Counts)
	}
	// Edges span the sample range exactly.
	if math.Abs(h.Edges[0]-0) > 1e-9 || math.Abs(h.Edges[len(h.Edges)-1]-2) > 1e-9 {
		t.Errorf("edges = %v, want [0..2]", h.Edges)
	}
}

func TestAutoscale(t *testing.T) {
	fig := New()
	ax := fig.Axes()
	ax.Plot([]float64{0, 10}, []float64{-5, 5})
	xmin, xmax, ymin, ymax := ax.dataRange()
	// x padded by 3% of span (10) => [-0.3, 10.3].
	if xmin > 0 || xmax < 10 {
		t.Errorf("x range [%v,%v] does not contain data", xmin, xmax)
	}
	if ymin > -5 || ymax < 5 {
		t.Errorf("y range [%v,%v] does not contain data", ymin, ymax)
	}
}

func TestXYLimOverride(t *testing.T) {
	ax := New().Axes()
	ax.Plot([]float64{0, 100}, []float64{0, 100})
	ax.XLim(1, 2).YLim(3, 4)
	xmin, xmax, ymin, ymax := ax.dataRange()
	if xmin != 1 || xmax != 2 || ymin != 3 || ymax != 4 {
		t.Errorf("limits not honored: [%v,%v]x[%v,%v]", xmin, xmax, ymin, ymax)
	}
}

func TestColorCycle(t *testing.T) {
	if ColorCycle(0) != Blue {
		t.Error("ColorCycle(0) should be Blue")
	}
	if ColorCycle(len(DefaultColors)) != Blue {
		t.Error("ColorCycle should wrap around")
	}
	if ColorCycle(-1) != Cyan {
		t.Error("ColorCycle(-1) should wrap to last color")
	}
}

func TestSeriesAutoColor(t *testing.T) {
	ax := New().Axes()
	l0 := ax.Plot([]float64{0}, []float64{0})
	l1 := ax.Plot([]float64{0}, []float64{0})
	if l0.Color != DefaultColors[0] || l1.Color != DefaultColors[1] {
		t.Error("successive series should get successive palette colors")
	}
}

func TestEmptyFigureRenders(t *testing.T) {
	fig := New()
	b, err := fig.RenderPNG()
	if err != nil {
		t.Fatalf("RenderPNG on empty figure: %v", err)
	}
	if !bytes.Equal(b[:8], pngSignature) {
		t.Fatal("empty figure did not produce a PNG")
	}
	if !strings.Contains(fig.RenderSVG(), "<svg") {
		t.Fatal("empty figure did not produce SVG")
	}
}

func TestNonFiniteBreaksLine(t *testing.T) {
	// A NaN in the middle should not crash rendering.
	fig := New()
	fig.Axes().Plot([]float64{0, 1, 2, 3}, []float64{0, math.NaN(), 2, 3})
	if _, err := fig.RenderPNG(); err != nil {
		t.Fatalf("render with NaN: %v", err)
	}
}

func TestSaveFiles(t *testing.T) {
	dir := t.TempDir()
	fig := sampleFigure()
	if err := fig.SavePNG(dir + "/out.png"); err != nil {
		t.Fatalf("SavePNG: %v", err)
	}
	if err := fig.SaveSVG(dir + "/out.svg"); err != nil {
		t.Fatalf("SaveSVG: %v", err)
	}
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}
