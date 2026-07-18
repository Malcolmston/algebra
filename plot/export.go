package plot

import (
	"fmt"
	"image/color"
	"strconv"
	"strings"
)

// ExportMatplotlib returns a runnable Python script that reproduces the figure
// using the real Matplotlib library. The script imports matplotlib, rebuilds
// every series with the equivalent pyplot call, applies the title, axis labels,
// limits, grid and legend, and finally shows the figure.
//
// This lets the zero-dependency Go renderer act as a drop-in preview while
// still offering a path to pixel-perfect Matplotlib output: write the script to
// a .py file and run it with a Python interpreter that has matplotlib
// installed. The returned string always contains the line "import matplotlib".
func (f *Figure) ExportMatplotlib() string {
	a := f.ax
	var b strings.Builder

	b.WriteString("import matplotlib\n")
	b.WriteString("import matplotlib.pyplot as plt\n\n")
	fmt.Fprintf(&b, "fig, ax = plt.subplots(figsize=(%s, %s), dpi=100)\n\n",
		pyf(float64(f.Width)/100), pyf(float64(f.Height)/100))

	for i, s := range a.series {
		writeSeries(&b, s, i)
	}
	b.WriteString("\n")

	if a.title != "" {
		fmt.Fprintf(&b, "ax.set_title(%s)\n", pystr(a.title))
	}
	if a.xlabel != "" {
		fmt.Fprintf(&b, "ax.set_xlabel(%s)\n", pystr(a.xlabel))
	}
	if a.ylabel != "" {
		fmt.Fprintf(&b, "ax.set_ylabel(%s)\n", pystr(a.ylabel))
	}
	if a.xlimSet {
		fmt.Fprintf(&b, "ax.set_xlim(%s, %s)\n", pyf(a.xmin), pyf(a.xmax))
	}
	if a.ylimSet {
		fmt.Fprintf(&b, "ax.set_ylim(%s, %s)\n", pyf(a.ymin), pyf(a.ymax))
	}
	if a.grid {
		b.WriteString("ax.grid(True)\n")
	}
	if a.legend && hasLabels(a) {
		b.WriteString("ax.legend()\n")
	}
	b.WriteString("\nplt.tight_layout()\nplt.show()\n")
	return b.String()
}

// writeSeries emits the pyplot call for one series.
func writeSeries(b *strings.Builder, s series, idx int) {
	switch v := s.(type) {
	case *LineSeries:
		fmt.Fprintf(b, "ax.plot(%s, %s, color=%s%s)\n",
			pyList(v.Xs), pyList(v.Ys), pyColor(v.Color), labelKW(v.Label))
	case *ScatterSeries:
		fmt.Fprintf(b, "ax.scatter(%s, %s, color=%s, s=%d%s)\n",
			pyList(v.Xs), pyList(v.Ys), pyColor(v.Color), v.Size*v.Size*4, labelKW(v.Label))
	case *BarSeries:
		fmt.Fprintf(b, "ax.bar(%s, %s, width=%s, bottom=%s, color=%s%s)\n",
			pyList(v.Xs), pyList(v.Heights), pyf(v.Width), pyf(v.Baseline), pyColor(v.Color), labelKW(v.Label))
	case *HistSeries:
		// Reconstruct bar positions from bin edges and counts.
		centers := make([]float64, len(v.Counts))
		widths := make([]float64, len(v.Counts))
		for i := range v.Counts {
			centers[i] = (v.Edges[i] + v.Edges[i+1]) / 2
			widths[i] = v.Edges[i+1] - v.Edges[i]
		}
		w := 0.0
		if len(widths) > 0 {
			w = widths[0]
		}
		fmt.Fprintf(b, "ax.bar(%s, %s, width=%s, color=%s%s)  # histogram\n",
			pyList(centers), pyList(v.Counts), pyf(w), pyColor(v.Color), labelKW(v.Label))
	case *StepSeries:
		fmt.Fprintf(b, "ax.step(%s, %s, where='post', color=%s%s)\n",
			pyList(v.Xs), pyList(v.Ys), pyColor(v.Color), labelKW(v.Label))
	case *FillSeries:
		fmt.Fprintf(b, "ax.fill_between(%s, %s, %s, color=%s, alpha=%s%s)\n",
			pyList(v.Xs), pyList(v.Y1), pyList(v.Y2), pyColor(v.Color), pyf(float64(v.Color.A)/255), labelKW(v.Label))
	}
}

// hasLabels reports whether any series carries a legend label.
func hasLabels(a *Axes) bool {
	for _, s := range a.series {
		if s.legendLabel() != "" {
			return true
		}
	}
	return false
}

// labelKW returns a ", label=..." keyword fragment, or "" when label is empty.
func labelKW(label string) string {
	if label == "" {
		return ""
	}
	return ", label=" + pystr(label)
}

// pyColor renders an RGBA color as a Python hex string literal like "#1f77b4".
func pyColor(c color.RGBA) string { return pystr(hexColor(c)) }

// pyf formats a float64 as a Python numeric literal.
func pyf(v float64) string { return strconv.FormatFloat(v, 'g', -1, 64) }

// pyList formats a float slice as a Python list literal.
func pyList(xs []float64) string {
	parts := make([]string, len(xs))
	for i, v := range xs {
		parts[i] = pyf(v)
	}
	return "[" + strings.Join(parts, ", ") + "]"
}

// pystr formats a string as a single-quoted Python string literal with the
// necessary escapes.
func pystr(s string) string {
	r := strings.NewReplacer("\\", "\\\\", "'", "\\'", "\n", "\\n")
	return "'" + r.Replace(s) + "'"
}
