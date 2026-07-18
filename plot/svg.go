package plot

import (
	"fmt"
	"math"
	"os"
	"strings"
)

// svgBuilder accumulates SVG markup.
type svgBuilder struct {
	b strings.Builder
}

func (s *svgBuilder) writef(format string, args ...interface{}) {
	fmt.Fprintf(&s.b, format, args...)
}

// escape replaces the XML metacharacters in text so it is safe inside SVG
// element content and attribute values.
func escape(text string) string {
	r := strings.NewReplacer(
		"&", "&amp;", "<", "&lt;", ">", "&gt;", "\"", "&quot;", "'", "&apos;",
	)
	return r.Replace(text)
}

// f formats a float for SVG coordinates with modest precision and no trailing
// zeros, keeping the output compact and deterministic.
func fnum(v float64) string {
	if math.IsNaN(v) || math.IsInf(v, 0) {
		v = 0
	}
	s := fmt.Sprintf("%.2f", v)
	s = strings.TrimRight(s, "0")
	s = strings.TrimRight(s, ".")
	if s == "" || s == "-0" {
		s = "0"
	}
	return s
}

// RenderSVG renders the figure and returns the SVG document as a string. The
// output is deterministic for a given figure.
func (f *Figure) RenderSVG() string {
	w, h := f.Width, f.Height
	s := &svgBuilder{}
	s.writef(`<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d" font-family="sans-serif">`, w, h, w, h)
	s.writef(`<rect x="0" y="0" width="%d" height="%d" fill="%s"/>`, w, h, hexColor(f.Background))

	a := f.ax
	tr := transform{
		px0: marginLeft, py0: marginTop,
		px1: float64(w - marginRight), py1: float64(h - marginBottom),
	}
	tr.xmin, tr.xmax, tr.ymin, tr.ymax = a.dataRange()

	xticks, xstep := ticksFor(tr.xmin, tr.xmax, 8)
	yticks, ystep := ticksFor(tr.ymin, tr.ymax, 6)

	if a.grid {
		for _, xv := range xticks {
			px := tr.X(xv)
			s.writef(`<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="1"/>`,
				fnum(px), fnum(tr.py0), fnum(px), fnum(tr.py1), hexColor(LightGray))
		}
		for _, yv := range yticks {
			py := tr.Y(yv)
			s.writef(`<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="1"/>`,
				fnum(tr.px0), fnum(py), fnum(tr.px1), fnum(py), hexColor(LightGray))
		}
	}

	// Clip series to the plot rectangle.
	s.writef(`<clipPath id="plotclip"><rect x="%s" y="%s" width="%s" height="%s"/></clipPath>`,
		fnum(tr.px0), fnum(tr.py0), fnum(tr.px1-tr.px0), fnum(tr.py1-tr.py0))
	s.writef(`<g clip-path="url(#plotclip)">`)
	for _, ser := range a.series {
		svgSeries(s, ser, tr)
	}
	s.writef(`</g>`)

	// Border.
	s.writef(`<rect x="%s" y="%s" width="%s" height="%s" fill="none" stroke="%s" stroke-width="1"/>`,
		fnum(tr.px0), fnum(tr.py0), fnum(tr.px1-tr.px0), fnum(tr.py1-tr.py0), hexColor(Black))

	// Ticks and labels.
	for _, xv := range xticks {
		px := tr.X(xv)
		s.writef(`<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="1"/>`,
			fnum(px), fnum(tr.py1), fnum(px), fnum(tr.py1+4), hexColor(Black))
		s.writef(`<text x="%s" y="%s" font-size="10" text-anchor="middle" fill="%s">%s</text>`,
			fnum(px), fnum(tr.py1+16), hexColor(Black), escape(formatTick(xv, xstep)))
	}
	for _, yv := range yticks {
		py := tr.Y(yv)
		s.writef(`<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="1"/>`,
			fnum(tr.px0-4), fnum(py), fnum(tr.px0), fnum(py), hexColor(Black))
		s.writef(`<text x="%s" y="%s" font-size="10" text-anchor="end" fill="%s">%s</text>`,
			fnum(tr.px0-7), fnum(py+3), hexColor(Black), escape(formatTick(yv, ystep)))
	}

	// Title and axis labels.
	if a.title != "" {
		s.writef(`<text x="%s" y="22" font-size="16" font-weight="bold" text-anchor="middle" fill="%s">%s</text>`,
			fnum(float64(w)/2), hexColor(Black), escape(a.title))
	}
	if a.xlabel != "" {
		s.writef(`<text x="%s" y="%d" font-size="12" text-anchor="middle" fill="%s">%s</text>`,
			fnum(float64(w)/2), h-8, hexColor(Black), escape(a.xlabel))
	}
	if a.ylabel != "" {
		cy := (tr.py0 + tr.py1) / 2
		s.writef(`<text x="16" y="%s" font-size="12" text-anchor="middle" fill="%s" transform="rotate(-90 16 %s)">%s</text>`,
			fnum(cy), hexColor(Black), fnum(cy), escape(a.ylabel))
	}

	if a.legend {
		svgLegend(s, a, tr)
	}

	s.writef(`</svg>`)
	return s.b.String()
}

// SaveSVG renders the figure and writes the SVG document to the file at path.
func (f *Figure) SaveSVG(path string) error {
	return os.WriteFile(path, []byte(f.RenderSVG()), 0o644)
}

// svgSeries appends the markup for one series to the builder.
func svgSeries(s *svgBuilder, ser series, tr transform) {
	switch v := ser.(type) {
	case *LineSeries:
		svgPolyline(s, v.Xs, v.Ys, v.Width, v.Color, tr, false)
	case *ScatterSeries:
		for i := range v.Xs {
			if i >= len(v.Ys) || isBad(v.Xs[i]) || isBad(v.Ys[i]) {
				continue
			}
			s.writef(`<circle cx="%s" cy="%s" r="%d" fill="%s"/>`,
				fnum(tr.X(v.Xs[i])), fnum(tr.Y(v.Ys[i])), v.Size, hexColor(v.Color))
		}
	case *BarSeries:
		for i := range v.Xs {
			if i >= len(v.Heights) || isBad(v.Xs[i]) || isBad(v.Heights[i]) {
				continue
			}
			svgBar(s, v.Xs[i]-v.Width/2, v.Xs[i]+v.Width/2, v.Baseline, v.Baseline+v.Heights[i], v.Color, tr)
		}
	case *HistSeries:
		for i := range v.Counts {
			svgBar(s, v.Edges[i], v.Edges[i+1], 0, v.Counts[i], v.Color, tr)
		}
	case *StepSeries:
		svgPolyline(s, v.Xs, v.Ys, v.Width, v.Color, tr, true)
	case *FillSeries:
		svgFill(s, v, tr)
	}
}

// svgPolyline appends a stroked path through the finite runs of the points. If
// step is true the path is drawn as a post-step (horizontal then vertical).
func svgPolyline(s *svgBuilder, xs, ys []float64, width int, col colorRGBA, tr transform, step bool) {
	if width < 1 {
		width = 1
	}
	n := len(xs)
	if len(ys) < n {
		n = len(ys)
	}
	var path strings.Builder
	prevOK := false
	var ppy float64
	for i := 0; i < n; i++ {
		if isBad(xs[i]) || isBad(ys[i]) {
			prevOK = false
			continue
		}
		px, py := tr.X(xs[i]), tr.Y(ys[i])
		if !prevOK {
			fmt.Fprintf(&path, "M%s %s", fnum(px), fnum(py))
		} else if step {
			fmt.Fprintf(&path, "L%s %s", fnum(px), fnum(ppy))
			fmt.Fprintf(&path, "L%s %s", fnum(px), fnum(py))
		} else {
			fmt.Fprintf(&path, "L%s %s", fnum(px), fnum(py))
		}
		ppy = py
		prevOK = true
	}
	if path.Len() == 0 {
		return
	}
	s.writef(`<path d="%s" fill="none" stroke="%s" stroke-width="%d" stroke-linejoin="round" stroke-linecap="round"/>`,
		path.String(), hexColor(col), width)
}

// svgBar appends a filled, outlined rectangle for one bar.
func svgBar(s *svgBuilder, xl, xr, yb, yt float64, col colorRGBA, tr transform) {
	x0, x1 := tr.X(xl), tr.X(xr)
	y0, y1 := tr.Y(yt), tr.Y(yb)
	if x1 < x0 {
		x0, x1 = x1, x0
	}
	if y1 < y0 {
		y0, y1 = y1, y0
	}
	edge := colorRGBA{R: col.R / 2, G: col.G / 2, B: col.B / 2, A: 0xff}
	s.writef(`<rect x="%s" y="%s" width="%s" height="%s" fill="%s" stroke="%s" stroke-width="1"/>`,
		fnum(x0), fnum(y0), fnum(x1-x0), fnum(y1-y0), hexColor(col), hexColor(edge))
}

// svgFill appends a filled polygon for the band between the two curves.
func svgFill(s *svgBuilder, v *FillSeries, tr transform) {
	var top, bottom []string
	for i := range v.Xs {
		if isBad(v.Xs[i]) || isBad(v.Y1[i]) || isBad(v.Y2[i]) {
			continue
		}
		top = append(top, fnum(tr.X(v.Xs[i]))+","+fnum(tr.Y(v.Y2[i])))
		bottom = append(bottom, fnum(tr.X(v.Xs[i]))+","+fnum(tr.Y(v.Y1[i])))
	}
	if len(top) == 0 {
		return
	}
	// Reverse the lower boundary to close the polygon.
	for i, j := 0, len(bottom)-1; i < j; i, j = i+1, j-1 {
		bottom[i], bottom[j] = bottom[j], bottom[i]
	}
	pts := strings.Join(append(top, bottom...), " ")
	op := float64(v.Color.A) / 255
	s.writef(`<polygon points="%s" fill="%s" fill-opacity="%s"/>`,
		pts, hexColor(colorRGBA{R: v.Color.R, G: v.Color.G, B: v.Color.B, A: 0xff}), fnum(op))
}

// svgLegend appends a boxed legend in the top-right corner.
func svgLegend(s *svgBuilder, a *Axes, tr transform) {
	type entry struct {
		label string
		col   colorRGBA
	}
	var entries []entry
	maxLbl := 0
	for _, ser := range a.series {
		lbl := ser.legendLabel()
		if lbl == "" {
			continue
		}
		entries = append(entries, entry{lbl, ser.styleColor()})
		if l := len([]rune(lbl)); l > maxLbl {
			maxLbl = l
		}
	}
	if len(entries) == 0 {
		return
	}
	rowH := 18.0
	boxW := 40.0 + float64(maxLbl)*7
	boxH := rowH*float64(len(entries)) + 8
	x1 := tr.px1 - 8
	x0 := x1 - boxW
	y0 := tr.py0 + 8
	s.writef(`<rect x="%s" y="%s" width="%s" height="%s" fill="%s" stroke="%s" stroke-width="1"/>`,
		fnum(x0), fnum(y0), fnum(boxW), fnum(boxH), hexColor(White), hexColor(Black))
	for i, e := range entries {
		ly := y0 + 8 + float64(i)*rowH
		col := e.col
		col.A = 0xff
		s.writef(`<rect x="%s" y="%s" width="20" height="4" fill="%s"/>`,
			fnum(x0+8), fnum(ly+4), hexColor(col))
		s.writef(`<text x="%s" y="%s" font-size="11" fill="%s">%s</text>`,
			fnum(x0+34), fnum(ly+9), hexColor(Black), escape(e.label))
	}
}
